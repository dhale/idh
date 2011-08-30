/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

import dnp.LocalSlopeFinder;
import het.RecursiveExponentialFilter;

/*
Algorithm:
given image f
find slopes p2,p3
for all phis
  rotate images f,p2,p3 to align with axis 2
  for all i3
    use slopes p2,p3 to make plus-minus images fp and fm
    compute fmp, fmm, and fpp
    smooth2 fmp, fmm and fpp
  for all i2
    copy fmp, fmm, and fpp to 2D arrays with axes 13
    for all thetas
      shear3 fmp, fmm, and fpp to align with axis 1
      smooth1 products fmp, fmm and fpp
      c = (<fmp>*<fmp>)/(<fmm><fpp>), if <fmp> > 0; 0, otherwise
      map c to fault likelihood c = 1-c
      unshear3 image c
      remember cmax and corresponding tmax (theta max)
  unrotate cmax and tmax
  remember cmax, tmax and pmax (phi max)
thin cmax by smoothing laterally and picking peaks
use p2,p3,cmax for structure-oriented smoothing
for all fault surfaces
  gather image samples alongside fault
  cross-correlate to find displacements
*/

/**
 * Finds faults in 3D seismic images.
 * <p>
 * This version computes fault likelihoods from correlation coefficients,
 * which are ratios of smoothed correlation products, and it uses rotation 
 * and shearing to enable that smoothing to be performed with axis-aligned
 * recursive filters.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.06.24
 */
public class FaultFinder3 {

  /**
   * Constructs a fault finder with specified parameters.
   * @param slopeMax maximum slope of seismic reflections.
   * @param shiftMax maximum fault shift, in samples.
   * @param thetaMax maximum fault angle theta, in degrees.
   */
  public FaultFinder3(double slopeMax, double shiftMax, double thetaMax) {
    _slopeMax = (float)slopeMax;
    _shiftMax = (float)shiftMax;
    _thetaMax = (float)thetaMax;
    _faultLengthMin = 2.0f*_shiftMax;
    _sigma = 2.0f*_shiftMax;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    _st  = makeThetaSampling(-_thetaMax,_thetaMax);
    _sp = makePhiSampling(-90.0,90.0);
  }

  public void setPhiSampling(Sampling sp) {
    _sp = sp;
  }

  public void setThetaSampling(Sampling st) {
    _st = st;
  }

  public Sampling makePhiSampling(double phiMin, double phiMax) {
    return angleSampling(SIGMA2C,phiMin,phiMax);
  }

  public Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma,thetaMin,thetaMax);
  }

  public float[][][][] scan(float[][][] p2, float[][][] p3, float[][][] f) {
    final int n1 = n1(f);
    final int n2 = n2(f);
    final int n3 = n3(f);
    final float[][][] c = new float[n3][n2][n1];
    final float[][][] p = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    int nphi = _sp.getCount();
    for (int iphi=0; iphi<nphi; ++iphi) {
      final float phi = (float)_sp.getValue(iphi);
      Rotator r = new Rotator(phi,n1,n2,n3);
      float[][][] rf = r.rotate(f);
      float[][][] rp2 = r.rotate(p2);
      float[][][] rp3 = r.rotate(p3);
      float[][][][] fmp = align(phi,rp2,rp3,rf);
      final int n2r = n2(rf);
      final int n3r = n3(rf);
      final float[][][] fm = fmp[0];
      final float[][][] fp = fmp[1];
      final float[][][] cmp = new float[n3r][n2r][];
      final float[][][] cmm = fm;
      final float[][][] cpp = fp;
      loop(n3r,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2r; ++i2) {
          if (fm[i3][i2]!=null) {
            cmp[i3][i2] = new float[n1];
            float[] fm32 = fm[i3][i2];
            float[] fp32 = fp[i3][i2];
            float[] cmp32 = cmp[i3][i2];
            float[] cmm32 = cmm[i3][i2];
            float[] cpp32 = cpp[i3][i2];
            for (int i1=0; i1<n1; ++i1) {
              float fmi = fm32[i1];
              float fpi = fp32[i1];
              cmp32[i1] = fmi*fpi;
              cmm32[i1] = fmi*fmi;
              cpp32[i1] = fpi*fpi;
            }
          }
        }
      }});
      smooth2(SIGMA2C,cmp);
      smooth2(SIGMA2C,cmm);
      smooth2(SIGMA2C,cpp);
      float[][][][] cp = {cmp,cmm,cpp};
      float[][][][] ct = ctScan(cp);
      float[][][] cr = ct[0];
      float[][][] tr = ct[1];
      final float[][][] ci = r.unrotate(cr);
      final float[][][] ti = r.unrotate(tr);
      loop(n3,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] c32 = c[i3][i2];
          float[] p32 = p[i3][i2];
          float[] t32 = t[i3][i2];
          float[] ci32 = ci[i3][i2];
          float[] ti32 = ti[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float cii = ci32[i1];
            float tii = ti32[i1];
            if (cii>c32[i1]) {
              c32[i1] = cii;
              t32[i1] = tii;
              p32[i1] = phi;
            }
          }
        }
      }});
    }
    return new float[][][][]{c,p,t};
  }

  public float[][][][] findSlopes(float[][][] f) {
    int n1 = n1(f);
    int n2 = n2(f);
    int n3 = n3(f);
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(SIGMA1P,_slopeMax);
    lsf.findSlopes(f,p2,p3,null);
    return new float[][][][]{p2,p3};
  }

  public static float[][][] taper(int m, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] g = copy(f);
    float[] t = new float[m];
    for (int i=0; i<m; ++i) {
      t[i] = (float)(0.54+0.46*cos(PI*(m-i)/m));
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0,j1=n1-1; i1<m; ++i1,--j1) {
          float ti = t[i1];
          g[i3][i2][i1] *= ti;
          g[i3][i2][j1] *= ti;
        }
      }
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0,j2=n2-1; i2<m; ++i2,--j2) {
        float ti = t[i2];
        for (int i1=0; i1<n1; ++i1) {
          g[i3][i2][i1] *= ti;
          g[i3][j2][i1] *= ti;
        }
      }
    }
    for (int i3=0,j3=n3-1; i3<m; ++i3,--j3) {
      float ti = t[i3];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          g[i3][i2][i1] *= ti;
          g[j3][i2][i1] *= ti;
        }
      }
    }
    return g;
  }

  public static void nullsToZeros(float[][][] f) {
    int n1 = n1(f);
    int n2 = n2(f);
    int n3 = n3(f);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        if (f[i3][i2]==null)
          f[i3][i2] = new float[n1];
      }
    }
  }

  public static void zerosToNulls(float[][][] f) {
    int n1 = n1(f);
    int n2 = n2(f);
    int n3 = n3(f);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] f32 = f[i3][i2];
        boolean zeros = true;
        for (int i1=0; i1<n1 && zeros; ++i1) {
          if (f32[i1]!=0.0f)
            zeros = false;
        }
        if (zeros)
          f[i3][i2] = null;
      }
    }
  }

  public float[][][][] align(
    double phi, float[][][] p2, float[][][] p3, float[][][] f) 
  {
    int n1 = n1(f);
    int n2 = n2(f);
    int n3 = n3(f);
    float[][][] fm = new float[n3][n2][];
    float[][][] fp = new float[n3][n2][];
    float[] xm = new float[n1];
    float[] xp = new float[n1];
    float phir = (float)toRadians(phi);
    float cosp = cos(phir);
    float sinp = sin(phir);
    _si.setUniformSampling(n1,1.0,0.0);
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        float[] p2m = p2[i3m][i2];
        float[] p2p = p2[i3p][i2];
        float[] p3m = p3[i3m][i2];
        float[] p3p = p3[i3p][i2];
        float[] f2m = f[i3m][i2];
        float[] f2p = f[i3p][i2];
        if (p2m!=null && p2p!=null && 
            p3m!=null && p3p!=null &&
            f2m!=null && f2p!=null) {
          for (int i1=0; i1<n1; ++i1) {
            // for phi =   0, p3 =  p3 (no change)
            // for phi =  90, p3 = -p2
            // for phi = -90, p3 =  p2
            float p3mi = cosp*p3m[i1]-sinp*p2m[i1];
            float p3pi = cosp*p3p[i1]-sinp*p2p[i1];
            xm[i1] = i1-p3mi;
            xp[i1] = i1+p3pi;
          }
          fm[i3][i2] = new float[n1];
          fp[i3][i2] = new float[n1];
          _si.setUniformSamples(f2m);
          _si.interpolate(n1,xm,fm[i3][i2]);
          _si.setUniformSamples(f2p);
          _si.interpolate(n1,xp,fp[i3][i2]);
        }
      }
    }
    return new float[][][][]{fm,fp};
  }
  public float[][][] xcor(float theta, float[][][][] fmp) {
    // This version performs only one unshear, but computes correlation
    // coefficients in the sheared space that has more samples.
    theta = toRadians(theta);
    float shear = tan(theta);
    float sigma = _sigma*cos(theta);
    final RecursiveExponentialFilter ref = 
      new RecursiveExponentialFilter(sigma,SIGMA2C);
    final int n1 = n1(fmp[0]);
    final int n2 = n2(fmp[0]);
    final int n3 = n3(fmp[0]);
    float[][][] fm = fmp[0];
    float[][][] fp = fmp[1];
    final float[][][] sm = shear(shear,fm);
    final float[][][] sp = shear(shear,fp);
    int n3s = sm.length;
    final float[][][] c = sm;
    loop(n3s,new LoopInt() {
      public void compute(int i3) {
        float[][] cmp = new float[n2][n1];
        float[][] cmm = new float[n2][n1];
        float[][] cpp = new float[n2][n1];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float smi = sm[i3][i2][i1];
            float spi = sp[i3][i2][i1];
            cmp[i2][i1] = smi*spi;
            cmm[i2][i1] = smi*smi;
            cpp[i2][i1] = spi*spi;
          }
        }
        ref.apply(cmp,cmp);
        ref.apply(cmm,cmm);
        ref.apply(cpp,cpp);
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float cmpi = cmp[i2][i1];
            float cmmi = cmm[i2][i1];
            float cppi = cpp[i2][i1];
            if (cmpi<0.0f) {
              c[i3][i2][i1] = 0.0f;
            } else {
              float cnum = cmpi*cmpi;
              float cden = cmmi*cppi;
              c[i3][i2][i1] = (cnum<cden)?cnum/cden:1.0f;
            }
          }
        }
      }
    });

    float[][][] cc = unshear(shear,c);
    clip(0.0f,1.0f,cc,cc);
    return cc;
  }

  // public for development only
  public static class Rotator {

    public Rotator(double phi, int n1, int n2, int n3) {
      _n1 = n1;

      // angle phi in radians, cosine and sine
      _phir = toRadians(phi);
      _cosp = cos(_phir);
      _sinp = sin(_phir);

      // center of rotation
      _x2c = 0.5*(n2-1.0);
      _x3c = 0.5*(n3-1.0);

      // input sampling
      _s2p = new Sampling(n2,1.0,0.0);
      _s3p = new Sampling(n3,1.0,0.0);

      // corners of input sampling rectangle
      double[] x2s = { 0.0, 0.0,n2-1,n2-1};
      double[] x3s = { 0.0,n3-1,n3-1, 0.0};

      // bounds after rotation
      double x2min =  Double.MAX_VALUE;
      double x3min =  Double.MAX_VALUE;
      double x2max = -Double.MAX_VALUE;
      double x3max = -Double.MAX_VALUE;
      for (int i=0; i<4; ++i) {
        double x2q = x2q(x2s[i],x3s[i]);
        double x3q = x3q(x2s[i],x3s[i]);
        if (x2q<x2min) x2min = x2q;
        if (x2q>x2max) x2max = x2q;
        if (x3q<x3min) x3min = x3q;
        if (x3q>x3max) x3max = x3q;
      }
      x2min = floor(x2min);
      x2max = ceil(x2max);
      x3min = floor(x3min);
      x3max = ceil(x3max);

      // sampling after rotation
      int n2q = max(2,1+(int)(x2max-x2min+0.5));
      int n3q = max(2,1+(int)(x3max-x3min+0.5));
      double d2q = 1.0;
      double d3q = 1.0;
      double f2q = x2min;
      double f3q = x3min;
      _s2q = new Sampling(n2q,d2q,f2q);
      _s3q = new Sampling(n3q,d3q,f3q);
      //System.out.println("s2p: n2p="+n2);
      //System.out.println("s3p: n3p="+n3);
      //System.out.println("s2q: n2q="+n2q+" d2q="+d2q+" f2q="+f2q);
      //System.out.println("s3q: n3q="+n3q+" d3q="+d3q+" f3q="+f3q);
    }

    public float[][][] rotate(float[][][] p) {
      final float[][][] fp = p;
      final float[][] siTable = _siTable;
      final int nsinc = siTable.length;
      final int lsinc = siTable[0].length;
      final Sampling s2p = _s2p;
      final Sampling s3p = _s3p;
      final Sampling s2q = _s2q;
      final Sampling s3q = _s3q;
      final int n1 = _n1;
      final int n2p = _s2p.getCount();
      final int n3p = _s3p.getCount();
      final int n2q = _s2q.getCount();
      final int n3q = _s3q.getCount();
      final float[][][] q = new float[n3q][n2q][];
      loop(n3q,new LoopInt() {
        public void compute(int i3) {
          double x3q = s3q.getValue(i3);
          for (int i2=0; i2<n2q; ++i2) {
            double x2q = s2q.getValue(i2);
            double x2p = x2p(x2q,x3q);
            double x3p = x3p(x2q,x3q);
            if (inBounds(x2p,x3p)) {
              float[] q32 = q[i3][i2] = new float[n1];
              int i2p = (int)floor(x2p);
              int i3p = (int)floor(x3p);
              double f2p = x2p-i2p;
              double f3p = x3p-i3p;
              int k2p = (int)(f2p*(nsinc-1)+0.5);
              int k3p = (int)(f3p*(nsinc-1)+0.5);
              for (int k3s=0; k3s<lsinc; ++k3s) {
                float s3 = siTable[k3p][k3s];
                int j3p = i3p+k3s-lsinc/2+1;
                if (j3p<   0) j3p = 0;
                if (j3p>=n3p) j3p = n3p-1;
                for (int k2s=0; k2s<lsinc; ++k2s) {
                  float s2 = siTable[k2p][k2s];
                  int j2p = i2p+k2s-lsinc/2+1;
                  if (j2p<   0) j2p = 0;
                  if (j2p>=n2p) j2p = n2p-1;
                  float[] p32 = fp[j3p][j2p];
                  float s32 = s3*s2;
                  for (int i1=0; i1<n1; ++i1)
                    q32[i1] += p32[i1]*s32;
                }
              }
            }
          }
        }
      });
      return q;
    }

    public float[][][] unrotate(float[][][] q) {
      final float[][][] fq = q;
      final float[][] siTable = _siTable;
      final int nsinc = siTable.length;
      final int lsinc = siTable[0].length;
      final Sampling s2p = _s2p;
      final Sampling s3p = _s3p;
      final Sampling s2q = _s2q;
      final Sampling s3q = _s3q;
      final int n1 = _n1;
      final int n2p = s2p.getCount();
      final int n3p = s3p.getCount();
      final int n2q = s2q.getCount();
      final int n3q = s3q.getCount();
      //System.out.println("n2p="+n2p+" n3p="+n3p+" n2q="+n2q+" n3q="+n3q);
      final double d2q = s2q.getDelta();
      final double d3q = s3q.getDelta();
      final double f2q = s2q.getFirst();
      final double f3q = s3q.getFirst();
      final float[][][] p = new float[n3p][n2p][n1];
      loop(n3p,new LoopInt() {
        public void compute(int i3) {
          double x3p = s3p.getValue(i3);
          for (int i2=0; i2<n2p; ++i2) {
            float[] p32 = p[i3][i2];
            double x2p = s2p.getValue(i2);
            double x2q = x2q(x2p,x3p);
            double x3q = x3q(x2p,x3p);
            double y2q = (x2q-f2q)/d2q;
            double y3q = (x3q-f3q)/d3q;
            int i2q = (int)floor(y2q);
            int i3q = (int)floor(y3q);
            double e2q = y2q-i2q;
            double e3q = y3q-i3q;
            int k2q = (int)(e2q*(nsinc-1)+0.5);
            int k3q = (int)(e3q*(nsinc-1)+0.5);
            for (int k3s=0; k3s<lsinc; ++k3s) {
              float s3 = siTable[k3q][k3s];
              int j3q = i3q+k3s-lsinc/2+1;
              if (j3q<   0) j3q = 0;
              if (j3q>=n3q) j3q = n3q-1;
              for (int k2s=0; k2s<lsinc; ++k2s) {
                float s2 = siTable[k2q][k2s];
                int j2q = i2q+k2s-lsinc/2+1;
                if (j2q<   0) j2q = 0;
                if (j2q>=n2q) j2q = n2q-1;
                float[] q32 = fq[j3q][j2q];
                if (q32!=null) {
                  float s32 = s3*s2;
                  for (int i1=0; i1<n1; ++i1)
                    p32[i1] += q32[i1]*s32;
                }
              }
            }
          }
        }
      });
      return p;
    }

    private int _n1; // number of samples in 1st dimension
    private double _phir,_cosp,_sinp; // angle phi in radians, cosine, sine
    private double _x2c,_x3c; // coordinates of center of rotation
    private Sampling _s2p,_s3p; // samplings in original coordinates
    private Sampling _s2q,_s3q; // samplings in rotated coordinates
    private static float[][] _siTable; // sinc interpolation coefficients
    static {
      SincInterpolator si = new SincInterpolator();
      _siTable = si.getTable();
    }
    private double x2p(double x2q, double x3q) {
      return _x2c+(x2q-_x2c)*_cosp-(x3q-_x3c)*_sinp;
    }
    private double x3p(double x2q, double x3q) {
      return _x3c+(x2q-_x2c)*_sinp+(x3q-_x3c)*_cosp;
    }
    private double x2q(double x2p, double x3p) {
      return _x2c+(x2p-_x2c)*_cosp+(x3p-_x3c)*_sinp;
    }
    private double x3q(double x2p, double x3p) {
      return _x3c-(x2p-_x2c)*_sinp+(x3p-_x3c)*_cosp;
    }
    private boolean inBounds(double x2p, double x3p) {
      return _s2p.getFirst()<=x2p && x2p<=_s2p.getLast() &&
             _s3p.getFirst()<=x3p && x3p<=_s3p.getLast();
    }
  }

  /**
   * Shears an image horizontally with q(i1,i2,i3) = p(i1,i2,i3+s*i1).
   * For non-zero shears, the number of samples n2q in the 2nd dimension 
   * of the output image will exceed the number of input samples n2p.
   * @param s the shear.
   * @param p the input image to shear.
   * @return the output sheared image q.
   */
  public float[][][] shear(final double s, final float[][][] p) {
    final int n1 = n1(p);
    final int n2 = n2(p);
    final int n3p = n3(p);
    final int n3q = n3p+(int)(abs(s)*n1);
    final double dqp = n3q-n3p;
    final float[][][] q = new float[n3q][n2][n1]; 
    final Unsafe<SincInterpolator> usi = new Unsafe<SincInterpolator>();
    loop(n2,new LoopInt() {
      public void compute(int i2) {
        SincInterpolator si = usi.get();
        if (si==null) {
          si = new SincInterpolator();
          si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
          usi.set(si);
        }
        float[] pp = new float[n3p];
        float[] qq = new float[n3q];
        si.setUniform(n3p,1.0,0.0,pp);
        for (int i1=0; i1<n1; ++i1) {
          for (int i3=0; i3<n3p; ++i3)
            pp[i3] = p[i3][i2][i1];
          double f3q = (s<0.0f)?s*i1:s*i1-dqp;
          si.interpolate(n3q,1.0f,f3q,qq);
          for (int i3=0; i3<n3q; ++i3)
            q[i3][i2][i1] = qq[i3];
        }
      }
    });
    return q;
  }
 
  /**
   * Unshears an image horizontally with p(i1,i2,i3) = q(i1,i2,i3-s*i1).
   * Except for interpolation errors, this method is the inverse of the
   * method {@link #shear(double,float[][][])}.
   * @param s the shear.
   * @param q the input image to unshear.
   * @return the output unsheared image p.
   */
  public float[][][] unshear(final double s, final float[][][] q) {
    final int n1 = n1(q);
    final int n2 = n2(q);
    final int n3q = n3(q);
    final int n3p = n3q-(int)(abs(s)*n1);
    final double dqp = n3q-n3p;
    final float[][][] p = new float[n3p][n2][n1]; 
    final Unsafe<SincInterpolator> usi = new Unsafe<SincInterpolator>();
    loop(n2,new LoopInt() {
      public void compute(int i2) {
        SincInterpolator si = usi.get();
        if (si==null) {
          si = new SincInterpolator();
          si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
          usi.set(si);
        }
        float[] pp = new float[n3p];
        float[] qq = new float[n3q];
        si.setUniform(n3q,1.0,0.0,qq);
        for (int i1=0; i1<n1; ++i1) {
          for (int i3=0; i3<n3q; ++i3)
            qq[i3] = q[i3][i2][i1];
          double f3p = (s<0.0f)?-s*i1:-s*i1+dqp;
          si.interpolate(n3p,1.0f,f3p,pp);
          for (int i3=0; i3<n3p; ++i3)
            p[i3][i2][i1] = pp[i3];
        }
      }
    });
    return p;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final double SIGMA1P = 8.0; // for estimating slopes
  private static final double SIGMA2C = 4.0; // for correlation coeff

  private float _sigma; // for smoothing correlation products vertically 
  private float _slopeMax; // maximum slope used in alignment
  private float _thetaMax; // maximum fault angle theta
  private float _shiftMax; // maximum vertical throw for any fault
  private float _faultLengthMin; // min number of samples in a fault
  private SincInterpolator _si; // for shifting, squeezing and stretching
  private Sampling _st; // sampling of fault angle theta
  private Sampling _sp; // sampling of fault strike phi

  // A range of non-null arrays within a 3D array.
  private static class Range {
    public int first; // index of first non-null array
    public int count; // number of non-null arrays
    private Range(int first, int count) {
      this.first = first;
      this.count = count;
    }
    public static Range get2(int i3, float[][][] f) {
      int n2 = n2(f);
      int j2 = -1;
      for (int i2=0; i2<n2 && j2<0; ++i2) {
        if (f[i3][i2]!=null)
          j2 = i2;
      }
      int k2 = -1;
      for (int i2=n2-1; i2>=0 && k2<0; --i2) {
        if (f[i3][i2]!=null)
          k2 = i2;
      }
      int m2 = (j2>=0)?1+k2-j2:0;
      return new Range(j2,m2);
    }
    public static Range get3(int i2, float[][][] f) {
      int n3 = n3(f);
      int j3 = -1;
      for (int i3=0; i3<n3 && j3<0; ++i3) {
        if (f[i3][i2]!=null)
          j3 = i3;
      }
      int k3 = -1;
      for (int i3=n3-1; i3>=0 && k3<0; --i3) {
        if (f[i3][i2]!=null)
          k3 = i3;
      }
      int m3 = (j3>=0)?1+k3-j3:0;
      return new Range(j3,m3);
    }
  }

  // Samples on opposite sides of a fault.
  private static class FaultSamples {
    public static final int KA = 2; // offset used to find fault shifts
    public static final int KB = KA+1; // offset used with KA for slopes
    public int[] k1,k2; // sample indices of fault
    public float[] fmb,fma,fpa,fpb; // image samples at -ka-1, -ka, ka, ka+1
    public FaultSamples(
      int[] k1, int[] k2,
      float[] fmb, float[] fma, float[] fpa, float[] fpb)
    {
      this.k1 = k1;
      this.k2 = k2;
      this.fmb = fmb;
      this.fma = fma;
      this.fpa = fpa;
      this.fpb = fpb;
    }
  }

  private FaultSamples findFaultSamples(
    int i1, int i2, float[][] c, float[][] f, boolean[][] b) 
  {
    if (c[i2][i1]==0.0f)
      return null;
    int n1 = c[0].length;
    int n2 = c.length;
    int k2min = 0;
    int k2max = n2-1;
    int[] k2s = new int[n1];
    int n = 0;
    boolean done = false;
    for (int k1=i1,k2=i2; !done && k1<n1; ++k1) {
      if (c[k2][k1]>0.0f) {
        b[k2][k1] = true;
        k2s[n++] = k2;
      } else if (k2<k2max && c[k2+1][k1]>0.0f) {
        ++k2;
        b[k2][k1] = true;
        k2s[n++] = k2;
      } else if (k2>k2min && c[k2-1][k1]>0.0f) {
        --k2;
        b[k2][k1] = true;
        k2s[n++] = k2;
      } else {
        done = true;
      }
    }
    if (n<_faultLengthMin)
      return null;
    int[] k1 = new int[n];
    int[] k2 = new int[n];
    float[] fmb = new float[n];
    float[] fma = new float[n];
    float[] fpa = new float[n];
    float[] fpb = new float[n];
    int ka = FaultSamples.KA;
    int kb = FaultSamples.KB;
    for (int k=0; k<n; ++k) {
      int k1k = i1+k;
      int k2k = k2s[k];
      int kmb = max(k2min,k2k-kb);
      int kma = max(k2min,k2k-ka);
      int kpa = min(k2max,k2k+ka);
      int kpb = min(k2max,k2k+kb);
      fmb[k] = f[kmb][k1k];
      fma[k] = f[kma][k1k];
      fpa[k] = f[kpa][k1k];
      fpb[k] = f[kpb][k1k];
      k1[k] = k1k;
      k2[k] = k2k;
    }
    return new FaultSamples(k1,k2,fmb,fma,fpa,fpb);
  }

  // Get numbers of samples in 3D arrays.
  private static int n1(float[][][] f) {
    int n1 = 0;
    int n2 = f[0].length;
    int n3 = f.length;
    for (int i3=0; i3<n3 && n1==0; ++i3) {
      for (int i2=0; i2<n2 && n1==0; ++i2) {
        if (f[i3][i2]!=null)
          n1 = f[i3][i2].length;
      }
    }
    return n1;
  }
  private static int n2(float[][][] f) {
    return f[0].length;
  }
  private static int n3(float[][][] f) {
    return f.length;
  }

  // Smoothing along 2nd dimension of the specified 3D array.
  private static void smooth2(double sigma, final float[][][] f) {
    final int n1 = n1(f);
    final int n2 = n2(f);
    final int n3 = n3(f);
    final RecursiveExponentialFilter ref = 
      new RecursiveExponentialFilter(sigma);
    loop(n3,new LoopInt() { 
    public void compute(int i3) {
      Range r = Range.get2(i3,f);
      int j2 = r.first;
      int m2 = r.count;
      if (m2>0) {
        float[][] f3 = new float[m2][];
        for (int i2=0; i2<m2; ++i2)
          f3[i2] = f[i3][j2+i2];
        ref.apply2(f3,f3);
      }
    }});
  }

  // Scan for maximum fault likelihood c and corresponding theta.
  private float[][][][] ctScan(float[][][][] cp) {
    final float[][][] cmp = cp[0];
    final float[][][] cmm = cp[1];
    final float[][][] cpp = cp[2];
    final int n1 = n1(cmp);
    final int n2 = n2(cmp);
    final int n3 = n3(cmp);
    final float[][][] c = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    final Unsafe<SincInterpolator> usi = new Unsafe<SincInterpolator>();
    for (int i2=0; i2<n2; ++i2) {
    //loop(n2,new LoopInt() {
    //public void compute(int i2) {
      SincInterpolator si = usi.get();
      if (si==null) {
        si = new SincInterpolator();
        si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
        usi.set(si);
        System.out.println("si="+si);
      }
      Range r = Range.get3(i2,cmp);
      int j3 = r.first;
      int m3 = r.count;
      if (m3>0) {
        float[][] cmp2 = new float[m3][];
        float[][] cmm2 = new float[m3][];
        float[][] cpp2 = new float[m3][];
        for (int i3=0; i3<m3; ++i3) {
          cmp2[i3] = cmp[i3+j3][i2];
          cmm2[i3] = cmm[i3+j3][i2];
          cpp2[i3] = cpp[i3+j3][i2];
        }
        float[][][] cp2 = {cmp2,cmm2,cpp2};
        int nt = _st.getCount();
        for (int it=0; it<nt; ++it) {
          float ti = (float)_st.getValue(it);
          float cscale = (it==0 || it==nt-1)?0.5f:1.0f;
          float[][] cx = xcorForOneTheta(si,ti,cp2);
          for (int i3=0; i3<m3; ++i3) {
            float[] cx3 = cx[i3];
            float[] c32 = c[i3+j3][i2];
            float[] t32 = t[i3+j3][i2];
            for (int i1=0; i1<n1; ++i1) {
              float ci = (1.0f-cx3[i1])*cscale;
              if (ci>c32[i1]) {
                c32[i1] = ci;
                t32[i1] = ti;
              }
            }
          }
        }
      }
    //}});
    }
    return new float[][][][]{c,t};
  }

  // Cross-correlation for specified theta and correlation products.
  private float[][] xcorForOneTheta(
    SincInterpolator si, float theta, float[][][] cp) 
  {
    theta = toRadians(theta);
    float shear = tan(theta);
    float sigma = _sigma*cos(theta);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    int n1 = cp[0][0].length;
    int n2 = cp[0].length;
    float[][][] sp = new float[3][][];
    for (int i=0; i<3; ++i)
      sp[i] = shear(si,shear,cp[i]);
    int n2s = sp[0].length;
    float[][] c = new float[n2s][n1];
    float[] cmp = new float[n1];
    float[] cmm = new float[n1];
    float[] cpp = new float[n1];
    for (int i2=0; i2<n2s; ++i2) {
      float[] smp = sp[0][i2];
      float[] smm = sp[1][i2];
      float[] spp = sp[2][i2];
      ref.apply(smp,cmp);
      ref.apply(smm,cmm);
      ref.apply(spp,cpp);
      for (int i1=0; i1<n1; ++i1) {
        float cmpi = cmp[i1];
        float cmmi = cmm[i1];
        float cppi = cpp[i1];
        if (cmpi<0.0f) {
          c[i2][i1] = 0.0f;
        } else {
          float cnum = cmpi*cmpi;
          float cden = cmmi*cppi;
          c[i2][i1] = (cnum<cden)?cnum/cden:1.0f;
        }
      }
    }
    c = unshear(si,shear,c);
    clip(0.0f,1.0f,c,c);
    return c;
  }

  private Sampling angleSampling(double sigma, double amin, double amax) {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  private static float[][] shear(
    SincInterpolator si, double s, float[][] p) 
  {
    int n1 = p[0].length;
    int n2p = p.length;
    int n2q = n2p+(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] q = new float[n2q][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    si.setUniform(n2p,1.0,0.0,pp);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2p; ++i2)
        pp[i2] = p[i2][i1];
      double f2q = (s<0.0f)?s*i1:s*i1-dqp;
      si.interpolate(n2q,1.0f,f2q,qq);
      for (int i2=0; i2<n2q; ++i2)
        q[i2][i1] = qq[i2];
    }
    return q;
  }
  private static float[][] unshear(
    SincInterpolator si, double s, float[][] q) 
  {
    int n1 = q[0].length;
    int n2q = q.length;
    int n2p = n2q-(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] p = new float[n2p][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    si.setUniform(n2q,1.0,0.0,qq);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2q; ++i2)
        qq[i2] = q[i2][i1];
      double f2p = (s<0.0f)?-s*i1:-s*i1+dqp;
      si.interpolate(n2p,1.0f,f2p,pp);
      for (int i2=0; i2<n2p; ++i2)
        p[i2][i1] = pp[i2];
    }
    return p;
  }
}
