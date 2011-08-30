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

import edu.mines.jtk.util.Stopwatch;
import dnp.LocalSlopeFinder;

import edu.mines.jtk.mosaic.*;

/*
Algorithm:
given image f
find slopes p2,p3
compute semblance snum,sden (2)
initialize c,p,t (3)

for all phi:
  rotate snum,sden (2)
  smooth snum,sden
  compute cphi,tphi (2)
  unrotate cphi,tphi
  update c,p,t
Total memory is 9 3D volumes

For FFT smoothing we need
FFT of snum,sden (2)
smoothed snum,sden (2)
c,p,t (3)
Total memory is 7 3D volumes

for all phi:
  rotate (parallel i3) snum,sden to align fault plane with axis 2
  for all i3 (parallel i3)
    extract non-null i3 slices from rotated snum,sden
    smooth snum,sden slices horizontally along axis 2
    restore non-null i3 slices into rotated snum,sden
  initialize cphi,tphi in rotated coordinates
  for all i2 (parallel i2):
    extract non-null i2 slices of rotated snum,sden
    for all theta:
      shear snum,sden slices to align fault plane with axis 1
      smooth snum,sden slices vertically along axis 1
      compute semblance s
      unshear semblance s
      compute fault likelihood c from s
      use c update cphi,tphi to maximize cphi
  unrotate (parallel i3) cphi,tphi
  for all i3 (parallel i3)
    update c,p,t
thin c,p,t by smoothing laterally and picking peaks
use p2,p3,c for structure-oriented smoothing
for all fault surfaces
  gather image samples alongside fault
  cross-correlate to find displacements
*/

/**
 * Finds faults in 3D seismic images.
 * <p>
 * This version computes fault likelihoods from semblances, using
 * image rotation and shearing to perform axis-aligned smoothing of 
 * semblance numerators and denominators for different fault angles.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.08.12
 */
public class FaultFinder3S {

  /**
   * Constructs a fault finder with specified parameters.
   * @param slopeMax maximum slope of image features.
   * @param shiftMax maximum fault shift, in samples.
   * @param thetaMax maximum fault angle, in degrees.
   */
  public FaultFinder3S(double slopeMax, double shiftMax, double thetaMax) {
    _slopeMax = (float)slopeMax;
    _shiftMax = (float)shiftMax;
    _thetaMax = (float)thetaMax;
    _faultLengthMin = 2.0f*_shiftMax;
    _sigma = 2.0f*_shiftMax;
    _rgf1 = new RecursiveGaussianFilter(1.0f);
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
    return angleSampling(4.0,phiMin,phiMax);
  }

  public Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma,thetaMin,thetaMax);
  }

  /**
   * Returns slopes of locally planar features in the specified image.
   * @param f the input image.
   * @return array {p2,p3} of slopes.
   */
  public float[][][][] findSlopes(float[][][] f) {
    int n1 = n1(f), n2 = n2(f), n3 = n3(f);
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,_slopeMax);
    lsf.findSlopes(f,p2,p3,null);
    return new float[][][][]{p2,p3};
  }

  /**
   * Returns an image with specified samples taper to zero at edges.
   * Tapering enables simple zero-value boundary conditions to be
   * used in fault detection without artifacts caused by abrupt 
   * truncations of strong image features near edges.
   * @param m width of the tapered band of samples at each edge.
   * @param f input image.
   * @return the tapered image.
   */
  public static float[][][] taper(int m, float[][][] f) {
    int n1 = n1(f), n2 = n2(f), n3 = n3(f);
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

  /**
   * Returns semblance numerators and denominators.
   * Each numerator is a squared average of image values, and each
   * denominator is an average of squared values. Specified slopes 
   * are used to align image samples before this averaging.
   * @param p input array of slopes.
   * @param f input image.
   * @return array {snum,sden} of semblance numerators and denominators.
   */
  public float[][][][] semblanceNumDen(
    final float[][][][] p, final float[][][] f) 
  {
    final int n1 = n1(f), n2 = n2(f), n3 = n3(f);
    final float[][][] sn = new float[n3][n2][n1];
    final float[][][] sd = new float[n3][n2][n1];
    final float[][][] p2 = p[0];
    final float[][][] p3 = p[1];
    final Unsafe<SincInterpolator> usi = new Unsafe<SincInterpolator>();
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      SincInterpolator si = usi.get();
      if (si==null) {
        si = new SincInterpolator();
        si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
        usi.set(si);
      }
      si.setUniformSampling(n1,1.0,0.0);
      float[] xmm = new float[n1];
      float[] xm0 = new float[n1];
      float[] xmp = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xpm = new float[n1];
      float[] xp0 = new float[n1];
      float[] xpp = new float[n1];
      float[] gmm = new float[n1];
      float[] gm0 = new float[n1];
      float[] gmp = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gpm = new float[n1];
      float[] gp0 = new float[n1];
      float[] gpp = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p2mm = p2[i3m][i2m];
        float[] p2m0 = p2[i3m][i2 ];
        float[] p2mp = p2[i3m][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p200 = p2[i3 ][i2 ];
        float[] p20p = p2[i3 ][i2p];
        float[] p2pm = p2[i3p][i2m];
        float[] p2p0 = p2[i3p][i2 ];
        float[] p2pp = p2[i3p][i2p];
        float[] p3mm = p3[i3m][i2m];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3mp = p3[i3m][i2p];
        float[] p30m = p3[i3 ][i2m];
        float[] p300 = p3[i3 ][i2 ];
        float[] p30p = p3[i3 ][i2p];
        float[] p3pm = p3[i3p][i2m];
        float[] p3p0 = p3[i3p][i2 ];
        float[] p3pp = p3[i3p][i2p];
        float[] sn32 = sn[i3][i2];
        float[] sd32 = sd[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xmm[i1] = i1-p3mm[i1]-p2mm[i1];
          xm0[i1] = i1-p3m0[i1]         ;
          xmp[i1] = i1-p3mp[i1]+p2mp[i1];
          x0m[i1] = i1         -p20m[i1];
          x0p[i1] = i1         +p20p[i1];
          xpm[i1] = i1+p3pm[i1]-p2pm[i1];
          xp0[i1] = i1+p3p0[i1]         ;
          xpp[i1] = i1+p3pp[i1]+p2pp[i1];
        }
        si.setUniformSamples(fmm); si.interpolate(n1,xmm,gmm);
        si.setUniformSamples(fm0); si.interpolate(n1,xm0,gm0);
        si.setUniformSamples(fmp); si.interpolate(n1,xmp,gmp);
        si.setUniformSamples(f0m); si.interpolate(n1,x0m,g0m);
        si.setUniformSamples(f0p); si.interpolate(n1,x0p,g0p);
        si.setUniformSamples(fpm); si.interpolate(n1,xpm,gpm);
        si.setUniformSamples(fp0); si.interpolate(n1,xp0,gp0);
        si.setUniformSamples(fpp); si.interpolate(n1,xpp,gpp);
        for (int i1=0; i1<n1; ++i1) {
          float gmmi = gmm[i1];
          float gm0i = gm0[i1];
          float gmpi = gmp[i1];
          float g0mi = g0m[i1];
          float g00i = f00[i1];
          float g0pi = g0p[i1];
          float gpmi = gpm[i1];
          float gp0i = gp0[i1];
          float gppi = gpp[i1];
          float sumn = gmmi+gm0i+gmpi+
                       g0mi+g00i+g0pi+
                       gpmi+gp0i+gppi;
          float sumd = gmmi*gmmi+gm0i*gm0i+gmpi*gmpi+
                       g0mi*g0mi+g00i*g00i+g0pi*g0pi+
                       gpmi*gpmi+gp0i*gp0i+gppi*gppi;
          sn32[i1] = sumn*sumn;
          sd32[i1] = 9.0f*sumd;
        }
      }
    }});
    return new float[][][][]{sn,sd};
  }

  /**
   * Returns fault likelihoods and corresponding fault angles.
   * @param snd array {snum,sden} of semblance numerators and denominators.
   * @return array {c,p,t} of fault likelihoods c and fault angles p and t.
   */
  public float[][][][] faultPhiThetaScan(float[][][][] snd) {
    final int n1 = n1(snd), n2 = n2(snd), n3 = n3(snd);
    final float[][][] c = new float[n3][n2][n1];
    final float[][][] p = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    int np = _sp.getCount();
    for (int ip=0; ip<np; ++ip) {
      System.out.println("ip="+ip);
      final float phi = (float)_sp.getValue(ip);
      Rotator r = new Rotator(phi,n1,n2,n3);
      float[][][][] rsnd = r.rotate(snd);
      smooth2(rsnd);
      float[][][][] rctp = scanTheta(rsnd);
      float[][][][] ctp = r.unrotate(rctp);
      final float[][][] cp = ctp[0];
      final float[][][] tp = ctp[1];
      loop(n3,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] c32 = c[i3][i2];
          float[] p32 = p[i3][i2];
          float[] t32 = t[i3][i2];
          float[] cp32 = cp[i3][i2];
          float[] tp32 = tp[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float cpi = cp32[i1];
            if (cpi<0.0f) cpi = 0.0f;
            if (cpi>1.0f) cpi = 1.0f;
            if (cpi>c32[i1]) {
              c32[i1] = cpi;
              p32[i1] = phi;
              t32[i1] = tp32[i1];
            }
          }
        }
      }});
    }
    return new float[][][][]{c,p,t};
  }

  // Horizontal smoothing of rotated snum,sden along axis 2.
  private void smooth2(float[][][][] snd) {
    final float[][][] snum = snd[0];
    final float[][][] sden = snd[1];
    final int n1 = n1(snum);
    final int n2 = n2(snum);
    final int n3 = n3(snum);
    final RecursiveExponentialFilter ref = new RecursiveExponentialFilter(4.0);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      float[][] s3;
      s3 = getSlice3(i3,snum); 
      if (s3!=null) {
        ref.apply2(s3,s3); 
        setSlice3(i3,snum,s3);
      }
      s3 = getSlice3(i3,sden); 
      if (s3!=null) {
        ref.apply2(s3,s3); 
        setSlice3(i3,sden,s3);
      }
    }});
  }

  // Scan over fault angles theta to maximize fault likelihood.
  private float[][][][] scanTheta(float[][][][] snd) {
    final float[][][][] ct = likeRotated(snd);
    final float[][][] sn = snd[0];
    final float[][][] sd = snd[1];
    final float[][][] c = ct[0];
    final float[][][] t = ct[1];
    final int n1 = n1(c);
    final int n2 = n2(c);
    final int n3 = n3(c);
    final Unsafe<SincInterpolator> usi = new Unsafe<SincInterpolator>();
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      SincInterpolator si = usi.get();
      if (si==null) {
        usi.set(si=new SincInterpolator());
        si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
      }
      float[][] sn2 = getSlice2(i2,sn);
      float[][] sd2 = getSlice2(i2,sd);
      if (sn2==null || sd2==null)
        return;
      int n3 = sn2.length;
      int nt = _st.getCount();
      for (int it=0; it<nt; ++it) {
        float ti = (float)_st.getValue(it);
        float theta = toRadians(ti);
        float shear = tan(theta);
        float sigma = _sigma*cos(theta);
        float[][] sns = shear(si,shear,sn2);
        float[][] sds = shear(si,shear,sd2);
        RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
        ref.apply1(sns,sns);
        ref.apply1(sds,sds);
        float[][] ss = semblance(sns,sds);
        float[][] s2 = unshear(si,shear,ss);
        for (int i3=0,j3=i3lo(i2,c); i3<n3; ++i3,++j3) {
          float[] s32 = s2[i3];
          float[] c32 = c[j3][i2];
          float[] t32 = t[j3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float st = s32[i1]; // semblance
            st = st*st; // semblance^2
            st = st*st; // semblance^4
            st = st*st; // semblance^8
            float ci = 1.0f-st;
            if (ci>c32[i1]) {
              c32[i1] = ci;
              t32[i1] = ti;
            }
          }
        }
      }
    }});
    return ct;
  }
  private float[][][][] likeRotated(float[][][][] p) {
    int n1 = n1(p[0]);
    int n2 = n2(p[0]);
    int n3 = n3(p[0]);
    int np = p.length;
    float[][][][] q = new float[np][n3][n2][];
    for (int ip=0; ip<np; ++ip) {
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          q[ip][i3][i2] = (p[ip][i3][i2]!=null)?new float[n1]:null;
        }
      }
    }
    return q;
  }
  private float[][] semblance(float[][] sn, float[][] sd) {
    int n1 = sn[0].length;
    int n2 = sn.length;
    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] s2 = s[i2];
      float[] sn2 = sn[i2];
      float[] sd2 = sd[i2];
      for (int i1=0; i1<n1; ++i1) {
        float sni = sn2[i1];
        float sdi = sd2[i1];
        if (sdi<=0.0f || sni<0.0f) {
          s2[i1] = 0.0f;
        } else if (sdi<sni) {
          s2[i1] = 1.0f;
        } else {
          s2[i1] = sni/sdi;
        }
      }
    }
    return s;
  }

  /**
   * Returns thinned fault likelihoods and corresponding fault angles.
   * @param ct array {c,t} of input fault likelihoods c and fault angles t.
   * @return array {c,t} of thinned fault likelihoods c and fault angles t.
   */
  /*
  public float[][][] faultThetaThin(float[][][] ct) {
    float[][] c = ct[0];
    float[][] t = ct[1];
    int n1 = c[0].length;
    int n2 = c.length;
    //pow(c,4.00f,c);
    //_rgf1.apply00(c,c);
    //pow(c,0.25f,c);
    float dt = (float)_st.getDelta();
    float ft = (float)_st.getFirst();
    float lt = (float)_st.getLast();
    float tmin = ft-0.5f*dt;
    float tmax = lt+0.5f*dt;
    float[][] cc = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
      float[] ti = t[i2  ];
      float[] ci = c[i2  ];
      float[] cm = c[i2-1];
      float[] cp = c[i2+1];
      for (int i1=0; i1<n1; ++i1) {
        float cii = ci[i1];
        float tii = ti[i1];
        if (cm[i1]<cii && cp[i1]<cii && tmin<tii && tii<tmax) {
          cc[i2][i1] = cii;
          tt[i2][i1] = tii;
        }
      }
    }
    return new float[][][]{cc,tt};
  }
  */

  /**
   * Returns a smoothed image for specified slopes and fault likelihoods.
   * The smoothing is structure-oriented and limited near fault locations.
   * @param sigma half-width of smoothing filter where no faults are present.
   * @param p array of slopes.
   * @param c array of fault likelihoods.
   * @param f array of input image samples.
   * @return array smoothed output image samples.
   */
  /*
  public float[][] smooth(
    double sigma, float[][] p, float[][] c, float[][] f) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    EigenTensors2 d = new EigenTensors2(n1,n2);
    d.setEigenvalues(0.001f,1.00f); // smooth mostly along structure
    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        s[i2][i1] = 1.0f-pow(c[i2][i1],0.1f); // almost a binary image
        float dip = atan(p[i2][i1]);
        float u1 =  cos(dip);
        float u2 = -sin(dip);
        d.setEigenvectorU(i1,i2,u1,u2);
      }
    }
    float a = (float)(0.5*sigma*sigma);
    float[][] g = new float[n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,a,s,f,g);
    return g;
  }
  */

  /**
   * Returns vertical fault displacements for specified (thinned) faults.
   * @param ct array of thinned fault likelihoods and angles.
   * @param f array of input image samples.
   * @return array of vertical fault displacements.
   */
  /*
  public float[][] findShifts(float[][][] ct, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] c = ct[0];
    //float[][] t = ct[1]; // currently not used
    float[][] u = new float[n2][n1];
    boolean[][] b = new boolean[n2][n1];
    int ka = FaultSamples.KA;
    int kb = FaultSamples.KB;
    int uMinShift = -(int)_shiftMax;
    int uMaxShift =  (int)_shiftMax;
    int pMinShift = -(int)_slopeMax*kb;
    int pMaxShift =  (int)_slopeMax*kb;
    LocalShiftFinderX lsfu = new LocalShiftFinderX(_sigma);
    LocalShiftFinderX lsfp = new LocalShiftFinderX(8.0);
    lsfu.setSmoothShifts(true);
    lsfp.setSmoothShifts(false);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (!b[i2][i1] && c[i2][i1]>0.0f) {
          FaultSamples fs = findFaultSamples(i1,i2,c,f,b);
          if (fs!=null) {
            int[] k1 = fs.k1;
            int[] k2 = fs.k2;
            float[] fmb = fs.fmb;
            float[] fma = fs.fma;
            float[] fpa = fs.fpa;
            float[] fpb = fs.fpb;
            int n = k1.length;
            float[] pm = new float[n];
            float[] pp = new float[n];
            float[] uk = new float[n];
            lsfp.find1(pMinShift,pMaxShift,fmb,fma,pm);
            lsfp.find1(pMinShift,pMaxShift,fpa,fpb,pp);
            lsfu.find1(uMinShift,uMaxShift,fma,fpa,uk);
            if (n>1500) {
              System.out.println("i1="+i1+" i2="+i2+" n="+n);
              SimplePlot sp = new SimplePlot();
              sp.setTitle("i1="+i1+" i2="+i2);
              //sp.addPoints(fm2).setLineColor(Color.GREEN);
              sp.addPoints(fma).setLineColor(Color.RED);
              sp.addPoints(fpa).setLineColor(Color.BLUE);
              sp = new SimplePlot();
              sp.setTitle("i1="+i1+" i2="+i2);
              sp.addPoints(uk).setLineColor(Color.RED);
              //sp.addPoints(fp2).setLineColor(Color.MAGENTA);
            }
            float scale = (float)ka/(float)(kb-ka);
            for (int k=0; k<n; ++k)
              u[k2[k]][k1[k]] = uk[k]-scale*(pm[k]+pp[k]);
          }
        }
      }
    }
    return u;
  }
  */

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma;
  private float _slopeMax;
  private float _thetaMax;
  private float _shiftMax;
  private float _faultLengthMin;
  private RecursiveGaussianFilter _rgf1;
  private SincInterpolator _si;
  private Sampling _sp,_st;

  /**
   * Returns a sampling to avoid aliasing for a range of angles.
   * @param sigma half-width of Gaussian smoothing filter.
   * @param amin minimum angle to be sampled.
   * @param amax maximum angle to be sampled.
   * @return the angle sampling.
   */
  private Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  // Get numbers of samples in 3D and 4D arrays that may have nulls.
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
  private static int n1(float[][][][] f) {
    return n1(f[0]);
  }
  private static int n2(float[][][][] f) {
    return n2(f[0]);
  }
  private static int n3(float[][][][] f) {
    return n2(f[0]);
  }

  // Get/set non-null slices of rotated 3D arrays
  private static float[][] getSlice2(int i2, float[][][] x) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i3lo = i3lo(i2,x);
    int i3hi = i3hi(i2,x);
    int m3 = 1+i3hi-i3lo;
    float[][] x2 = (m3>0)?new float[m3][n1]:null;
    for (int i3=0; i3<m3; ++i3)
      copy(x[i3+i3lo][i2],x2[i3]);
    return x2;
  }
  private static void setSlice2(int i2, float[][][] x, float[][] x2) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i3lo = i3lo(i2,x);
    int i3hi = i3hi(i2,x);
    int m3 = 1+i3hi-i3lo;
    assert x2.length==m3:"x2 length is correct";
    for (int i3=0; i3<m3; ++i3)
      copy(x2[i3],x[i3+i3lo][i2]);
  }
  private static float[][] getSlice3(int i3, float[][][] x) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i2lo = i2lo(i3,x);
    int i2hi = i2hi(i3,x);
    int m2 = 1+i2hi-i2lo;
    float[][] x3 = (m2>0)?new float[m2][n1]:null;
    for (int i2=0; i2<m2; ++i2)
      copy(x[i3][i2+i2lo],x3[i2]);
    return x3;
  }
  private static void setSlice3(int i3, float[][][] x, float[][] x3) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i2lo = i2lo(i3,x);
    int i2hi = i2hi(i3,x);
    int m2 = 1+i2hi-i2lo;
    assert x3.length==m2:"x3 length is correct";
    for (int i2=0; i2<m2; ++i2)
      copy(x3[i2],x[i3][i2+i2lo]);
  }
  private static int i2lo(int i3, float[][][] x) {
    int n2 = x[0].length;
    int i2lo = 0;
    while (i2lo<n2 && x[i3][i2lo]==null)
      ++i2lo;
    return i2lo;
  }
  private static int i2hi(int i3, float[][][] x) {
    int n2 = x[0].length;
    int i2hi = n2-1;
    while (i2hi>=0 && x[i3][i2hi]==null)
      --i2hi;
    return i2hi;
  }
  private static int i3lo(int i2, float[][][] x) {
    int n3 = x.length;
    int i3lo = 0;
    while (i3lo<n3 && x[i3lo][i2]==null)
      ++i3lo;
    return i3lo;
  }
  private static int i3hi(int i2, float[][][] x) {
    int n3 = x.length;
    int i3hi = n3-1;
    while (i3hi>=0 && x[i3hi][i2]==null)
      --i3hi;
    return i3hi;
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

    public float[][][][] rotate(float[][][][] p) {
      int n = p.length;
      float[][][][] q = new float[n][][][];
      for (int i=0; i<n; ++i)
        q[i] = rotate(p[i]);
      return q;
    }

    public float[][][][] unrotate(float[][][][] p) {
      int n = p.length;
      float[][][][] q = new float[n][][][];
      for (int i=0; i<n; ++i)
        q[i] = unrotate(p[i]);
      return q;
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

  // Shear horizontally such that q(i1,i2) = p(i1,i2+s*i1).
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

  // Unshear horizontally such that p(i1,i2) = q(i1,i2-s*i1).
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

  /**
   * Samples collected from both sides of a fault.
   */
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
}
