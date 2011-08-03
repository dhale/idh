/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

import dnp.LocalSlopeFinder;
import het.RecursiveExponentialFilter;

/**
 * Finds faults in 3D seismic images.
 * <p>
 * This version computes fault likelihoods from semblances, and uses 
 * rotated Gaussian filters implemented with FFTs to perform smoothing
 * of semblance numerators and denominators for different fault angles.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.06.26
 */
public class FaultFinder3SB {

  /**
   * Constructs a fault finder with specified parameters.
   * @param slopeMax maximum slope of image features.
   * @param shiftMax maximum fault shift, in samples.
   * @param thetaMax maximum fault angle, in degrees.
   */
  public FaultFinder3SB(double slopeMax, double shiftMax, double thetaMax) {
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
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
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
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
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
   * Returns semblances computed for specified fault strike and angle.
   * @param phi fault strike, in degrees.
   * @param theta fault angle, in degrees.
   * @param snd array {snum,sden} of semblance numerators and denominators.
   * @return array of semblances.
   */
  public float[][][] semblance(double phi, double theta, float[][][][] snd) {
    final int n1 = snd[0][0][0].length;
    final int n2 = snd[0][0].length;
    final int n3 = snd[0].length;
    //FaultPlaneSmoother fps = new FaultPlaneSmoother(_sigma,4.0,0.0,snd);
    FaultPlaneSmoother fps = new FaultPlaneSmoother(4.0,0.0,0.0,snd);
    snd = new float[2][n3][n2][n1];
    fps.apply(phi,theta,snd);
    final float[][][] sn = snd[0];
    final float[][][] sd = snd[1];
    final float[][][] s = sn;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] s32 = s[i3][i2];
        float[] sn32 = sn[i3][i2];
        float[] sd32 = sd[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float sni = sn32[i1];
          float sdi = sd32[i1];
          if (sdi<=0.0f || sni<0.0f) {
            s32[i1] = 0.0f;
          } else if (sdi<sni) {
            s32[i1] = 1.0f;
          } else {
            s32[i1] = sni/sdi;
          }
        }
      }
    }});
    return s;
  }

  /**
   * Returns fault likelihoods and corresponding fault angles.
   * @param snd array {snum,sden} of semblance numerators and denominators.
   * @return array {c,p,t} of fault likelihoods c and fault angles p and t.
   */
  public float[][][][] faultPhiThetaScan(float[][][][] snd) {
    final int n1 = snd[0][0][0].length;
    final int n2 = snd[0][0].length;
    final int n3 = snd[0].length;
    final float[][][] c = new float[n3][n2][n1];
    final float[][][] p = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    int np = _sp.getCount();
    int nt = _st.getCount();
    for (int ip=0; ip<np; ++ip) {
      final float pi = (float)_sp.getValue(ip);
      for (int it=0; it<nt; ++it) {
        final float ti = (float)_st.getValue(it);
        final float[][][] s = semblance(pi,ti,snd);
        loop(n3,new LoopInt() {
        public void compute(int i3) {
          for (int i2=0; i2<n2; ++i2) {
            float[] s32 = s[i3][i2];
            float[] c32 = c[i3][i2];
            float[] p32 = p[i3][i2];
            float[] t32 = t[i3][i2];
            for (int i1=0; i1<n1; ++i1) {
              float si = s32[i1]; // semblance
              si = si*si; // semblance^2
              si = si*si; // semblance^4
              si = si*si; // semblance^8
              float ci = 1.0f-si;
              if (ci>c32[i1]) {
                c32[i1] = ci;
                p32[i1] = pi;
                t32[i1] = ti;
              }
            }
          }
        }});
      }
    }
    return new float[][][][]{c,p,t};
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
