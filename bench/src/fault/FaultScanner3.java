/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import java.util.ArrayList;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

import edu.mines.jtk.util.Stopwatch;
import dnp.LocalSlopeFinder;

import edu.mines.jtk.mosaic.*;

/**
 * Computes fault likelihoods by scanning over fault orientations.
 * Fault likelihoods are in the range [0,1], where 0 and 1 denote 
 * lowest and highest likelihoods, respectively.
 * <p>
 * A fault scanner is constructed from precomputed semblance numerators 
 * and denominators. During a scan over fault orientations, these images 
 * are smoothed within fault planes having different orientations, before 
 * computing the semblance ratios. Fault likelihoods are computed from 
 * semblances. For each sample location, the fault likelihood equals the
 * maximum likelihood computed at that location during a scan over all
 * orientations.
 * <p>
 * Scan results include the fault strike and dip angle for which the
 * maximum likelihoods occured. Both angles are measured in degrees.
 * The fault strike angle in [-90,90] is measured relative to the 
 * horizontal 2nd image axis, and is positive when the horizontal 
 * strike vector has a positive component in the direction of the 
 * 3rd image axis. The fault dip angle, also in [-90,90], is measured 
 * from the vertical 1st image axis, and is positive when the horizontal 
 * projection of the dip vector has a positive component in the direction
 * of the 2nd image axis. Note that fault dip is here defined such that 
 * vertical faults have zero dip.
 * <p>
 * EXPERIMENTAL: this class provides two methods for smoothing
 * semblance numerators and denominators along fault planes.
 * One method first rotates and shears those factors before 
 * applying fast recursive axis-aligned smoothing filters.
 * The other method uses rotated Gaussian filters implemented
 * with FFTs.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.08.13
 */
public class FaultScanner3 {

  /**
   * Method used to smooth semblance numerators and denominators.
   */
  public enum Smoother {
    ROTATE_AND_SHEAR,
    FFT_GAUSSIAN
  }

  /**
   * Constructs a fault scanner with specified parameters.
   * @param sigmaPhi half-width for smoothing along strike of fault planes.
   * @param sigmaTheta half-width for smoothing up-down dip of fault planes.
   * @param snd array {snum,sden} of semblance numerators/denominators.
   */
  public FaultScanner3(
    double sigmaPhi, double sigmaTheta, float[][][][] snd)
  {
    this(sigmaPhi,sigmaTheta,snd,Smoother.ROTATE_AND_SHEAR);
  }

  /**
   * Constructs a fault scanner with specified parameters.
   * @param sigmaPhi half-width for smoothing along strike of fault planes.
   * @param sigmaTheta half-width for smoothing up-down dip of fault planes.
   * @param snd array {snum,sden} of semblance numerators/denominators.
   * @param s smoothing method.
   */
  public FaultScanner3(
    double sigmaPhi, double sigmaTheta, float[][][][] snd, Smoother s)
  {
    _sigmaPhi = sigmaPhi;
    _sigmaTheta = sigmaTheta;
    if (s==Smoother.ROTATE_AND_SHEAR) {
      _snd = snd;
    } else {
      _fps = new FaultPlaneSmoother(_sigmaTheta,_sigmaPhi,0.0,snd);
    }
    _n1 = snd[0][0][0].length;
    _n2 = snd[0][0].length;
    _n3 = snd[0].length;
  }

  /**
   * Computes and returns fault likelihood for only one strike and dip.
   * @param phi fault strike, in degrees.
   * @param theta fault dip, in degrees.
   * @return array of fault likelihoods.
   */
  public float[][][] likelihood(double phi, double theta) {
    float[][][][] fpt = scan(phi,phi,theta,theta);
    return fpt[0];
  }

  /**
   * Scans within specified bounds for fault strikes and dips.
   * @param phiMin minimum fault strike, in degrees.
   * @param phiMax maximum fault strike, in degrees.
   * @param thetaMin minimum fault dip, in degrees.
   * @param thetaMax maximum fault dip, in degrees.
   * @return array {f,p,t} of fault likelihoods, strikes, and dips.
   */
  public float[][][][] scan(
    double phiMin, double phiMax,
    double thetaMin, double thetaMax)
  {
    Sampling sp = makePhiSampling(phiMin,phiMax);
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scan(sp,st);
  }

  /**
   * Scans with the specified sampling of fault strikes and dips.
   * @param phiSampling sampling of fault strikes, in degrees.
   * @param thetaSampling sampling of fault dip angles, in degrees.
   * @return array {f,p,t} of fault likelihoods, strikes, and dips.
   */
  public float[][][][] scan(Sampling phiSampling, Sampling thetaSampling) {
    if (_snd!=null) {
      return scanS(phiSampling,thetaSampling);
    } else {
      return scanF(phiSampling,thetaSampling);
    }
  }

  public static float[][][][] thin(float[][][][] fpt) {
    int n1 = n1(fpt);
    int n2 = n2(fpt);
    int n3 = n3(fpt);
    float[][][] f = fpt[0];
    float[][][] p = fpt[1];
    float[][][] t = fpt[2];
    float[][][] ff = new float[n3][n2][n1];
    float[][][] pp = new float[n3][n2][n1];
    float[][][] tt = new float[n3][n2][n1];
    for (int i3=1; i3<n3-1; ++i3) {
      for (int i2=1; i2<n2-1; ++i2) {
        float[] fmm = f[i3-1][i2-1];
        float[] fm0 = f[i3-1][i2  ];
        float[] fmp = f[i3-1][i2+1];
        float[] f0m = f[i3  ][i2-1];
        float[] f00 = f[i3  ][i2  ];
        float[] f0p = f[i3  ][i2+1];
        float[] fpm = f[i3+1][i2-1];
        float[] fp0 = f[i3+1][i2  ];
        float[] fpp = f[i3+1][i2+1];
        float[] p00 = p[i3  ][i2  ];
        float[] t00 = t[i3  ][i2  ];
        for (int i1=0; i1<n1; ++i1) {
          float f000 = f00[i1];
          float p000 = p00[i1];
          float t000 = t00[i1];
          if (
            (-90.0f<=p000 && p000<=-67.5f && f0m[i1]<f000 && f0p[i1]<f000) ||
            (-67.5f<=p000 && p000<=-22.5f && fmm[i1]<f000 && fpp[i1]<f000) ||
            (-22.5f<=p000 && p000<= 22.5f && fm0[i1]<f000 && fp0[i1]<f000) ||
            ( 22.5f<=p000 && p000<= 67.5f && fmp[i1]<f000 && fpm[i1]<f000) ||
            ( 67.5f<=p000 && p000<= 90.0f && f0m[i1]<f000 && f0p[i1]<f000)) {
            ff[i3][i2][i1] = f000;
            pp[i3][i2][i1] = p000;
            tt[i3][i2][i1] = t000;
          }
        }
      }
    }
    return new float[][][][]{ff,pp,tt};
  }

  public static float[][][] smooth(
    double sigma, float[][][] p2, float[][][] p3, 
    float[][][] fl, float[][][] g)
  {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    EigenTensors3 d = new EigenTensors3(n1,n2,n3,true);
    d.setEigenvalues(0.001f,1.00f,1.00f); // smooth mostly along structure
    float[][][] s = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          //s[i3][i2][i1] = 1.0f-pow(fl[i3][i2][i1],0.1f); // almost binary
          s[i3][i2][i1] = (fl[i3][i2][i1]<0.1f)?1.0f:0.0f; // binary
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float u1i = 1.0f/sqrt(1.0f+p2i*p2i+p3i*p3i);
          float u2i = -p2i*u1i;
          float u3i = -p3i*u1i;
          float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
          float w1i = -u2i*usi;
          float w2i =  u1i*usi;
          float w3i = 0.0f;
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
        }
      }
    }
    float c = (float)(0.5*sigma*sigma);
    float[][][] h = new float[n3][n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,s,g,h);
    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][][] _snd;
  private FaultPlaneSmoother _fps;
  private double _sigmaPhi,_sigmaTheta;
  private int _n1,_n2,_n3;

  // This scan smooths semblance numerators and denominators along fault
  // planes using Gaussian filters implemented with FFTs.
  private float[][][][] scanF(Sampling phiSampling, Sampling thetaSampling) {
    // Algorithm: given snum,sden (semblance numerators and denominators)
    // construct fault plane smoother for snum,sden
    // initialize f,p,t (fault likelihood, phi, and theta)
    // for all phi:
    //   for all theta:
    //   use fault plane smoother to compute semblance
    //   compute fault likelihood from semblance
    //   update f,p,t for maximum likelihood
    Sampling sp = phiSampling;
    Sampling st = thetaSampling;
    FaultSemblance fs = new FaultSemblance();
    final int n1 = _n1;
    final int n2 = _n2;
    final int n3 = _n3;
    final float[][][] f = new float[n3][n2][n1];
    final float[][][] p = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    int np = sp.getCount();
    int nt = st.getCount();
    for (int ip=0; ip<np; ++ip) {
      System.out.println("FaultScanner3.scanF: ip="+ip);
      final float phi = (float)sp.getValue(ip);
      for (int it=0; it<nt; ++it) {
        final float theta = (float)st.getValue(it);
        float[][][][] snd = _fps.apply(phi,theta);
        final float[][][] s = fs.semblanceFromNumDen(snd);
        loop(n3,new LoopInt() {
        public void compute(int i3) {
          for (int i2=0; i2<n2; ++i2) {
            float[] s32 = s[i3][i2];
            float[] f32 = f[i3][i2];
            float[] p32 = p[i3][i2];
            float[] t32 = t[i3][i2];
            for (int i1=0; i1<n1; ++i1) {
              float si = s32[i1]; // semblance
              si = si*si; // semblance^2
              si = si*si; // semblance^4
              si = si*si; // semblance^8
              float fi = 1.0f-si;
              if (fi>f32[i1]) {
                f32[i1] = fi;
                p32[i1] = phi;
                t32[i1] = theta;
              }
            }
          }
        }});
      }
    }
    return new float[][][][]{f,p,t};
  }

  // This scan smooths semblance numerators and denominators along fault
  // planes by first rotating and shearing those factors before applying
  // fast recursive axis-aligned smoothing filters.
  private float[][][][] scanS(Sampling phiSampling, Sampling thetaSampling) {
    // Algorithm: given snum,sden (semblance numerators and denominators)
    // initialize f,p,t (fault likelihood, phi, and theta)
    // for all phi:
    //   rotate snum,sden so that strike vector is aligned with axis 2
    //   smooth snum,sden along fault strike
    //   compute fphi,tphi (fault likelihood and theta, as for 2D faults)
    //   unrotate fphi,tphi to original coordinates
    //   update f,p,t for maximum likelihood
    Sampling sp = phiSampling;
    Sampling st = thetaSampling;
    final int n1 = _n1, n2 = _n2, n3 = _n3;
    final float[][][] f = new float[n3][n2][n1];
    final float[][][] p = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    int np = sp.getCount();
    for (int ip=0; ip<np; ++ip) {
      System.out.println("FaultScanner3.scanS: ip="+ip);
      final float phi = (float)sp.getValue(ip);
      Rotator r = new Rotator(phi,n1,n2,n3);
      float[][][][] rsnd = r.rotate(_snd);
      smooth2(rsnd);
      float[][][][] rftp = scanTheta(st,rsnd); rsnd = null;
      float[][][][] ftp = r.unrotate(rftp);
      final float[][][] fp = ftp[0];
      final float[][][] tp = ftp[1];
      loop(n3,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] f32 = f[i3][i2];
          float[] p32 = p[i3][i2];
          float[] t32 = t[i3][i2];
          float[] fp32 = fp[i3][i2];
          float[] tp32 = tp[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float fpi = fp32[i1];
            float tpi = tp32[i1];
            if (fpi<0.0f) fpi = 0.0f;
            if (fpi>1.0f) fpi = 1.0f;
            if (fpi>f32[i1]) {
              f32[i1] = fpi;
              p32[i1] = phi;
              t32[i1] = tpi;
            }
          }
        }
      }});
    }
    return new float[][][][]{f,p,t};
  }

  // Sampling of angles depends on extent of smoothing.
  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }
  private Sampling makePhiSampling(double phiMin, double phiMax) {
    return angleSampling(_sigmaPhi,phiMin,phiMax);
  }
  private Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigmaTheta,thetaMin,thetaMax);
  }

  // Numbers of samples in 3D arrays (arrays of arrays of arrays)
  // that contain null arrays after horizontal rotation.
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
    return n3(f[0]);
  }

  // Get/set non-null slices of rotated 3D arrays
  private static float[][] extractSlice2(int i2, float[][][] x) {
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
  private static void restoreSlice2(int i2, float[][][] x, float[][] x2) {
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
  private static float[][] extractSlice3(int i3, float[][][] x) {
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
  private static void restoreSlice3(int i3, float[][][] x, float[][] x3) {
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

  // Horizontal smoothing of rotated snum,sden along axis 2.
  private void smooth2(final float[][][][] snd) {
    final int n1 = n1(snd), n2 = n2(snd), n3 = n3(snd);
    final RecursiveExponentialFilter ref = makeRef(_sigmaPhi);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int is=0; is<2; ++is) {
        float[][] s3 = extractSlice3(i3,snd[is]);
        if (s3==null)
          continue;
        ref.apply2(s3,s3); 
        restoreSlice3(i3,snd[is],s3);
      }
    }});
  }

  // Smoothing filter
  private static RecursiveExponentialFilter makeRef(double sigma) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    return ref;
  }

  // Scan over fault angles theta to maximize fault likelihood.
  private float[][][][] scanTheta(Sampling thetaSampling, float[][][][] snd) {
    final int n1 = n1(snd), n2 = n2(snd), n3 = n3(snd);
    final Sampling st = thetaSampling;
    final float[][][][] ft = like(snd);
    final float[][][] sn = snd[0];
    final float[][][] sd = snd[1];
    final float[][][] f = ft[0];
    final float[][][] t = ft[1];
    final Unsafe<SincInterpolator> usi = new Unsafe<SincInterpolator>();
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      SincInterpolator si = usi.get();
      if (si==null) {
        usi.set(si=new SincInterpolator());
        si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
      }
      float[][] sn2 = extractSlice2(i2,sn);
      float[][] sd2 = extractSlice2(i2,sd);
      if (sn2==null)
        return;
      int n3 = sn2.length;
      int nt = st.getCount();
      for (int it=0; it<nt; ++it) {
        float ti = (float)st.getValue(it);
        float theta = toRadians(ti);
        float shear = tan(theta);
        float sigma = (float)_sigmaTheta*cos(theta);
        float[][] sns = shear(si,shear,sn2);
        float[][] sds = shear(si,shear,sd2);
        RecursiveExponentialFilter ref = makeRef(_sigmaTheta);
        ref.apply1(sns,sns);
        ref.apply1(sds,sds);
        float[][] ss = semblance(sns,sds);
        float[][] s2 = unshear(si,shear,ss);
        for (int i3=0,j3=i3lo(i2,f); i3<n3; ++i3,++j3) {
          float[] s32 = s2[i3];
          float[] f32 = f[j3][i2];
          float[] t32 = t[j3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float st = s32[i1]; // semblance
            st = st*st; // semblance^2
            st = st*st; // semblance^4
            st = st*st; // semblance^8
            float fi = 1.0f-st;
            if (fi>f32[i1]) {
              f32[i1] = fi;
              t32[i1] = ti;
            }
          }
        }
      }
    }});
    return ft;
  }

  // Makes an array like that specified, which may contain null arrays.
  private float[][][][] like(float[][][][] p) {
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

  // Computes semblance from specified numerators and denominators.
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
    private static int HALF_LSINC; // half length of sinc interpolator
    static {
      SincInterpolator si = new SincInterpolator();
      _siTable = si.getTable();
      HALF_LSINC = _siTable[0].length/2;
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
      return _s2p.getFirst()-HALF_LSINC<=x2p && 
             _s3p.getFirst()-HALF_LSINC<=x3p && 
              x2p<=_s2p.getLast()+HALF_LSINC &&
              x3p<=_s3p.getLast()+HALF_LSINC;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // fault analysis

  /*
  public Faults findFaults(float[][][][] fpt, int minLength) {
    return new Faults(fpt,minLength);
  }

  public static class Faults {
    Faults(float[][][] ft, int minLength) {
      _n1 = ft[0][0].length;
      _n2 = ft[0].length;
      FaultNode[][] fns = findFaultNodes(ft);
      _fls = findFaultLines(fns,minLength);
    }

    public void findShifts(float[][] g, float[][] p, int smin, int smax) {
      for (FaultLine fl:_fls)
        findShiftsFl(fl,g,p,smin,smax);
    }

    public void clean() {
      ArrayList<FaultLine> fls = new ArrayList<FaultLine>();
      for (FaultLine fl:_fls)
        if (shiftsValid(fl))
          fls.add(fl);
      _fls = fls;
    }

    public int getCount() {
      return _fls.size();
    }

    public float[][] getShifts() {
      float[][] s = new float[_n2][_n1];
      for (FaultLine fl:_fls) {
        FaultNode[] fns = fl.getNodes();
        int nn = fns.length;
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns[in];
          int i1 = fn.i1;
          int i2m = fn.i2m;
          int i2p = fn.i2p;
          s[i2m][i1] = fn.smp;
          s[i2p][i1] = fn.spm;
        }
      }
      return s;
    }

    public float[][] getLikelihoods() {
      float[][] f = new float[_n2][_n1];
      for (FaultLine fl:_fls) {
        FaultNode[] fns = fl.getNodes();
        int nn = fns.length;
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns[in];
          int i1 = fn.i1;
          int i2 = fn.i2;
          f[i2][i1] = fn.fl;
        }
      }
      return f;
    }

    public float[][] getDips() {
      float[][] t = new float[_n2][_n1];
      for (FaultLine fl:_fls) {
        FaultNode[] fns = fl.getNodes();
        int nn = fns.length;
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns[in];
          int i1 = fn.i1;
          int i2 = fn.i2;
          t[i2][i1] = fn.ft;
        }
      }
      return t;
    }

    private int _n1,_n2;
    private ArrayList<FaultLine> _fls;
  }
  */

  // One node in a fault corresponds to one image sample.
  private static class FaultNode {
    FaultNode(
      int i1, int i2, int i3, 
      int i2m, int i2p, int i3m, int i3p, float x2, float x3, 
      float fl, float fp, float ft) 
    {
      this.i1 = i1; this.i2 = i2; this.i3 = i3;
      this.i2m = i2m; this.i2p = i2p;
      this.i3m = i3m; this.i3p = i3p;
      this.x2 = x2; this.x3 = x3;
      this.fl = fl; this.fp = fp; this.ft = ft;
    }
    int i1,i2,i3; // sample indices (i1,i2,i3) for this node
    int i2m,i2p,i3m,i3p; // indices i2,i3 for minus and plus sides of fault
    float x2,x3; // horizontal fault location in [i2m,i2p],[i3m,i3p]
    float fl,fp,ft; // fault likelihood, phi and theta
    float smp,spm; // fault shifts for both mp and pm sides of fault
    FaultNode above,below,left,right; // node neighbors
  }

  //private static class FaultSurf {
  //  ArrayList<FaultLineH> flh;
  //  ArrayList<FaultLineV> flv;
  //}

  // A vertical fault line is a linked top-to-bottom list of fault nodes.
  private static class FaultLineV {
    FaultLineV(FaultNode top) {
      this.top = top;
      this.n = 1;
      for (FaultNode fn=top.below; fn!=null; fn=fn.below)
        ++n;
    }
    int n; // number of nodes
    FaultNode top; // top node
    FaultNode[] getNodes() {
      FaultNode[] fns = new FaultNode[n];
      int in = 0;
      for (FaultNode fn=top; fn!=null; fn=fn.below,++in)
        fns[in] = fn;
      return fns;
    }
  }

  // Returns index for specified fault strike:
  //   0 for phi ~ -90 (or phi ~ 90) (fault strike aligned with axis 3)
  //   1 for phi ~ -45
  //   2 for phi ~   0 (fault strike aligned with axis 2)
  //   3 for phi ~  45
  private static int indexPhi(float phi) {
    return (int)((phi+112.5f)/45.0f)%4;
  }

  private static FaultNode[][][] findFaultNodes(float[][][][] fpt) {
    int n1 = n1(fpt);
    int n2 = n2(fpt);
    int n3 = n3(fpt);
    float[][][] f = fpt[0];
    float[][][] p = fpt[1];
    float[][][] t = fpt[2];
    float[][][] fs = new float[n3][n2][n1];
    new RecursiveGaussianFilter(1.0).apply000(f,fs);
    FaultNode[][][] fns = new FaultNode[n3][n2][n1];
    for (int i3=1; i3<n3-1; ++i3) {
      for (int i2=1; i2<n2-1; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          int ifp = indexPhi(p[i3][i2][i1]);
          if (ifp==0) {
            fns[i3][i2][i1] = makeFaultNode0(i1,i2,i3,fs,f,p,t);
          } else if (ifp==1) {
            fns[i3][i2][i1] = makeFaultNode1(i1,i2,i3,fs,f,p,t);
          } else if (ifp==2) {
            fns[i3][i2][i1] = makeFaultNode2(i1,i2,i3,fs,f,p,t);
          } else if (ifp==3) {
            fns[i3][i2][i1] = makeFaultNode3(i1,i2,i3,fs,f,p,t);
          }
        }
      }
    }
    return fns;
  }
  private static FaultNode makeFaultNode0(
    int i1, int i2, int i3, 
    float[][][] fs, float[][][] f, float[][][] p, float[][][] t) 
  {
    float f0m = fs[i3][i2-1][i1]; //   o
    float f00 = fs[i3][i2  ][i1]; // - o -  (phi ~ 90 or phi ~ 90)
    float f0p = fs[i3][i2+1][i1]; //   o
    FaultNode fn = null;
    if (f0m<f00 && f0p<f00) {
      float u = peak(f0m,f00,f0p);
      float x2 = i2+u;
      float x3 = i3;
      int i2m = i2, i2p = i2;
      int i3m = i3, i3p = i3;
      if (u>=0.0f) {
        ++i2p;
      } else {
        --i2p; 
      }
      fn = new FaultNode(i1,i2,i3,i2m,i2p,i3m,i3p,x2,x3,
        f[i3][i2][i1],p[i3][i2][i1],t[i3][i2][i1]);
    }
    return fn;
  }
  private static FaultNode makeFaultNode1(
    int i1, int i2, int i3, 
    float[][][] fs, float[][][] f, float[][][] p, float[][][] t) 
  {
    float fmm = fs[i3-1][i2-1][i1]; // o  / 
    float f00 = fs[i3  ][i2  ][i1]; //   o    (phi ~ -45)
    float fpp = fs[i3+1][i2+1][i1]; // /   o
    FaultNode fn = null;
    if (fmm<f00 && fpp<f00) {
      float u = peak(fmm,f00,fpp);
      float x2 = i2+u;
      float x3 = i3+u;
      int i2m = i2, i2p = i2;
      int i3m = i3, i3p = i3;
      if (u>=0.0f) {
        ++i2p; ++i3p; 
      } else {
        --i2p; --i3m;
      }
      fn = new FaultNode(i1,i2,i3,i2m,i2p,i3m,i3p,x2,x3,
        f[i3][i2][i1],p[i3][i2][i1],t[i3][i2][i1]);
    }
    return fn;
  }
  private static FaultNode makeFaultNode2(
    int i1, int i2, int i3, 
    float[][][] fs, float[][][] f, float[][][] p, float[][][] t) 
  {
    float fm0 = fs[i3-1][i2][i1]; //    |
    float f00 = fs[i3  ][i2][i1]; // o  o  o (phi ~ 0) 
    float fp0 = fs[i3+1][i2][i1]; //    |
    FaultNode fn = null;
    if (fm0<f00 && fp0<f00) {
      float u = peak(fm0,f00,fp0);
      float x2 = i2;
      float x3 = i3+u;
      int i2m = i2, i2p = i2;
      int i3m = i3, i3p = i3;
      if (u>=0.0f) {
        ++i3p;
      } else {
        --i3m;
      }
      fn = new FaultNode(i1,i2,i3,i2m,i2p,i3m,i3p,x2,x3,
        f[i3][i2][i1],p[i3][i2][i1],t[i3][i2][i1]);
    }
    return fn;
  }
  private static FaultNode makeFaultNode3(
    int i1, int i2, int i3, 
    float[][][] fs, float[][][] f, float[][][] p, float[][][] t) 
  {
    float fpm = fs[i3+1][i2-1][i1]; //  \  o
    float f00 = fs[i3  ][i2  ][i1]; //   o    (phi ~ 45)
    float fmp = fs[i3-1][i2+1][i1]; // o  \
    FaultNode fn = null;
    if (fpm<f00 && fmp<f00) {
      float u = peak(fpm,f00,fmp);
      float x2 = i2-u;
      float x3 = i3+u;
      int i2m = i2, i2p = i2;
      int i3m = i3, i3p = i3;
      if (u>=0.0f) {
        --i2p; ++i3p; 
      } else {
        ++i2p; --i3m;
      }
      fn = new FaultNode(i1,i2,i3,i2m,i2p,i3m,i3p,x2,x3,
        f[i3][i2][i1],p[i3][i2][i1],t[i3][i2][i1]);
    }
    return fn;
  }
  private static float peak(float fa, float fb, float fc) {
    float c0 = fb;
    float c1 = 0.5f*(fc-fa);
    float c2 = 0.5f*(fc+fa)-fb;
    float u = (c2<0.0f)?-0.5f*c1/c2:0.0f;
    return u;
  }

  /*
  private static ArrayList<FaultLine> findFaultLines(
    FaultNode[][] fns, int minNodesInLine) 
  {
    int n1 = fns[0].length;
    int n2 = fns.length;

    // List of fault lines.
    ArrayList<FaultLine> fls = new ArrayList<FaultLine>();

    // For all fault nodes in the array of nodes, ...
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        FaultNode fni = fns[i2][i1];

        // If a node exists and has not yet been linked to others, ...
        if (fni!=null && fni.above==null && fni.below==null) {

          // Count the number of nodes.
          int nn = 1;

          // Search below for nabors to link to this node.
          for (FaultNode fnt=fni; fnt!=null; fnt=fnt.below) {
            FaultNode fnb = naborBelow(fnt,fns);
            if (fnb!=null && fnt==naborAbove(fnb,fns)) {
              fnb.above = fnt;
              fnt.below = fnb;
              ++nn;
            }
          }

          // Search above for nabors to link to this node.
          // Also remember the top node, which may become the
          // head node in the list of nodes for a fault line.
          FaultNode fnr = fni;
          for (FaultNode fnt=fni; fnt!=null; fnt=fnt.above) {
            FaultNode fna = naborAbove(fnt,fns);
            if (fna!=null && fnt==naborBelow(fna,fns)) {
              fna.below = fnt;
              fnt.above = fna;
              ++nn;
            }
            fnr = fnt;
          }

          // If sufficient number of nodes, construct new fault line.
          if (nn>=minNodesInLine)
            fls.add(new FaultLine(fnr));
        }
      }
    }
    return fls;
  }

  private static float[][] getImageSamples(FaultLine fl, int k, float[][] g) {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    int n1 = g[0].length;
    int n2 = g.length;
    int n2m = n2-1;
    float[] gm = new float[nn];
    float[] gp = new float[nn];
    for (int in=0; in<nn; ++in) {
      FaultNode fn = fns[in];
      int i1 = fn.i1;
      int i2 = fn.i2;
      int i2m = max(  0,i2-k);
      int i2p = min(n2m,i2+k);
      gm[in] = g[i2m][i1];
      gp[in] = g[i2p][i1];
    }
    return new float[][]{gm,gp};
  }

  private static boolean shiftsValid(FaultLine fl) {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    int nv = 0;
    for (int in=0; in<nn; ++in) {
      FaultNode fn = fns[in];
      if (fn.smp*fn.spm<=0.0f)
        ++nv;
    }
    return 10*nv>8*nn; // 80% of shifts must have opposite signs
  }

  private static void cleanShifts(FaultLine fl) {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    for (int in=0; in<nn; ++in) {
      FaultNode fn = fns[in];
      if (fn.smp*fn.spm>=0.0f) {
        fn.smp = 0.0f;
        fn.spm = 0.0f;
      }
    }
  }

  private static void findShiftsFl(
    FaultLine fl, float[][] g, float[][] p, int smin, int smax) 
  {
    findShiftsFl( 1,fl,g,p,smin,smax);
    findShiftsFl(-1,fl,g,p,smin,smax);
  }

  private static void findShiftsFl(
    int sign, FaultLine fl, float[][] g, float[][] p, int smin, int smax) 
  {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    int ssmin = (nn+smin>10)?smin:10-nn;
    int ssmax = (nn-smax>10)?smax:nn-10;
    DynamicWarping dw = new DynamicWarping(ssmin,ssmax);
    dw.setStretchMax(0.2);
    float[][] gmp = getImageSamples(fl,2,g);
    float[][] pmp = getImageSamples(fl,2,p);
    float[] g1,g2,p1,p2;
    float ps;
    if (sign>0) {
      g1 = gmp[0]; g2 = gmp[1];
      p1 = pmp[0]; p2 = pmp[1];
      ps = 2.0f;
    } else {
      g1 = gmp[1]; g2 = gmp[0];
      p1 = pmp[1]; p2 = pmp[0];
      ps = -2.0f;
    }
    float[][] se = dw.computeErrors(g1,g2);
    float[][] sd = dw.accumulateForward(se);
    float[] s = dw.findShiftsReverse(sd,se);
    for (int in=0; in<nn; ++in)
      s[in] -= ps*(p1[in]+p2[in]);
    if (sign>0) {
      for (int in=0; in<nn; ++in)
        fns[in].smp = s[in];
    } else {
      for (int in=0; in<nn; ++in)
        fns[in].spm = s[in];
    }
  }
  */

  // Minimum acceptable dot product of two aligned fault dip vectors.
  // Nodes can be neighbors only if their dip vectors are so aligned.
  private static final float AMIN = 0.9f;

  // Returns dot product of dip vectors for two specified nodes.
  private static float dotDipVectors(FaultNode fn1, FaultNode fn2) {
    float t1 = toRadians(fn1.ft);
    float t2 = toRadians(fn2.ft);
    float c1 = cos(t1);
    float c2 = cos(t2);
    float s1 = sin(t1);
    float s2 = sin(t2);
    return c1*c2+s1*s2;
  }

  // Determines which (if any) of nodes f1,f2,f3 is aligned with fn.
  private static FaultNode alignDip(
    FaultNode fn, FaultNode f1, FaultNode f2, FaultNode f3) 
  {
    float a1 = (f1!=null)?dotDipVectors(fn,f1):0.0f;
    float a2 = (f2!=null)?dotDipVectors(fn,f2):0.0f;
    float a3 = (f3!=null)?dotDipVectors(fn,f3):0.0f;
    float am = max(a1,a2,a3);
    FaultNode fm = null;
    if (am>=AMIN) 
      fm = (a1==am)?f1:(a2==am)?f2:f3;
    return fm;
  }

  // Returns the best nabor (if any) above the specified node.
  private static FaultNode naborAbove(FaultNode fn, FaultNode[][] fns) {
    int n1 = fns[0].length;
    int n2 = fns.length;
    int i1 = fn.i1;
    int i2 = fn.i2;
    FaultNode fm = null;
    if (i1>0) {
      FaultNode f1 = fns[i2][i1-1];
      FaultNode f2 = (i2>0)?fns[i2-1][i1-1]:null;
      FaultNode f3 = (i2<n2-1)?fns[i2+1][i1-1]:null;
      fm = alignDip(fn,f1,f2,f3);
    }
    return fm;
  }

  // Returns the best nabor (if any) below the specified node.
  private static FaultNode naborBelow(FaultNode fn, FaultNode[][] fns) {
    int n1 = fns[0].length;
    int n2 = fns.length;
    int i1 = fn.i1;
    int i2 = fn.i2;
    FaultNode fm = null;
    if (i1<n1-1) {
      FaultNode f1 = fns[i2][i1+1];
      FaultNode f2 = (i2>0)?fns[i2-1][i1+1]:null;
      FaultNode f3 = (i2<n2-1)?fns[i2+1][i1+1]:null;
      fm = alignDip(fn,f1,f2,f3);
    }
    return fm;
  }
}
