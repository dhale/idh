/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import java.util.*;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

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
      trace("FaultScanner3.scanF: ip/np="+ip+"/"+np);
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
      trace("FaultScanner3.scanS: ip/np="+ip+"/"+np);
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
      //trace("s2p: n2p="+n2);
      //trace("s3p: n3p="+n3);
      //trace("s2q: n2q="+n2q+" d2q="+d2q+" f2q="+f2q);
      //trace("s3q: n3q="+n3q+" d3q="+d3q+" f3q="+f3q);
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

  public static Faults findFaults(float[][][][] fpt, int minSize) {
    return new Faults(fpt,minSize);
  }

  public static class Faults {
    Faults(float[][][][] fpt, int minSize) {
      _n1 = fpt[0][0][0].length;
      _n2 = fpt[0][0].length;
      _n3 = fpt[0].length;
      trace("find nodes ...");
      FaultNode[][][] fns = findFaultNodes(fpt);
      trace("find surfs ...");
      _fsl = findFaultSurfs(fns,minSize);
    }

    public void findShifts(
      float[][][] g, float[][][][] p, int smin, int smax) 
    {
      //for (FaultSurf fs:_fsl)
      //  findShiftsFs(fs,g,p,smin,smax);
    }

    public int getCount() {
      return _fsl.size();
    }

    public float[][][] getShifts() {
      float[][][] s = new float[_n3][_n2][_n1];
      for (FaultSurf fs:_fsl) {
        ArrayList<FaultNode> fns = fs.getNodes();
        int nn = fns.size();
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns.get(in);
          int i1 = fn.i1;
          int i2m = fn.i2m;
          int i2p = fn.i2p;
          int i3m = fn.i3m;
          int i3p = fn.i3p;
          s[i3m][i2m][i1] = -1.0f;
          s[i3p][i2p][i1] =  1.0f;
        }
      }
      return s;
    }

    public float[][][] getLikelihoods() {
      float[][][] f = new float[_n3][_n2][_n1];
      for (FaultSurf fs:_fsl) {
        ArrayList<FaultNode> fns = fs.getNodes();
        int nn = fns.size();
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns.get(in);
          int i1 = fn.i1;
          int i2 = fn.i2;
          int i3 = fn.i3;
          f[i3][i2][i1] = fn.fl;
        }
      }
      return f;
    }

    private int _n1,_n2,_n3;
    private ArrayList<FaultSurf> _fsl;
  }

  // One node in a fault corresponds to one image sample.
  private static class FaultNode {
    int i1,i2,i3; // sample indices (i1,i2,i3) for this node
    int i2m,i2p,i3m,i3p; // indices i2,i3 for minus and plus sides of fault
    float x2,x3; // horizontal fault location in [i2m,i2p],[i3m,i3p]
    float fl,fp,ft; // fault likelihood, phi and theta
    float smp,spm; // fault shifts for both mp and pm sides of fault
    FaultNode above,below,left,right; // node neighbors
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
    void flipOrientation() {
      int i2t = i2m; i2m = i2p; i2p = i2t;
      int i3t = i3m; i3m = i3p; i3p = i3t;
      float st = smp; smp = spm; spm = st;
      FaultNode fnt = left; left = right; right = fnt;
    }
    public String toString() {
      return "("+i1+","+i2+","+i3+")";
    }
  }

  private static class FaultSurf {
    public FaultSurf(FaultNode node, FaultNode[][][] fns) {
      //System.out.print("linking ... ");
      linkNodeNabors(node,fns);
      //System.out.print("cleaning ... ");
      //cleanNodeNabors();
      checkNodeNabors();
      //System.out.println("count="+countNodes());
    }
    public int countNodes() {
      return _nodes.size();
    }
    public ArrayList<FaultNode> getNodes() {
      int n = countNodes();
      ArrayList<FaultNode> fnl = new ArrayList<FaultNode>(n);
      for (FaultNode fn:_nodes)
        fnl.add(fn);
      return fnl;
    }
    HashSet<FaultNode> _nodes = new HashSet<FaultNode>();

    // Recursively link nodes with nabors above, below, left and right.
    private void linkNodeNabors(FaultNode fn, FaultNode[][][] fns) {
      ArrayDeque<FaultNode> fnq = new ArrayDeque<FaultNode>();
      fnq.add(fn);
      while (!fnq.isEmpty()) {
        fn = fnq.removeLast();
        FaultNode fna = null;
        FaultNode fnb = null;
        FaultNode fnl = null;
        FaultNode fnr = null;
        if (fn.above==null) {
          fna = naborAbove(fn,fns);
          if (fna!=null && fn==naborBelow(fna,fns)) {
            fn.above = fna;
            fna.below = fn;
          } else {
            fna = null;
          }
        }
        if (fn.below==null) {
          fnb = naborBelow(fn,fns);
          if (fnb!=null && fn==naborAbove(fnb,fns)) {
            fn.below = fnb;
            fnb.above = fn;
          } else {
            fnb = null;
          }
        }
        if (fn.left==null) {
          fnl = naborLeft(fn,fns);
          if (fnl!=null && fn==naborLeft(fnl,fns))
            fnl.flipOrientation();
          if (fnl!=null && fn==naborRight(fnl,fns)) {
            fn.left = fnl;
            fnl.right = fn;
          } else {
            fnl = null;
          }
        }
        if (fn.right==null) {
          fnr = naborRight(fn,fns);
          if (fnr!=null && fn==naborRight(fnr,fns))
            fnr.flipOrientation();
          if (fnr!=null && fn==naborLeft(fnr,fns)) {
            fn.right = fnr;
            fnr.left = fn;
          } else {
            fnr = null;
          }
        }

        if ((fn.above!=null || fn.below!=null) &&
            (fn.left !=null || fn.right!=null)) {
          _nodes.add(fn);
          /*
          if (fn.below!=null && fn.right!=null) {
            if (linkNaborsLeftRight(fn.below,naborBelow(fn.right,fns),fns))
              linkNaborsAboveBelow(fn.right,fn.below.right,fns);
          }
          if (fn.right!=null && fn.above!=null) {
            if (linkNaborsLeftRight(fn.above,naborAbove(fn.right,fns),fns))
              linkNaborsAboveBelow(fn.above.right,fn.right,fns);
          }
          if (fn.above!=null && fn.left!=null) {
            if (linkNaborsLeftRight(naborAbove(fn.left,fns),fn.above,fns))
              linkNaborsAboveBelow(fn.above.left,fn.left,fns);
          }
          if (fn.left!=null && fn.below!=null) {
            if (linkNaborsLeftRight(naborBelow(fn.left,fns),fn.below,fns))
              linkNaborsAboveBelow(fn.left,fn.below.left,fns);
          }
          */
        }
        if (fna!=null) fnq.add(fna);
        if (fnb!=null) fnq.add(fnb);
        if (fnl!=null) fnq.add(fnl);
        if (fnr!=null) fnq.add(fnr);
      }
    }

    // Attempts to link above-below nodes. Returns true, if the
    // nodes are successfully linked (or were already linked) as 
    // above-below nabors; false, otherwise.
    private static boolean linkNaborsAboveBelow(
      FaultNode fna, FaultNode fnb, FaultNode[][][] fns) 
    {
      if (fna==null || fnb==null)
        return false;
      if (fna.below==fnb.above)
        return true;
      FaultNode fnab = naborBelow(fna,fns);
      FaultNode fnba = naborAbove(fnb,fns);
      if (fna==fnba && fnb==fnab) {
        fna.below = fnb;
        fnb.above = fna;
        return true;
      }
      return false;
    }
    
    // Attempts to link specified left-right nodes. If necessary,
    // modifies the orientations of those nodes to make them nabors. 
    // Returns true, if the nodes are successfully linked (or were 
    // already linked) as left-right nabors; false, otherwise.
    private static boolean linkNaborsLeftRight(
      FaultNode fnl, FaultNode fnr, FaultNode[][][] fns) 
    {
      if (fnl==null || fnr==null)
        return false;
      if (fnl.right==fnr.left)
        return true;
      FaultNode fnll = naborLeft(fnl,fns);
      FaultNode fnlr = naborRight(fnl,fns);
      FaultNode fnrl = naborLeft(fnr,fns);
      FaultNode fnrr = naborRight(fnr,fns);
      if ((fnl==fnrl || fnl==fnrr) && (fnr==fnlr || fnr==fnll)) {
        if (fnl==fnrr) {
          fnr.flipOrientation();
          fnrl = fnrr;
        }
        if (fnr==fnll) {
          fnl.flipOrientation();
          fnlr = fnll;
        }
        fnl.right = fnr;
        fnr.left = fnl;
        return true;
      }
      return false;
    }

    private void checkNodeNabors() {
      for (FaultNode fn:_nodes) {
        if (fn.right!=null && fn.right.left!=fn)
          trace("fn.right.left!=fn for fn="+fn);
        if (fn.left!=null && fn.left.right!=fn)
          trace("fn.left.right!=fn for fn="+fn);
        /*
        if (fn.right!=null && fn.above!=null) {
          FaultNode fnra = fn.right.above;
          FaultNode fnal = fn.above.left;
          if (fnra!=null && fnal!=null && fnra==fnal) {
            trace("fnra="+fnra+" fnal="+fnal);
            trace("fnra: phi="+fnra.fp+" theta="+fnra.ft);
            trace("fnal: phi="+fnal.fp+" theta="+fnal.ft);
            //assert fnra==null || fnal==null || fnra!=fnal;
          }
        } else if (fn.right!=null && fn.below!=null) {
          FaultNode fnrb = fn.right.below;
          FaultNode fnbl = fn.below.left;
          //assert fnrb==null || fnbl==null || fnrb!=fnbl;
        } else if (fn.left!=null && fn.below!=null) {
          FaultNode fnlb = fn.left.below;
          FaultNode fnbr = fn.below.right;
          //assert fnlb==null || fnbr==null || fnlb!=fnbr;
        } else if (fn.left!=null && fn.above!=null) {
          FaultNode fnla = fn.left.above;
          FaultNode fnar = fn.above.right;
          //assert fnla==null || fnar==null || fnla!=fnar;
        }
        */
      }
    }

    // Ensure nodes belong to at least one quad loop of nodes.
    // Unlinks any nodes for which this condition is not satisfied. 
    private void cleanNodeNabors() {
      HashMap<FaultNode,Integer> fnm = new HashMap<FaultNode,Integer>();
      for (FaultNode fni:_nodes) {
        FaultNode fnr = fni.right;
        if (fnr!=null) {
          FaultNode fna = fnr.above;
          if (fna!=null) {
            FaultNode fnl = fna.left;
            if (fnl!=null) {
              FaultNode fnb = fnl.below;
              if (fni==fnb) {
                fnm.put(fnb,(fnm.containsKey(fnb)?fnm.get(fnb):0)+1);
                fnm.put(fnr,(fnm.containsKey(fnr)?fnm.get(fnr):0)+1);
                fnm.put(fna,(fnm.containsKey(fna)?fnm.get(fna):0)+1);
                fnm.put(fnl,(fnm.containsKey(fnl)?fnm.get(fnl):0)+1);
              }
            }
          }
        }
      }
      int nclean = 0;
      for (FaultNode fni:fnm.keySet()) {
        int count = fnm.get(fni);
        if (count==0) {
          _nodes.remove(fni);
          ++nclean;
        }
      }
      if (nclean>0)
        trace("cnn: nclean="+nclean);
    }
  }

  private static ArrayList<FaultSurf> findFaultSurfs(
    FaultNode[][][] fns, int minSize) 
  {
    int n1 = fns[0][0].length;
    int n2 = fns[0].length;
    int n3 = fns.length;
    ArrayList<FaultSurf> fsl = new ArrayList<FaultSurf>();
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          //if (i1!=209 || i2!=99 || i3!=29) continue;
          FaultNode fn = fns[i3][i2][i1];
          if (fn!=null &&
              fn.above==null && fn.below==null &&
               fn.left==null && fn.right==null) {
            FaultSurf fs = new FaultSurf(fn,fns);
            if (fs.countNodes()>minSize)
              fsl.add(fs);
          }
        }
      }
    }
    return fsl;
  }

  // Returns index for specified fault strike phi:
  //   0 for phi ~ -90 (strike aligned with axis 3)
  //   1 for phi ~ -45
  //   2 for phi ~   0 (strike aligned with axis 2)
  //   3 for phi ~  45
  //   0 for phi ~  90 (strike aligned with axis 3)
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
    dumpSlice(9,5,80,120,99,fs);
    dumpSlice(9,5,80,120,99,p);
    for (int i3=1; i3<n3-1; ++i3) {
    //for (int i3=131; i3<165; ++i3) {
    //for (int i3=90; i3<110; ++i3) {
      for (int i2=1; i2<n2-1; ++i2) {
      //for (int i2=28; i2<59; ++i2) {
      //for (int i2=115; i2<135; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
        //for (int i1=80; i1<85; ++i1) {
          if (f[i3][i2][i1]>FMIN) {
            float fpi = p[i3][i2][i1];
            int ifp = indexPhi(fpi);
            FaultNode fn = null;
            if (ifp==0) {
              fn = makeFaultNode0(i1,i2,i3,fs,f,p,t);
              //if (fn==null) fn = (fpi<=0.0f) ?
              //  makeFaultNode1(i1,i2,i3,fs,f,p,t) :
              //  makeFaultNode3(i1,i2,i3,fs,f,p,t);
            } else if (ifp==1) {
              fn = makeFaultNode1(i1,i2,i3,fs,f,p,t);
              //if (fn==null) fn = (fpi<=-45.0f) ?
              //  makeFaultNode0(i1,i2,i3,fs,f,p,t) :
              //  makeFaultNode2(i1,i2,i3,fs,f,p,t);
            } else if (ifp==2) {
              fn = makeFaultNode2(i1,i2,i3,fs,f,p,t);
              //if (fn==null) fn = (fpi<=0.0f) ?
              //  makeFaultNode1(i1,i2,i3,fs,f,p,t) :
              //  makeFaultNode3(i1,i2,i3,fs,f,p,t);
            } else if (ifp==3) {
              fn = makeFaultNode3(i1,i2,i3,fs,f,p,t);
              //if (fn==null) fn = (fpi<=45.0f) ?
              //  makeFaultNode2(i1,i2,i3,fs,f,p,t) :
              //  makeFaultNode0(i1,i2,i3,fs,f,p,t);
            }
            fns[i3][i2][i1] = fn;
          }
        }
      }
    }
    return fns;
  }
  private static void dumpSlice(
    int n2, int n3, int i1, int i2, int i3, float[][][] x)
  {
    float[][] s = new float[n3][n2];
    for (int j3=i3; j3<i3+n3; ++j3) {
      for (int j2=i2; j2<i2+n2; ++j2) {
        s[j3-i3][j2-i2] = x[j3][j2][i1];
      }
    }
    dump(transpose(s));
  }
  private static FaultNode makeFaultNode0(
    int i1, int i2, int i3, 
    float[][][] fs, float[][][] f, float[][][] p, float[][][] t) 
  {
    float f0m = fs[i3][i2-1][i1]; //   o-
    float f00 = fs[i3][i2  ][i1]; // - o -  (phi ~ +-90)
    float f0p = fs[i3][i2+1][i1]; //   o+
    FaultNode fn = null;
    if (i1==80 && i2==121 && i3==100)
      trace("f0m="+f0m+" f00="+f00+" f0p="+f0p);
    if (f0m<f00 && f0p<f00) {
      float u = peak(f0m,f00,f0p);
      float x2 = i2+u;
      float x3 = i3;
      int i2m = i2, i2p = i2;
      int i3m = i3, i3p = i3;
      if (u>=0.0f) {
        ++i2p;
      } else {
        --i2m; 
      }
      fn = new FaultNode(i1,i2,i3,i2m,i2p,i3m,i3p,x2,x3,
        f[i3][i2][i1],p[i3][i2][i1],t[i3][i2][i1]);
      if (i1==80 && 120<i2 && i2<130 && 95<i3 && i3<105)
        trace("fn="+fn+" i2m="+i2m+" i2p="+i2p+" i3m="+i3m+" i3p="+i3p+
              " p="+p[i3][i2][i1]);
    }
    return fn;
  }
  private static FaultNode makeFaultNode1(
    int i1, int i2, int i3, 
    float[][][] fs, float[][][] f, float[][][] p, float[][][] t) 
  {
    float fmm = fs[i3-1][i2-1][i1]; // o- / 
    float f00 = fs[i3  ][i2  ][i1]; //   o    (phi ~ -45)
    float fpp = fs[i3+1][i2+1][i1]; // /   o+
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
        --i2m; --i3m;
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
    float f00 = fs[i3  ][i2][i1]; // o- o  o+ (phi ~ 0) 
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
    float fpm = fs[i3+1][i2-1][i1]; //  \  o+
    float f00 = fs[i3  ][i2  ][i1]; //   o    (phi ~ 45)
    float fmp = fs[i3-1][i2+1][i1]; // o- \
    if (i1==80 && i2==126 && i3==102)
      trace("fpm="+fpm+" f00="+f00+" fmp="+fmp);
    FaultNode fn = null;
    if (fpm<f00 && fmp<f00) {
      float u = peak(fmp,f00,fpm);
      float x2 = i2-u;
      float x3 = i3+u;
      int i2m = i2, i2p = i2;
      int i3m = i3, i3p = i3;
      if (u>=0.0f) {
        --i2p; ++i3p;
      } else {
        ++i2m; --i3m;
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

  // Minimum acceptable fault likelihood.
  private static final float FMIN = 0.1f;

  // Minimum acceptable dot product of two aligned vectors. Nodes can 
  // be neighbors only if their dip/strike vectors are so aligned.
  private static final float AMIN = 0.9f;

  // Returns dot product of dip vectors for two specified nodes.
  private static float dotDipVectors(FaultNode fna, FaultNode fnb) {
    float pa = toRadians(fna.fp);
    float ta = toRadians(fna.ft);
    float pb = toRadians(fnb.fp);
    float tb = toRadians(fnb.ft);
    float cpa = cos(pa), cta = cos(ta);
    float spa = sin(pa), sta = sin(ta);
    float cpb = cos(pb), ctb = cos(tb);
    float spb = sin(pb), stb = sin(tb);
    float wa1 = cta, wa2 = -spa*sta, wa3 = cpa*sta;
    float wb1 = ctb, wb2 = -spb*stb, wb3 = cpb*stb;
    return wa1*wb1+wa2*wb2+wa3*wb3;
  }

  // Returns dot product of strike vectors for two specified nodes.
  private static float dotStrikeVectors(FaultNode fna, FaultNode fnb) {
    float pa = toRadians(fna.fp);
    float pb = toRadians(fnb.fp);
    float cpa = cos(pa);
    float spa = sin(pa);
    float cpb = cos(pb);
    float spb = sin(pb);
    return cpa*cpb+spa*spb;
  }

  // Determines which (if any) of nodes fa,fb,fc is dip-aligned with fn.
  private static FaultNode alignDip(
    FaultNode fn, FaultNode fa, FaultNode fb, FaultNode fc) 
  {
    float da = (fa!=null)?dotDipVectors(fn,fa):0.0f;
    float db = (fb!=null)?dotDipVectors(fn,fb):0.0f;
    float dc = (fc!=null)?dotDipVectors(fn,fc):0.0f;
    float dm = max(da,db,dc);
    FaultNode fm = null;
    if (dm>=AMIN) 
      fm = (db==dm)?fb:(da==dm)?fa:fc;
    return fm;
  }

  // Determines which (if any) of nodes fa,fb,fc is strike-aligned with fn.
  private static FaultNode alignStrike(
    FaultNode fn, FaultNode fa, FaultNode fb, FaultNode fc) 
  {
    float da = (fa!=null)?dotStrikeVectors(fn,fa):0.0f;
    float db = (fb!=null)?dotStrikeVectors(fn,fb):0.0f;
    float dc = (fc!=null)?dotStrikeVectors(fn,fc):0.0f;
    if (da<0.0f) da = -da;
    if (db<0.0f) db = -db;
    if (dc<0.0f) dc = -dc;
    float dm = max(da,db,dc);
    FaultNode fm = null;
    if (dm>=AMIN) 
      fm = (db==dm)?fb:(da==dm)?fa:fc;
    return fm;
  }

  // Returns node with specified indices, or null if out of bounds.
  private static FaultNode node(int i1, int i2, int i3, FaultNode[][][] fns) {
    int n1 = fns[0][0].length;
    int n2 = fns[0].length;
    int n3 = fns.length;
    return (0<=i1 && i1<n1 && 0<=i2 && i2<n2 && 0<=i3 && i3<n3) ?
      fns[i3][i2][i1]:
      null;
  }

  // Returns the best nabor (if any) above the specified node.
  private static FaultNode naborAbove(FaultNode fn, FaultNode[][][] fns) {
    return naborAboveBelow(-1,fn,fns);
  }
  // Returns the best nabor (if any) below the specified node.
  private static FaultNode naborBelow(FaultNode fn, FaultNode[][][] fns) {
    return naborAboveBelow(1,fn,fns);
  }
  private static FaultNode naborAboveBelow(
    int k, FaultNode fn, FaultNode[][][] fns) 
  {
    int n1 = fns[0][0].length;
    int n2 = fns[0].length;
    int n3 = fns.length;
    int i1 = fn.i1;
    int i2 = fn.i2;
    int i3 = fn.i3;
    FaultNode fa = null;
    FaultNode fb = node(i1+k,i2,i3,fns);
    FaultNode fc = null;
    int ifp = indexPhi(fn.fp);
    if (ifp==0) {
      fa = node(i1+k,i2-1,i3,fns);
      fc = node(i1+k,i2+1,i3,fns);
    } else if (ifp==1) {
      fa = node(i1+k,i2-1,i3-1,fns);
      fc = node(i1+k,i2+1,i3+1,fns);
    } else if (ifp==2) {
      fa = node(i1+k,i2,i3-1,fns);
      fc = node(i1+k,i2,i3+1,fns);
    } else if (ifp==3) {
      fa = node(i1+k,i2+1,i3-1,fns);
      fc = node(i1+k,i2-1,i3+1,fns);
    }
    return alignDip(fn,fa,fb,fc);
  }

  // Returns the best nabor (if any) left of the specified node.
  private static FaultNode naborLeft(FaultNode fn, FaultNode[][][] fns) {
    return naborLeftRight(-1,fn,fns);
  }
  private static FaultNode naborRight(FaultNode fn, FaultNode[][][] fns) {
    return naborLeftRight(1,fn,fns);
  }
  private static FaultNode naborLeftRight(
    int k, FaultNode fn, FaultNode[][][] fns) 
  {
    int n1 = fns[0][0].length;
    int n2 = fns[0].length;
    int n3 = fns.length;
    int i1 = fn.i1;
    int i2 = fn.i2;
    int i3 = fn.i3;
    int i2m = fn.i2m, i2p = fn.i2p;
    int i3m = fn.i3m, i3p = fn.i3p;
    FaultNode fa = null;
    FaultNode fb = null;
    FaultNode fc = null;
    if (i3m==i3p) {
      int k2 = (k<0)?i2p-i2m:i2m-i2p;
      assert k2==-1 || k2==1;
      fa = node(i1,i2-k2,i3-k2,fns);
      fb = node(i1,i2   ,i3-k2,fns);
      fc = node(i1,i2+k2,i3-k2,fns);
    } else if (i2m<i2p && i3m<i3p) {
      int k2 = (k<0)?i2p-i2m:i2m-i2p;
      int k3 = (k<0)?i3p-i3m:i3m-i3p;
      assert k2==-1 || k2==1;
      assert k3==-1 || k3==1;
      fa = node(i1,i2   ,i3-k3,fns);
      fb = node(i1,i2+k2,i3-k3,fns);
      fc = node(i1,i2+k2,i3   ,fns);
    } else if (i2m==i2p) {
      int k3 = (k<0)?i3p-i3m:i3m-i3p;
      assert k3==-1 || k3==1;
      fa = node(i1,i2+k3,i3-k3,fns);
      fb = node(i1,i2+k3,i3   ,fns);
      fc = node(i1,i2+k3,i3+k3,fns);
    } else {
      int k2 = (k<0)?i2p-i2m:i2m-i2p;
      int k3 = (k<0)?i3p-i3m:i3m-i3p;
      assert k2==-1 || k2==1;
      assert k3==-1 || k3==1;
      fa = node(i1,i2-k2,i3   ,fns);
      fb = node(i1,i2-k2,i3+k3,fns);
      fc = node(i1,i2   ,i3+k3,fns);
    }
    return alignStrike(fn,fa,fb,fc);
  }

  public static float[][][][] fakeFpt(int n1, int n2, int n3) {
    Random r = new Random(314159);
    float[][][] f = new float[n3][n2][n1];
    float[][][] p = new float[n3][n2][n1];
    float[][][] t = new float[n3][n2][n1];
    float c2 = n2/2+0.5f;
    float c3 = n3/2+0.5f;
    float rc = min(n2,n3)/3;
    for (int i1=0; i1<n1; ++i1) {
      float rs = 0.5f+0.5f*i1/n1;
      for (int i3=0; i3<n3; ++i3) {
        float x3 = i3-c3;
        for (int i2=0; i2<n2; ++i2) {
          float x2 = i2-c2;
          float y2 = x2+0.0f*r.nextFloat();
          float y3 = x3+0.0f*r.nextFloat();
          float rx = sqrt(y2*y2+y3*y3);
          float dr = rx-rs*rc;
          float fx = exp(-0.125f*dr*dr);
          float px = toDegrees(atan(-x2/x3));
          float tx = 0.0f;
          f[i3][i2][i1] = fx;
          p[i3][i2][i1] = px;
          t[i3][i2][i1] = tx;
        }
      }
    }
    return new float[][][][]{f,p,t};
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
