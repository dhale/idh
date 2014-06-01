/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates two wavelets from alignment by warping of sequences or images.
 * The two sequences or images are assumed to have been convolved with the
 * different wavelets. Warping of one sequence or image to align with the
 * other will cause the convolved wavelet to be stretched or squeezed, and
 * this distortion enables us to estimate both wavelets.
 * <p>
 * For images, convolution with the wavelet is assumed to be in only the 1st
 * dimension. For definiteness, this 1st dimension is assumed to be time in
 * the documentation below.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.02.15
 */
public class WaveletWarping2 {

  /**
   * Sets the min-max range of times used to estimate wavelet.
   * @param itmin minimum time, in samples.
   * @param itmax maximum time, in samples.
   */
  public void setTimeRange(int itmin, int itmax) {
    Check.argument(itmin<=itmax,"itmin<=itmax");
    _itmin = itmin;
    _itmax = itmax;
  }

  /**
   * Sets the min-max range of frequencies used to estimate the wavelet.
   * If the specified min-max bounds on frequency are not a subset of the
   * zero-Nyquist range [0,0.5], then no bandpass filter is used. The default
   * is to use no bandpass filter.
   * @param fmin minimum frequency, in cycles/sample.
   * @param fmax maximum frequency, in cycles/sample.
   */
  public void setFrequencyRange(double fmin, double fmax) {
    if (fmin<fmax && (0.0<fmin || fmax<0.5)) {
      _bpf = new BandPassFilter(fmin,fmax,0.05,0.01);
    } else {
      _bpf = null;
    }
  }

  /**
   * Sets the stability factor by which to scale zero-lag of correlations.
   * A factor slightly greater than one may stabilize estimates of
   * inverse wavelets A.
   * @param sfac stability factor.
   */
  public void setStabilityFactor(double sfac) {
    _sfac = sfac;
  }

  /**
   * Returns inverse wavelets a and b estimated by warping one sequence to
   * another. The sequences are related by warping such that f[t] ~ g[u[t]].
   * @param nab number of samples in the inverse wavelets a and b.
   * @param kab the sample index for a[0] and b[0].
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t]
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[][] getInverseAB(
    int nab, int kab, float[] u, float[] f, float[] g)
  {
    Check.argument(-nab<kab,"-nab<kab");
    Check.argument(kab<=0,"kab<=0");

    // Sequences WF and WSLG.
    float[][] s = new float[2*nab][];
    float[] lg = applyL(u,g); // L commutes with delay applied in loop
    for (int iab=0,ia=0,ib=nab,lag=kab; iab<nab; ++iab,++ia,++ib,++lag) {
      float[] df = delay(lag,f);
      float[] ldg = delay(lag,lg);
      float[] sldg = applyS(u,ldg); 
      s[ia] = applyW(df);
      s[ib] = applyW(sldg);
    }

    // The matrix Q and right-hand-side vector r, for Qe = r. For zero lag, we
    // have a0 = a[-ka] = 1, so that only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags in a
    // and all lags in b.
    int mab = 2*nab-1;
    DMatrix q = new DMatrix(mab,mab);
    DMatrix r = new DMatrix(mab,1);
    for (int iab=0,iq=0,ir=0; iab<2*nab; ++iab) {
      if (iab==-kab) continue; // skip lag zero, because a0 = 1
      for (int jab=0,jq=0; jab<2*nab; ++jab) {
        if (jab==-kab) continue; // skip lag zero, because a0 = 1
        double qij = dot(s[iab],s[jab]);
        q.set(iq,jq,qij);
        ++jq;
      }
      q.set(iq,iq,q.get(iq,iq)*_sfac);
      double ri = dot(s[iab],s[-kab]);
      r.set(ir,0,ri);
      ++iq;
      ++ir;
    }
    //System.out.println("q=\n"+q);
    //System.out.println("r=\n"+r);

    // Solve for inverse filters a and b using Cholesky decomposition of Q.
    DMatrixChd chd = new DMatrixChd(q);
    DMatrix e = chd.solve(r);
    float[][] ab = new float[2][nab];
    float[] a = ab[0];
    for (int iab=0,ie=0; iab<nab; ++iab) {
      if (iab==-kab) {
        a[iab] = 1.0f; // lag 0, so a0 = 1
      } else {
        a[iab] = -(float)e.get(ie,0);
        ++ie;
      }
    }
    float[] b = ab[1];
    for (int iab=0,ie=nab-1; iab<nab; ++iab,++ie) {
      b[iab] = (float)e.get(ie,0);
    }
    return ab;
  }

  /**
   * Estimates wavelets c and d from inverse wavelets a and b.
   * @param nab number of samples in inverse wavelets a and b.
   * @param kab the sample index for a[0] and b[0].
   * @param ab array {a,b} of coefficients for inverse wavelets a and b.
   * @param ncd number of samples in the wavelets c and d.
   * @param kcd the sample index for c[0] and d[0].
   * @return array {c,d} of wavelets c and d.
   */
  public float[][] getWaveletCD(
    int nab, int kab, float[][] ab, int ncd, int kcd) 
  {
    float[][] cd = new float[2][];
    cd[0] = inverse(nab,kab,ab[0],ncd,kcd);
    cd[1] = inverse(nab,kab,ab[1],ncd,kcd);
    return cd;
  }

  /**
   * Applies the specified inverse wavelet A.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyA(int na, int ka, float[] a, float[] f) {
    return convolve(na,ka,a,f);
  }

  /**
   * Applies the specified inverse wavelet B.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyB(int nb, int kb, float[] b, float[] f) {
    return convolve(nb,kb,b,f);
  }

  /**
   * Applies the specified wavelet C.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyC(int nc, int kc, float[] c, float[] f) {
    return convolve(nc,kc,c,f);
  }

  /**
   * Applies the specified wavelet D.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyD(int nd, int kd, float[] d, float[] f) {
    return convolve(nd,kd,d,f);
  }

  /**
   * Applies the weighting filter W, if any was specified. If none specified,
   * then this method simply returns a copy of the specified input sequence.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyW(float[] f) {
    int nt = f.length;
    float[] g = new float[nt];
    if (_bpf!=null) {
      _bpf.apply(f,g);
    } else {
      copy(f,g);
    }
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    for (int it=0; it<itlo; ++it)
      g[it] = 0.0f;
    for (int it=ithi+1; it<nt; ++it)
      g[it] = 0.0f;
    return g;
  }

  /**
   * Applies the low-pass anti-alias filter L.
   * If the specified warping includes squeezing, then this method attenuates
   * high frequencies that could be aliased during warping.
   * @param u array of warping times u(t).
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyL(float[] u, float[] f) {
    return aaf(RMAX,u,f);
  }

  /**
   * Applies the warping operator S.
   * Does not apply an anti-alias low-pass filter.
   * @param u array of warping times u(t).
   * @param f array with input sequence f(t).
   * @return array with warped output sequence.
   */
  public float[] applyS(float[] u, float[] f) {
    return warp(u,f);
  }

  /**
   * Applies the composite linear operator CSLB.
   * The sequence of operations is (1) convolution with the inverse wavelet b,
   * (2) anti-alias filtering (if squeezing), (3) warping, and (4) convolution
   * with the wavelet c.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index for b[0].
   * @param b array of coefficients for the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index for c[0].
   * @param c array of coefficients for the wavelet c.
   * @param u array[nt] of warping times u(t).
   * @param g array[nt] with input sequence.
   * @return array[nt] with output sequence.
   */
  public float[] applyCSLB(
    int nb, int kb, float[] b,
    int nc, int kc, float[] c,
    float[] u, float[] g) 
  {
    int nt = g.length;
    float[] bg = applyB(nb,kb,b,g);
    float[] lbg = applyL(u,bg);
    float[] sbg = applyS(u,lbg);
    float[] csbg = applyC(nc,kc,c,sbg);
    return csbg;
  }

  /**
   * Applies the composite linear operator WSLB.
   * The sequence of operations is (1) convolution with the inverse wavelet b,
   * (2) anti-alias filtering (if squeezing), (3) warping, and (4) application
   * of the weighting filter w.
   * @param nb number of samples in inverse wavelet b.
   * @param kb the sample index for b[0].
   * @param b array of coefficients for inverse wavelet b.
   * @param u array[nt] of warping times u(t).
   * @param g array[nt] with input sequence.
   * @return array[nt] with output sequence.
   */
  public float[] applyWSLB(
    int nb, int kb, float[] b, float[] u, float[] g) 
  {
    int nt = g.length;
    float[] bg = applyB(nb,kb,b,g);
    float[] lbg = applyL(u,bg);
    float[] sbg = applyS(u,lbg);
    float[] wsbg = applyW(sbg);
    return wsbg;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float RMAX = 10.0f; // limits anti-alias filter
  private static final SincInterpolator _si =
    SincInterpolator.fromErrorAndFrequency(0.01,0.40);

  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private BandPassFilter _bpf;

  /**
   * Returns the array of lagged WF and WSLG.
   */
  private float[][] computeLagged(
    int nab, int kab, float[] u, float[] f, float[] g)
  {
    g = applyL(u,g);
    int nt = u.length;
    float[][] d = new float[2*nab][];
    for (int iab=0,ia=0,ib=nab,lag=kab; iab<nab; ++iab,++ia,++ib,++lag) {
      float[] df = delay(lag,f);
      float[] dg = delay(lag,g);
      float[] sdg = applyS(u,dg);
      d[ia] = applyW(df);
      d[ib] = applyW(sdg);
    }
    return d;
  }

  private double dot(float[] x, float[] y) {
    int nt = x.length;
    double sum = 0.0;
    for (int it=0; it<nt; ++it) 
      sum += x[it]*y[it];
    return sum;
  }

  /**
   * Returns the largest squeezing r(t) = u'(t) not greater than rmax.
   * If less than or equal to one, then no squeezing is implied by u(t).
   */
  private float squeezing(float rmax, float[] u) {
    int nt = u.length;
    int itlo = max(1,_itmin);
    int ithi = min(_itmax,nt-1);
    float r = 0.0f;
    for (int it=itlo; it<=ithi; ++it) {
      float du = u[it]-u[it-1];
      if (r<du)
        r = du;
    }
    return min(r,rmax);
  }

  /**
   * If necessary, applies an anti-alias filter to the sequence x(t).
   * An anti-alias filter is necessary if the warping includes squeezing.
   */
  private float[] aaf(float rmax, float[] u, float[] x) {
    int nt = x.length;
    float r = squeezing(RMAX,u);
    if (r>1.0) {
      float[] y = new float[nt];
      BandPassFilter aaf = new BandPassFilter(0.0,0.5/r,0.10/r,0.01);
      aaf.apply(x,y);
      return y;
    } else {
      return copy(x);
    }
  }

  /**
   * Returns y(t) = x(t-lag).
   */
  private static float[] delay(int lag, float[] x) {
    int nt = x.length;
    int itlo = max(0,lag);   // 0 <= it-lag
    int ithi = min(nt,nt+lag); // it-lag < nt
    float[] y = new float[nt];
    for (int it=0; it<itlo; ++it)
      y[it] = 0.0f;
    for (int it=itlo; it<ithi; ++it)
      y[it] = x[it-lag];
    for (int it=ithi; it<nt; ++it)
      y[it] = 0.0f;
    return y;
  }

  /**
   * Returns y(t) = x(u(t)).
   */
  private static float[] warp(float[] u, float[] x) {
    int nt = u.length;
    float[] y = new float[nt];
    _si.interpolate(nt,1.0,0.0,x,nt,u,y);
    y[0] *= u[1]-u[0];
    for (int it=1; it<nt; ++it)
      y[it] *= u[it]-u[it-1];
    return y;
  }

  /**
   * Returns y(t) = h(t)*x(t), where * denotes convolution.
   */
  private static float[] convolve(int nh, int kh, float[] h, float[] x) {
    int nt = x.length;
    float[] y = new float[nt];
    convolve(nh,kh,h,x,y);
    return y;
  }
  private static void convolve(
    int nh, int kh, float[] h, float[] f,  float[] g)
  {
    int nt = f.length;
    conv(nh,kh,h,nt,0,f,nt,0,g);
  }

  /**
   * Returns inverse c(t) of a(t), such that c(t)*a(t) = delta(t).
   */
  private float[] inverse(int na, int ka, float[] a, int nc, int kc) {
    float[] one = {1.0f};
    float[] ca1 = new float[nc];
    float[] caa = new float[nc];
    xcor(na,ka,a,1,0,one,nc,kc,ca1);
    xcor(na,ka,a,na,ka,a,nc, 0,caa);
    caa[0] *= _sfac;
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }

  private float rms(float[] x) {
    int nt = x.length;
    return (float)sqrt(dot(x,x)/nt);
  }
}
