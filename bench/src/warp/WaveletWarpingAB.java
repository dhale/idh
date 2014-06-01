/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;

import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates a wavelet from alignment by warping of sequences or images.
 * The two sequences or images are assumed to have been convolved with the
 * same wavelet. Warping of one sequence or image to align with the other will
 * cause the wavelet to be stretched or squeezed, and this distortion enables
 * us to estimate the wavelet.
 * <p>
 * For images, convolution with the wavelet is assumed to be in only the 1st
 * dimension. For definiteness, this 1st dimension is assumed to be time in
 * the documentation below.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.04.19
 */
public class WaveletWarpingAB {

  /**
   * Sets the min-max range of times used to estimate wavelet.
   * @param itmin minimum time, in samples.
   * @param itmax maximum time, in samples.
   */
  public void setTimeRange(int itmin, int itmax) {
    _itmin = itmin;
    _itmax = itmax;
  }

  /**
   * Returns inverse wavelets a and b estimated by warping.
   * The two sequences f and g are related by f[t] ~ g[u[t]].
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index for b[0].
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t]
   * @return array {a,b} of coefficients for inverse wavelets a and b.
   */
  public float[][] getInverseAB(
    int na, int ka, int nb, int kb, 
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;

    // Matrix Q = [F  -SLG]'[F  -SLG].
    int nc = na+nb;
    DMatrix q = new DMatrix(nc,nc);
    for (int ia=0,ic=0,ilag=ka; ia<na; ++ia,++ic,++ilag) {
      float[] fi = delay(ilag,f);
      for (int ja=0,jc=0,jlag=ka; ja<na; ++ja,++jc,++jlag) {
        float[] fj = delay(jlag,f);
        q.set(ic,jc,dot(fi,fj));
      }
      for (int jb=0,jc=na,jlag=kb; jb<nb; ++jb,++jc,++jlag) {
        float[] gj = applyS(u,delay(jlag,g));
        q.set(ic,jc,-dot(fi,gj));
      }
    }
    for (int ib=0,ic=na,ilag=kb; ib<nb; ++ib,++ic,++ilag) {
      float[] gi = applyS(u,delay(ilag,g));
      for (int ja=0,jc=0,jlag=ka; ja<na; ++ja,++jc,++jlag) {
        float[] fj = delay(jlag,f);
        q.set(ic,jc,-dot(gi,fj));
      }
      for (int jb=0,jc=na,jlag=kb; jb<nb; ++jb,++jc,++jlag) {
        float[] gj = applyS(u,delay(jlag,g));
        q.set(ic,jc,dot(gi,gj));
      }
    }

    // Get coefficients a and b from eigenvector for smallest eigenvalue.
    DMatrixEvd evd = new DMatrixEvd(q);
    DMatrix v = evd.getV().get(0,nc-1,0,0);
    float[] a = new float[na];
    float[] b = new float[nb];
    for (int ia=0,ic=0; ia<na; ++ia,++ic)
      a[ia] = (float)v.get(ic,0);
    for (int ib=0,ic=na; ib<nb; ++ib,++ic)
      b[ib] = (float)v.get(ic,0);

    return new float[][]{a,b};
  }

  /**
   * Estimates the wavelet h from the inverse wavelet a.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   */
  public float[] getWaveletH(int na, int ka, float[] a, int nh, int kh) {
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);
    caa[0] *= _sfac;
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }

  /**
   * Applies the specified inverse wavelet A.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyA(int na, int ka, float[] a, float[] f) {
    return convolve(na,ka,a,f);
  }
  public float[][] applyA(int na, int ka, float[] a, float[][] f) {
    return convolve(na,ka,a,f);
  }

  /**
   * Applies the specified wavelet H.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyH(int nh, int kh, float[] h, float[] f) {
    return convolve(nh,kh,h,f);
  }
  public float[][] applyH(int nh, int kh, float[] h, float[][] f) {
    return convolve(nh,kh,h,f);
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
  public float[][] applyS(float[][] u, float[][] f) {
    return warp(u,f);
  }

  /**
   * Applies the composite linear operator HSA.
   * The sequence of operations is (1) convolution with the inverse wavelet a,
   * (2) warping (with anti-alias filtering), and (3) convolution with the
   * wavelet h.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param u array[nt] of warping times u(t).
   * @param f array[nt] with input sequence.
   * @return array[nt] with output sequence.
   */
  public float[] applyHSA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[] u, float[] f) 
  {
    int nt = f.length;
    float[] af = applyA(na,ka,a,f);
    float[] saf = applyS(u,af);
    float[] hsaf = applyH(nh,kh,h,saf);
    return hsaf;
  }
  public float[][] applyHSA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[][] u, float[][] f) 
  {
    int nt = f.length;
    float[][] af = applyA(na,ka,a,f);
    float[][] saf = applyS(u,af);
    float[][] hsaf = applyH(nh,kh,h,saf);
    return hsaf;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float RMAX = 10.0f; // limits anti-alias filter
  private static final SincInterpolator _si =
    SincInterpolator.fromErrorAndFrequency(0.01,0.40);
  private static final WarpingFilter _wf = new WarpingFilter();

  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;

  private double dot(float[] x, float[] y) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it];
    return sum;
  }
  private double dot(float[][] x, float[][] y) {
    int n = x.length;
    double sum = 0.0;
    for (int i=0; i<n; ++i) 
      sum += dot(x[i],y[i]);
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
  private float squeezing(float rmax, float[][] u) {
    int n = u.length;
    float r = 0.0f;
    for (int i=0; i<n; ++i)
      r = max(r,squeezing(rmax,u[i]));
    return r;
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
  private float[][] aaf(float rmax, float[][] u, float[][] x) {
    float r = squeezing(RMAX,u);
    if (r>1.0) {
      int nx = x.length;
      int nt = x[0].length;
      float[][] y = new float[nx][nt];
      BandPassFilter aaf = new BandPassFilter(0.0,0.5/r,0.10/r,0.01);
      for (int ix=0; ix<nx; ++ix)
        aaf.apply(x[ix],y[ix]);
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
  private static float[][] delay(int lag, float[][] x) {
    int n = x.length;
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = delay(lag,x[i]);
    return y;
  }

  /**
   * Returns y(t) = x(u(t)).
   */
  private float[] warp(float[] u, float[] x) {
    return _wf.apply(u,x);
  }
  private float[][] warp(float[][] u, float[][] x) {
    int n = u.length;
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = warp(u[i],x[i]);
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
  private static float[][] convolve(int nh, int kh, float[] h, float[][] x) {
    int n = x.length;
    int nt = x[0].length;
    float[][] y = new float[n][nt];
    for (int i=0; i<n; ++i)
      convolve(nh,kh,h,x[i],y[i]);
    return y;
  }

  private float rms(float[] x) {
    return (float)sqrt(dot(x,x)/x.length);
  }
  public float rms(float[][] x) {
    return (float)sqrt(dot(x,x)/x.length/x[0].length);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
