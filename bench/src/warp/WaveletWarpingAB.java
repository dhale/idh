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

  public void setStabilityFactor(double sfac) {
    _sfac = (float)sfac;
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
  public float[][] getInverseABQ(
    int na, int ka, int nb, int kb, 
    float[] u, float[] f, float[] g)
  {
    // Matrix Q = [F | -SG]'[F | -SG].
    int nc = na+nb;
    DMatrix q = new DMatrix(nc,nc);
    for (int ia=0,ic=0,ilag=ka; ia<na; ++ia,++ic,++ilag) {
      float[] fi = delay(ilag,f);
      for (int ja=0,jc=0,jlag=ka; ja<na; ++ja,++jc,++jlag) {
        float[] fj = delay(jlag,f);
        q.set(ic,jc,dot(fi,fj));
      }
      for (int jb=0,jc=na,jlag=kb; jb<nb; ++jb,++jc,++jlag) {
        float[] gj = warp(u,delay(jlag,g));
        q.set(ic,jc,-dot(fi,gj));
      }
    }
    for (int ib=0,ic=na,ilag=kb; ib<nb; ++ib,++ic,++ilag) {
      float[] gi = warp(u,delay(ilag,g));
      for (int ja=0,jc=0,jlag=ka; ja<na; ++ja,++jc,++jlag) {
        float[] fj = delay(jlag,f);
        q.set(ic,jc,-dot(gi,fj));
      }
      for (int jb=0,jc=na,jlag=kb; jb<nb; ++jb,++jc,++jlag) {
        float[] gj = warp(u,delay(jlag,g));
        q.set(ic,jc,dot(gi,gj));
      }
    }

    // Get coefficients a and b from eigenvector for smallest eigenvalue.
    return abFromQ(na,nb,q);
  }
  private static float[][] abFromQ(int na, int nb, DMatrix q) {
    int nc = na+nb;
    DMatrixEvd evd = new DMatrixEvd(q);
    DMatrix v = evd.getV().get(0,nc-1,0,0);
    double[] ev = evd.getRealEigenvalues();
    trace("edelta = "+(ev[1]-ev[0])/(FLT_EPSILON*ev[nc-1]));
    //dump(evd.getRealEigenvalues());
    float[] a = new float[na];
    float[] b = new float[nb];
    for (int ia=0,ic=0; ia<na; ++ia,++ic)
      a[ia] = (float)v.get(ic,0);
    for (int ib=0,ic=na; ib<nb; ++ib,++ic)
      b[ib] = (float)v.get(ic,0);
    return new float[][]{a,b};
  }

  public float[][] getInverseAB(
    int na, int ka, int nb, int kb, 
    float[] u, float[] f, float[] g)
  {
    trace("getInverseAB: begin");
    DMatrix aa = new DMatrix(na,na);
    for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
      float[] fi = delay(ilag,f);
      for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
        if (ia<=ja) {
          float[] fj = delay(jlag,f);
          aa.set(ia,ja,dot(fi,fj));
        } else {
          aa.set(ia,ja,aa.get(ja,ia));
        }
      }
    }
    DMatrix bb = new DMatrix(nb,nb);
    for (int ib=0,ilag=kb; ib<nb; ++ib,++ilag) {
      float[] gi = warp(u,delay(ilag,g));
      for (int jb=0,jlag=kb; jb<nb; ++jb,++jlag) {
        if (ib<=jb) {
          float[] gj = warp(u,delay(jlag,g));
          bb.set(ib,jb,dot(gi,gj));
        } else {
          bb.set(ib,jb,bb.get(jb,ib));
        }
      }
    }
    DMatrix ab = new DMatrix(na,nb);
    DMatrix ba = new DMatrix(nb,na);
    for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
      float[] fi = delay(ilag,f);
      for (int jb=0,jlag=kb; jb<nb; ++jb,++jlag) {
        float[] gj = warp(u,delay(jlag,g));
        double abij = -dot(fi,gj);
        ab.set(ia,jb,abij);
        ba.set(jb,ia,abij);
      }
    }
    float[][] abs = abSolve(aa,ab,ba,bb);
    trace("getInverseAB: done");
    return abs;
  }
  public float[][] getInverseAB(
    int na, int ka, int nb, int kb, 
    float[][] u, float[][] f, float[][] g)
  {
    trace("getInverseAB: begin");
    DMatrix aa = new DMatrix(na,na);
    for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
      float[][] fi = delay(ilag,f);
      for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
        if (ia<=ja) {
          float[][] fj = delay(jlag,f);
          aa.set(ia,ja,dot(fi,fj));
        } else {
          aa.set(ia,ja,aa.get(ja,ia));
        }
      }
    }
    DMatrix bb = new DMatrix(nb,nb);
    for (int ib=0,ilag=kb; ib<nb; ++ib,++ilag) {
      float[][] gi = warp(u,delay(ilag,g));
      for (int jb=0,jlag=kb; jb<nb; ++jb,++jlag) {
        if (ib<=jb) {
          float[][] gj = warp(u,delay(jlag,g));
          bb.set(ib,jb,dot(gi,gj));
        } else {
          bb.set(ib,jb,bb.get(jb,ib));
        }
      }
    }
    DMatrix ab = new DMatrix(na,nb);
    DMatrix ba = new DMatrix(nb,na);
    for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
      float[][] fi = delay(ilag,f);
      for (int jb=0,jlag=kb; jb<nb; ++jb,++jlag) {
        float[][] gj = warp(u,delay(jlag,g));
        double abij = -dot(fi,gj);
        ab.set(ia,jb,abij);
        ba.set(jb,ia,abij);
      }
    }
    float[][] abc = abSolve(aa,ab,ba,bb);
    trace("getInverseAB: done");
    return abc;
  }
  private static float[][] abSolve(
    DMatrix aa, DMatrix ab, DMatrix ba, DMatrix bb)
  {
    assert aa.isSymmetric();
    assert bb.isSymmetric();
    assert ab.equals(ba.transpose());
    aa.set(0,0,4.0);
    bb.set(0,0,9.0);
    ab.set(0,0,-6.0);
    ba.set(0,0,-6.0);
    int na = aa.getN();
    int nb = bb.getN();
    trace("na="+na+" nb="+nb+" aa,bb,ab:");
    trace(aa.toString());
    trace(bb.toString());
    trace(ab.toString());
    DMatrix ai = DMatrix.identity(na,na);
    DMatrix bi = DMatrix.identity(nb,nb);
    DMatrix aam,bbm;
    DMatrixEvd aevd = null;
    DMatrixEvd bevd = null;
    double ea = 0.0; // eigenvalue associated with a
    double eb = 0.0; // eigenvalue associated with b
    double ebmax = bb.normF()/nb;
    double ebmin = -ebmax;
    int neb = 101;
    double deb = (ebmax-ebmin)/(neb-1);
    for (int ieb=0; ieb<neb; ++ieb) {
      double ebi = ebmin+ieb*deb;
      bbm = bb.minus(bi.times(ebi)).inverse();
      aam = aa.minus(ab.times(bbm.times(ba)));
      aevd = new DMatrixEvd(aam);
      double[] eas = aevd.getRealEigenvalues();
      for (int ia=0; ia<na; ++ia) {
        double eai = eas[ia];
        aam = aa.minus(ai.times(eai)).inverse();
        bbm = bb.minus(ba.times(aam.times(ab)));
        bevd = new DMatrixEvd(bbm);
        double[] ebs = bevd.getRealEigenvalues();
        for (int ib=0; ib<nb; ++ib) {
        }
      }
      bbm = bb.minus(bi.times(ebi)).inverse();
      aam = aa.minus(ab.times(bbm.times(ba)));
      aevd = new DMatrixEvd(aam);
      double eai = aevd.getRealEigenvalues()[0];
      trace("ea="+eai+" eb="+ebi+" sum="+(eai+ebi));
      if (eai+ebi<ea+eb) {
        ea = eai;
        eb = ebi;
      }
    }
    bbm = bb.minus(bi.times(eb)).inverse();
    aam = aa.minus(ab.times(bbm.times(ba)));
    aevd = new DMatrixEvd(aam);
    aam = aa.minus(ai.times(ea)).inverse();
    bbm = bb.minus(ba.times(aam.times(ab)));
    bevd = new DMatrixEvd(bbm);
    ea = aevd.getRealEigenvalues()[0];
    eb = bevd.getRealEigenvalues()[0];
    DMatrix av = aevd.getV().get(0,na-1,0,0);
    DMatrix bv = bevd.getV().get(0,nb-1,0,0);
    DMatrix aerr = aa.times(av).minus(ab.times(bv).minus(av.times(ea)));
    DMatrix berr = bb.times(bv).minus(ba.times(av).minus(bv.times(eb)));
    trace("aerr=\n"+aerr);
    trace("berr=\n"+berr);
    float[] a = new float[na];
    float[] b = new float[nb];
    for (int ia=0; ia<na; ++ia)
      a[ia] = (float)av.get(ia,0);
    for (int ib=0; ib<nb; ++ib)
      b[ib] = (float)bv.get(ib,0);
    return new float[][]{a,b};
  }

  /**
   * Estimates the inverse h of a specified wavelet a.
   * @param na number of samples in the wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   */
  public float[] getInverse(int na, int ka, float[] a, int nh, int kh) {
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);
    caa[0] *= 1.0+_sfac;
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }

  /**
   * Convolves with the specified filter h.
   * @param nh number of samples in h.
   * @param kh the sample index for h[0].
   * @param h input array of coefficients in h.
   * @param f input array with input sequence f.
   * @return output array with filtered sequence h*f.
   */
  public static float[] convolve(int nh, int kh, float[] h, float[] f) {
    int n = f.length;
    float[] g = new float[n];
    convolve(nh,kh,h,f,g);
    return g;
  }
  public static float[][] convolve(int nh, int kh, float[] h, float[][] f) {
    int n = f.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      g[i] = convolve(nh,kh,h,f[i]);
    return g;
  }

  /**
   * Applies the warping operator S.
   * Does not apply an anti-alias low-pass filter.
   * @param u array of warping times u(t).
   * @param f array with input sequence f(t).
   * @return array with warped output sequence.
   */
  public float[] warp(float[] u, float[] x) {
    return _wf.apply(u,x);
  }
  public float[][] warp(float[][] u, float[][] x) {
    int n = u.length;
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = warp(u[i],x[i]);
    return y;
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
    float[] af = convolve(na,ka,a,f);
    float[] saf = warp(u,af);
    float[] hsaf = convolve(nh,kh,h,saf);
    return hsaf;
  }
  public float[][] applyHSA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[][] u, float[][] f) 
  {
    int nt = f.length;
    float[][] af = convolve(na,ka,a,f);
    float[][] saf = warp(u,af);
    float[][] hsaf = convolve(nh,kh,h,saf);
    return hsaf;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float RMAX = 10.0f; // limits anti-alias filter
  private static final SincInterpolator _si =
    SincInterpolator.fromErrorAndFrequency(0.01,0.40);
  private static final WarpingFilter _wf = new WarpingFilter();
  static {
    _wf.setOversamplingLimit(3.0);
  }

  private float _sfac = 0.0f;
  private int _itmin = -1;
  private int _itmax = -1;

  private double dot(float[] x, float[] y) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it];
    return sum/(1+ithi-itlo);
  }
  private double dot(float[][] x, float[][] y) {
    int n = x.length;
    double sum = 0.0;
    for (int i=0; i<n; ++i) 
      sum += dot(x[i],y[i]);
    return sum/n;
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
   * Computes y(t) = h(t)*x(t), where * denotes convolution.
   */
  private static void convolve(
    int nh, int kh, float[] h, float[] f,  float[] g)
  {
    int nt = f.length;
    conv(nh,kh,h,nt,0,f,nt,0,g);
  }

  private float rms(float[] x) {
    return (float)sqrt(dot(x,x));
  }
  public float rms(float[][] x) {
    return (float)sqrt(dot(x,x));
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
