/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.opt.BrentMinFinder;

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
 * @version 2015.02.03
 */
public class WaveletWarpingAH {

  /**
   * Sets the min-max range of times used to estimate wavelet.
   * @param itmin minimum time, in samples.
   * @param itmax maximum time, in samples.
   */
  public void setTimeRange(int itmin, int itmax) {
    _itmin = itmin;
    _itmax = itmax;
  }

  public void setStabilityFactor(double factor) {
    _stabilityFactor = (float)factor;
  }

  public void setMaxIterations(int maxiter) {
    _maxiter = maxiter;
  }

  /**
   * Return inverse wavelet a estimated by warping.
   * The two sequences f and g are related by f[t] ~ g[u[t]].
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t]
   * @return array a of coefficients for inverse wavelet a.
   */
  public float[] getInverse(
    int na, int ka, int nh, int kh, float[] h,
    float[] u, float[] f, float[] g)
  {
    DMatrix q = new DMatrix(na,na);
    DMatrix r = new DMatrix(na,1);
    for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
      float[] gi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
      for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
        float[] gj = convolve(nh,kh,h,warp(u,delay(jlag,g)));
        if (ilag==0 && jlag==0)
          q.set(ia,ja,dot(gi,gj));
        else if (ilag==0 || jlag==0)
          q.set(ia,ja,0.0);
        else
          q.set(ia,ja,dot(gi,gj));
      }
      double s = 0.0001*abs(ilag); // for stability, penalize larger lags
      q.set(ia,ia,q.get(ia,ia)*(1.0+s*s));
      if (ilag==0)
        r.set(ia,0,q.get(ia,ia));
      else
        r.set(ia,0,dot(gi,f));
    }
    return solveChd(q,r);
  }
  public float[] getInverse(
    int na, int ka, int nh, int kh, float[] h,
    float[][] u, float[][] f, float[][] g)
  {
    DMatrix q = new DMatrix(na,na);
    DMatrix r = new DMatrix(na,1);
    for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
      float[][] gi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
      for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
        float[][] gj = convolve(nh,kh,h,warp(u,delay(jlag,g)));
        q.set(ia,ja,dot(gi,gj));
      }
      r.set(ia,0,dot(gi,f));
    }
    return solveChd(q,r);
  }

  public float[] getWavelet(
    int nh, int kh, int na, int ka, float[] a,
    float[] u, float[] f, float[] g)
  {
    float[] p = warp(u,convolve(na,ka,a,g));
    DMatrix q = new DMatrix(nh,nh);
    DMatrix r = new DMatrix(nh,1);
    for (int ih=0,ilag=kh; ih<nh; ++ih,++ilag) {
      float[] pi = delay(ilag,p);
      for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
        float[] pj = delay(jlag,p);
        q.set(ih,jh,dot(pi,pj));
      }
      q.set(ih,ih,q.get(ih,ih)*(1.0+_stabilityFactor));
      r.set(ih,0,dot(pi,f));
    }
    return solveChd(q,r);
  }
  public float[] getWavelet(
    int nh, int kh, int na, int ka, float[] a,
    float[][] u, float[][] f, float[][] g)
  {
    float[][] p = warp(u,convolve(na,ka,a,g));
    DMatrix q = new DMatrix(nh,nh);
    DMatrix r = new DMatrix(nh,1);
    for (int ih=0,ilag=kh; ih<nh; ++ih,++ilag) {
      float[][] pi = delay(ilag,p);
      for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
        float[][] pj = delay(jlag,p);
        q.set(ih,jh,dot(pi,pj));
      }
      r.set(ih,0,dot(pi,f));
    }
    return solveChd(q,r);
  }

  /**
   * Estimates the wavelet h from a specified inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   */
  public float[] getWavelet(int nh, int kh, int na, int ka, float[] a) {
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);
    caa[0] *= 1.0+_stabilityFactor;
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }
  public float[] getInverse(int na, int ka, int nh, int kh, float[] h) {
    return getWavelet(na,ka,nh,kh,h);
  }

  public float[][] getAHSimple(
    int na, int ka, int nh, int kh,
    float[] u, float[] f, float[] g)
  {
    trace("getAHSimple: begin");
    int namax = na;
    int kamin = ka;
    na = 1;
    ka = 0;
    float[] a = new float[na]; a[0] = 1.0f;
    float[] h = getWavelet(nh,kh,na,ka,a,u,f,g);
    float[] r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
    double rmsr = rms(r);
    double rmsrOld = 2.0*rmsr;
    trace("getAHSimple: rmsr="+rmsr);
    for (int jouter=0; jouter<_maxiter && rmsr<rmsrOld-0.01*rmsr; ++jouter) {
      na += 2;
      ka -= 1;
      if (na>namax) na = namax;
      if (ka<kamin) ka = kamin;
      for (int jinner=0; jinner<4; ++jinner) {
        a = getInverse(na,ka,nh,kh,h,u,f,g);
        h = getWavelet(nh,kh,na,ka,a,u,f,g);
      }
      r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
      rmsrOld = rmsr;
      rmsr = rms(r);
      trace("getAHSimple: jouter="+jouter+" na="+na+" rmsr="+rmsr);
    }
    trace("getAHSimple: end");
    return new float[][]{a,h};
  }

  public float[][] getAHNewton(
    int na, int ka, int nh, int kh,
    float[] u, float[] f, float[] g)
  {
    trace("getAHNewton: begin");
    float[] a = new float[na];
    a[-ka] = 1.0f;
    float[] h = getWavelet(nh,kh,na,ka,a,u,f,g);
    float[] r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
    double rmsr = rms(r);
    double rmsrOld = 2.0*rmsr;
    trace("getAHNewton: rmsr="+rmsr);

    // Newton iterations.
    for (int iter=0; iter<_maxiter && rmsr<rmsrOld-0.01*rmsr; ++iter) {

      // Components of matrix and right-hand-side.
      DMatrix pp = new DMatrix(na,na);
      DMatrix pr = new DMatrix(na,1);
      for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
        float[] pi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
        for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
          if (ia<=ja) {
            float[] pj = convolve(nh,kh,h,warp(u,delay(jlag,g)));
            pp.set(ia,ja,dot(pi,pj));
          } else {
            pp.set(ia,ja,pp.get(ja,ia));
          }
        }
        pr.set(ia,0,dot(pi,r));
      }
      float[] q0 = warp(u,convolve(na,ka,a,g));
      DMatrix qq = new DMatrix(nh,nh);
      DMatrix qr = new DMatrix(nh,1);
      for (int ih=0,ilag=kh; ih<nh; ++ih,++ilag) {
        float[] qi = delay(ilag,q0);
        for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
          if (ih<=jh) {
            float[] qj = delay(jlag,q0);
            qq.set(ih,jh,dot(qi,qj));
          } else {
            qq.set(ih,jh,qq.get(jh,ih));
          }
        }
        qr.set(ih,0,dot(qi,r));
      }
      DMatrix pq = new DMatrix(na,nh);
      DMatrix qp = new DMatrix(nh,na);
      for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
        float[] sgi = warp(u,delay(ilag,g));
        float[] pi = convolve(nh,kh,h,sgi);
        for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
          float[] qj = delay(jlag,q0);
          float[] sgij = delay(jlag,sgi);
          double pqij = dot(pi,qj); //+dot(sgij,r);
          pq.set(ia,jh,pqij);
          qp.set(jh,ia,pqij);
        }
      }

      // Composite matrix and right-hand-side.
      int nm = na+nh;
      DMatrix m = new DMatrix(nm,nm);
      DMatrix b = new DMatrix(nm,1);
      m.set( 0,na-1,0,na-1,pp); m.set( 0,na-1,na,nm-1,pq);
      m.set(na,nm-1,0,na-1,qp); m.set(na,nm-1,na,nm-1,qq);
      b.set( 0,na-1,0,0,pr);
      b.set(na,nm-1,0,0,qr);

      // Ensure that a = 1 for lag = 0.
      for (int im=0; im<nm; ++im)
        if (im!=-ka)
          m.set(im,-ka,0.0);
      for (int im=0; im<nm; ++im)
        if (im!=-ka)
          m.set(-ka,im,0.0);
      b.set(-ka,0,0.0);

      // Ensure positive-definite.
      for (int im=0; im<nm; ++im)
        m.set(im,im,m.get(im,im)*(1.0+_stabilityFactor));

      // Solve SPD system of equations.
      DMatrix x = m.solve(b);
      float[] da = new float[na];
      float[] dh = new float[nh];
      for (int ia=0; ia<na; ++ia)
        da[ia] = (float)x.get(ia,0);
      for (int ih=0; ih<nh; ++ih)
        dh[ih] = (float)x.get(ih+na,0);

      // Update a and h after line search for good step size.
      float step = findStep(na,ka,a,da,nh,kh,h,dh,u,f,g);
      //float step = 1.0f;
      a = sub(a,mul(step,da));
      h = sub(h,mul(step,dh));

      // Update rms of residual r.
      r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
      rmsrOld = rmsr;
      rmsr = rms(r);
      trace("getAHNewton: step="+step+" rmsr="+rmsr);
    }
    trace("getAHNewton: end");
    return new float[][]{a,h};
  }

  public float[][] getAHGaussNewton(
    int na, int ka, int nh, int kh,
    float[] u, float[] f, float[] g)
  {
    trace("getAHGaussNewton: begin");

    // Begin with a = unit impulse and h = shaping filter.
    float[] a = new float[na];
    a[-ka] = 1.0f;
    float[] h = getWavelet(nh,kh,na,ka,a,u,f,g);
    //trace("getAHGaussNewton: h="); dump(h);

    // Initial residual.
    float[] r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
    float rmsr = (float)rms(r);
    float rmsrOld = 2.0f*rmsr;
    trace("getAHGaussNewton: rmsr="+rmsr);

    // Gauss-Newton iterations.
    for (int iter=0; iter<_maxiter && rmsr<=rmsrOld*0.99999; ++iter) {
    //for (int iter=0; iter<_maxiter; ++iter) {

      // Components of matrix and right-hand-side.
      DMatrix pp = new DMatrix(na,na);
      DMatrix pr = new DMatrix(na,1);
      for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
        float[] pi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
        for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
          if (ia<=ja) {
            float[] pj = convolve(nh,kh,h,warp(u,delay(jlag,g)));
            pp.set(ia,ja,dot(pi,pj));
          } else {
            pp.set(ia,ja,pp.get(ja,ia));
          }
        }
        pr.set(ia,0,dot(pi,r));
      }
      float[] q0 = warp(u,convolve(na,ka,a,g));
      DMatrix qq = new DMatrix(nh,nh);
      DMatrix qr = new DMatrix(nh,1);
      for (int ih=0,ilag=kh; ih<nh; ++ih,++ilag) {
        float[] qi = delay(ilag,q0);
        for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
          if (ih<=jh) {
            float[] qj = delay(jlag,q0);
            qq.set(ih,jh,dot(qi,qj));
          } else {
            qq.set(ih,jh,qq.get(jh,ih));
          }
        }
        qr.set(ih,0,dot(qi,r));
      }
      DMatrix pq = new DMatrix(na,nh);
      DMatrix qp = new DMatrix(nh,na);
      for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
        float[] pi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
        for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
          float[] qj = delay(jlag,q0);
          double piqj = dot(pi,qj);
          pq.set(ia,jh,piqj);
          qp.set(jh,ia,piqj);
        }
      }

      // Composite matrix and right-hand-side.
      int nm = na+nh;
      DMatrix m = new DMatrix(nm,nm);
      DMatrix b = new DMatrix(nm,1);
      m.set( 0,na-1,0,na-1,pp); m.set( 0,na-1,na,nm-1,pq);
      m.set(na,nm-1,0,na-1,qp); m.set(na,nm-1,na,nm-1,qq);
      b.set( 0,na-1,0,0,pr);
      b.set(na,nm-1,0,0,qr);

      // Ensure that a = 1 for lag = 0.
      for (int im=0; im<nm; ++im)
        if (im!=-ka)
          m.set(im,-ka,0.0);
      for (int im=0; im<nm; ++im)
        if (im!=-ka)
          m.set(-ka,im,0.0);
      b.set(-ka,0,0.0);

      // Ensure positive-definite.
      for (int im=0; im<nm; ++im)
        m.set(im,im,m.get(im,im)*(1.0+_stabilityFactor));
      //trace("m=\n"+m.times(301));
      //trace("b=\n"+b.times(301));

      // Solve SPD system of equations.
      DMatrix x = m.chd().solve(b);
      float[] da = new float[na];
      float[] dh = new float[nh];
      for (int ia=0; ia<na; ++ia)
        da[ia] = (float)x.get(ia,0);
      for (int ih=0; ih<nh; ++ih)
        dh[ih] = (float)x.get(ih+na,0);

      // Update a and h after line search for good step size.
      float step = findStep(na,ka,a,da,nh,kh,h,dh,u,f,g);
      a = sub(a,mul(step,da));
      h = sub(h,mul(step,dh));

      // Update rms of residual r.
      r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
      rmsrOld = rmsr;
      rmsr = (float)rms(r);
      float bnrm = (float)b.norm2();
      trace("  iter="+iter+" step="+step+" rmsr="+rmsr+" bnrm="+bnrm);
    }
    trace("getAHGaussNewton: end");
    return new float[][]{a,h};
  }
  private float findStep(
    final int na, final int ka, final float[] a, final float[] da,
    final int nh, final int kh, final float[] h, final float[] dh,
    final float[] u, final float[] f, final float[] g)
  {
    BrentMinFinder bmf = new BrentMinFinder(
      new BrentMinFinder.Function() {
        public double evaluate(double step) {
          //trace("bmf: evaluate step="+step);
          float s = (float)step;
          float[] an = sub(a,mul(s,da));
          float[] hn = sub(h,mul(s,dh));
          float[] rn = sub(convolve(nh,kh,hn,warp(u,convolve(na,ka,an,g))),f);
          return rms(rn);
        }
      });
    return (float)bmf.findMin(0.0,1.0,0.0001);
  }

  public float[][] getAHGaussNewton(
    int na, int ka, int nh, int kh,
    float[][] u, float[][] f, float[][] g)
  {
    trace("getAHGaussNewton: begin");
    float[] a = new float[na];
    a[-ka] = 1.0f;
    float[] h = getWavelet(nh,kh,na,ka,a,u,f,g);
    float[][] r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
    double rmsr = rms(r);
    double rmsrOld = 2.0*rmsr;
    trace("getAHGaussNewton: rmsr="+rmsr);

    // Gauss-Newton iterations.
    for (int iter=0; iter<_maxiter && rmsr<rmsrOld-0.001*rmsr; ++iter) {

      // Components of matrix and right-hand-side.
      DMatrix pp = new DMatrix(na,na);
      DMatrix pr = new DMatrix(na,1);
      for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
        float[][] pi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
        for (int ja=0,jlag=ka; ja<na; ++ja,++jlag) {
          if (ia<=ja) {
            float[][] pj = convolve(nh,kh,h,warp(u,delay(jlag,g)));
            pp.set(ia,ja,dot(pi,pj));
          } else {
            pp.set(ia,ja,pp.get(ja,ia));
          }
        }
        pr.set(ia,0,dot(pi,r));
      }
      float[][] q0 = warp(u,convolve(na,ka,a,g));
      DMatrix qq = new DMatrix(nh,nh);
      DMatrix qr = new DMatrix(nh,1);
      for (int ih=0,ilag=kh; ih<nh; ++ih,++ilag) {
        float[][] qi = delay(ilag,q0);
        for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
          if (ih<=jh) {
            float[][] qj = delay(jlag,q0);
            qq.set(ih,jh,dot(qi,qj));
          } else {
            qq.set(ih,jh,qq.get(jh,ih));
          }
        }
        qr.set(ih,0,dot(qi,r));
      }
      DMatrix pq = new DMatrix(na,nh);
      DMatrix qp = new DMatrix(nh,na);
      for (int ia=0,ilag=ka; ia<na; ++ia,++ilag) {
        float[][] pi = convolve(nh,kh,h,warp(u,delay(ilag,g)));
        for (int jh=0,jlag=kh; jh<nh; ++jh,++jlag) {
          float[][] qj = delay(jlag,q0);
          double piqj = dot(pi,qj);
          pq.set(ia,jh,piqj);
          qp.set(jh,ia,piqj);
        }
      }

      // Composite matrix and right-hand-side.
      int nm = na+nh;
      DMatrix m = new DMatrix(nm,nm);
      DMatrix b = new DMatrix(nm,1);
      m.set( 0,na-1,0,na-1,pp); m.set( 0,na-1,na,nm-1,pq);
      m.set(na,nm-1,0,na-1,qp); m.set(na,nm-1,na,nm-1,qq);
      b.set( 0,na-1,0,0,pr);
      b.set(na,nm-1,0,0,qr);

      // Ensure that a = 1 for lag = 0.
      for (int im=0; im<nm; ++im)
        if (im!=-ka)
          m.set(im,-ka,0.0);
      for (int im=0; im<nm; ++im)
        if (im!=-ka)
          m.set(-ka,im,0.0);
      b.set(-ka,0,0.0);

      // Ensure positive-definite.
      for (int im=0; im<nm; ++im)
        m.set(im,im,m.get(im,im)*(1.0+_stabilityFactor));

      // Solve SPD system of equations.
      DMatrix x = m.chd().solve(b);
      float[] da = new float[na];
      float[] dh = new float[nh];
      for (int ia=0; ia<na; ++ia)
        da[ia] = (float)x.get(ia,0);
      for (int ih=0; ih<nh; ++ih)
        dh[ih] = (float)x.get(ih+na,0);

      // Update a and h after line search for good step size.
      float step = findStep(na,ka,a,da,nh,kh,h,dh,u,f,g);
      a = sub(a,mul(step,da));
      h = sub(h,mul(step,dh));

      // Update rms of residual r.
      r = sub(convolve(nh,kh,h,warp(u,convolve(na,ka,a,g))),f);
      rmsrOld = rmsr;
      rmsr = rms(r);
      trace("getAHGaussNewton: step="+step+" rmsr="+rmsr);
    }
    trace("getAHGaussNewton: end");
    return new float[][]{a,h};
  }
  private float findStep(
    final int na, final int ka, final float[] a, final float[] da,
    final int nh, final int kh, final float[] h, final float[] dh,
    final float[][] u, final float[][] f, final float[][] g)
  {
    BrentMinFinder bmf = new BrentMinFinder(
      new BrentMinFinder.Function() {
        public double evaluate(double step) {
          //trace("bmf: evaluate step="+step);
          float s = (float)step;
          float[] an = sub(a,mul(s,da));
          float[] hn = sub(h,mul(s,dh));
          float[][] rn = sub(convolve(nh,kh,hn,warp(u,convolve(na,ka,an,g))),f);
          return rms(rn);
        }
      });
    return (float)bmf.findMin(0.01,1.00,0.01);
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

  /**
   * Some useful array processing.
   */
  public double msq(float[] f) {
    return dot(f,f);
  }
  public double msq(float[][] f) {
    return dot(f,f);
  }
  public double rms(float[] f) {
    return sqrt(msq(f));
  }
  public double rms(float[][] f) {
    return sqrt(msq(f));
  }
  public float[] msqNormalize(float[] f) {
    return mul(1.0f/(float)msq(f),f);
  }
  public float[][] msqNormalize(float[][] f) {
    return mul(1.0f/(float)msq(f),f);
  }
  public float[] rmsNormalize(float[] f) {
    return mul(1.0f/(float)rms(f),f);
  }
  public float[][] rmsNormalize(float[][] f) {
    return mul(1.0f/(float)rms(f),f);
  }

  private static float[] solveChd(DMatrix q, DMatrix r) {
    int nx = r.getM();
    DMatrixChd chd = new DMatrixChd(q);
    DMatrix v = chd.solve(r);
    float[] x = new float[nx];
    for (int ix=0; ix<nx; ++ix)
      x[ix] = (float)v.get(ix,0);
    return x;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final SincInterpolator _si =
    SincInterpolator.fromErrorAndFrequency(0.01,0.40);
  private static final WarpingFilter _wf = new WarpingFilter();
  static {
    _wf.setOversamplingLimit(4.0);
  }

  private int _itmin = -1;
  private int _itmax = -1;
  private float _stabilityFactor = 0.0f;
  private int _maxiter = 0;

  // Dot products normalized by the number of products summed.
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

  private static void trace(String s) {
    System.out.println(s);
  }
}
