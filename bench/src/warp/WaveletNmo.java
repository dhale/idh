/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
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
 * Estimates a wavelet via NMO correction of a CMP gather.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.27
 */
public class WaveletNmo {

  /**
   * Constructs a wavelet estimator for specified sampling and NMO velocity.
   * @param st sampling of time.
   * @param sx sampling of offset.
   * @param vnmo NMO velocity.
   */
  public WaveletNmo(Sampling st, Sampling sx, double vnmo) {
    this(st,sx,fillfloat((float)vnmo,st.getCount()));
  }

  /**
   * Constructs a wavelet estimator for specified sampling and NMO velocities.
   * @param st sampling of time.
   * @param sx sampling of offset.
   * @param vnmo array[nt] of NMO velocities.
   */
  public WaveletNmo(Sampling st, Sampling sx, float[] vnmo) {
    _st = st;
    _sx = sx;
    _nmo = new Nmo(st,sx,vnmo);
  }

  /**
   * Sets the min-max range of frequencies for wavelet.
   * @param fmin minimum frequency.
   * @param fmax maximum frequency.
   */
  public void setFrequencyRange(double fmin, double fmax) {
    double dt = _st.getDelta();
    _bpf = new BandPassFilter(fmin*dt,fmax*dt,0.1*dt,0.01);
  }

  /**
   * Returns inverse wavelet a estimated via spiking decon of CMP gather.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param f array[nx][nt] for CMP gather.
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getSpikingDeconA(int na, int ka, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;

    // CMP gather for different time shifts (band-pass filtered?).
    float[][][] d = new float[na][nx][nt];
    for (int ia=0; ia<na; ++ia) {
      d[ia] = delay(ka+ia,f);
      //for (int ix=0; ix<nx; ++ix)
      //  _bpf.apply(d[ia][ix],d[ia][ix]);
    }

    // The matrix C and right-hand-side vector b, for Ca = b. For zero lag, we
    // have a0 = a[-ka] = 1, so that only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*1.0000);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    System.out.println("c=\n"+c);
    System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }

  /**
   * Returns inverse wavelet a estimated via NMO correction of CMP gather.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param st sampling of time.
   * @param sx sampling of offset.
   * @param vnmo array[nt] of nmo velocities.
   * @param f array[nx][nt] for CMP gather.
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getInverseA(int na, int ka, float[][] f) {
    return estimateInverseWavelet(na,ka,_nmo,_bpf,_st,_sx,f);
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
    float[] caa = new float[nh];
    float[] ca1 = new float[nh];
    xcor(na,ka,a,na,ka,a,nh, 0,caa); //caa[0] *= 1.0001;
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }

  /**
   * Applies nmo(f).
   * @param f array[nx][nt] with input CMP gather.
   * @return array[nx][nt] with output CMP gather.
   * 
   */
  public float[][] applyNmo(float[][] f) {
    return _nmo.apply(f);
  }

  /**
   * Returns differences between NMO-corrected gathers and stacks.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param f array[nx][nt] with input CMP gather.
   * @return array[na][nx][nt] of difference gathers.
   */
  public float[][][] getDifferenceGathers(int na,int ka, float[][] f) {
    float[][][] d = computeDifferences(na,ka,_nmo,_bpf,_st,_sx,f);
    for (int ia=0,lag=ka; ia<na; ++ia,++lag)
      d[ia] = delay(-lag,d[ia]);
    //for (int ia=1; ia<na; ++ia)
    //  d[ia] = sub(d[ia],d[0]);
    return d;
  }

  /**
   * Applies the composite h*nmo(a*f), where * denotes convolution.
   * The sequence of operations is (1) convolution with the inverse wavelet a,
   * (2) NMO correction, and (3) convolution with the wavelet h.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param f array[nx][nt] with input CMP gather.
   * @return array[nx][nt] with output CMP gather.
   */
  public float[][] applyHNmoA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[][] f)
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][] g = new float[nx][nt];
    float[][] af = g;
    for (int ix=0; ix<nx; ++ix)
      conv(na,ka,a,nt,0,f[ix],nt,0,af[ix]);
    float[][] saf = _nmo.apply(af);
    for (int ix=0; ix<nx; ++ix)
      conv(nh,kh,h,nt,0,saf[ix],nt,0,g[ix]);
    return g;
  }

  public float[][] applyBNmoA(int na, int ka, float[] a, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    float[][] g = new float[nx][nt];
    float[][] af = g;
    for (int ix=0; ix<nx; ++ix)
      conv(na,ka,a,nt,0,f[ix],nt,0,af[ix]);
    float[][] saf = _nmo.apply(af);
    for (int ix=0; ix<nx; ++ix)
      _bpf.apply(saf[ix],g[ix]);
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _st;
  private Sampling _sx;
  private Nmo _nmo;
  private BandPassFilter _bpf;

  private static class Nmo {
    public Nmo(Sampling st, Sampling sx, double vnmo) {
      this(st,sx,fillfloat((float)vnmo,st.getCount()));
    }
    public Nmo(Sampling st, Sampling sx, float[] vnmo) {
      int nt = st.getCount();
      double dt = st.getDelta();
      double ft = st.getFirst();
      int nx = sx.getCount();
      _st = st;
      _sx = sx;
      _si = SincInterp.fromErrorAndFrequency(0.01,0.45);
      System.out.println("si: max length = "+_si.getMaximumLength());
      _t = new float[nx][nt];
      _a = new float[nx][nt];
      for (int ix=0; ix<nx; ++ix) {
        double x = sx.getValue(ix);
        for (int it=0; it<nt; ++it) {
          double ti = ft+it*dt;
          _t[ix][it] = (float)sqrt(ti*ti+(x*x)/(vnmo[it]*vnmo[it]));
        }
        float odt = 1.0f/(float)dt;
        _a[ix][0] = (_t[ix][1]-_t[ix][0])*odt;
        for (int it=1; it<nt; ++it)
          _a[ix][it] = (_t[ix][it]-_t[ix][it-1])*odt;
      }
    }
    public float[][] apply(float[][] f) {
      int nt = _st.getCount();
      double dt = _st.getDelta();
      double ft = _st.getFirst();
      int nx = _sx.getCount();
      float[][] g = new float[nx][nt];
      for (int ix=0; ix<nx; ++ix) {
        _si.interpolate(nt,dt,ft,f[ix],nt,_t[ix],g[ix]);
        for (int it=0; it<nt; ++it) {
          g[ix][it] *= _a[ix][it];
        }
        //for (int it=nt-200; it<nt; ++it)
        //  g[ix][it] *= 0.5f*(1.0f+cos(FLT_PI*0.005f*(it-nt+200)));
      }
      return g;
    }
    public float[][] getTimes() {
      return copy(_t);
    }
    public float[][] getAmplitudes() {
      return copy(_a);
    }
    private Sampling _st,_sx;
    private SincInterp _si;
    private float[][] _t;
    private float[][] _a;
  }

  /**
   * Estimates the inverse a of the wavelet contained in CMP gather f.
   */
  private static float[] estimateInverseWavelet(
    int na, int ka, Nmo nmo, BandPassFilter bpf,
    Sampling st, Sampling sx, float[][] f) 
  {
    int nt = f[0].length;
    int nx = f.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Differences d for all lags of inverse wavelet a.
    float[][][] d = computeDifferences(na,ka,nmo,bpf,st,sx,f);

    // The matrix C and right-hand-side vector b, for Ca = b. For zero lag, we
    // have a0 = a[-ka] = 1, so that only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*1.0000);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    System.out.println("c=\n"+c);
    System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }

  /**
   * For each lag of the inverse wavelet, computes differences between
   * NMO-corrected gathers and the stacked-and-replicated versions of those
   * gathers.
   */
  private static float[][][] computeDifferences(
    int na, int ka, Nmo nmo, BandPassFilter bpf,
    Sampling st, Sampling sx, float[][] f)
  {
    final float smax = 3.0f;
    final float wmin = 1.0f/(smax-1.0f);
    int nt = f[0].length;
    int nx = f.length;
    float[][] w = nmo.getAmplitudes();
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        w[ix][it] = (w[ix][it]>wmin)?1.0f/max(wmin,w[ix][it])-1.0f:0.0f;
        //w[ix][it] = (ix==0 || ix==nx-1)?1.0f:0.0f;
        //w[ix][it] = 1.0f;
      }
    }
    float[][][] d = new float[na][nx][nt];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[][] df = delay(lag,f);
      float[][] dg = nmo.apply(df);
      for (int ix=0; ix<nx; ++ix)
        bpf.apply(dg[ix],dg[ix]);
      float[][] rdg = stackAndReplicate(dg);
      for (int ix=0; ix<nx; ++ix) {
        for (int it=0; it<nt; ++it) {
          d[ia][ix][it] = w[ix][it]*(dg[ix][it]-rdg[ix][it]);
        }
      }
    }
    return d;
  }

  /**
   * Delays the CMP gather f by specified lag (which may be negative).
   */
  private static float[][] delay(int lag, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    int itlo = max(0,lag);   // 0 <= it-lag
    int ithi = min(nt,nt+lag); // it-lag < nt
    float[][] g = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<itlo; ++it)
        g[ix][it] = 0.0f;
      for (int it=itlo; it<ithi; ++it)
        g[ix][it] = f[ix][it-lag];
      for (int it=ithi; it<nt; ++it)
        g[ix][it] = 0.0f;
    }
    return g;
  }

  /**
   * Stacks the input CMP gather f, and then replicates the stack into a
   * returned CMP gather g. If any traces in f contain only zeros, the
   * corresponding traces in g will contain only zeros.
   */
  private static float[][] stackAndReplicate(float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    int nnz = 0;
    boolean[] nz = new boolean[nx];
    float[][] g = new float[nx][nt];
    float[] s = g[0];

    // For all traces in input gather, ...
    for (int ix=0; ix<nx; ++ix) {
      
      // Look for at least one non-zero sample.
      boolean nonzero = false;
      for (int it=0; it<nt && !nonzero; ++it) {
        if (f[ix][it]!=0.0f)
          nonzero = true;
      }

      // If at least one sample is non-zero, accumulate.
      if (nonzero) {
        ++nnz;
        nz[ix] = true;
        for (int it=0; it<nt; ++it)
          s[it] += f[ix][it];
      }
    }

    // Scale stack by number of non-zero samples.
    float scale = 1.0f/nnz;
    for (int it=0; it<nt; ++it)
      s[it] *= scale;

    // RESEARCH (use zero-offset trace from input gather)!
    //s = f[0];

    // Replicate stack into output gather.
    for (int ix=0; ix<nx; ++ix) {
      if (nz[ix]) {
        for (int it=0; it<nt; ++it) {
          g[ix][it] = s[it];
        }
      }
    }
    return g;
  }

  private static double dot(float[][] f, float[][] g) {
    int nt = f[0].length;
    int nx = f.length;
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix) 
      for (int it=0; it<nt; ++it) 
        sum += f[ix][it]*g[ix][it];
    return sum;
  }
}
