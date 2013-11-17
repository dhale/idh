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
   * @param smax maximum NMO stretch.
   * @param vnmo NMO velocity.
   */
  public WaveletNmo(Sampling st, Sampling sx, double smax, double vnmo) {
    this(st,sx,smax,fillfloat((float)vnmo,st.getCount()));
  }

  /**
   * Constructs a wavelet estimator for specified sampling and NMO velocities.
   * @param st sampling of time.
   * @param sx sampling of offset.
   * @param smax maximum NMO stretch.
   * @param vnmo array[nt] of NMO velocities.
   */
  public WaveletNmo(Sampling st, Sampling sx, double smax, float[] vnmo) {
    _st = st;
    _sx = sx;
    _nmo = new Nmo(st,sx,smax,vnmo);
  }

  /**
   * Sets the min-max range of frequencies for wavelet.
   * @param fmin minimum frequency.
   * @param fmax maximum frequency.
   */
  public void setFrequencyRange(double fmin, double fmax) {
    double dt = _st.getDelta();
    _bpf = new BandPassFilter(fmin*dt,fmax*dt,0.05*dt,0.01);
  }

  /**
   * Sets the stability factor by which to scale zero-lag correlation.
   * @param sfac stability factor.
   */
  public void setStabilityFactor(double sfac) {
    _sfac = sfac;
  }

  public float getVariancePef(int na, int ka, float[] a, float[][] f) {
    float[][] g = applyFilter(na,ka,a,f);
    return pow(rms(g),2.0f);
  }
  public float getVariance(int na, int ka, float[] a, float[][] f) {
    float[][] g = applyBNmoA(na,ka,a,f);
    float[][] r = stackAndReplicate(_nmo,g);
    return pow(rms(sub(g,r)),2.0f);
  }
  public float getNormalizedVariance(int na, int ka, float[] a, float[][] f) {
    float[][] g = applyBNmoA(na,ka,a,f);
    float[][] r = stackAndReplicate(_nmo,g);
    return pow(rms(sub(g,r))/rms(g),2.0f);
  }

  /**
   * Returns inverse wavelet a estimated via PEF of CMP gather.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param f array[nx][nt] for CMP gather.
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getInverseAPef(int na, int ka, float[][] f) {
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
      c.set(ic,ic,c.get(ic,ic)*_sfac);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

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
   * @param f array[nx][nt] for CMP gather.
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getInverseA(int na, int ka, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Differences d for all lags of inverse wavelet a.
    float[][][] d = computeDifferences(na,ka,_nmo,_bpf,_st,_sx,f);

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
      c.set(ic,ic,c.get(ic,ic)*_sfac);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

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
    float[][] af = applyFilter(na,ka,a,f);
    float[][] saf = _nmo.apply(af);
    return applyFilter(nh,kh,h,saf);
  }

  public float[][] applyBNmoA(int na, int ka, float[] a, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    float[][] af = applyFilter(na,ka,a,f);
    float[][] saf = _nmo.apply(af);
    float[][] g = saf;
    if (_bpf!=null) {
      for (int ix=0; ix<nx; ++ix)
        _bpf.apply(saf[ix],g[ix]);
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _st;
  private Sampling _sx;
  private Nmo _nmo;
  private BandPassFilter _bpf;
  private double _sfac = 1.0;

  private static class Nmo {
    public Nmo(Sampling st, Sampling sx, double smax, double vnmo) {
      this(st,sx,smax,fillfloat((float)vnmo,st.getCount()));
    }
    public Nmo(Sampling st, Sampling sx, double smax, float[] vnmo) {
      int nt = st.getCount();
      double dt = st.getDelta();
      double ft = st.getFirst();
      int nx = sx.getCount();
      float amin = 1.0f/(float)smax;
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
        for (int it=1; it<nt; ++it) {
          _a[ix][it] = (_t[ix][it]-_t[ix][it-1])*odt;
          if (_a[ix][it]<amin)
            _a[ix][it] = 0.0f;
        }
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
   * For each lag of the inverse wavelet, computes differences between
   * NMO-corrected gathers and the stacked-and-replicated versions of those
   * gathers.
   */
  private static float[][][] computeDifferences(
    int na, int ka, Nmo nmo, BandPassFilter bpf,
    Sampling st, Sampling sx, float[][] f)
  {
    final float smax = 3.0f;
    final float wmin = 1.0f/smax;
    int nt = f[0].length;
    int nx = f.length;
    float[][] w = nmo.getAmplitudes();
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        w[ix][it] = 1.0f/max(wmin,w[ix][it])-1.0f;
        //w[ix][it] = 1.0f;
      }
    }
    float[][][] d = new float[na][nx][nt];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[][] df = delay(lag,f);
      float[][] dg = nmo.apply(df);
      if (bpf!=null) {
        for (int ix=0; ix<nx; ++ix)
          bpf.apply(dg[ix],dg[ix]);
      }
      float[][] rdg = stackAndReplicate(nmo,dg);
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
   * Returns true if array contains non-zero samples.
   */
  private static boolean nonZero(float[] f) {
    int n = f.length;
    boolean nonzero = false;
    for (int i=0; i<n && !nonzero; ++i) {
      if (f[i]!=0.0f)
        nonzero = true;
    }
    return nonzero;
  }

  /**
   * Returns the number of non-zero traces in gather.
   */
  private static int countNonZero(float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    int nnz = 0;
    for (int ix=0; ix<nx; ++ix) {
      if (nonZero(f[ix]))
        ++nnz;
    }
    return nnz;
  }

  /**
   * Stacks the input CMP gather f, and then replicates the stack into a
   * returned CMP gather g. If any traces in f contain only zeros, the
   * corresponding traces in g will contain only zeros.
   */
  private static float[][] stackAndReplicate(Nmo nmo, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    float[] s = new float[nt];
    for (int ix=0; ix<nx; ++ix)
      add(f[ix],s,s);
    float[] ss = getStackScaling(nmo,f);
    boolean[][] nz = getNonZeroMask(nmo,f);
    float[][] g = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      mul(ss,s,g[ix]);
      for (int it=0; it<nt; ++it) {
        if (!nz[ix][it])
          g[ix][it] = 0.0f;
      }
    }
    return g;
  }
  private static float[] getStackScaling(Nmo nmo, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    boolean[][] mask = getNonZeroMask(nmo,f);
    float[] s = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        if (mask[ix][it])
          s[it] += 1.0f;
      }
    }
    for (int it=0; it<nt; ++it)
      s[it] = 1.0f/max(s[it],1.0f);
    return s;
  }
  private static boolean[][] getNonZeroMask(Nmo nmo, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    float[][] a = nmo.getAmplitudes();
    boolean[][] mask = new boolean[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      boolean nz = nonZero(f[ix]);
      if (nz) {
        for (int it=0; it<nt; ++it) {
          mask[ix][it] = a[ix][it]>0.0f;
        }
      }
    }
    return mask;
  }

  private static float[][] applyFilter(int nh, int kh, float[] h, float[][] f)
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][] g = new float[nx][nt];
    applyFilter(nh,kh,h,f,g);
    return g;
  }

  private static void applyFilter(
    int nh, int kh, float[] h, float[][] f,  float[][] g)
  {
    int nt = f[0].length;
    int nx = f.length;
    for (int ix=0; ix<nx; ++ix)
      conv(nh,kh,h,nt,0,f[ix],nt,0,g[ix]);
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

  private static float rms(float[][] f) {
    int nt = f[0].length;
    int nxnz = countNonZero(f);
    return (float)(sqrt(dot(f,f)/nxnz/nt));
  }
}
