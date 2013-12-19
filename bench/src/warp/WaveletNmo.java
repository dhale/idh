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
 * @version 2013.12.18
 */
public class WaveletNmo {

  /**
   * Constructs a wavelet estimator for specified NMO.
   * @param nmo the NMO correction.
   */
  public WaveletNmo(NormalMoveout nmo) {
    _nmo = nmo;
  }

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
   * Sets the min-max range of frequencies in wavelet.
   * @param fmin minimum frequency, in cycles/sample.
   * @param fmax maximum frequency, in cycles/sample.
   */
  public void setFrequencyRange(double fmin, double fmax) {
    _bpf = new BandPassFilter(fmin,fmax,0.05,0.01);
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

  public float getVariancePef(
    int na, int ka, float[] a, float[][] f) 
  {
    float[][] g = applyFilter(na,ka,a,f);
    return pow(rms(g),2.0f);
  }
  public float getVarianceNmo(
    int na, int ka, float[] a, 
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    float[][] g = applyBNmoA(na,ka,a,st,sx,vnmo,f);
    float[][] r = _nmo.stackAndReplicate(g);
    return pow(rms(sub(g,r)),2.0f);
  }
  public float getNormalizedVarianceNmo(
    int na, int ka, float[] a, 
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    float[][] g = applyBNmoA(na,ka,a,st,sx,vnmo,f);
    float[][] r = _nmo.stackAndReplicate(g);
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
      /* band-pass filter causes unstable estimates of wavelet
      if (_bpf!=null) {
        for (int ix=0; ix<nx; ++ix)
          _bpf.apply(d[ia][ix],d[ia][ix]);
      }
      */
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
   * @param st time sampling.
   * @param sx offset sampling.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] for CMP gather.
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getInverseANmo(
    int na, int ka, 
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    int nt = f[0].length;
    int nx = f.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Differences d for all lags of inverse wavelet a.
    float[][][] d = computeDifferences(na,ka,_nmo,_bpf,st,sx,vnmo,f);

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
   * @param st time sampling.
   * @param sx offset sampling.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] with input CMP gather.
   * @return array[nx][nt] with output CMP gather.
   * 
   */
  public float[][] applyNmo(
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    return _nmo.apply(st,sx,vnmo,f);
  }

  /**
   * Returns differences between NMO-corrected gathers and stacks.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param st time sampling.
   * @param sx offset sampling.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] with input CMP gather.
   * @return array[na][nx][nt] of difference gathers.
   */
  public float[][][] getDifferenceGathers(
    int na, int ka, 
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    float[][][] d = computeDifferences(na,ka,_nmo,_bpf,st,sx,vnmo,f);
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
   * @param st time sampling.
   * @param sx offset sampling.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] with input CMP gather.
   * @return array[nx][nt] with output CMP gather.
   */
  public float[][] applyHNmoA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][] af = applyFilter(na,ka,a,f);
    float[][] saf = _nmo.apply(st,sx,vnmo,af);
    return applyFilter(nh,kh,h,saf);
  }

  public float[][] applyBNmoA(
    int na, int ka, float[] a, 
    Sampling st, Sampling sx, float[] vnmo, float[][] f) 
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][] af = applyFilter(na,ka,a,f);
    float[][] saf = _nmo.apply(st,sx,vnmo,af);
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
  private NormalMoveout _nmo;
  private BandPassFilter _bpf;
  private double _sfac = 1.0;
  private int _itmin,_itmax;

  /**
   * For each lag of the inverse wavelet, computes differences between
   * NMO-corrected gathers and the stacked-and-replicated versions of those
   * gathers.
   */
  private static float[][][] computeDifferences(
    int na, int ka, NormalMoveout nmo, BandPassFilter bpf,
    Sampling st, Sampling sx, float[] vnmo, float[][] f)
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][][] ta = nmo.timesAndAmplitudes(st,sx,vnmo,f);
    float[][] t = ta[0];
    float[][] a = ta[1];
    float[][][] d = new float[na][nx][nt];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[][] df = delay(lag,f);
      float[][] dg = nmo.apply(st,t,a,df);
      if (bpf!=null) {
        for (int ix=0; ix<nx; ++ix) {
          bpf.apply(dg[ix],dg[ix]);
          for (int it=0; it<nt; ++it) {
            if (a[ix][it]==0.0f)
              dg[ix][it] = 0.0f;
          }
        }
      }
      float[][] rdg = nmo.stackAndReplicate(dg);
      for (int ix=0; ix<nx; ++ix) {
        for (int it=0; it<nt; ++it) {
          d[ia][ix][it] = dg[ix][it]-rdg[ix][it];
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

  private static float[][] applyFilter(
    int nh, int kh, float[] h, float[][] f)
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
    preserveLeadingZeros(f,g);
  }

  private static void preserveLeadingZeros(float[][] f, float[][] g) {
    int nx = f.length;
    for (int ix=0; ix<nx; ++ix)
      preserveLeadingZeros(f[ix],g[ix]);
  }

  private static void preserveLeadingZeros(float[] f, float[] g) {
    int nt = f.length;
    int nz = countLeadingZeros(f);
    for (int it=0; it<nz; ++it)
      g[it] = 0.0f;
  }

  private static int countLeadingZeros(float[] f) {
    int n = f.length;
    int nz = 0;
    for (int i=0; i<n && f[i]==0.0f; ++i)
      ++nz;
    return nz;
  }

  private double dot(float[][] f, float[][] g) {
    int nt = f[0].length;
    int nx = f.length;
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix) 
      for (int it=_itmin; it<=_itmax; ++it) 
        sum += f[ix][it]*g[ix][it];
    return sum;
  }

  private float rms(float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    return (float)(sqrt(dot(f,f)/nx/nt));
  }
}
