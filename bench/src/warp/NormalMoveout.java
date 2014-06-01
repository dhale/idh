/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Normal-moveout (NMO) correction.
 * For a gather of input seismic traces f(t,x), NMO correction is defined by
 * the transformation g(u,x) = a(u,x)*f(t(u,x),x), where t(u,x) are times at
 * which to evaluate f(t,x), a(u,x) are amplitude scale factors, and g(u,x) is
 * an NMO-corrected gather of output traces.
 * <p>
 * The times t(u,x) may be specified directly, or computed from specified NMO
 * velocities for hyperbolic moveout. Amplitudes a(u,x) may also be specified,
 * or computed to preserve mutes and/or limit the maximum NMO stretch. Note
 * that times t(u,x), amplitudes a(u,x), and the output gather g(u,x) are
 * functions of output time u, whereas f(t,x) is a function of input time t.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.12.17
 */
public class NormalMoveout {

  /**
   * Sets the maximum NMO stretch factor.
   * Used to zero amplitudes for which NMO stretch is excessive.
   * @param smax the maximum stretch factor.
   */
  public void setStretchMax(double smax) {
    _smax = (float)smax;
  }

  /**
   * Returns arrays of times and amplitudes for NMO correction.
   * Sets to zero any amplitudes corresponding to (1) leading zeros in the
   * input gather or (2) samples for which NMO stretch is excessive. For all
   * other samples, non-zero amplitudes are simply the inverse of NMO stretch.
   * @param st uniform time sampling.
   * @param sx offset sampling; need not be uniform.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] for input gather f(t,x).
   * @return array {t,a} of times and amplitudes.
   */
  public float[][][] timesAndAmplitudes(
    Sampling st, Sampling sx, float[] vnmo, float[][] f)
  {
    float[][] t = getTimes(st,sx,vnmo);
    float[][] a = getAmplitudes(st,_smax,f,t);
    return new float[][][]{t,a};
  }

  /**
   * Applies this correction for specified times and amplitudes.
   * @param st uniform time sampling.
   * @param t array[nx][nt] of times t(u,x).
   * @param a array[nx][nt] of amplitudes a(u,x).
   * @param f array[nx][nt] for input gather f(t,x).
   * @return array[nx][nt] for output NMO-corrected gather g(u,x).
   */
  public float[][] apply(Sampling st, float[][] t, float[][] a, float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    double dt = st.getDelta();
    double ft = st.getFirst();
    float[][] g = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      _si.interpolate(nt,dt,ft,f[ix],nt,t[ix],g[ix]);
      for (int it=0; it<nt; ++it) {
        g[ix][it] *= a[ix][it];
      }
    }
    return g;
  }

  /**
   * Applies this correction for specified NMO velocities.
   * @param st uniform time sampling.
   * @param sx offset sampling; need not be uniform.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] for input gather.
   * @return array[nx][nt] NMO-corrected output gather.
   */
  public float[][] apply(Sampling st, Sampling sx, float[] vnmo, float[][] f) {
    float[][][] ta = timesAndAmplitudes(st,sx,vnmo,f);
    float[][] t = ta[0];
    float[][] a = ta[1];
    return apply(st,t,a,f);
  }

  /**
   * Applies this correction for specified constant NMO velocity.
   * @param st uniform time sampling.
   * @param sx offset sampling; need not be uniform.
   * @param vnmo NMO velocity.
   * @param f array[nx][nt] for input gather.
   * @return array[nx][nt] NMO-corrected output gather.
   */
  public float[][] apply(Sampling st, Sampling sx, double vnmo, float[][] f) {
    return apply(st,sx,fillfloat((float)vnmo,st.getCount()),f);
  }

  /**
   * For each offset, counts the number of leading zeros in a trace.
   * @param f array[nx][nt] in which to count leading zeros.
   * @return array[nx] counts of leading zeros, one for each trace.
   */
  public static int[] countLeadingZeros(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    int[] ifnz = new int[nx];
    for (int ix=0; ix<nx; ++ix)
      ifnz[ix] = countLeadingZeros(f[ix]);
    return ifnz;
  }

  /**
   * For each time sample, counts the number of non-zero values.
   * Omits only leading zeros from the counts. In other words, for every
   * offset, any zero values that occur after the first non-zero value are
   * included in the returned counts. These counts are typically used to
   * normalize a stack over offsets, after NMO correction.
   * @param f array[nx][nt] in which to count non-zero values.
   * @return array[nt] of counts, one for each time sample.
   */
  public static int[] countNonZero(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    int[] nnz = new int[nt];
    for (int ix=0; ix<nx; ++ix) {
      int nz = countLeadingZeros(f[ix]);
      for (int it=nz; it<nt; ++it)
        nnz[it] += 1.0f;
    }
    return nnz;
  }

  /**
   * Stacks the specified gather over offset.
   * For each time sample, normalizes the stack by the number of values
   * summed, not counting any leading zeros.
   * @param f array[nx][nt] for the gather to be stacked.
   * @return array[nt] for the stack.
   */
  public float[] stack(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    float[] s = new float[nt];
    float[] c = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      int nz = countLeadingZeros(f[ix]);
      for (int it=nz; it<nt; ++it) {
        c[it] += 1.0f;
        s[it] += f[ix][it];
      }
    }
    for (int it=0; it<nt; ++it)
      s[it] /= max(c[it],1.0f);
    return s;
  }

  /**
   * Stacks and replicates the specified gather over offset.
   * For each time sample, normalizes the stack by the number of values
   * summed, not counting any leading zeros. Leading zeros in the input
   * will be zero in the output gather.
   * @param f array[nx][nt] for the gather to be stacked and replicated.
   * @return array[nt] for the stacked and replicated gather.
   */
  public float[][] stackAndReplicate(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    float[] s = new float[nt];
    float[] c = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      int nz = countLeadingZeros(f[ix]);
      for (int it=nz; it<nt; ++it) {
        c[it] += 1.0f;
        s[it] += f[ix][it];
      }
    }
    for (int it=0; it<nt; ++it)
      s[it] /= max(c[it],1.0f);
    float[][] g = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      int nz = countLeadingZeros(f[ix]);
      for (int it=nz; it<nt; ++it)
        g[ix][it] = s[it];
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private SincInterpolator
    _si = SincInterpolator.fromErrorAndFrequency(0.01,0.40);
  private float _smax = 0.1f*Float.MAX_VALUE;

  private static int countLeadingZeros(float[] f) {
    int n = f.length;
    int nz = 0;
    for (int i=0; i<n && f[i]==0.0f; ++i)
      ++nz;
    return nz;
  }

  private static float[][] getTimes(Sampling st, Sampling sx, float[] vnmo) {
    int nt = st.getCount();
    int nx = sx.getCount();
    float[][] t = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
      for (int it=0; it<nt; ++it) {
        double tnmo = st.getValue(it);
        t[ix][it] = (float)sqrt(tnmo*tnmo+(x*x)/(vnmo[it]*vnmo[it]));
      }
    }
    return t;
  }

  private static float[][] getAmplitudes(
    Sampling st, float smax, float[][] f, float[][] t) 
  {
    int nx = f.length;
    int nt = f[0].length;
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
    float odt = 1.0f/dt;
    float dtmin = dt/smax;
    float[][] a = new float[nx][nt];

    // For all offsets, ...
    for (int ix=0; ix<nx; ++ix) {

      // Time of first non-zero input sample.
      int nz = countLeadingZeros(f[ix]);
      float tnz = ft+nz*dt;

      // Number of leading zeros in output. A leading output sample is zero 
      // if either (1) the corresponding input samples and all prior input 
      // samples are zero, or (2) NMO stretch would exceed the maximum.
      nz = 0;
      if (t[ix][0]<tnz || t[ix][1]-t[ix][0]<dtmin)
        ++nz;
      for (int it=1; it<nt; ++it) {
        if (t[ix][it]<tnz || t[ix][it]-t[ix][it-1]<dtmin)
          ++nz;
      }

      // Compute only the non-zero amplitudes. These amplitudes are simply the
      // inverse of NMO stretch.
      if (nz==0) {
        a[ix][0] = (t[ix][1]-t[ix][0])*odt;
        ++nz;
      }
      for (int it=nz; it<nt; ++it)
        a[ix][it] = (t[ix][it]-t[ix][it-1])*odt;
    }
    return a;
  }
}
