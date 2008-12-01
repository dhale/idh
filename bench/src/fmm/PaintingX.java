/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

// for testing
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;

/**
 * A 2D array of painted values, where most values are painted automatically.
 * Except for a relatively small number of fixed samples painted explicitly, 
 * most samples are painted by extrapolation and interpolation guided by
 * structure tensors. Intuitively, paint flows slowly through locations with
 * relatively high structure. Structure tensors may also cause paint to
 * flow anisotropically, in directions corresponding to relatively low
 * structure.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.12
 */
public class PaintingX {

  /**
   * An interface for classes of structure tensors. Each tensor is a
   * symmetric positive-definite 2-by-2 matrix {{s11,s12},{s12,s22}}.
   */
  public interface Tensors {

    /**
     * Gets structure tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param s array {s11,s12,s22} of tensor elements.
     */
    public void getTensor(int i1, int i2, float[] s);
  }

  /**
   * Constructs a painting with constant identity structure tensors.
   * In this case, time = distance, which is useful for testing.
   * Painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param nv number of values painted for each sample
   */
  public PaintingX(int n1, int n2, int nv) {
    this(n1,n2,nv,new IdentityTensors());
  }
  
  /**
   * Constructs a painting for the specified structure tensor field.
   * All painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param nv number of values painted for each sample
   * @param st structure tensors.
   */
  public PaintingX(int n1, int n2, int nv, Tensors st) {
    _n1 = n1;
    _n2 = n2;
    _nv = nv;
    _st = st;
    _dv = new float[nv];
    _k1 = new int[n2][n1];
    _k2 = new int[n2][n1];
    _tk = new float[n2][n1];
    _vk = new float[n2][n1][];
    _type = new byte[n2][n1];
    _mark = new int[n2][n1];
    _hmin = new TimeHeap2(TimeHeap2.Type.MIN,n1,n2);
    _hmax = new TimeHeap2(TimeHeap2.Type.MAX,n1,n2);
    clearAll();
  }

  /**
   * Sets the structure tensors used in this painting.
   * @param st structure tensors.
   */
  public void setTensors(Tensors st) {
    _st = st;
  }

  /**
   * Sets the default value for value index zero.
   * Default values are used for all clear (not painted) samples.
   * @param value the default value.
   */
  public void setDefaultValue(float value) {
    _dv[0] = value;
  }

  /**
   * Sets all default values.
   * Default values are used for all clear (not painted) samples.
   * @param value array[nv] of default values.
   */
  public void setDefaultValues(float[] values) {
    for (int kv=0; kv<_nv; ++kv)
      _dv[kv] = values[kv];
  }

  /**
   * Clears all painted values, including all fixed values.
   */
  public void clearAll() {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        clear(i1,i2);
      }
    }
  }

  /**
   * Clears painted values that are not fixed.
   */
  public void clearNotFixed() {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_type[i2][i1]!=FIXED) {
          clear(i1,i2);
        }
      }
    }
  }

  /**
   * Erases values for any fixed sample with specified indices.
   * If the specified sample is not fixed, this method does nothing.
   * @param i1 index in 1st dimension of sample to erase.
   * @param i2 index in 2nd dimension of sample to erase.
   */
  public void eraseFixedAt(int i1, int i2) {
    if (_type[i2][i1]==FIXED)
      clear(i1,i2);
  }
  private void clear(int i1, int i2) {
    _type[i2][i1] = CLEAR;
    _mark[i2][i1] = _known;
    _k1[i2][i1] = -1;
    _k2[i2][i1] = -1;
    _tk[i2][i1] = TIME_INVALID;
    _vk[i2][i1] = null;
  }

  /**
   * Paints the specified sample with one specified value at index zero.
   * Paints default values for indices greater than zero.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param v value at index zero for the painted sample.
   */
  public void paintAt(int i1, int i2, float v) {
    _type[i2][i1] = FIXED;
    _k1[i2][i1] = i1;
    _k2[i2][i1] = i2;
    _tk[i2][i1] = TIME_INVALID;
    _vk[i2][i1] = new float[_nv];
    _vk[i2][i1][0] = v;
    for (int iv=1; iv<_nv; ++iv)
      _vk[i2][i1][iv] = _dv[iv];
  }

  /**
   * Paints the specified sample with specified values.
   * After painting, the specified sample is fixed; its values will
   * not change in any subsequent extrapolation or interpolation.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param v array of values for painted sample; by copy, not by reference.
   */
  public void paintAt(int i1, int i2, float[] v) {
    _type[i2][i1] = FIXED;
    _k1[i2][i1] = i1;
    _k2[i2][i1] = i2;
    _tk[i2][i1] = TIME_INVALID;
    _vk[i2][i1] = Array.copy(v);
  }

  /**
   * Extrapolates values from all fixed (explicitly painted) samples.
   * After extrapolation, all samples are either fixed or extrapolated.
   */
  public void extrapolate() {

    // Clear all samples that are not fixed, and insert all fixed 
    // samples into the max-heap with huge (invalid) times that can
    // only get smaller. After the max-heap is built, any one of the
    // fixed samples could be at the top of the heap.
    _hmax.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        _tk[i2][i1] = TIME_INVALID;
        if (_type[i2][i1]==FIXED) {
          _hmax.insert(i1,i2,TIME_INVALID);
        } else {
          _type[i2][i1] = CLEAR;
        }
      }
    }

    // Extrapolate from all fixed samples, one at a time, in order of
    // decreasing time. The fixed sample with the largest time is at the
    // top of the max-heap. As we extrapolate from this fixed sample, times 
    // for other fixed samples that remain in the heap may be reduced, so 
    // that their order may change. We choose a decreasing order to reduce 
    // the number of times that must be reduced during extrapolation 
    // from each fixed sample.
    while (!_hmax.isEmpty()) {

      // Remove from the max-heap the fixed sample with largest time.
      TimeHeap2.Entry ef = _hmax.remove();
      int k1 = ef.i1;
      int k2 = ef.i2;

      // The values to be extrapolated.
      float[] vk = _vk[k2][k1];

      // Mark all samples as far, mark the fixed sample as known with 
      // time zero, and update its neighbors.
      clearMarks();
      _hmin.clear();
      _mark[k2][k1] = _known;
      _k1[k2][k1] = k1;
      _k2[k2][k1] = k2;
      _tk[k2][k1] = 0.0f;
      updateNabors(k1,k2,null);

      // Extrapolate from the fixed sample to all samples that are
      // nearer to the fixed sample than to any other fixed sample.
      while (!_hmin.isEmpty()) {
        TimeHeap2.Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        float t = e.t;
        _mark[i2][i1] = _known;
        if (_type[i2][i1]==FIXED) {
          _hmax.reduce(i1,i2,t);
        } else {
          _type[i2][i1] = EXTRA;
          _k1[i2][i1] = k1;
          _k2[i2][i1] = k2;
          _vk[i2][i1] = vk;
        }
        updateNabors(i1,i2,null);
      }
    }
  }

  /**
   * Interpolates values from all fixed samples, using extrapolated samples.
   * After interpolation, all samples are either fixed or interpolated.
   */
  public void interpolate() {

    // Interpolation occurs in two stages. Both stages compute times by 
    // fast marching away from the sample to be interpolated. In stage 1, 
    // times and values for samples reached while marching are modified,
    // so that each sample interpolated in this first stage will affect
    // samples interpolated later. In stage 2, times (but not values) 
    // are again modified during marching, but are restored after marching. 
    // Therefore, values interpolated in stage 2 do not affect other 
    // interpolated values. 
    boolean stage1 = true;
    boolean stage2 = false;

    // In stage 2, times that must be restored are saved in a list while 
    // marching, and values interpolated are saved in a separate array so
    // that they will not affect other interpolated values. At the end of
    // stage 2, we will merge the two arrays of values. In stage 1, both
    // the time list and values array are null and not used.
    TimeList tl = null;
    float[][][] va = null;

    // Minimum number of samples to interpolate in stage 2, a fraction 
    // of the the total number of samples. This is an important parameter. 
    // Higher fractions close to one yield smoother interpolations, but
    // can be much more costly than lower fractions.
    //int nstage2 = (int)(0.05*_n1*_n2);
    int nstage2 = (int)(0.9*_n1*_n2);

    // Insert all extrapolated samples into the max-heap with their
    // current times. After the max-heap is built, the extrapolated
    // sample with largest time is at the top of the heap.
    _hmax.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_type[i2][i1]==EXTRA) {
          _hmax.insert(i1,i2,_tk[i2][i1]);
        }
      }
    }

    // Interpolate all extrapolated (not-fixed) samples, one at a time, 
    // in order of decreasing time. The extrapolated sample with the 
    // largest time is at the top of the max-heap. In stage 1, as we 
    // interpolate this sample, times for other extrapolated samples 
    // that remain in the heap may be reduced, so that their order may 
    // change. We choose a decreasing order in stage 1 to reduce the 
    // number of extrapolated samples that are reached while marching.
    // In stage 2, the order will not matter.
    while (!_hmax.isEmpty()) {

      // Switch to stage 2 when number of samples left is small enough. 
      if (stage1 && _hmax.size()<nstage2) {
        stage1 = false;
        stage2 = true;
        tl = new TimeList();
        va = new float[_n2][_n1][];
      }

      // Remove from the max-heap the extrapolated sample with largest time.
      // This is the sample to be interpolated; the "interpolated sample".
      TimeHeap2.Entry te = _hmax.remove();
      int k1 = te.i1;
      int k2 = te.i2;
      float tk = te.t;

      // Initial weight and sum of weights in weighted sum of values.
      /*
      float[] ds = new float[3];
      _st.getTensor(k1,k2,ds);
      float d11 = ds[0], d12 = ds[1], d22 = ds[2];
      float wi = sqrt(d11*d22-d12*d12);
      */
      float wi = 1.0f;
      float ws = wi;

      // The values to be interpolated. In stage 1, this array will be
      // assigned to all extrapolated samples during the march away from 
      // the interpolated sample. This assignment is one reason that an 
      // array of values is so useful, for we will not actually know the 
      // interpolated values until the march is complete. Later, when we 
      // have completed the computation of the interpolated values, those 
      // values will already be referenced by all extrapolated samples 
      // nearest to the interpolated sample.
      // In stage 2, interpolated values are not extrapolated, and must not
      // affect other interpolated values, so we store them in a separate
      // array of values to be merged later.
      float[] vk = Array.mul(wi,_vk[k2][k1]);
      if (stage1) {
        _vk[k2][k1] = vk;
      } else {
        va[k2][k1] = vk;
      }

      // In stage 2, save the time for the interpolated sample.
      if (stage2) {
        tl.clear();
        tl.append(k1,k2,tk);
      }

      // Mark all samples as far, set the type of the interpolated sample,
      // mark the interpolated sample as known with time zero, and update 
      // its neighbors.
      clearMarks();
      _hmin.clear();
      _type[k2][k1] = INTER;
      _mark[k2][k1] = _known;
      _k1[k2][k1] = k1;
      _k2[k2][k1] = k2;
      _tk[k2][k1] = 0.0f;
      updateNabors(k1,k2,tl);

      // March away from the interpolated sample to all extrapolated
      // samples that are nearer to the interpolated sample than to any 
      // fixed samples or samples previously interpolated in stage 1.
      // While marching, accumulate values needed for interpolation.
      while (!_hmin.isEmpty()) {

        // Get the extrapolated sample with minimum time.
        TimeHeap2.Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        float ti = e.t;

        // Accumulate values for the extrapolated sample.
        /*
        _st.getTensor(i1,i2,ds);
        d11 = ds[0]; d12 = ds[1]; d22 = ds[2];
        wi = sqrt(d11*d22-d12*d12);
        */
        wi = 1.0f;
        ws += wi;
        float[] vki = _vk[i2][i1];
        for (int iv=0; iv<_nv; ++iv)
          vk[iv] += wi*vki[iv];

        // Mark the extrapolated sample known. In stage 1, reduce it's
        // time in the max-heap. Also, in stage 1, it's values will be 
        // those of the interpolated sample. Continue marching by updating 
        // the nabor samples.
        _mark[i2][i1] = _known;
        if (stage1) {
          _hmax.reduce(i1,i2,ti);
          _k1[i2][i1] = k1;
          _k2[i2][i1] = k2;
          _vk[i2][i1] = vk;
        }
        updateNabors(i1,i2,tl);
      }

      // The march is complete, and nearby values have been accumulated.
      // Now simply divide by the sum of the weights.
      for (int iv=0; iv<_nv; ++iv)
        vk[iv] /= ws;

      // In stage 2, restore any times saved during marching.
      if (stage2) {
        int nl = tl.n;
        int[] k1l = tl.k1List;
        int[] k2l = tl.k2List;
        float[] tkl = tl.tkList;
        for (int il=0; il<nl; ++il)
          _tk[k2l[il]][k1l[il]] = tkl[il];
      }
    }

    // Finally, merge any interpolated values that were saved in stage 2.
    if (va!=null) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          float[] vi = va[i2][i1];
          if (vi!=null)
            _vk[i2][i1] = vi;
        }
      }
    }
  }

  /**
   * Gets the array of times used in this painting.
   * These times are modified by extrapolation and interpolation.
   * @return array of times; by reference, not by copy.
   */
  public float[][] getTimes() {
    return _tk;
  }

  /**
   * Gets the painted value with index zero for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   */
  public float getValue(int i1, int i2) {
    return getValue(i1,i2,0);
  }

  /**
   * Gets the painted value for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param iv index of value to get.
   */
  public float getValue(int i1, int i2, int iv) {
    float[] vk = _vk[i2][i1];
    return (vk!=null)?vk[iv]:_dv[iv];
  }

  /**
   * Gets a copy of the painted values with index zero.
   * Returns default values for any samples not painted.
   * @return array of values; by copy, not by reference.
   */
  public float[][] getValues() {
    return getValues(0);
  }

  /**
   * Gets a copy of the painted values with specified index.
   * Returns default values for any samples not painted.
   * @param iv index of values to get.
   * @return array of values; by copy, not by reference.
   */
  public float[][] getValues(int iv) {
    float[][] v = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        v[i2][i1] = getValue(i1,i2,iv);
      }
    }
    return v;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // The value for times not yet computed. Also the value returned by
  // methods that compute times from nabor times when a valid time
  // cannot be computed. We use the maximum possible float so that
  // it will be larger than any valid times we compute.
  private static final float TIME_INVALID = Float.MAX_VALUE;

  // Type of paint.
  private static final byte CLEAR = 0; // values null (not painted)
  private static final byte FIXED = 1; // values painted explicitly
  private static final byte EXTRA = 2; // values painted by extrapolation
  private static final byte INTER = 3; // values painted by interpolation

  // Marks used during computation of times. For efficiency, we do not
  // loop over all the marks to clear them before beginning a fast 
  // marching loop. Instead, we simply modify the mark values.
  private int _far = 0; // samples with no time
  private int _trial = 1; // samples in min-heap with a proposed time
  private int _known = 2; // samples with a known time
  private void clearMarks() {
    if (_known+2>Integer.MAX_VALUE) { // if we must loop over all marks, ...
      _far = 0;
      _trial = 1;
      _known = 2;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _mark[i2][i1] = _far;
        }
      }
    } else {
      _far += 2; // all known samples instantly become far samples
      _trial +=2; // no samples are trial
      _known +=2; // no samples are known
    }
  }

  private Tensors _st; // the structure tensor field
  private int _n1,_n2; // painting dimensions
  private int _nv; // number of values associated with each sample
  private float[] _dv; // default values for clear samples
  private float[][][] _vk; // painted values; null for clear samples
  private float[][] _tk; // time to nearest fixed or interpolated sample
  private int[][] _k1,_k2; // indices of nearest fixed or interp sample
  private int[][] _mark; // samples are marked far, trial, or known
  private byte[][] _type; // fixed, extra, 
  private TimeHeap2 _hmin; // the min heap
  private TimeHeap2 _hmax; // the max heap

  // Times for each sample are computed from one of eight nabor triangles.
  // These triangles are indexed as follows:
  //       2 ^
  //   * - - * - - *
  //   | \ 2 | 1 / | 
  //   | 3 \ | / 0 |
  //   * - - X - - * >
  //   | 4 / | \ 7 | 1
  //   | / 5 | 6 \ | 
  //   * - - * - - *
  // The symbol X represents the vertex X0 shared by all eight triangles. 
  // The symbol * represents the other two triangle vertices X1 and X2, 
  // which are indexed in counter-clockwise order around X0.

  // Sample index offsets for vertices X1 of the eight nabor triangles.
  private static final int[] K11 = { 1, 1, 0,-1,-1,-1, 0, 1};
  private static final int[] K12 = { 0, 1, 1, 1, 0,-1,-1,-1};

  // Sample index offsets for vertices X2 of the eight nabor triangles.
  private static final int[] K21 = { 1, 0,-1,-1,-1, 0, 1, 1};
  private static final int[] K22 = { 1, 1, 1, 0,-1,-1,-1, 0};

  // Components of vectors Y1 = X1-X2 for the eight nabor triangles.
  private static final float[] Y11 =
    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,-1.0f,-1.0f, 0.0f};
  private static final float[] Y12 =
    {-1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,-1.0f};

  // Components of vectors Y2 = X0-X2 for the eight nabor triangles.
  private static final float[] Y21 =
    {-1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f,-1.0f,-1.0f};
  private static final float[] Y22 =
    {-1.0f,-1.0f,-1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f};

  // Structure tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, float[] s) {
      s[0] = 1.0f; // s11
      s[1] = 0.0f; // s12
      s[2] = 1.0f; // s22
    }
  }

  private void updateNabors(int i1, int i2, TimeList tl) {

    // For all eight nabors of specified sample at (i1,i2) ...
    for (int k=0; k<8; ++k) {
      int k1 = K11[k];
      int k2 = K12[k];

      // Sample indices (j1,j2) for this nabor; skip if out of bounds.
      int j1 = i1+k1;
      int j2 = i2+k2;
      if (j1<0 || j1>=_n1) continue;
      if (j2<0 || j2>=_n2) continue;

      // If time for nabor not already known, update it.
      if (_mark[j2][j1]!=_known)
        updateTime(i1,i2,j1,j2,tl);
    }
  }

  /**
   * Updates the time for one sample using times at eight nabors.
   * If not null, the time list is used to store times before they
   * are updated, so that they can later be restored.
   */
  private void updateTime(int i1, int i2, int j1, int j2, TimeList tl) {

    // Elements of structure tensor.
    float[] s = new float[3];
    _st.getTensor(j1,j2,s);
    float s11 = s[0];
    float s12 = s[1];
    float s22 = s[2];

    // The current minimum time.
    float tmin = _tk[j2][j1];

    // Initally assume that no computed time will be less than current min.
    boolean smallerTimeFound = false;

    // For all eight nabor triangles, ...
    for (int jt=0; jt<8; ++jt) {

      // Sample indices of vertices X1 and X2 of nabor triangle.
      int j11 = j1+K11[jt];
      int j12 = j2+K12[jt];
      int j21 = j1+K21[jt];
      int j22 = j2+K22[jt];

      // Either X1 or X2 must be newly known.
      if (j11!=i1 && j21!=i1 || j12!=i2 && j22!=i2) continue;

      // All indices must be in bounds.
      if (j11<0 || j11>=_n1) continue;
      if (j12<0 || j12>=_n2) continue;
      if (j21<0 || j21>=_n1) continue;
      if (j22<0 || j22>=_n2) continue;

      // Need at least one nabor with known time.
      int m1 = _mark[j12][j11];
      int m2 = _mark[j22][j21];
      if (m1!=_known && m2!=_known) continue;

      // Times T0, T1 and T2 at vertices X0, X1 and X2 of nabor triangle.
      float t0 = TIME_INVALID;
      float t1 = _tk[j12][j11];
      float t2 = _tk[j22][j21];

      // Components of vectors Y1 = X1-X2 and Y2 = X0-X2.
      float y11 = Y11[jt];
      float y12 = Y12[jt];
      float y21 = Y21[jt];
      float y22 = Y22[jt];

      // Inner products with respect to metric tensor S.
      float d11 = y11*s11*y11+y11*s12*y12+y12*s12*y11+y12*s22*y12;
      float d12 = y11*s11*y21+y11*s12*y22+y12*s12*y21+y12*s22*y22;
      float d22 = y21*s11*y21+y21*s12*y22+y22*s12*y21+y22*s22*y22;

      // Time T0 computed for one nabor triangle.
      if (m1!=_known) {
        t0 = t2+sqrt(d22); // a = 0
      } else if (m2!=_known) {
        t0 = t1+sqrt(d22-2.0f*d12+d11); // a = 1
      } else {
        float u1 = t1-t2;
        float u2 = t2;
        float dd = d11*d22-d12*d12;
        if (dd<0.0f) dd = 0.0f;
        float du = d11-u1*u1;
        if (du>0.0f) {
          float a = (d12-u1*sqrt(dd/du))/d11;
          if (a<=0.0f) { // a <= 0
            t0 = t2+sqrt(d22);
          } else if (a>=1.0f) { // a >= 1
            t0 = t1+sqrt(d22-2.0f*d12+d11);
          } else { // 0 < a < 1
            float da = d22-a*(2.0f*d12-a*d11);
            if (da<0.0f) da = 0.0f;
            t0 = u2+a*u1+sqrt(d22-2.0f*a*d12+a*a*d11);
          }
        }
      }

      // If computed time T0 is smaller than the min time, update the min time.
      if (t0<tmin) {
        tmin = t0;
        smallerTimeFound = true;
      }
    }

    // If a smaller time has been found, ...
    if (smallerTimeFound) {

      // If not been here before, so this sample not already in the min-heap, 
      // then insert this sample into the min-heap, and if the time list is
      // not null, save the current stored time so it can be restored later.
      if (_mark[j2][j1]!=_trial) {
        _mark[j2][j1] = _trial;
        _hmin.insert(j1,j2,tmin);
        if (tl!=null) 
          tl.append(j1,j2,_tk[j2][j1]);
      }

      // Else, simply reduce the time already stored in the min-heap.
      else {
        _hmin.reduce(j1,j2,tmin);
      }

      // Store the smaller time.
      _tk[j2][j1] = tmin;
    }
  }

  /**
   * For times that must be restored during interpolation.
   */
  private static class TimeList {
    private int n = 0;
    private int[] k1List = new int[1024];
    private int[] k2List = new int[1024];
    private float[] tkList = new float[1024];
    public void append(int k1, int k2, float tk) {
      if (n==tkList.length) {
        int[] k1New = new int[2*k1List.length];
        int[] k2New = new int[2*k2List.length];
        float[] tkNew = new float[2*tkList.length];
        Array.copy(n,k1List,k1New);
        Array.copy(n,k2List,k2New);
        Array.copy(n,tkList,tkNew);
        k1List = k1New;
        k2List = k2New;
        tkList = tkNew;
      }
      k1List[n] = k1;
      k2List[n] = k2;
      tkList[n] = tk;
      ++n;
    }
    public void clear() {
      n = 0;
    }
  }

  /**
   * Interpolates values from all fixed samples, using extrapolated samples.
   * After interpolation, all samples are either fixed or interpolated.
   */
  public void interpolateOldMethod() {

    //Plot plot = new Plot(getValues(),ColorMap.JET); // DEBUG

    // Insert all extrapolated samples into the max-heap with their
    // current times. After the max-heap is built, the extrpolated
    // sample with largest time is at the top of the heap.
    _hmax.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_type[i2][i1]==EXTRA) {
          _hmax.insert(i1,i2,_tk[i2][i1]);
        }
      }
    }

    // Interpolate all extrapolated (not-fixed) samples, one at a time, 
    // in order of decreasing time. The extrapolated sample with the 
    // largest time is at the top of the max-heap. As we interpolate 
    // this sample, times for other extrapolated samples that remain in 
    // the heap may be reduced, so that their order may change. We choose 
    // a decreasing order to reduce the number of extrapolated samples 
    // that must be modified during interpolation.
    while (!_hmax.isEmpty()) {

      // Remove from the max-heap the extrapolated sample with largest time.
      // This is the sample to be interpolated; the "interpolated sample".
      TimeHeap2.Entry ee = _hmax.remove();
      int k1 = ee.i1;
      int k2 = ee.i2;

      // The values to be interpolated. This array will be assigned to all
      // extrapolated samples during the march away from the interpolated
      // sample. This is one reason that an array of values is so useful, 
      // for we will not actually know the interpolated values until the 
      // march is complete. Then, when we compute the interpolated values, 
      // those values will already be referenced by all extrapolated samples
      // nearest to the interpolated sample.
      float[] vk = _vk[k2][k1] = Array.copy(_vk[k2][k1]);

      // Count of values accumulated for the interpolated sample.
      int nk = 1;

      // Mark all samples as far, set the type of the interpolated sample,
      // mark the interpolated sample as known with time zero, and update 
      // its neighbors.
      clearMarks();
      _hmin.clear();
      _type[k2][k1] = INTER;
      _mark[k2][k1] = _known;
      _k1[k2][k1] = k1;
      _k2[k2][k1] = k2;
      _tk[k2][k1] = 0.0f;
      updateNabors(k1,k2,null);

      // March away from the interpolated sample to all extrapolated
      // samples that are nearer to the interpolated sample than to any 
      // other fixed or interpolated samples. While marching, accumulate 
      // values needed for interpolation.
      while (!_hmin.isEmpty()) {

        // Get the extrapolated sample with minimum time.
        TimeHeap2.Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        float t = e.t;

        // Accumulate existing values for the extrapolated sample.
        float[] vki = _vk[i2][i1];
        for (int iv=0; iv<_nv; ++iv)
          vk[iv] += vki[iv];
        ++nk;

        // Mark the extrapolated sample known with reduced time.
        // It's values will be those of the interpolated sample.
        // Continue the march by updating the nabor samples.
        _mark[i2][i1] = _known;
        _hmax.reduce(i1,i2,t);
        _k1[i2][i1] = k1;
        _k2[i2][i1] = k2;
        _vk[i2][i1] = vk;
        updateNabors(i1,i2,null);
      }

      // The march is complete, and nearby values have been accumulated.
      // Now divide by the number of values accumulated.
      float vs = 1.0f/(float)nk;
      for (int iv=0; iv<_nv; ++iv)
        vk[iv] *= vs;
      
      //trace("k1="+k1+" k2="+k2+" nk="+nk);
      //plot(getValues()); // DEBUG
      //plot(_tk); // DEBUG
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void sleep(int ms) {
    try {
      Thread.currentThread().sleep(1000);
    } catch (InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static void plotImageTensors(float[][] x, EigenTensors2 et) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    //sp.setSize(650,600);
    //sp.setSize(950,900);
    sp.setSize(1130,820);
    PixelsView pv = sp.addPixels(x);
    //pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    int n1 = x[0].length;
    int n2 = x.length;
    float[][][] x12 = getTensorEllipses(n1,n2,10,et);
    float[][] x1 = x12[0];
    float[][] x2 = x12[1];
    PointsView ev = new PointsView(x1,x2);
    ev.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    ev.setLineColor(Color.YELLOW);
    sp.getPlotPanel().getTile(0,0).addTiledView(ev);
  }

  private static class Plot {
    Plot(float[][] f) {
      this(f,null);
    }
    Plot(float[][] f, IndexColorModel icm) {
      _sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
      //_sp.setSize(650,600);
      //_sp.setSize(950,900);
      _sp.setSize(1130,820);
      _pv = _sp.addPixels(f);
      if (icm==null) icm = ColorMap.JET;
      _pv.setColorModel(icm);
      //_pv.setInterpolation(PixelsView.Interpolation.NEAREST);
      _pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    }
    void set(final float[][] f) {
      sleep(1000);
      SwingUtilities.invokeLater(new Runnable() {
        public void run() {
          _pv.set(f);
        }
      });
    }
    final private SimplePlot _sp;
    final private PixelsView _pv;
  }

  private static void plot(float[][] f) {
    plot(f,null);
  }

  private static void plot(float[][] f, IndexColorModel icm) {
    new Plot(f,icm);
  }

  private static float[][] readImage(int n1, int n2, String fileName) {
    try {
      //java.nio.ByteOrder bo = java.nio.ByteOrder.LITTLE_ENDIAN;
      java.nio.ByteOrder bo = java.nio.ByteOrder.BIG_ENDIAN;
      ArrayInputStream ais = new ArrayInputStream(fileName,bo);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      return null;
    }
  }

  private static class StructureTensors 
    extends EigenTensors2
    implements PaintingX.Tensors 
  {
    StructureTensors(double sigma, float[][] x) {
      super(x[0].length,x.length);
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] u1 = new float[n2][n1];
      float[][] u2 = new float[n2][n1];
      float[][] su = new float[n2][n1];
      float[][] sv = new float[n2][n1];
      LocalOrientFilter lof = new LocalOrientFilter(sigma);
      lof.apply(x,null,u1,u2,null,null,su,sv,null);
      float[][] sc = Array.sub(1.0f,coherence(sigma,x));
      su = Array.mul(su,sc);
      sv = Array.mul(sv,sc);
      plot(su);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          setEigenvalues(i1,i2,su[i2][i1],sv[i2][i1]);
          setEigenvectorU(i1,i2,u1[i2][i1],u2[i2][i1]);
        }
      }
    }
  }

  private static float[][] coherence(double sigma, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    LocalOrientFilter lof1 = new LocalOrientFilter(sigma);
    LocalOrientFilter lof2 = new LocalOrientFilter(sigma*4);
    float[][] u11 = new float[n2][n1];
    float[][] u21 = new float[n2][n1];
    float[][] su1 = new float[n2][n1];
    float[][] sv1 = new float[n2][n1];
    float[][] u12 = new float[n2][n1];
    float[][] u22 = new float[n2][n1];
    float[][] su2 = new float[n2][n1];
    float[][] sv2 = new float[n2][n1];
    lof1.apply(x,null,u11,u21,null,null,su1,sv1,null);
    lof2.apply(x,null,u12,u22,null,null,su2,sv2,null);
    float[][] c = u11;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u11i = u11[i2][i1];
        float u21i = u21[i2][i1];
        float su1i = su1[i2][i1];
        float sv1i = sv1[i2][i1];
        float u12i = u12[i2][i1];
        float u22i = u22[i2][i1];
        float su2i = su2[i2][i1];
        float sv2i = sv2[i2][i1];
        float s111 = (su1i-sv1i)*u11i*u11i+sv1i;
        float s121 = (su1i-sv1i)*u11i*u21i     ;
        float s221 = (su1i-sv1i)*u21i*u21i+sv1i;
        float s112 = (su2i-sv2i)*u12i*u12i+sv2i;
        float s122 = (su2i-sv2i)*u12i*u22i     ;
        float s222 = (su2i-sv2i)*u22i*u22i+sv2i;
        float s113 = s111*s112+s121*s122;
        float s223 = s121*s122+s221*s222;
        float t1 = s111+s221;
        float t2 = s112+s222;
        float t3 = s113+s223;
        float t12 = t1*t2;
        c[i2][i1] = (t12>0.0f)?t3/t12:0.0f;
      }
    }
    return c;
  }

  private static float[][][] getTensorEllipses(
    int n1, int n2, int ns, EigenTensors2 et) 
  {
    int nt = 51;
    int m1 = (n1-1)/ns;
    int m2 = (n2-1)/ns;
    int j1 = (n1-1-(m1-1)*ns)/2;
    int j2 = (n2-1-(m2-1)*ns)/2;
    int nm = m1*m2;
    double r = 0.45*ns;
    float[][] x1 = new float[nm][nt];
    float[][] x2 = new float[nm][nt];
    double dt = 2.0*PI/(nt-1);
    double ft = 0.0f;
    for (int i2=j2,im=0; i2<n2; i2+=ns) {
      double y2 = i2+r;
      for (int i1=j1; i1<n1; i1+=ns,++im) {
        float[] u = et.getEigenvectorU(i1,i2);
        float[] s = et.getEigenvalues(i1,i2);
        double u1 = u[0];
        double u2 = u[1];
        double v1 = -u2;
        double v2 =  u1;
        double su = s[0];
        double sv = s[1];
        double a = r*sqrt(sv/su);
        double b = r;
        for (int it=0; it<nt; ++it) {
          double t = ft+it*dt;
          double cost = cos(t);
          double sint = sin(t);
          x1[im][it] = (float)(i1+a*cost*u1-b*sint*u2);
          x2[im][it] = (float)(i2+b*sint*u1+a*cost*u2);
        }
      }
    }
    return new float[][][]{x1,x2};
  }

  private static int[] K1_SEIS = {
      43,  87, 110, 147, 167, 188,
      39,  82, 105, 142, 160, 180,
      40,  86, 107, 144, 165, 185,
  };
  private static int[] K2_SEIS = {
     100, 100, 100, 100, 100, 100,
     170, 170, 170, 170, 170, 170,
     280, 280, 280, 280, 280, 280,
  };
  private static double[] VK_SEIS = {
       1,   2,   3,   4,   5,   6,
       1,   2,   3,   4,   5,   6,
       1,   2,   3,   4,   5,   6,
  };
  private static int[] K1_SEISB = {
      20,  75, 100, 130, 155, 180, 210,
      20,  70,  95, 142, 160, 180,
      40,  86, 107, 144, 165, 185,
  };
  private static int[] K2_SEISB = {
     100, 100, 100, 100, 100, 100,
     170, 170, 170, 170, 170, 170,
     280, 280, 280, 280, 280, 280,
  };
  private static double[] VK_SEISB = {
       1,   2,   3,   4,   5,   6,
       1,   2,   3,   4,   5,   6,
       1,   2,   3,   4,   5,   6,
  };
  private static int[] getK1SeisAboveBelow() {
    int nk = K1_SEIS.length;
    int[] k1 = new int[2*nk];
    for (int ik=0,jk=0; ik<nk; ++ik,jk+=2) {
      k1[jk  ] = K1_SEIS[ik]-4;
      k1[jk+1] = K1_SEIS[ik]+4;
    }
    return k1;
  }
  private static int[] getK2SeisAboveBelow() {
    int nk = K2_SEIS.length;
    int[] k2 = new int[2*nk];
    for (int ik=0,jk=0; ik<nk; ++ik,jk+=2) {
      k2[jk  ] = K2_SEIS[ik];
      k2[jk+1] = K2_SEIS[ik];
    }
    return k2;
  }
  private static float[] getVkSeisLayer() {
    int nk = VK_SEIS.length;
    float[] vk = new float[2*nk];
    for (int ik=0,jk=0; ik<nk; ++ik,jk+=2) {
      int vka = (int)(VK_SEIS[ik]-0.01);
      int vkb = (int)(VK_SEIS[ik]+0.01);
      vk[jk  ] = (float)vka;
      vk[jk+1] = (float)vkb;
    }
    return vk;
  }

  private static void insertFaults(StructureTensors st) {
    int[] k1 = {  25, 193,  27, 106, 107, 200};
    int[] k2 = { 135, 145, 232, 213, 213, 206};
    int nk = k1.length;
    float ehuge = 10000.0f;
    for (int ik=0; ik<nk; ik+=2) {
      int k1a = k1[ik  ];
      int k1b = k1[ik+1];
      int k2a = k2[ik  ];
      int k2b = k2[ik+1];
      double kscale = (double)(k2b-k2a)/(double)(k1b-k1a);
      for (int j1=k1a; j1<=k1b; ++j1) {
        double x2 = k2a+(j1-k1a)*kscale; 
        int j2l = (int)x2;
        int j2r = j2l+1;
        float[] el = st.getEigenvalues(j1,j2l);
        float[] er = st.getEigenvalues(j1,j2r);
        float e2r = (float)(x2-j2l);
        float e2l = (1.0f-e2r);
        el[0] *= e2l*ehuge;
        el[1] *= e2l*ehuge;
        er[0] *= e2r*ehuge;
        er[1] *= e2r*ehuge;
        st.setEigenvalues(j1,j2l,el);
        st.setEigenvalues(j1,j2r,er);
      }
    }
  }

  private static void testSeismic() {
    int n1 = 251;
    int n2 = 357;
    int nv = 1;
    float[][] x = readImage(n1,n2,"/data/seis/tp/tp73.dat");
    /*
    int m1 = 20;
    int nk = 1+(n1-1)/m1;
    int[] k1 = new int[nk];
    int[] k2 = new int[nk];
    float[] vk = new float[nk];
    for (int i1=0,ik=0; i1<n1; i1+=m1,++ik) {
      k1[ik] = i1;
      k2[ik] = n2/2;
      //vk[ik] = x[k2[ik]][k1[ik]];
      vk[ik] = (float)i1;
      //vk[ik] = (ik%2==0)?1.0f:2.0f;
    }
    */
    StructureTensors st = new StructureTensors(3,x);
    plotImageTensors(x,st);
    PaintingX p = new PaintingX(n1,n2,nv,st);
    int[] k1 = getK1SeisAboveBelow();
    int[] k2 = getK2SeisAboveBelow();
    float[] vk = getVkSeisLayer();
    int nk = k1.length;
    for (int ik=0; ik<nk; ++ik) {
      p.paintAt(k1[ik],k2[ik],vk[ik]);
    }
    plot(x,ColorMap.GRAY);
    p.extrapolate();
    plot(p.getTimes(),ColorMap.JET);
    plot(p.getValues(),ColorMap.JET);
    p.interpolate();
    plot(p.getValues(),ColorMap.JET);
  }

  private static void testChannels() {
    int n1 = 200;
    int n2 = 200;
    int nv = 1;
    float[][] x = readImage(n1,n2,"/data/seis/joe/x174.dat");
    plot(x,ColorMap.GRAY);
    StructureTensors st = new StructureTensors(8,x);

    /*
    int[] k1 =   {  92,  92,  92, 100, 100, 100,  60,  25,  20,  19};
    int[] k2 =   { 109, 102, 116, 132, 125, 139, 116, 124, 110, 117};
    float[] vk = {1.0f,2.0f,2.0f,1.0f,2.0f,2.0f,1.0f,1.0f,1.0f,2.0f};
    int nk = vk.length;
    */
    /*
    int[] k1 =    { 34,  92, 172,  27,  25,  12,  81, 117,  94,  14,  44};
    int[] k2 =    { 81, 109, 109, 111, 124, 138, 146,  82, 122,  99, 162};
    float[] vk = {1.0f,2.0f,2.0f,2.0f,2.0f,3.0f,3.0f,0.0f,0.0f,0.0f,0.0f};
    int nk = vk.length;
    */
    int m2 = 5;
    int nk = 1+(n2-1)/m2;
    int[] k1 = new int[nk];
    int[] k2 = new int[nk];
    float[] vk = new float[nk];
    for (int i2=0,ik=0; i2<n2; i2+=m2,++ik) {
      k1[ik] = 130;
      k2[ik] = i2;
      vk[ik] = (float)i2;
    }
    /*
    int[] k1 =   {  n1-1,  n1-1};
    int[] k2 =   {1*n2/4,3*n2/4};
    float[] vk = {  1.0f,  2.0f};
    int nk = vk.length;
    */
    PaintingX p = new PaintingX(n1,n2,nv,st);
    for (int ik=0; ik<nk; ++ik) {
      p.paintAt(k1[ik],k2[ik],vk[ik]);
    }
    float[][] v;
    p.extrapolate();
    v = p.getValues();
    plot(v,ColorMap.JET);
    p.interpolate();
    v = p.getValues();
    plot(v,ColorMap.JET);
  }

  private static class SimpleTensors 
    extends EigenTensors2 
    implements PaintingX.Tensors 
  {
    SimpleTensors(int n1, int n2, double su, double sv, double v1) {
      super(n1,n2);
      float u2 = -(float)v1;
      float u1 = sqrt(1.0f-u2*u2);
      float au = (float)su;
      float av = (float)sv;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float d1 = (float)(i1-n1/2);
          float d2 = (float)(i2-n2/2);
          float as = exp(-0.0001f*(d1*d1+d2*d2));
          if (i1==n1/3 && 1*n2/4<i2 && i2<3*n2/4)
            as = 1000.0f;
          else
            as = 1.0f;
          setEigenvalues(i1,i2,au*as,av*as);
          setEigenvectorU(i1,i2,u1,u2);
        }
      }
    }
  }

  private static void testIsotropic() {
    int n1 = 301;
    int n2 = 301;
    int nv = 1;
    //float su = 4.0f;
    float su = 1.0f;
    float sv = 1.0f;
    float v1 = sin(0.0f*FLT_PI/8.0f);
    SimpleTensors st = new SimpleTensors(n1,n2,su,sv,v1);
    PaintingX p = new PaintingX(n1,n2,nv,st);
    p.paintAt(1*n1/4,1*n2/4,1.0f);
    p.paintAt(3*n1/4,3*n2/4,2.0f);
    /*
    p.paintAt(   1,   1,1.0f);
    p.paintAt(n1-1,   1,1.0f);
    p.paintAt(   1,n2-1,1.0f);
    p.paintAt(n1-1,n2-1,1.0f);
    p.paintAt(n1/2,n2/2,2.0f);
    */
    p.extrapolate();
    plotImageTensors(p.getTimes(),st);
    plot(p.getValues());
    p.interpolate();
    plot(p.getValues());
  }

  private static float[][] makeTargetImage(int n1, int n2) {
    float k = 0.3f;
    float[][] x = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float d2 = (float)(i2-n2/2);
      for (int i1=0; i1<n1; ++i1) {
        float d1 = (float)(i1-n1/2);
        x[i2][i1] = 10.0f*sin(k*sqrt(d1*d1+d2*d2));
      }
    }
    return x;
  }

  private static void testTarget() {
    int n1 = 315;
    int n2 = 315;
    int nv = 1;
    float[][] x = makeTargetImage(n1,n2);
    StructureTensors st = new StructureTensors(8,x);
    PaintingX p = new PaintingX(n1,n2,nv,st);
    int m1 = 1;
    int m2 = 1;
    int nk = 1+(n2-1)/m2+1+(n1-1)/m1;
    int[] k1 = new int[nk];
    int[] k2 = new int[nk];
    float[] vk = new float[nk];
    int ik = 0;
    for (int i2=0; i2<n2; i2+=m2,++ik) {
      k1[ik] = n1/2;
      k2[ik] = i2;
      vk[ik] = x[k2[ik]][k1[ik]];
    }
    for (int i1=0; i1<n1; i1+=m1,++ik) {
      k1[ik] = i1;
      k2[ik] = n2/2;
      vk[ik] = x[k2[ik]][k1[ik]];
    }
    for (ik=0; ik<nk; ++ik)
      p.paintAt(k1[ik],k2[ik],vk[ik]);
    plotImageTensors(x,st);
    p.extrapolate();
    plot(p.getValues());
    p.interpolate();
    plot(p.getValues());
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //testSeismic();
        //testChannels();
        testIsotropic();
        //testTarget();
      }
    });
  }
}
