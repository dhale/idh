/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.concurrent.atomic.*;

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
 * A 3D array of painted values, where most values are painted automatically.
 * Except for a relatively small number of fixed samples painted explicitly, 
 * most samples are painted by extrapolation and interpolation guided by
 * structure tensors. Intuitively, paint flows slowly through locations with
 * relatively high structure. Structure tensors may also cause paint to
 * flow anisotropically, in directions corresponding to relatively low
 * structure.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.07
 */
public class Painting3 {

  /**
   * An interface for classes of structure tensors. Each tensor is a
   * symmetric positive-definite 3-by-3 matrix 
   * {{s11,s12,s13},{s12,s22,s23},{s13,s23,s33}}.
   */
  public interface Tensors {

    /**
     * Gets structure tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param i3 index for 3rd dimension.
     * @param s array {s11,s12,s13,s22,s23,s33} of tensor elements.
     */
    public void getTensor(int i1, int i2, int i3, float[] s);
  }

  /**
   * Constructs a painting with constant identity structure tensors.
   * In this case, time = distance, which is useful for testing.
   * Painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param nv number of values painted for each sample
   */
  public Painting3(int n1, int n2, int n3, int nv) {
    this(n1,n2,n3,nv,new IdentityTensors());
  }
  
  /**
   * Constructs a painting for the specified structure tensor field.
   * All painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param nv number of values painted for each sample
   * @param st structure tensors.
   */
  public Painting3(int n1, int n2, int n3, int nv, Tensors st) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _nv = nv;
    _st = st;
    _dv = new float[nv];
    _k1 = new int[n3][n2][n1];
    _k2 = new int[n3][n2][n1];
    _k3 = new int[n3][n2][n1];
    _tk = new float[n3][n2][n1];
    _vk = new float[n3][n2][n1][];
    _type = new byte[n3][n2][n1];
    _mark = new int[n3][n2][n1];
    _hmin = new TimeHeap3(TimeHeap3.Type.MIN,n1,n2,n3);
    _hmax = new TimeHeap3(TimeHeap3.Type.MAX,n1,n2,n3);
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
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          clear(i1,i2,i3);
        }
      }
    }
  }

  /**
   * Clears painted values that are not fixed.
   */
  public void clearNotFixed() {
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (_type[i3][i2][i1]!=FIXED) {
            clear(i1,i2,i3);
          }
        }
      }
    }
  }

  /**
   * Erases values for any fixed sample with specified indices.
   * If the specified sample is not fixed, this method does nothing.
   * @param i1 index in 1st dimension of sample to erase.
   * @param i2 index in 2nd dimension of sample to erase.
   * @param i3 index in 3rd dimension of sample to erase.
   */
  public void eraseFixedAt(int i1, int i2, int i3) {
    if (_type[i3][i2][i1]==FIXED)
      clear(i1,i2,i3);
  }
  private void clear(int i1, int i2, int i3) {
    _type[i3][i2][i1] = CLEAR;
    _mark[i3][i2][i1] = _known;
    _k1[i3][i2][i1] = -1;
    _k2[i3][i2][i1] = -1;
    _tk[i3][i2][i1] = TIME_INVALID;
    _vk[i3][i2][i1] = null;
  }

  /**
   * Paints the specified sample with one specified value at index zero.
   * Paints default values for indices greater than zero.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param i3 index in 3rd dimension of sample to paint.
   * @param v value at index zero for the painted sample.
   */
  public void paintAt(int i1, int i2, int i3, float v) {
    _type[i3][i2][i1] = FIXED;
    _k1[i3][i2][i1] = i1;
    _k2[i3][i2][i1] = i2;
    _tk[i3][i2][i1] = TIME_INVALID;
    _vk[i3][i2][i1] = new float[_nv];
    _vk[i3][i2][i1][0] = v;
    for (int iv=1; iv<_nv; ++iv)
      _vk[i3][i2][i1][iv] = _dv[iv];
  }

  /**
   * Paints the specified sample with specified values.
   * After painting, the specified sample is fixed; its values will
   * not change in any subsequent extrapolation or interpolation.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param i3 index in 3rd dimension of sample to paint.
   * @param v array of values for painted sample; by copy, not by reference.
   */
  public void paintAt(int i1, int i2, int i3, float[] v) {
    _type[i3][i2][i1] = FIXED;
    _k1[i3][i2][i1] = i1;
    _k2[i3][i2][i1] = i2;
    _tk[i3][i2][i1] = TIME_INVALID;
    _vk[i3][i2][i1] = Array.copy(v);
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
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _tk[i3][i2][i1] = TIME_INVALID;
          if (_type[i3][i2][i1]==FIXED) {
            _hmax.insert(i1,i2,i3,TIME_INVALID);
          } else {
            _type[i3][i2][i1] = CLEAR;
          }
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
      TimeHeap3.Entry ef = _hmax.remove();
      int k1 = ef.i1;
      int k2 = ef.i2;
      int k3 = ef.i3;

      // The values to be extrapolated.
      float[] vk = _vk[k3][k2][k1];

      // Mark all samples as far, mark the fixed sample as known with 
      // time zero, and update its neighbors.
      clearMarks();
      _hmin.clear();
      _mark[k3][k2][k1] = _known;
      _k1[k3][k2][k1] = k1;
      _k2[k3][k2][k1] = k2;
      _k3[k3][k2][k1] = k3;
      _tk[k3][k2][k1] = 0.0f;
      updateNabors(k1,k2,k3,null);

      // Extrapolate from the fixed sample to all samples that are
      // nearer to the fixed sample than to any other fixed sample.
      while (!_hmin.isEmpty()) {
        TimeHeap3.Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        int i3 = e.i3;
        float t = e.t;
        _mark[i3][i2][i1] = _known;
        if (_type[i3][i2][i1]==FIXED) {
          _hmax.reduce(i1,i2,i3,t);
        } else {
          _type[i3][i2][i1] = EXTRA;
          _k1[i3][i2][i1] = k1;
          _k2[i3][i2][i1] = k2;
          _vk[i3][i2][i1] = vk;
        }
        updateNabors(i1,i2,i3,null);
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
    float[][][][] va = null;

    // Minimum number of samples to interpolate in stage 2, a fraction 
    // of the the total number of samples. This is an important parameter. 
    // Higher fractions close to one yield smoother interpolations, but
    // can be much more costly than lower fractions.
    int nstage2 = (int)(0.05*_n1*_n2*_n3);

    // Insert all extrapolated samples into the max-heap with their
    // current times. After the max-heap is built, the extrapolated
    // sample with largest time is at the top of the heap.
    _hmax.clear();
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (_type[i3][i2][i1]==EXTRA) {
            _hmax.insert(i1,i2,i3,_tk[i3][i2][i1]);
          }
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
        va = new float[_n3][_n2][_n1][];
      }

      // Remove from the max-heap the extrapolated sample with largest time.
      // This is the sample to be interpolated; the "interpolated sample".
      TimeHeap3.Entry te = _hmax.remove();
      int k1 = te.i1;
      int k2 = te.i2;
      int k3 = te.i3;
      float tk = te.t;

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
      float[] vk = Array.copy(_vk[k3][k2][k1]);
      if (stage1) {
        _vk[k3][k2][k1] = vk;
      } else {
        va[k3][k2][k1] = vk;
      }

      // In stage 2, save the time for the interpolated sample.
      if (stage2) {
        tl.clear();
        tl.append(k1,k2,k3,tk);
      }

      // Count of values accumulated for the interpolated sample.
      int nk = 1;

      // Mark all samples as far, set the type of the interpolated sample,
      // mark the interpolated sample as known with time zero, and update 
      // its neighbors.
      clearMarks();
      _hmin.clear();
      _type[k3][k2][k1] = INTER;
      _mark[k3][k2][k1] = _known;
      _k1[k3][k2][k1] = k1;
      _k2[k3][k2][k1] = k2;
      _k3[k3][k2][k1] = k3;
      _tk[k3][k2][k1] = 0.0f;
      updateNabors(k1,k2,k3,tl);

      // March away from the interpolated sample to all extrapolated
      // samples that are nearer to the interpolated sample than to any 
      // fixed samples or samples previously interpolated in stage 1.
      // While marching, accumulate values needed for interpolation.
      while (!_hmin.isEmpty()) {

        // Get the extrapolated sample with minimum time.
        TimeHeap3.Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        int i3 = e.i3;
        float ti = e.t;

        // Accumulate values for the extrapolated sample.
        float[] vki = _vk[i3][i2][i1];
        for (int iv=0; iv<_nv; ++iv)
          vk[iv] += vki[iv];
        ++nk;

        // Mark the extrapolated sample known. In stage 1, reduce it's
        // time in the max-heap. Also, in stage 1, it's values will be 
        // those of the interpolated sample. Continue marching by updating 
        // the neighbor samples.
        _mark[i3][i2][i1] = _known;
        if (stage1) {
          _hmax.reduce(i1,i2,i3,ti);
          _k1[i3][i2][i1] = k1;
          _k2[i3][i2][i1] = k2;
          _k3[i3][i2][i1] = k3;
          _vk[i3][i2][i1] = vk;
        }
        updateNabors(i1,i2,i3,tl);
      }

      // The march is complete, and nearby values have been accumulated.
      // Now simply divide by the number of values accumulated.
      float vs = 1.0f/(float)nk;
      for (int iv=0; iv<_nv; ++iv)
        vk[iv] *= vs;

      // In stage 2, restore any times saved during marching.
      if (stage2) {
        int nl = tl.n;
        int[] k1l = tl.k1List;
        int[] k2l = tl.k2List;
        int[] k3l = tl.k3List;
        float[] tkl = tl.tkList;
        for (int il=0; il<nl; ++il)
          _tk[k3l[il]][k2l[il]][k1l[il]] = tkl[il];
      }
    }

    // Finally, merge any interpolated values that were saved in stage 2.
    if (va!=null) {
      for (int i3=0; i3<_n3; ++i3) {
        for (int i2=0; i2<_n2; ++i2) {
          for (int i1=0; i1<_n1; ++i1) {
            float[] vi = va[i3][i2][i1];
            if (vi!=null)
              _vk[i3][i2][i1] = vi;
          }
        }
      }
    }
  }

  /**
   * Gets the array of times used in this painting.
   * These times are modified by extrapolation and interpolation.
   * @return array of times; by reference, not by copy.
   */
  public float[][][] getTimes() {
    return _tk;
  }

  /**
   * Gets the painted value with index zero for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   */
  public float getValue(int i1, int i2, int i3) {
    return getValue(i1,i2,i3,0);
  }

  /**
   * Gets the painted value for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @param iv index of value to get.
   */
  public float getValue(int i1, int i2, int i3, int iv) {
    float[] vk = _vk[i3][i2][i1];
    return (vk!=null)?vk[iv]:_dv[iv];
  }

  /**
   * Gets a copy of the painted values with index zero.
   * Returns default values for any samples not painted.
   * @return array of values; by copy, not by reference.
   */
  public float[][][] getValues() {
    return getValues(0);
  }

  /**
   * Gets a copy of the painted values with specified index.
   * Returns default values for any samples not painted.
   * @param iv index of values to get.
   * @return array of values; by copy, not by reference.
   */
  public float[][][] getValues(int iv) {
    float[][][] v = new float[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          v[i3][i2][i1] = getValue(i1,i2,i3,iv);
        }
      }
    }
    return v;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // The value for times not yet computed. Also the value returned by
  // methods that compute times from neighbor times when a valid time
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
      for (int i3=0; i3<_n3; ++i3) {
        for (int i2=0; i2<_n2; ++i2) {
          for (int i1=0; i1<_n1; ++i1) {
            _mark[i3][i2][i1] = _far;
          }
        }
      }
    } else {
      _far += 2; // all known samples instantly become far samples
      _trial +=2; // no samples are trial
      _known +=2; // no samples are known
    }
  }

  private Tensors _st; // the structure tensor field
  private int _n1,_n2,_n3; // painting dimensions
  private int _nv; // number of values associated with each sample
  private float[] _dv; // default values for clear samples
  private float[][][][] _vk; // painted values; null for clear samples
  private float[][][] _tk; // time to nearest fixed or interpolated sample
  private int[][][] _k1,_k2,_k3; // indices of nearest fixed or interp sample
  private int[][][] _mark; // samples are marked far, trial, or known
  private byte[][][] _type; // sample types: clear, fixed, extra, inter
  private TimeHeap3 _hmin; // the min heap
  private TimeHeap3 _hmax; // the max heap

  // Sample index offsets for 26 neighbor samples.
  private static final int[] K1 = {
    -1, 0, 1, -1, 0, 1, -1, 0, 1,
    -1, 0, 1, -1,    1, -1, 0, 1,
    -1, 0, 1, -1, 0, 1, -1, 0, 1
  };
  private static final int[] K2 = {
    -1,-1,-1,  0, 0, 0,  1, 1, 1,
    -1,-1,-1,  0,    0,  1, 1, 1,
    -1,-1,-1,  0, 0, 0,  1, 1, 1
  };
  private static final int[] K3 = {
    -1,-1,-1, -1,-1,-1, -1,-1,-1,
     0, 0, 0,  0,    0,  0, 0, 0,
     1, 1, 1,  1, 1, 1,  1, 1, 1
  };

  // Sample index offsets for vertices X1 of the 48 neighbor tets.
  private static final int[] K11 = { 1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1};
  private static final int[] K12 = { 0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1};
  private static final int[] K13 = {-1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1};

  // Sample index offsets for vertices X2 of the 48 neighbor tets.
  private static final int[] K21 = { 1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1};
  private static final int[] K22 = { 1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1};
  private static final int[] K23 = {-1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0};

  // Sample index offsets for vertices X3 of the 48 neighbor tets.
  private static final int[] K31 = { 0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1};
  private static final int[] K32 = { 0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0};
  private static final int[] K33 = {-1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0};

  // Indices into offset arrays above for each of 26 neighbor samples.
  // These arrays enable us to loop over only the relevant tets for
  // each of the 26 neighbor samples. For example, for the neighbor with
  // offsets {-1,-1,-1}, the six relevant tets are { 4, 5,20,21,36,37},
  // and the other 42 tets can be ignored. Note that each tet index should
  // appear in these arrays exactly three times.
  private static final int[][] KT  = {
    { 4, 5,20,21,36,37},
    { 5, 6,21,22},
    { 6, 7,22,23,44,45},
    { 3, 4,37,38},
    { 0, 1, 2, 3, 4, 5, 6, 7},
    { 0, 7,45,46},
    { 2, 3,28,29,38,39},
    { 1, 2,29,30},
    { 0, 1,30,31,46,47},
    {19,20,35,36},
    {16,17,18,19,20,21,22,23},
    {16,23,43,44},
    {32,33,34,35,36,37,38,39},
    // no tets for the center of the 3x3x3 cube of 27 samples
    {40,41,42,43,44,45,46,47},
    {27,28,32,39},
    {24,25,26,27,28,29,30,31},
    {24,31,40,47},
    {12,13,18,19,34,35},
    {13,14,17,18},
    {14,15,16,17,42,43},
    {11,12,33,34},
    { 8, 9,10,11,12,13,14,15},
    { 8,15,41,42},
    {10,11,26,27,32,33},
    { 9,10,25,26},
    { 8, 9,24,25,40,41}
  };

  // Structure tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, int i3, float[] s) {
      s[0] = 1.0f; // s11
      s[1] = 0.0f; // s12
      s[2] = 0.0f; // s13
      s[3] = 1.0f; // s22
      s[4] = 0.0f; // s23
      s[5] = 1.0f; // s33
    }
  }

  /**
   * Updates times for all neighbors of sample with indices (i1,i2,i3).
   * If not null, the time list is used to store times before they are
   * updated, so that they can later be restored.
   */
  private void updateNabors(
    final int i1, final int i2, final int i3, final TimeList tl) 
  {
    Thread[] threads = Threads.makeArray();
    int nthread = threads.length;
    trace("updateNabors: nthread="+nthread+" i1="+i1+" i2="+i2+" i3="+i3);
    final AtomicInteger ak = new AtomicInteger();
    for (int ithread=0; ithread<nthread; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int k=ak.getAndIncrement(); k<26; k=ak.getAndIncrement()) {
            trace("  k="+k);
            int k1 = K1[k];
            int k2 = K2[k];
            int k3 = K3[k];

            // Neighbor sample indices (j1,j2,j3); skip if out of bounds.
            int j1 = i1+k1;
            int j2 = i2+k2;
            int j3 = i3+k3;
            if (j1<0 || j1>=_n1) continue;
            if (j2<0 || j2>=_n2) continue;
            if (j3<0 || j3>=_n3) continue;

            // Update time for the neighbor sample.
            updateTime(j1,j2,j3,KT[k],tl);
          }
        }
      });
    }
    Threads.startAndJoin(threads);
  }
  private void updateNaborsX(
    final int i1, final int i2, final int i3, final TimeList tl) 
  {
    for (int k=0; k<26; ++k) {
      int k1 = K1[k];
      int k2 = K2[k];
      int k3 = K3[k];

      // Neighbor sample indices (j1,j2,j3); skip if out of bounds.
      int j1 = i1+k1;
      int j2 = i2+k2;
      int j3 = i3+k3;
      if (j1<0 || j1>=_n1) continue;
      if (j2<0 || j2>=_n2) continue;
      if (j3<0 || j3>=_n3) continue;

      // Update time for the neighbor sample.
      updateTime(j1,j2,j3,KT[k],tl);
    }
  }

  /**
   * Updates the time for one sample using known times in neighbor tets.
   * The array kt contains indices in [0,47] of relevant neighbor tets.
   * If not null, the time list is used to store the time before it
   * is updated, so that it can later be restored.
   */
  private void updateTime(int j1, int j2, int j3, int[] kt, TimeList tl) {

    // Elements of structure tensor.
    float[] s = new float[6];
    _st.getTensor(j1,j2,j3,s);
    float s11 = s[0];
    float s12 = s[1];
    float s13 = s[2];
    float s22 = s[3];
    float s23 = s[4];
    float s33 = s[5];

    // The current minimum time.
    float tmin = _tk[j3][j2][j1];

    // Initally assume that no computed time will be less than current min.
    boolean smallerTimeFound = false;

    // For all relevant neighbor tets, ...
    for (int it=0; it<kt.length; ++it) {
      int jt = kt[it];

      // Sample indices of vertices X0, X1, X2, and X3 of neighbor tet.
      int j01 = j1;
      int j02 = j2;
      int j03 = j3;
      int j11 = j1-K11[jt];
      int j12 = j2-K12[jt];
      int j13 = j3-K13[jt];
      int j21 = j1-K21[jt];
      int j22 = j2-K22[jt];
      int j23 = j3-K23[jt];
      int j31 = j1-K31[jt];
      int j32 = j2-K32[jt];
      int j33 = j3-K33[jt];

      // All indices must be in bounds.
      if (j11<0 || j11>=_n1) continue;
      if (j12<0 || j12>=_n2) continue;
      if (j13<0 || j13>=_n3) continue;
      if (j21<0 || j21>=_n1) continue;
      if (j22<0 || j22>=_n2) continue;
      if (j23<0 || j23>=_n3) continue;
      if (j31<0 || j31>=_n1) continue;
      if (j32<0 || j32>=_n2) continue;
      if (j33<0 || j33>=_n3) continue;

      // Which neighbor vertices are known? (At least one must be!)
      int m1 = _mark[j13][j12][j11];
      int m2 = _mark[j23][j22][j21];
      int m3 = _mark[j33][j32][j31];

      // Times T0, T1, T2 and T3 at vertices X0, X1, X2 and X3 of tet.
      float t0 = TIME_INVALID;
      float t1 = _tk[j13][j12][j11];
      float t2 = _tk[j23][j22][j21];
      float t3 = _tk[j33][j32][j31];

      // Use only known times in {T1,T2,T3} to compute candidate time T0.
      if (m1==_known) {
        if (m2==_known) {
          if (m3==_known) { // use triangle 123
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1-t3,t2-t3,t3,
                             j11-j31,j12-j32,j13-j33,
                             j21-j31,j22-j32,j23-j33,
                             j01-j31,j02-j32,j03-j33);
          } else { // use segment 12
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1-t2,t2,
                             j11-j21,j12-j22,j13-j23,
                             j01-j21,j02-j22,j03-j23);
          }
        } else {
          if (m3==_known) { // use segment 13
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1-t3,t3,
                             j11-j31,j12-j32,j13-j33,
                             j01-j31,j02-j32,j03-j33);
          } else { // use vertex 1
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1,j01-j11,j02-j12,j03-j13);
          }
        }
      } else if (m2==_known) {
        if (m3==_known) { // use segment 23
          t0 = computeTime(s11,s12,s13,s22,s23,s33,
                           t2-t3,t3,
                           j21-j31,j22-j32,j23-j33,
                           j01-j31,j02-j32,j03-j33);
        } else { // use vertex 2
          t0 = computeTime(s11,s12,s13,s22,s23,s33,
                           t2,j01-j21,j02-j22,j03-j23);
        }
      } else if (m3==_known) { // use vertex 3
        t0 = computeTime(s11,s12,s13,s22,s23,s33,
                         t3,j01-j31,j02-j32,j03-j33);
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
      if (_mark[j3][j2][j1]!=_trial) {
        _mark[j3][j2][j1] = _trial;
        _hmin.insert(j1,j2,j3,tmin);
        if (tl!=null) 
          tl.append(j1,j2,j3,_tk[j3][j2][j1]);
      }

      // Else, simply reduce the time already stored in the min-heap.
      else {
        _hmin.reduce(j1,j2,j3,tmin);
      }

      // Store the smaller time.
      _tk[j3][j2][j1] = tmin;
    }
  }

  // Returns the parameter a (alpha) that minimizes
  // t(a) = a*u1+sqrt((y2-a*y1)'*S*(y2-a*y1))
  // subject to the constraint 0 <= a <= 1.
  private static float computeAlpha(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float u1,
    float y11, float y12, float y13,
    float y21, float y22, float y23)
  {
    float z11 = s11*y11+s12*y12+s13*y13;
    float z12 = s12*y11+s22*y12+s23*y13;
    float z13 = s13*y11+s23*y12+s33*y13;
    float z21 = s11*y21+s12*y22+s13*y23;
    float z22 = s12*y21+s22*y22+s23*y23;
    float z23 = s13*y21+s23*y22+s33*y23;
    float d11 = y11*z11+y12*z12+y13*z13;
    float d12 = y11*z21+y12*z22+y13*z23;
    float d22 = y21*z21+y22*z22+y23*z23;
    return computeAlpha(u1,d11,d12,d22);
  }

  // Returns the parameter a (alpha) that minimizes
  // t(a) = a*u1+sqrt(d22-2.0f*a*d12+a*a*d11)
  // subject to the constraint 0 <= a <= 1.
  private static float computeAlpha(
    float u1, float d11, float d12, float d22) 
  {
    float alpha;
    float dd = d11*d22-d12*d12;
    if (dd<0.0f)
      dd = 0.0f;
    float du = d11-u1*u1;
    if (du<=0.0f) {
      alpha = (u1>=0.0f)?0.0f:1.0f;
    } else {
      alpha = (d12-u1*sqrt(dd/du))/d11;
      if (alpha<=0.0f) {
        alpha = 0.0f;
      } else if (alpha>=1.0f) {
        alpha = 1.0f;
      }
    }
    return alpha;
  }

  // Returns the time t+sqrt(y'*S*y).
  private static float computeTime(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float t, float y1, float y2, float y3)
  {
    float z1 = s11*y1+s12*y2+s13*y3;
    float z2 = s12*y1+s22*y2+s23*y3;
    float z3 = s13*y1+s23*y2+s33*y3;
    float dd = y1*z1+y2*z2+y3*z3;
    return t+sqrt(dd);
  }

  // Returns the minimum time 
  // t(a) = u2+a*u1+sqrt((y2-a*y1)'*S*(y2-a*y1))
  // subject to the constraint 0 <= a <= 1.
  private static float computeTime(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float u1, float u2,
    float y11, float y12, float y13,
    float y21, float y22, float y23)
  {
    float z11 = s11*y11+s12*y12+s13*y13;
    float z12 = s12*y11+s22*y12+s23*y13;
    float z13 = s13*y11+s23*y12+s33*y13;
    float z21 = s11*y21+s12*y22+s13*y23;
    float z22 = s12*y21+s22*y22+s23*y23;
    float z23 = s13*y21+s23*y22+s33*y23;
    float d11 = y11*z11+y12*z12+y13*z13;
    float d12 = y11*z21+y12*z22+y13*z23;
    float d22 = y21*z21+y22*z22+y23*z23;
    float alpha = computeAlpha(u1,d11,d12,d22);
    return u2+alpha*u1+sqrt(d22-2.0f*alpha*d12+alpha*alpha*d11);
  }

  // Returns the minimum time
  // t(a1,a2) = u3+a1*u1+a2*u2+sqrt((y3-a1*y1-a2*y2)'*S*(y3-a1*y1-a2*y2))
  // subject to the constraints 0 <= a1, 0 <= a2, and 0 <= a3 = 1-a1-a2.
  private static float computeTime(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float u1, float u2, float u3,
    float y11, float y12, float y13,
    float y21, float y22, float y23,
    float y31, float y32, float y33)
  {
    // Inner products with respect to metric tensor S.
    float z11 = s11*y11+s12*y12+s13*y13;
    float z12 = s12*y11+s22*y12+s23*y13;
    float z13 = s13*y11+s23*y12+s33*y13;
    float z21 = s11*y21+s12*y22+s13*y23;
    float z22 = s12*y21+s22*y22+s23*y23;
    float z23 = s13*y21+s23*y22+s33*y23;
    float z31 = s11*y31+s12*y32+s13*y33;
    float z32 = s12*y31+s22*y32+s23*y33;
    float z33 = s13*y31+s23*y32+s33*y33;
    float d11 = y11*z11+y12*z12+y13*z13;
    float d12 = y11*z21+y12*z22+y13*z23;
    float d13 = y11*z31+y12*z32+y13*z33;
    float d22 = y21*z21+y22*z22+y23*z23;
    float d23 = y21*z31+y22*z32+y23*z33;
    float d33 = y31*z31+y32*z32+y33*z33;

    // The minimum lies on the line a1*alpha1 + a2*alpha2 + bb = 0.
    float a1 = u1*d12-u2*d11;
    float a2 = u1*d22-u2*d12;
    float bb = u2*d13-u1*d23;
    float aa1 = (a1>=0.0f)?a1:-a1;
    float aa2 = (a2>=0.0f)?a2:-a2;
    float alpha1,alpha2;

    // If abs(a1) < abs(a2), solve for alpha1 first, then alpha2.
    if (aa1<aa2) {
      float aoa = a1/a2;
      float boa = bb/a2;
      float v1 = u1-aoa*u2;
      float w21 = y31+boa*y21;
      float w22 = y32+boa*y22;
      float w23 = y33+boa*y23;
      float w11 = y11-aoa*y21;
      float w12 = y12-aoa*y22;
      float w13 = y13-aoa*y23;
      alpha1 = computeAlpha(s11,s12,s13,s22,s23,s33,
                            v1,w11,w12,w13,w21,w22,w23);
      alpha2 = -boa-aoa*alpha1;
    }

    // Else if abs(a2) < abs(a1), solve for alpha2 first, then alpha1.
    else if (aa2<aa1) {
      float aoa = a2/a1;
      float boa = bb/a1;
      float v1 = u2-aoa*u1;
      float w21 = y31+boa*y11;
      float w22 = y32+boa*y12;
      float w23 = y33+boa*y13;
      float w11 = y21-aoa*y11;
      float w12 = y22-aoa*y12;
      float w13 = y23-aoa*y13;
      alpha2 = computeAlpha(s11,s12,s13,s22,s23,s33,
                            v1,w11,w12,w13,w21,w22,w23);
      alpha1 = -boa-aoa*alpha2;
    } 
    
    // Else if a1 = a2 = bb = 0, solve easily for both alpha1 and alpha2.
    else {
      float dd = d11*d22-d12*d12;
      alpha1 = (d12*d23-d13*d22)/dd;
      alpha2 = (d13*d12-d23*d11)/dd;
    }

    // The sum alpha1 + alpha2 + alpha3 = 1.
    float alpha3 = 1.0f-alpha1-alpha2;

    // Initial time is huge.
    float t = TIME_INVALID;

    // If minimum is strictly inside the triangle 123, ...
    if (alpha1>0.0f && alpha2>0.0f && alpha3>0.0f) {
      t = u3+alpha1*u1+alpha2*u2 + 
        sqrt(d33+alpha1*alpha1*d11+alpha2*alpha2*d22 +
             2.0f*alpha1*alpha2*d12-2.0f*alpha1*d13-2.0f*alpha2*d23);

    // Else if the minimum is not strictly inside the triangle 123,
    // search for the minimum on the edges of the triangle.
    } else {

      // If minimum could be on the edge 23 (alpha1=0), ...
      if (alpha1<=0.0f) {
        alpha2 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              u2,y21,y22,y23,y31,y32,y33);
        float t23 = u3+alpha2*u2+sqrt(d33+alpha2*alpha2*d22-2.0f*alpha2*d23);
        if (t23<t) t = t23;
      }

      // If minimum could be on the edge 13 (alpha2=0), ...
      if (alpha2<=0.0f) {
        alpha1 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              u1,y11,y12,y13,y31,y32,y33);
        float t13 = u3+alpha1*u1+sqrt(d33+alpha1*alpha1*d11-2.0f*alpha1*d13);
        if (t13<t) t = t13;
      }

      // If minimum could be on the edge 12 (alpha3=0), ...
      if (alpha3<=0.0f) {
        alpha1 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              u1-u2,
                              y11-y21,y12-y22,y13-y23,
                              y31-y21,y32-y22,y33-y23);
        alpha2 = 1.0f-alpha1;
        float t12 = u3+alpha1*u1+alpha2*u2 + 
          sqrt(d33+alpha1*alpha1*d11+alpha2*alpha2*d22 +
               2.0f*alpha1*alpha2*d12-2.0f*alpha1*d13-2.0f*alpha2*d23);
        if (t12<t) t = t12;
      }
    }
    return t;
  }

  /**
   * For times that must be restored during interpolation.
   */
  private static class TimeList {
    private int n = 0;
    private int[] k1List = new int[1024];
    private int[] k2List = new int[1024];
    private int[] k3List = new int[1024];
    private float[] tkList = new float[1024];
    public synchronized void append(int k1, int k2, int k3, float tk) {
      if (n==tkList.length) {
        int[] k1New = new int[2*k1List.length];
        int[] k2New = new int[2*k2List.length];
        int[] k3New = new int[2*k3List.length];
        float[] tkNew = new float[2*tkList.length];
        Array.copy(n,k1List,k1New);
        Array.copy(n,k2List,k2New);
        Array.copy(n,k3List,k3New);
        Array.copy(n,tkList,tkNew);
        k1List = k1New;
        k2List = k2New;
        k3List = k3New;
        tkList = tkNew;
      }
      k1List[n] = k1;
      k2List[n] = k2;
      k3List[n] = k3;
      tkList[n] = tk;
      ++n;
    }
    public synchronized void clear() {
      n = 0;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static float[][][] readImage(
    int n1, int n2, int n3, String fileName) 
  {
    try {
      //java.nio.ByteOrder bo = java.nio.ByteOrder.LITTLE_ENDIAN;
      java.nio.ByteOrder bo = java.nio.ByteOrder.BIG_ENDIAN;
      ArrayInputStream ais = new ArrayInputStream(fileName,bo);
      float[][][] x = new float[n3][n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      return null;
    }
  }

  private static void plot(float[][][] x) {
  }

  private static class StructureTensors 
    extends EigenTensors3
    implements Painting3.Tensors 
  {
    StructureTensors(double sigma, float[][][] x) {
      super(x[0][0].length,x[0].length,x.length,false);
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      float[][][] u1 = new float[n3][n2][n1];
      float[][][] u2 = new float[n3][n2][n1];
      float[][][] u3 = new float[n3][n2][n1];
      float[][][] w1 = new float[n3][n2][n1];
      float[][][] w2 = new float[n3][n2][n1];
      float[][][] w3 = new float[n3][n2][n1];
      float[][][] su = new float[n3][n2][n1];
      float[][][] sv = new float[n3][n2][n1];
      float[][][] sw = new float[n3][n2][n1];
      LocalOrientFilter lof = new LocalOrientFilter(sigma);
      lof.apply(x,null,null,
                u1,u2,u3,
                null,null,null,
                w1,w2,w3,
                su,sv,sw,
                null,null);
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float sui = su[i3][i2][i1];
            float svi = sv[i3][i2][i1];
            float swi = sw[i3][i2][i1];
            float u1i = u1[i3][i2][i1];
            float u2i = u2[i3][i2][i1];
            float u3i = u3[i3][i2][i1];
            float w1i = w1[i3][i2][i1];
            float w2i = w2[i3][i2][i1];
            float w3i = w3[i3][i2][i1];
            setEigenvalues(i1,i2,i3,sui,svi,swi);
            setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
            setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
          }
        }
      }
    }
  }

  private static class ConstantTensors 
    extends EigenTensors3
    implements Painting3.Tensors 
  {
    ConstantTensors(
      int n1, int n2, int n3,
      float su, float sv, float sw) 
    {
      super(n1,n2,n3,false);
      float u1 = 1.0f, u2 = 0.0f, u3 = 0.0f;
      float v1 = 0.0f, v2 = 1.0f, v3 = 0.0f;
      float w1 = 0.0f, w2 = 0.0f, w3 = 1.0f;
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            setEigenvalues(i1,i2,i3,su,sv,sw);
            setEigenvectorU(i1,i2,i3,u1,u2,u3);
            setEigenvectorW(i1,i2,i3,w1,w2,w3);
          }
        }
      }
    }
  }

  private static void testConstant() {
    int n1 = 25;
    int n2 = 25;
    int n3 = 25;
    int nv = 1;
    float su = 1.0f;
    float sv = 4.0f;
    float sw = 1.0f;
    ConstantTensors ct = new ConstantTensors(n1,n2,n3,su,sv,sw);
    Painting3 p = new Painting3(n1,n2,n3,nv,ct);
    p.paintAt(n1/2,n2/2,n3/2,1.0f);
    //p.paintAt(   0,   0,   0,1.0f);
    //p.paintAt(n1-1,n2-1,n3-1,1.0f);
    trace("extrapolate ...");
    p.extrapolate();
    trace("done");
    //Array.dump(p.getTimes());
    //plot(p.getValues());
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testConstant();
      }
    });
  }
}
