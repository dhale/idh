/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * 1D interpolation using fast marching methods.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.23
 */
public class ImagePainter1 {

  public ImagePainter1(float[] d) {
    _d = Array.copy(d);
    _nv = 1;
    _n1 = d.length;
    _dv = new float[_nv];
    _k1 = new int[_n1];
    _tk = new float[_n1];
    _vk = new float[_n1][];
    _type = new byte[_n1];
    _mark = new int[_n1];
    _hmin = new TimeHeap1(TimeHeap1.Type.MIN,_n1);
    _hmax = new TimeHeap1(TimeHeap1.Type.MAX,_n1);
    clearAll();
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
   * Sets normalization used in interpolation. With normalization,
   * interpolation weights are computed from distance divided by
   * velocity.
   * @param normalize true, for normalization; false, otherwise.
   */
  public void setNormalize(boolean normalize) {
    _normalize = normalize;
  }

  /**
   * Clears all painted values, including all fixed values.
   */
  public void clearAll() {
    for (int i1=0; i1<_n1; ++i1) {
      clear(i1);
    }
  }

  /**
   * Clears painted values that are not fixed.
   */
  public void clearNotFixed() {
    for (int i1=0; i1<_n1; ++i1) {
      if (_type[i1]!=FIXED) {
        clear(i1);
      }
    }
  }

  /**
   * Erases values for any fixed sample with specified indices.
   * If the specified sample is not fixed, this method does nothing.
   * @param i1 index in 1st dimension of sample to erase.
   */
  public void eraseFixedAt(int i1) {
    if (_type[i1]==FIXED)
      clear(i1);
  }
  private void clear(int i1) {
    _type[i1] = CLEAR;
    _k1[i1] = -1;
    _tk[i1] = INFINITY;
    _vk[i1] = null;
  }

  /**
   * Paints the specified sample with one specified value at index zero.
   * Paints default values for indices greater than zero.
   * @param i1 index in 1st dimension of sample to paint.
   * @param v value at index zero for the painted sample.
   */
  public void paintAt(int i1, float v) {
    _type[i1] = FIXED;
    _k1[i1] = i1;
    _tk[i1] = INFINITY;
    _vk[i1] = new float[_nv];
    _vk[i1][0] = v;
    for (int iv=1; iv<_nv; ++iv)
      _vk[i1][iv] = _dv[iv];
  }

  /**
   * Paints the specified sample with specified values.
   * After painting, the specified sample is fixed; its values will
   * not change in any subsequent extrapolation or interpolation.
   * @param i1 index in 1st dimension of sample to paint.
   * @param v array of values for painted sample; by copy, not by reference.
   */
  public void paintAt(int i1, float[] v) {
    _type[i1] = FIXED;
    _k1[i1] = i1;
    _tk[i1] = INFINITY;
    _vk[i1] = Array.copy(v);
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
    for (int i1=0; i1<_n1; ++i1) {
      _tk[i1] = INFINITY;
      if (_type[i1]==FIXED) {
        _hmax.insert(i1,INFINITY);
      } else {
        _type[i1] = CLEAR;
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
      TimeHeap1.Entry ef = _hmax.remove();
      int k1 = ef.i1;

      // The values to be extrapolated.
      float[] vk = _vk[k1];

      // Mark all samples as far, mark the fixed sample as known with 
      // time zero, and update its neighbors.
      clearMarks();
      _hmin.clear();
      _mark[k1] = _known;
      _k1[k1] = k1;
      _tk[k1] = 0.0f;
      updateNabors(k1);

      // Extrapolate from the fixed sample to all samples that are
      // nearer to the fixed sample than to any other fixed sample.
      while (!_hmin.isEmpty()) {
        TimeHeap1.Entry e = _hmin.remove();
        int i1 = e.i1;
        float t = e.t;
        _mark[i1] = _known;
        if (_type[i1]==FIXED) {
          _hmax.reduce(i1,t);
        } else {
          _type[i1] = EXTRA;
          _k1[i1] = k1;
          _vk[i1] = vk;
        }
        updateNabors(i1);
      }
    }
  }

  /**
   * Interpolates values from all fixed samples, using extrapolated samples.
   * After interpolation, all samples are either fixed or interpolated.
   */
  public void interpolate() {

    // Insert all extrapolated samples into the max-heap with their
    // current times. After the max-heap is built, the extrapolated
    // sample with largest time is at the top of the heap.
    _hmax.clear();
    for (int i1=0; i1<_n1; ++i1) {
      if (_type[i1]==EXTRA) {
        _hmax.insert(i1,_tk[i1]);
      }
    }

    // Interpolate all extrapolated (not-fixed) samples, one at a time, 
    // in order of decreasing time. The extrapolated sample with the 
    // largest time is at the top of the max-heap. As we interpolate 
    // this sample, times for other extrapolated samples that remain 
    // in the heap may be reduced, so that their order may change. 
    while (!_hmax.isEmpty()) {

      // Remove from the max-heap the extrapolated sample with largest time.
      // This is the sample to be interpolated; the "interpolated sample".
      TimeHeap1.Entry te = _hmax.remove();
      int k1 = te.i1;
      float tk = te.t;

      // Sum of weights in weighted sum of values.
      float wi = 1.0f/_d[k1];
      if (_normalize)
        wi /= _d[k1];
      float ws = wi;

      // The values to be interpolated. This array will be assigned to 
      // all extrapolated samples during the march away from the 
      // interpolated sample. This assignment is one reason that an 
      // array of values is so useful, for we will not actually know the 
      // interpolated values until the march is complete. Later, when we 
      // have completed the computation of the interpolated values, those 
      // values will already be referenced by all extrapolated samples 
      // nearest to the interpolated sample.
      float[] vk = Array.mul(wi,_vk[k1]);
      _vk[k1] = vk;

      // Mark all samples as far, set the type of the interpolated sample,
      // mark the interpolated sample as known with time zero, and update 
      // its neighbors.
      clearMarks();
      _hmin.clear();
      _type[k1] = INTER;
      _mark[k1] = _known;
      _k1[k1] = k1;
      _tk[k1] = 0.0f;
      updateNabors(k1);

      // March away from the interpolated sample to all extrapolated
      // samples that are nearer to the interpolated sample than to any 
      // fixed samples or samples previously interpolated.
      // While marching, accumulate values needed for interpolation.
      while (!_hmin.isEmpty()) {

        // Get the extrapolated sample with minimum time.
        TimeHeap1.Entry e = _hmin.remove();
        int i1 = e.i1;
        float ti = e.t;

        // Weighted sum of values for the extrapolated sample.
        float[] vki = _vk[i1];
        wi = 1.0f;
        if (_normalize)
          wi /= _d[i1];
        ws += wi;
        for (int iv=0; iv<_nv; ++iv)
          vk[iv] += wi*vki[iv];

        // Mark the extrapolated sample known. Reduce it's time in the 
        // max-heap. Its values will be those of the interpolated sample.
        // Continue marching by updating the nabor samples.
        _mark[i1] = _known;
        _hmax.reduce(i1,ti);
        _k1[i1] = k1;
        _vk[i1] = vk;
        updateNabors(i1);
      }

      // The march is complete, and we have a weighted sum of values.
      // Now simply divide by sum of the weights.
      for (int iv=0; iv<_nv; ++iv)
        vk[iv] /= ws;
    }
  }

  /**
   * Gets the array of times used in this painting.
   * These times are modified by extrapolation and interpolation.
   * @return array of times; by reference, not by copy.
   */
  public float[] getTimes() {
    return _tk;
  }

  /**
   * Gets the painted value with index zero for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   */
  public float getValue(int i1) {
    return getValue(i1,0);
  }

  /**
   * Gets the painted value for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param iv index of value to get.
   */
  public float getValue(int i1, int iv) {
    float[] vk = _vk[i1];
    return (vk!=null)?vk[iv]:_dv[iv];
  }

  /**
   * Gets a copy of the painted values with index zero.
   * Returns default values for any samples not painted.
   * @return array of values; by copy, not by reference.
   */
  public float[] getValues() {
    return getValues(0);
  }

  /**
   * Gets a copy of the painted values with specified index.
   * Returns default values for any samples not painted.
   * @param iv index of values to get.
   * @return array of values; by copy, not by reference.
   */
  public float[] getValues(int iv) {
    float[] v = new float[_n1];
    for (int i1=0; i1<_n1; ++i1) {
      v[i1] = getValue(i1,iv);
    }
    return v;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // The value for times not yet computed. Also the value returned by
  // methods that compute times from nabor times when a valid time
  // cannot be computed. We use the maximum possible float so that
  // it will be larger than any valid times we compute.
  private static final float INFINITY = Float.MAX_VALUE;

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
      for (int i1=0; i1<_n1; ++i1) {
        _mark[i1] = _far;
      }
    } else {
      _far += 2; // all known samples instantly become far samples
      _trial +=2; // no samples are trial
      _known +=2; // no samples are known
    }
  }

  private int _n1; // painting dimensions
  private int _nv; // number of values associated with each sample
  private float[] _d; // diffusion coefficients (like velocities)
  private float[] _dv; // default values for clear samples
  private float[][] _vk; // painted values; null for clear samples
  private float[] _tk; // time to nearest fixed or interpolated sample
  private int[] _k1; // indices of nearest fixed or interp sample
  private int[] _mark; // samples are marked far, trial, or known
  private byte[] _type; // fixed, extra, 
  private TimeHeap1 _hmin; // the min heap
  private TimeHeap1 _hmax; // the max heap
  private boolean _normalize; // true, for normalization by velocity.

  private void updateNabors(int i1) {
    for (int k1=-1; k1<=1; k1+=2) {
      int j1 = i1+k1;
      if (j1<0 || j1>=_n1) continue;
      if (_mark[j1]!=_known) {
        float tj = _tk[i1]+1.0f/_d[j1];
        if (tj<_tk[j1]) {
          if (_mark[j1]!=_trial) {
            _mark[j1] = _trial;
            _hmin.insert(j1,tj);
          } else {
            _hmin.reduce(j1,tj);
          }
          _tk[j1] = tj;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void trace(String s) {
    System.out.println(s);
  }

  private static void plotSequence(Sampling sx, float[] f) {
    SimplePlot sp = SimplePlot.asSequence(sx,f);
    double xmin = sx.getFirst();
    double xmax = sx.getLast();
    double fmin = Array.min(f);
    double fmax = Array.max(f);
    double xpad = 0.05*(xmax-xmin);
    double fpad = 0.05*(fmax-fmin);
    sp.setHLimits(xmin-xpad,xmax+xpad);
    sp.setVLimits(fmin-fpad,fmax+fpad);
  }

  private static void plotPoints(Sampling sx, float[] f, String file) {
    double xmin = sx.getFirst();
    double xmax = sx.getLast();
    double fmin = Array.min(f);
    double fmax = Array.max(f);
    double xpad = 0.05*(xmax-xmin);
    double fpad = 0.05*(fmax-fmin);
    SimplePlot sp = new SimplePlot();
    sp.setVisible(false);
    sp.setHLimits(xmin-xpad,xmax+xpad);
    sp.setVLimits(fmin-fpad,fmax+fpad);
    sp.setSize(700,500);
    sp.addPoints(sx,f);
    sp.setVisible(true);
    if (file!=null)
      sp.paintToPng(300,6,"png/"+file+".png");
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        int n1 = 4001;
        float d1 = 1.0f/(float)(n1-1);
        //float d1 = 1.0f;
        float f1 = 0.0f;
        Sampling s1 = new Sampling(n1,d1,f1);
        float[] d = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          d[i1] = i1<4*(n1-1)/10?1.0f:0.5f;
          //float x1 = f1+i1*d1;
          //d[i1] = pow((1.00001f-x1*x1*(3.0f-2.0f*x1)),0.25f);
        }
        plotPoints(s1,d,"ip1v");
        ImagePainter1 ip = new ImagePainter1(d);
        ip.setNormalize(true);
        ip.paintAt(   0,0.0f);
        ip.paintAt(n1-1,1.0f);
        plotPoints(s1,ip.getValues(),null);
        ip.extrapolate();
        plotPoints(s1,ip.getValues(),"ip1fe");
        plotPoints(s1,Array.mul(d1,ip.getTimes()),"ip1t");
        ip.interpolate();
        plotPoints(s1,ip.getValues(),"ip1fi");
      }
    });
  }
}
