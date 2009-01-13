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
 * @version 2009.01.11
 */
public class PaintingM {

  /**
   * Constructs a painting for the specified structure tensor field.
   * All painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param pt painting tensors.
   */
  public PaintingM(int n1, int n2, Tensors2 pt) {
    _n1 = n1;
    _n2 = n2;
    _pt = pt;
    _dv = 0.0f;
    _tm = new TimeMarker2(n1,n2,pt);
    _tm.setConcurrency(TimeMarker2.Concurrency.PARALLEL);
    _types = new byte[n2][n1];
    _marks = new int[n2][n1];
    _times = new float[n2][n1];
    _valus = new float[n2][n1];
    clearAll();
  }

  /**
   * Sets the tensors used in this painting.
   * @param pt painting tensors.
   */
  public void setTensors(Tensors2 pt) {
    _pt = pt;
    _tm.setTensors(pt);
  }

  /**
   * Sets the default paint value.
   * Default values are used for all clear (not painted) samples.
   * @param value the default value.
   */
  public void setDefaultValue(float value) {
    _dv = value;
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
        if (_types[i2][i1]!=FIXED) {
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
    if (_types[i2][i1]==FIXED)
      clear(i1,i2);
  }

  /**
   * Paints the specified sample with the specified value.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param v value for the painted sample.
   */
  public void paintAt(int i1, int i2, float v) {
    _types[i2][i1] = FIXED;
    _marks[i2][i1] = i1+i2*_n1;
    _times[i2][i1] = 0.0f;
    _valus[i2][i1] = v;
  }

  /**
   * Extrapolates values from all fixed (explicitly painted) samples.
   * After extrapolation, all samples are either fixed or extrapolated.
   */
  public void extrapolate() {

    // First compute times and marks.
    _tm.apply(_times,_marks);

    // Use marks to compute values.
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_types[i2][i1]!=FIXED) {
          int mark = _marks[i2][i1];
          int j1 = mark%_n1;
          int j2 = mark/_n1;
          _types[i2][i1] = EXTRA;
          _valus[i2][i1] = _valus[j2][j1];
        }
      }
    }
  }

  /**
   * Interpolates values from all fixed samples, using extrapolated samples.
   * After interpolation, all samples are either fixed or interpolated.
   */
  public void interpolate() {

    // Prepare for local smoothing filter.
    float c = 0.25f;
    float[][] s = Array.mul(_times,_times);
    float[][] v = new float[_n2][_n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.01,1000);

    // For all values, interpolate by smoothing extrapolated values.
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        v[i2][i1] = _valus[i2][i1];
      }
    }
    lsf.apply(_pt,c,s,v,_valus);
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_types[i2][i1]!=FIXED) {
          _types[i2][i1] = INTER;
        } else {
          _valus[i2][i1] = v[i2][i1];
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
    return _times;
  }

  /**
   * Gets the painted value for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   */
  public float getValue(int i1, int i2) {
    return _valus[i2][i1];
  }

  /**
   * Gets the array of painted values with index zero.
   * Returns default values for any samples not painted.
   * @return array of values; by reference, not by copy.
   */
  public float[][] getValues() {
    return _valus;
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

  private Tensors2 _pt; // painting tensors`
  private int _n1,_n2; // painting dimensions
  private float _dv; // default values for clear samples
  private byte[][] _types; // clear, fixed, extra, or inter
  private int[][] _marks; // marks used to extrapolate values
  private float[][] _times; // times to nearest fixed sample
  private float[][] _valus; // painted values; default for clear samples
  private TimeMarker2 _tm; // computes times and marks

  private void clear(int i1, int i2) {
    _types[i2][i1] = CLEAR;
    _marks[i2][i1] = -1;
    _times[i2][i1] = INFINITY;
    _valus[i2][i1] = _dv;
  }
}
