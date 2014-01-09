/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.util.MedianFinder;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Running median (and weighted median) filtering of sequences.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.02
 */
public class MedianFilter {

  /**
   * Constructs a median filter with specified half-width m.
   * @param hw the half-width. The total width is 1+2*hw.
   */
  public MedianFilter(int hw) {
    this(null,hw);
  }

  /**
   * Constructs a weighted-median filter with specified weights.
   * A specified index k and the length m of the array of weights determine
   * which input samples are used to compute each output sample. Except near
   * the ends of arrays, each output sample y[i] will be a weighted median of
   * input samples x[i-k:i-k+m-1], with corresponding weights w[0:m-1].
   * @param w array of weights; by copy, not by reference.
   * @param k index in input x[i-k] used with w[0] to compute output y[i].
   */
  public MedianFilter(float[] w, int k) {
    _w = copy(w);
    _k = k;
    _m = (w!=null)?w.length:1+2*k;
    _mf = new MedianFinder(_m);
  }

  /** 
   * Applies this filter.
   * @param x input array.
   * @return output array.
   */
  public float[] apply(float[] x) {
    float[] y = new float[x.length];
    apply(x,y);
    return y;
  }

  /** 
   * Applies this filter.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float[] x, float[] y) {
    int n = x.length;
    int m = _w.length;
    float[] t = new float[m];
    for (int i=0; i<n; ++i) {
      int j0 = i-_k;
      int jm = j0+m;
      int jf = max(j0,0);
      int jl = min(jm,n);
      int jt = 0;
      for (int jx=j0; jx<jf; ++jx,++jt)
        t[jt] = x[0];
      for (int jx=jf; jx<jl; ++jx,++jt)
        t[jt] = x[jx];
      for (int jx=jl; jx<jm; ++jx,++jt)
        t[jt] = x[n-1];
      y[i] = _mf.findMedian(t);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
 
  private float[] _w;
  private int _k;
  private int _m; 
  private MedianFinder _mf;
}
