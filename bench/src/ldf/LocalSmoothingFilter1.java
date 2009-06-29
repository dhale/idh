/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import static edu.mines.jtk.util.MathPlus.sqrt;

/**
 * Local anisotropic smoothing filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.10.11
 */
public class LocalSmoothingFilter1 {

  /**
   * Constructs a local smoothing filter.
   * @param sigma the maximum filter half-width.
   * @param n1 the number of samples in sequences to be filtered.
   * @param ds array of diffusivity scale factors
   */
  public LocalSmoothingFilter1(double sigma, int n1, float[] ds) {
    _sigma = (float)sigma;
    _n1 = n1;
    init(sigma,n1,ds);
  }

  /**
   * Applies this filter.
   * The input and output arrays may be the same arrays.
   * @param x array of input samples.
   * @param y array of output samples.
   */
  public void apply(float[] x, float[] y) {

    // y = D * x.
    for (int i1=0; i1<_n1; ++i1) {
      y[i1] = x[i1]*_d[i1];
    }

    // y = inv(L') * D * x.
    for (int i1=_n1-1; i1>0; --i1) {
      y[i1-1] -= _e[i1]*y[i1];
    }
  }

  /**
   * Applies the transpose of this filter.
   * The input and output arrays may be the same arrays.
   * @param x array of input samples.
   * @param y array of output samples.
   */
  public void applyTranspose(float[] x, float[] y) {

    // y = inv(L) * x.
    y[0] = x[0];
    for (int i1=1; i1<_n1; ++i1) {
      y[i1] = x[i1]-_e[i1]*y[i1-1];
    }

    // y = D * inv(L) * x.
    for (int i1=0; i1<_n1; ++i1) {
      y[i1] *= _d[i1];
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // maximum half-width of filter
  private int _n1; // number of samples in 1st dimension
  private float[] _d; // diagonal of matrix D in A = L*inv(D*D)*L'
  private float[] _e; // sub-diagonal of matrix L in A = L*inv(D*D)*L'

  private void init(double sigma, int n1, float[] ds) {
    _sigma = (float)sigma;
    _n1 = n1;

    // Compute sub-diagonal of SPD tridiagonal matrix A in array e.
    _e = new float[n1];
    float ss = 0.50f*_sigma*_sigma;
    for (int i1=1; i1<n1; ++i1)
      _e[i1] = -ss;
    if (ds!=null) {
      for (int i1=1; i1<n1; ++i1) {
        float dsi = 0.50f*(ds[i1-1]+ds[i1]);
        _e[i1] *= dsi*dsi;
      }
    }

    // Compute diagonal of SPD tridiagonal matrix A in array d.
    _d = new float[n1];
    _d[0] = 1.0f-_e[1];
    for (int i1=1; i1<n1-1; ++i1) {
      _d[i1] = 1.0f-_e[i1]-_e[i1+1];
    }
    _d[n1-1] = 1.0f-_e[n1-1];

    // A = L*inv(D*D)*L', where L is lower unit bidiagonal, and D is diagonal.
    // Lower sub-diagonal of L goes in array e; diagonal of D goes in array d.
    _d[0] = 1.0f/_d[0];
    for (int i1=1; i1<n1; ++i1) {
      float t = _e[i1];
      _e[i1] = t*_d[i1-1];
      _d[i1] = 1.0f/(_d[i1]-t*_e[i1]);
    }
    for (int i1=0; i1<n1; ++i1)
      _d[i1] = sqrt(_d[i1]);
  }
}
