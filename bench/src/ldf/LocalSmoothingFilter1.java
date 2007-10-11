/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic smoothing filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.21
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

  public void apply(float[] x, float[] y) {

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

  public void applyTranspose(float[] x, float[] y) {

    // y = D * inv(L) * x.
    for (int i1=0; i1<_n1; ++i1) {
      y[i1] = x[i1]*_d[i1];
    }

    // y = inv(L') * D * inv(L) * x.
    for (int i1=_n1-1; i1>0; --i1) {
      y[i1-1] -= _e[i1]*y[i1];
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma;
  private int _n1;
  private float[] _d,_e;

  private void init(double sigma, int n1, float[] ds) {
    _sigma = (float)sigma;
    _n1 = n1;

    // Sub-diagonal of SPD tridiagonal matrix A in array e.
    float[] e = new float[n1];
    float ss = 0.50f*_sigma*_sigma;
    for (int i1=1; i1<n1; ++i1)
      e[i1] = -ss;
    if (ds!=null) {
      for (int i1=1; i1<n1; ++i1)
        e[i1] *= 0.50f*(ds[i1-1]+ds[i1]);
    }

    // Diagonal of SPD tridiagonal matrix A in array d.
    float[] d = new float[n1];
    d[0] = 1.0f-e[1];
    for (int i1=1; i1<n1-1; ++i1) {
      d[i1] = 1.0f-e[i1]-e[i1+1];
    }
    d[n1-1] = 1.0f-e[n1-1];

    // A = L*inv(D*D)*L', where L is lower unit bidiagonal, and D is diagonal.
    // Lower sub-diagonal of L goes in array e; diagonal of D in array d.
    d[0] = 1.0f/d[0];
    for (int i1=1; i1<n1; ++i1) {
      float t = e[i1];
      e[i1] = t*d[i1-1];
      d[i1] = 1.0f/(d[i1]-t*e[i1]);
    }
    for (int i1=0; i1<n1; ++i1)
      d[i1] = sqrt(d[i1]);
    
    // Save the decomposition used to apply this filter or its transpose.
    _d = d;
    _e = e;
  }
}
