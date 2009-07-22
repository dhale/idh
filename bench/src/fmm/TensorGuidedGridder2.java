/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.dsp.LocalSmoothingFilter;
import edu.mines.jtk.dsp.Tensors2;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Tensor-guided gridding of 2D data.
 * Gridding is interpolation of a set of known sample values on a 
 * uniformly sampled grid. Here the interpolation is performed by 
 * a two-step process described by Hale (2009).
 * <p>
 * The first step is to compute for all samples the distance to the 
 * nearest known sample and the value of that known sample. This first
 * step produces a distance map and a nearest-neighbor interpolant.
 * The second step is to blend (smooth) the nearest-neighbor interpolant,
 * where the extent of smoothing varies spatially and is proportional to 
 * distances in the distance map.
 * <p>
 * In tensor-guided gridding, we replace distance with time. Time is a
 * simple term for non-Euclidean distance measured in a metric-tensor
 * field. So "nearest" now means nearest in time. In the first step we 
 * compute a time map by solving an eikonal equation with coefficients 
 * that may be both anisotropic and spatially varying. In the second 
 * step, we blend the nearest-neighbor interpolant with an anisotropic 
 * and spatially varying smoothing filter.
 * <p>
 * The default tensor field is homogeneous and isotropic. In this
 * special case, time is equivalent to distance, and tensor-guided
 * gridding is similar to gridding with Sibson's natural neighbor 
 * interpolant.
 * <p>
 * Reference: 
 * <a href="http://www.mines.edu/papers/Hale09ImageGuidedBlendedNeighborInterpolation.pdf">
 * Hale, D., 2009, Image-guided blended neighbor interpolation, CWP-634</a>
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.07.21
 */
public class TensorGuidedGridder2 {

  /**
   * Constructs a gridder for a homogeneous and isotropic tensor field.
   */
  public TensorGuidedGridder2() {
    this(null);
  }

  /**
   * Constructs a gridder for the specified tensors.
   * @param tensors the tensors.
   */
  public TensorGuidedGridder2(Tensors2 tensors) {
    setTensors(tensors);
  }

  /**
   * Sets the tensor field used by this gridder.
   * The default is a homogeneous and isotropic tensor field.
   * @param tensors the tensors; null for default tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = tensors;
    if (_tensors==null) {
      _tensors = new Tensors2() {
        public void getTensor(int i1, int i2, float[] d) {
          d[0] = 1.0f;
          d[1] = 0.0f;
          d[2] = 1.0f;
        }
      };
    }
  }

  /**
   * Computes gridded values using nearest neighbors.
   * Gridded values in the array p are computed for only unknown 
   * samples denoted by corresponding non-zero times in the array t. 
   * This method does not change known values in p, which correspond
   * to zero times in t.
   * @param t array of times to nearest known samples.
   * @param p array of nearest-neighbor gridded values.
   */
  public void gridNearest(float[][] t, float[][] p) {
    int n1 = t[0].length;
    int n2 = t.length;

    // First count the known samples, the number of marks.
    int nmark = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]==0.0f)
          ++nmark;
      }
    }

    // Make an array for marks, while storing values of known samples
    // in an array of values indexed by mark.
    float[] pmark = new float[nmark];
    int[][] m = new int[n2][n1];
    int mark = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]==0.0f) {
          pmark[mark] = p[i2][i1];
          m[i2][i1] = mark;
          ++mark;
        }
      }
    }

    // Use the time marker to compute both times and marks.
    TimeMarker2 tm = new TimeMarker2(n1,n2,_tensors);
    tm.apply(t,m);

    // Use the marks to compute the nearest-neighbor interpolant.
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]!=0.0f)
          p[i2][i1] = pmark[m[i2][i1]];
      }
    }
  }

  /**
   * Computes gridded values using blended neighbors. 
   * Note that blended-neighbor gridding can be performed only 
   * after nearest-neighbor gridding. Blending does not change
   * the values of known samples for which times are zero.
   * @param t array of times to nearest known samples.
   * @param p array of nearest-neighbor gridded values.
   * @param q array of blended-neighbor gridded values.
   */
  public void gridBlended(float[][] t, float[][] p, float[][] q) {
    int n1 = t[0].length;
    int n2 = t.length;

    // Compute time squared, shifted to account for the shift in the
    // finite-difference stencil usd in the local smoothing filter.
    float[][] s = mul(t,t);
    for (int i2=n2-1; i2>0; --i2) {
      for (int i1=n1-1; i1>0; --i1) {
        s[i2][i1] = 0.25f*(s[i2  ][i1  ] +
                           s[i2  ][i1-1] +
                           s[i2-1][i1  ] +
                           s[i2-1][i1-1]);
      }
    }

    // Construct an apply a local smoothing filter.
    float c = 0.5f; // constant for linear precision
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.01,10000);
    lsf.apply(_tensors,c,s,p,q);

    // Restore the known sample values. Due to errors in finite-difference
    // approximations, these values may have changed during smoothing.
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]==0.0f) {
          q[i2][i1] = p[i2][i1];
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Tensors2 _tensors;
}
