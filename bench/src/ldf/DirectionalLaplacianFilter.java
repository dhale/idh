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
 * Directional Laplacian filter. The Laplacian operator appears in 
 * diffusion equations, and this filter is useful in the context of 
 * anisotropic diffusion filtering.
 * <p>
 * This filter computes y = y+G'DGx where G is the gradient operator, 
 * G' is its adjoint, and D is a local diffusion tensor that determines 
 * for each image sample the direction of the Laplacian filter.
 * <p>
 * For example, if D = dvv' for local diffusivities d and unit vectors v,
 * then G'DGx is zero for image features that are constant in the direction 
 * of v. The diffusivities d depend on local scale factors ds multiplied by 
 * a nominal filter half-width sigma. Specifically, for each sample, 
 * diffusivities d = 0.5*(ds*sigma)*(ds*sigma).
 * <p>
 * Alternatively, if D = d(I-UU'), then the right-hand side G'DGx is zero
 * for image features that are constant in all directions orthogonal to the 
 * unit vectors u. 
 * <p>
 * Directional Laplacian filters are rarely used alone. While zeroing some 
 * features in images, they tend to attenuate many other features as well. 
 * Therefore, these filters are typically used in combinations with others.
 * <p>
 * For example, the filter implied by (I+G'DG)y = G'DGx acts as a notch
 * filter. It attenuates features for which G'DGx is zero while preserving 
 * other features. Diffusivities d (inside D) control the width of the notch.
 * Note that application of this filter requires solving a sparse symmetric 
 * positive-definite system of equations.
 * <p>
 * An even simpler example is the filter implied by (I+G'DG)y = x. This
 * filter smooths features in the directions implied by the tensors D.
 * Again, application of this filter requires solving a sparse symmetric 
 * positive-definite system of equations.
 * <p>
 * The accumulation of the filter output in y = y+G'DGx is useful when
 * constructing such combination filters. Given y = 0, this filter 
 * computes y = G'DGx. Given y = x, this filter computes y = (I+G'DG)x.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.05
 */
public class DirectionalLaplacianFilter {

  /**
   * Constructs a directional Laplacian filter with nominal half-width sigma.
   * @param sigma the nominal half-width for this filter.
   */
  public DirectionalLaplacianFilter(double sigma) {
    _sigma = (float)sigma;
  }

  /**
   * Computes y = y+G'DGx, where D = dvv' and G is the gradient operator. 
   * The right-hand-side G'DGx is zero in the direction of the unit vectors v.
   * Diffusivities d are scale factors that multiply the nominal half-width 
   * sigma for this filter.
   * <p>
   * Only components v1 of the inline vectors v are specified; all components 
   * v2 = sqrt(1-v1*v1) are assumed to be non-negative.
   * @param ds scale factors for diffusivity in direction of unit vectors v;
   *  if null, this method uses constant ds = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float dsi = (ds!=null)?_sigma*ds[i2][i1]:_sigma;
        float svi = 0.5f*dsi*dsi;
        float v1i = v1[i2][i1];
        float v2i = sqrt(1.0f-v1i*v1i);
        float d11 = svi*v1i*v1i;
        float d12 = svi*v1i*v2i;
        float d22 = svi*v2i*v2i;
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = d(I-uu') and G is the gradient operator. 
   * The right-hand-side G'DGx is zero in all directions orthogonal to the 
   * unit vectors u. Diffusivities d are scale factors that multiply the 
   * nominal half-width sigma for this filter.
   * <p>
   * Only components u2 of the normal vectors u are specified; all components 
   * u1 = sqrt(1-u2*u2) are assumed to be non-negative.
   * @param ds scale factors for diffusivity in directions orthogonal to
   *  unit vectors u; if null, this method uses constant ds = 1.
   * @param u2 array of 2nd components of normal unit vectors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyNormal(
    float[][] ds, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float dsi = (ds!=null)?_sigma*ds[i2][i1]:_sigma;
        float sui = 0.5f*dsi*dsi;
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float d11 = 1.0f-sui*u1i*u1i;
        float d12 =     -sui*u1i*u2i;
        float d22 = 1.0f-sui*u2i*u2i;
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // nominal filter half-width

  // Computes y = y+G'DGx for one sample.
  private void apply(
   float d11, float d12, float d22,
   int i1, int i2, float[][] x, float[][] y) 
  {
    float x00 = x[i2  ][i1  ];
    float x01 = x[i2  ][i1-1];
    float x10 = x[i2-1][i1  ];
    float x11 = x[i2-1][i1-1];
    float xa = x00-x11;
    float xb = x01-x10;
    float x1 = 0.5f*(xa-xb);
    float x2 = 0.5f*(xa+xb);
    float y1 = d11*x1+d12*x2;
    float y2 = d12*x1+d22*x2;
    float ya = 0.5f*(y1+y2);
    float yb = 0.5f*(y1-y2);
    y[i2  ][i1  ] += ya;
    y[i2  ][i1-1] -= yb;
    y[i2-1][i1  ] += yb;
    y[i2-1][i1-1] -= ya;
  }
} 
