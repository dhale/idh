/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Laplacian diffusion filter. The Laplacian operator appears in 
 * diffusion equations, and this filter is useful in the context of 
 * anisotropic diffusion filtering.
 * <p>
 * This filter computes y = y+G'DGx where G is the gradient operator, 
 * G' is its adjoint, and D is a local diffusion tensor that determines 
 * for each image sample the filter coefficients.
 * <p>
 * Laplacian diffusion filters are rarely used alone. While zeroing some 
 * features in images, they tend to attenuate many other features as well. 
 * Therefore, these filters are typically used in combinations with others.
 * <p>
 * For example, the filter implied by (I+G'DG)y = G'DGx acts as a notch
 * filter. It attenuates features for which G'DGx is zero while preserving 
 * other features. The diffusion tensors D control the width, orientation,
 * and anisotropy of the spectral notch. Note that application of this filter 
 * requires solving a sparse symmetric positive-definite system of equations.
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
 * @version 2007.11.10
 */
public class LaplacianDiffusionFilter {

  /**
   * Constructs a Laplacian diffusion filter.
   */
  public LaplacianDiffusionFilter() {
    this(1.0f/12.0f);
  }

  /**
   * Constructs a Laplacian diffusion filter with experimental factor.
   * @param rs the experimental factor for tuning gradient filters.
   */
  public LaplacianDiffusionFilter(double rs) {
    _rs = (float)rs;
  }

  ///////////////////////////////////////////////////////////////////////////
  // 2-D

  /**
   * Computes y = y+G'DGx, for specified local diffusion tensors D.
   * @param ldt local diffusion tensors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void apply(LocalDiffusionTensors2 ldt, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] d = new float[3];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        ldt.getTensor(i1,i2,d);
        apply(_rs,d[0],d[1],d[2],i1,i2,x,y);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // 3-D


  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _rs; // experimental parameter for gradient filters

  // Computes y = y+G'DGx for one sample.
  private static void apply(
   float rs, float d11, float d12, float d22,
   int i1, int i2, float[][] x, float[][] y) 
  {
    float dd = rs*(d11+d22);
    float aa = (0.5f*d11-dd);
    float bm = (0.5f*d12-dd);
    float bp = (0.5f*d12+dd);
    float cc = (0.5f*d22-dd);
    float x00 = x[i2  ][i1  ];
    float x01 = x[i2  ][i1-1];
    float x10 = x[i2-1][i1  ];
    float x11 = x[i2-1][i1-1];
    float a0 = aa*(x00-x01);
    float a1 = aa*(x10-x11);
    float b0 = bp*(x00-x11);
    float b1 = bm*(x01-x10);
    float c0 = cc*(x00-x10);
    float c1 = cc*(x01-x11);
    y[i2  ][i1  ] += a0+b0+c0;
    y[i2  ][i1-1] -= a0+b1-c1;
    y[i2-1][i1  ] += a1+b1-c0;
    y[i2-1][i1-1] -= a1+b0+c1;
  }

  // Computes y = y+G'DGx for one sample.
  // This old method was designed for rs = 0.25.
  private static void applyOld(
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

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void testApplyMethods() {
    int n1 = 3;
    int n2 = 3;
    float[][] x = Array.randfloat(n1,n2);
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    float d11 = 0.5f;
    float d12 = 0.5f;
    float d22 = 0.5f;
    float rs = 0.25f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        apply(rs,d11,d12,d22,i1,i2,x,y);
        applyOld(d11,d12,d22,i1,i2,x,z);
      }
    }
    Array.dump(y);
    Array.dump(z);
  }

  public static void main(String[] args) {
    testApplyMethods();
  }
}
