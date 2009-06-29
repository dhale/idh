/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import static edu.mines.jtk.util.MathPlus.sqrt;

/**
 * Local 2-dimensional diffusion tensors. These tensors are matrices 
 * <pre><code>
 * D = |d11 d12| = (s0*d0)^2 * I + (s1*d1)^2 * v * v',
 *     |d12 d22|
 * </code></pre>
 * where s0 and s1 are constant (non-local) scale factors, d0 is the local 
 * isotropic diffusivity, d1 is the local linear diffusivity, and v is a 
 * unit inline diffusion vector. The symbol I denotes the 2-by-2 identity 
 * matrix, and v * v' denotes the 2-by-2 outer product of a vector v.
 * <p>
 * These diffusion tensors are local in the sense that d0, d1, and v may
 * vary as a function of sample indices (i1,i2).
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.10
 */
public class LocalDiffusionTensors2 {

  /**
   * Constructs tensors for specified diffusivities and inline vectors.
   * Local isotropic and linear diffusivities are assumed equal to one.
   * <p>
   * The array v1 is passed by reference, not by copy.
   * @param s0 constant scale factor for isotropic diffusivity.
   * @param s1 constant scale factor for linear diffusivity.
   * @param v1 array of 1st components of inline unit vectors.
   */
  public LocalDiffusionTensors2(double s0, double s1, float[][] v1) {
    this(s0,s1,null,null,v1);
  }

  /**
   * Constructs tensors for specified diffusivities and inline vectors.
   * Initially sets anisotropy to zero, so that anisotropy depends
   * entirely on the specified diffusivities.
   * <p>
   * All arrays are passed by reference, not by copy.
   * @param s0 constant scale factor for isotropic diffusivity.
   * @param s1 constant scale factor for linear diffusivity.
   * @param d0 array of isotropic diffusivities; if null, assumed to be one.
   * @param d1 array of linear diffusivities; if null, assumed to be one.
   * @param v1 array of 1st components of inline unit vectors.
   */
  public LocalDiffusionTensors2(
    double s0, double s1, float[][] d0, float[][] d1, float[][] v1) {
    _n1 = v1[0].length;
    _n2 = v1.length;
    _s0 = (float)s0;
    _s1 = (float)s1;
    _d0 = d0;
    _d1 = d1;
    _v1 = v1;
  }

  /**
   * Gets tensor elements {d11, d12, d22} for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param d array of tensor elements {d11,d12,d22}.
   */
  public void getTensor(int i1, int i2, float[] d) {
    float d0i = (_d0!=null)?_s0*_d0[i2][i1]:_s0;
    float d1i = (_d1!=null)?_s1*_d1[i2][i1]:_s1;
    float s0i = d0i*d0i;
    float s1i = d1i*d1i;
    float v1i = _v1[i2][i1];
    float v2i = sqrt(1.0f-v1i*v1i);
    float d11 = s0i+s1i*v1i*v1i;
    float d12 =     s1i*v1i*v2i;
    float d22 = s0i+s1i*v2i*v2i;
    d[0] = d11;
    d[1] = d12;
    d[2] = d22;
  }

  /**
   * Gets the number of tensors in the 1st dimension.
   * @return the number of tensors in the 1st dimension.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of tensors in the 2nd dimension.
   * @return the number of tensors in the 2nd dimension.
   */
  public int getN2() {
    return _n2;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2;
  private float _s0,_s1;
  private float[][] _d0,_d1,_v1;
} 
