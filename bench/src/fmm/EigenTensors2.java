/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.dsp.Eigen;

/**
 * An array of eigen-decompositions of tensors for 2D image processing.
 * Each tensor is a symmetric positive-semidefinite 2-by-2 matrix:
 * <pre><code>
 * A = |a11 a12|
 *     |a12 a22|
 * </code></pre>
 * Such tensors can be used to parameterize anisotropic image processing.
 * <p>
 * The eigen-decomposition of the matrix A is
 * <pre><code>
 * A = au*u*u' + av*v*v'; au &gt;= av &gt;= 0
 * </code></pre>
 * where u and v are orthogonal unit eigenvectors of A. (The notation 
 * u' denotes the transpose of u.) The outer products of eigenvectors are
 * scaled by the corresponding eigenvalues au and av.
 * <p>
 * The ordering among eigenvalues is easily ensured by the equivalent 
 * representation in terms of a1 = au-av and a2 = av:
 * <pre><code>
 * A = (a1+a2)*u*u' + a2*v*v'
 *     a1*u*u' + a2*(u*u'+v*v')
 *   = a1*u*u' + a2*I
 * </code></pre>
 * With this representation, the non-negative coefficients (a1,a2)
 * are unordered, and the redundancy of the eigenvector v is apparent.
 * <p>
 * An intuitive interpretation of the coefficients (a1,a2) is that a1 
 * corresponds to 1D scaling in the direction of u and a2 corresponds 
 * to 2D isotropic scaling.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.09
 */
public class EigenTensors2 {

  /**
   * Constructs tensors for specified array dimensions. All coefficients 
   * (a1,a2) and eigenvectors u are not set and are initially zero.
   * @param n1 number of tensors in 1st dimension.
   * @param n2 number of tensors in 2nd dimension.
   */
  public EigenTensors2(int n1, int n2) {
    _n1 = n1;
    _n2 = n2;
    _a1 = new float[n2][n1];
    _a2 = new float[n2][n1];
    _u1 = new float[n2][n1];
    _u2 = new float[n2][n1];
  }

  /**
   * Constructs tensors for specified array dimensions and coefficients.
   * The 3rd components of eigenvectors u and v are computed from the 1st 
   * and 2nd components and are assumed to be non-negative.
   * @param u1 array of 1st components of u.
   * @param u2 array of 2nd components of u.
   * @param a1 array of 1D coefficients.
   * @param a2 array of 2D coefficients.
   * @param compressed true, for compressed tensors; false, otherwise.
   */
  public EigenTensors2(
    float[][] u1, float[][] u2,
    float[][] a1, float[][] a2)
  {
    this(u1[0].length,u1.length);
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        float a1i = a1[i2][i1];
        float a2i = a2[i2][i1];
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        setCoefficients(i1,i2,a1i,a2i);
        setEigenvectorU(i1,i2,u1i,u2i);
      }
    }
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

  /**
   * Gets tensor elements for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param a array {a11,a12,a22} of tensor elements.
   */
  public void getTensor(int i1, int i2, float[] a) {
    float a1 = _a1[i2][i1];
    float a2 = _a2[i2][i1];
    float u1 = _u1[i2][i1];
    float u2 = _u2[i2][i1];
    a[0] = a1*u1*u1+a2; // a11
    a[1] = a1*u1*u2   ; // a12
    a[2] = a1*u2*u2+a2; // a22
  }

  /**
   * Gets tensor elements for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @return a array {a11,a12,a13,a22,a23,a33} of tensor elements.
   */
  public float[] getTensor(int i1, int i2) {
    float[] a = new float[3];
    getTensor(i1,i2,a);
    return a;
  }


  /**
   * Gets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param a array {a1,a2} of coefficients.
   */
  public void getCoefficients(int i1, int i2, float[] a) {
    a[0] = _a1[i2][i1];
    a[1] = _a2[i2][i1];
  }

  /**
   * Gets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @return array {a1,a2} of coefficients.
   */
  public float[] getCoefficients(int i1, int i2) {
    float[] a = new float[2];
    getCoefficients(i1,i2,a);
    return a;
  }

  /**
   * Gets the eigenvector u for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param u array {u1,u2} of eigenvector components.
   */
  public void getEigenvectorU(int i1, int i2, float[] u) {
    u[0] = _u1[i2][i1];
    u[1] = _u2[i2][i1];
  }

  /**
   * Gets the eigenvector u for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @return array {u1,u2} of eigenvector components.
   */
  public float[] getEigenvectorU(int i1, int i2) {
    float[] u = new float[2];
    getEigenvectorU(i1,i2,u);
    return u;
  }

  /**
   * Gets the eigenvector v for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param v array {v1,v2} of eigenvector components.
   */
  public void getEigenvectorV(int i1, int i2, float[] v) {
    v[0] =  _u2[i2][i1];
    v[1] = -_u1[i2][i1];
  }

  /**
   * Gets the eigenvector v for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @return array {v1,v2} of eigenvector components.
   */
  public float[] getEigenvectorV(int i1, int i2) {
    float[] v = new float[2];
    getEigenvectorV(i1,i2,v);
    return v;
  }

  /**
   * Sets tensor elements for specified indices.
   * This method first computes an eigen-decomposition of the specified
   * tensor, and then stores the computed eigenvectors and coefficients.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param a array {a11,a12,a22} of tensor elements.
   */
  public void setTensor(int i1, int i2, float[] a) {
    setTensor(i1,i2,a[0],a[1],a[2]);
  }

  /**
   * Sets tensor elements for specified indices.
   * This method first computes an eigen-decomposition of the specified
   * tensor, and then stores the computed eigenvectors and coefficients.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param a11 tensor element a11.
   * @param a12 tensor element a12.
   * @param a22 tensor element a22.
   */
  public void setTensor(int i1, int i2, float a11, float a12, float a22) {
    float[][] aa = {
      {a11,a12},
      {a12,a22}
    };
    float[][] vv = new float[2][2];
    float[] ev = new float[2];
    Eigen.solveSymmetric22(aa,vv,ev);
    float[] u = vv[0];
    float au = ev[0]; if (au<0.0f) au = 0.0f;
    float av = ev[1]; if (av<0.0f) av = 0.0f;
    setEigenvectorU(i1,i2,u);
    setCoefficients(i1,i2,au-av,av);
  }

  /**
   * Sets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param a1 1D coefficient.
   * @param a2 2D coefficient.
   */
  public void setCoefficients(int i1, int i2, float a1, float a2) {
    _a1[i2][i1] = a1;
    _a2[i2][i1] = a2;
  }

  /**
   * Sets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param a array {a1,a2} of coefficients.
   */
  public void setCoefficients(int i1, int i2, float[] a) {
    setCoefficients(i1,i2,a[0],a[1]);
  }

  /**
   * Sets the eigenvector u for the tensor with specified indices.
   * The specified vector is assumed to have length one.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param u1 1st component of u.
   * @param u2 2nd component of u.
   */
  public void setEigenvectorU(int i1, int i2, float u1, float u2) {
    _u1[i2][i1] = u1;
    _u2[i2][i1] = u2;
  }

  /**
   * Sets the eigenvector u for the tensor with specified indices.
   * The specified vector is assumed to have length one.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param array {u1,u2} of eigenvector components.
   */
  public void setEigenvectorU(int i1, int i2, float[] u) {
    setEigenvectorU(i1,i2,u[0],u[1]);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2;
  private float[][] _a1;
  private float[][] _a2;
  private float[][] _u1;
  private float[][] _u2;
}
