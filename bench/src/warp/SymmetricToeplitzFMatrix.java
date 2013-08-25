/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

/**
 * A symmetric Toeplitz matrix is a square matrix specified by one row.
 * Elements of a Toeplitz matrix are Aij = a[i-j]. In other words, all
 * elements on each diagonal of a Toeplitz matrix are equal. This class
 * therefore enables solution of equations of the form
 * <pre><code>
 *  |a[0]    a[1]    a[2]    a[3]| |x[0]|     |b[0]|
 *  |a[1]    a[0]    a[1]    a[2]| |x[1]|  =  |b[1]|
 *  |a[2]    a[1]    a[0]    a[1]| |x[2]|     |b[2]|
 *  |a[3]    a[2]    a[1]    a[0]| |x[3]|     |b[3]|
 * </code></pre>
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.25
 */
public class SymmetricToeplitzFMatrix {

  /**
   * Constructs a symmetric Toeplitz matrix with specified elements.
   * Stores the specified array by reference, not by copy.
   * @param a array of elements for the top row of the matrix.
   */
  public SymmetricToeplitzFMatrix(float[] a) {
    _a = a;
    _t = new float[a.length];
  }

  /**
   * Solves this symmetric Toeplitz system for specified right-hand-side.
   * @param b input array containing the right-hand-side column vector.
   * @return array containing the left-hand-side solution vector.
   */
  public float[] solve(float[] b) {
    float[] x = new float[b.length];
    solve(b,x);
    return x;
  }

  /**
   * Solves this symmetric Toeplitz system for specified right-hand-side.
   * @param b input array containing the right-hand-side column vector.
   * @param x output array containing the left-hand-side solution vector.
   */
  public void solve(float[] b, float[] x) {
    solve(_t,_a,b,x);
  }

  /**
   * Solves a symmetric Toeplitz system for specified right-hand-side.
   * @param a input array of elements for top row of matrix.
   * @param b input array containing the right-hand-side column vector.
   * @param x output array containing the left-hand-side vector of unknowns.
   */
  public static void solve(float[] a, float[] b, float[] x) {
    float[] t = new float[a.length];
    solve(t,a,b,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[] _a; // top row of matrix
  private float[] _t; // work array

  // Solution of Ax = b using a work array t and Levinson recursion.
  private static void solve(float[] t, float[] a, float[] b, float[] x) {
    int n = a.length;
    t[0] = 1.0f;
    float v = a[0];
    x[0] = b[0]/a[0];
    for (int i=1; i<n; ++i) {
      t[i] = 0.0f;
      x[i] = 0.0f;
      float e = 0.0f;
      for (int j=0; j<i; ++j)
        e += t[j]*a[i-j];
      float c = e/v;
      v -= c*e;
      for (int j=0; j<=i/2; ++j) {
        float timj = t[i-j]-c*t[j];
        t[j] -= c*t[i-j];
        t[i-j] = timj;
      }
      float w = 0.0f;
      for (int j=0; j<i; ++j)
        w += x[j]*a[i-j];
      c = (w-b[i])/v;
      for (int j=0; j<=i; ++j)
        x[j] -= c*t[i-j];
    }
  }
}

