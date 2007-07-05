/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

/**
 * Interface for 3-D diffusion tensors. A diffusion tensor D is a 3 x 3
 * symmetric positive-semidefinite matrix with six unique elements:
 * <pre><code>
 *      [ D11 D12 D13 ]
 *  D = [ D12 D22 D23 ]
 *      [ D13 D23 D33 ]
 * </code></pre>
 * Because diffusion tensors D are symmetric and positive-semidefinite,
 * their eigenvalues are real and non-negative. We let u, v, and w
 * denote the eigenvectors corresponding to the eigenvalues eu, ev, ew,
 * where eu &gt;= ev &gt;= ew. The eigenvectors are unit vectors with
 * indefinite signs; e.g., both u and -u are eigenvectors of a tensor D.
 * <p>
 * Eigenvectors and eigenvalues of diffusion tensors are often used in 
 * image processing. Large 3-D images may consume large amounts of memory,
 * and their diffusion tensors may consume even more. Therefore, classes 
 * that implement this interface may use sophisticated data structures and 
 * compression techniques to reduce memory requirements.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.26
 */
public interface DiffusionTensor3 {

  /**
   * Gets array of elements {D11, D12, D13, D22, D23, D33}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getElements(int i1, int i2, int i3);

  /**
   * Gets the eigenvector u corresponding to the largest eigenvalue.
   * The array contains the vector components {u1, u2, u3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvectorU(int i1, int i2, int i3);

  /**
   * Gets the eigenvector v corresponding to the second largest eigenvalue.
   * The array contains the vector components {v1, v2, v3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvectorV(int i1, int i2, int i3);

  /**
   * Gets the eigenvector u corresponding to the smallest eigenvalue.
   * The array contains the vector components {w1, w2, w3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvectorW(int i1, int i2, int i3);

  /**
   * Gets the three eigenvalues {eu, ev, ew}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvalues(int i1, int i2, int i3);
}
