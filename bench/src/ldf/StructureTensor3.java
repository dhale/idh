/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

/**
 * Interface for 3-D structure tensors. A structure tensor S is a 3 x 3
 * symmetric positive-semidefinite matrix of averaged products of components
 * of gradients. (This matrix is also called the gradient-squared tensor.)
 * <pre><code>
 *      [ S11 S12 S13 ]
 *  S = [ S12 S22 S23 ]
 *      [ S13 S23 S33 ]
 * </code></pre>
 * For example, the element S12 is an average of products of g1 and g2,
 * where g1 and g2 denote the first and second components of gradient
 * vectors. In three dimensions, gradients have three components g1, g2,
 * and g3, and so each structure tensor has three unique elements.
 * <p>
 * Because the structure tensor S is symmetric positive-semidefinite,
 * it's eigenvalues are real and non-negative. We let u, v, and w
 * denote the eigenvectors corresponding to the eigenvalues eu, ev, ew,
 * where eu &gt;= ev &gt;= ew. The eigenvectors are unit vectors with
 * indefinite signs; e.g., both u and -u are eigenvectors of a tensor S.
 * <p>
 * Eigenvectors and eigenvalues of structure tensors are often used in 
 * image processing. Large 3-D images may consume large amounts of memory,
 * and their structure tensors may consume even more. Therefore, classes 
 * that implement this interface may use sophisticated data structures and 
 * compression techniques to reduce memory requirements.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.26
 */
public interface StructureTensor3 {

  /**
   * Gets array of elements {S11, S12, S13, S22, S23, S33}.
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
  public float[] getVectorU(int i1, int i2, int i3);

  /**
   * Gets the eigenvector v corresponding to the second largest eigenvalue.
   * The array contains the vector components {v1, v2, v3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getVectorV(int i1, int i2, int i3);

  /**
   * Gets the eigenvector u corresponding to the smallest eigenvalue.
   * The array contains the vector components {w1, w2, w3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getVectorW(int i1, int i2, int i3);

  /**
   * Gets the three eigenvalues {eu, ev, ew}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getValues(int i1, int i2, int i3);
}
