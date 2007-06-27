
/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Subsampled structure tensors for 3-D images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.26
 */
public class SubsampledStructureTensor3 implements StructureTensor3 {

  public SubsampledStructureTensor3(double sigma, float[][][] x) {
    LocalOrientFilter lof = new LocalOrientFilter(sigma);
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] eu = new float[n3][n2][n1];
    float[][][] ev = new float[n3][n2][n1];
    float[][][] ew = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    lof.apply(x,null,null,
      null,u2,u3,
      null,null,null,
      null,w2,w2,
      eu,ev,ew,
      null,null);
  }
  private void subsample(int ms, float[][][] x, float[][][] y) {
  }

  public SubsampledStructureTensor3(
    float[][][] s11, float[][][] s12, float[][][] s13,
    float[][][] s22, float[][][] s23, float[][][] s33)
  {
  }

  /**
   * Gets array of elements {S11, S12, S13, S22, S23, S33}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getElements(int i1, int i2, int i3) {
    float[] u = getVectorU(i1,i2,i3);
    float[] v = getVectorV(i1,i2,i3);
    float[] w = getVectorW(i1,i2,i3);
    float[] e = getValues(i1,i2,i3);
    float u1 = u[0], u2 = u[1], u3 = u[2];
    float v1 = v[0], v2 = v[1], v3 = v[2];
    float w1 = w[0], w2 = w[1], w3 = w[2];
    float eu = e[0], ev = e[1], ew = e[2];
    float s11 = eu*u1*u1+ev*v1*v1+ew*w1*w1;
    float s12 = eu*u1*u2+ev*v1*v2+ew*w1*w2;
    float s13 = eu*u1*u3+ev*v1*v3+ew*w1*w3;
    float s22 = eu*u2*u2+ev*v2*v2+ew*w2*w2;
    float s23 = eu*u2*u3+ev*v2*v3+ew*w2*w3;
    float s33 = eu*u3*u3+ev*v3*v3+ew*w3*w3;
    return new float[]{s11,s12,s13,s22,s23,s33};
  }

  /**
   * Gets the eigenvector u corresponding to the largest eigenvalue.
   * The array contains the vector components {u1, u2, u3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getVectorU(int i1, int i2, int i3) {
    int j1 = index1(i1);
    int j2 = index2(i2);
    int j3 = index3(i3);
    float u2 = _u2[j3][j2][j1];
    float u3 = _u3[j3][j2][j1];
    float us = u2*u2+u3*u3;
    float u1 = (us<1.0f)?sqrt(1.0f-us):0.0f;
    return new float[]{u1,u2,u3};
  }

  /**
   * Gets the eigenvector v corresponding to the second largest eigenvalue.
   * The array contains the vector components {v1, v2, v3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getVectorV(int i1, int i2, int i3) {
    int j1 = index1(i1);
    int j2 = index2(i2);
    int j3 = index3(i3);
    float u2 = _u2[j3][j2][j1];
    float u3 = _u3[j3][j2][j1];
    float us = u2*u2+u3*u3;
    float u1 = (us<1.0f)?sqrt(1.0f-us):0.0f;
    float w2 = _w2[j3][j2][j1];
    float w3 = _w3[j3][j2][j1];
    float ws = w2*w2+w3*w3;
    float w1 = (ws<1.0f)?sqrt(1.0f-ws):0.0f;
    float v1 = u2*w3-u3*w2;
    float v2 = u3*w1-u1*w3;
    float v3 = u1*w2-u2*w1;
    return new float[]{v1,v2,v3};
  }

  /**
   * Gets the eigenvector u corresponding to the smallest eigenvalue.
   * The array contains the vector components {w1, w2, w3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getVectorW(int i1, int i2, int i3) {
    int j1 = index1(i1);
    int j2 = index2(i2);
    int j3 = index3(i3);
    float w2 = _w2[j3][j2][j1];
    float w3 = _w3[j3][j2][j1];
    float ws = w2*w2+w3*w3;
    float w1 = (ws<1.0f)?sqrt(1.0f-ws):0.0f;
    return new float[]{w1,w2,w3};
  }

  /**
   * Gets the three eigenvalues {eu, ev, ew}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getValues(int i1, int i2, int i3) {
    int j1 = index1(i1);
    int j2 = index2(i2);
    int j3 = index3(i3);
    float eu = _eu[j3][j2][j1];
    float ev = _ev[j3][j2][j1];
    float ew = _ew[j3][j2][j1];
    return new float[]{eu,ev,ew};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _ms;
  private int _m1,_m2,_m3;
  private int _n1,_n2,_n3;
  private float[][][] _eu,_ev,_ew;
  private float[][][] _u2,_u3;
  private float[][][] _w2,_w3;

  private int index1(int i1) {
    int j1 = i1/_ms;
    return (j1<_m1)?j1:_m1-1;
  }
  private int index2(int i2) {
    int j2 = i2/_ms;
    return (j2<_m2)?j2:_m2-1;
  }
  private int index3(int i3) {
    int j3 = i3/_ms;
    return (j3<_m3)?j3:_m3-1;
  }
}
