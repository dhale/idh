/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;

/**
 * A 3-D array of structure tensors. Each structure tensor is a symmetric 
 * positive-semidefinite 3-by-3 matrix:
 * <pre><code>
 *     |s11 s12 s13|
 * S = |s12 s22 s23| = eu*u*u' + ev*v*v' + ew*w*w'
 *     |s13 s23 s33|
 * </code></pre>
 * where eu, ev, and ew are non-negative eigenvalues such that
 * eu %gt;= ev &gt;= ew and u, v, and w are corresponding orthogonal
 * unit eigenvectors. (The notation u' denotes the transpose of the
 * vector u.) In such an eigen-decomposition, the eigenvalues eu, ev, 
 * and ew serve as weights for the 3-by-3 outer outer products of the
 * corresponding eigenvectors u, v, and w.
 * <p>
 * Eigen-decompositions of structure tensors are useful because the
 * eigenvectors represent the orientations of image features, while 
 * the eigenvalues contain information related to the dimensionality 
 * of those features. For example, planar (2-D) features have large 
 * eu and small ev and ew.
 * <p>
 * Computation of eigenvalues and eigenvectors is straightforward
 * but can be costly if performed repeatedly as needed. And caching
 * of eigenvalues and eigenvectors computed once can require more
 * memory than the six 3-D arrays required to store tensor elements.
 * <p> 
 * To balance computational efficiency with memory requirements, the
 * eigen-decomposition of each structure tensor is computed once and 
 * stored approximately in only 10 bytes per tensor. Eigenvalues eu 
 * are stored in single-precision floats. Errors in the approximations 
 * of smaller eigenvalues ev and ew are less than one percent of eu.
 * Angles between approximate eigenvectors and exact eigenvectors u, v, 
 * and w are less than one degree.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.02.12
 */
public class StructureTensors3 {

  /**
   * Constructs tensors for specified array dimensions.
   * @param n1 number of tensors in 1st dimension.
   * @param n2 number of tensors in 2nd dimension.
   * @param n3 number of tensors in 3rd dimension.
   */
  public StructureTensors3(int n1, int n2, int n3) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _eu = new float[n3][n2][n1];
    _bv = new byte[n3][n2][n1];
    _bw = new byte[n3][n2][n1];
    _iu = new short[n3][n2][n1];
    _iw = new short[n3][n2][n1];
  }

  /**
   * Constructs tensors with specified elements.
   * @param s array[6][n3][n2][n1] of tensor elements {s11,s12,s13,s22,s23,s33}.
   */
  public StructureTensors3(float[][][][] s) {
    _n1 = s[0][0][0].length;
    _n2 = s[0][0].length;
    _n3 = s[0].length;
    float[][][] s11 = s[0];
    float[][][] s12 = s[1];
    float[][][] s13 = s[2];
    float[][][] s22 = s[3];
    float[][][] s23 = s[4];
    float[][][] s33 = s[5];
    float[] si = new float[6];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          si[0] = s11[i3][i2][i1];
          si[1] = s12[i3][i2][i1];
          si[2] = s13[i3][i2][i1];
          si[3] = s22[i3][i2][i1];
          si[4] = s23[i3][i2][i1];
          si[5] = s33[i3][i2][i1];
          setTensor(i1,i2,i3,si);
        }
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
   * Gets the number of tensors in the 3rd dimension.
   * @return the number of tensors in the 3rd dimension.
   */
  public int getN3() {
    return _n3;
  }

  /**
   * Gets tensor elements {s11,s12,s13,s22,s23,s33} for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param s array {s11,s12,s13,s22,s23,s33} of tensor elements.
   */
  public void getTensor(int i1, int i2, int i3, float[] s) {
    float bvi = _bv[i3][i2][i1];
    float bwi = _bw[i3][i2][i1];
    if (bvi<0.0) bvi += 256.0f;
    if (bwi<0.0) bwi += 256.0f;
    float eui = _eu[i3][i2][i1];
    float esi = eui*ES_GET;
    float evi = bvi*esi;
    float ewi = bwi*esi;
    float[] u = _uss.getPoint(_iu[i3][i2][i1]);
    float[] w = _uss.getPoint(_iw[i3][i2][i1]);
    float u1i = u[0];
    float u2i = u[1];
    float u3i = u[2];
    float w1i = w[0];
    float w2i = w[1];
    float w3i = w[2];
    float v1i = w2i*u3i-w3i*u2i; // v = w cross u
    float v2i = w3i*u1i-w1i*u3i;
    float v3i = w1i*u2i-w2i*u1i;
    float eu1i = eui*u1i;
    float eu2i = eui*u2i;
    float eu3i = eui*u3i;
    float ev1i = evi*v1i;
    float ev2i = evi*v2i;
    float ev3i = evi*v3i;
    float ew1i = ewi*w1i;
    float ew2i = ewi*w2i;
    float ew3i = ewi*w3i;
    s[0] = eu1i*u1i+ev1i*v1i+ew1i*w1i; // s11
    s[1] = eu1i*u2i+ev1i*v2i+ew1i*w2i; // s12
    s[2] = eu1i*u3i+ev1i*v3i+ew1i*w3i; // s13
    s[3] = eu2i*u2i+ev2i*v2i+ew2i*w2i; // s22
    s[4] = eu2i*u3i+ev2i*v3i+ew2i*w3i; // s23
    s[5] = eu3i*u3i+ev3i*v3i+ew3i*w3i; // s33
  }

  /**
   * Gets the eigenvalues for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param e array {eu,ev,ew} of eigenvalues.
   */
  public void getEigenvalues(int i1, int i2, int i3, float[] e) {
    float bvi = _bv[i3][i2][i1];
    float bwi = _bw[i3][i2][i1];
    if (bvi<0.0) bvi += 256.0f;
    if (bwi<0.0) bwi += 256.0f;
    float eui = _eu[i3][i2][i1];
    float esi = eui*ES_GET;
    e[0] = eui;
    e[1] = bvi*esi;
    e[2] = bwi*esi;
  }

  /**
   * Gets the eigenvalues for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @return array {eu,ev,ew} of eigenvalues.
   */
  public float[] getEigenvalues(int i1, int i2, int i3) {
    float[] e = new float[3];
    getEigenvalues(i1,i2,i3,e);
    return e;
  }

  /**
   * Gets the eigenvector u for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param u array {u1,u2,u3} of eigenvector components.
   */
  public void getEigenvectorU(int i1, int i2, int i3, float[] u) {
    float[] ui = _uss.getPoint(_iu[i3][i2][i1]);
    u[0] = ui[0];
    u[1] = ui[1];
    u[2] = ui[2];
  }

  /**
   * Gets the eigenvector u for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @return array {u1,u2,u3} of eigenvector components.
   */
  public float[] getEigenvectorU(int i1, int i2, int i3) {
    float[] u = new float[3];
    getEigenvectorU(i1,i2,i3,u);
    return u;
  }

  /**
   * Gets the eigenvector v for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param v array {v1,v2,v3} of eigenvector components.
   */
  public void getEigenvectorV(int i1, int i2, int i3, float[] v) {
    float[] u = getEigenvectorU(i1,i2,i3);
    float[] w = getEigenvectorW(i1,i2,i3);
    v[0] = w[1]*u[2]-w[2]*u[1]; // v = w cross u
    v[1] = w[2]*u[0]-w[0]*u[2];
    v[2] = w[0]*u[1]-w[1]*u[1];
  }

  /**
   * Gets the eigenvector v for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @return array {v1,v2,v3} of eigenvector components.
   */
  public float[] getEigenvectorV(int i1, int i2, int i3) {
    float[] v = new float[3];
    getEigenvectorV(i1,i2,i3,v);
    return v;
  }

  /**
   * Gets the eigenvector w for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param w array {w1,w2,w3} of eigenvector components.
   */
  public void getEigenvectorW(int i1, int i2, int i3, float[] w) {
    float[] wi = _uss.getPoint(_iw[i3][i2][i1]);
    w[0] = wi[0];
    w[1] = wi[1];
    w[2] = wi[2];
  }

  /**
   * Gets the eigenvector w for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @return array {w1,w2,w3} of eigenvector components.
   */
  public float[] getEigenvectorW(int i1, int i2, int i3) {
    float[] w = new float[3];
    getEigenvectorW(i1,i2,i3,w);
    return w;
  }

  /**
   * Sets tensor elements {s11,s12,s13,s22,s23,s33} for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param s array of tensor elements {s11,s12,s13,s22,s23,s33}.
   */
  public void setTensor(int i1, int i2, int i3, float[] s) {
    float[][] a = new float[][]{
      {s[0],s[1],s[2]},
      {s[1],s[3],s[4]},
      {s[2],s[4],s[5]}
    };
    float[][] uvw = new float[3][3];
    float[] e = new float[3];
    Eigen.solveSymmetric33(a,uvw,e);
    setEigenvalues(i1,i2,i3,e);
    setEigenvectorU(i1,i2,i3,uvw[0]);
    setEigenvectorW(i1,i2,i3,uvw[2]);
  }

  /**
   * Sets the eigenvalues for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param e array of eigenvalues {eu,ev,ew}.
   */
  public void setEigenvalues(int i1, int i2, int i3, float[] e) {
    Check.argument(e[0]>=e[1],"e[0]>=e[1]");
    Check.argument(e[1]>=e[2],"e[1]>=e[2]");
    Check.argument(e[2]>=0.0f,"e[2]>=0.0");
    float eu = e[0];
    float ev = e[1];
    float ew = e[2];
    float es = (eu>0.0f)?ES_SET/eu:0.0f;
    _eu[i3][i2][i1] = eu;
    _bv[i3][i2][i1] = (byte)(ev*es+0.5f);
    _bw[i3][i2][i1] = (byte)(ew*es+0.5f);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float ES_SET = 255.0f;
  private static final float ES_GET = 1.0f/ES_SET;
  private static UnitSphereSampling _uss = new UnitSphereSampling(16);

  private int _n1,_n2,_n3;
  private float[][][] _eu;
  private byte[][][] _bv;
  private byte[][][] _bw;
  private short[][][] _iu;
  private short[][][] _iw;

  private void setEigenvectorU(int i1, int i2, int i3, float[] u) {
    _iu[i3][i2][i1] = (short)_uss.getIndex(u);
  }

  private void setEigenvectorW(int i1, int i2, int i3, float[] w) {
    _iw[i3][i2][i1] = (short)_uss.getIndex(w);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    //int n1 = 3, n2 = 4, n3 = 5;
    int n1 = 1, n2 = 1, n3 = 1;
    StructureTensors3 st = new StructureTensors3(n1,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] s = makeRandomTensor();
          float[] t = new float[6];
          st.setTensor(i1,i2,i3,s);
          st.getTensor(i1,i2,i3,t);
          checkTensors(s,t);
        }
      }
    }
  }

  private static void checkVectors(float[] u, float[] v) {
    float uv = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
    float uu = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
    float vv = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    double angle = Math.toDegrees(Math.acos(uv));
    System.out.println("angle="+angle+" uu="+uu+" vv="+vv);
  }

  private static void checkTensors(float[] s, float[] t) {
    edu.mines.jtk.util.Array.dump(s);
    edu.mines.jtk.util.Array.dump(t);
    float e = 0.0f, d = 0.0f;
    e += Math.abs(s[0]-t[0]); d += Math.abs(s[0]);
    e += Math.abs(s[1]-t[1]); d += Math.abs(s[1]);
    e += Math.abs(s[2]-t[2]); d += Math.abs(s[2]);
    e += Math.abs(s[3]-t[3]); d += Math.abs(s[3]);
    e += Math.abs(s[4]-t[4]); d += Math.abs(s[4]);
    e += Math.abs(s[5]-t[5]); d += Math.abs(s[5]);
    System.out.println("e/d="+(e/d));
  }

  // A random positive-semidefinite tensor.
  private static float[] makeRandomTensor() {
    java.util.Random r = new java.util.Random();
    float a11 = r.nextFloat(), a12 = r.nextFloat(), a13 = r.nextFloat();
    float a21 = r.nextFloat(), a22 = r.nextFloat(), a23 = r.nextFloat();
    float a31 = r.nextFloat(), a32 = r.nextFloat(), a33 = r.nextFloat();
    float b11 = a11, b12 = a21, b13 = a31;
    float b21 = a12, b22 = a22, b23 = a32;
    float b31 = a13, b32 = a23, b33 = a33;
    float c11 = b11*a11+b12*a21+b13*a31;
    float c12 = b11*a12+b12*a22+b13*a32;
    float c13 = b11*a13+b12*a23+b13*a33;
    float c22 = b21*a12+b22*a22+b23*a32;
    float c23 = b21*a13+b22*a23+b23*a33;
    float c33 = b31*a13+b32*a23+b33*a33;
    return new float[]{c11,c12,c13,c22,c23,c33};
  }
}
