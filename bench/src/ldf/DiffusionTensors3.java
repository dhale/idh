/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

/**
 * A 3-D array of diffusion tensors. Each diffusion tensor is a symmetric 
 * positive-semidefinite 3-by-3 matrix:
 * <pre><code>
 *     |d11 d12 d13|
 * D = |d12 d22 d23|
 *     |d13 d23 d33|
 * </code></pre>
 * The six independent elements in this matrix may be chosen to cause
 * anisotropic diffusion in three dimensions. This anisotropy is more 
 * easily described by an eigen-decomposition of the diffusion tensors
 * <pre><code>
 * D = eu*u*u' + ev*v*v' + ew*w*w'
 * </code></pre>
 * where u, v, and w are orthogonal unit eigenvectors of D. (The notation 
 * u' denotes the transpose of u.) The outer products of eigenvectors are
 * scaled by the corresponding eigenvalues eu, ev, and ew. Here, these
 * eigenvalues are normalized such that 0 &lt;= eu &lt;= ev &lt;=ew = 1.
 * <p>
 * These relationships among eigenvalues are more easily enforced by an 
 * alternative representation:
 * <pre><code>
 * D =  d3*u*u' + (d2+d3)*v*v' + (d1+d2+d3)*w*w'
 *   =  d3*u*u' + (d2+d3)*v*v' + w*w'
 *   =  d1*w*w' + d2*(I-u*u') + d3*I
 * </code></pre>
 * where d1 controls the amount of 1-D diffusion along lines parallel
 * to the eigenvector w, d2 controls the amount of 2-D diffusion within
 * planes perpendicular to the eigenvector u, and d3 controls the amount
 * of 3-D isotropic diffusion. (The symbol I denotes a 3-by-3 identity
 * matrix.) All three diffusion coefficients d1, d2, and d3 are in the 
 * range [0,1], and are normalized so that their sum d1+d2+d3 = 1.
 * <p>
 * Normalization implies that d1, d2, and d3 are not independent. While
 * they may be specified for each sample, only two of them must be stored; 
 * the third coefficient can be readily computed from the other two.
 * <p>
 * Normalization reduces storage requirements but implies a constant
 * amount of diffusion (d1+d2+d3 = 1) along the eigenvector w. To 
 * facilitate more general diffusions, global scale factors s1, s2, 
 * and s3 may be may be specified to obtain diffusion tensors
 * <pre><code>
 * D =  s1*s1*d1*w*w' + s2*s2*d2*(I-u*u') + s3*s3*(1-d1-d2)*I
 * </code></pre>
 * Theese scale factors are squared so that they more closely correspond
 * to diffusion distances.
 * <p>
 * Computation of eigenvalues and eigenvectors is straightforward
 * but can be costly if performed repeatedly as needed. And caching
 * of eigenvalues and eigenvectors computed once can require more
 * memory than the six 3-D arrays required to store tensor elements.
 * <p> 
 * To balance computational efficiency with memory requirements, the
 * eigen-decomposition of each diffusion tensor is computed once and 
 * stored approximately in only 6 bytes per tensor. Coefficients d2 
 * and d3 are approximated with less than 0.5% error, and angles 
 * between approximate and exact eigenvectors u and w are less than 
 * one degree. Eigenvectors v are computed by cross products v = w x u, 
 * and coefficients d1 = 1-d2-d3.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.02.18
 */
public class DiffusionTensors3 {

  /**
   * Constructs tensors for specified array dimensions.
   * @param n1 number of tensors in 1st dimension.
   * @param n2 number of tensors in 2nd dimension.
   * @param n3 number of tensors in 3rd dimension.
   * @param s1 scale factor for 1-D diffusion.
   * @param s2 scale factor for 2-D diffusion.
   * @param s3 scale factor for 3-D diffusion.
   */
  public DiffusionTensors3(
    int n1, int n2, int n3,
    double s1, double s2, double s3) 
  {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _ss1 = (float)(s1*s1);
    _ss2 = (float)(s2*s2);
    _ss3 = (float)(s3*s3);
    _b2 = new byte[n3][n2][n1];
    _b3 = new byte[n3][n2][n1];
    _iu = new short[n3][n2][n1];
    _iw = new short[n3][n2][n1];
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
   * Gets tensor elements {d11,d12,d13,d22,d23,d33} for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param d array {d11,d12,d13,d22,d23,d33} of tensor elements.
   */
  public void getTensor(int i1, int i2, int i3, float[] d) {
    float b2i = _b2[i3][i2][i1];
    float b3i = _b3[i3][i2][i1];
    if (b2i<0.0) b2i += 256.0f;
    if (b3i<0.0) b3i += 256.0f;
    float d2i = b2i*DS_GET;
    float d3i = b3i*DS_GET;
    float d1i = 1.0f-d2i-d3i;
    float[] u = _uss.getPoint(_iu[i3][i2][i1]);
    float[] w = _uss.getPoint(_iw[i3][i2][i1]);
    float u1i = u[0];
    float u2i = u[1];
    float u3i = u[2];
    float w1i = w[0];
    float w2i = w[1];
    float w3i = w[2];
    float e1i = d1i*_ss1;
    float e2i = d2i*_ss2;
    float e3i = d3i*_ss3;
    d[0] = e1i*w1i*w1i+e2i*(1.0f-u1i*u1i)+e3i; // d11
    d[1] = e1i*w1i*w2i+e2i*(    -u1i*u2i)    ; // d12
    d[2] = e1i*w1i*w3i+e2i*(    -u1i*u3i)    ; // d13
    d[3] = e1i*w2i*w2i+e2i*(1.0f-u2i*u2i)+e3i; // d22
    d[4] = e1i*w2i*w3i+e2i*(    -u2i*u3i)    ; // d23
    d[5] = e1i*w3i*w3i+e2i*(1.0f-u3i*u3i)+e3i; // d33
  }

  /**
   * Gets diffusion coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param d array {d1,d2,d3} of diffusion coefficients.
   */
  public void getCoefficients(int i1, int i2, int i3, float[] d) {
    float b2i = _b2[i3][i2][i1];
    float b3i = _b3[i3][i2][i1];
    if (b2i<0.0) b2i += 256.0f;
    if (b3i<0.0) b3i += 256.0f;
    float d2 = b2i*DS_GET;
    float d3 = b3i*DS_GET;
    float d1 = 1.0f-d2-d3;
    d[0] = d1;
    d[1] = d2;
    d[2] = d3;
  }

  /**
   * Gets diffusion coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @return array {d1,d2,d3} of diffusion coefficients.
   */
  public float[] getCoefficients(int i1, int i2, int i3) {
    float[] d = new float[3];
    getCoefficients(i1,i2,i3,d);
    return d;
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
   * Sets diffusion coefficients for the tensor with specified indices.
   * If necessary, the specified coefficients are scaled so that their 
   * sum d1+d2+d3 = 1.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param d array {d1,d2,d3} of diffusion coefficients.
   */
  public void setCoefficients(int i1, int i2, int i3, float[] d) {
    float d1 = d[0];
    float d2 = d[1];
    float d3 = d[2];
    float ds = d1+d2+d3;
    ds = (ds>0.0f)?DS_SET/ds:0.0f;
    _b2[i3][i2][i1] = (byte)(d2*ds+0.5f);
    _b3[i3][i2][i1] = (byte)(d3*ds+0.5f);
  }

  /**
   * Sets the eigenvector u for the tensor with specified indices.
   * The unit vector u is orthogonal to the plane of 2-D diffusion.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param array {u1,u2,u3} of eigenvector components.
   */
  public void setEigenvectorU(int i1, int i2, int i3, float[] u) {
    _iu[i3][i2][i1] = (short)_uss.getIndex(u);
  }

  /**
   * Sets the eigenvector w for the tensor with specified indices.
   * The unit vector w is parallel to the line of 1-D diffusion.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param array {w1,w2,w3} of eigenvector components.
   */
  public void setEigenvectorW(int i1, int i2, int i3, float[] w) {
    _iw[i3][i2][i1] = (short)_uss.getIndex(w);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float DS_SET = 255.0f;
  private static final float DS_GET = 1.0f/DS_SET;
  private static UnitSphereSampling _uss = new UnitSphereSampling(16);

  private int _n1,_n2,_n3;
  private float _ss1,_ss2,_ss3;
  private byte[][][] _b2;
  private byte[][][] _b3;
  private short[][][] _iu;
  private short[][][] _iw;

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    //int n1 = 3, n2 = 4, n3 = 5;
    int n1 = 1, n2 = 1, n3 = 1;
    double s1 = 1.0, s2 = 1.0, s3 = 1.0;
    DiffusionTensors3 dt = new DiffusionTensors3(n1,n2,n3,s1,s2,s3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] d = makeRandomCoefficients();
          float[] u = makeRandomVector();
          float[] w = makeOrthogonalVector(u);
          dt.setCoefficients(i1,i2,i3,d);
          dt.setEigenvectorU(i1,i2,i3,u);
          dt.setEigenvectorW(i1,i2,i3,w);
          float[] c;
          c = dt.getEigenvectorU(i1,i2,i3); checkVectors(u,c);
          c = dt.getEigenvectorW(i1,i2,i3); checkVectors(w,c);
          c = dt.getCoefficients(i1,i2,i3); checkCoefficients(c,d);
        }
      }
    }
  }

  private static void checkCoefficients(float[] c, float[] d) {
    float e1 = c[0]-d[0];
    float e2 = c[1]-d[1];
    float e3 = c[2]-d[2];
    System.out.println("e1="+e1+" e2="+e2+" e3="+e3);
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

  private static java.util.Random r = new java.util.Random();

  // Random diffusion coefficients.
  private static float[] makeRandomCoefficients() {
    float d1 = r.nextFloat();
    float d2 = r.nextFloat();
    float d3 = r.nextFloat();
    float ds = 1.0f/(d1+d2+d3);
    return new float[]{d1*ds,d2*ds,d3*ds};
  }

  // Random unit vector.
  private static float[] makeRandomVector() {
    float a = r.nextFloat()-0.5f;
    float b = r.nextFloat()-0.5f;
    float c = r.nextFloat()-0.5f;
    float s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }

  // Random unit vector orthogonal to specified vector.
  private static float[] makeOrthogonalVector(float[] v1) {
    float a1 = v1[0];
    float b1 = v1[1];
    float c1 = v1[2];
    float a2 = r.nextFloat()-0.5f;
    float b2 = r.nextFloat()-0.5f;
    float c2 = r.nextFloat()-0.5f;
    float d11 = a1*a1+b1*b1+c1*c1;
    float d12 = a1*a2+b1*b2+c1*c2;
    float s = d12/d11;
    float a = a2-s*a1;
    float b = b2-s*b1;
    float c = c2-s*c1;
    s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }
}
