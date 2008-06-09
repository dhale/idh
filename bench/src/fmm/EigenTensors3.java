/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import ldf.UnitSphereSampling;

/**
 * An array of eigen-decompositions of tensors for 3D image processing. 
 * Each tensor is a symmetric positive-semidefinite 3-by-3 matrix:
 * <pre><code>
 *     |a11 a12 a13|
 * A = |a12 a22 a23|
 *     |a13 a23 a33|
 * </code></pre>
 * Such tensors facilitate anisotropic image processing.
 * <p>
 * The eigen-decomposition of the matrix A is
 * <pre><code>
 * A = au*u*u' + av*v*v' + aw*w*w'; au &gt;= av &gt;= aw &gt;= 0
 * </code></pre>
 * where u, v, and w are orthogonal unit eigenvectors of A. (The notation 
 * u' denotes the transpose of u.) The outer products of eigenvectors are
 * scaled by the corresponding eigenvalues au, av, and aw.
 * <p>
 * The ordering among eigenvalues is easily ensured by the equivalent 
 * representation
 * <pre><code>
 * A = (a1+a2+a3)*u*u' + (a2+a3)*v*v' + a3*w*w' 
 *     a1*u*u' + a2*(u*u'+v*v') + a3*(u*u'+v*v'+w*w')
 *   = a1*u*u' + a2*(I-w*w') + a3*I
 * </code></pre>
 * With this representation, the non-negative coefficients (a1,a2,a3)
 * are unordered, and the redundancy of the eigenvector v is apparent.
 * <p>
 * An intuitive interpretation of the coefficients (a1,a2,a3) is that
 * a1 corresponds to 1D scaling in the direction of u, a2 corresponds 
 * to 2D scaling in a plane orthogonal to w, and a3 corresponds to 3D 
 * isotropic scaling.
 * <p>
 * The coefficients (a1,a2,a3) are local; they may vary from sample to
 * sample. To provide additional control for 1D, 2D and 3D scaling, the 
 * coefficients (a1,a2,a3) may be multiplied by non-negative constant 
 * factors (s1,s2,s3), so that
 * <pre><code>
 * A = s1*a1*u*u' + s2*a2*(I-w*w') + s3*a3*I
 * </code></pre>
 * Only the 1st and 2nd components of the eigenvectors u and w are stored. 
 * Except for a sign, the 3rd components may be computed from the 1st and 
 * 2nd. Because the tensors are independent of the choice of sign, the 
 * eigenvectors u and v are stored with an implied non-negative 3rd 
 * component.
 * <p>
 * Storage may be further reduced by compression, whereby coefficients
 * and vectors are quantized. Quantization errors for coefficients
 * (a1,a2,a3) are less than 0.001*(a1+a2+a3). Quantization errors for 
 * eigenvectors are less than one degree of arc on the unit sphere.
 * Memory required to store each tensor is 12 bytes if compressed, and
 * 28 bytes if not compressed.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.07
 */
public class EigenTensors3 {

  /**
   * Constructs tensors for specified array dimensions. All coefficients 
   * (a1,a2,a3) and eigenvectors u and w are not set and are initially zero.
   * @param n1 number of tensors in 1st dimension.
   * @param n2 number of tensors in 2nd dimension.
   * @param n3 number of tensors in 3rd dimension.
   * @param s1 global scale factor for 1D coefficients.
   * @param s2 global scale factor for 2D coefficients.
   * @param s3 global scale factor for 3D coefficients.
   * @param compressed true, for compressed tensors; false, otherwise.
   */
  public EigenTensors3(
    int n1, int n2, int n3,
    double s1, double s2, double s3,
    boolean compressed)
  {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _s1 = (float)s1;
    _s2 = (float)s2;
    _s3 = (float)s3;
    _compressed = compressed;
    if (compressed) {
      _b1 = new short[n3][n2][n1];
      _b2 = new short[n3][n2][n1];
      _iu = new short[n3][n2][n1];
      _iw = new short[n3][n2][n1];
    } else {
      _a1 = new float[n3][n2][n1];
      _a2 = new float[n3][n2][n1];
      _u1 = new float[n3][n2][n1];
      _u2 = new float[n3][n2][n1];
      _w1 = new float[n3][n2][n1];
      _w2 = new float[n3][n2][n1];
    }
    _as = new float[n3][n2][n1];
  }

  /**
   * Constructs tensors for specified array dimensions and coefficients.
   * The 3rd components of eigenvectors u and v are computed from the 1st 
   * and 2nd components and are assumed to be non-negative.
   * @param u1 array of 1st components of u.
   * @param u2 array of 2nd components of u.
   * @param w1 array of 1st components of w.
   * @param w2 array of 2nd components of w.
   * @param a1 array of 1D coefficients.
   * @param a2 array of 2D coefficients.
   * @param a3 array of 3D coefficients.
   * @param s1 global scale factor for 1D coefficients.
   * @param s2 global scale factor for 2D coefficients.
   * @param s3 global scale factor for 3D coefficients.
   * @param compressed true, for compressed tensors; false, otherwise.
   */
  public EigenTensors3(
    float[][][] u1, float[][][] u2,
    float[][][] w1, float[][][] w2,
    float[][][] a1, float[][][] a2, float[][][] a3,
    double s1, double s2, double s3,
    boolean compressed)
  {
    this(a1[0][0].length,a1[0].length,a1.length,s1,s2,s3,compressed);
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          float a1i = a1[i3][i2][i1];
          float a2i = a2[i3][i2][i1];
          float a3i = a3[i3][i2][i1];
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = c3(u1i,u2i);
          float w1i = w1[i3][i2][i1];
          float w2i = w2[i3][i2][i1];
          float w3i = c3(w1i,w2i);
          setCoefficients(i1,i2,i3,a1i,a2i,a3i);
          setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
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
   * Gets tensor elements {a11,a12,a13,a22,a23,a33} for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param a array {a11,a12,a13,a22,a23,a33} of tensor elements.
   */
  public void getTensor(int i1, int i2, int i3, float[] a) {
    float asum = _as[i3][i2][i1];
    float a1,a2,a3,u1,u2,u3,w1,w2,w3;
    if (_compressed) {
      float ascale = asum*AS_GET;
      a1 = ascale*_b1[i3][i2][i1];
      a2 = ascale*_b2[i3][i2][i1];
      float[] u = _uss.getPoint(_iu[i3][i2][i1]);
      u1 = u[0]; u2 = u[1]; u3 = u[2];
      float[] w = _uss.getPoint(_iw[i3][i2][i1]);
      w1 = w[0]; w2 = w[1]; w3 = w[2];
    } else {
      a1 = _a1[i3][i2][i1];
      a2 = _a2[i3][i2][i1];
      u1 = _u1[i3][i2][i1];
      u2 = _u2[i3][i2][i1];
      u3 = c3(u1,u2);
      w1 = _w1[i3][i2][i1];
      w2 = _w2[i3][i2][i1];
      w3 = c3(w1,w2);
    }
    a3 = asum-a1-a2;
    a1 *= _s1;
    a2 *= _s2;
    a3 *= _s3;
    a[0] = a1*u1*u1+a2*(1.0f-w1*w1)+a3; // a11
    a[1] = a1*u1*u2+a2*(    -w1*w2)   ; // a12
    a[2] = a1*u1*u3+a2*(    -w1*w3)   ; // a13
    a[3] = a1*u2*u2+a2*(1.0f-w2*w2)+a3; // a22
    a[4] = a1*u2*u3+a2*(    -w2*w3)   ; // a23
    a[5] = a1*u3*u3+a2*(1.0f-w3*w3)+a3; // a33
  }


  /**
   * Gets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param a array {a1,a2,a3} of coefficients.
   */
  public void getCoefficients(int i1, int i2, int i3, float[] a) {
    float asum = _as[i3][i2][i1];
    float a1,a2;
    if (_compressed) {
      float ascale = asum*AS_GET;
      a1 = ascale*_b1[i3][i2][i1];
      a2 = ascale*_b2[i3][i2][i1];
    } else {
      a1 = _a1[i3][i2][i1];
      a2 = _a2[i3][i2][i1];
    }
    a[0] = a1; 
    a[1] = a2; 
    a[2] = asum-a1-a2;
  }

  /**
   * Gets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @return array {a1,a2,a3} of coefficients.
   */
  public float[] getCoefficients(int i1, int i2, int i3) {
    float[] a = new float[3];
    getCoefficients(i1,i2,i3,a);
    return a;
  }

  /**
   * Gets the eigenvector u for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param u array {u1,u2,u3} of eigenvector components.
   */
  public void getEigenvectorU(int i1, int i2, int i3, float[] u) {
    if (_compressed) {
      float[] ui = _uss.getPoint(_iu[i3][i2][i1]);
      u[0] = ui[0];
      u[1] = ui[1];
      u[2] = ui[2];
    } else {
      u[0] = _u1[i3][i2][i1];
      u[1] = _u2[i3][i2][i1];
      u[2] = c3(u[0],u[1]);
    }
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
    if (_compressed) {
      float[] wi = _uss.getPoint(_iw[i3][i2][i1]);
      w[0] = wi[0];
      w[1] = wi[1];
      w[2] = wi[2];
    } else {
      w[0] = _w1[i3][i2][i1];
      w[1] = _w2[i3][i2][i1];
      w[2] = c3(w[0],w[1]);
    }
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
   * Sets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param a1 1D coefficient.
   * @param a2 2D coefficient.
   * @param a3 3D coefficient.
   */
  public void setCoefficients(
    int i1, int i2, int i3, float a1, float a2, float a3)
  {
    float asum = a1+a2+a3;
    if (_compressed) {
      float ascale = (asum>0.0f)?AS_SET/asum:0.0f;
      _b1[i3][i2][i1] = (short)(a1*AS_SET+0.5f);
      _b2[i3][i2][i1] = (short)(a2*AS_SET+0.5f);
    } else {
      _a1[i3][i2][i1] = a1;
      _a2[i3][i2][i1] = a2;
    }
    _as[i3][i2][i1] = asum;
  }

  /**
   * Sets coefficients for the tensor with specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param a array {a1,a2,a3} of coefficients.
   */
  public void setCoefficients(int i1, int i2, int i3, float[] a) {
    setCoefficients(i1,i2,i3,a[0],a[1],a[2]);
  }

  /**
   * Sets the eigenvector u for the tensor with specified indices.
   * The specified vector is assumed to have length one. If the 3rd 
   * component is negative, this method stores the negative of the 
   * specified vector, so that the 3rd component is positive.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param u1 1st component of u.
   * @param u2 2nd component of u.
   * @param u3 3nd component of u.
   */
  public void setEigenvectorU(
    int i1, int i2, int i3, float u1, float u2, float u3)
  {
    if (u3<0.0f) {
      u1 = -u1;
      u2 = -u2;
      u3 = -u3;
    }
    if (_compressed) {
      _iu[i3][i2][i1] = (short)_uss.getIndex(u1,u2,u3);
    } else {
      _u1[i3][i2][i1] = u1;
      _u2[i3][i2][i1] = u2;
    }
  }

  /**
   * Sets the eigenvector u for the tensor with specified indices.
   * The specified vector is assumed to have length one. If the 3rd 
   * component is negative, this method stores the negative of the 
   * specified vector, so that the 3rd component is positive.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param array {u1,u2,u3} of eigenvector components.
   */
  public void setEigenvectorU(int i1, int i2, int i3, float[] u) {
    setEigenvectorU(i1,i2,i3,u[0],u[1],u[2]);
  }

  /**
   * Sets the eigenvector w for the tensor with specified indices.
   * The specified vector is assumed to have length one. If the 3rd 
   * component is negative, this method stores the negative of the 
   * specified vector, so that the 3rd component is positive.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param w1 1st component of w.
   * @param w2 2nd component of w.
   * @param w3 3nd component of w.
   */
  public void setEigenvectorW(
    int i1, int i2, int i3, float w1, float w2, float w3)
  {
    if (w3<0.0f) {
      w1 = -w1;
      w2 = -w2;
      w3 = -w3;
    }
    if (_compressed) {
      _iw[i3][i2][i1] = (short)_uss.getIndex(w1,w2,w3);
    } else {
      _w1[i3][i2][i1] = w1;
      _w2[i3][i2][i1] = w2;
    }
  }

  /**
   * Sets the eigenvector w for the tensor with specified indices.
   * The specified vector is assumed to have length one. If the 3rd 
   * component is negative, this method stores the negative of the 
   * specified vector, so that the 3rd component is positive.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param array {w1,w2,w3} of eigenvector components.
   */
  public void setEigenvectorW(int i1, int i2, int i3, float[] w) {
    setEigenvectorW(i1,i2,i3,w[0],w[1],w[2]);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float AS_SET = (float)Short.MAX_VALUE;
  private static final float AS_GET = 1.0f/AS_SET;
  private static UnitSphereSampling _uss = new UnitSphereSampling(16);

  private boolean _compressed; // true if tensors compressed
  private int _n1,_n2,_n3; // array dimensions
  private float _s1,_s2,_s3; // multipliers for (a1,a2,a3)
  private float[][][] _as; // sum a1+a2+a3
  private short[][][] _b1; // a1 compressed
  private short[][][] _b2; // a2 compressed
  private short[][][] _iu; // (u1,u2,u3) compressed
  private short[][][] _iw; // (w1,w2,w3) compressed
  private float[][][] _a1; // a1 not compressed
  private float[][][] _a2; // a2 not compressed
  private float[][][] _u1; // u1 not compressed
  private float[][][] _u2; // u2 not compressed
  private float[][][] _w1; // w1 not compressed
  private float[][][] _w2; // w2 not compressed

  private static float c3(float c1, float c2) {
    float c3s = 1.0f-c1*c1-c2*c2;
    return (c3s>0.0f)?(float)Math.sqrt(c3s):0.0f;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    testRandom();
  }

  private static void testRandom() {
    testRandom(false,0.1,1.0e-6);
    testRandom(true,1.0,1.0e-4);
  }

  private static void testRandom(
    boolean compressed, double errorAngle, double errorCoeff) 
  {
    int n1 = 3, n2 = 4, n3 = 5;
    double s1 = 1.1, s2 = 1.2, s3 = 1.3;
    EigenTensors3 dt = new EigenTensors3(n1,n2,n3,s1,s2,s3,false);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] a = makeRandomCoefficients();
          float[] u = makeRandomVector();
          float[] w = makeOrthogonalVector(u);
          dt.setCoefficients(i1,i2,i3,a);
          dt.setEigenvectorU(i1,i2,i3,u);
          dt.setEigenvectorW(i1,i2,i3,w);
          float[] c;
          c = dt.getEigenvectorU(i1,i2,i3); checkVectors(u,c,errorAngle);
          c = dt.getEigenvectorW(i1,i2,i3); checkVectors(w,c,errorAngle);
          c = dt.getCoefficients(i1,i2,i3); checkCoefficients(c,a,errorCoeff);
        }
      }
    }
  }

  private static void checkCoefficients(float[] c, float[] a, double error) {
    float e1 = Math.abs(c[0]-a[0]);
    float e2 = Math.abs(c[1]-a[1]);
    float e3 = Math.abs(c[2]-a[2]);
    //System.out.println("e1="+e1+" e2="+e2+" e3="+e3);
    assert e1<error:"error in a1 less than "+error;
    assert e2<error:"error in a2 less than "+error;
    assert e3<error:"error in a3 less than "+error;
  }

  private static void checkVectors(float[] u, float[] v, double error) {
    float uv = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
    float uu = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
    float vv = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    double ca = Math.min(uv,1.0f);
    double angle = Math.toDegrees(Math.acos(ca));
    //System.out.println("angle="+angle+" uv="+uv+" uu="+uu+" vv="+vv);
    assert angle<error:"angle between u and v less than "+error+" degrees";
  }

  private static void checkTensors(float[] s, float[] t, double error) {
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

  private static float[] normalize(float[] a) {
    float s = 1.0f/(a[0]+a[1]+a[2]);
    return new float[]{a[0]*s,a[1]*s,a[2]*s};
  }

  private static java.util.Random r = new java.util.Random();

  // Random coefficients, not normalized.
  private static float[] makeRandomCoefficients() {
    float d1 = r.nextFloat();
    float d2 = r.nextFloat();
    float d3 = r.nextFloat();
    //float ds = 1.0f/(d1+d2+d3);
    //return new float[]{d1*ds,d2*ds,d3*ds};
    return new float[]{d1,d2,d3};
  }

  // Random unit vector.
  private static float[] makeRandomVector() {
    float a = r.nextFloat()-0.5f;
    float b = r.nextFloat()-0.5f;
    float c = r.nextFloat()-0.5f;
    if (c<0.0f) {
      a = -a;
      b = -b;
      c = -c;
    }
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
    if (c<0.0f) {
      a = -a;
      b = -b;
      c = -c;
    }
    s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }
}
