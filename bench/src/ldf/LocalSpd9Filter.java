/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Local symmetric positive-definite (SPD) filter with a 2-D 9-point stencil.
 * This filter is local in the sense that filter coefficients may differ for 
 * each sample. These coefficients form a 9-point stencil:
 * <pre><code>
 * y[i2][i1] = smm*x[i2-1][i1-1] + s0m*x[i2  ][i1-1] + spm*x[i2+1][i1-1] +
 *             sm0*x[i2-1][i1  ] + s00*x[i2  ][i1  ] + sp0*x[i2+1][i1  ] +
 *             smp*x[i2-1][i1+1] + s0p*x[i2  ][i1+1] + spp*x[i2+1][i1+1]
 * </code></pre>
 * The suffixes m, 0, and p denote minus, zero, and plus, respectively.
 * For example, smp is the coefficient of x[i2-1][i1+1].
 * <p>
 * For symmetric filters with constant coefficients, this stencil is 
 * symmetric about the central coefficient s00. For example, smm = spp.
 * Therefore, only five of the nine coefficients need be specified.
 * <p>
 * For symmetric filters with variable coefficients, this stencil is 
 * <em>not</em> symmetric. That is, smm[i2][i1] does not equal spp[i2][i1];
 * rather, smm[i2][i1] = spp[i2-1][i1-1]. Still, only five filter 
 * coefficients need be specified for each sample. If we choose those five
 * coefficients to be s00, s0p, spm, sp0, and spp, then
 * <pre><code>
 * y[i2][i1] = spp[i2-1][i1-1]*x[i2-1][i1-1] +
 *             sp0[i2-1][i1  ]*x[i2-1][i1  ] +
 *             spm[i2-1][i1+1]*x[i2-1][i1+1] +
 *             s0p[i2  ][i1-1]*x[i2  ][i1-1] +
 *             s00[i2  ][i1  ]*x[i2  ][i1  ] +
 *             s0p[i2  ][i1  ]*x[i2  ][i1+1] +
 *             spm[i2  ][i1  ]*x[i2+1][i1-1] +
 *             sp0[i2  ][i1  ]*x[i2+1][i1  ] +
 *             spp[i2  ][i1  ]*x[i2+1][i1+1]
 * </code></pre>
 * <p>
 * Becouse this filter is SPD, it may in theory be factored exactly with 
 * Cholesky decomposition. However, the factors seldom fit in a 9-point
 * stencil. Therefore, only approximate factors are typically computed 
 * using incomplete Cholesky (IC) decomposition. The factors may then be 
 * used to apply an approximate inverse of this filter. This approximate 
 * inverse is especially useful as a pre-conditioner in the method of 
 * conjugate gradients.
 * <p>
 * Unfortunately, IC decomposition may fail with non-positive pivots for 
 * filters that are not diagonally-dominant. To enable IC decomposition
 * to succeed, filter coefficients s00*(1+bias) may be used instead of s00. 
 * (Any bias is used only during IC decomposition; the specified s00 are 
 * not changed.) For filters known to be diagonally dominant, zero bias 
 * should be specified.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.11
 */
public class LocalSpd9Filter {

  /**
   * Constructs a filter with specified coefficients.
   * Any approximate inverse filter (when required) will be computed with an 
   * initial bias of zero.
   * @param s arrays[5][n2][n1] of coefficients; by reference, not by copy.
   *  The elements of the array s are {s00,s0p,spm,sp0,spp}.
   */
  public LocalSpd9Filter(float[][][] s) {
    this(s,0.0);
  }

  /**
   * Constructs a filter with specified coefficients.
   * @param s arrays[5][n2][n1] of coefficients; by reference, not by copy.
   *  The elements of the array s are {s00,s0p,spm,sp0,spp}.
   * @param bias the initial non-negative amount by which to perturb the 
   *  coefficients s00 during computation of an approximate inverse filter.
   */
  public LocalSpd9Filter(float[][][] s, double bias) {
    _s = s;
    _b = (float)bias;
  }

  /**
   * Applies this filter by computing y = A*x.
   * @param x input array. Must be distinct from y.
   * @param y output array. Must be distinct from x.
   */
  public void apply(float[][] x, float[][] y) {
    applyFilter(x,y);
  }

  /**
   * Applies this filter by computing y = A*x using an incomplete Cholesky 
   * decomposition of A. In effect, this method applies an approximation of 
   * this filter.
   * @param x input array. Must be distinct from y.
   * @param y output array. Must be distinct from x.
   */
  public void applyApproximate(float[][] x, float[][] y) {
    applyFactors(x,y);
  }

  /**
   * Solves A*x = y using an incomplete Cholesky decomposition of A.
   * In effect, this method applies an approximate inverse of this filter.
   * @param y the input right-hand side array.
   * @param x the output solution array.
   */
  public void applyApproximateInverse(float[][] y, float[][] x) {
    solveWithFactors(y,x);
  }

  /**
   * Gets the sparse matrix A equivalent to this filter.
   * Most elements in this matrix will be zero. For small numbers of
   * samples, it may be useful for visualization of matrix sparsity.
   * @return an array[n][n] representing A, where n = n1*n2.
   */
  public float[][] getMatrix() {
    float[][] s00 = _s[0];
    float[][] s0p = _s[1];
    float[][] spm = _s[2];
    float[][] sp0 = _s[3];
    float[][] spp = _s[4];
    int n1 = s00[0].length;
    int n2 = s00.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n = n1*n2;
    float[][] a = new float[n][n];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
                                     a[i     ][i] = s00[i2][i1];
        if (i1<n1m)   a[i][i+1   ] = a[i+1   ][i] = s0p[i2][i1];
        if (i2<n2m) {
          if (0<i1)   a[i][i+n1-1] = a[i+n1-1][i] = spm[i2][i1];
                      a[i][i+n1  ] = a[i+n1  ][i] = sp0[i2][i1];
          if (i1<n1m) a[i][i+n1+1] = a[i+n1+1][i] = spp[i2][i1];
        }
      }
    }
    return a;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _s; // specified SPD filter coefficients
  private float[][][] _l; // coefficients of IC(0) decomposition
  private float _b; // initial bias for IC(0) decomposition.

  /**
   * Makes the IC(0) factors, if not already made.
   */
  private void ensureFactors() {
    if (_l==null)
      _l = factorIC0(_s,_b);
    Check.state(_l!=null,"incomplete Cholesky decomposition successful");
  }

  /**
   * Computes y = A*x.
   */
  private void applyFilter(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] s00 = _s[0];
    float[][] s0p = _s[1];
    float[][] spm = _s[2];
    float[][] sp0 = _s[3];
    float[][] spp = _s[4];
    int i1,i2,i1m,i1p;
    i1 = n1m;
    i2 = n2m;
    y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
    for (i1=n1m-1; i1>=0; --i1) {
      y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += s0p[i2][i1]*x[i2  ][i1+1];
      y[i2  ][i1+1] += s0p[i2][i1]*x[i2  ][i1  ];
    }
    for (i2=n2m-1; i2>=0; --i2) {
      float[] s00i2 = s00[i2];
      float[] s0pi2 = s0p[i2];
      float[] spmi2 = spm[i2];
      float[] sp0i2 = sp0[i2];
      float[] sppi2 = spp[i2];
      float[] yi2p0 = y[i2  ];
      float[] yi2p1 = y[i2+1];
      float[] xi2p0 = x[i2  ];
      float[] xi2p1 = x[i2+1];
      i1 = n1m;
      yi2p0[i1  ]  = s00i2[i1]*xi2p0[i1  ];
      yi2p0[i1  ] += spmi2[i1]*xi2p1[i1-1];
      yi2p1[i1-1] += spmi2[i1]*xi2p0[i1  ];
      yi2p0[i1  ] += sp0i2[i1]*xi2p1[i1  ];
      yi2p1[i1  ] += sp0i2[i1]*xi2p0[i1  ];
      for (i1=n1m-1,i1m=i1-1,i1p=i1+1; i1>=1; --i1,--i1m,--i1p) {
        yi2p0[i1 ]  = s00i2[i1]*xi2p0[i1 ];
        yi2p0[i1 ] += s0pi2[i1]*xi2p0[i1p];
        yi2p0[i1p] += s0pi2[i1]*xi2p0[i1 ];
        yi2p0[i1 ] += spmi2[i1]*xi2p1[i1m];
        yi2p1[i1m] += spmi2[i1]*xi2p0[i1 ];
        yi2p0[i1 ] += sp0i2[i1]*xi2p1[i1 ];
        yi2p1[i1 ] += sp0i2[i1]*xi2p0[i1 ];
        yi2p0[i1 ] += sppi2[i1]*xi2p1[i1p];
        yi2p1[i1p] += sppi2[i1]*xi2p0[i1 ];
      }
      yi2p0[i1  ]  = s00i2[i1]*xi2p0[i1  ];
      yi2p0[i1  ] += s0pi2[i1]*xi2p0[i1+1];
      yi2p0[i1+1] += s0pi2[i1]*xi2p0[i1  ];
      yi2p0[i1  ] += sp0i2[i1]*xi2p1[i1  ];
      yi2p1[i1  ] += sp0i2[i1]*xi2p0[i1  ];
      yi2p0[i1  ] += sppi2[i1]*xi2p1[i1+1];
      yi2p1[i1+1] += sppi2[i1]*xi2p0[i1  ];
    }
  }

  /**
   * Solves L*D*L'*x = b.
   */
  private void solveWithFactors(float[][] b, float[][] x) {
    ensureFactors();
    int n1 = b[0].length;
    int n2 = b.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] d00 = _l[0];
    float[][] l0p = _l[1];
    float[][] lpm = _l[2];
    float[][] lp0 = _l[3];
    float[][] lpp = _l[4];
    int i1,i2;

    // Solve L*z = b.
    zero(x);
    for (i2=0; i2<n2m; ++i2) {
      float[] bi2p0 = b[i2  ];
      float[] xi2p0 = x[i2  ];
      float[] xi2p1 = x[i2+1];
      float[] l0pi2 = l0p[i2];
      float[] lpmi2 = lpm[i2];
      float[] lp0i2 = lp0[i2];
      float[] lppi2 = lpp[i2];
      i1 = 0;
      xi2p0[i1  ] += bi2p0[i1];
      xi2p0[i1+1] -= l0pi2[i1]*xi2p0[i1];
      xi2p1[i1  ] -= lp0i2[i1]*xi2p0[i1];
      xi2p1[i1+1] -= lppi2[i1]*xi2p0[i1];
      for (i1=1; i1<n1m; ++i1) {
        xi2p0[i1  ] += bi2p0[i1];
        xi2p0[i1+1] -= l0pi2[i1]*xi2p0[i1];
        xi2p1[i1-1] -= lpmi2[i1]*xi2p0[i1];
        xi2p1[i1  ] -= lp0i2[i1]*xi2p0[i1];
        xi2p1[i1+1] -= lppi2[i1]*xi2p0[i1];
      }
      xi2p0[i1  ] += bi2p0[i1];
      xi2p1[i1-1] -= lpmi2[i1]*xi2p0[i1];
      xi2p1[i1  ] -= lp0i2[i1]*xi2p0[i1];
    }
    i1 = 0;
    x[i2  ][i1  ] += b[i2][i1];
    x[i2  ][i1+1] -= l0p[i2][i1]*x[i2][i1];
    for (i1=1; i1<n1m; ++i1) {
      x[i2  ][i1  ] += b[i2][i1];
      x[i2  ][i1+1] -= l0p[i2][i1]*x[i2][i1];
    }
    x[i2  ][i1  ] += b[i2][i1];

    // Solve D*y = z and L'*x = y.
    i2 = n2m;
    i1 = n1m;
    x[i2][i1] = d00[i2][i1]*x[i2  ][i1  ];
    for (i1=n1m-1; i1>=0; --i1) {
      x[i2][i1] = d00[i2][i1]*x[i2  ][i1  ] -
                  l0p[i2][i1]*x[i2  ][i1+1];
    }
    for (i2=n2m-1; i2>=0; --i2) {
      float[] xi2p0 = x[i2  ];
      float[] xi2p1 = x[i2+1];
      float[] d00i2 = d00[i2];
      float[] l0pi2 = l0p[i2];
      float[] lpmi2 = lpm[i2];
      float[] lp0i2 = lp0[i2];
      float[] lppi2 = lpp[i2];
      i1 = n1m;
      xi2p0[i1] = d00[i2][i1]*xi2p0[i1  ] -
                  lp0[i2][i1]*xi2p1[i1  ] -
                  lpm[i2][i1]*xi2p1[i1-1];
      for (i1=n1m-1; i1>=1; --i1) {
        xi2p0[i1] = d00i2[i1]*xi2p0[i1  ] -
                    lppi2[i1]*xi2p1[i1+1] -
                    lp0i2[i1]*xi2p1[i1  ] -
                    lpmi2[i1]*xi2p1[i1-1] -
                    l0pi2[i1]*xi2p0[i1+1];
      }
      xi2p0[i1] = d00i2[i1]*xi2p0[i1  ] -
                  lppi2[i1]*xi2p1[i1+1] -
                  lp0i2[i1]*xi2p1[i1  ] -
                  l0pi2[i1]*xi2p0[i1+1];
    }
  }
  
  /**
   * Computes y = L*D*L'*x. For testing, only.
   */
  private void applyFactors(float[][] x, float[][] y) {
    ensureFactors();
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] d00 = _l[0];
    float[][] l0p = _l[1];
    float[][] lpm = _l[2];
    float[][] lp0 = _l[3];
    float[][] lpp = _l[4];
    int i1,i2;

    // y = L'*x
    for (i2=0; i2<n2m; ++i2) {
      i1 = 0;
      y[i2][i1] = x[i2  ][i1  ] +
                  x[i2  ][i1+1]*l0p[i2][i1] +
                  x[i2+1][i1  ]*lp0[i2][i1] +
                  x[i2+1][i1+1]*lpp[i2][i1];
      for (i1=1; i1<n1m; ++i1) {
        y[i2][i1] = x[i2  ][i1  ] +
                    x[i2  ][i1+1]*l0p[i2][i1] +
                    x[i2+1][i1-1]*lpm[i2][i1] +
                    x[i2+1][i1  ]*lp0[i2][i1] +
                    x[i2+1][i1+1]*lpp[i2][i1];
      }
      y[i2][i1] = x[i2  ][i1  ] +
                  x[i2+1][i1-1]*lpm[i2][i1] +
                  x[i2+1][i1  ]*lp0[i2][i1];
    }
    for (i1=0; i1<n1m; ++i1) {
      y[i2][i1] = x[i2  ][i1  ] +
                  x[i2  ][i1+1]*l0p[i2][i1];
    }
    y[i2][i1] = x[i2  ][i1  ];

    // y = L*D*y
    i2 = n2m;
    i1 = n1m;
    y[i2  ][i1  ] /= d00[i2][i1];
    for (i1=n1m-1; i1>=0; --i1) {
      y[i2  ][i1  ] /= d00[i2][i1];
      y[i2  ][i1+1] += l0p[i2][i1]*y[i2][i1];
    }
    for (i2=n2m-1; i2>=0; --i2) {
      i1 = n1m;
      y[i2  ][i1  ] /= d00[i2][i1];
      y[i2+1][i1-1] += lpm[i2][i1]*y[i2][i1];
      y[i2+1][i1  ] += lp0[i2][i1]*y[i2][i1];
      for (i1=n1m-1; i1>=1; --i1) {
        y[i2  ][i1  ] /= d00[i2][i1];
        y[i2  ][i1+1] += l0p[i2][i1]*y[i2][i1];
        y[i2+1][i1-1] += lpm[i2][i1]*y[i2][i1];
        y[i2+1][i1  ] += lp0[i2][i1]*y[i2][i1];
        y[i2+1][i1+1] += lpp[i2][i1]*y[i2][i1];
      }
      y[i2  ][i1  ] /= d00[i2][i1];
      y[i2  ][i1+1] += l0p[i2][i1]*y[i2][i1];
      y[i2+1][i1  ] += lp0[i2][i1]*y[i2][i1];
      y[i2+1][i1+1] += lpp[i2][i1]*y[i2][i1];
    }
  }

  /**
   * Factors this filter with incomplete Cholesky decomposition IC(0).
   * The approximate factorization is A ~ L*D*L', where L is a 
   * unit-lower-triangular matrix, and D is a diagonal matrix.
   * @param bias the amount by which to perturb the diagonal of A.
   * @return The coefficients {d00,l0p,lpm,lp0,lpp}. The elements in the
   *  array d00 are the inverse of the elements of the diagonal matrix D.
   */
  private static float[][][] factorIC0(float[][][] a, float bias) {
    float[][][] l = null;
    float bmin = (bias>0.0f)?bias:0.001f;
    for (float b=bias; l==null; b=max(bmin,2*b)) {
      l = attemptIC0(a,b);
      if (l==null)
        trace("factorIC0: failed for bias="+b);
      else
        trace("factorIC0: success for bias="+b);
    }
    return l;
  }
  private static float[][][] attemptIC0(float[][][] a, float bias) {
    float[][][] l = copy(a);
    if (bias>0.0f)
      mul(1.0f+bias,l[0],l[0]);
    float[][] l00 = l[0], l0p = l[1], lpm = l[2], lp0 = l[3], lpp = l[4];
    float[][] d00 = l00; // will contain inverse of diagonal matrix D

    // Incomplete Cholesky decomposition, in-place.
    int i1m,i1p,i2m;
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    int i1 = 0;
    int i2 = 0;
    d00[i2][i1] = 1.0f/l00[i2][i1];
    for (i1=1,i1m=i1-1,i1p=i1+1; i1<n1; ++i1,++i1m,++i1p) {
      l00[i2][i1] -= d00[i2 ][i1m]*l0p[i2 ][i1m]*l0p[i2 ][i1m];
      lpm[i2][i1] -= d00[i2 ][i1m]*lp0[i2 ][i1m]*l0p[i2 ][i1m];
      lp0[i2][i1] -= d00[i2 ][i1m]*lpp[i2 ][i1m]*l0p[i2 ][i1m];
      if (l00[i2][i1]<=0.0f) 
        return null;
      d00[i2][i1] = 1.0f/l00[i2][i1];
    }
    for (i2=1,i2m=i2-1; i2<n2; ++i2,++i2m) {
      i1 = 0; i1p = i1+1;
      l00[i2][i1] -= d00[i2m][i1p]*lpm[i2m][i1p]*lpm[i2m][i1p] +
                     d00[i2m][i1 ]*lp0[i2m][i1 ]*lp0[i2m][i1 ];
      l0p[i2][i1] -= d00[i2m][i1p]*lp0[i2m][i1p]*lpm[i2m][i1p] +
                     d00[i2m][i1 ]*lpp[i2m][i1 ]*lp0[i2m][i1 ];
      if (l00[i2][i1]<=0.0f) 
        return null;
      d00[i2][i1] = 1.0f/l00[i2][i1];
      for (i1=1,i1m=i1-1,i1p=i1+1; i1<n1-1; ++i1,++i1m,++i1p) {
        l00[i2][i1] -= d00[i2 ][i1m]*l0p[i2 ][i1m]*l0p[i2 ][i1m] +
                       d00[i2m][i1p]*lpm[i2m][i1p]*lpm[i2m][i1p] +
                       d00[i2m][i1 ]*lp0[i2m][i1 ]*lp0[i2m][i1 ] +
                       d00[i2m][i1m]*lpp[i2m][i1m]*lpp[i2m][i1m];
        l0p[i2][i1] -= d00[i2m][i1p]*lp0[i2m][i1p]*lpm[i2m][i1p] +
                       d00[i2m][i1 ]*lpp[i2m][i1 ]*lp0[i2m][i1 ];
        lpm[i2][i1] -= d00[i2 ][i1m]*lp0[i2 ][i1m]*l0p[i2 ][i1m];
        lp0[i2][i1] -= d00[i2 ][i1m]*lpp[i2 ][i1m]*l0p[i2 ][i1m];
        if (l00[i2][i1]<=0.0f) 
          return null;
        d00[i2][i1] = 1.0f/l00[i2][i1];
      }
      l00[i2][i1] -= d00[i2 ][i1m]*l0p[i2 ][i1m]*l0p[i2 ][i1m] +
                     d00[i2m][i1 ]*lp0[i2m][i1 ]*lp0[i2m][i1 ] +
                     d00[i2m][i1m]*lpp[i2m][i1m]*lpp[i2m][i1m];
      l0p[i2][i1] -= d00[i2m][i1 ]*lpp[i2m][i1 ]*lp0[i2m][i1 ];
      lpm[i2][i1] -= d00[i2 ][i1m]*lp0[i2 ][i1m]*l0p[i2 ][i1m];
      lp0[i2][i1] -= d00[i2 ][i1m]*lpp[i2 ][i1m]*l0p[i2 ][i1m];
      if (l00[i2][i1]<=0.0f) 
        return null;
      d00[i2][i1] = 1.0f/l00[i2][i1];
    }

    // At this point, L is lower-triangular. Now make L have unit diagonal. 
    // This will enable application of inverse of L*D*L' without division.
    for (i2=0; i2<n2; ++i2) {
      for (i1=0; i1<n1; ++i1) {
        l0p[i2][i1] *= d00[i2][i1];
        lpm[i2][i1] *= d00[i2][i1];
        lp0[i2][i1] *= d00[i2][i1];
        lpp[i2][i1] *= d00[i2][i1];
      }
    }
    return l;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  private static void testFactor() {
    int n1 = 5;
    int n2 = 7;
    //float[][] x = zerofloat(n1,n2);
    //x[0][0] = x[n2-1][0] = x[0][n1-1] = x[n2-1][n1-1] = 1.0f;
    //x[2][2] = 1.0f;
    float[][] x = randfloat(n1,n2);
    float[][] y = randfloat(n1,n2);
    float[][] z = randfloat(n1,n2);
    float[][] w = randfloat(n1,n2);
    float theta = FLT_PI*2.0f/8.0f;
    float[][] d0 = fillfloat(1.0f,n1,n2);
    float[][] d1 = fillfloat(1.0f,n1,n2);
    float[][] v1 = fillfloat(sin(theta),n1,n2);
    LocalDiffusionTensors2 ldt = 
      new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    LocalDiffusionKernelX ldk = new LocalDiffusionKernelX();
    float[][][] s = ldk.getCoefficients(ldt);
    float[][] s00 = s[0], s0p = s[1], spm = s[2], sp0 = s[3], spp = s[4];
    int k1 = 1;
    int k2 = 1;
    s00[k2][k1] = 1.0f;
    s0p[k2][k1] = s0p[k2  ][k1-1] = 0.0f;
    spm[k2][k1] = spm[k2-1][k1+1] = 0.0f;
    sp0[k2][k1] = sp0[k2-1][k1  ] = 0.0f;
    spp[k2][k1] = spp[k2-1][k1-1] = 0.0f;
    LocalSpd9Filter lsf = new LocalSpd9Filter(s);
    //float[][] a = lsf.getMatrix();
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
    lsf.apply(x,y);
    lsf.applyApproximate(x,z);
    lsf.applyApproximateInverse(z,w);
    float[][] ldl = factorIC0(lsf,0.0f);
    float[][] v = factorMul(ldl,x);
    float[][] ez = sub(z,v);
    System.out.println("ez: error="+sum(abs(ez)));
    dump(z);
    dump(v);
    dump(ez);
    /*
    dump(x);
    dump(y);
    dump(z);
    dump(w);
    dump(sub(w,x));
    */
  }

  private static void testMatrix() {
    int n1 = 4;
    int n2 = 4;
    float[][] d0 = fillfloat(1.0f,n1,n2);
    float[][] d1 = fillfloat(1.0f,n1,n2);
    float[][] v1 = fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = 
      new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    LocalDiffusionKernelX ldk = new LocalDiffusionKernelX();
    float[][][] s = ldk.getCoefficients(ldt);
    add(0.1f,s[0],s[0]);
    LocalSpd9Filter lsf = new LocalSpd9Filter(s);
    float[][] a = lsf.getMatrix();
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
    float[][] x = randfloat(n1,n2);
    float[][] y = randfloat(n1,n2);
    lsf.apply(x,y);
    float[][] z = matrixMul(a,x);
    float[][] e = sub(z,y);
    System.out.println("error="+ sum(abs(e)));
    dump(y);
    dump(z);
    dump(e);
  }
  private static float[][] factorIC0(LocalSpd9Filter lsf, float bias) {
    float[][] a = lsf.getMatrix();
    int n = a.length;
    float scale = 1.0f+bias;
    for (int k=0; k<n; ++k) {
      a[k][k] = sqrt(a[k][k]*scale);
      for (int i=k+1; i<n; ++i) {
        if (a[k][i]!=0.0f)
          a[k][i] /= a[k][k];
      }
      for (int j=k+1; j<n; ++j) {
        for (int i=j; i<n; ++i) {
          if (a[j][i]!=0.0f)
            a[j][i] -= a[k][i]*a[k][j];
        }
      }
      for (int i=0; i<k; ++i)
        a[k][i] = 0.0f;
    }
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
    return a;
  }
  private static float[][] factorMul(float[][] a, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        float ti = 0.0f;
        for (int j2=0,j=0; j2<n2; ++j2) {
          for (int j1=0; j1<n1; ++j1,++j) {
            if (i<=j)
              ti += a[i][j]*x[j2][j1];
          }
        }
        for (int j2=0,j=0; j2<n2; ++j2) {
          for (int j1=0; j1<n1; ++j1,++j) {
            if (i<=j)
              y[j2][j1] += a[i][j]*ti;
          }
        }
      }
    }
    return y;
  }
  private static float[][] matrixMul(float[][] a, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        for (int j2=0,j=0; j2<n2; ++j2) {
          for (int j1=0; j1<n1; ++j1,++j) {
            y[i2][i1] += a[i][j]*x[j2][j1];
          }
        }
      }
    }
    return y;
  }

  public static void main(String[] args) {
    testFactor();
    //testMatrix();
  }

  ///////////////////////////////////////////////////////////////////////////
  // currently unused

  /**
   * Makes the specified matrix an M-matrix.
   */
  private static void mmatrix(float[][][] a) {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (a0p[i2][i1]>0.0f) {
          //trace("a0p>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += a0p[i2][i1];
          if (i1<n1-1)
            a00[i2][i1+1] += a0p[i2][i1];
          a0p[i2][i1] = 0.0f;
        }
        if (apm[i2][i1]>0.0f) {
          //trace("apm>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += apm[i2][i1];
          if (0<i1 && i2<n2-1)
            a00[i2+1][i1-1] += apm[i2][i1];
          apm[i2][i1] = 0.0f;
        }
        if (ap0[i2][i1]>0.0f) {
          //trace("ap0>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += ap0[i2][i1];
          if (i2<n2-1)
            a00[i2+1][i1] += ap0[i2][i1];
          ap0[i2][i1] = 0.0f;
        }
        if (app[i2][i1]>0.0f) {
          //trace("app>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += app[i2][i1];
          if (i1<n1-1 && i2<n2-1)
            a00[i2+1][i1+1] += app[i2][i1];
          app[i2][i1] = 0.0f;
        }
      }
    }
  }

  /**
   * Makes boundary elements of the specified matrix an M-matrix.
   */
  private static void mmatrixBoundaries(float[][][] a) {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (i1!=0 && i2!=0 && i1!=n1-1 && i2!=n2-1)
          continue;
        if (a0p[i2][i1]>0.0f) {
          //trace("a0p>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += a0p[i2][i1];
          if (i1<n1-1)
            a00[i2][i1+1] += a0p[i2][i1];
          a0p[i2][i1] = 0.0f;
        }
        if (apm[i2][i1]>0.0f) {
          //trace("apm>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += apm[i2][i1];
          if (0<i1 && i2<n2-1)
            a00[i2+1][i1-1] += apm[i2][i1];
          apm[i2][i1] = 0.0f;
        }
        if (ap0[i2][i1]>0.0f) {
          //trace("ap0>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += ap0[i2][i1];
          if (i2<n2-1)
            a00[i2+1][i1] += ap0[i2][i1];
          ap0[i2][i1] = 0.0f;
        }
        if (app[i2][i1]>0.0f) {
          //trace("app>0: i1="+i1+" i2="+i2);
          a00[i2][i1] += app[i2][i1];
          if (i1<n1-1 && i2<n2-1)
            a00[i2+1][i1+1] += app[i2][i1];
          app[i2][i1] = 0.0f;
        }
      }
    }
  }

  private static void fixPivot(float[][][] a, int i1, int i2) {
    trace("fixPivot: i1="+i1+" i2="+i2);
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    float sum = 0.0f;
    sum += abs(a0p[i2][i1]);
    sum += abs(apm[i2][i1]);
    sum += abs(ap0[i2][i1]);
    sum += abs(app[i2][i1]);
    if (0<i1)
      sum += abs(a0p[i2][i1-1]);
    if (0<i2) {
        sum += abs(ap0[i2-1][i1]);
      if (i1<n1-1)
        sum += abs(apm[i2-1][i1+1]);
      if (0<i1)
        sum += abs(app[i2-1][i1-1]);
    }
    a00[i2][i1] = sum;
  }
}
