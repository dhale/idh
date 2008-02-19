/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local symmetric positive-definite (SPD) filter with a 3-D 19-point stencil.
 * This filter is local in the sense that filter coefficients may differ for 
 * each sample. These coefficients form a 19-point stencil:
 * <pre><code>
 * y[i3][i2][i1] =             sm0m*xm0m +
 *                 smm0*xmm0 + sm00*xm00 + smp0*xmp0 +
 *                             sm0p*xm0p +
 *                 s0mm*x0mm + s00m*x00m + s0pm*x0pm +
 *                 s0m0*x0m0 + s000*x000 + s0p0*x0p0 +
 *                 s0mp*x0mp + s00p*x00p + s0pp*x0pp +
 *                             sp0m*xp0m +
 *                 spm0*xpm0 + sp00*xp00 + spp0*xpp0 +
 *                             sp0p*xp0p
 * </code></pre>
 * The suffixes m, 0, and p denote minus, zero, and plus, respectively.
 * For example sm0p is the coefficient of xm0p = x[i3-1][i2  ][i1+1].
 * <p>
 * For symmetric filters with constant coefficients, this stencil is 
 * symmetric about the central coefficient s000. For example, smmm = sppp.
 * Therefore, only ten of the nineteen coefficients need be specified.
 * <p>
 * For symmetric filters with variable coefficients, this stencil is 
 * <em>not</em> symmetric. That is, smmm[i3][i2][i1] does not equal 
 * sppp[i3][i2][i1]; rather, smmm[i3][i2][i1] = sppp[i3-1][i2-1][i1-1]. 
 * Still, only ten filter coefficients need be specified for each sample. 
 * If we choose those ten coefficients to be s000, s00p, s0pm, s0p0, s0pp, 
 * spm0, sp0m, sp00, sp0p, and spp0, then
 * <pre><code>
 * y[i3][i2][i1] = spp0[i3-1][i2-1][i1  ]*x[i3-1][i2-1][i1  ] +
 *                 sp0p[i3-1][i2  ][i1-1]*x[i3-1][i2  ][i1+1] +
 *                 sp00[i3-1][i2  ][i1  ]*x[i3-1][i2  ][i1  ] +
 *                 sp0m[i3-1][i2  ][i1+1]*x[i3-1][i2  ][i1-1] +
 *                 spm0[i3-1][i2+1][i1  ]*x[i3-1][i2+1][i1  ] +
 *                 s0pp[i3  ][i2-1][i1-1]*x[i3  ][i2-1][i1-1] +
 *                 s0p0[i3  ][i2-1][i1  ]*x[i3  ][i2-1][i1  ] +
 *                 s0pm[i3  ][i2-1][i1+1]*x[i3  ][i2-1][i1+1] +
 *                 s00p[i3  ][i2  ][i1-1]*x[i3  ][i2  ][i1-1] +
 *                 s000[i3  ][i2  ][i1  ]*x[i3  ][i2  ][i1  ] +
 *                 s00p[i3  ][i2  ][i1  ]*x[i3  ][i2  ][i1+1] +
 *                 s0pm[i3  ][i2  ][i1  ]*x[i3  ][i2+1][i1-1] +
 *                 s0p0[i3  ][i2  ][i1  ]*x[i3  ][i2+1][i1  ] +
 *                 s0pp[i3  ][i2  ][i1  ]*x[i3  ][i2+1][i1+1] +
 *                 spm0[i3  ][i2  ][i1  ]*x[i3+1][i2-1][i1  ] +
 *                 sp0m[i3  ][i2  ][i1  ]*x[i3+1][i2  ][i1-1] +
 *                 sp00[i3  ][i2  ][i1  ]*x[i3+1][i2  ][i1  ] +
 *                 sp0p[i3  ][i2  ][i1  ]*x[i3+1][i2  ][i1+1] +
 *                 spp0[i3  ][i2  ][i1  ]*x[i3+1][i2+1][i1  ]
 * </code></pre>
 * <p>
 * Becouse this filter is SPD, it may in theory be factored exactly with 
 * Cholesky decomposition. However, the factors seldom fit in a 19-point
 * stencil. Therefore, only approximate factors are typically computed 
 * using incomplete Cholesky (IC) decomposition. The factors may then be 
 * used to apply an approximate inverse of this filter. This approximate 
 * inverse is especially useful as a pre-conditioner in the method of 
 * conjugate gradients.
 * <p>
 * Unfortunately, IC decomposition may fail with non-positive pivots for 
 * filters that are not diagonally-dominant. To enable IC decomposition
 * to succeed, filter coefficients s000*(1+bias) may be used instead of 
 * s000. (Any bias is used only during IC decomposition; the specified 
 * s000 are not changed.) For filters known to be diagonally dominant, 
 * zero bias should be specified.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.02.19
 */
public class LocalSpd19Filter {

  /**
   * Constructs a filter with specified coefficients.
   * Any approximate inverse filter (when required) will be computed with an 
   * initial bias of zero.
   * @param s arrays[10][n3][n2][n1] of coefficients; by reference, not copy.
   *  The elements of the array s are 
   *  {s000,s00p,s0pm,s0p0,s0pp,spm0,sp0m,sp00,sp0p,spp0}.
   */
  public LocalSpd19Filter(float[][][][] s) {
    this(s,0.0);
  }

  /**
   * Constructs a filter with specified coefficients.
   * @param s arrays[10][n3][n2][n1] of coefficients; by reference, not copy.
   *  The elements of the array s are 
   *  {s000,s00p,s0pm,s0p0,s0pp,spm0,sp0m,sp00,sp0p,spp0}.
   * @param bias the initial non-negative amount by which to perturb the 
   *  coefficients s000 during computation of an approximate inverse filter.
   */
  public LocalSpd19Filter(float[][][][] s, double bias) {
    _s = s;
    _b = (float)bias;
  }

  /**
   * Applies this filter by computing y = A*x.
   * @param x input array. Must be distinct from y.
   * @param y output array. Must be distinct from x.
   */
  public void apply(float[][][] x, float[][][] y) {
    applyFilter(x,y);
  }

  /**
   * Applies this filter by computing y = A*x using an incomplete Cholesky 
   * decomposition of A. In effect, this method applies an approximation of 
   * this filter.
   * @param x input array. Must be distinct from y.
   * @param y output array. Must be distinct from x.
   */
  public void applyApproximate(float[][][] x, float[][][] y) {
    //applyFactors(x,y);
  }

  /**
   * Solves A*x = y using an incomplete Cholesky decomposition of A.
   * In effect, this method applies an approximate inverse of this filter.
   * @param y the input right-hand side array.
   * @param x the output solution array.
   */
  public void applyApproximateInverse(float[][][] y, float[][][] x) {
    //solveWithFactors(y,x);
  }

  /**
   * Gets the sparse matrix A equivalent to this filter.
   * Most elements in this matrix will be zero. For small numbers of
   * samples, it may be useful for visualization of matrix sparsity.
   * @return an array[n][n] representing A, where n = n1*n2*n3.
   */
  public float[][] getMatrix() {
    float[][][] s000 = _s[0];
    float[][][] s00p = _s[1];
    float[][][] s0pm = _s[2];
    float[][][] s0p0 = _s[3];
    float[][][] s0pp = _s[4];
    float[][][] spm0 = _s[5];
    float[][][] sp0m = _s[6];
    float[][][] sp00 = _s[7];
    float[][][] sp0p = _s[8];
    float[][][] spp0 = _s[9];
    int n1 = s000[0][0].length;
    int n2 = s000[0].length;
    int n3 = s000.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    int n = n1*n2*n3;
    float[][] a = new float[n][n];
    for (int i3=0,i=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1,++i) {
          int j = i+n1;
          int k = i+n1*n2;
                                     a[i   ][i] = s000[i3][i2][i1];
          if (i1<n1m)   a[i][i+1 ] = a[i+1 ][i] = s00p[i3][i2][i1];
          if (i2<n2m) {
            if (0<i1)   a[i][j-1 ] = a[j-1 ][i] = s0pm[i3][i2][i1];
                        a[i][j   ] = a[j   ][i] = s0p0[i3][i2][i1];
            if (i1<n1m) a[i][j+1 ] = a[j+1 ][i] = s0pp[i3][i2][i1];
          }
          if (i3<n3m) {
            if (0<i2)   a[i][k-n1] = a[k-n1][i] = spm0[i3][i2][i1];
            if (0<i1)   a[i][k-1 ] = a[k-1 ][i] = sp0m[i3][i2][i1];
                        a[i][k   ] = a[k   ][i] = sp00[i3][i2][i1];
            if (i1<n1m) a[i][k+1 ] = a[k+1 ][i] = sp0p[i3][i2][i1];
            if (i2<n2m) a[i][k+n1] = a[k+n1][i] = spp0[i3][i2][i1];
          }
        }
      }
    }
    return a;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][][] _s; // specified SPD filter coefficients
  private float[][][][] _l; // coefficients of IC(0) decomposition
  private float _b; // initial bias for IC(0) decomposition.

  /**
   * Makes the IC(0) factors, if not already made.
   */
  private void ensureFactors() {
    //if (_l==null)
    //  _l = factorIC0(_s,_b);
    Check.state(_l!=null,"incomplete Cholesky decomposition successful");
  }

  /**
   * Computes y = A*x.
   */
  private void applyFilter(float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float[][][] s000 = _s[0];
    float[][][] s00p = _s[1];
    float[][][] s0pm = _s[2];
    float[][][] s0p0 = _s[3];
    float[][][] s0pp = _s[4];
    float[][][] spm0 = _s[5];
    float[][][] sp0m = _s[6];
    float[][][] sp00 = _s[7];
    float[][][] sp0p = _s[8];
    float[][][] spp0 = _s[9];
    int i1,i2,i3,i1m,i1p;
    i1 = n1m;
    i2 = n2m;
    i3 = n3m;
    y[i3  ][i2  ][i1  ]  = s000[i3][i2][i1]*x[i3  ][i2  ][i1  ];
    for (i1=n1m-1; i1>=0; --i1) {
      y[i3  ][i2  ][i1  ]  = s000[i3][i2][i1]*x[i3  ][i2  ][i1  ];
      y[i3  ][i2  ][i1  ] += s00p[i3][i2][i1]*x[i3  ][i2  ][i1+1];
      y[i3  ][i2  ][i1+1] += s00p[i3][i2][i1]*x[i3  ][i2  ][i1  ];
    }


    for (i3=n3m-1; i3>=0; --i3) {
      for (i2=n2m-1; i2>=1; --i2) {
        float[] s000i = s000[i3][i2];
        float[] s00pi = s00p[i3][i2];
        float[] s0pmi = s0pm[i3][i2];
        float[] s0p0i = s0p0[i3][i2];
        float[] s0ppi = s0pp[i3][i2];
        float[] spm0i = spm0[i3][i2];
        float[] sp0mi = sp0m[i3][i2];
        float[] sp00i = sp00[i3][i2];
        float[] sp0pi = sp0p[i3][i2];
        float[] spp0i = spp0[i3][i2];
        float[] y00 = y[i3  ][i2  ];
        float[] y0p = y[i3  ][i2+1];
        float[] ypm = y[i3+1][i2-1];
        float[] yp0 = y[i3+1][i2  ];
        float[] ypp = y[i3+1][i2+1];
        float[] x00 = x[i3  ][i2  ];
        float[] x0p = x[i3  ][i2+1];
        float[] xpm = x[i3+1][i2-1];
        float[] xp0 = x[i3+1][i2  ];
        float[] xpp = x[i3+1][i2+1];
        for (i1=n1m-1,i1m=i1-1,i1p=i1+1; i1>=1; --i1,--i1m,--i1p) {
          y00[i1 ]  = s000i[i1]*x00[i1 ];
          y00[i1 ] += s00pi[i1]*x00[i1p];
          y00[i1p] += s00pi[i1]*x00[i1 ];
          y00[i1 ] += s0pmi[i1]*x0p[i1m];
          y0p[i1m] += s0pmi[i1]*x00[i1 ];
          y00[i1 ] += s0p0i[i1]*x0p[i1 ];
          y0p[i1 ] += s0p0i[i1]*x00[i1 ];
          y00[i1 ] += s0ppi[i1]*x0p[i1p];
          y0p[i1p] += s0ppi[i1]*x00[i1 ];
          y00[i1 ] += spm0i[i1]*xpm[i1 ];
          ypm[i1 ] += spm0i[i1]*x00[i1 ];
          y00[i1 ] += sp0mi[i1]*xp0[i1m];
          yp0[i1m] += sp0mi[i1]*x00[i1 ];
          y00[i1 ] += sp00i[i1]*xp0[i1 ];
          yp0[i1 ] += sp00i[i1]*x00[i1 ];
          y00[i1 ] += sp0pi[i1]*xp0[i1p];
          yp0[i1p] += sp0pi[i1]*x00[i1 ];
          y00[i1 ] += spp0i[i1]*xpp[i1 ];
          ypp[i1 ] += spp0i[i1]*x00[i1 ];
        }
      }
    }
  }

  /**
   * Solves L*D*L'*x = b.
   */
/*
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
    Array.zero(x);
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
*/
  
  /**
   * Computes y = L*D*L'*x. For testing, only.
   */
/*
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
*/

  /**
   * Factors this filter with incomplete Cholesky decomposition IC(0).
   * The approximate factorization is A ~ L*D*L', where L is a 
   * unit-lower-triangular matrix, and D is a diagonal matrix.
   * @param bias the amount by which to perturb the diagonal of A.
   * @return The coefficients {d00,l0p,lpm,lp0,lpp}. The elements in the
   *  array d00 are the inverse of the elements of the diagonal matrix D.
   */
/*
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
    float[][][] l = Array.copy(a);
    if (bias>0.0f)
      Array.mul(1.0f+bias,l[0],l[0]);
    float[][] l00 = l[0], l0p = l[1], lpm = l[2], lp0 = l[3], lpp = l[4];
    float[][] d00 = l00; // will contain inverse of diagonal matrix D

    // Incomplete Cholesky decomposition, in-place.
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    int i1 = 0;
    int i2 = 0;
    d00[i2][i1] = 1.0f/l00[i2][i1];
    for (i1=1; i1<n1; ++i1) {
      l00[i2][i1] -= d00[i2  ][i1-1]*l0p[i2  ][i1-1]*l0p[i2  ][i1-1];
      lpm[i2][i1] -= d00[i2  ][i1-1]*lp0[i2  ][i1-1]*l0p[i2  ][i1-1];
      lp0[i2][i1] -= d00[i2  ][i1-1]*lpp[i2  ][i1-1]*l0p[i2  ][i1-1];
      if (l00[i2][i1]<=0.0f) 
        return null;
      d00[i2][i1] = 1.0f/l00[i2][i1];
    }
    for (i2=1; i2<n2; ++i2) {
      i1 = 0;
      l00[i2][i1] -= d00[i2-1][i1+1]*lpm[i2-1][i1+1]*lpm[i2-1][i1+1] +
                     d00[i2-1][i1  ]*lp0[i2-1][i1  ]*lp0[i2-1][i1  ];
      l0p[i2][i1] -= d00[i2-1][i1  ]*lpp[i2-1][i1  ]*lp0[i2-1][i1  ];
      if (l00[i2][i1]<=0.0f) 
        return null;
      d00[i2][i1] = 1.0f/l00[i2][i1];
      for (i1=1; i1<n1-1; ++i1) {
        l00[i2][i1] -= d00[i2  ][i1-1]*l0p[i2  ][i1-1]*l0p[i2  ][i1-1] +
                       d00[i2-1][i1+1]*lpm[i2-1][i1+1]*lpm[i2-1][i1+1] +
                       d00[i2-1][i1  ]*lp0[i2-1][i1  ]*lp0[i2-1][i1  ] +
                       d00[i2-1][i1-1]*lpp[i2-1][i1-1]*lpp[i2-1][i1-1];
        l0p[i2][i1] -= d00[i2-1][i1+1]*lp0[i2-1][i1+1]*lpm[i2-1][i1+1] +
                       d00[i2-1][i1  ]*lpp[i2-1][i1  ]*lp0[i2-1][i1  ];
        lpm[i2][i1] -= d00[i2  ][i1-1]*lp0[i2  ][i1-1]*l0p[i2  ][i1-1];
        lp0[i2][i1] -= d00[i2  ][i1-1]*lpp[i2  ][i1-1]*l0p[i2  ][i1-1];
        if (l00[i2][i1]<=0.0f) 
          return null;
        d00[i2][i1] = 1.0f/l00[i2][i1];
      }
      l00[i2][i1] -= d00[i2  ][i1-1]*l0p[i2  ][i1-1]*l0p[i2  ][i1-1] +
                     d00[i2-1][i1  ]*lp0[i2-1][i1  ]*lp0[i2-1][i1  ] +
                     d00[i2-1][i1-1]*lpp[i2-1][i1-1]*lpp[i2-1][i1-1];
      l0p[i2][i1] -= d00[i2-1][i1  ]*lpp[i2-1][i1  ]*lp0[i2-1][i1  ];
      lpm[i2][i1] -= d00[i2  ][i1-1]*lp0[i2  ][i1-1]*l0p[i2  ][i1-1];
      lp0[i2][i1] -= d00[i2  ][i1-1]*lpp[i2  ][i1-1]*l0p[i2  ][i1-1];
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
*/

  ///////////////////////////////////////////////////////////////////////////
  // testing

/*
  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  private static void testFactor() {
    int n1 = 5;
    int n2 = 5;
    //float[][] x = Array.zerofloat(n1,n2);
    //x[0][0] = x[n2-1][0] = x[0][n1-1] = x[n2-1][n1-1] = 1.0f;
    //x[2][2] = 1.0f;
    float[][] x = Array.randfloat(n1,n2);
    float[][] y = Array.randfloat(n1,n2);
    float[][] z = Array.randfloat(n1,n2);
    float[][] w = Array.randfloat(n1,n2);
    float theta = FLT_PI*2.0f/8.0f;
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sin(theta),n1,n2);
    LocalDiffusionTensors2 ldt = 
      new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
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
    float[][] a = lsf.getMatrix();
    edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
    lsf.apply(x,y);
    lsf.applyApproximate(x,z);
    lsf.applyApproximateInverse(z,w);
    Array.dump(x);
    Array.dump(y);
    Array.dump(z);
    Array.dump(w);
    Array.dump(Array.sub(w,x));
  }

  private static void testMatrix() {
    int n1 = 4;
    int n2 = 4;
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = 
      new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    float[][][] s = ldk.getCoefficients(ldt);
    Array.add(0.1f,s[0],s[0]);
    LocalSpd9Filter lsf = new LocalSpd9Filter(s);
    float[][] a = lsf.getMatrix();
    edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
  }

  public static void main(String[] args) {
    testFactor();
    //testMatrix();
  }
*/
}
