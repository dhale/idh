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
 * Local symmetric positive-definite (SPD) filter with a 3-D 27-point stencil.
 * This filter is local in the sense that filter coefficients may differ for 
 * each sample. These coefficients form a 27-point stencil:
 * <pre><code>
 * y[i3][i2][i1] = smmm*xmmm + sm0m*xm0m + smpm*xmpm +
 *                 smm0*xmm0 + sm00*xm00 + smp0*xmp0 +
 *                 smmp*xmmp + sm0p*xm0p + smpp*xmpp +
 *
 *                 s0mm*x0mm + s00m*x00m + s0pm*x0pm +
 *                 s0m0*x0m0 + s000*x000 + s0p0*x0p0 +
 *                 s0mp*x0mp + s00p*x00p + s0pp*x0pp +
 *
 *                 spmm*xpmm + sp0m*xp0m + sppm*xppm +
 *                 spm0*xpm0 + sp00*xp00 + spp0*xpp0 +
 *                 spmp*xpmp + sp0p*xp0p + sppp*xppp
 * </code></pre>
 * The suffixes m, 0, and p denote minus, zero, and plus, respectively.
 * For example sm0p is the coefficient of xm0p = x[i3-1][i2  ][i1+1].
 * <p>
 * For symmetric filters with constant coefficients, this stencil is 
 * symmetric about the central coefficient s000. For example, smmm = sppp.
 * Therefore, only 14 of the 27 coefficients need be specified.
 * <p>
 * For symmetric filters with variable coefficients, this stencil is 
 * <em>not</em> symmetric. That is, smmm[i3][i2][i1] does not equal 
 * sppp[i3][i2][i1]; rather, smmm[i3][i2][i1] = sppp[i3-1][i2-1][i1-1]. 
 * Still, only 14 filter coefficients need be specified for each sample. 
 * If we choose those 14 coefficients to be s000, s00p, s0pm, s0p0, s0pp, 
 * spmm, spm0, spmp, sp0m, sp00, sp0p, sppm, spp0 and sppp, then
 * <pre><code>
 * y[i3][i2][i1] = sppp[i3-1][i2-1][i1-1]*x[i3-1][i2-1][i1-1] +
 *                 spp0[i3-1][i2-1][i1  ]*x[i3-1][i2-1][i1  ] +
 *                 sppm[i3-1][i2-1][i1+1]*x[i3-1][i2-1][i1+1] +
 *                 sp0p[i3-1][i2  ][i1-1]*x[i3-1][i2  ][i1-1] +
 *                 sp00[i3-1][i2  ][i1  ]*x[i3-1][i2  ][i1  ] +
 *                 sp0m[i3-1][i2  ][i1+1]*x[i3-1][i2  ][i1+1] +
 *                 spmp[i3-1][i2+1][i1-1]*x[i3-1][i2+1][i1-1] +
 *                 spm0[i3-1][i2+1][i1  ]*x[i3-1][i2+1][i1  ] +
 *                 spmm[i3-1][i2+1][i1+1]*x[i3-1][i2+1][i1+1] +
 *
 *                 s0pp[i3  ][i2-1][i1-1]*x[i3  ][i2-1][i1-1] +
 *                 s0p0[i3  ][i2-1][i1  ]*x[i3  ][i2-1][i1  ] +
 *                 s0pm[i3  ][i2-1][i1+1]*x[i3  ][i2-1][i1+1] +
 *                 s00p[i3  ][i2  ][i1-1]*x[i3  ][i2  ][i1-1] +
 *                 s000[i3  ][i2  ][i1  ]*x[i3  ][i2  ][i1  ] +
 *                 s00p[i3  ][i2  ][i1  ]*x[i3  ][i2  ][i1+1] +
 *                 s0pm[i3  ][i2  ][i1  ]*x[i3  ][i2+1][i1-1] +
 *                 s0p0[i3  ][i2  ][i1  ]*x[i3  ][i2+1][i1  ] +
 *                 s0pp[i3  ][i2  ][i1  ]*x[i3  ][i2+1][i1+1] +
 *
 *                 spmm[i3  ][i2  ][i1  ]*x[i3+1][i2-1][i1-1] +
 *                 spm0[i3  ][i2  ][i1  ]*x[i3+1][i2-1][i1  ] +
 *                 spmp[i3  ][i2  ][i1  ]*x[i3+1][i2-1][i1+1] +
 *                 sp0m[i3  ][i2  ][i1  ]*x[i3+1][i2  ][i1-1] +
 *                 sp00[i3  ][i2  ][i1  ]*x[i3+1][i2  ][i1  ] +
 *                 sp0p[i3  ][i2  ][i1  ]*x[i3+1][i2  ][i1+1] +
 *                 sppm[i3  ][i2  ][i1  ]*x[i3+1][i2+1][i1-1] +
 *                 spp0[i3  ][i2  ][i1  ]*x[i3+1][i2+1][i1  ] +
 *                 sppp[i3  ][i2  ][i1  ]*x[i3+1][i2+1][i1+1]
 * </code></pre>
 * <p>
 * Becouse this filter is SPD, it may in theory be factored exactly with 
 * Cholesky decomposition. However, the factors seldom fit in a 27-point
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
 * @version 2008.05.08
 */
public class LocalSpd27Filter {

  /**
   * Constructs a filter with specified coefficients.
   * Any approximate inverse filter (when required) will be computed with an 
   * initial bias of zero.
   * @param s arrays[14][n3][n2][n1] of coefficients; by reference, not copy.
   *  The elements of the array s are 
   *  {s000,s00p,s0pm,s0p0,s0pp,spm0,sp0m,sp00,sp0p,spp0}.
   */
  public LocalSpd27Filter(float[][][][] s) {
    this(s,0.0);
  }

  /**
   * Constructs a filter with specified coefficients.
   * @param s arrays[14][n3][n2][n1] of coefficients; by reference, not copy.
   *  The elements of the array s are 
   *  {s000,s00p,s0pm,s0p0,s0pp,spmm,spm0,spmp,sp0m,sp00,sp0p,sppm,spp0,sppp}.
   * @param bias the initial non-negative amount by which to perturb the 
   *  coefficients s000 during computation of an approximate inverse filter.
   */
  public LocalSpd27Filter(float[][][][] s, double bias) {
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
    applyFactors(x,y);
  }

  /**
   * Solves A*x = y using an incomplete Cholesky decomposition of A.
   * In effect, this method applies an approximate inverse of this filter.
   * @param y the input right-hand side array.
   * @param x the output solution array.
   */
  public void applyApproximateInverse(float[][][] y, float[][][] x) {
    solveWithFactors(y,x);
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
    float[][][] spmm = _s[5];
    float[][][] spm0 = _s[6];
    float[][][] spmp = _s[7];
    float[][][] sp0m = _s[8];
    float[][][] sp00 = _s[9];
    float[][][] sp0p = _s[10];
    float[][][] sppm = _s[11];
    float[][][] spp0 = _s[12];
    float[][][] sppp = _s[13];
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
          int k0 = i+n1*n2;
          int km = k0-n1;
          int kp = k0+n1;
                                       a[i   ][i] = s000[i3][i2][i1];
          if (i1<n1m)     a[i][i+1 ] = a[i+1 ][i] = s00p[i3][i2][i1];
          if (i2<n2m) {
            if (0<i1)     a[i][j-1 ] = a[j-1 ][i] = s0pm[i3][i2][i1];
                          a[i][j   ] = a[j   ][i] = s0p0[i3][i2][i1];
            if (i1<n1m)   a[i][j+1 ] = a[j+1 ][i] = s0pp[i3][i2][i1];
          }
          if (i3<n3m) {
            if (0<i2) {
              if (0<i1)   a[i][km-1] = a[km-1][i] = spmm[i3][i2][i1];
                          a[i][km  ] = a[km  ][i] = spm0[i3][i2][i1];
              if (i1<n1m) a[i][km+1] = a[km+1][i] = spmp[i3][i2][i1]; 
            }
            if (0<i1)     a[i][k0-1] = a[k0-1][i] = sp0m[i3][i2][i1];
                          a[i][k0  ] = a[k0  ][i] = sp00[i3][i2][i1];
            if (i1<n1m)   a[i][k0+1] = a[k0+1][i] = sp0p[i3][i2][i1];
            if (i2<n2m) {
              if (0<i1)   a[i][kp-1] = a[kp-1][i] = sppm[i3][i2][i1];
                          a[i][kp  ] = a[kp  ][i] = spp0[i3][i2][i1];
              if (i1<n1m) a[i][kp+1] = a[kp+1][i] = sppp[i3][i2][i1];
            }
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
    if (_l==null)
      _l = factorIC0(_s,_b);
    Check.state(_l!=null,"incomplete Cholesky decomposition successful");
  }

  /**
   * Computes y = A*x.
   */
  private void applyFilter(float[][][] x, float[][][] y) {
    // 27-point stencil (elements marked X are determined by symmetry):
    //  0     X     0  |   X     X    s0pm  |  spmm  sp0m  sppm
    //  X     X     X  |   X    s000  s0p0  |  spm0  sp00  spp0
    //  0     X     0  |   X    s00p  s0pp  |  spmp  sp0p  sppp
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
    float[][][] spmm = _s[5];
    float[][][] spm0 = _s[6];
    float[][][] spmp = _s[7];
    float[][][] sp0m = _s[8];
    float[][][] sp00 = _s[9];
    float[][][] sp0p = _s[10];
    float[][][] sppm = _s[11];
    float[][][] spp0 = _s[12];
    float[][][] sppp = _s[13];
    int i1,i2,i3,i1m,i2m,i1p,i2p,i3p;
    for (i3=n3m,i3p=i3+1; i3>=0; --i3,--i3p) {
      for (i2=n2m,i2m=i2-1,i2p=i2+1; i2>=0; --i2,--i2m,--i2p) {
        float[] s000i = s000[i3][i2];
        float[] s00pi = s00p[i3][i2];
        float[] s0pmi = s0pm[i3][i2];
        float[] s0p0i = s0p0[i3][i2];
        float[] s0ppi = s0pp[i3][i2];
        float[] spmmi = spmm[i3][i2];
        float[] spm0i = spm0[i3][i2];
        float[] spmpi = spmp[i3][i2];
        float[] sp0mi = sp0m[i3][i2];
        float[] sp00i = sp00[i3][i2];
        float[] sp0pi = sp0p[i3][i2];
        float[] sppmi = sppm[i3][i2];
        float[] spp0i = spp0[i3][i2];
        float[] spppi = sppp[i3][i2];
        if (0<=i2m && i2p<n2 && i3p<n3) {
          float[] x00 = x[i3 ][i2 ];
          float[] x0p = x[i3 ][i2p];
          float[] xpm = x[i3p][i2m];
          float[] xp0 = x[i3p][i2 ];
          float[] xpp = x[i3p][i2p];
          float[] y00 = y[i3 ][i2 ];
          float[] y0p = y[i3 ][i2p];
          float[] ypm = y[i3p][i2m];
          float[] yp0 = y[i3p][i2 ];
          float[] ypp = y[i3p][i2p];
          i1 = n1m;
          i1m = i1-1;
          y00[i1 ]  = s000i[i1]*x00[i1 ];
          y00[i1 ] += s0pmi[i1]*x0p[i1m];
          y0p[i1m] += s0pmi[i1]*x00[i1 ];
          y00[i1 ] += s0p0i[i1]*x0p[i1 ];
          y0p[i1 ] += s0p0i[i1]*x00[i1 ];
          y00[i1 ] += spmmi[i1]*xpm[i1m];
          ypm[i1m] += spmmi[i1]*x00[i1 ];
          y00[i1 ] += spm0i[i1]*xpm[i1 ];
          ypm[i1 ] += spm0i[i1]*x00[i1 ];
          y00[i1 ] += sp0mi[i1]*xp0[i1m];
          yp0[i1m] += sp0mi[i1]*x00[i1 ];
          y00[i1 ] += sp00i[i1]*xp0[i1 ];
          yp0[i1 ] += sp00i[i1]*x00[i1 ];
          y00[i1 ] += sppmi[i1]*xpp[i1m];
          ypp[i1m] += sppmi[i1]*x00[i1 ];
          y00[i1 ] += spp0i[i1]*xpp[i1 ];
          ypp[i1 ] += spp0i[i1]*x00[i1 ];
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
            y00[i1 ] += spmmi[i1]*xpm[i1m];
            ypm[i1m] += spmmi[i1]*x00[i1 ];
            y00[i1 ] += spm0i[i1]*xpm[i1 ];
            ypm[i1 ] += spm0i[i1]*x00[i1 ];
            y00[i1 ] += spmpi[i1]*xpm[i1p];
            ypm[i1p] += spmpi[i1]*x00[i1 ];
            y00[i1 ] += sp0mi[i1]*xp0[i1m];
            yp0[i1m] += sp0mi[i1]*x00[i1 ];
            y00[i1 ] += sp00i[i1]*xp0[i1 ];
            yp0[i1 ] += sp00i[i1]*x00[i1 ];
            y00[i1 ] += sp0pi[i1]*xp0[i1p];
            yp0[i1p] += sp0pi[i1]*x00[i1 ];
            y00[i1 ] += sppmi[i1]*xpp[i1m];
            ypp[i1m] += sppmi[i1]*x00[i1 ];
            y00[i1 ] += spp0i[i1]*xpp[i1 ];
            ypp[i1 ] += spp0i[i1]*x00[i1 ];
            y00[i1 ] += spppi[i1]*xpp[i1p];
            ypp[i1p] += spppi[i1]*x00[i1 ];
          }
          y00[i1 ]  = s000i[i1]*x00[i1 ];
          y00[i1 ] += s00pi[i1]*x00[i1p];
          y00[i1p] += s00pi[i1]*x00[i1 ];
          y00[i1 ] += s0p0i[i1]*x0p[i1 ];
          y0p[i1 ] += s0p0i[i1]*x00[i1 ];
          y00[i1 ] += s0ppi[i1]*x0p[i1p];
          y0p[i1p] += s0ppi[i1]*x00[i1 ];
          y00[i1 ] += spm0i[i1]*xpm[i1 ];
          ypm[i1 ] += spm0i[i1]*x00[i1 ];
          y00[i1 ] += spmpi[i1]*xpm[i1p];
          ypm[i1p] += spmpi[i1]*x00[i1 ];
          y00[i1 ] += sp00i[i1]*xp0[i1 ];
          yp0[i1 ] += sp00i[i1]*x00[i1 ];
          y00[i1 ] += sp0pi[i1]*xp0[i1p];
          yp0[i1p] += sp0pi[i1]*x00[i1 ];
          y00[i1 ] += spp0i[i1]*xpp[i1 ];
          ypp[i1 ] += spp0i[i1]*x00[i1 ];
          y00[i1 ] += spppi[i1]*xpp[i1p];
          ypp[i1p] += spppi[i1]*x00[i1 ];
        } else {
          for (i1=n1m,i1m=i1-1,i1p=i1+1; i1>=0; --i1,--i1m,--i1p) {
            y[i3 ][i2 ][i1 ]  = s000i[i1]*x[i3 ][i2 ][i1 ];
            if (i1p<n1) {
              y[i3 ][i2 ][i1 ] += s00pi[i1]*x[i3 ][i2 ][i1p];
              y[i3 ][i2 ][i1p] += s00pi[i1]*x[i3 ][i2 ][i1 ];
            }
            if (i2p<n2) {
              if (0<=i1m) {
                y[i3 ][i2 ][i1 ] += s0pmi[i1]*x[i3 ][i2p][i1m];
                y[i3 ][i2p][i1m] += s0pmi[i1]*x[i3 ][i2 ][i1 ];
              }
              y[i3 ][i2 ][i1 ] += s0p0i[i1]*x[i3 ][i2p][i1 ];
              y[i3 ][i2p][i1 ] += s0p0i[i1]*x[i3 ][i2 ][i1 ];
              if (i1p<n1) {
                y[i3 ][i2 ][i1 ] += s0ppi[i1]*x[i3 ][i2p][i1p];
                y[i3 ][i2p][i1p] += s0ppi[i1]*x[i3 ][i2 ][i1 ];
              }
            }
            if (i3p<n3) {
              if (0<=i2m) {
                if (0<=i1m) {
                  y[i3 ][i2 ][i1 ] += spmmi[i1]*x[i3p][i2m][i1m];
                  y[i3p][i2m][i1m] += spmmi[i1]*x[i3 ][i2 ][i1 ];
                }
                y[i3 ][i2 ][i1 ] += spm0i[i1]*x[i3p][i2m][i1 ];
                y[i3p][i2m][i1 ] += spm0i[i1]*x[i3 ][i2 ][i1 ];
                if (i1p<n1) {
                  y[i3 ][i2 ][i1 ] += spmpi[i1]*x[i3p][i2m][i1p];
                  y[i3p][i2m][i1p] += spmpi[i1]*x[i3 ][i2 ][i1 ];
                }
              }
              if (0<=i1m) {
                y[i3 ][i2 ][i1 ] += sp0mi[i1]*x[i3p][i2 ][i1m];
                y[i3p][i2 ][i1m] += sp0mi[i1]*x[i3 ][i2 ][i1 ];
              }
              y[i3 ][i2 ][i1 ] += sp00i[i1]*x[i3p][i2 ][i1 ];
              y[i3p][i2 ][i1 ] += sp00i[i1]*x[i3 ][i2 ][i1 ];
              if (i1p<n1) {
                y[i3 ][i2 ][i1 ] += sp0pi[i1]*x[i3p][i2 ][i1p];
                y[i3p][i2 ][i1p] += sp0pi[i1]*x[i3 ][i2 ][i1 ];
              }
              if (i2p<n2) {
                if (0<=i1m) {
                  y[i3 ][i2 ][i1 ] += sppmi[i1]*x[i3p][i2p][i1m];
                  y[i3p][i2p][i1m] += sppmi[i1]*x[i3 ][i2 ][i1 ];
                }
                y[i3 ][i2 ][i1 ] += spp0i[i1]*x[i3p][i2p][i1 ];
                y[i3p][i2p][i1 ] += spp0i[i1]*x[i3 ][i2 ][i1 ];
                if (i1p<n1) {
                  y[i3 ][i2 ][i1 ] += spppi[i1]*x[i3p][i2p][i1p];
                  y[i3p][i2p][i1p] += spppi[i1]*x[i3 ][i2 ][i1 ];
                }
              }
            }
          }
        }
      }
    }
  }

  /**
   * Solves L*D*L'*x = b.
   */
  private void solveWithFactors(float[][][] b, float[][][] x) {
    ensureFactors();
    int n1 = b[0][0].length;
    int n2 = b[0].length;
    int n3 = b.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float[][][] d000 = _l[0];
    float[][][] l00p = _l[1];
    float[][][] l0pm = _l[2];
    float[][][] l0p0 = _l[3];
    float[][][] l0pp = _l[4];
    float[][][] lpmm = _l[5];
    float[][][] lpm0 = _l[6];
    float[][][] lpmp = _l[7];
    float[][][] lp0m = _l[8];
    float[][][] lp00 = _l[9];
    float[][][] lp0p = _l[10];
    float[][][] lppm = _l[11];
    float[][][] lpp0 = _l[12];
    float[][][] lppp = _l[13];

    int i1,i2,i3,i1m,i2m,i1p,i2p,i3p;
    float xi;

    // Solve L*z = b.
    Array.zero(x);
    for (i3=0,i3p=i3+1; i3<n3; ++i3,++i3p) {
      for (i2=0,i2m=i2-1,i2p=i2+1; i2<n2; ++i2,++i2m,++i2p) {
        if (0<=i2m && i2p<n2 && i3p<n3) {
          float[] b00 = b[i3 ][i2 ];
          float[] x0p = x[i3 ][i2p];
          float[] x00 = x[i3 ][i2 ];
          float[] x0m = x[i3 ][i2m];
          float[] xpm = x[i3p][i2m];
          float[] xp0 = x[i3p][i2 ];
          float[] xpp = x[i3p][i2p];
          float[] l00pi = l00p[i3][i2];
          float[] l0pmi = l0pm[i3][i2];
          float[] l0p0i = l0p0[i3][i2];
          float[] l0ppi = l0pp[i3][i2];
          float[] lpmmi = lpmm[i3][i2];
          float[] lpm0i = lpm0[i3][i2];
          float[] lpmpi = lpmp[i3][i2];
          float[] lp0mi = lp0m[i3][i2];
          float[] lp00i = lp00[i3][i2];
          float[] lp0pi = lp0p[i3][i2];
          float[] lppmi = lppm[i3][i2];
          float[] lpp0i = lpp0[i3][i2];
          float[] lpppi = lppp[i3][i2];
          i1 = 0;
          i1p = i1+1;
          x00[i1] += b00[i1];
          xi = x00[i1];
          x00[i1p] -= l00pi[i1]*xi;
          x0p[i1 ] -= l0p0i[i1]*xi;
          x0p[i1p] -= l0ppi[i1]*xi;
          xpm[i1 ] -= lpm0i[i1]*xi;
          xpm[i1p] -= lpmpi[i1]*xi;
          xp0[i1 ] -= lp00i[i1]*xi;
          xp0[i1p] -= lp0pi[i1]*xi;
          xpp[i1 ] -= lpp0i[i1]*xi;
          xpp[i1p] -= lpppi[i1]*xi;
          for (i1=1,i1m=i1-1,i1p=i1+1; i1<n1m; ++i1,++i1m,++i1p) {
            x00[i1] += b00[i1];
            xi = x00[i1];
            x00[i1p] -= l00pi[i1]*xi;
            x0p[i1m] -= l0pmi[i1]*xi;
            x0p[i1 ] -= l0p0i[i1]*xi;
            x0p[i1p] -= l0ppi[i1]*xi;
            xpm[i1m] -= lpmmi[i1]*xi;
            xpm[i1 ] -= lpm0i[i1]*xi;
            xpm[i1p] -= lpmpi[i1]*xi;
            xp0[i1m] -= lp0mi[i1]*xi;
            xp0[i1 ] -= lp00i[i1]*xi;
            xp0[i1p] -= lp0pi[i1]*xi;
            xpp[i1m] -= lppmi[i1]*xi;
            xpp[i1 ] -= lpp0i[i1]*xi;
            xpp[i1p] -= lpppi[i1]*xi;
          }
          x00[i1] += b00[i1];
          xi = x00[i1];
          x0p[i1m] -= l0pmi[i1]*xi;
          x0p[i1 ] -= l0p0i[i1]*xi;
          xpm[i1m] -= lpmmi[i1]*xi;
          xpm[i1 ] -= lpm0i[i1]*xi;
          xp0[i1m] -= lp0mi[i1]*xi;
          xp0[i1 ] -= lp00i[i1]*xi;
          xpp[i1m] -= lppmi[i1]*xi;
          xpp[i1 ] -= lpp0i[i1]*xi;
        } else {
          for (i1=0,i1m=i1-1,i1p=i1+1; i1<n1; ++i1,++i1m,++i1p) {
            x[i3][i2][i1] += b[i3][i2][i1];
            xi = x[i3][i2][i1];
            if (i1p<n1)     x[i3 ][i2 ][i1p] -= l00p[i3][i2][i1]*xi;
            if (i2p<n2) {
              if (0<=i1m)   x[i3 ][i2p][i1m] -= l0pm[i3][i2][i1]*xi;
                            x[i3 ][i2p][i1 ] -= l0p0[i3][i2][i1]*xi;
              if (i1p<n1)   x[i3 ][i2p][i1p] -= l0pp[i3][i2][i1]*xi;
            }
            if (i3p<n3) {
              if (0<=i2m) {
                if (0<=i1m) x[i3p][i2m][i1m] -= lpmm[i3][i2][i1]*xi;
                            x[i3p][i2m][i1 ] -= lpm0[i3][i2][i1]*xi;
                if (i1p<n1) x[i3p][i2m][i1p] -= lpmp[i3][i2][i1]*xi;
              }
              if (0<=i1m)   x[i3p][i2 ][i1m] -= lp0m[i3][i2][i1]*xi;
                            x[i3p][i2 ][i1 ] -= lp00[i3][i2][i1]*xi;
              if (i1p<n1)   x[i3p][i2 ][i1p] -= lp0p[i3][i2][i1]*xi;
              if (i2p<n2) {
                if (0<=i1m) x[i3p][i2p][i1m] -= lppm[i3][i2][i1]*xi;
                            x[i3p][i2p][i1 ] -= lpp0[i3][i2][i1]*xi;
                if (i1p<n1) x[i3p][i2p][i1p] -= lppp[i3][i2][i1]*xi;
              }
            }
          }
        }
      }
    }

    // Solve D*y = z and L'*x = y.
    for (i3=n3m,i3p=i3+1; i3>=0; --i3,--i3p) {
      for (i2=n2m,i2m=i2-1,i2p=i2+1; i2>=0; --i2,--i2m,--i2p) {
        if (0<=i2m && i2p<n2 && i3p<n3) {
          i1 = n1m;
          i1m = i1-1;
          x[i3][i2][i1] = d000[i3][i2][i1]*x[i3 ][i2 ][i1 ] -
                          l0pm[i3][i2][i1]*x[i3 ][i2p][i1m] -
                          l0p0[i3][i2][i1]*x[i3 ][i2p][i1 ] -
                          lpmm[i3][i2][i1]*x[i3p][i2m][i1m] -
                          lpm0[i3][i2][i1]*x[i3p][i2m][i1 ] -
                          lp0m[i3][i2][i1]*x[i3p][i2 ][i1m] -
                          lp00[i3][i2][i1]*x[i3p][i2 ][i1 ] -
                          lppm[i3][i2][i1]*x[i3p][i2p][i1m] -
                          lpp0[i3][i2][i1]*x[i3p][i2p][i1 ];
          for (i1=n1m-1,i1m=i1-1,i1p=i1+1; i1>=1; --i1,--i1m,--i1p) {
            x[i3][i2][i1] = d000[i3][i2][i1]*x[i3 ][i2 ][i1 ] -
                            l00p[i3][i2][i1]*x[i3 ][i2 ][i1p] -
                            l0pm[i3][i2][i1]*x[i3 ][i2p][i1m] -
                            l0p0[i3][i2][i1]*x[i3 ][i2p][i1 ] -
                            l0pp[i3][i2][i1]*x[i3 ][i2p][i1p] -
                            lpmm[i3][i2][i1]*x[i3p][i2m][i1m] -
                            lpm0[i3][i2][i1]*x[i3p][i2m][i1 ] -
                            lpmp[i3][i2][i1]*x[i3p][i2m][i1p] -
                            lp0m[i3][i2][i1]*x[i3p][i2 ][i1m] -
                            lp00[i3][i2][i1]*x[i3p][i2 ][i1 ] -
                            lp0p[i3][i2][i1]*x[i3p][i2 ][i1p] -
                            lppm[i3][i2][i1]*x[i3p][i2p][i1m] -
                            lpp0[i3][i2][i1]*x[i3p][i2p][i1 ] -
                            lppp[i3][i2][i1]*x[i3p][i2p][i1p];
          }
          x[i3][i2][i1] = d000[i3][i2][i1]*x[i3 ][i2 ][i1 ] -
                          l00p[i3][i2][i1]*x[i3 ][i2 ][i1p] -
                          l0p0[i3][i2][i1]*x[i3 ][i2p][i1 ] -
                          l0pp[i3][i2][i1]*x[i3 ][i2p][i1p] -
                          lpm0[i3][i2][i1]*x[i3p][i2m][i1 ] -
                          lpmp[i3][i2][i1]*x[i3p][i2m][i1p] -
                          lp00[i3][i2][i1]*x[i3p][i2 ][i1 ] -
                          lp0p[i3][i2][i1]*x[i3p][i2 ][i1p] -
                          lpp0[i3][i2][i1]*x[i3p][i2p][i1 ] -
                          lppp[i3][i2][i1]*x[i3p][i2p][i1p];
        } else {
          for (i1=n1m,i1m=i1-1,i1p=i1+1; i1>=0; --i1,--i1m,--i1p) {
                            xi  = d000[i3][i2][i1]*x[i3 ][i2 ][i1 ];
            if (i1p<n1)     xi -= l00p[i3][i2][i1]*x[i3 ][i2 ][i1p];
            if (i2p<n2) {
              if (0<=i1m)   xi -= l0pm[i3][i2][i1]*x[i3 ][i2p][i1m];
                            xi -= l0p0[i3][i2][i1]*x[i3 ][i2p][i1 ];
              if (i1p<n1)   xi -= l0pp[i3][i2][i1]*x[i3 ][i2p][i1p];
            }
            if (i3p<n3) {
              if (0<=i2m) {
                if (0<=i1m) xi -= lpmm[i3][i2][i1]*x[i3p][i2m][i1m];
                            xi -= lpm0[i3][i2][i1]*x[i3p][i2m][i1 ];
                if (i1p<n1) xi -= lpmp[i3][i2][i1]*x[i3p][i2m][i1p];
              }
              if (0<=i1m)   xi -= lp0m[i3][i2][i1]*x[i3p][i2 ][i1m];
                            xi -= lp00[i3][i2][i1]*x[i3p][i2 ][i1 ];
              if (i1p<n1)   xi -= lp0p[i3][i2][i1]*x[i3p][i2 ][i1p];
              if (i2p<n2) {
                if (0<=i1m) xi -= lppm[i3][i2][i1]*x[i3p][i2p][i1m];
                            xi -= lpp0[i3][i2][i1]*x[i3p][i2p][i1 ];
                if (i1p<n1) xi -= lppp[i3][i2][i1]*x[i3p][i2p][i1p];
              }
            }
            x[i3][i2][i1] = xi;
          }
        }
      }
    }
  }
  
  /**
   * Computes y = L*D*L'*x. For testing, only.
   */
  private void applyFactors(float[][][] x, float[][][] y) {
    ensureFactors();
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float[][][] d000 = _l[0];
    float[][][] l00p = _l[1];
    float[][][] l0pm = _l[2];
    float[][][] l0p0 = _l[3];
    float[][][] l0pp = _l[4];
    float[][][] lpmm = _l[5];
    float[][][] lpm0 = _l[6];
    float[][][] lpmp = _l[7];
    float[][][] lp0m = _l[8];
    float[][][] lp00 = _l[9];
    float[][][] lp0p = _l[10];
    float[][][] lppm = _l[11];
    float[][][] lpp0 = _l[12];
    float[][][] lppp = _l[13];

    // y = L'*x
    for (int i3=0,i3p=i3+1; i3<n3; ++i3,++i3p) {
      for (int i2=0,i2m=i2-1,i2p=i2+1; i2<n2; ++i2,++i2m,++i2p) {
        for (int i1=0,i1m=i1-1,i1p=i1+1; i1<n1; ++i1,++i1m,++i1p) {
          float yi = x[i3][i2][i1];
          if (i1p<n1)     yi += x[i3 ][i2 ][i1p]*l00p[i3][i2][i1];
          if (i2p<n2) {
            if (0<=i1m)   yi += x[i3 ][i2p][i1m]*l0pm[i3][i2][i1];
                          yi += x[i3 ][i2p][i1 ]*l0p0[i3][i2][i1];
            if (i1p<n1)   yi += x[i3 ][i2p][i1p]*l0pp[i3][i2][i1];
          }
          if (i3p<n3) {
            if (0<=i2m) {
              if (0<=i1m) yi += x[i3p][i2m][i1m]*lpmm[i3][i2][i1];
                          yi += x[i3p][i2m][i1 ]*lpm0[i3][i2][i1];
              if (i1p<n1) yi += x[i3p][i2m][i1p]*lpmp[i3][i2][i1];
            }
            if (0<=i1m)   yi += x[i3p][i2 ][i1m]*lp0m[i3][i2][i1];
                          yi += x[i3p][i2 ][i1 ]*lp00[i3][i2][i1];
            if (i1p<n1)   yi += x[i3p][i2 ][i1p]*lp0p[i3][i2][i1];
            if (i2p<n2) {
              if (0<=i1m) yi += x[i3p][i2p][i1m]*lppm[i3][i2][i1];
                          yi += x[i3p][i2p][i1 ]*lpp0[i3][i2][i1];
              if (i1p<n1) yi += x[i3p][i2p][i1p]*lppp[i3][i2][i1];
            }
          }
          y[i3][i2][i1] = yi;
        }
      }
    }

    // y = L*D*y
    for (int i3=n3m,i3p=i3+1; i3>=0; --i3,--i3p) {
      for (int i2=n2m,i2m=i2-1,i2p=i2+1; i2>=0; --i2,--i2m,--i2p) {
        for (int i1=n1m,i1m=i1-1,i1p=i1+1; i1>=0; --i1,--i1m,--i1p) {
          y[i3][i2][i1] /= d000[i3][i2][i1];
          float yi = y[i3][i2][i1];
          if (i1p<n1)     y[i3 ][i2 ][i1p] += l00p[i3][i2][i1]*yi;
          if (i2p<n2) {
            if (0<=i1m)   y[i3 ][i2p][i1m] += l0pm[i3][i2][i1]*yi;
                          y[i3 ][i2p][i1 ] += l0p0[i3][i2][i1]*yi;
            if (i1p<n1)   y[i3 ][i2p][i1p] += l0pp[i3][i2][i1]*yi;
          }
          if (i3p<n3) {
            if (0<=i2m) {
              if (0<=i1m) y[i3p][i2m][i1m] += lpmm[i3][i2][i1]*yi;
                          y[i3p][i2m][i1 ] += lpm0[i3][i2][i1]*yi;
              if (i1p<n1) y[i3p][i2m][i1p] += lpmp[i3][i2][i1]*yi;
            }
            if (0<=i1m)   y[i3p][i2 ][i1m] += lp0m[i3][i2][i1]*yi;
                          y[i3p][i2 ][i1 ] += lp00[i3][i2][i1]*yi;
            if (i1p<n1)   y[i3p][i2 ][i1p] += lp0p[i3][i2][i1]*yi;
            if (i2p<n2) {
              if (0<=i1m) y[i3p][i2p][i1m] += lppm[i3][i2][i1]*yi;
                          y[i3p][i2p][i1 ] += lpp0[i3][i2][i1]*yi;
              if (i1p<n1) y[i3p][i2p][i1p] += lppp[i3][i2][i1]*yi;
            }
          }
        }
      }
    }
  }

  /**
   * Factors this filter with incomplete Cholesky decomposition IC(0).
   * The approximate factorization is A ~ L*D*L', where L is a 
   * unit-lower-triangular matrix, and D is a diagonal matrix.
   * @param bias the amount by which to perturb the diagonal of A.
   * @return coefficients {d000,l00p,l0pm,l0p0,l0pp,lpm0,lp0m,lp00,lp0p,lpp0}.
   *  The elements in the array d000 are the inverse of the elements of the 
   *  diagonal matrix D.
   */
  private static float[][][][] factorIC0(float[][][][] a, float bias) {
    float[][][][] l = null;
    float bmin = (bias>0.0f)?bias:0.001f;
    for (float b=bias; l==null && b<1000.0f; b=max(bmin,2*b)) {
      l = attemptIC0(a,b);
      if (l==null)
        trace("factorIC0: failed for bias="+b);
      else
        trace("factorIC0: success for bias="+b);
    }
    return l;
  }
  private static float[][][][] attemptIC0(float[][][][] a, float bias) {
    float[][][] l000 = Array.copy(a[0]);
    float[][][] l00p = Array.copy(a[1]);
    float[][][] l0pm = Array.copy(a[2]);
    float[][][] l0p0 = Array.copy(a[3]);
    float[][][] l0pp = Array.copy(a[4]);
    float[][][] lpmm = Array.copy(a[5]);
    float[][][] lpm0 = Array.copy(a[6]);
    float[][][] lpmp = Array.copy(a[7]);
    float[][][] lp0m = Array.copy(a[8]);
    float[][][] lp00 = Array.copy(a[9]);
    float[][][] lp0p = Array.copy(a[10]);
    float[][][] lppm = Array.copy(a[11]);
    float[][][] lpp0 = Array.copy(a[12]);
    float[][][] lppp = Array.copy(a[13]);
    float[][][] d000 = l000; // will contain inverse of diagonal matrix D
    float[][][][] l = {l000,l00p,l0pm,l0p0,l0pp,
                       lpmm,lpm0,lpmp,lp0m,lp00,lp0p,lppm,lpp0,lppp};
    if (bias>0.0f)
      Array.mul(1.0f+bias,l000,l000);

    // Incomplete Cholesky decomposition, in-place.
    int n1 = a[0][0][0].length;
    int n2 = a[0][0].length;
    int n3 = a[0].length;
    int i1,i2,i3,i1m,i2m,i3m,i1p,i2p;
    for (i3=0,i3m=i3-1; i3<n3; ++i3,++i3m) {
      for (i2=0,i2m=i2-1,i2p=i2+1; i2<n2; ++i2,++i2m,++i2p) {
        if (0<=i2m && i2p<n2 && 0<=i3m) {
          i1 = 0;
          i1p = i1+1;
          l000[i3][i2][i1] -= 
            d000[i3 ][i2m][i1p]*l0pm[i3 ][i2m][i1p]*l0pm[i3 ][i2m][i1p] +
            d000[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ] +
            d000[i3m][i2p][i1p]*lpmm[i3m][i2p][i1p]*lpmm[i3m][i2p][i1p] +
            d000[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1p]*lp0m[i3m][i2 ][i1p]*lp0m[i3m][i2 ][i1p] +
            d000[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ] +
            d000[i3m][i2m][i1p]*lppm[i3m][i2m][i1p]*lppm[i3m][i2m][i1p] +
            d000[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ];
          l00p[i3][i2][i1] -= 
            d000[i3 ][i2m][i1p]*l0p0[i3 ][i2m][i1p]*l0p0[i3 ][i2m][i1p] +
            d000[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ] +
            d000[i3m][i2p][i1p]*lpm0[i3m][i2p][i1p]*lpm0[i3m][i2p][i1p] +
            d000[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1p]*lp00[i3m][i2 ][i1p]*lp00[i3m][i2 ][i1p] +
            d000[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ] +
            d000[i3m][i2m][i1p]*lpp0[i3m][i2m][i1p]*lpp0[i3m][i2m][i1p] +
            d000[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ];
          l0pm[i3][i2][i1] -= 
            d000[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ];
          l0p0[i3][i2][i1] -= 
            d000[i3m][i2p][i1p]*lp0m[i3m][i2p][i1p]*lp0m[i3m][i2p][i1p] +
            d000[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1p]*lppm[i3m][i2 ][i1p]*lppm[i3m][i2 ][i1p] +
            d000[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ];
          l0pp[i3][i2][i1] -= 
            d000[i3m][i2p][i1p]*lp00[i3m][i2p][i1p]*lp00[i3m][i2p][i1p] +
            d000[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1p]*lpp0[i3m][i2 ][i1p]*lpp0[i3m][i2 ][i1p] +
            d000[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ];
          lpmm[i3][i2][i1] -= 
            d000[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ];
          lpm0[i3][i2][i1] -= 
            d000[i3 ][i2m][i1p]*lp0m[i3 ][i2m][i1p]*lp0m[i3 ][i2m][i1p] +
            d000[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ];
          lpmp[i3][i2][i1] -= 
            d000[i3 ][i2m][i1p]*lp00[i3 ][i2m][i1p]*lp00[i3 ][i2m][i1p] +
            d000[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ];
          lp0m[i3][i2][i1] -= 
            d000[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ];
          lp00[i3][i2][i1] -= 
            d000[i3 ][i2m][i1p]*lppm[i3 ][i2m][i1p]*lppm[i3 ][i2m][i1p] +
            d000[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ];
          lp0p[i3][i2][i1] -= 
            d000[i3 ][i2m][i1p]*lpp0[i3 ][i2m][i1p]*lpp0[i3 ][i2m][i1p] +
            d000[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ];
          if (l000[i3][i2][i1]<=0.0f) 
            return null;
          d000[i3][i2][i1] = 1.0f/l000[i3][i2][i1];
          for (i1=1,i1m=i1-1,i1p=i1+1; i1<n1-1; ++i1,++i1m,++i1p) {
            //       l0pm  |  lpmm lp0m lppm
            //  l000 l0p0  |  lpm0 lp00 lpp0
            //  l00p l0pp  |  lpmp lp0p lppp // 57 updates
            l000[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*l00p[i3 ][i2 ][i1m]*l00p[i3 ][i2 ][i1m] +
              d000[i3 ][i2m][i1p]*l0pm[i3 ][i2m][i1p]*l0pm[i3 ][i2m][i1p] +
              d000[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ] +
              d000[i3 ][i2m][i1m]*l0pp[i3 ][i2m][i1m]*l0pp[i3 ][i2m][i1m] +
              d000[i3m][i2p][i1p]*lpmm[i3m][i2p][i1p]*lpmm[i3m][i2p][i1p] +
              d000[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ] +
              d000[i3m][i2p][i1m]*lpmp[i3m][i2p][i1m]*lpmp[i3m][i2p][i1m] +
              d000[i3m][i2 ][i1p]*lp0m[i3m][i2 ][i1p]*lp0m[i3m][i2 ][i1p] +
              d000[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ] +
              d000[i3m][i2 ][i1m]*lp0p[i3m][i2 ][i1m]*lp0p[i3m][i2 ][i1m] +
              d000[i3m][i2m][i1p]*lppm[i3m][i2m][i1p]*lppm[i3m][i2m][i1p] +
              d000[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ] +
              d000[i3m][i2m][i1m]*lppp[i3m][i2m][i1m]*lppp[i3m][i2m][i1m];
            l00p[i3][i2][i1] -= 
              d000[i3 ][i2m][i1p]*l0p0[i3 ][i2m][i1p]*l0p0[i3 ][i2m][i1p] +
              d000[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ] +
              d000[i3m][i2p][i1p]*lpm0[i3m][i2p][i1p]*lpm0[i3m][i2p][i1p] +
              d000[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ] +
              d000[i3m][i2 ][i1p]*lp00[i3m][i2 ][i1p]*lp00[i3m][i2 ][i1p] +
              d000[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ] +
              d000[i3m][i2m][i1p]*lpp0[i3m][i2m][i1p]*lpp0[i3m][i2m][i1p] +
              d000[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ];
            l0pm[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*l0p0[i3 ][i2 ][i1m]*l0p0[i3 ][i2 ][i1m] +
              d000[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ] +
              d000[i3m][i2p][i1m]*lp00[i3m][i2p][i1m]*lp00[i3m][i2p][i1m] +
              d000[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ] +
              d000[i3m][i2 ][i1m]*lpp0[i3m][i2 ][i1m]*lpp0[i3m][i2 ][i1m];
            l0p0[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*l0pp[i3 ][i2 ][i1m]*l0pp[i3 ][i2 ][i1m] +
              d000[i3m][i2p][i1p]*lp0m[i3m][i2p][i1p]*lp0m[i3m][i2p][i1p] +
              d000[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ] +
              d000[i3m][i2p][i1m]*lp0p[i3m][i2p][i1m]*lp0p[i3m][i2p][i1m] +
              d000[i3m][i2 ][i1p]*lppm[i3m][i2 ][i1p]*lppm[i3m][i2 ][i1p] +
              d000[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ] +
              d000[i3m][i2 ][i1m]*lppp[i3m][i2 ][i1m]*lppp[i3m][i2 ][i1m];
            l0pp[i3][i2][i1] -= 
              d000[i3m][i2p][i1p]*lp00[i3m][i2p][i1p]*lp00[i3m][i2p][i1p] +
              d000[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ] +
              d000[i3m][i2 ][i1p]*lpp0[i3m][i2 ][i1p]*lpp0[i3m][i2 ][i1p] +
              d000[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ];
            lpmm[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*lpm0[i3 ][i2 ][i1m]*lpm0[i3 ][i2 ][i1m] +
              d000[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ] +
              d000[i3 ][i2m][i1m]*lp00[i3 ][i2m][i1m]*lp00[i3 ][i2m][i1m];
            lpm0[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*lpmp[i3 ][i2 ][i1m]*lpmp[i3 ][i2 ][i1m] +
              d000[i3 ][i2m][i1p]*lp0m[i3 ][i2m][i1p]*lp0m[i3 ][i2m][i1p] +
              d000[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ] +
              d000[i3 ][i2m][i1m]*lp0p[i3 ][i2m][i1m]*lp0p[i3 ][i2m][i1m];
            lpmp[i3][i2][i1] -= 
              d000[i3 ][i2m][i1p]*lp00[i3 ][i2m][i1p]*lp00[i3 ][i2m][i1p] +
              d000[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ];
            lp0m[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*lp00[i3 ][i2 ][i1m]*lp00[i3 ][i2 ][i1m] +
              d000[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ] +
              d000[i3 ][i2m][i1m]*lpp0[i3 ][i2m][i1m]*lpp0[i3 ][i2m][i1m];
            lp00[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*lp0p[i3 ][i2 ][i1m]*lp0p[i3 ][i2 ][i1m] +
              d000[i3 ][i2m][i1p]*lppm[i3 ][i2m][i1p]*lppm[i3 ][i2m][i1p] +
              d000[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ] +
              d000[i3 ][i2m][i1m]*lppp[i3 ][i2m][i1m]*lppp[i3 ][i2m][i1m];
            lp0p[i3][i2][i1] -= 
              d000[i3 ][i2m][i1p]*lpp0[i3 ][i2m][i1p]*lpp0[i3 ][i2m][i1p] +
              d000[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ];
            lppm[i3][i2][i1] -= 
              d000[i3 ][i2 ][i1m]*lpp0[i3 ][i2 ][i1m]*lpp0[i3 ][i2 ][i1m];
            lpp0[i3][i2][i1] -= 
              d000[i3 ][i2m][i1m]*lppp[i3 ][i2 ][i1m]*lppp[i3 ][i2 ][i1m];
            if (l000[i3][i2][i1]<=0.0f) 
              return null;
            d000[i3][i2][i1] = 1.0f/l000[i3][i2][i1];
          }
          l000[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*l00p[i3 ][i2 ][i1m]*l00p[i3 ][i2 ][i1m] +
            d000[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ] +
            d000[i3 ][i2m][i1m]*l0pp[i3 ][i2m][i1m]*l0pp[i3 ][i2m][i1m] +
            d000[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ] +
            d000[i3m][i2p][i1m]*lpmp[i3m][i2p][i1m]*lpmp[i3m][i2p][i1m] +
            d000[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ] +
            d000[i3m][i2 ][i1m]*lp0p[i3m][i2 ][i1m]*lp0p[i3m][i2 ][i1m] +
            d000[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ] +
            d000[i3m][i2m][i1m]*lppp[i3m][i2m][i1m]*lppp[i3m][i2m][i1m];
          l00p[i3][i2][i1] -= 
            d000[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ] +
            d000[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ] +
            d000[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ];
          l0pm[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*l0p0[i3 ][i2 ][i1m]*l0p0[i3 ][i2 ][i1m] +
            d000[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ] +
            d000[i3m][i2p][i1m]*lp00[i3m][i2p][i1m]*lp00[i3m][i2p][i1m] +
            d000[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ] +
            d000[i3m][i2 ][i1m]*lpp0[i3m][i2 ][i1m]*lpp0[i3m][i2 ][i1m];
          l0p0[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*l0pp[i3 ][i2 ][i1m]*l0pp[i3 ][i2 ][i1m] +
            d000[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ] +
            d000[i3m][i2p][i1m]*lp0p[i3m][i2p][i1m]*lp0p[i3m][i2p][i1m] +
            d000[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ] +
            d000[i3m][i2 ][i1m]*lppp[i3m][i2 ][i1m]*lppp[i3m][i2 ][i1m];
          l0pp[i3][i2][i1] -= 
            d000[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ] +
            d000[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ];
          lpmm[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*lpm0[i3 ][i2 ][i1m]*lpm0[i3 ][i2 ][i1m] +
            d000[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ] +
            d000[i3 ][i2m][i1m]*lp00[i3 ][i2m][i1m]*lp00[i3 ][i2m][i1m];
          lpm0[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*lpmp[i3 ][i2 ][i1m]*lpmp[i3 ][i2 ][i1m] +
            d000[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ] +
            d000[i3 ][i2m][i1m]*lp0p[i3 ][i2m][i1m]*lp0p[i3 ][i2m][i1m];
          lpmp[i3][i2][i1] -= 
            d000[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ];
          lp0m[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*lp00[i3 ][i2 ][i1m]*lp00[i3 ][i2 ][i1m] +
            d000[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ] +
            d000[i3 ][i2m][i1m]*lpp0[i3 ][i2m][i1m]*lpp0[i3 ][i2m][i1m];
          lp00[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*lp0p[i3 ][i2 ][i1m]*lp0p[i3 ][i2 ][i1m] +
            d000[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ] +
            d000[i3 ][i2m][i1m]*lppp[i3 ][i2m][i1m]*lppp[i3 ][i2m][i1m];
          lp0p[i3][i2][i1] -= 
            d000[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ];
          lppm[i3][i2][i1] -= 
            d000[i3 ][i2 ][i1m]*lpp0[i3 ][i2 ][i1m]*lpp0[i3 ][i2 ][i1m];
          lpp0[i3][i2][i1] -= 
            d000[i3 ][i2m][i1m]*lppp[i3 ][i2 ][i1m]*lppp[i3 ][i2 ][i1m];
          if (l000[i3][i2][i1]<=0.0f) 
            return null;
          d000[i3][i2][i1] = 1.0f/l000[i3][i2][i1];
        } else {
          for (i1=0,i1m=i1-1,i1p=i1+1; i1<n1; ++i1,++i1m,++i1p) {
            //       l0pm  |  lpmm lp0m lppm
            //  l000 l0p0  |  lpm0 lp00 lpp0
            //  l00p l0pp  |  lpmp lp0p lppp // 57 updates
            if (0<=i1m) {
              l000[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*l00p[i3 ][i2 ][i1m]*l00p[i3 ][i2 ][i1m];
              l0pm[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*l0p0[i3 ][i2 ][i1m]*l0p0[i3 ][i2 ][i1m];
              l0p0[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*l0pp[i3 ][i2 ][i1m]*l0pp[i3 ][i2 ][i1m];
              lpmm[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*lpm0[i3 ][i2 ][i1m]*lpm0[i3 ][i2 ][i1m];
              lpm0[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*lpmp[i3 ][i2 ][i1m]*lpmp[i3 ][i2 ][i1m];
              lp0m[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*lp00[i3 ][i2 ][i1m]*lp00[i3 ][i2 ][i1m];
              lp00[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*lp0p[i3 ][i2 ][i1m]*lp0p[i3 ][i2 ][i1m];
              lppm[i3][i2][i1] -= 
                d000[i3 ][i2 ][i1m]*lpp0[i3 ][i2 ][i1m]*lpp0[i3 ][i2 ][i1m];
              lpp0[i3][i2][i1] -= 
                d000[i3 ][i2m][i1m]*lppp[i3 ][i2 ][i1m]*lppp[i3 ][i2 ][i1m];
            }
            if (0<=i2m) {
              if (0<=i1m) {
                l000[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1m]*l0pp[i3 ][i2m][i1m]*l0pp[i3 ][i2m][i1m];
                lpmm[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1m]*lp00[i3 ][i2m][i1m]*lp00[i3 ][i2m][i1m];
                lpm0[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1m]*lp0p[i3 ][i2m][i1m]*lp0p[i3 ][i2m][i1m];
                lp0m[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1m]*lpp0[i3 ][i2m][i1m]*lpp0[i3 ][i2m][i1m];
                lp00[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1m]*lppp[i3 ][i2m][i1m]*lppp[i3 ][i2m][i1m];
              }
              l000[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ]*l0p0[i3 ][i2m][i1 ];
              l00p[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ]*l0pp[i3 ][i2m][i1 ];
              lpmm[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ]*lp0m[i3 ][i2m][i1 ];
              lpm0[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ]*lp00[i3 ][i2m][i1 ];
              lpmp[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ]*lp0p[i3 ][i2m][i1 ];
              lp0m[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ]*lppm[i3 ][i2m][i1 ];
              lp00[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ]*lpp0[i3 ][i2m][i1 ];
              lp0p[i3][i2][i1] -= 
                d000[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ]*lppp[i3 ][i2m][i1 ];
              if (i1p<n1) {
                l000[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1p]*l0pm[i3 ][i2m][i1p]*l0pm[i3 ][i2m][i1p];
                l00p[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1p]*l0p0[i3 ][i2m][i1p]*l0p0[i3 ][i2m][i1p];
                lpm0[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1p]*lp0m[i3 ][i2m][i1p]*lp0m[i3 ][i2m][i1p];
                lpmp[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1p]*lp00[i3 ][i2m][i1p]*lp00[i3 ][i2m][i1p];
                lp00[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1p]*lppm[i3 ][i2m][i1p]*lppm[i3 ][i2m][i1p];
                lp0p[i3][i2][i1] -= 
                  d000[i3 ][i2m][i1p]*lpp0[i3 ][i2m][i1p]*lpp0[i3 ][i2m][i1p];
              }
            }
            if (0<=i3m) {
              if (0<=i2m && 0<=i1m) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2m][i1m]*lppp[i3m][i2m][i1m]*lppp[i3m][i2m][i1m];
              }
              if (0<=i2m) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ]*lpp0[i3m][i2m][i1 ];
                l00p[i3][i2][i1] -= 
                  d000[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ]*lppp[i3m][i2m][i1 ];
              }
              if (0<=i2m && i1p<n1) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2m][i1p]*lppm[i3m][i2m][i1p]*lppm[i3m][i2m][i1p];
                l00p[i3][i2][i1] -= 
                  d000[i3m][i2m][i1p]*lpp0[i3m][i2m][i1p]*lpp0[i3m][i2m][i1p];
              }
              if (0<=i1m) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1m]*lp0p[i3m][i2 ][i1m]*lp0p[i3m][i2 ][i1m];
                l0pm[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1m]*lpp0[i3m][i2 ][i1m]*lpp0[i3m][i2 ][i1m];
                l0p0[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1m]*lppp[i3m][i2 ][i1m]*lppp[i3m][i2 ][i1m];
              }
              l000[i3][i2][i1] -= 
                d000[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ]*lp00[i3m][i2 ][i1 ];
              l00p[i3][i2][i1] -= 
                d000[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ]*lp0p[i3m][i2 ][i1 ];
              l0pm[i3][i2][i1] -= 
                d000[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ]*lppm[i3m][i2 ][i1 ];
              l0p0[i3][i2][i1] -= 
                d000[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ]*lpp0[i3m][i2 ][i1 ];
              l0pp[i3][i2][i1] -= 
                d000[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ]*lppp[i3m][i2 ][i1 ];
              if (i1p<n1) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1p]*lp0m[i3m][i2 ][i1p]*lp0m[i3m][i2 ][i1p];
                l00p[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1p]*lp00[i3m][i2 ][i1p]*lp00[i3m][i2 ][i1p];
                l0p0[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1p]*lppm[i3m][i2 ][i1p]*lppm[i3m][i2 ][i1p];
                l0pp[i3][i2][i1] -= 
                  d000[i3m][i2 ][i1p]*lpp0[i3m][i2 ][i1p]*lpp0[i3m][i2 ][i1p];
              }
              if (i2p<n2 && 0<=i1m) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2p][i1m]*lpmp[i3m][i2p][i1m]*lpmp[i3m][i2p][i1m];
                l0pm[i3][i2][i1] -= 
                  d000[i3m][i2p][i1m]*lp00[i3m][i2p][i1m]*lp00[i3m][i2p][i1m];
                l0p0[i3][i2][i1] -= 
                  d000[i3m][i2p][i1m]*lp0p[i3m][i2p][i1m]*lp0p[i3m][i2p][i1m];
              }
              if (i2p<n2) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ]*lpm0[i3m][i2p][i1 ];
                l00p[i3][i2][i1] -= 
                  d000[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ]*lpmp[i3m][i2p][i1 ];
                l0pm[i3][i2][i1] -= 
                  d000[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ]*lp0m[i3m][i2p][i1 ];
                l0p0[i3][i2][i1] -= 
                  d000[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ]*lp00[i3m][i2p][i1 ];
                l0pp[i3][i2][i1] -= 
                  d000[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ]*lp0p[i3m][i2p][i1 ];
              }
              if (i2p<n2 && i1p<n1) {
                l000[i3][i2][i1] -= 
                  d000[i3m][i2p][i1p]*lpmm[i3m][i2p][i1p]*lpmm[i3m][i2p][i1p];
                l00p[i3][i2][i1] -= 
                  d000[i3m][i2p][i1p]*lpm0[i3m][i2p][i1p]*lpm0[i3m][i2p][i1p];
                l0p0[i3][i2][i1] -= 
                  d000[i3m][i2p][i1p]*lp0m[i3m][i2p][i1p]*lp0m[i3m][i2p][i1p];
                l0pp[i3][i2][i1] -= 
                  d000[i3m][i2p][i1p]*lp00[i3m][i2p][i1p]*lp00[i3m][i2p][i1p];
              }
            }
            if (l000[i3][i2][i1]<=0.0f) 
              return null;
            d000[i3][i2][i1] = 1.0f/l000[i3][i2][i1];
          }
        }
      }
    }

    // At this point, L is lower-triangular. Now make L have unit diagonal. 
    // This will enable application of inverse of L*D*L' without division.
    for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
        for (i1=0; i1<n1; ++i1) {
          float d000i = d000[i3][i2][i1];
          l00p[i3][i2][i1] *= d000i;
          l0pm[i3][i2][i1] *= d000i;
          l0p0[i3][i2][i1] *= d000i;
          l0pp[i3][i2][i1] *= d000i;
          lpmm[i3][i2][i1] *= d000i;
          lpm0[i3][i2][i1] *= d000i;
          lpmp[i3][i2][i1] *= d000i;
          lp0m[i3][i2][i1] *= d000i;
          lp00[i3][i2][i1] *= d000i;
          lp0p[i3][i2][i1] *= d000i;
          lppm[i3][i2][i1] *= d000i;
          lpp0[i3][i2][i1] *= d000i;
          lppp[i3][i2][i1] *= d000i;

        }
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
   //{s000,s00p,s0pm,s0p0,s0pp,spm0,sp0m,sp00,sp0p,spp0}.
  private static void testFactor() {
    java.util.Random r = new java.util.Random();
    int n1 = 5;
    int n2 = 4;
    int n3 = 3;
    float[][][][] s = new float[14][n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float s0 = 0.0f;
          for (int i=1; i<14; ++i) {
            s[i][i3][i2][i1] = -1.0f; //-r.nextFloat();
          }
          s[0][i3][i2][i1] = 27.0f;
        }
      }
    }
    LocalSpd27Filter lsf = new LocalSpd27Filter(s);
    //float[][] a = lsf.getMatrix();
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
    float[][][] x = Array.randfloat(n1,n2,n3);
    float[][][] y = Array.randfloat(n1,n2,n3);
    float[][][] z = Array.randfloat(n1,n2,n3);
    float[][][] w = Array.randfloat(n1,n2,n3);
    lsf.apply(x,y);
    lsf.applyApproximate(x,z);
    lsf.applyApproximateInverse(z,w);
    Array.dump(x);
    Array.dump(y);
    Array.dump(z);
    Array.dump(Array.sub(z,y));
    Array.dump(w);
    Array.dump(Array.sub(w,x));
  }

  public static void main(String[] args) {
    testFactor();
  }
}
