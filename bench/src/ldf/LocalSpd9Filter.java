/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local symmetric positive-definite (SPD) filter with a 2-D 9-point stencil.
 * This filter is local in the sense that filter coefficients may differ for 
 * each sample. These coefficients form a 9-point stencil:
 * <pre><code>
 * y[i2][i1] = smm*x[i2-1][i1-1] + s0m*x[i2  ][i1-1] + spm*x[i2+1][i1-1] +
 *             sm0*x[i2-1][i1  ] + s00*x[i2  ][i1  ] + sp0*x[i2+1][i1  ] +
 *             smp*x[i2-1][i1+1] + s0p*x[i2  ][i1+1] + spp*x[i2+1][i1+1]
 * </code></pre>
 * (The suffixes m, 0, and p denote minus, zero, and plus, respectively.)
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
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.11
 */
public class LocalSpd9Filter {

  /**
   * Constructs a local SPD filter with nine-point stencil.
   * @param s00 stencil coefficient for lag [ 0][ 0].
   * @param s0p stencil coefficient for lag [ 0][ 1].
   * @param spm stencil coefficient for lag [ 1][-1].
   * @param sp0 stencil coefficient for lag [ 1][ 0].
   * @param spp stencil coefficient for lag [ 1][ 1].
   */
  public LocalSpd9Filter(
    float[][] s00, float[][] s0p, 
    float[][] spm, float[][] sp0, float[][] spp) 
  {
    _s = new float[][][]{s00,s0p,spm,sp0,spp};
  }

  /**
   * Applies this filter by computing y = A*x.
   * @param x input array. Must be distinct from y.
   * @param y output array. Must be distinct from x.
   */
  public void apply(float[][] x, float[][] y) {
    if (_factored) {
      applyFactors(_s,x,y);
    } else {
      apply(_s,x,y);
    }
  }

  /**
   * Solves A*x = b, assuming that this filter has been factored.
   * @param b the input right-hand-side.
   * @param x the output app
   */
  public void solve(float[][] b, float[][] x) {
    Check.state(_factored,"filter has already been factored");
    solveWithFactors(_s,b,x);
  }

  /**
   * Factors this filter approximately via incomplete Cholesky decomposition.
   * This method uses zero fill-in.
   */
  public void factorIC0() {
    Check.state(!_factored,"filter has not already been factored");
    factorIC0(_s);
    _factored = true;
  }

  /**
   * Gets the sparse matrix A equivalent to this filter.
   * Most elements in this matrix will be zero. For small numbers of
   * samples, it may be useful for visualization of matrix sparsity.
   * <p>
   * If this filter is factored, the returned matrix represents the
   * decomposition A = L*D*L', and the diagonal contains the inverse
   * of the elements of the diagonal matrix D.
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
    int i1,i2,i;
    for (i2=0,i=0; i2<n2m; ++i2) {
      i1 = 0;
      a[i     ][i     ] = s00[i2][i1];
      a[i     ][i+1   ] = s0p[i2][i1];
      a[i     ][i+n1  ] = sp0[i2][i1];
      a[i     ][i+n1+1] = spp[i2][i1];
      a[i+1   ][i     ] = a[i][i+1];
      a[i+n1  ][i     ] = a[i][i+n1  ];
      a[i+n1+1][i     ] = a[i][i+n1+1];
      for (i1=1,++i; i1<n1m; ++i1,++i) {
        a[i     ][i     ] = s00[i2][i1];
        a[i     ][i+1   ] = s0p[i2][i1];
        a[i     ][i+n1-1] = spm[i2][i1];
        a[i     ][i+n1  ] = sp0[i2][i1];
        a[i     ][i+n1+1] = spp[i2][i1];
        a[i+1   ][i     ] = a[i][i+1];
        a[i+n1-1][i     ] = a[i][i+n1-1];
        a[i+n1  ][i     ] = a[i][i+n1  ];
        a[i+n1+1][i     ] = a[i][i+n1+1];
      }
      a[i     ][i     ] = s00[i2][i1];
      a[i     ][i+n1-1] = spm[i2][i1];
      a[i     ][i+n1  ] = sp0[i2][i1];
      a[i+n1-1][i     ] = a[i][i+n1-1];
      a[i+n1  ][i     ] = a[i][i+n1  ];
      ++i;
    }
    for (i1=0; i1<n1m; ++i1,++i) {
      a[i     ][i     ] = s00[i2][i1];
      a[i     ][i+1   ] = s0p[i2][i1];
      a[i+1   ][i     ] = a[i][i+1];
    }
    a[i     ][i     ] = s00[i2][i1];
    if (_factored) {
      for (i=0; i<n; ++i)
        a[i][i] = 1.0f/a[i][i];
    }
    return a;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _s;
  private boolean _factored;

  /**
   * Computes y = A*x, with non-zero elements of A stored in the array s.
   */
  private void apply(float[][][] s, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] s00 = s[0];
    float[][] s0p = s[1];
    float[][] spm = s[2];
    float[][] sp0 = s[3];
    float[][] spp = s[4];
    int i1,i2;
    i1 = n1m;
    i2 = n2m;
    y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
    for (i1=n1m-1; i1>=0; --i1) {
      y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += s0p[i2][i1]*x[i2  ][i1+1];
      y[i2  ][i1+1] += s0p[i2][i1]*x[i2  ][i1  ];
    }
    for (i2=n2m-1; i2>=0; --i2) {
      i1 = n1m;
      y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += spm[i2][i1]*x[i2+1][i1-1];
      y[i2+1][i1-1] += spm[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += sp0[i2][i1]*x[i2+1][i1  ];
      y[i2+1][i1  ] += sp0[i2][i1]*x[i2  ][i1  ];
      for (i1=n1m-1; i1>=1; --i1) {
        y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
        y[i2  ][i1  ] += s0p[i2][i1]*x[i2  ][i1+1];
        y[i2  ][i1+1] += s0p[i2][i1]*x[i2  ][i1  ];
        y[i2  ][i1  ] += spm[i2][i1]*x[i2+1][i1-1];
        y[i2+1][i1-1] += spm[i2][i1]*x[i2  ][i1  ];
        y[i2  ][i1  ] += sp0[i2][i1]*x[i2+1][i1  ];
        y[i2+1][i1  ] += sp0[i2][i1]*x[i2  ][i1  ];
        y[i2  ][i1  ] += spp[i2][i1]*x[i2+1][i1+1];
        y[i2+1][i1+1] += spp[i2][i1]*x[i2  ][i1  ];
      }
      y[i2  ][i1  ]  = s00[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += s0p[i2][i1]*x[i2  ][i1+1];
      y[i2  ][i1+1] += s0p[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += sp0[i2][i1]*x[i2+1][i1  ];
      y[i2+1][i1  ] += sp0[i2][i1]*x[i2  ][i1  ];
      y[i2  ][i1  ] += spp[i2][i1]*x[i2+1][i1+1];
      y[i2+1][i1+1] += spp[i2][i1]*x[i2  ][i1  ];
    }
  }

  /**
   * Factors this filter with incomplete Cholesky decomposition IC(0).
   * Elements of the diagonal matrix d00 (inverse of a00) are stored in a[0].
   */
  private static void factorIC0(float[][][] a) {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    float[][] d00 = a[0];
    int i1 = 0;
    int i2 = 0;
    d00[i2][i1] = 1.0f/a00[i2][i1];
    for (i1=1; i1<n1; ++i1) {
      a00[i2][i1] -= d00[i2  ][i1-1]*a0p[i2  ][i1-1]*a0p[i2  ][i1-1];
      apm[i2][i1] -= d00[i2  ][i1-1]*ap0[i2  ][i1-1]*a0p[i2  ][i1-1];
      ap0[i2][i1] -= d00[i2  ][i1-1]*app[i2  ][i1-1]*a0p[i2  ][i1-1];
      d00[i2][i1] = 1.0f/a00[i2][i1];
    }
    for (i2=1; i2<n2; ++i2) {
      i1 = 0;
      a00[i2][i1] -= d00[i2-1][i1+1]*apm[i2-1][i1+1]*apm[i2-1][i1+1] +
                     d00[i2-1][i1  ]*ap0[i2-1][i1  ]*ap0[i2-1][i1  ];
      a0p[i2][i1] -= d00[i2-1][i1  ]*app[i2-1][i1  ]*ap0[i2-1][i1  ];
      d00[i2][i1] = 1.0f/a00[i2][i1];
      for (i1=1; i1<n1-1; ++i1) {
        a00[i2][i1] -= d00[i2  ][i1-1]*a0p[i2  ][i1-1]*a0p[i2  ][i1-1] +
                       d00[i2-1][i1+1]*apm[i2-1][i1+1]*apm[i2-1][i1+1] +
                       d00[i2-1][i1  ]*ap0[i2-1][i1  ]*ap0[i2-1][i1  ] +
                       d00[i2-1][i1-1]*app[i2-1][i1-1]*app[i2-1][i1-1];
        a0p[i2][i1] -= d00[i2-1][i1+1]*ap0[i2-1][i1+1]*apm[i2-1][i1+1] +
                       d00[i2-1][i1  ]*app[i2-1][i1  ]*ap0[i2-1][i1  ];
        apm[i2][i1] -= d00[i2  ][i1-1]*ap0[i2  ][i1-1]*a0p[i2  ][i1-1];
        ap0[i2][i1] -= d00[i2  ][i1-1]*app[i2  ][i1-1]*a0p[i2  ][i1-1];
        d00[i2][i1] = 1.0f/a00[i2][i1];
      }
      a00[i2][i1] -= d00[i2  ][i1-1]*a0p[i2  ][i1-1]*a0p[i2  ][i1-1] +
                     d00[i2-1][i1  ]*ap0[i2-1][i1  ]*ap0[i2-1][i1  ] +
                     d00[i2-1][i1-1]*app[i2-1][i1-1]*app[i2-1][i1-1];
      a0p[i2][i1] -= d00[i2-1][i1  ]*app[i2-1][i1  ]*ap0[i2-1][i1  ];
      apm[i2][i1] -= d00[i2  ][i1-1]*ap0[i2  ][i1-1]*a0p[i2  ][i1-1];
      ap0[i2][i1] -= d00[i2  ][i1-1]*app[i2  ][i1-1]*a0p[i2  ][i1-1];
      d00[i2][i1] = 1.0f/a00[i2][i1];
    }
  }
  
  /**
   * Computes y = L*D*L'*x, with factors D and L' stored in the array s.
   */
  private static void applyFactors(float[][][] s, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] s00 = s[0];
    float[][] s0p = s[1];
    float[][] spm = s[2];
    float[][] sp0 = s[3];
    float[][] spp = s[4];
    int i1,i2;
    float y00;

    // y = L'*x
    for (i2=0; i2<n2m; ++i2) {
      i1 = 0;
      y00  = x[i2  ][i1  ]/s00[i2][i1];
      y00 += s0p[i2][i1]*x[i2  ][i1+1];
      y00 += sp0[i2][i1]*x[i2+1][i1  ];
      y00 += spp[i2][i1]*x[i2+1][i1+1];
      y[i2  ][i1  ]  = y00;
      for (i1=1; i1<n1m; ++i1) {
        y00  = x[i2  ][i1  ]/s00[i2][i1];
        y00 += s0p[i2][i1]*x[i2  ][i1+1];
        y00 += spm[i2][i1]*x[i2+1][i1-1];
        y00 += sp0[i2][i1]*x[i2+1][i1  ];
        y00 += spp[i2][i1]*x[i2+1][i1+1];
        y[i2  ][i1  ]  = y00;
      }
      y00  = x[i2  ][i1  ]/s00[i2][i1];
      y00 += spm[i2][i1]*x[i2+1][i1-1];
      y00 += sp0[i2][i1]*x[i2+1][i1  ];
      y[i2  ][i1  ]  = y00;
    }
    for (i1=0; i1<n1m; ++i1) {
      y00  = x[i2  ][i1  ]/s00[i2][i1];
      y00 += s0p[i2][i1]*x[i2  ][i1+1];
      y[i2  ][i1  ]  = y00;
    }
    y[i2  ][i1  ]  = x[i2  ][i1  ]/s00[i2][i1];

    // y = L*D*y
    i2 = n2m;
    i1 = n1m;
    for (i1=n1m-1; i1>=0; --i1) {
      y00 = s00[i2][i1]*y[i2  ][i1  ];
      y[i2  ][i1+1] += s0p[i2][i1]*y00;
    }
    for (i2=n2m-1; i2>=0; --i2) {
      i1 = n1m;
      y00 = s00[i2][i1]*y[i2  ][i1  ];
      y[i2+1][i1-1] += spm[i2][i1]*y00;
      y[i2+1][i1  ] += sp0[i2][i1]*y00;
      for (i1=n1m-1; i1>=1; --i1) {
        y00 = s00[i2][i1]*y[i2  ][i1  ];
        y[i2  ][i1+1] += s0p[i2][i1]*y00;
        y[i2+1][i1-1] += spm[i2][i1]*y00;
        y[i2+1][i1  ] += sp0[i2][i1]*y00;
        y[i2+1][i1+1] += spp[i2][i1]*y00;
      }
      y00 = s00[i2][i1]*y[i2  ][i1  ];
      y[i2  ][i1+1] += s0p[i2][i1]*y00;
      y[i2+1][i1  ] += sp0[i2][i1]*y00;
      y[i2+1][i1+1] += spp[i2][i1]*y00;
    }
  }

  /**
   * Solves L*D*L'*x = b, for elements of factors L, D, and L' stored in a.
   */
  private void solveWithFactors(float[][][] a, float[][] b, float[][] x) {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    float[][] d00 = a[0];
    int i1,i2;
    float yi,zi;

    // Solve L*z = b.
    float[][] z = x;
    Array.zero(z);
    for (i2=0; i2<n2m; ++i2) {
      i1 = 0;
      zi = d00[i2][i1]*(b[i2][i1]+z[i2][i1]);
      z[i2  ][i1  ] = zi;
      z[i2  ][i1+1] -= a0p[i2][i1]*zi;
      z[i2+1][i1  ] -= ap0[i2][i1]*zi;
      z[i2+1][i1+1] -= app[i2][i1]*zi;
      for (i1=1; i1<n1m; ++i1) {
        zi = d00[i2][i1]*(b[i2][i1]+z[i2][i1]);
        z[i2  ][i1  ] = zi;
        z[i2  ][i1+1] -= a0p[i2][i1]*zi;
        z[i2+1][i1-1] -= apm[i2][i1]*zi;
        z[i2+1][i1  ] -= ap0[i2][i1]*zi;
        z[i2+1][i1+1] -= app[i2][i1]*zi;
      }
      zi = d00[i2][i1]*(b[i2][i1]+z[i2][i1]);
      z[i2  ][i1  ] = zi;
      z[i2+1][i1-1] -= apm[i2][i1]*zi;
      z[i2+1][i1  ] -= ap0[i2][i1]*zi;
    }
    i1 = 0;
    zi = d00[i2][i1]*(b[i2][i1]+z[i2][i1]);
    z[i2  ][i1  ] = zi;
    z[i2  ][i1+1] -= a0p[i2][i1]*zi;
    for (i1=1; i1<n1m; ++i1) {
      zi = d00[i2][i1]*(b[i2][i1]+z[i2][i1]);
      z[i2  ][i1  ] = zi;
      z[i2  ][i1+1] -= a0p[i2][i1]*zi;
    }
    zi = d00[i2][i1]*b[i2][i1];
    z[i2  ][i1  ] = zi;

    // Solve D*y = z and L'*x = y.
    i2 = n2m;
    i1 = n1m;
    yi = z[i2][i1]/d00[i2][i1];
    x[i2][i1] = d00[i2][i1]*yi;
    for (i1=n1m-1; i1>=0; --i1) {
      yi = z[i2][i1]/d00[i2][i1];
      yi -= a0p[i2][i1]*x[i2  ][i1+1];
      x[i2][i1] = d00[i2][i1]*yi;
    }
    for (i2=n2m-1; i2>=0; --i2) {
      i1 = n1m;
      yi = z[i2][i1]/d00[i2][i1];
      yi -= ap0[i2][i1]*x[i2+1][i1  ];
      yi -= apm[i2][i1]*x[i2+1][i1-1];
      x[i2][i1] = d00[i2][i1]*yi;
      for (i1=n1m-1; i1>=1; --i1) {
        yi = z[i2][i1]/d00[i2][i1];
        yi -= app[i2][i1]*x[i2+1][i1+1];
        yi -= ap0[i2][i1]*x[i2+1][i1  ];
        yi -= apm[i2][i1]*x[i2+1][i1-1];
        yi -= a0p[i2][i1]*x[i2  ][i1+1];
        x[i2][i1] = d00[i2][i1]*yi;
      }
      yi = z[i2][i1]/d00[i2][i1];
      yi -= app[i2][i1]*x[i2+1][i1+1];
      yi -= ap0[i2][i1]*x[i2+1][i1  ];
      yi -= a0p[i2][i1]*x[i2  ][i1+1];
      x[i2][i1] = d00[i2][i1]*yi;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static float[][] mul(float[][] a, float[][] b) {
    int n = a.length;
    float[][] c = new float[n][n];
    for (int i=0; i<n; ++i) {
      for (int j=0; j<n; ++j) {
        for (int k=0; k<n; ++k) {
          c[i][j] += a[i][k]*b[k][i];
        }
      }
    }
    return c;
  }

  private static void testFactor() {
    int n1 = 5;
    int n2 = 5;
    float[][] x = Array.zerofloat(n1,n2);  x[2][2] = 1.0f;
    float[][] y = Array.randfloat(n1,n2);
    float[][] z = Array.randfloat(n1,n2);
    float[][] w = Array.randfloat(n1,n2);
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = 
      new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    float[][][] s = ldk.getStencils(ldt);
    Array.add(0.1f,s[0],s[0]);
    LocalSpd9Filter lsf = new LocalSpd9Filter(s[0],s[1],s[2],s[3],s[4]);
    lsf.apply(x,y);
    lsf.factorIC0();
    lsf.apply(x,z);
    lsf.solve(z,w);
    Array.dump(y);
    Array.dump(z);
    Array.dump(w);
  }

  private static void testMatrix() {
    int n1 = 4;
    int n2 = 4;
    //float[][] d0 = Array.randfloat(n1,n2);
    //float[][] d1 = Array.randfloat(n1,n2);
    //float[][] v1 = Array.randfloat(n1,n2);
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = 
      new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    float[][][] s = ldk.getStencils(ldt);
    LocalSpd9Filter lsf = new LocalSpd9Filter(s[0],s[1],s[2],s[3],s[4]);
    float[][] a = lsf.getMatrix();
    edu.mines.jtk.mosaic.SimplePlot.asPixels(a);
    lsf.factorIC0();
    float[][] b = lsf.getMatrix();
    edu.mines.jtk.mosaic.SimplePlot.asPixels(b);
  }

  public static void main(String[] args) {
    testFactor();
    //testMatrix();
  }
}
