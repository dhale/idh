/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * A multigrid solver for 2-D grids. Solves Ax = b, where A is a linear
 * operator, b is a known function sampled on a 2-D grid, and x is an
 * unknown sampled function to be computed.
 * <p>
 * This solver currently supports operators with 3x3 stencils. For such 
 * stencils, each sample of b is coupled to nine samples nearest to the 
 * corresponding sample of x. That is,
 * <pre><code>
 * b[i2][i1] = 
 *   a[0]*x[i2-1][i1-1] + a[1]*x[i2-1][i1  ] + a[2]*x[i2-1][i1+1] +
 *   a[3]*x[i2  ][i1-1] + a[4]*x[i2  ][i1  ] + a[5]*x[i2  ][i1+1] +
 *   a[6]*x[i2+1][i1-1] + a[7]*x[i2+1][i1  ] + a[8]*x[i2+1][i1+1];
 * </code></pre>
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.17
 */
public class Multigrid2 {

  /**
   * A linear operator represented by a 3x3 nine-point stencil. 
   */
  public interface A33 {

    /**
     * Gets nine stencil coefficients for specified indices.
     * @param i1 index in 1st dimension.
     * @param i2 index in 2nd dimension.
     * @param a array of stencil coefficients.
     */
    public void get(int i1, int i2, float[] a);
  }

  /**
   * A simple implementation for 3x3 stencils.
   */
  public static class SimpleA33 implements A33 {

    /**
     * Constructs a simple 3x3 stencil. See description above.
     * @param a array[n2][n1][9] of stencil coefficients.
     */
    public SimpleA33(float[][][] a) {
      _a = a;
    }
    public void get(int i1, int i2, float[] a) {
      float[] ai = _a[i2][i1];
      a[0] = ai[0];  a[1] = ai[1];  a[2] = ai[2];
      a[3] = ai[3];  a[4] = ai[4];  a[5] = ai[5];
      a[6] = ai[6];  a[7] = ai[7];  a[8] = ai[8];
    }
    private float[][][] _a;
  }

  /**
   * Constructs a multigrid solver.
   * @param n1 number of samples in 1st dimension of grid.
   * @param n2 number of samples in 2nd dimension of grid.
   * @param a33 a 3x3 stencil.
   * @param nbefore number of smoothings before downsampling.
   * @param ncycle number of recursive cycles at each coarse grid level.
   *  In the terminology of multigrid methods, ncycle=1 yields a V cycle, 
   *  and ncycle=2 yields a W cycle. Values greater than 2 are unusual.
   * @param nafter number of smoothings after upsampling.
   */
  public Multigrid2(
    int n1, int n2, A33 a33, 
    int nbefore, int ncycle, int nafter) 
  {
    _n1 = n1;
    _n2 = n2;
    _nlevel = nlevel(n1,n2);
    _nbefore = nbefore;
    _ncycle = ncycle;
    _nafter = nafter;

    // Array of operators, for the specified grid and coarser grids.
    _a33s = new A33[_nlevel];
    _a33s[_nlevel-1] = a33;
    for (int ilevel=_nlevel-2; ilevel>=0; --ilevel) {
      _a33s[ilevel] = coarsen(n1,n2,_a33s[ilevel+1]);
      n1 = (n1+1)/2;
      n2 = (n2+1)/2;
    }
  }

  /**
   * Updates the multigrid solution x of Ax = b with one cycle for the
   * finest grid level. Typically, this update corresponds to one V or
   * W cycle. If a good initial guess is not available, the solution x 
   * may be initially zero.
   * @param b array[n2][n1] for the right-hand-side.
   * @param x array[n2][n1] for the solution to be updated.
   */
  public void update(float[][] b, float[][] x) {
    cycleDownUp(_nlevel-1,b,x);
  }

  /**
   * Returns the sum of squared residuals r = b-Ax.
   * @param b array[n2][n1] for the right-hand-side.
   * @param x array[n2][n1] for the solution.
   */
  public float normResidual(float[][] b, float[][] x) {
    float[][] c = Array.copy(x);
    apply(_a33s[_nlevel-1],x,c);
    return normError(b,c);
  }

  /**
   * Returns the sum of squared errors e = x-y, where y is a known solution.
   * Useful only for test problems for which the solution is known.
   * @param x array[n2][n1] for the estimated solution.
   * @param y array[n2][n1] for the known solution.
   */
  public static float normError(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    double sum = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float e = x[i2][i1]-y[i2][i1];
        sum += e*e;
      }
    }
    return (float)sum;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2; // grid dimensions
  private int _nlevel; // number of multigrid levels
  private int _nbefore; // number of smoothings before downsampling
  private int _ncycle; // number of recursive cycles on coarse grids
  private int _nafter; // number of smoothings after upsampling
  private A33[] _a33s; // array of matrices, one for each level

  /**
   * Facilitates loops for 3x3 stencils. In effect, this class pads 2-D 
   * arrays with zeros on all sides, without making a complete copy of the 
   * entire array. Because typically only three consecutive columns of the 
   * array are needed at any time, this class maintains three columns that
   * are padded with one extra zero sample at each end.
   */
  private static class Buffer33 {

    /**
     * Constructs a buffer for the specified array.
     */
    Buffer33(float[][] a) {
      _n1 = a[0].length;
      _n2 = a.length;
      _a = a;
      _i = new int[]{-2,-2,-2};
      _b = new float[3][1+_n1+1];
    }

    /**
     * Gets an array for the specified column indexed by i2.
     * The returned array is padded with one extra zero sample at each end.
     */
    float[] get(int i2) {
      Check.argument(-1<=i2 && i2<=_n2,"index i2 is in bounds");
      int j2 = (i2+1)%3;
      if (_i[j2]!=i2) {
        if (0<=i2 && i2<_n2) {
          Array.copy(_n1,0,_a[i2],1,_b[j2]);
        } else {
          Array.zero(_b[j2]);
        }
        _i[j2] = i2;
      }
      return _b[j2];
    }

    private int _n1,_n2;
    private float[][] _a;
    private int[] _i;
    private float[][] _b;
  }

  /**
   * Computes the number of levels such that neither n1 nor n2 on the
   * coarsest grid is less than 3.
   */
  private static int nlevel(int n1, int n2) {
    int nlevel = 1;
    while ((n1>=5 && n2>=5)) {
      n1 = (n1+1)/2;
      n2 = (n2+1)/2;
      ++nlevel;
    }
    return nlevel;
  }

  
  /**
   * Applies the specified operator. Computes y = Ax.
   */
  private static void apply(A33 a33, float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    float[] a = new float[9];
    Buffer33 xb = new Buffer33(x);
    for (int i2=0; i2<n2; ++i2) {
      float[] xi2m = xb.get(i2-1);
      float[] xi20 = xb.get(i2  );
      float[] xi2p = xb.get(i2+1);
      for (int i1=0,j1=1; i1<n1; ++i1,++j1) {
        a33.get(i1,i2,a);
        y[i2][i1] = a[0]*xi2m[j1-1] + a[3]*xi20[j1-1] + a[6]*xi2p[j1-1] +
                    a[1]*xi2m[j1  ] + a[4]*xi20[j1  ] + a[7]*xi2p[j1  ] +
                    a[2]*xi2m[j1+1] + a[5]*xi20[j1+1] + a[8]*xi2p[j1+1];
      }
    }
  }

  /**
   * Computes the solution x of the system Ax = b. Computes x by LU
   * decomposition of A, which is costly unless the specified system is small.
   */
  private static void solve(A33 a33, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    int n = n1*n2;
    DMatrix am = new DMatrix(n,n);
    DMatrix bm = new DMatrix(n,1);
    float[] ai = new float[9];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        a33.get(i1,i2,ai);
        for (int k=0; k<9; ++k) {
          int k1 = k%3-1;
          int k2 = k/3-1;
          int j1 = i1+k1;
          int j2 = i2+k2;
          if (0<=j1 && j1<n1 && 0<=j2 && j2<n2) {
            int j = j1+j2*n1;
            am.set(i,j,ai[k]);
          }
        }
        bm.set(i,0,b[i2][i1]);
      }
    }
    DMatrixLud lud = new DMatrixLud(am);
    DMatrix xm = lud.solve(bm);
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        x[i2][i1] = (float)xm.get(i,0);
      }
    }
  }

  /**
   * Performs one cycle of the multigrid for the specified grid level.
   * This cycle may recursively perform multiple cycles for lower levels
   * corresponding to coarser grids. If the specified level is zero,
   * corresponding to the coarsest grid, then this method solves this
   * smallest system directly by LU decomposition.
   * 
   */
  private void cycleDownUp(int ilevel, float[][] b, float[][] x) {
    A33 a33 = _a33s[ilevel];
    int n1 = x[0].length;
    int n2 = x.length;

    // If coarsest grid, solve the coarsest system exactly.
    if (ilevel==0) {
      solve(a33,b,x);
    } 
    
    // Else, cycle recursively on coarser grids.
    else {

      // Smooth the solution x.
      for (int ibefore=0; ibefore<_nbefore; ++ibefore)
        smooth(a33,b,x);

      // Apply operator and compute residual.
      float[][] r = new float[n2][n1];
      apply(a33,x,r);
      Array.sub(b,r,r);

      // Downsample the residual.
      int m1 = (n1+1)/2;
      int m2 = (n2+1)/2;
      float[][] rc = new float[m2][m1];
      downsample(1.0f,r,rc);

      // Estimate error from downsampled residual on coarser grid.
      float[][] ec = new float[m2][m1];
      for (int icycle=0; icycle<_ncycle; ++icycle)
        cycleDownUp(ilevel-1,rc,ec);

      // Upsample the estimated error and accumulate in solution x.
      upsample(1.0f,ec,x);

      // Smooth the solution x.
      for (int iafter=0; iafter<_nafter; ++iafter)
        smooth(a33,b,x);
    }
  }

  /**
   * Weighted Jacobi relaxation with weight w = 2/3.
   */
  private static void smoothJacobi(A33 a33, float[][] b, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float w = 2.0f/3.0f;
    float omw = 1.0f-w;
    float[] a = new float[9];
    Buffer33 xb = new Buffer33(x);
    for (int i2=0; i2<n2; ++i2) {
      float[] xi2m = xb.get(i2-1);
      float[] xi20 = xb.get(i2  );
      float[] xi2p = xb.get(i2+1);
      for (int i1=0,j1=1; i1<n1; ++i1,++j1) {
        a33.get(i1,i2,a);
        float si = w/a[4];
        float ti = a[0]*xi2m[j1-1] + a[3]*xi20[j1-1] + a[6]*xi2p[j1-1] +
                   a[1]*xi2m[j1  ] +                   a[7]*xi2p[j1  ] +
                   a[2]*xi2m[j1+1] + a[5]*xi20[j1+1] + a[8]*xi2p[j1+1];
        x[i2][i1] = omw*x[i2][i1]+si*(b[i2][i1]-ti);
      }
    }
  }

  /**
   * Four-color Gauss-Seidel relaxation.
   */
  private static void smoothGaussSeidel4(A33 a33, float[][] b, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] a = new float[9];
    Buffer33 xb = new Buffer33(x);
    int[] k1s = {0,1,0,1};
    int[] k2s = {0,0,1,1};
    for (int k=0; k<4; ++k) {
      int k1 = k1s[k];
      int k2 = k2s[k];
      for (int i2=k2; i2<n2; i2+=2) {
        float[] xi2m = xb.get(i2-1);
        float[] xi20 = xb.get(i2  );
        float[] xi2p = xb.get(i2+1);
        for (int i1=k1,j1=1+k1; i1<n1; i1+=2,j1+=2) {
          a33.get(i1,i2,a);
          float ti = a[0]*xi2m[j1-1] + a[3]*xi20[j1-1] + a[6]*xi2p[j1-1] +
                     a[1]*xi2m[j1  ] +                   a[7]*xi2p[j1  ] +
                     a[2]*xi2m[j1+1] + a[5]*xi20[j1+1] + a[8]*xi2p[j1+1];
          float si = 1.0f/a[4];
          x[i2][i1] = si*(b[i2][i1]-ti);
        }
      }
    }
  }

  private static void smooth(A33 a33, float[][] b, float[][] x) {
    smoothJacobi(a33,b,x); // requires only pass
    //smoothGaussSeidel4(a33,b,x); // requires four passes
  }

  /**
   * Returns a coarsened version of the specified operator.
   * The specified operator has indices for an array[n2][n1].
   * The coarse operator has indices for an array[(n2+1)/2][(n1+1)/2].
   */
  private A33 coarsen(int n1, int n2, A33 a33) {
    int m1 = (n1+1)/2;
    int m2 = (n2+1)/2;
    float[][][] b = new float[m2][m1][9];
    float[] s = {0.25f,0.50f,0.25f,
                 0.50f,1.00f,0.50f,
                 0.25f,0.50f,0.25f};
    int[] index1 = {-1, 0, 1,
                    -1, 0, 1,
                    -1, 0, 1};
    int[] index2 = {-1,-1,-1,
                     0, 0, 0,
                     1, 1, 1};
    float[] aj = new float[9];
    for (int i2=0; i2<m2; ++i2) {
      for (int i1=0; i1<m1; ++i1) {
        float[] bi = b[i2][i1];
        for (int e=0; e<9; ++e) {
          int e1 = index1[e];
          int e2 = index2[e];
          int j1 = 2*i1+e1;
          int j2 = 2*i2+e2;
          if (0<=j1 && j1<n1 && 0<=j2 && j2<n2) {
            float se = 0.25f*s[e];
            a33.get(j1,j2,aj);
            for (int g=0; g<9; ++g) {
              int g1 = index1[g];
              int g2 = index2[g];
              for (int d=0; d<9; ++d) {
                int d1 = index1[d];
                int d2 = index2[d];
                int f1 = e1+g1-2*d1;
                int f2 = e2+g2-2*d2;
                if (-1<=f1 && f1<=1 && -1<=f2 && f2<=1) {
                  int f = 4+f1+3*f2;
                  bi[d] += se*s[f]*aj[g];
                }
              }
            }
          }
        }
      }
    }
    return new SimpleA33(b);
  }
 
  // Downsamples from [n1] x samples to [(n1+1)/2] y samples.
  // The gathering stencil is scale*[1/4,1/2,1/4].
  // Accumulates into the output array y.
  private static void downsample(float scale, float[] x, float[] y) {
    int n1 = x.length;
    int i1 = 0;
    int j1 = 0;
    float s = scale/4.0f;
    float t = x[i1];
    y[i1] += s*t;
    for (i1=1,j1=0; i1<n1; ++i1,j1=i1/2) {
      t = x[i1  ] +
          x[i1-1];
      y[j1] += s*t;
    }
  }

  // Downsample from [n2][n1] x samples to [(n2+1)/2][(n1+1)/2] y samples.
  //                                [1/16, 1/8, 1/16]
  // The gathering stencil is scale*[1/8,  1/4, 1/8 ]
  //                                [1/16, 1/8, 1/16].
  // Accumulates into the output array y.
  private static void xdownsample(float scale, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int m1 = y[0].length;
    int m2 = y.length;
    Check.argument(m1==(n1+1)/2,"m1=(n1+1)/2");
    Check.argument(m2==(n2+1)/2,"m2=(n2+1)/2");
    float s0 = scale*1.0f/4.0f;
    float s1 = scale*1.0f/8.0f;
    float s2 = scale*1.0f/16.0f;
    Buffer33 xb = new Buffer33(x);
    for (int i2=0,j2=0; i2<m2; ++i2,j2+=2) {
      float[] xj2m = xb.get(j2-1);
      float[] xj20 = xb.get(j2  );
      float[] xj2p = xb.get(j2+1);
      for (int i1=0,j1=1; i1<m1; ++i1,j1+=2) {
        y[i2][i1] = s0*(xj20[j1  ]) +
                    s1*(xj2m[j1  ]+xj2p[j1  ]+xj20[j1-1]+xj20[j1+1]) +
                    s2*(xj2m[j1-1]+xj2m[j1+1]+xj2p[j1-1]+xj2p[j1+1]);
      }
    }
  }

  // Downsample from [n2][n1] x samples to [(n2+1)/2][(n1+1)/2] y samples.
  //                                [1/16, 1/8, 1/16]
  // The gathering stencil is scale*[1/8,  1/4, 1/8 ]
  //                                [1/16, 1/8, 1/16].
  // Accumulates into the output array y.
  private static void downsample(float scale, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int m1 = y[0].length;
    int m2 = y.length;
    float s = scale/16.0f;

    // Rolling on.
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float t = x[i2][i1];
    y[j2][j1] += s*t;
    for (i1=1,j1=0; i1<n1; ++i1,j1=i1/2) {
      t = x[i2][i1  ] +
          x[i2][i1-1];
      y[j2][j1] += s*t;
    }
    if (j1<m1) {
      t = x[i2][i1-1];
      y[j2][j1] += s*t;
    }

    // Interior.
    for (i2=1,j2=0; i2<n2; ++i2,j2=i2/2) {
      i1 = j1 = 0;
      t = x[i2  ][i1  ] +
          x[i2-1][i1  ];
      y[j2][j1] += s*t;
      for (i1=1,j1=0; i1<n1; ++i1,j1=i1/2) {
        t = x[i2  ][i1  ] +
            x[i2  ][i1-1] +
            x[i2-1][i1  ] +
            x[i2-1][i1-1];
        y[j2][j1] += s*t;
      }
      if (j1<m1) {
        t = x[i2  ][i1-1] +
            x[i2-1][i1-1];
        y[j2][j1] += s*t;
      }
    }

    // Rolling off.
    if (j2<m2) {
      i1 = j1 = 0;
      t = x[i2-1][i1  ];
      y[j2][j1] += s*t;
      for (i1=1,j1=0; i1<n1; ++i1,j1=i1/2) {
        t = x[i2-1][i1  ] +
            x[i2-1][i1-1];
        y[j2][j1] += s*t;
      }
      if (j1<m1) {
        t = x[i2-1][i1-1];
        y[j2][j1] += s*t;
      }
    }
  }

  // Upsample from [(n1+1)/2] x samples to [n1] y samples.
  // The scattering stencil is scale*[1/2,1/1,1/2].
  // Accumulates into the output array y.
  private static void upsample(float scale, float[] x, float[] y) {
    int n1 = y.length;
    int i1 = 0;
    int j1 = 0;
    float s = scale/2.0f;
    float t = s*x[i1];
    y[j1] += t;
    for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
      t = s*x[i1];
      y[j1  ] += t;
      y[j1-1] += t;
    }
  }

  // Upsample from [(n2+1)/2][(n1+1)/2] x samples to [n2][n1] y samples.
  //                                 [1/4, 1/2, 1/4]
  // The scattering stencil is scale*[1/2, 1/1  1/2]
  //                                 [1/4, 1/2, 1/4].
  // Accumulates into the output array y.
  private static void upsample(float scale, float[][] x, float[][] y) {
    int m1 = x[0].length;
    int m2 = x.length;
    int n1 = y[0].length;
    int n2 = y.length;
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float s = scale/4.0f;
    float t = s*x[i2][i1];
    y[j2][j1] += t;
    for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
      t = s*x[i2][i1];
      y[j2][j1  ] += t;
      y[j2][j1-1] += t;
    }
    if (i1<m1) {
      t = s*x[i2][i1];
      y[j2][j1-1] += t;
    }
    for (j2=1,i2=0; j2<n2; ++j2,i2=j2/2) {
      i1 = j1 = 0;
      t = s*x[i2][i1];
      y[j2  ][j1] += t;
      y[j2-1][j1] += t;
      for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
        t = s*x[i2][i1];
        y[j2  ][j1  ] += t;
        y[j2  ][j1-1] += t;
        y[j2-1][j1  ] += t;
        y[j2-1][j1-1] += t;
      }
      if (i1<m1) {
        t = s*x[i2][i1];
        y[j2  ][j1-1] += t;
        y[j2-1][j1-1] += t;
      }
    }
    if (i2<m2) {
      i1 = j1 = 0;
      t = s*x[i2][i1];
      y[j2-1][j1] += t;
      for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
        t = s*x[i2][i1];
        y[j2-1][j1  ] += t;
        y[j2-1][j1-1] += t;
      }
      if (i1<m1) {
        t = s*x[i2][i1];
        y[j2-1][j1-1] += t;
      }
    }
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }
  private static void traceSequence(float[] x) {
    if (TRACE)
      edu.mines.jtk.mosaic.SimplePlot.asSequence(x);
  }
  private static void tracePixels(float[][] x) {
    if (TRACE) {
      trace("tracePixels: min="+Array.min(x)+" max="+Array.max(x));
      edu.mines.jtk.mosaic.SimplePlot.asPixels(x);
      /*
      edu.mines.jtk.mosaic.SimplePlot sp =
        new edu.mines.jtk.mosaic.SimplePlot();
      edu.mines.jtk.mosaic.PixelsView pv = sp.addPixels(x);
      pv.setInterpolation(
        edu.mines.jtk.mosaic.PixelsView.Interpolation.NEAREST);
      pv.setClips(-0.05f,0.05f);
      */
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  // Test code ensures that upsampling is the transpose of downsampling,
  // to within a known scale factor.
  private static void testDownUpSampling() {
    testDownUpSampling(107);
    testDownUpSampling(108);
    testDownUpSampling(107,107);
    testDownUpSampling(107,108);
    testDownUpSampling(108,107);
    testDownUpSampling(108,108);
  }
  private static void testDownUpSampling(int n) {
    System.out.println("testDownUpSampling: n="+n);
    int nx = n;
    int ny = (nx+1)/2;
    float[] x = Array.randfloat(nx);
    float[] y = Array.randfloat(ny);
    float[] ax = Array.zerofloat(ny);
    float[] ay = Array.zerofloat(nx);
    downsample(3.1f,x,ax);
    upsample(3.1f,y,ay);
    double xay = 0.0;
    for (int ix=0; ix<nx; ++ix)
      xay += x[ix]*ay[ix];
    double yax = 0.0;
    for (int iy=0; iy<ny; ++iy)
      yax += y[iy]*ax[iy];
    yax *= 2.0f;
    System.out.println("  xay="+xay);
    System.out.println("  yax="+yax);
  }
  private static void testDownUpSampling(int n1, int n2) {
    System.out.println("testDownUpSampling: n1="+n1+" n2="+n2);
    int n1x = n1;
    int n1y = (n1x+1)/2;
    int n2x = n2;
    int n2y = (n2x+1)/2;
    float[][] x = Array.randfloat(n1x,n2x);
    float[][] y = Array.randfloat(n1y,n2y);
    float[][] ax = Array.zerofloat(n1y,n2y);
    float[][] ay = Array.zerofloat(n1x,n2x);
    downsample(2.3f,x,ax);
    upsample(2.3f,y,ay);
    double xay = 0.0;
    for (int i2x=0; i2x<n2x; ++i2x)
      for (int i1x=0; i1x<n1x; ++i1x)
        xay += x[i2x][i1x]*ay[i2x][i1x];
    double yax = 0.0;
    for (int i2y=0; i2y<n2y; ++i2y)
      for (int i1y=0; i1y<n1y; ++i1y)
        yax += y[i2y][i1y]*ax[i2y][i1y];
    yax *= 4.0f;
    System.out.println("  xay="+xay);
    System.out.println("  yax="+yax);
  }

  private static void testSolve(int n) {
    int n1 = n;
    int n2 = n;
    float[][][] a = new float[n2][n1][9];
    float[][] b = new float[n2][n1]; // rhs b
    float[][] x = new float[n2][n1]; // multigrid solution
    float[][] y = new float[n2][n1]; // exact solution
    //loadBriggs(a,y,b);
    loadSimple(a,y,b);
    trace("y: min="+Array.min(y)+" max="+Array.max(y));
    tracePixels(y);

    A33 a33 = new Multigrid2.SimpleA33(a);
    Multigrid2 m2 = new Multigrid2(n1,n2,a33,0,2,2);
    float enew = m2.normError(x,y);
    float rnew = m2.normResidual(b,x);
    trace("initial e="+enew);
    trace("initial r="+rnew);
    int ncycle = 2;
    for (int icycle=0; icycle<ncycle; ++icycle) {
      m2.update(b,x);
      tracePixels(x);
      float eold = enew;
      float rold = rnew;
      enew = m2.normError(x,y);
      rnew = m2.normResidual(b,x);
      trace("  e="+enew+" ratio="+enew/eold);
      trace("  r="+rnew+" ratio="+rnew/rold);
    }
  }

  /**
   * Ax = b from Brigg's tutorial.
   */
  private static void loadBriggs(float[][][] a, float[][] x, float[][] b) {
    int n1 = x[0].length;
    int n2 = x.length;
    float d1 = 1.0f/(n1+1);
    float d2 = 1.0f/(n2+1);
    float d1s = d1*d1;
    float d2s = d2*d2;
    float f1 = d1;
    float f2 = d2;
    float[] ai = { 0.0f,    -1.0f/d2s,           0.0f,
                  -1.0f/d1s, 2.0f/d1s+2.0f/d2s, -1.0f/d1s,
                   0.0f,    -1.0f/d2s,           0.0f};
    for (int i2=0; i2<n2; ++i2) {
      float x2 = f2+i2*d2;
      float x2s = x2*x2;
      for (int i1=0; i1<n1; ++i1) {
        float x1 = f1+i1*d1;
        float x1s = x1*x1;
        for (int k=0; k<9; ++k)
          a[i2][i1][k] = ai[k];
        b[i2][i1] = 2.0f*((1.0f-6.0f*x1s)*x2s*(1.0f-x2s) +
                          (1.0f-6.0f*x2s)*x1s*(1.0f-x1s));
        x[i2][i1] = (x1s-x1s*x1s)*(x2s*x2s-x2s);
      }
    }
  }

  /**
   * Simple Ax = b for which x = s*b, where s is a constant scale factor.
   */
  private static void loadSimple(float[][][] a, float[][] x, float[][] b) {
    int n1 = x[0].length;
    int n2 = x.length;
    float d1 = 1.0f/(n1+1);
    float d2 = 1.0f/(n2+1);
    float d1s = d1*d1;
    float d2s = d2*d2;
    float f1 = d1;
    float f2 = d2;
    float[] ai = { 0.0f,    -1.0f/d2s,           0.0f,
                  -1.0f/d1s, 2.0f/d1s+2.0f/d2s, -1.0f/d1s,
                   0.0f,    -1.0f/d2s,           0.0f};
    for (int i2=0; i2<n2; ++i2) {
      float x2 = f2+i2*d2;
      float x2s = x2*x2;
      for (int i1=0; i1<n1; ++i1) {
        float x1 = f1+i1*d1;
        float x1s = x1*x1;
        for (int k=0; k<9; ++k)
          a[i2][i1][k] = ai[k];
        x[i2][i1] = sin(FLT_PI*x1)*sin(FLT_PI*x2);
        b[i2][i1] = 2.0f*FLT_PI*FLT_PI*x[i2][i1];
      }
    }
  }

  // Test code for multigrid.
  public static void main(String[] args) {
    // testDownUpSampling();
    testSolve(127);
  }
}
