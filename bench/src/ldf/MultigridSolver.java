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
 * Multigrid solver.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.17
 */
public class MultigridSolver {

  public interface A33 {
    public int getN1();
    public int getN2();
    public void getA(int i1, int i2, float[] a);
  }

  public static class SimpleA33 implements A33 {
    public SimpleA33(float[][][] a) {
      _n1 = a[0].length;
      _n2 = a.length;
      _a = a;
    }
    public int getN1() {
      return _n1;
    }
    public int getN2() {
      return _n2;
    }
    public void getA(int i1, int i2, float[] a) {
      float[] ai = _a[i2][i1];
      a[0] = ai[0];  a[1] = ai[1];  a[2] = ai[2];
      a[3] = ai[3];  a[4] = ai[4];  a[5] = ai[5];
      a[6] = ai[6];  a[7] = ai[7];  a[8] = ai[8];
    }
    private int _n1,_n2;
    private float[][][] _a;
  }

  /**
   * Constructs a multigrid solver for 3x3 stencil.
   * @param a33 the linear operator.
   * @param nbefore number of smoothings before downsampling
   * @param ncycle number of recursive cycles for coarse grids
   * @param nafter number smoothings after upsampling
   */
  public MultigridSolver(A33 a33, int nbefore, int ncycle, int nafter) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
    _nlevel = nlevel(n1,n2);
    _nbefore = nbefore;
    _ncycle = ncycle;
    _nafter = nafter;
    _n1 = n1;
    _n2 = n2;
    _a33s = new A33[_nlevel];
    _a33s[_nlevel-1] = a33;
    for (int ilevel=_nlevel-2; ilevel>=0; --ilevel)
      _a33s[ilevel] = coarsen(_a33s[ilevel+1]);
  }

  /**
   * Updates the multigrid solution x of Ax = b.
   * @param b the right-hand-side.
   * @param x the updated solution.
   */
  public void solve(float[][] b, float[][] x) {
    cycleDownUp(_nlevel-1,b,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nlevel; // number of multigrid levels
  private int _nbefore; // number of smoothings before downsampling
  private int _ncycle; // number of recursive cycles on coarse grids
  private int _nafter; // number of smoothings after upsampling
  private int _n1,_n2; // problem dimensions
  private A33[] _a33s; // array of matrices, one for each level

  private static class Buffer33 {
    Buffer33(float[][] a) {
      _n1 = a[0].length;
      _n2 = a.length;
      _a = a;
      _i = new int[]{-2,-2,-2};
      _b = new float[3][1+_n1+1];
    }
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

  private static int nlevel(int n1, int n2) {
    int nlevel = 1;
    while ((n1==1 || n1>=5) && 
           (n2==1 || n2>=5)) {
      n1 = (n1+1)/2;
      n2 = (n2+1)/2;
      ++nlevel;
    }
    return nlevel;
  }

  private static void apply(A33 a33, float[][] x, float[][] y) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
    float[] a = new float[9];
    Buffer33 xb = new Buffer33(x);
    for (int i2=0; i2<n2; ++i2) {
      float[] xi2m = xb.get(i2-1);
      float[] xi20 = xb.get(i2  );
      float[] xi2p = xb.get(i2+1);
      for (int i1=0,j1=1; i1<n1; ++i1,++j1) {
        a33.getA(i1,i2,a);
        y[i2][i1] = a[0]*xi2m[j1-1] + a[3]*xi20[j1-1] + a[6]*xi2p[j1-1] +
                    a[1]*xi2m[j1  ] + a[4]*xi20[j1  ] + a[7]*xi2p[j1  ] +
                    a[2]*xi2m[j1+1] + a[5]*xi20[j1+1] + a[8]*xi2p[j1+1];
      }
    }
  }

  private static void solve(A33 a33, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    int n = n1*n2;
    DMatrix am = new DMatrix(n,n);
    DMatrix bm = new DMatrix(n,1);
    float[] ai = new float[9];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        a33.getA(i1,i2,ai);
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

  private void cycleDownUp(int ilevel, float[][] b, float[][] x) {
    A33 a33 = _a33s[ilevel];
    int n1 = a33.getN1();
    int n2 = a33.getN2();

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

  private static void smooth(A33 a33, float[][] b, float[][] x) {
    smoothJacobi(a33,b,x); // requires only pass
    //smoothGaussSeidel4(a33,b,x); // requires four passes
  }

  private static float norm2(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        sum += xi*xi;
      }
    }
    return (float)sum;
  }

  private static void traceResidual(A33 a33, float[][] b, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1];
    apply(a33,x,r);
    Array.sub(b,r,r);
    tracePixels(r);
  }

  private static float residual(A33 a33, float[][] b, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1];
    apply(a33,x,r);
    Array.sub(b,r,r);
    return norm2(r);
  }

  // Weighted Jacobi relaxation with weight w = 2/3.
  private static void smoothJacobi(A33 a33, float[][] b, float[][] x) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
    float w = 2.0f/3.0f;
    float omw = 1.0f-w;
    float[] a = new float[9];
    Buffer33 xb = new Buffer33(x);
    for (int i2=0; i2<n2; ++i2) {
      float[] xi2m = xb.get(i2-1);
      float[] xi20 = xb.get(i2  );
      float[] xi2p = xb.get(i2+1);
      for (int i1=0,j1=1; i1<n1; ++i1,++j1) {
        a33.getA(i1,i2,a);
        float si = w/a[4];
        float ti = a[0]*xi2m[j1-1] + a[3]*xi20[j1-1] + a[6]*xi2p[j1-1] +
                   a[1]*xi2m[j1  ] +                   a[7]*xi2p[j1  ] +
                   a[2]*xi2m[j1+1] + a[5]*xi20[j1+1] + a[8]*xi2p[j1+1];
        x[i2][i1] = omw*x[i2][i1]+si*(b[i2][i1]-ti);
      }
    }
  }

  // Four-color Gauss-Seidel relaxation.
  private static void smoothGaussSeidel4(A33 a33, float[][] b, float[][] x) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
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
          a33.getA(i1,i2,a);
          float ti = a[0]*xi2m[j1-1] + a[3]*xi20[j1-1] + a[6]*xi2p[j1-1] +
                     a[1]*xi2m[j1  ] +                   a[7]*xi2p[j1  ] +
                     a[2]*xi2m[j1+1] + a[5]*xi20[j1+1] + a[8]*xi2p[j1+1];
          float si = 1.0f/a[4];
          x[i2][i1] = si*(b[i2][i1]-ti);
        }
      }
    }
  }

  // Returns a coarsened version of the specified operator.
  private A33 coarsen(A33 a33) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
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
            a33.getA(j1,j2,aj);
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
    //trace("coarsen: n1="+n1+" n2="+n2);
    //a33.getA(n1/2,n2/2,aj);
    //Array.dump(aj);
    //Array.dump(b[m2/2][m1/2]);
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
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float s = scale/16.0f;
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
      trace("x: min="+Array.min(x)+" max="+Array.max(x));
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
    float d1 = 1.0f/(n1+1);
    float d2 = 1.0f/(n2+1);
    float d1s = d1*d1;
    float d2s = d2*d2;
    float f1 = d1;
    float f2 = d2;
    float[] ai = { 0.0f,    -1.0f/d2s,           0.0f,
                  -1.0f/d1s, 2.0f/d1s+2.0f/d2s, -1.0f/d1s,
                   0.0f,    -1.0f/d2s,           0.0f};
    //float[] ai = { 0.0f, 0.0f, 0.0f,
    //               0.0f, 1.0f, 0.0f,
    //               0.0f, 0.0f, 0.0f};
    float[][][] a = new float[n2][n1][9];
    float[][] b = new float[n2][n1]; // rhs b
    float[][] c = new float[n2][n1]; // c = Ax
    float[][] d = new float[n2][n1]; // d = Ay
    float[][] x = new float[n2][n1]; // multigrid solution
    float[][] y = new float[n2][n1]; // exact solution

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
        y[i2][i1] = (x1s-x1s*x1s)*(x2s*x2s-x2s);
      }
    }
    /*
    for (int i2=0; i2<n2; ++i2) {
      float x2 = f2+i2*d2;
      float x2s = x2*x2;
      for (int i1=0; i1<n1; ++i1) {
        float x1 = f1+i1*d1;
        float x1s = x1*x1;
        for (int k=0; k<9; ++k)
          a[i2][i1][k] = ai[k];
        b[i2][i1] = 2.0f*FLT_PI*FLT_PI*sin(FLT_PI*x1)*sin(FLT_PI*x2);
        y[i2][i1] = sin(FLT_PI*x1)*sin(FLT_PI*x2);
      }
    }
    */
    trace("y: min="+Array.min(y)+" max="+Array.max(y));
    tracePixels(y);

    A33 a33 = new MultigridSolver.SimpleA33(a);
    MultigridSolver ms = new MultigridSolver(a33,0,2,2);
    int ncycle = 2;
    float rnew = residual(a33,b,x);
    trace("initial r="+rnew);
    trace("  |x-y|^2 = "+norm2(Array.sub(x,y)));
    for (int icycle=0; icycle<ncycle; ++icycle) {
      ms.solve(b,x);
      tracePixels(x);
      float rold = rnew;
      rnew = residual(a33,b,x);
      trace("  r="+rnew+" ratio="+rnew/rold);
      //trace("  min="+Array.min(x)+" max="+Array.max(x));
      trace("  |x-y|^2 = "+norm2(Array.sub(x,y)));
    }
  }

  // Test code for multigrid.
  public static void main(String[] args) {
    // testDownUpSampling();
    testSolve(127);
  }
}
