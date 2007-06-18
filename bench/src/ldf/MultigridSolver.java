/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
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
      a[0] = ai[0]; a[1] = ai[1]; a[2] = ai[2];
      a[3] = ai[3]; a[4] = ai[4]; a[5] = ai[5];
      a[6] = ai[6]; a[7] = ai[7]; a[8] = ai[8];
    }
    private int _n1,_n2;
    private float[][][] _a;
  }

  public MultigridSolver(A33 a33) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
    _nlevel = nlevel(n1,n2);
    _n1 = n1;
    _n2 = n2;
    _a33s = new A33[_nlevel];
    _a33s[0] = a33;
    for (int ilevel=1; ilevel<_nlevel; ++ilevel)
      _a33s[ilevel] = coarsen(_a33s[ilevel-1]);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nlevel; // number of multigrid levels
  private int _n1,_n2; // problem dimensions
  private A33[] _a33s; // array of matrices, one for each level

  private static class Buffer33 {
    Buffer33(float[][] a) {
      _n1 = a[0].length;
      _n2 = a.length;
      _a = a;
      _i = new int[]{EMPTY,EMPTY,EMPTY};
      _b = new float[3][1+_n1+1];
    }
    float[] get(int i2) {
      Check.argument(-1<=i2 && i2<=_n2,"index i2 is in bounds");
      int j2 = (i2+1)%3;
      if (_i[j2]==EMPTY) {
        if (0<=i2 && i2<_n2) {
          Array.copy(_n1,0,_a[i2],1,_b[j2]);
        } else {
          Array.zero(_b[j2]);
        }
        _i[j2] = i2;
      }
      return _b[j2];
    }
    private static final int EMPTY = -2;
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

  private void vcycle(
    int ilevel, int nrelax1, int nrelax2, 
    float[][] b, float[][] x) 
  {
    int n1 = b[0].length;
    int n2 = b.length;
    int m1 = (n1+1)/2;
    int m2 = (n2+1)/2;
    if (ilevel==_nlevel) {
      // Solve coarsest system exactly.
      
    } else {
      A33 a33 = _a33s[ilevel];

      // Pre-smooth.
      for (int irelax=0; irelax<nrelax1; ++irelax)
        relaxJacobi(a33,b,x);

      // Compute and restrict residual.
      float[][] r = new float[n2][n1];
      apply(a33,x,r);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          r[i2][i1] = b[i2][i1]-r[i2][i1];
        }
      }
      float[][] bc = new float[m2][m1];
      downsample(1.0f,r,bc);

      // Recurse.
      float[][] xc = new float[m2][m1];
      vcycle(ilevel+1,nrelax1,nrelax2,bc,xc);

      // Interpolate and correct.
      upsample(1.0f,xc,r);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          x[i2][i1] += r[i2][i1];
        }
      }

      // Post-smooth.
      for (int irelax=0; irelax<nrelax1; ++irelax)
        relaxJacobi(a33,b,x);
    }
  }

  private void apply(A33 a33, float[][] x, float[][] y) {
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

  private void relaxJacobi(A33 a33, float[][] b, float[][] v) {
    int n1 = a33.getN1();
    int n2 = a33.getN2();
    float w = 2.0f/3.0f;
    float omw = 1.0f-w;
    float[] a = new float[9];
    Buffer33 vb = new Buffer33(v);
    for (int i2=0; i2<n2; ++i2) {
      float[] vi2m = vb.get(i2-1);
      float[] vi20 = vb.get(i2  );
      float[] vi2p = vb.get(i2+1);
      for (int i1=0,j1=1; i1<n1; ++i1,++j1) {
        a33.getA(i1,i2,a);
        float si = w/a[4];
        float ti = a[0]*vi2m[j1-1] + a[3]*vi20[j1-1] + a[6]*vi2p[j1-1] +
                   a[1]*vi2m[j1  ] +                   a[7]*vi2p[j1  ] +
                   a[2]*vi2m[j1+1] + a[5]*vi20[j1+1] + a[8]*vi2p[j1+1];
        v[i2][i1] = omw*v[i2][i1]+si*(b[i2][i1]-ti);
      }
    }
  }

  private A33 coarsen(A33 a) {
    int n1 = a.getN1();
    int n2 = a.getN2();
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
        for (int d=0; d<9; ++d) {
          int d1 = index1[d];
          int d2 = index2[d];
          int j1 = 2*i1+d1;
          int j2 = 2*i2+d2;
          if (0<=j1 && j1<n1 && 0<=j2 && j2<n2) {
            a.getA(j1,j2,aj);
            for (int e=0; e<9; ++e) {
              int e1 = index1[e];
              int e2 = index2[e];
              float se = s[e];
              for (int g=0; g<9; ++g) {
                int g1 = index1[g];
                int g2 = index2[g];
                int f1 = e1+g1-2*d1;
                int f2 = e2+g2-2*d2;
                int f = 4+f1+3*f2;
                if (0<=f && f<9)
                  bi[d] += se*s[f]*aj[g];
              }
            }
          }
        }
        for (int i=0; i<9; ++i)
          bi[i] *= 0.25f;
      }
    }
    return new SimpleA33(b);
  }
 
  // Downsamples from [n1] x samples to [(n1+1)/2] y samples.
  // The gathering stencil is scale*[1/4,1/2,1/4].
  private static void downsample(float scale, float[] x, float[] y) {
    Array.zero(y);
    int n1 = x.length;
    int i1 = 0;
    int j1 = 0;
    float s = scale/4.0f;
    float t = x[i1];
    y[i1] = s*t;
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
  private static void downsample(float scale, float[][] x, float[][] y) {
    Array.zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float s = scale/16.0f;
    float t = x[i2][i1];
    y[j2][j1] = s*t;
    for (i1=1,j1=0; i1<n1; ++i1,j1=i1/2) {
      t = x[i2][i1  ] +
          x[i2][i1-1];
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
    }
  }

  // Upsample from [(n1+1)/2] x samples to [n1] y samples.
  // The scattering stencil is scale*[1/2,1/1,1/2].
  private static void upsample(float scale, float[] x, float[] y) {
    int n1 = y.length;
    int i1 = 0;
    int j1 = 0;
    float s = scale/2.0f;
    float t = s*x[i1];
    y[j1]  = t;
    for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
      t = s*x[i1];
      y[j1  ]  = t;
      y[j1-1] += t;
    }
  }

  // Upsample from [(n2+1)/2][(n1+1)/2] x samples to [n2][n1] y samples.
  //                                 [1/4, 1/2, 1/4]
  // The scattering stencil is scale*[1/2, 1/1  1/2]
  //                                 [1/4, 1/2, 1/4].
  private static void upsample(float scale, float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float s = scale/4.0f;
    float t = s*x[i2][i1];
    y[j2][j1]  = t;
    for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
      t = s*x[i2][i1];
      y[j2][j1  ]  = t;
      y[j2][j1-1] += t;
    }
    for (j2=1,i2=0; j2<n2; ++j2,i2=j2/2) {
      i1 = j1 = 0;
      t = s*x[i2][i1];
      y[j2  ][j1]  = t;
      y[j2-1][j1] += t;
      for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
        t = s*x[i2][i1];
        y[j2  ][j1  ]  = t;
        y[j2  ][j1-1] += t;
        y[j2-1][j1  ] += t;
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
    if (TRACE)
      edu.mines.jtk.mosaic.SimplePlot.asPixels(x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  // Test code ensures that upsampling is the transpose of downsampling,
  // to within a known scale factor.
  private static void mainTestDownUpSampling(String[] args) {
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

  // Test code for multigrid.
  public static void main(String[] args) {
    testMg(200,1,1000);
    testMg(16,1,1000);
    testMg(5,4,1000);
  }
  private static void testMg(int niter, int nlevel, int n1) {
    float[] d = new float[n1];
    float[] x = new float[n1];
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      d[i1] = 1.0f+256.0f*sin(FLT_PI*(float)i1/(float)(n1-1));
      if (i1>0 && i1%100==0)
        x[i1] = 1.0f;
    }
    float sigma = 1.0f;
    float small = 0.0001f;
    LocalDiffusionFilter ldf = 
      new LocalDiffusionFilter(sigma,small,niter,nlevel);
    ldf.applyMg(d,x,y);
    traceSequence(y);
  }
} 
