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
 * Local anisotropic diffusion filtering via a multigrid method.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.05
 */
public class LocalDiffusionFilterMg {

  public LocalDiffusionFilterMg(double sigma) {
    this(sigma,0.000001,100,2,2,2);
  }

  public LocalDiffusionFilterMg(
    double sigma, double small, int niter,
    int nbefore, int ncycle, int nafter) 
  {
    _dlf = new DirectionalLaplacianFilter(sigma);
    _sigma = (float)sigma;
    _small = (float)small;
    _niter = niter;
    _nbefore = nbefore;
    _ncycle = ncycle;
    _nafter = nafter;
  }

  /**
   * Applies a line smoothing filter using a multigrid method.
   * Diffusivities parallel to linear features are one.
   * Diffusivities perpendicular to linear features are zero.
   * @param v1 array of 1st como
   * @param u2 array of 2nd components of vectors normal to linear features.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyLineSmoothing(float[][] u2, float[][] x, float[][] y) {
    int n1 = u2[0].length;
    int n2 = u2.length;
    float[][] d11 = new float[n2][n1];
    float[][] d12 = new float[n2][n1];
    float[][] d22 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float v2i =  u1i;
        float v1i = -u2i;
        d11[i2][i1] = v1i*v1i;
        d12[i2][i1] = v1i*v2i;
        d22[i2][i1] = v2i*v2i;
      }
    }
    float[][][] d = {d11,d12,d22};
    Array.copy(x,y);
    solveMg(d,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // maximum filter half-width
  private float _small; // stop iterations when rr decreases by this factor
  private int _niter; // number of iterations
  private int _nbefore; // number of mg smoothings before downsampling
  private int _ncycle; // number of mg recursive cycles on coarse grids
  private int _nafter; // number of mg smoothings after upsampling

  /**
   * Solves the diffusion system via multigrid iterations.
   */
  private void solveMg(float[][][] d, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    Multigrid2.A33 a33 = makeOperator(d);
    Multigrid2 m2 = new Multigrid2(n1,n2,a33,_nbefore,_ncycle,_nafter);
    float rr = m2.normResidual(x,y);
    trace("solveMg: r="+rr);
    int miter;
    double rrsmall = rr*_small;
    for (miter=0; miter<_niter && rr>rrsmall; ++miter) {
      m2.update(x,y);
      rr = m2.normResidual(x,y);
      //trace("  miter="+miter+" rr="+rr);
    }
    trace("  miter="+miter+" rr="+rr);
  }

  /**
   * Makes the diffusion operator for the multigrid method.
   */
  private Multigrid2.A33 makeOperator(float[][][] d) {
    int n1 = d[0][0].length;
    int n2 = d[0].length;

    // Coefficients for finite-difference approximation.
    //float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    //float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    float r = 0.5f;
    float s = 0.5f;

    // Array of stencil coefficients.
    float[][][] a = new float[n2][n1][9];

    // Arrays of ones and zeros to pick out stencil coefficients.
    // Array x00s picks out the coefficient of x[i2  ][i1  ].
    // Array x01s picks out the coefficient of x[i2  ][i1-1].
    // Array x10s picks out the coefficient of x[i2-1][i1  ].
    // Array x11s picks out the coefficient of x[i2-1][i1-1].
    float[][] x00s = {{1.0f,0.0f},{0.0f,0.0f}};
    float[][] x01s = {{0.0f,1.0f},{0.0f,0.0f}};
    float[][] x10s = {{0.0f,0.0f},{1.0f,0.0f}};
    float[][] x11s = {{0.0f,0.0f},{0.0f,1.0f}};

    // Arrays of components of diffusion tensor.
    float[][] d11s = d[0];
    float[][] d12s = d[1];
    float[][] d22s = d[2];

    // Scale factor applied to diffusion tensor.
    float ss = 0.5f*_sigma*_sigma;

    // Identity matrix.
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[i2][i1][4] = 1.0f;
      }
    }

    // For all samples, ...
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float[] a00s = a[i2  ][i1  ];
        float[] a01s = a[i2  ][i1-1];
        float[] a10s = a[i2-1][i1  ];
        float[] a11s = a[i2-1][i1-1];

        // Average diffusion tensor.
        float d11 = 0.25f*ss*(d11s[i2  ][i1  ]+d11s[i2  ][i1-1] +
                              d11s[i2-1][i1  ]+d11s[i2-1][i1-1]);
        float d12 = 0.25f*ss*(d12s[i2  ][i1  ]+d12s[i2  ][i1-1] +
                              d12s[i2-1][i1  ]+d12s[i2-1][i1-1]);
        float d22 = 0.25f*ss*(d22s[i2  ][i1  ]+d22s[i2  ][i1-1] +
                              d22s[i2-1][i1  ]+d22s[i2-1][i1-1]);

        // Accumulate stencil coefficients.
        for (int k2=0; k2<2; ++k2) {
          for (int k1=0; k1<2; ++k1) {
            float x00 = x00s[k2][k1];
            float x01 = x01s[k2][k1];
            float x10 = x10s[k2][k1];
            float x11 = x11s[k2][k1];
            float x1 = r*(x00-x01)+s*(x10-x11);
            float x2 = r*(x00-x10)+s*(x01-x11);
            float y1 = d11*x1+d12*x2;
            float y2 = d12*x1+d22*x2;
            a00s[4-k1-3*k2] += r*y1+r*y2; // updates a[4],a[3],a[1],a[0]
            a01s[5-k1-3*k2] -= r*y1-s*y2; // updates a[5],a[4],a[2],a[1]
            a10s[7-k1-3*k2] += s*y1-r*y2; // updates a[7],a[6],a[4],a[3]
            a11s[8-k1-3*k2] -= s*y1+s*y2; // updates a[8],a[7],a[5],a[4]
          }
        }
      }
    }
    //Array.dump(a[n2/2][n1/2]);
    return new Multigrid2.SimpleA33(a);
  }

  /**
   * Solves the diffusion system via conjugate gradient iterations.
   */
  private void solveCg(float[][][] d, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    float[][] t = new float[n2][n1];
    applyOperator(d,y,t);
    double rr = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] r2 = r[i2];
      float[] s2 = s[i2];
      float[] t2 = t[i2];
      for (int i1=0; i1<n1; ++i1) {
        float ri = x2[i1]-t2[i1];
        r2[i1] = ri;
        s2[i1] = ri;
        rr += ri*ri;
      }
    }
    trace("solveCg: r="+rr);
    int miter;
    double rrsmall = rr*_small;
    for (miter=0; miter<_niter && rr>rrsmall; ++miter) {
      applyOperator(d,s,t);
      double st = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] s2 = s[i2];
        float[] t2 = t[i2];
        for (int i1=0; i1<n1; ++i1)
          st += s2[i1]*t2[i1];
      }
      float alpha = (float)(rr/st);
      double rrold = rr;
      rr = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] y2 = y[i2];
        float[] r2 = r[i2];
        float[] s2 = s[i2];
        float[] t2 = t[i2];
        for (int i1=0; i1<n1; ++i1) {
          y2[i1] += alpha*s2[i1];
          r2[i1] -= alpha*t2[i1];
          rr += r2[i1]*r2[i1];
        }
      }
      if (rr<=rrsmall)
        break;
      float beta = (float)(rr/rrold);
      for (int i2=0; i2<n2; ++i2) {
        float[] r2 = r[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = r2[i1]+beta*s2[i1];
      }
    }
    trace("  miter="+miter+" rr="+rr);
  }

  private void applyOperator(float[][][] d, float[][] x, float[][] y) {
    Array.copy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float ss = 0.5f*_sigma*_sigma;
    //float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    //float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    float r = 0.5f;
    float s = 0.5f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float d11 = ss*d[0][i2][i1];
        float d12 = ss*d[1][i2][i1];
        float d22 = ss*d[2][i2][i1];
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float x1 = r*(x00-x01)+s*(x10-x11);
        float x2 = r*(x00-x10)+s*(x01-x11);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        y[i2  ][i1  ] += r*y1+r*y2;
        y[i2  ][i1-1] -= r*y1-s*y2;
        y[i2-1][i1  ] += s*y1-r*y2;
        y[i2-1][i1-1] -= s*y1+s*y2;
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
} 
