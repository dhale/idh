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
 * Local anisotropic diffusion filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.04.08
 */
public class LocalDiffusionFilter {

  public LocalDiffusionFilter(double sigma) {
    this(sigma,0.0001,5,4);
  }

  public LocalDiffusionFilter(
    double sigma, double small, int niter, int nlevel) 
  {
    _sigma = (float)sigma;
    _small = (float)small;
    _niter = niter;
    _nlevel = nlevel;
  }

  /**
   * Applies this local anisotropic diffusion filter.
   * Input and output arrays must be distinct.
   * @param su diffusivities in direction of normal vector u.
   * @param sv diffusivities in direction of normal vector v.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void apply(
    float[][] su, float[][] sv, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    int ns = 1+(int)(_sigma*_sigma);
    float ss = 0.5f*(_sigma*_sigma)/(float)ns;
    float[] xi20 = new float[n1];
    float[] xi21 = new float[n1];
    for (int is=0; is<ns; ++is,x=y) {
      for (int i2=0; i2<n2; ++i2) {
        float[] xtmp = xi20;  xi20 = xi21;  xi21 = xtmp;
        for (int i1=0; i1<n1; ++i1)
          xi20[i1] = x[i2][i1];
        for (int i1=0; i1<n1; ++i1) {
          float sui = su[i2][i1];
          float svi = sv[i2][i1];
          float u2i = u2[i2][i1];
          float u1i = sqrt(1.0f-u2i*u2i);
          float v2i =  u1i;
          float v1i = -u2i;
          float a11 = ss*(sui*u1i*u1i+svi*v1i*v1i);
          float a12 = ss*(sui*u1i*u2i+svi*v1i*v2i);
          float a22 = ss*(sui*u2i*u2i+svi*v2i*v2i);
          float x00 = xi20[i1];
          float x10 = xi21[i1];
          float x01 = (i1>0)?xi20[i1-1]:0.0f;
          float x11 = (i1>0)?xi21[i1-1]:0.0f;
          float xa = x00-x11;
          float xb = x01-x10;
          float x1 = 0.5f*(xa-xb);
          float x2 = 0.5f*(xa+xb);
          float y1 = a11*x1+a12*x2;
          float y2 = a12*x1+a22*x2;
          float ya = 0.5f*(y1+y2);
          float yb = 0.5f*(y1-y2);
          y[i2][i1] = x00-ya;
          if (i1>0) y[i2][i1-1] += yb;
          if (i2>0) y[i2-1][i1] -= yb;
          if (i2>0 && i1>0) y[i2-1][i1-1] += ya;
        }
      }
    }
  }

  /**
   * Applies a dip smoothing filter using a multigrid method.
   * @param u2 array of 2nd components of vectors normal to dip.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyDipSmoothing(float[][] u2, float[][] x, float[][] y) {
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
    //applyMg(d,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // filter half-width
  private float _small; // small value used to terminate CG iterations
  private int _niter; // number of CG iterations per multigrid level
  private int _nlevel; // number of CG iterations per multigrid level

  // Solves the diffusion system via conjugate gradient iterations.
  private void solveCg(float[] d, float[] x, float[] y) {
    int n1 = x.length;
    float[] r = new float[n1];
    float[] s = new float[n1];
    float[] t = new float[n1];
    applyOperator(d,y,t);
    double rr = 0.0;
    for (int i1=0; i1<n1; ++i1) {
      float ri = x[i1]-t[i1];
      r[i1] = ri;
      s[i1] = ri;
      rr += ri*ri;
    }
    trace("solveCg: n1="+n1+" rr="+rr);
    double rrsmall = rr*_small;
    int miter;
    for (miter=0; miter<_niter && rr>rrsmall; ++miter) {
      applyOperator(d,s,t);
      double st = 0.0;
      for (int i1=0; i1<n1; ++i1)
        st += s[i1]*t[i1];
      float alpha = (float)(rr/st);
      double rrold = rr;
      rr = 0.0;
      for (int i1=0; i1<n1; ++i1) {
        y[i1] += alpha*s[i1];
        r[i1] -= alpha*t[i1];
        rr += r[i1]*r[i1];
      }
      if (rr<=rrsmall)
        break;
      float beta = (float)(rr/rrold);
      for (int i1=0; i1<n1; ++i1)
        s[i1] = r[i1]+beta*s[i1];
    }
    trace("  miter="+miter+" rr="+rr);
  }

  // Solves the diffusion system via conjugate gradient iterations.
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
    trace("solveCg: n1="+n1+" rr="+rr);
    double rrsmall = rr*_small;
    int miter;
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

  // Applies the forward operator. Diffusion filtering is the inverse of
  // this operator. The operator is I + 0.5*sigma*sigma*G'DD'G, where G 
  // approximates the gradient operator and D is the diffusion tensor.
  private void applyOperator(float[] d, float[] x, float[] y) {
    Array.copy(x,y);
    int n1 = x.length;
    float ss = 0.5f*_sigma*_sigma;
    for (int i1=0; i1<n1; ++i1) {
      float d11 = ss*d[i1];
      float xi0 = x[i1];
      float xi1 = (i1>0)?x[i1-1]:0.0f;
      float x1 = xi0-xi1;
      float y1 = d11*x1;
      y[i1] += y1;
      if (i1>0) y[i1-1] -= y1;
    }
  }

  // Applies the forward operator. Diffusion filtering is the inverse of
  // this operator. The operator is I + 0.5*sigma*sigma*G'DD'G, where G 
  // approximates the gradient operator and D is the diffusion tensor.
  private void xapplyOperator(float[][][] d, float[][] x, float[][] y) {
    Array.copy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float ss = 0.5f*_sigma*_sigma;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float d11 = ss*d[0][i2][i1];
        float d12 = ss*d[1][i2][i1];
        float d22 = ss*d[2][i2][i1];
        float x00 = x[i2][i1];
        float x01 = (i1>0)?x[i2][i1-1]:0.0f;
        float x10 = (i2>0)?x[i2-1][i1]:0.0f;
        float x11 = (i2>0 && i1>0)?x[i2-1][i1-1]:0.0f;
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.5f*(xa-xb);
        float x2 = 0.5f*(xa+xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2][i1] += ya;
        if (i1>0) y[i2][i1-1] -= yb;
        if (i2>0) {
          y[i2-1][i1] += yb;
          if (i1>0) y[i2-1][i1-1] -= ya;
        }
      }
    }
  }
  private void applyOperator(float[][][] d, float[][] x, float[][] y) {
    Array.copy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float ss = 0.5f*_sigma*_sigma;
    float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float d11 = ss*d[0][i2][i1];
        float d12 = ss*d[1][i2][i1];
        float d22 = ss*d[2][i2][i1];
        float x00 = x[i2][i1];
        float x01 = (i1>0)?x[i2][i1-1]:0.0f;
        float x10 = (i2>0)?x[i2-1][i1]:0.0f;
        float x11 = (i2>0 && i1>0)?x[i2-1][i1-1]:0.0f;
        float x1 = r*(x00-x01)+s*(x10-x11);
        float x2 = r*(x00-x10)+s*(x01-x11);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        y[i2][i1] += r*y1+r*y2;
        if (i1>0) y[i2][i1-1] -= r*y1-s*y2;
        if (i2>0) y[i2-1][i1  ] += s*y1-r*y2;
        if (i2>0 && i1>0) y[i2-1][i1-1] -= s*y1+s*y2;
      }
    }
  }

  private static float dot(float[] x, float[] y) {
    int n1 = x.length;
    double s = 0.0;
    for (int i1=0; i1<n1; ++i1)
      s += x[i1]*y[i1];
    return (float)s;
  }
  private static float dot(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    double s = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        s += x2[i1]*y2[i1];
    }
    return (float)s;
  }
  private static float dot(float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    double s = 0.0;
    for (int i3=0; i3<n3; ++i3) {
      float[][] x3 = x[i3];
      float[][] y3 = y[i3];
      for (int i2=0; i2<n2; ++i2) {
        float[] x32 = x3[i2];
        float[] y32 = y3[i2];
        for (int i1=0; i1<n1; ++i1)
          s += x32[i1]*y32[i1];
      }
    }
    return (float)s;
  }

  private static void saxpy(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += a*x[i1];
  }
  private static void saxpy(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        y2[i1] += a*x2[i1];
    }
  }
  private static void saxpy(float a, float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      float[][] x3 = x[i3];
      float[][] y3 = y[i3];
      for (int i2=0; i2<n2; ++i2) {
        float[] x32 = x3[i2];
        float[] y32 = y3[i2];
        for (int i1=0; i1<n1; ++i1)
          y32[i1] += a*x32[i1];
      }
    }
  }

  private static void upsampleComplicated(float[] x, float[] y) {
    int n1x = x.length;
    int n1y = y.length;
    int i1xmax = min(n1x-1,n1y/2-1);
    int i1x = 0;
    int i1y = 0;
    float xi = x[i1x];
    float xio2 = 0.5f*xi;
    y[i1y  ] = xi;
    y[i1y+1] = xio2;
    for (i1x=1,i1y=2; i1x<=i1xmax; ++i1x,i1y+=2) {
      xi = x[i1x];
      xio2 = 0.5f*xi;
      y[i1y-1] += xio2;
      y[i1y  ]  = xi;
      y[i1y+1]  = xio2;
    }
    if (i1x<n1x && i1y<n1y) {
      xi = x[i1x];
      xio2 = 0.5f*xi;
      y[i1y-1] += xio2;
      y[i1y  ]  = xi;
    }
  }

  private static void upsampleComplicated(float[][] x, float[][] y) {
    int n1x = x[0].length;
    int n1y = y[0].length;
    int n2x = x.length;
    int n2y = y.length;
    int i1xmax = min(n1x-1,n1y/2-1);
    int i2xmax = min(n2x-1,n2y/2-1);
    int i1x = 0;
    int i1y = 0;
    int i2x = 0;
    int i2y = 0;

    // i2x = i2y = 0
    float[] y2m = null;
    float[] y20 = y[i2y  ];
    float[] y2p = y[i2y+1];
    float xi = x[i2x][i1x];
    float xio2 = 0.50f*xi;
    float xio4 = 0.25f*xi;
    y20[i1y  ] = xi;
    y20[i1y+1] = xio2;
    y2p[i1y  ] = xio2;
    y2p[i1y+1] = xio4;
    for (i1x=1,i1y=2; i1x<=i1xmax; ++i1x,i1y+=2) {
      xi = x[i2x][i1x];
      xio2 = 0.50f*xi;
      xio4 = 0.25f*xi;
      y20[i1y-1] += xio2;
      y20[i1y  ]  = xi;
      y20[i1y+1]  = xio2;
      y2p[i1y-1] += xio4;
      y2p[i1y  ]  = xio2;
      y2p[i1y+1]  = xio4;
    }
    if (i1x<n1x && i1y<n1y) {
      xi = x[i2x][i1x];
      xio2 = 0.50f*xi;
      xio4 = 0.25f*xi;
      y20[i1y-1] += xio2;
      y20[i1y  ]  = xi;
      y2p[i1y-1] += xio4;
      y2p[i1y  ]  = xio2;
    }

    // i2x < n2x && i2y < n2y-1
    for (i2x=1,i2y=2; i2x<=i2xmax; ++i2x,i2y+=2) {
      i1x = i1y = 0;
      y2m = y20;
      y20 = y2p;
      y2p = y[i2y+1];
      xi = x[i2x][i1x];
      xio2 = 0.50f*xi;
      xio4 = 0.25f*xi;
      y2m[i1y  ] += xio2;
      y2m[i1y+1] += xio4;
      y20[i1y  ]  = xi;
      y20[i1y+1]  = xio2;
      y2p[i1y  ]  = xio2;
      y2p[i1y+1]  = xio4;
      for (i1x=1,i1y=2; i1x<=i1xmax; ++i1x,i1y+=2) {
        xi = x[i2x][i1x];
        xio2 = 0.50f*xi;
        xio4 = 0.25f*xi;
        y2m[i1y-1] += xio4;
        y2m[i1y  ] += xio2;
        y2m[i1y+1] += xio4;
        y20[i1y-1] += xio2;
        y20[i1y  ]  = xi;
        y20[i1y+1]  = xio2;
        y2p[i1y-1] += xio4;
        y2p[i1y  ]  = xio2;
        y2p[i1y+1]  = xio4;
      }
      if (i1x<n1x && i1y<n1y) {
        xi = x[i2x][i1x];
        xio2 = 0.50f*xi;
        xio4 = 0.25f*xi;
        y2m[i1y-1] += xio4;
        y2m[i1y  ] += xio2;
        y20[i1y-1] += xio2;
        y20[i1y  ]  = xi;
        y2p[i1y-1] += xio4;
        y2p[i1y  ]  = xio2;
      }
    }

    // i2x==n2x-1 && i2y==n2y-1
    if (i2x<n2x && i2y<n2y) {
      i1x = i1y = 0;
      y2m = y20;
      y20 = y2p;
      y2p = null;
      xi = x[i2x][i1x];
      xio2 = 0.50f*xi;
      xio4 = 0.25f*xi;
      y2m[i1y  ] += xio2;
      y2m[i1y+1] += xio4;
      y20[i1y  ]  = xi;
      y20[i1y+1]  = xio2;
      for (i1x=1,i1y=2; i1x<=i1xmax; ++i1x,i1y+=2) {
        xi = x[i2x][i1x];
        xio2 = 0.50f*xi;
        xio4 = 0.25f*xi;
        y2m[i1y-1] += xio4;
        y2m[i1y  ] += xio2;
        y2m[i1y+1] += xio4;
        y20[i1y-1] += xio2;
        y20[i1y  ]  = xi;
        y20[i1y+1]  = xio2;
      }
      if (i1x<n1x && i1y<n1y) {
        xi = x[i2x][i1x];
        xio2 = 0.50f*xi;
        xio4 = 0.25f*xi;
        y2m[i1y-1] += xio4;
        y2m[i1y  ] += xio2;
        y20[i1y-1] += xio2;
        y20[i1y  ]  = xi;
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

  /*
  private static Multigrid2.A33 makeOperator(float[][][] d) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][][] a = new float[n2][n1][9];
    float ss = 0.5f*_sigma*_sigma;
    //float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    //float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    float r = 0.5f;
    float s = 0.5f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float d11 = ss*d[0][i2][i1];
        float d12 = ss*d[1][i2][i1];
        float d22 = ss*d[2][i2][i1];
        float x00 = x[i2][i1];
        float x01 = (i1>0)?x[i2][i1-1]:0.0f;
        float x10 = (i2>0)?x[i2-1][i1]:0.0f;
        float x11 = (i2>0 && i1>0)?x[i2-1][i1-1]:0.0f;
        float x1 = r*(x00-x01)+s*(x10-x11);
        float x2 = r*(x00-x10)+s*(x01-x11);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        y[i2][i1] += r*y1+r*y2;
        if (i1>0) y[i2][i1-1] -= r*y1-s*y2;
        if (i2>0) y[i2-1][i1  ] += s*y1-r*y2;
        if (i2>0 && i1>0) y[i2-1][i1-1] -= s*y1+s*y2;
      }
    }
    t = 0.5*(d[i]+d[i-1])*(x[i]-x[i-1]);
    y[i  ]  = t;
    y[i-1] -= t;
    
    y[i] = 0.5*(d[i]+d[i-1])*(x[i]-x[i-1]) -
           0.5*(d[i+1]+d[i])*(x[i+1]-x[i]);

    am = -0.5*(d[i]-d[i-1])
    a0 =  0.5*(2.0*d[i]-d[i-1]-d[i+1])
    ap = -0.5*(d[i+1]-d[i])

  }
  */
} 
