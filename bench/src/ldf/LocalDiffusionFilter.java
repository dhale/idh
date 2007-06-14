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
    _sigma = (float)sigma;
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

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float CG_SMALL = 0.000001f;

  private float _sigma;


  // Test code for multigrid.
  public static void main(String[] args) {
    testMg(1,1001);
    testMg(4,1001);
  }
  private static void testMg(int nlevel, int n1) {
    float[] d = new float[n1];
    float[] x = new float[n1];
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      d[i1] = 1.0f+64.0f*(float)i1/(float)(n1-1);
      if (i1%(n1/10)==1)
        x[i1] = 1.0f;
    }
    LocalDiffusionFilter ldf = new LocalDiffusionFilter(1.0);
    ldf.applyMg(nlevel,d,x,y);
    //edu.mines.jtk.mosaic.SimplePlot.asSequence(x);
    edu.mines.jtk.mosaic.SimplePlot.asSequence(y);
  }
  private void applyMg(int nlevel, float[] d, float[] x, float[] y) {
    float[][] dd = makeGaussianPyramid(nlevel,d);
    float[][] xx = makeGaussianPyramid(nlevel,x);
    float[] yi = new float[xx[nlevel-1].length];
    for (int ilevel=nlevel-1; ilevel>=0; --ilevel) {
      float[] di = dd[ilevel];
      float[] xi = xx[ilevel];
      if (ilevel>0)
        solveCg(di,xi,yi);
      if (ilevel>0) {
        int m1 = xx[ilevel-1].length;
        float[] yt = new float[m1];
        upsample(yi,yt);
        yi = yt;
      } else {
        Array.copy(yi,y);
      }
    }
  }

  private float[][] makeGaussianPyramid(int nlevel, float[] x) {
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    int n1 = x.length;
    float[][] y = new float[nlevel][];
    y[0] = x;
    for (int ilevel=1; ilevel<nlevel; ++ilevel) {
      float[] t = new float[n1];
      rgf.apply0(y[ilevel-1],t);
      n1 = (n1+1)/2;
      y[ilevel] = Array.copy(n1,0,2,t);
    }
    return y;
  }

  // Solve via conjugate gradient iterations.
  private void solveCg(float[] d, float[] x, float[] y) {
    int n1 = x.length;
    float[] r = new float[n1];
    float[] s = new float[n1];
    float[] t = new float[n1];
    Array.copy(x,r);
    applyForward(d,y,t);
    saxpy(-1.0f,t,r);
    Array.copy(r,s);
    float rr = dot(r,r);
    float stop = rr*CG_SMALL;
    trace("solveCg: n1="+n1+" stop="+stop+" rr="+rr);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForward(d,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      float rrold = rr;
      rr = dot(r,r);
      float beta = rr/rrold;
      for (int i1=0; i1<n1; ++i1)
        s[i1] = r[i1]+beta*s[i1];
    }
    trace("  niter="+niter+" rr="+rr);
  }

  // Apply diffusion operator.
  private void applyForward(float[] d, float[] x, float[] y) {
    int n1 = x.length;
    float ss = 0.5f*_sigma*_sigma;
    for (int i1=0; i1<n1; ++i1) {
      float d11 = ss*d[i1];
      float xi0 = x[i1];
      float xi1 = (i1>0)?x[i1-1]:0.0f;
      float x1 = xi0-xi1;
      float y1 = d11*x1;
      y[i1] = xi0+y1;
      if (i1>0) y[i1-1] -= y1;
    }
  }

  // Downsample from n1 x samples to (n1+1)/2 y samples.
  private static void downsample(float[] x, float[] y) {
    Array.zero(y);
    int n1 = x.length;
    int i1 = 0;
    int j1 = 0;
    float s = 1.0f/4.0f;
    float t = x[i1];
    y[i1] = s*t;
    for (i1=1,j1=0; i1<n1; ++i1,j1=i1/2) {
      t = x[i1  ] +
          x[i1-1];
      y[j1] += s*t;
    }
  }

  // Downsample from (n1,n2) x samples to ((n1+1)/2,(n2+1)/2) y samples.
  private static void downsample(float[][] x, float[][] y) {
    Array.zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float s = 1.0f/16.0f;
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

  // Upsample from (n1+1)/2 x samples to n1 y samples.
  private static void upsample(float[] x, float[] y) {
    int n1 = y.length;
    int i1 = 0;
    int j1 = 0;
    float s = 1.0f/2.0f;
    float t = s*x[i1];
    y[j1]  = t;
    for (j1=1,i1=0; j1<n1; ++j1,i1=j1/2) {
      t = s*x[i1];
      y[j1  ]  = t;
      y[j1-1] += t;
    }
  }

  // Upsample from ((n1+1)/2,(n2+1)/2) x samples to (n1,n2) y samples.
  private static void upsample(float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    float s = 1.0f/4.0f;
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



  // Test code ensures that upsampling is the transpose of downsampling,
  // to within a known scale factor. If in doubt, make main public and 
  // run these tests.
  private static void xmain(String[] args) {
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
    downsample(x,ax);
    upsample(y,ay);
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
    downsample(x,ax);
    upsample(y,ay);
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
}
