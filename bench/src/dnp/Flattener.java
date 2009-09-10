/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates and applies shifts to flatten features in 2D and 3D images.
 * In 2D, the shifts (in samples) are functions s(x1,x2) that flatten
 * image features. Specifically, for a 2D input image f(x1,x2), flattening 
 * is performed by computing g(x1,x2) = f(x1-s(x1,x2),x2). The shifts 
 * s(x1,x2) added to x1 in this expression are computed from slopes
 * measured in the input image f(x1,x2) so that features in the output
 * image g(x1,x2) vary as little as possible with x2.
 *
 * Flattens features in 3D images in the same way, but using shifts
 * s(x1,x2,x3) to compute g(x1,x2,x3) = f(x1-s(x1,x2,x3),x2,x3).
 * Features in the flattened image g(x1,x2,x3 vary little with x2 and x3.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.08
 */
public class Flattener {

  public Flattener(double sigma, double eps) {
    _sigma = (float)sigma;
    _eps = (float)eps;
  }

  public float[][] findShifts(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] u2 = new float[n2][n1]; // 2nd component of normal vectors u
    float[][] el = new float[n2][n1]; // linearities
    float[][] r = new float[n2][n1]; // right-hand side
    float[][] s = new float[n2][n1]; // the shifts
    makeNormals(_sigma,f,u2,el);
    makeRhs(_eps,u2,el,r);
    LhsOperator2 a = new LhsOperator2(_eps,u2,el);
    CgLinearSolver cls = new CgLinearSolver(_niter,_small);
    cls.solve(a,r,s);
    invertShifts(s);
    return s;
  }

  public float[][] applyShifts(float[][] f, float[][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      sub(r,s[i2],t);
      si.setUniformSamples(f[i2]);
      si.interpolate(n1,t,g[i2]);
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean PARALLEL = true; // false for single-threaded

  private float _sigma = 8.0f; // smoothing half-width to estimate slopes
  private float _eps = 0.1f; // penalty for ds/dx1 term in least-squares
  private float _small = 0.01f; // stop iterations when residuals are small
  private int _niter = 2000; // maximum number of iterations

  private static class LhsOperator2 
    implements CgLinearSolver.Operator2 {
    LhsOperator2(float eps, float[][] u2, float[][] el) {
      _eps = eps;
      _u2 = u2;
      _el = el;
    }
    public void apply(float[][] x, float[][] y) {
      applyLhs(_eps,_u2,_el,x,y);
    }
    private float _eps;
    private float[][] _u2;
    private float[][] _el;
  }
 
  private static void makeNormals(
    float sigma, float[][] f, float[][] u2, float[][] el) 
  {
    int n1 = f[0].length;
    int n2 = f.length;

    // Normal vectors and linearities.
    float[][] u1 = new float[n2][n1];
    float[][] us = u1;
    float sigma1 = sigma;
    float sigma2 = 1.0f;
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    lof.applyForNormalLinear(f,u1,u2,el);
    float usmax = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float usi = sqrt(u1i*u1i+u2i*u2i);
        if (usi>usmax)
          usmax = usi;
        us[i2][i1] = usi;
      }
    }

    // Normalize.
    float tiny = 0.001f;
    float utiny = usmax*tiny;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float uscale = 1.0f/(us[i2][i1]+utiny);
        u2[i2][i1] *= uscale;
      }
    }
    u1 = us = null;

    // Use tiny values in large areas where image is zero.
    ZeroMask zm = new ZeroMask(0.1,sigma,1,f);
    zm.apply(tiny,u2);
    zm.apply(tiny,el);
  }

  private static void invertShifts(float[][] s) {
    cleanShifts(s);
    int n1 = s[0].length;
    int n2 = s.length;
    float[] t1 = rampfloat(0.0f,1.0f,n1);
    float[] t2 = new float[n1];
    InverseInterpolator ii = new InverseInterpolator(n1,n1);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1)
        s[i2][i1] += t1[i1];
      ii.invert(s[i2],t2);
      float tmin = -5.0f;
      float tmax = n1-1+5.0f;
      for (int i1=0; i1<n1; ++i1) {
        if (t2[i1]<tmin) t2[i1] = tmin;
        if (t2[i1]>tmax) t2[i1] = tmax;
        s[i2][i1] = t1[i1]-t2[i1];
      }
    }
  }

  private static void cleanShifts(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if (s[i2][i1]<=s[i2][i1-1]-1.0f)
          s[i2][i1] = s[i2][i1-1]-0.99f;
      }
    }
  }

  private static void makeRhs(
    float eps, float[][] u2, float[][] el, float[][] y) 
  {
    int n1 = y[0].length;
    int n2 = y.length;
    zero(y);
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float b12 = -u2i*eli;
        float b22 =  u1i*eli;
        // float x1 = 0.0f;
        float x2 = 0.5f*u2i;
        float y1 = b12*x2;
        float y2 = b22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void applyLhs(
    float eps, float[][] u2, float[][] el, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    zero(y);
    float epss = eps*eps;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float els = eli*eli;
        float u2i = u2[i2][i1];
        float u2s = u2i*u2i;
        float u1s = 1.0f-u2s;
        float u1i = sqrt(u1s);
        float d11 = epss+u2s*els;
        float d12 = -u2i*u1i*els;
        float d22 = u1s*els;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.25f*(xa-xb);
        float x2 = 0.25f*(xa+xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // TODO
  public static void normals(
    float[][][] f, float[][][] u2, float[][][] u3, float[][][] ep) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] us = u1;
    LocalOrientFilter lof = new LocalOrientFilter(8,1,1);
    lof.applyForNormalPlanar(f,u1,u2,u3,ep);
    float usmax = 0.0f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float usi = sqrt(u1i*u1i+u2i*u2i+u3i*u3i);
          if (usi>usmax)
            usmax = usi;
          us[i3][i2][i1] = usi;
        }
      }
    }
    float utiny = 0.001f*usmax;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float uscale = 1.0f/(us[i3][i2][i1]+utiny);
          u2[i3][i2][i1] *= uscale;
          u3[i3][i2][i1] *= uscale;
        }
      }
    }
  }
} 

/*
  //lhs

  private static void applyLhs(
    float eps, float[][][] u2, float[][] el, float[][] x, float[][] y) 
    Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    if (SMOOTH) {
      szero(y);
    } else {
      scopy(x,y);
    }
    if (PARALLEL) {
      applyLhsParallel(d,c,s,x,y);
    } else {
      applyLhsSerial(d,c,s,x,y);
    }
  }

  private static void applyLhsSerial(
    Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3)
      applyLhsSlice3(i3,d,c,s,x,y);
  }

  private static void applyLhsParallel(
    final Tensors3 d, final float c, final float[][][] s, 
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;

    // i3 = 1, 3, 5, ...
    final AtomicInteger a1 = new AtomicInteger(1);
    Thread[] thread1 = Threads.makeArray();
    for (int ithread=0; ithread<thread1.length; ++ithread) {
      thread1[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a1.getAndAdd(2); i3<n3; i3=a1.getAndAdd(2))
            applyLhsSlice3(i3,d,c,s,x,y);
        }
      });
    }
    Threads.startAndJoin(thread1);

    // i3 = 2, 4, 6, ...
    final AtomicInteger a2 = new AtomicInteger(2);
    Thread[] thread2 = Threads.makeArray();
    for (int ithread=0; ithread<thread2.length; ++ithread) {
      thread2[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a2.getAndAdd(2); i3<n3; i3=a2.getAndAdd(2))
            applyLhsSlice3(i3,d,c,s,x,y);
        }
      });
    }
    Threads.startAndJoin(thread2);
  }
*/


  /**
   * Computes y = (S'S+D'TD)x for one constant-i3 slice.
   */
   /*
  private static void applyLhsSlice3(
    int i3, Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    float[] di = new float[6];
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        d.getTensor(i1,i2,i3,di);
        float csi = (s!=null)?c*s[i3][i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d13 = di[2]*csi;
        float d22 = di[3]*csi;
        float d23 = di[4]*csi;
        float d33 = di[5]*csi;
        applyLhs(i1,d11,d12,d13,d22,d23,d33,x00,x01,x10,x11,y00,y01,y10,y11);
      }
    }
  }
   */

  /**
   * Computes y = (S'S+D'TD)x for one sample.
   */
   /*
  private static void applyLhs(int i1,
   float d11, float d12, float d13, float d22, float d23, float d33,
   float[] x00, float[] x01, float[] x10, float[] x11,
   float[] y00, float[] y01, float[] y10, float[] y11)
  {
    int i1m = i1-1;
    float x000 = x00[i1 ];
    float x001 = x00[i1m];
    float x010 = x01[i1 ];
    float x011 = x01[i1m];
    float x100 = x10[i1 ];
    float x101 = x10[i1m];
    float x110 = x11[i1 ];
    float x111 = x11[i1m];
    //float x1 = 0.0625f*(x000+x010+x100+x110-x001-x011-x101-x111);
    //float x2 = 0.0625f*(x000+x001+x100+x101-x010-x011-x110-x111);
    //float x3 = 0.0625f*(x000+x001+x010+x011-x100-x101-x110-x111);
    float xa = x000-x111;
    float xb = x001-x110;
    float xc = x010-x101;
    float xd = x100-x011;
    float x1 = 0.0625f*(xa-xb+xc+xd);
    float x2 = 0.0625f*(xa+xb-xc+xd);
    float x3 = 0.0625f*(xa+xb+xc-xd);
    float y1 = d11*x1+d12*x2+d13*x3;
    float y2 = d12*x1+d22*x2+d23*x3;
    float y3 = d13*x1+d23*x2+d33*x3;
    float ya = y1+y2+y3;
    float yb = y1-y2+y3;
    float yc = y1+y2-y3;
    float yd = y1-y2-y3;
    if (SMOOTH) {
      float xs = 0.015625f*(x000+x001+x010+x011+x100+x101+x110+x111);
      y00[i1 ] += ya+xs;
      y00[i1m] -= yd-xs;
      y01[i1 ] += yb+xs;
      y01[i1m] -= yc-xs;
      y10[i1 ] += yc+xs;
      y10[i1m] -= yb-xs;
      y11[i1 ] += yd+xs;
      y11[i1m] -= ya-xs;
    } else {
      y00[i1 ] += ya;
      y00[i1m] -= yd;
      y01[i1 ] += yb;
      y01[i1m] -= yc;
      y10[i1 ] += yc;
      y10[i1m] -= yb;
      y11[i1 ] += yd;
      y11[i1m] -= ya;
    }
  }
   */
