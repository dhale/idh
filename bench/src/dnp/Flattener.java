/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.logging.Logger;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates and applies shifts to flatten features in 2D and 3D images.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.03
 */
public class Flattener {

  public Flattener() {
  }

  public Flattener(double sigma) {
    _sigma = (float)sigma;
  }

  public float[][] findShifts(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] u2 = new float[n2][n1]; // 2nd component of normal vectors u
    float[][] el = new float[n2][n1]; // linearities
    float[][] r = new float[n2][n1]; // right-hand side
    float[][] s = new float[n2][n1]; // the shifts
    makeNormals(_sigma,f,u2,el);
    //fill(1.0f,el);
    float eps = 0.1f;
    makeRhs(eps,u2,el,r);
    Operator2 a = new LhsOperator2(eps,u2,el);
    u2 = el = null;
    solve(a,r,s);
    cleanShifts(s);
    float[] t1 = rampfloat(0.0f,1.0f,n1);
    float[] t2 = new float[n1];
    float[][] t = s;
    InverseInterpolator ii = new InverseInterpolator(n1,n1);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        t[i2][i1] = s[i2][i1]+t1[i1];
        if (i1>0 && t[i2][i1]<t[i2][i1-1])
          System.out.println("bad time = "+t[i2][i1]+", i1,i2 ="+i1+","+i2);
      }
      ii.invert(t[i2],t2);
      for (int i1=0; i1<n1; ++i1)
        s[i2][i1] = t1[i1]-t2[i1];
    }
    return s;
  }
  private static void cleanShifts(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if (s[i2][i1]<=s[i2][i1-1]-1.0f) {
          System.out.println("bad shift = "+s[i2][i1]+", i1,i2 ="+i1+","+i2);
          s[i2][i1] = s[i2][i1-1]-0.99f;
        }
      }
    }
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

  private static Logger log = 
    Logger.getLogger(Flattener.class.getName());

  private static final boolean PARALLEL = true; // false for single-threaded

  private float _sigma = 8.0f; // smoothing half-width in 1st dimension
  private float _small = 0.01f; // stop iterations when residuals are small
  private int _niter = 2000; // number of iterations

  private static interface Operator2 {
    public void apply(float[][] x, float[][] y);
  }

  private static class LhsOperator2 implements Operator2 {
    LhsOperator2(float sigma, float[][] u2, float[][] el) {
      _sigma = sigma;
      _u2 = u2;
      _el = el;
    }
    public void apply(float[][] x, float[][] y) {
      applyLhs(_sigma,_u2,_el,x,y);
    }
    private float _sigma;
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

  private static void xmakeRhs(
    float eps, float[][] u2, float[][] el, float[][] y) 
  {
    int n1 = y[0].length;
    int n2 = y.length;
    szero(y);
    float d11Constant = eps;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float p2i;
        if (u2i<0.0f) {
          p2i = (u2i>-10.0f*u1i)?u2i/u1i:-10.0f;
        } else {
          p2i = (u2i<10.0f*u1i)?u2i/u1i:10.0f;
        }
        float d11 = d11Constant;
        float d12 = 0.0f; // axis-aligned anisotropy
        float d22 = 1.0f;
        float x1 = 0.0f;
        float x2 = 0.5f*p2i;
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

  private static void xapplyLhs(
    float eps, float[][] u2, float[][] el, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    szero(y);
    float d11Constant = eps*eps;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float d11 = d11Constant;
        float d12 = 0.0f; // axis-aligned anisotropy
        float d22 = 1.0f;
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

  private static void makeRhs(
    float eps, float[][] u2, float[][] el, float[][] y) 
  {
    int n1 = y[0].length;
    int n2 = y.length;
    szero(y);
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float b11 = eps;
        float b12 = -u2i*eli;
        float b22 =  u1i*eli;
        float x1 = 0.0f;
        float x2 = 0.5f*u2i;
        float y1 = b11*x1+b12*x2;
        float y2 =        b22*x2;
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
    szero(y);
    eps = eps*eps;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float els = eli*eli;
        float u2i = u2[i2][i1];
        float u2s = u2i*u2i;
        float u1s = 1.0f-u2s;
        float u1i = sqrt(u1s);
        float d11 = eps+u2s*els;
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

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   */
  private void solve(Operator2 a, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] d = new float[n2][n1];
    float[][] q = new float[n2][n1];
    float[][] r = new float[n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    scopy(r,d);
    float delta = sdot(r,r);
    float deltaBegin = delta;
    float deltaSmall = sdot(b,b)*_small*_small;
    log.fine("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      log.finer("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
      a.apply(d,q);
      float dq = sdot(d,q);
      float alpha = delta/dq;
      saxpy( alpha,d,x);
      saxpy(-alpha,q,r);
      float deltaOld = delta;
      delta = sdot(r,r);
      float beta = delta/deltaOld;
      sxpay(beta,r,d);
    }
    log.fine("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }

  // Zeros array x.
  private static void szero(float[][] x) {
    zero(x);
  }
  private static void szero(float[][][] x) {
    if (PARALLEL) {
      szeroP(x);
    } else {
      szeroS(x);
    }
  }
  private static void szeroS(float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      szero(x[i3]);
  }
  private static void szeroP(final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            szero(x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Copys array x to array y.
  private static void scopy(float[][] x, float[][] y) {
    copy(x,y);
  }
  private static void scopy(float[][][] x, float[][][] y) {
    if (PARALLEL) {
      scopyP(x,y);
    } else {
      scopyS(x,y);
    }
  }
  private static void scopyS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      scopy(x[i3],y[i3]);
  }
  private static void scopyP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            scopy(x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Returns the dot product x'y.
  private static float sdot(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float d = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        d += x2[i1]*y2[i1];
      }
    }
    return d;
  }
  private static float sdot(float[][][] x, float[][][] y) {
    if (PARALLEL) {
      return sdotP(x,y);
    } else {
      return sdotS(x,y);
    }
  }
  private static float sdotS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    float d = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      d += sdot(x[i3],y[i3]);
    return d;
  }
  private static float sdotP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicFloat ad = new AtomicFloat(0.0f);
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          float d = 0.0f;
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            d += sdot(x[i3],y[i3]);
          ad.getAndAdd(d);
        }
      });
    }
    Threads.startAndJoin(threads);
    return ad.get();
  }

  // Computes y = y + a*x.
  private static void saxpy(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        y2[i1] += a*x2[i1];
      }
    }
  }
  private static void saxpy(float a, float[][][] x, float[][][] y) {
    if (PARALLEL) {
      saxpyP(a,x,y);
    } else {
      saxpyS(a,x,y);
    }
  }
  private static void saxpyS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpy(a,x[i3],y[i3]);
  }
  private static void saxpyP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            saxpy(a,x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Computes y = x + a*y.
  private static void sxpay(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        y2[i1] = a*y2[i1]+x2[i1];
      }
    }
  }
  private static void sxpay(float a, float[][][] x, float[][][] y) {
    if (PARALLEL) {
      sxpayP(a,x,y);
    } else {
      sxpayS(a,x,y);
    }
  }
  private static void sxpayS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sxpay(a,x[i3],y[i3]);
  }
  private static void sxpayP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            sxpay(a,x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
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
