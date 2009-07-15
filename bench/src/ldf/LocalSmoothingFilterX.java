/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.Tensors2;
import edu.mines.jtk.dsp.Tensors3;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Local smoothing of images with tensor filter coefficients.
 * Smoothing is performed by solving a sparse symmetric positive-definite
 * system of equations: (S'S+D'TD)y = S'Sx, where D is a matrix of 
 * derivative (difference) operators, S is a matrix of smoothing (sum)
 * operators, T is a matrix of tensor filter coefficients, x is an input 
 * image and y is an output image.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.10.13
 */
public class LocalSmoothingFilterX {

  /**
   * Constructs a local smoothing filter.
   * @param scale scale factor for all tensor coefficients.
   * @param small stop when L2 norm of residuals decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalSmoothingFilterX(double scale, double small, int niter) {
    _scale = (float)scale;
    _small = (float)small;
    _niter = niter;
  }

  /**
   * Applies this filter for specified tensor coefficients.
   * @param t tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(Tensors2 t, float[][] x, float[][] y) {
    Operator2 a = new LhsOperator2(_scale,t);
    float[][] r = SMOOTH?applyRhs(x):x;
    solve(a,r,y);
  }

  /**
   * Applies this filter for specified tensor coefficients.
   * @param t tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(Tensors3 t, float[][][] x, float[][][] y) {
    Operator3 a = new LhsOperator3(_scale,t);
    float[][][] r = SMOOTH?applyRhs(x):x;
    solve(a,r,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean PARALLEL = true;
  private static final boolean SMOOTH = true;

  private float _scale; // scale factor for tensor coefficients
  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator2 {
    public void apply(float[][] x, float[][] y);
  }
  private static interface Operator3 {
    public void apply(float[][][] x, float[][][] y);
  }

  private static class LhsOperator2 implements Operator2 {
    LhsOperator2(float s, Tensors2 t) {
      _s = s;
      _t = t;
    }
    public void apply(float[][] x, float[][] y) {
      applyLhs(_s,_t,x,y);
    }
    private float _s;
    private Tensors2 _t;
  }

  private static class LhsOperator3 implements Operator3 {
    LhsOperator3(float s, Tensors3 t) {
      _s = s;
      _t = t;
    }
    public void apply(float[][][] x, float[][][] y) {
      applyLhs(_s,_t,x,y);
    }
    private float _s;
    private Tensors3 _t;
  }

  /**
   * Returns y = S'Sx.
   */
  private static float[][] applyRhs(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        //         0.0625 = 1/16
        float xs = 0.0625f*(x00+x01+x10+x11);
        y[i2  ][i1  ] += xs;
        y[i2  ][i1-1] += xs;
        y[i2-1][i1  ] += xs;
        y[i2-1][i1-1] += xs;
      }
    }
    return y;
  }

  /**
   * Returns y = S'Sx.
   */
  private static float[][][] applyRhs(float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] y = new float[n3][n2][n1];
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        float[] x00 = x[i3  ][i2  ];
        float[] x01 = x[i3  ][i2-1];
        float[] x10 = x[i3-1][i2  ];
        float[] x11 = x[i3-1][i2-1];
        float[] y00 = y[i3  ][i2  ];
        float[] y01 = y[i3  ][i2-1];
        float[] y10 = y[i3-1][i2  ];
        float[] y11 = y[i3-1][i2-1];
        for (int i1=1; i1<n1; ++i1) {
          int i1m = i1-1;
          float x000 = x00[i1 ];
          float x001 = x00[i1m];
          float x010 = x01[i1 ];
          float x011 = x01[i1m];
          float x100 = x10[i1 ];
          float x101 = x10[i1m];
          float x110 = x11[i1 ];
          float x111 = x11[i1m];
          //         0.015625 = 1/64
          float xs = 0.015625f*(x000+x001+x010+x011+x100+x101+x110+x111);
          y00[i1 ] += xs;
          y00[i1m] += xs;
          y01[i1 ] += xs;
          y01[i1m] += xs;
          y10[i1 ] += xs;
          y10[i1m] += xs;
          y11[i1 ] += xs;
          y11[i1m] += xs;
        }
      }
    }
    return y;
  }

  /**
   * Computes y = (S'S+D'TD)x. Arrays x and y must be distinct.
   */
  private static void applyLhs(
    float s, Tensors2 t, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    if (SMOOTH) {
      szero(y);
    } else {
      scopy(x,y);
    }
    float[] ti = new float[3];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        t.getTensor(i1,i2,ti);
        float t11 = ti[0]*s;
        float t12 = ti[1]*s;
        float t22 = ti[2]*s;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.2500f*(xa-xb);
        float x2 = 0.2500f*(xa+xb);
        float y1 = t11*x1+t12*x2;
        float y2 = t12*x1+t22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        if (SMOOTH) {
          float xs = 0.0625f*(x00+x01+x10+x11);
          y[i2  ][i1  ] += ya+xs;
          y[i2  ][i1-1] -= yb-xs;
          y[i2-1][i1  ] += yb+xs;
          y[i2-1][i1-1] -= ya-xs;
        } else {
          y[i2  ][i1  ] += ya;
          y[i2  ][i1-1] -= yb;
          y[i2-1][i1  ] += yb;
          y[i2-1][i1-1] -= ya;
        }
      }
    }
  }

  private static void applyLhs(
    float s, Tensors3 t, float[][][] x, float[][][] y) 
  {
    if (SMOOTH) {
      szero(y);
    } else {
      scopy(x,y);
    }
    if (PARALLEL) {
      applyLhsParallel(s,t,x,y);
    } else {
      applyLhsSerial(s,t,x,y);
    }
  }

  private static void applyLhsSerial(
    float s, Tensors3 t, float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3)
      applyLhsSlice3(i3,s,t,x,y);
  }

  private static void applyLhsParallel(
    final float s, final Tensors3 t, final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;

    // i3 = 1, 3, 5, ...
    final AtomicInteger a1 = new AtomicInteger(1);
    Thread[] thread1 = Threads.makeArray();
    for (int ithread=0; ithread<thread1.length; ++ithread) {
      thread1[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a1.getAndAdd(2); i3<n3; i3=a1.getAndAdd(2))
            applyLhsSlice3(i3,s,t,x,y);
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
            applyLhsSlice3(i3,s,t,x,y);
        }
      });
    }
    Threads.startAndJoin(thread2);
  }


  /**
   * Computes y = (S'S+D'TD)x for one constant-i3 slice.
   */
  private static void applyLhsSlice3(
    int i3, float s, Tensors3 t, float[][][] x, float[][][] y) 
  {
    float[] ti = new float[6];
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
        t.getTensor(i1,i2,i3,ti);
        float t11 = ti[0]*s;
        float t12 = ti[1]*s;
        float t13 = ti[2]*s;
        float t22 = ti[3]*s;
        float t23 = ti[4]*s;
        float t33 = ti[5]*s;
        applyLhs(i1,t11,t12,t13,t22,t23,t33,x00,x01,x10,x11,y00,y01,y10,y11);
      }
    }
  }

  /**
   * Computes y = (S'S+D'TD)x for one sample.
   */
  private static void applyLhs(int i1,
   float t11, float t12, float t13, float t22, float t23, float t33,
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
    float x1 = 0.062500f*(xa-xb+xc+xd);
    float x2 = 0.062500f*(xa+xb-xc+xd);
    float x3 = 0.062500f*(xa+xb+xc-xd);
    float y1 = t11*x1+t12*x2+t13*x3;
    float y2 = t12*x1+t22*x2+t23*x3;
    float y3 = t13*x1+t23*x2+t33*x3;
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
    trace("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
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
    trace("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }
  private void solve(Operator3 a, float[][][] b, float[][][] x) {
    int n1 = b[0][0].length;
    int n2 = b[0].length;
    int n3 = b.length;
    float[][][] d = new float[n3][n2][n1];
    float[][][] q = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    scopy(r,d);
    float delta = sdot(r,r);
    float deltaBegin = delta;
    float deltaSmall = sdot(b,b)*_small*_small;
    trace("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      //trace("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
      //trace("  r min="+min(r)+" max="+max(r));
      //trace("  x min="+min(x)+" max="+max(x));
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
    trace("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
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
} 
