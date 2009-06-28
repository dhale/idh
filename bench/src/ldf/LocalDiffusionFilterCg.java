/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic diffusion filter via conjugate gradient iterations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.05
 */
public class LocalDiffusionFilterCg extends LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   * @param small stop when L2 norm of residuals decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalDiffusionFilterCg(double sigma, double small, int niter) {
    super(sigma);
    _dlf = new DirectionalLaplacianFilter(sigma);
    _sigma = (float)sigma;
    _small = (float)small;
    _niter = niter;
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void solveLinear(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveLinearSimple(ds,v1,x,y);
    //solveLinearSsor(ds,v1,x,y);
  }

  protected void solveLinear(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    solveLinearSimple(is,iw,x,y);
  }

  protected void solvePlanar(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    solvePlanarSimple(is,iu,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean PARALLEL = true;

  private float _sigma; // maximum filter half-width
  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations
  private DirectionalLaplacianFilter _dlf;

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator2 {
    public void apply(float[][] x, float[][] y);
  }
  private static interface Operator3 {
    public void apply(float[][][] x, float[][][] y);
  }

  private static class LinearOperator2 implements Operator2 {
    LinearOperator2(
      DirectionalLaplacianFilter dlf, float[][] ds, float[][] v1) 
    {
      _dlf = dlf;
      _ds = ds;
      _v1 = v1;
    }
    public void apply(float[][] x, float[][] y) {
      scopy(x,y);
      _dlf.applyLinear(_ds,_v1,x,y);
    }
    private float[][] _ds,_v1;
    private DirectionalLaplacianFilter _dlf;
  }
  private static class LinearOperator3 implements Operator3 {
    LinearOperator3(
      DirectionalLaplacianFilter dlf, byte[][][] is, short[][][] iw) 
    {
      _dlf = dlf;
      _is = is;
      _iw = iw;
    }
    public void apply(float[][][] x, float[][][] y) {
      scopy(x,y);
      _dlf.applyLinear(_is,_iw,x,y);
    }
    private byte[][][] _is;
    private short[][][] _iw;
    private DirectionalLaplacianFilter _dlf;
  }
  private static class PlanarOperator3 implements Operator3 {
    PlanarOperator3(
      DirectionalLaplacianFilter dlf, byte[][][] is, short[][][] iu) 
    {
      _dlf = dlf;
      _is = is;
      _iu = iu;
    }
    public void apply(float[][][] x, float[][][] y) {
      scopy(x,y);
      _dlf.applyPlanar(_is,_iu,x,y);
    }
    private byte[][][] _is;
    private short[][][] _iu;
    private DirectionalLaplacianFilter _dlf;
  }

  private static class LinearSsorOperator2 implements Operator2 {
    LinearSsorOperator2(
      DirectionalLaplacianFilter dlf, float[][] ds, float[][] v1) 
    {
      _s33 = dlf.makeLinearStencil33(ds,v1);
    }
    public void apply(float[][] x, float[][] y) {
      // y = (2-w) * inv(D/w+L') * D/w * inv(D/w+L) * x
      int n1 = x[0].length;
      int n2 = x.length;
      int n1m = n1-1;
      int n2m = n2-1;
      float w1 = 1.02f;    // optimal choice seems to be near 1.00, so
      float w2 = 2.00f-w1; // could possibly eliminate the w1, w2, ...
      float ow1 = 1.0f/w1;
      float[][] a = new float[n1][9];
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2m);
        float[] x0 = x[i2 ];
        float[] ym = y[i2m];
        float[] y0 = y[i2 ];
        float[] yp = y[i2p];
        _s33.get(i2,a);
        for (int i1=0; i1<n1; ++i1) {
          int i1m = max(i1-1,0);
          int i1p = min(i1+1,n1m);
          float[] ai = a[i1];
          float dow = (1.0f+ai[4])*ow1;
          y0[i1] = (x0[i1]-(ai[0]*ym[i1m] +
                            ai[1]*ym[i1 ] +
                            ai[2]*ym[i1p] +
                            ai[3]*y0[i1m]))/dow;
        }
      }
      for (int i2=n2-1; i2>=0; --i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2m);
        float[] x0 = x[i2 ];
        float[] ym = y[i2m];
        float[] y0 = y[i2 ];
        float[] yp = y[i2p];
        _s33.get(i2,a);
        for (int i1=n1-1; i1>=0; --i1) {
          int i1m = max(i1-1,0);
          int i1p = min(i1+1,n1m);
          float[] ai = a[i1];
          float d = (1.0f+ai[4])*ow1;
          y0[i1] *= w2;
          y0[i1] =  y0[i1]-(ai[5]*y0[i1p] +
                            ai[6]*yp[i1m] +
                            ai[7]*yp[i1 ] +
                            ai[8]*yp[i1p])/d;
        }
      }
    }
    private DirectionalLaplacianFilter.Stencil33 _s33;
  }

  private void solveLinearSimple(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    Operator2 a = new LinearOperator2(_dlf,ds,v1);
    solve(a,x,y);
  }
  private void solveLinearSsor(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    Operator2 a = new LinearOperator2(_dlf,ds,v1);
    Operator2 m = new LinearSsorOperator2(_dlf,ds,v1);
    solve(a,m,x,y);
  }
  private void solveLinearSimple(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    Operator3 a = new LinearOperator3(_dlf,is,iw);
    solve(a,x,y);
  }
  private void solvePlanarSimple(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    Operator3 a = new PlanarOperator3(_dlf,is,iu);
    solve(a,x,y);
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

  /**
   * Solves Ax = b via conjugate gradient iterations with preconditioner M.
   * Uses the initial values of x; does not assume they are zero.
   */
  private void solve(Operator2 a, Operator2 m, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] d = new float[n2][n1];
    float[][] q = new float[n2][n1];
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    m.apply(r,d);
    m.apply(b,s);
    float delta = sdot(r,d);
    float deltaBegin = delta;
    float deltaSmall = sdot(s,s)*_small*_small;
    trace("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      a.apply(d,q);
      float dq = sdot(d,q);
      float alpha = delta/dq;
      saxpy( alpha,d,x);
      saxpy(-alpha,q,r);
      m.apply(r,s);
      float deltaOld = delta;
      delta = sdot(r,s);
      float beta = delta/deltaOld;
      sxpay(beta,s,d);
    }
    trace("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }
  private void solve(Operator3 a, Operator3 m, float[][][] b, float[][][] x) {
    int n1 = b[0][0].length;
    int n2 = b[0].length;
    int n3 = b.length;
    float[][][] d = new float[n3][n2][n1];
    float[][][] q = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    float[][][] s = new float[n3][n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    m.apply(r,d);
    m.apply(b,s);
    float delta = sdot(r,d);
    float deltaBegin = delta;
    float deltaSmall = sdot(s,s)*_small*_small;
    trace("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      a.apply(d,q);
      float dq = sdot(d,q);
      float alpha = delta/dq;
      saxpy( alpha,d,x);
      saxpy(-alpha,q,r);
      m.apply(r,s);
      float deltaOld = delta;
      delta = sdot(r,s);
      float beta = delta/deltaOld;
      sxpay(beta,s,d);
    }
    trace("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }

  private static void scopy(float[][] x, float[][] y) {
    ArrayMath.copy(x,y);
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

  // Computes y = y + ax.
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

  // Computes y = x + ay.
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
  // testing

  public static void main(String[] args) {
    testOperators();
  }

  private static void testOperators() {
    int n1 = 100;
    int n2 = 101;
    DirectionalLaplacianFilter dlf = new DirectionalLaplacianFilter(1.0);
    float[][] ds = ArrayMath.randfloat(n1,n2);
    float[][] v1 = ArrayMath.randfloat(n1,n2);
    Operator2 a = new LinearOperator2(dlf,ds,v1);
    Operator2 m = new LinearSsorOperator2(dlf,ds,v1);
    testSpd(n1,n2,a);
    testSpd(n1,n2,m);
  }

  private static void testSpd(int n1, int n2, Operator2 a) {
    float[][] x = ArrayMath.sub(ArrayMath.randfloat(n1,n2),0.5f);
    float[][] y = ArrayMath.sub(ArrayMath.randfloat(n1,n2),0.5f);
    float[][] ax = ArrayMath.zerofloat(n1,n2);
    float[][] ay = ArrayMath.zerofloat(n1,n2);
    a.apply(x,ax);
    a.apply(y,ay);
    float xax = sdot(x,ax);
    float yay = sdot(y,ay);
    float yax = sdot(y,ax);
    float xay = sdot(x,ay);
    System.out.println("xax="+xax+" yay="+yay+" (should be positive)");
    System.out.println("yax="+yax+" xay="+xay+" (should be equal)");
  }
} 
