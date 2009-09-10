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
 * Conjugate-gradient solver for symmetric positive-semidefinite problems.
 * The problem solved is a linear system of equations Ax = b, where
 * multiplication by the matrix A is specified by an abstract operator.
 * Elements of the matrix A are not required by this solver. The solution 
 * column vector x and the right-hand-column vector b can be stored in 
 * 1D, 2D, or 3D arrays of floats.
 *
 * The operator A must be symmetric positive-semidefinite (SPSD).
 * Because the conjugate-gradient algorithm computes the solution 
 * vector x as a linear combination of products Ar, AAr, AAAr, etc., 
 * that solution vector x will include no components from the null 
 * space of A (except those due to rounding errors).
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.09
 */
public class CgLinearSolver {

  /**
   * An abstract linear operator A for 1D arrays.
   */
  public interface A1 {

    /**
     * Computes the matrix-vector product y = Ax.
     * @param x the input vector x.
     * @param y the output vector y.
     */
    public void apply(float[] x, float[] y);
  }

  /**
   * An abstract linear operator A for 2D arrays.
   */
  public interface A2 {

    /**
     * Computes the matrix-vector product y = Ax.
     * @param x the input vector x.
     * @param y the output vector y.
     */
    public void apply(float[][] x, float[][] y);
  }

  /**
   * An abstract linear operator A for 3D arrays.
   */
  public interface A3 {

    /**
     * Computes the matrix-vector product y = Ax.
     * @param x the input vector x.
     * @param y the output vector y.
     */
    public void apply(float[][][] x, float[][][] y);
  }

  /**
   * Constructs a solver with specified parameters.
   * Iterations end after either the specified number of iterations
   * or after the L2 norm of the residuals r = b-Ax is less than the 
   * specified small fraction of the L2 norm of the initial residuals.
   * @param niter the maximum number of iterations.
   * @param small the small fraction used to define convergence.
   */
  public CgLinearSolver(int niter, double small) {
    _niter = niter;
    _small = (float)small;
    setParallel(1);
  }

  /**
   * Sets the number of parallel threads used by this solver.
   * The number of threads used equals the specified factor times the 
   * number of available processors, or 1, whichever is greater. The
   * default factor is 1, so that the number of available processors
   * is used.
   * @param factor the factor.
   */
  public void setParallel(double factor) {
    int nprocessors = Threads.getAvailableProcessors();
    _nthread = Math.max(1,(int)factor*nprocessors);
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   * @param a the linear operator that represents the matrix A.
   * @param b the right-hand-side column vector.
   * @param x the solution column vector.
   */
  public void solve(A1 a, float[] b, float[] x) {
    int n1 = b.length;
    float[] d = new float[n1];
    float[] q = new float[n1];
    float[] r = new float[n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    scopy(r,d);
    float delta = sdot(r,r);
    float deltaBegin = delta;
    float deltaSmall = sdot(b,b)*_small*_small;
    int iter;
    logInit(delta);
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      logIter(iter,delta,deltaBegin);
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
    logDone(iter,delta,deltaBegin);
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   * @param a the linear operator that represents the matrix A.
   * @param b the right-hand-side column vector.
   * @param x the solution column vector.
   */
  public void solve(A2 a, float[][] b, float[][] x) {
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
    int iter;
    logInit(delta);
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      logIter(iter,delta,deltaBegin);
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
    logDone(iter,delta,deltaBegin);
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   * @param a the linear operator that represents the matrix A.
   * @param b the right-hand-side column vector.
   * @param x the solution column vector.
   */
  public void solve(A3 a, float[][][] b, float[][][] x) {
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
    int iter;
    logInit(delta);
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      logIter(iter,delta,deltaBegin);
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
    logDone(iter,delta,deltaBegin);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _niter; // maximum number of CG iterations
  private float _small; // stop iterations when residuals r = b-Ax are small
  private float _nthread; // number of threads to use for parallel solver

  // Logging.
  private static Logger _log = 
    Logger.getLogger(CgLinearSolver.class.getName());
  private static void logInit(float delta) {
    _log.fine("solve: delta="+delta);
  }
  private static void logIter(int iter, float delta, float deltaBegin) {
    _log.finer("solve: iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }
  private static void logDone(int iter, float delta, float deltaBegin) {
    _log.fine("solve: iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }

  // Zeros array x.
  private void szero(float[] x) {
    zero(x);
  }
  private void szero(float[][] x) {
    zero(x);
  }
  private void szero(float[][][] x) {
    if (_nthread>1) {
      szeroP(x);
    } else {
      szeroS(x);
    }
  }
  private void szeroS(float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      szero(x[i3]);
  }
  private void szeroP(final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
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
  private void scopy(float[] x, float[] y) {
    copy(x,y);
  }
  private void scopy(float[][] x, float[][] y) {
    copy(x,y);
  }
  private void scopy(float[][][] x, float[][][] y) {
    if (_nthread>1) {
      scopyP(x,y);
    } else {
      scopyS(x,y);
    }
  }
  private void scopyS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      scopy(x[i3],y[i3]);
  }
  private void scopyP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
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
  private float sdot(float[] x, float[] y) {
    int n1 = x.length;
    float d = 0.0f;
    for (int i1=0; i1<n1; ++i1)
      d += x[i1]*y[i1];
    return d;
  }
  private float sdot(float[][] x, float[][] y) {
    int n2 = x.length;
    float d = 0.0f;
    for (int i2=0; i2<n2; ++i2)
      d += sdot(x[i2],y[i2]);
    return d;
  }
  private float sdot(float[][][] x, float[][][] y) {
    if (_nthread>1) {
      return sdotP(x,y);
    } else {
      return sdotS(x,y);
    }
  }
  private float sdotS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    float d = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      d += sdot(x[i3],y[i3]);
    return d;
  }
  private float sdotP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicFloat ad = new AtomicFloat(0.0f);
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
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
  private void saxpy(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += a*x[i1];
  }
  private void saxpy(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      saxpy(a,x[i2],y[i2]);
  }
  private void saxpy(float a, float[][][] x, float[][][] y) {
    if (_nthread>1) {
      saxpyP(a,x,y);
    } else {
      saxpyS(a,x,y);
    }
  }
  private void saxpyS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpy(a,x[i3],y[i3]);
  }
  private void saxpyP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
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
  private void sxpay(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = a*y[i1]+x[i1];
  }
  private void sxpay(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      sxpay(a,x[i2],y[i2]);
  }
  private void sxpay(float a, float[][][] x, float[][][] y) {
    if (_nthread>1) {
      sxpayP(a,x,y);
    } else {
      sxpayS(a,x,y);
    }
  }
  private void sxpayS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sxpay(a,x[i3],y[i3]);
  }
  private void sxpayP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
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
}
