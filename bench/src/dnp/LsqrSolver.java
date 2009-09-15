/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.logging.Logger;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * An implementation of the LSQR algorithm of Paige and Saunders.
 * Let A denote a sparse matrix with m rows and n columns, b a column 
 * vector with m rows, and damp a scalar. This solver iteratively
 * computes a solution x to one of the following problems:
 * <pre><code>
 * (1) unsymmetric equations -- A*x = b
 * (2) linear least squares  -- A*x = b, in the least-squares sense
 * (3) damped least squares  -- Abar*x = bbar, in the least-squares sense,
 *     where   Abar = | A      |  and  bbar = |b| 
 *                    | damp*I |              |0|
 * </code></pre>
 * This solver does not require that the elements of the matrix A be 
 * specified explicitly. Rather, it requires only multiplication by A 
 * and its transpose A'. Specifically, when solving any of the three 
 * problems listed above, the user must provide an abstract linear 
 * operator A that updates vectors x and y via y = y+A*x and x = x+A'*y.
 * The column vector x has n rows; the column vector y has m rows.
 * <p>
 * The solver is iterative, and iterations typically end before any of the 
 * problems listed above are solved exactly. Iterations may, for example, 
 * end when errors in the residual r = A*x - b are less than estimated 
 * errors in the elements of the column vector b. When constructing a 
 * solver, the user must specify the following parameters:
 * <dl>
 * <dt>atol</dt>
 * <dd>estimate of relative error in elements of the matrix A.</dd>
 * <dt>btol</dt>
 * <dd>estimate of relative error in the right-hand-side vector b.</dd>
 * <dt>ctol</dt>
 * <dd>estimate of the reciprocal of the condition number of matrix Abar.</dd>
 * <dt>maxi</dt>
 * <dd>maximum number of iterations to perform.</dd>
 * </dl>
 * For example, if b is accurate to about six digits, then set btol = 1.0e-6.
 * <p>
 * The number of iterations required to achieve a certain accuracy depends
 * strongly on the scaling of the problem. For example, in problem (1) above
 * the solution is unaltered by row scaling. If a row of A is very small or
 * large compared to the other rows of A, then the small or large row of A
 * should be scaled up or down.
 * <p>
 * In problems 1 and 2, the solution x is easily recovered following column
 * scaling. Unless better information is available, and where feasible, the 
 * non-zero columns of A should be scaled so that they all have the same
 * Euclidean norm.
 * <p>
 * In problem 3, there is no freedom to scale if damp is nonzero. However,
 * the value of damp should be assigned only after scaling A appropriately.
 * The parameter damp exists to help regularize ill-conditioned systems, by
 * preventing the true solution from being very large.
 * <p>
 * Adapted from Paige, C.C, and Saunders, M.A., 1982, Algorithm 583,
 * Collected algorithms of the ACM. Any differences here are likely 
 * errors in programming.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.13
 */
public class LsqrSolver {

  /**
   * Iterations stop when one of these conditions is met.
   */
  public enum Stop {
    /**
     * x = 0 is the exact solution. No iterations were performed.
     */
    ZERO,
    /**
     * The equations A*x = b are probably compatible. The L2 norm ||b-A*x|| 
     * of residuals is sufficiently small, given the values of atol and btol.
     */
    RTOL,
    /**
     * The equations A*x = b are probably not compatible. The least-squares
     * solution x is sufficiently accurate, given the value of atol.
     */
    ATOL,
    /**
     * An estimate of 1/cond(Abar) is less than ctol. The system A*x = b
     * appears to be ill-conditioned. Or, the implementations of A and A'
     * may be inconsistent.
     */
    CTOL,
    /**
     * The equations A*x = b are probably compatible. Further iterations
     * likely cannot further reduce the L2 norm ||b-A*x|| of residuals.
     */
    RTOL_EPSILON,
    /**
     * The equations A*x = b are probably not compatible. Further iterations
     * likely cannot increase the accuracy of the least-squares solution.
     */
    ATOL_EPSILON,
    /**
     * An estimate of 1/cond(Abar) is too small for further iterations.
     * The implementations of A and A' may be inconsistent.
     */
    CTOL_EPSILON,
    /**
     * The maximum number of iterations was performed.
     */
    MAXI
  }

  /**
   * Information returned by this iterative solver.
   * Fields in this information are described below using the following
   * definitions and relationships:
   * <pre><code>
   *     Abar = | A      |  ,  bbar = |b| 
   *            | damp*I |            |0|
   *
   *        r = b - A*x,       rbar = bbar - Abar*x
   *
   *    norm(r) = the L2 norm ||r||
   *
   *    rnorm = norm(rbar) = sqrt(norm(r)^2 + damp^2 * norm(x)^2)
   * </code></pre>
   * This solver minimizes the function rnorm with respect to x.
   */
  public static class Info {
    private Info(
      Stop stop, int niter, 
      float anorm, float acond, float rnorm, float arnorm, float xnorm)
    {
      this.stop = stop;
      this.niter = niter;
      this.anorm = anorm;
      this.acond = acond;
      this.rnorm = rnorm;
      this.arnorm = arnorm;
      this.xnorm = xnorm;
    }
    /** 
     * The condition that caused iterations to stop. 
     */
    public Stop stop;
    /**
     * The number of iterations performed.
     */
    public int niter;
    /**
     * An estimate of the Frobenius norm of Abar.
     * This is the square-root of the sum of squares of the elements of Abar.
     * If damp is small and if the columns of A have all been scaled to have
     * length 1.0, then anorm should equal roughly sqrt(n). A radically
     * different value for anorm could indicate that implementations of A
     * and A' are inconsistent.
     */
    public float anorm;
    /**
     * An estimate of cond(Abar), the condition number of Abar. A very
     * high value may indicate that the implementations of A and A' are
     * inconsistent.
     */
    public float acond;
    /**
     * An estimate of the final value of norm(rbar). This is the function
     * that is being minimized. Its value will be small if A*x = b has a
     * solution.
     */
    public float rnorm;
    /**
     * An estimate of the final value of norm(Abar'*rbar). This is the 
     * norm of the residual for the usual normal equations. This value
     * should be small in all cases; it will often be smaller than the
     * true value computed from the solution vector x.
     */
    public float arnorm;
    /**
     * An estimate of the L2 norm of the final solution vector x.
     */
    public float xnorm;
  }

  /**
   * Abstract linear operators A and A' for 1D arrays.
   */
  public interface A1 {

    /**
     * Accumulates the matrix-vector product y += Ax.
     * @param x the input vector x[n].
     * @param y the output vector y[m].
     */
    public void apply(float[] x, float[] y);

    /**
     * Accumulates the matrix-vector product x += A'y.
     * @param x the input vector y[m].
     * @param y the output vector x[n].
     */
    public void applyTranspose(float[] y, float[] x);
  }

  /**
   * Abstract linear operators A and A' for 2D arrays.
   */
  public interface A2 {

    /**
     * Accumulates the matrix-vector product y += Ax.
     * @param x the input vector x.
     * @param y the output vector y.
     */
    public void apply(float[][] x, float[][] y);

    /**
     * Accumulates the matrix-vector product x += A'y.
     * @param x the input vector y.
     * @param y the output vector x.
     */
    public void applyTranspose(float[][] y, float[][] x);
  }

  /**
   * Abstract linear operators A and A' for 3D arrays.
   */
  public interface A3 {

    /**
     * Accumulates the matrix-vector product y += Ax.
     * @param x the input vector x.
     * @param y the output vector y.
     */
    public void apply(float[][][] x, float[][][] y);

    /**
     * Accumulates the matrix-vector product x += A'y.
     * @param x the input vector y.
     * @param y the output vector x.
     */
    public void applyTranspose(float[][][] y, float[][][] x);
  }

  /**
   * Constructs a solver with specified parameters.
   * @param atol estimate of relative error in the matrix A.
   * @param btol estimate of relative error in the right-hand-side vector b.
   * @param ctol estimate of reciprocal of condition number of matrix Abar.
   * @param maxi maximum number of iterations to perform.
   */
  public LsqrSolver(double atol, double btol, double ctol, int maxi) {
    _atol = (float)atol;
    _btol = (float)btol;
    _ctol = (float)ctol;
    _maxi = maxi;
    setParallel(1);
  }

  /**
   * Sets the number of parallel threads used by this solver.
   * The number of threads used equals the specified factor times the 
   * number of available processors, or 1, whichever is greater. The
   * default factor is 1, so that the number of threads equals the
   * number of available processors.
   * @param factor the factor.
   */
  public void setParallel(double factor) {
    int nprocessors = Threads.getAvailableProcessors();
    _nthread = Math.max(1,(int)factor*nprocessors);
  }

  /**
   * Solves the least-squares problem.
   * Uses the initial values of x; does not assume they are zero.
   * @param a the linear operator that represents the matrix A.
   * @param b array[m] for the right-hand-side column vector.
   * @param x array[n] for the solution column vector.
   */
  public Info solve(A1 a, float damp, float[] b, float[] x) {
    int m = b.length;
    int n = x.length;
    float zero = 0.0f;
    float one = 1.0f;

    // Initialize.
    Stop stop = Stop.ZERO;
    int niter = 0;
    int nstop = 0;
    float anorm = zero;
    float acond = zero;
    float bbnorm = zero;
    float dampsq = damp*damp;
    float ddnorm = zero;
    float res1 = zero;
    float res2 = zero;
    float rnorm = zero;
    float arnorm = zero;
    float xnorm = zero;
    float xxnorm = zero;
    float cs2 = -one;
    float sn2 = zero;
    float z = zero;
    float[] u = new float[m]; // or save space with u = b?
    float[] v = new float[n];
    float[] w = new float[n];
    scopy(b,u);
    szero(x);

    // Set up the first vectors u and v for the bidiagonalization.
    // These vectors satisfy beta*u = b and alfa*v = A'u.
    float alfa = zero;
    float beta = snorm(u);
    if (beta>zero) {
      sscal(one/beta,u);
      a.applyTranspose(u,v);
      alfa = snorm(v);
    }
    if (alfa>zero) {
      sscal(one/alfa,v);
      scopy(v,w);
    }
    arnorm = alfa*beta;

    // If trivial solution x = 0, simply return.
    if (arnorm==zero)
      return new Info(stop,niter,anorm,acond,rnorm,arnorm,xnorm);

    // Otherwise, continue set up.
    rnorm = beta;
    float bnorm = beta;
    float phibar = beta;
    float rhobar = alfa;

    // While the number of iterations is less than the maximum, ...
    while (niter<_maxi) {
      ++niter;

      // Perform the next step of the bidiagonalization to obtain
      // the next beta, u, alfa, and v. These satisfy the relations
      //   beta*u = A v - alfa*u 
      //   alfa*v = A'u - beta*v 
      sscal(-alfa,u);
      a.apply(v,u);
      beta = snorm(u);
      bbnorm += alfa*alfa+beta*beta+dampsq;
      if (beta>zero) {
        sscal(one/beta,u);
        sscal(-beta,v);
        a.applyTranspose(u,v);
        alfa = snorm(v);
        if (alfa>zero)
          sscal(one/alfa,v);
      }

      // Use a plane rotation to eliminate the damping parameter. This
      // alters the diagonal (rhobar) of the lower-bidiagonal matrix.
      float rhbar2 = rhobar*rhobar+dampsq;
      float rhbar1 = sqrt(rhbar2);
      float cs1 = rhobar/rhbar1;
      float sn1 = damp/rhbar1;
      float psi = sn1*phibar;
      phibar = cs1*phibar;

      // Use a plane rotation to eliminate the subdiagonal element (beta)
      // of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
      float rho = sqrt(rhbar2+beta*beta);
      float cs = rhbar1/rho;
      float sn = beta/rho;
      float theta = sn*alfa;
      float phi = cs*phibar;
      float tau = sn*phi;
      phibar = sn*phibar;
      rhobar = -cs*alfa;

      // Update x, w, and the standard error estimates.
      float t1 = phi/rho;
      float t2 = -theta/rho;
      float t3 = one/rho;
      for (int i=0; i<n; ++i) {
        float wi = w[i];
        float wt1 = wi*t1;
        float wt2 = wi*t2;
        float wt3 = wi*t3;
        x[i] += wt1;
        w[i] = wt2+v[i];
        ddnorm += wt3*wt3;
      }

      // Use a plane rotation to eliminate the super-diagonal element
      // (theta) of the upper-bidiagonal matrix. The use the result to
      // estimate norm(x).
      float delta = sn2*rho;
      float gambar = -cs2*rho;
      float rhs = phi-delta*z;
      float zbar = rhs/gambar;
      float gamma = sqrt(gambar*gambar+theta*theta);
      cs2 = gambar/gamma;
      sn2 = theta/gamma;
      z = rhs/gamma;
      xnorm = sqrt(xxnorm+zbar*zbar);
      xxnorm += z*z;

      // Test for convergence. First estimate the norm and condition
      // number of the matrix Abar, and the norms of rbar and Abar'.
      anorm = sqrt(bbnorm);
      acond = anorm*sqrt(ddnorm);
      res1 = phibar*phibar;
      res2 += psi*psi;
      rnorm = sqrt(res1+res2);
      arnorm = alfa*abs(tau);

      // Now use these norms to estimate certain other quantities,
      // some of which will be small near a solution.
      float test1 = rnorm/bnorm;
      float test2 = (rnorm>zero)?arnorm/(anorm*rnorm):zero;
      float test3 = one/acond;
      float rtol = _btol+_atol*anorm*xnorm/bnorm;

      // The following tests guard against extremely small values of
      // atol, btol, or ctol. (The user may have set any or all of
      // the parameters atol, btol, and ctol to zero.) The effect is
      // equivalent to the user specifying atol = btol = ctol = epsilon.
      t1 = one+test1/(one+anorm*xnorm/bnorm);
      t2 = one+test2;
      t3 = one+test3;
      if (niter>=_maxi) stop = Stop.MAXI;
      if (t3<=one) stop = Stop.CTOL_EPSILON;
      if (t2<=one) stop = Stop.ATOL_EPSILON;
      if (t1<=one) stop = Stop.RTOL_EPSILON;

      // Check tolerances set by the user.
      if (test3<=_ctol) stop = Stop.CTOL;
      if (test2<=_atol) stop = Stop.ATOL;
      if (test1<=rtol) stop = Stop.RTOL;

      // Stop if appropriate. Convergence criteria are required to
      // be met on nconv consecutive iterations, where nconv is a
      // constant defined here. Suggested value: nconv = 1, 2, or 3.
      final int nconv = 1;
      if (stop==Stop.ZERO) {
        nstop = 0;
      } else {
        ++nstop;
        if (nstop<nconv && niter<_maxi) {
          stop = Stop.ZERO;
        } else {
          break;
        }
      }
    }
    return new Info(stop,niter,anorm,acond,rnorm,arnorm,xnorm);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _atol; // relative error in A
  private float _btol; // relative error in B
  private float _ctol; // reciprocal of upper limit on cond number of Abar
  private int _maxi; // upper limit on number of iterations
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
  private static void szero(float[] x) {
    zero(x);
  }
  private static void szero(float[][] x) {
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
  private static void scopy(float[] x, float[] y) {
    copy(x,y);
  }
  private static void scopy(float[][] x, float[][] y) {
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

  // Scales array x.
  private static void sscal(float a, float[] x) {
    mul(a,x,x);
  }
  private static void sscal(float a, float[][] x) {
    sscal(a,x);
  }
  private void sscal(float a, float[][][] x) {
    if (_nthread>1) {
      sscalP(a,x);
    } else {
      sscalS(a,x);
    }
  }
  private void sscalS(float a, float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sscal(a,x[i3]);
  }
  private void sscalP(final float a, final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            sscal(a,x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Returns the sum of squares of elements of x.
  private static float ssums(float[] x) {
    int n1 = x.length;
    float sum = 0.0f;
    for (int i1=0; i1<n1; ++i1) {
      float xi = x[i1];
      sum += xi*xi;
    }
    return sum;
  }
  private static float ssums(float[][] x) {
    int n2 = x.length;
    float sum = 0.0f;
    for (int i2=0; i2<n2; ++i2)
      sum += ssums(x[i2]);
    return sum;
  }
  private float ssums(float[][][] x) {
    if (_nthread>1) {
      return ssumsP(x);
    } else {
      return ssumsS(x);
    }
  }
  private float ssumsS(float[][][] x) {
    int n3 = x.length;
    float sum = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      sum += ssums(x[i3]);
    return sum;
  }
  private float ssumsP(final float[][][] x) {
    final int n3 = x.length;
    final AtomicFloat as = new AtomicFloat(0.0f);
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          float sum = 0.0f;
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            sum += ssums(x[i3]);
          as.getAndAdd(sum);
        }
      });
    }
    Threads.startAndJoin(threads);
    return as.get();
  }

  // Returns L2 norm of x.
  private static float snorm(float[] x) {
    return sqrt(ssums(x));
  }
  private static float snorm(float[][] x) {
    return sqrt(ssums(x));
  }
  private float snorm(float[][][] x) {
    return sqrt(ssums(x));
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing
  
  private static class TestH implements A1 {
    public int m,n;
    public int ndup,npow;
    public float damp,acond,rnorm;
    public float[] x,b,d,hy,hz,w;
    public TestH(int m, int n, int ndup, int npow, float damp, float[] x) {
      this.m = m;
      this.n = n;
      this.ndup = ndup;
      this.npow = npow;
      this.damp = damp;
      this.x = x;
      float[] b = this.b = new float[m];
      float[] d = this.d = new float[n];
      float[] hy = this.hy = new float[m];
      float[] hz = this.hz = new float[n];
      float[] w = this.w = new float[m];

      // Make two vectors with norm 1 for Householder transformations.
      float dampsq = damp*damp;
      float fourpi = 4.0f*FLT_PI;
      float alfa = fourpi/m;
      float beta = fourpi/n;
      for (int i=0; i<m; ++i)
        hy[i] = sin((i+1)*alfa);
      for (int i=0; i<n; ++i)
        hz[i] = cos((i+1)*beta);
      alfa = snorm(hy);
      beta = snorm(hz);
      sscal(-1.0f/alfa,hy);
      sscal(-1.0f/beta,hz);

      // Set the diagonal matrix D, which contains the singular values of A.
      for (int i=0; i<n; ++i) {
        int j = (i+ndup)/ndup;
        float t = j*ndup;
        t /= n;
        d[i] = pow(t,npow);
      }
      acond = sqrt((d[n-1]*d[n-1]+dampsq)/(d[0]*d[0]+dampsq));

      // Compute the residual vector, storing it in b.
      hprod(hz,x,b);
      for (int i=0; i<n; ++i)
        b[i] *= dampsq/d[i];
      float t = 1.0f;
      for (int i=n; i<m; ++i) {
        int j = i-n;
        b[i] = (t*j)/m;
        t = -t;
      }
      hprod(hy,b,b);
      float bnorm = snorm(b);
      float xnorm = snorm(x);
      rnorm = sqrt(bnorm*bnorm+dampsq*xnorm*xnorm);

      // Now compute the true b = r + A*x.
      apply(x,b);
    }
    public void apply(float[] x, float[] y) {
      // y += A*x, where A = Hy*D*Hz
      hprod(hz,x,w);
      for (int i=0; i<n; ++i)
        w[i] *= d[i];
      for (int i=n; i<m; ++i)
        w[i] = 0.0f;
      hprod(hy,w,w);
      for (int i=0; i<m; ++i)
        y[i] += w[i];
    }
    public void applyTranspose(float[] y, float[] x) {
      // x += A'*y, where A = Hy*D*Hz
      hprod(hy,y,w);
      for (int i=0; i<n; ++i)
        w[i] *= d[i];
      hprod(hz,w,w);
      for (int i=0; i<n; ++i)
        x[i] += w[i];
    }
    private void hprod(float[] hz, float[] x, float[] y) {
      // Householder transformation: y = (I-2*Hz*Hz')*x
      int n = hz.length;
      float s = 0.0f;
      for (int i=0; i<n; ++i)
         s += hz[i]*x[i];
      s += s;
      for (int i=0; i<n; ++i)
        y[i] = x[i]-s*hz[i];
    }
  }
  public static void testH(int m, int n, int ndup, int npow, float damp) {
    float[] xtrue = new float[n];
    for (int j=0; j<n; ++j)
      xtrue[j] = n-j-1;
    TestH p = new TestH(m,n,ndup,npow,damp,xtrue);
    float acond = p.acond;
    float rnorm = p.rnorm;
    float[] b = p.b;
    float atol = 1.0e-6f;
    float btol = atol;
    float ctol = 0.1f/acond;
    int maxi = m+n+50;
    LsqrSolver ls = new LsqrSolver(atol,btol,ctol,maxi);
    float[] x = new float[n];
    LsqrSolver.Info info = ls.solve(p,damp,b,x);
    println("x="); dump(x);
    println("xtrue="); dump(x);
    printInfo(info);
  }

  public static class Test1 implements LsqrSolver.A1 {
    public void apply(float[] x, float[] y) {
      y[0] += 2.0f*x[0]+1.0f*x[1];
      y[1] += 3.0f*x[0]+4.0f*x[1];
    }
    public void applyTranspose(float[] y, float[] x) {
      x[0] += 2.0f*y[0]+3.0f*y[1];
      x[1] += 1.0f*y[0]+4.0f*y[1];
    }
  }
  public static void test1() {
    float atol = 0.0f;
    float btol = 0.0f;
    float ctol = 0.0f;
    int maxi = 4;
    LsqrSolver ls = new LsqrSolver(atol,btol,ctol,maxi);
    float[] b = {4.0f,11.0f};
    float[] x = {1.0f,2.0f};
    float[] xtrue = copy(x);
    float damp = 0.0f;
    LsqrSolver.A1 a = new Test1();
    LsqrSolver.Info info = ls.solve(a,damp,b,x);
    println("x="); dump(x);
    println("xtrue="); dump(xtrue);
    printInfo(info);
  }

  private static void println(String s) {
    System.out.println(s);
  }
  private static void printInfo(LsqrSolver.Info info) {
    println("info:");
    println("    stop="+info.stop);
    println("   niter="+info.niter);
    println("   anorm="+info.anorm);
    println("   acond="+info.acond);
    println("   rnorm="+info.rnorm);
    println("  arnorm="+info.arnorm);
    println("   xnorm="+info.xnorm);
  }

  public static void main(String[] args) {
    //test1();
    testH(1,1,1,1,0.0f);
    testH(2,1,1,1,0.0f);
    testH(40,40,4,4,0.0f);
    testH(40,40,4,4,0.01f);
    testH(80,40,4,4,0.01f);
    testH(20,10,1,1,0.001f);
  }
}
