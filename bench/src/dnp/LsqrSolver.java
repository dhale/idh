/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.logging.Logger;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

/**
 * An implementation of the LSQR algorithm of Paige and Saunders.
 * Let A denote a sparse matrix with m rows and n columns, let b denote
 * a column vector with m rows, and let damp denote a scalar. This solver 
 * iteratively computes a solution x to one of the following problems:
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
 * problems listed above are solved exactly. For example, iterations may
 * end when errors in the residual r = A*x - b are less than the estimated 
 * error in the elements of the column vector b. When constructing a solver, 
 * the user must specify the following parameters:
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
      double anorm, double acond, double rnorm, double arnorm, double xnorm)
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
    public double anorm;
    /**
     * An estimate of cond(Abar), the condition number of Abar. A very
     * high value may indicate that the implementations of A and A' are
     * inconsistent.
     */
    public double acond;
    /**
     * An estimate of the final value of norm(rbar). This is the function
     * that is being minimized. Its value will be small if A*x = b has a
     * solution.
     */
    public double rnorm;
    /**
     * An estimate of the final value of norm(Abar'*rbar). This is the 
     * norm of the residual for the usual normal equations. This value
     * should be small in all cases; it will often be smaller than the
     * true value computed from the solution vector x.
     */
    public double arnorm;
    /**
     * An estimate of the L2 norm of the final solution vector x.
     */
    public double xnorm;
  }

  /**
   * Abstract linear operators A and A' applied to abstract vectors.
   */
  public interface A {

    /**
     * Accumulates the matrix-vector product y += Ax.
     * @param x the input n-vector.
     * @param y the output m-vector.
     */
    public void apply(Vec x, Vec y);

    /**
     * Accumulates the matrix-vector product x += A'y.
     * @param x the input m-vector y.
     * @param y the output n-vector x.
     */
    public void applyTranspose(Vec y, Vec x);
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
  }

  /**
   * Solves the least-squares problem.
   * @param a the linear operator that represents the m x n matrix A.
   * @param b the right-hand-side m-vector.
   * @param x the solution n-vector.
   */
  public Info solve(A a, double damp, Vec b, Vec x) {
    double zero = 0.0;
    double one = 1.0;
    double epsilon = x.epsilon();

    // Initialize.
    Stop stop = Stop.ZERO;
    int niter = 0;
    int nstop = 0;
    double anorm = zero;
    double acond = zero;
    double bbnorm = zero;
    double dampsq = damp*damp;
    double ddnorm = zero;
    double res1 = zero;
    double res2 = zero;
    double rnorm = zero;
    double arnorm = zero;
    double xnorm = zero;
    double xxnorm = zero;
    double cs2 = -one;
    double sn2 = zero;
    double z = zero;
    x.zero();
    Vec u = b.clone(); // or save space with u = b?
    Vec v = x.clone();
    Vec w = null;

    // Set up the first vectors u and v for the bidiagonalization.
    // These vectors satisfy beta*u = b and alfa*v = A'u.
    double alfa = zero;
    double beta = u.norm2();
    if (beta>zero) {
      u.scale(one/beta);
      a.applyTranspose(u,v);
      alfa = v.norm2();
    }
    if (alfa>zero) {
      v.scale(one/alfa);
      w = v.clone();
    }
    arnorm = alfa*beta;

    // If trivial solution x = 0, simply return.
    if (arnorm==zero)
      return new Info(stop,niter,anorm,acond,rnorm,arnorm,xnorm);

    // Otherwise, complete preparation for iterations.
    rnorm = beta;
    double bnorm = beta;
    double phibar = beta;
    double rhobar = alfa;

    // Logging.
    logInit("niter="+niter+" rnorm="+rnorm+" anorm="+anorm);

    // While the number of iterations is less than the maximum, ...
    while (niter<_maxi) {
      ++niter;

      // Perform the next step of the bidiagonalization to obtain
      // the next beta, u, alfa, and v, such that
      //   beta*u = A v - alfa*u 
      //   alfa*v = A'u - beta*v 
      u.scale(-alfa);
      a.apply(v,u);
      beta = u.norm2();
      bbnorm += alfa*alfa+beta*beta+dampsq;
      if (beta>zero) {
        u.scale(one/beta);
        v.scale(-beta);
        a.applyTranspose(u,v);
        alfa = v.norm2();
        if (alfa>zero)
          v.scale(one/alfa);
      }

      // Use a plane rotation to eliminate the damping parameter. This
      // alters the diagonal (rhobar) of the lower-bidiagonal matrix.
      double rhbar2 = rhobar*rhobar+dampsq;
      double rhbar1 = sqrt(rhbar2);
      double cs1 = rhobar/rhbar1;
      double sn1 = damp/rhbar1;
      double psi = sn1*phibar;
      phibar = cs1*phibar;

      // Use a plane rotation to eliminate the subdiagonal element (beta)
      // of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
      double rho = sqrt(rhbar2+beta*beta);
      double cs = rhbar1/rho;
      double sn = beta/rho;
      double theta = sn*alfa;
      double phi = cs*phibar;
      double tau = sn*phi;
      phibar = sn*phibar;
      rhobar = -cs*alfa;

      // Update x and w.
      double dnorm = w.norm2()/rho;
      ddnorm += dnorm*dnorm;
      x.add(1.0,w,phi/rho);
      w.add(-theta/rho,v,1.0);

      // Use a plane rotation to eliminate the super-diagonal element
      // (theta) of the upper-bidiagonal matrix. The use the result to
      // estimate norm(x).
      double delta = sn2*rho;
      double gambar = -cs2*rho;
      double rhs = phi-delta*z;
      double zbar = rhs/gambar;
      double gamma = sqrt(gambar*gambar+theta*theta);
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
      double testr = rnorm/bnorm;
      double testa = (rnorm>zero)?arnorm/(anorm*rnorm):zero;
      double testc = one/acond;
      double rtol = _btol+_atol*anorm*xnorm/bnorm;

      // Logging.
      logIter("niter="+niter+" rnorm="+rnorm+" anorm="+anorm);
      logIter("     testr="+testr+" testa="+testa);
      logIter("     anorm="+anorm+" acond="+acond);

      // The following tests guard against extremely small values of
      // atol, btol, or ctol. (The user may have set any or all of
      // the parameters atol, btol, and ctol to zero.) The effect is
      // equivalent to the user specifying atol = btol = ctol = epsilon.
      if (niter>=_maxi) stop = Stop.MAXI;
      if (testc<epsilon) 
        stop = Stop.CTOL_EPSILON;
      if (testa<epsilon) 
        stop = Stop.ATOL_EPSILON;
      if (testr/(one+anorm*xnorm/bnorm)<epsilon) 
        stop = Stop.RTOL_EPSILON;

      // Check tolerances set by the user.
      if (testc<=_ctol) 
        stop = Stop.CTOL;
      if (testa<=_atol) 
        stop = Stop.ATOL;
      if (testr<=rtol) 
        stop = Stop.RTOL;

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

    // Logging.
    logDone("niter="+niter+" rnorm="+rnorm);
    logDone("     anorm="+anorm+" acond="+acond);
    logDone("     stop="+stop);

    return new Info(stop,niter,anorm,acond,rnorm,arnorm,xnorm);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _atol; // relative error in A
  private double _btol; // relative error in B
  private double _ctol; // reciprocal of upper limit on cond number of Abar
  private int _maxi; // upper limit on number of iterations

  // Logging.
  private static Logger _log = 
    Logger.getLogger(LsqrSolver.class.getName());
  private static void logInit(String s) {
    _log.fine(s);
  }
  private static void logIter(String s) {
    _log.finer(s);
  }
  private static void logDone(String s) {
    _log.fine(s);
  }
}
