/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package util;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * The simplex method for a linear objective function with linear constraints.
 * The simplex algorithm solves linear programs. The standard form of a linear
 * program is to find a set of non-negative variables x that maximize the
 * linear function c'x, subject to a set of linear constraints Ax&le;b. Here,
 * c is an array of n coefficients, b is a set of m upper bounds, and A is an
 * m-by-n matrix of coefficients that together with b define the inequality
 * constraints.
 * <p>
 * This implementation of the simplex algorithm follows that described by
 * Cormen, T.H., C.E. Leiserson, R.L. Rivest, and C. Stein, 2009,
 * Introduction to Algorithms, 3rd edition, MIT Press.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.15
 */
public class SimplexSolver {

  /**
   * Exception thrown if the objective function is unbounded.
   * This exception implies that the linear program is underconstrained, so
   * that the maximized value of the objective function is infinite.
   */
  public static class UnboundedException extends IllegalArgumentException {
    UnboundedException() {
      super("objective function is unbounded");
    }
  }

  /**
   * Exception thrown if no feasible solution exists.
   * This exception implies that the linear program is overconstrained, so
   * that no solution can satisfy simultaneously all of the constraints.
   */
  public static class InfeasibleException extends IllegalArgumentException {
    InfeasibleException() {
      super("no feasible solution exists");
    }
  }

  /**
   * Constructs a solver.
   */
  public SimplexSolver() {
  }

  /**
   * Solves a linear program specified in standard form.
   * @param a input array[m][n] of coefficients in the constraints Ax&le;b.
   * @param b input array[m] of bounds in the constraints Ax&le;b.
   * @param c input array[n] of coefficients in the objective function c'x.
   * @return array[n] containing the solution x.
   */
  public double[] solve(double[][] a, double[] b, double[] c) {
    double[] x = new double[c.length];
    SlackForm sf = new SlackForm(a,b,c);
    sf.solve(x);
    return x;
  }

  /**
   * Solves a linear program specified in standard form.
   * @param a input array[m][n] of coefficients in the constraints Ax&le;b.
   * @param b input array[m] of bounds in the constraints Ax&le;b.
   * @param c input array[n] of coefficients in the objective function c'x.
   * @param x output array[n] containing the solution x.
   * @return value of the objective function c'x for solution x.
   */
  public double solve(double[][] a, double[] b, double[] c, double[] x) {
    SlackForm sf = new SlackForm(a,b,c);
    return sf.solve(x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  /**
   * The slack form of equations to be solved.
   * <p>
   * In CLRS, integers i, j, l and e are indices of variables. Here, i and l
   * denote row indices of basic variables, and j and e are column indices of
   * non-basic variables. This difference enables use of i, j, l and e
   * directly as indices in arrays of coefficients a, b and c. Here, the
   * integer k is the variable index in the range [0:m+n).
   */
  private static class SlackForm {

    // The coefficients for the equations represented by this slack form. Rows
    // and columns correspond to basic and non-basic variables, respectively.
    int m; // number of rows = number of basic variables
    int n; // number of columns = number of non-basic variables
    double[][] a; // m-by-n matrix of coefficients for constraints
    double[] b; // m right-hand-sides for constraints
    double[] c; // n coefficients of objective function
    double v; // value of objective function

    // Mapping from a variable index k to a row index i or column index j.
    // If the index ij[k] is less than n, then ij[k] is a column index j
    // and the variable k is in the non-basic set. Otherwise, ij[k]-n is 
    // a row index i, and the variable k is in the basic set.
    int[] ij; // ij[k] = column index j, or ij[k]-n = row index i

    // Mappings from row and column indices i and j to variable indices k.
    // These mappings must satisfy the following invariants:
    // ij[km[i]] = n+i, for i in [0,m)
    // ij[kn[j]] =   j, for j in [0,n)
    int[] km; // km[i] is index k of variable with row index i
    int[] kn; // kn[i] is index k of variable with column index j

    /**
     * Constructs a slack form from specified standard form.
     */
    SlackForm(double[][] a, double[] b, double[] c) {
      int m = b.length;
      int n = c.length;
      this.m = m;
      this.n = n;
      this.a = copy(a);
      this.b = copy(b);
      this.c = copy(c);
      this.v = 0.0;
      this.ij = new int[n+m];
      this.km = new int[m];
      this.kn = new int[n];
      for (int i=0; i<m; ++i) {
        km[i] = n+i;
        ij[km[i]] = n+i;
      }
      for (int j=0; j<n; ++j) {
        kn[j] = j;
        ij[kn[j]] = j;
      }
    }

    void dumpState() {
      trace("m="+m+" n="+n+" v="+v);
      trace("a="); dump(a);
      trace("b="); dump(b);
      trace("c="); dump(c);
      trace("km="); dump(km);
      trace("kn="); dump(kn);
    }

    /**
     * Pivots this slack form by exchanging basic and non-basic variables.
     * @param l row index in [0:m) of basic variable that is leaving.
     * @param e column index in [0:n) of non-basic variable that is entering.
     */
    void pivot(int l, int e) {

      // Swap indices for leaving and entering variables.
      int kl = km[l];
      int ke = kn[e];
      kn[e] = kl;
      km[l] = ke;

      // Update row and column indices of variables.
      ij[ke] = n+l; // new row index of entering variable
      ij[kl] = e; // new column index of leaving variable

      // Compute coefficients of constraint for the entering variable.
      double oale = 1.0/a[l][e];
      b[l] *= oale;
      for (int j=0; j<n; ++j) {
        if (j!=e) {
          a[l][j] *= oale;
        }
      }
      a[l][e] = oale;

      // Update coefficients of constraints for the other basic variables.
      for (int i=0; i<m; ++i) {
        if (i!=l) {
          b[i] -= a[i][e]*b[l];
          for (int j=0; j<n; ++j) {
            if (j!=e) {
              a[i][j] -= a[i][e]*a[l][j];
            }
          }
          a[i][e] *= -a[l][e];
        }
      }

      // Update coefficients for the objective function.
      v += c[e]*b[l];
      for (int j=0; j<n; ++j) {
        if (j!=e) {
          c[j] -= c[e]*a[l][j];
        }
      }
      c[e] *= -a[l][e];
    }

    /**
     * Solves the equations represented by this slack form.
     * @param x output array[n] containing the solution x.
     * @return the value of the objective function.
     */
    double solve(double[] x) {

      // Initialize this slack form to have a feasible basic solution.
      boolean feasible = init();
      if (!feasible)
        throw new InfeasibleException();
      //dumpState();

      // While we have a column index e for which c[e] > 0, ...
      for (int e=nexte(); e<n; e=nexte()) {

        // The column index e is for the entering variable.
        // Now find the row index l of the leaving variable.
        double dmin = Double.MAX_VALUE;
        int l = m;
        for (int i=0; i<m; ++i) {
          if (a[i][e]>0.0) {
            double d = b[i]/a[i][e];
            if (d<dmin) {
              dmin = d;
              l = i;
            }
          }
        }

        // If a leaving variable was found, do the pivot.
        if (l<m) {
          pivot(l,e);
        } 
        
        // Otherwise, objective function is unbounded.
        else {
          throw new UnboundedException();
        }
        //dumpState();
      }

      // All c[e]>0, so objective function is bounded; get solution x.
      for (int k=0; k<n; ++k) {
        int i = ij[k]-n; // row index, if non-negative
        x[k] = (i>=0)?b[i]:0.0;
      }
      return v;
    }

    /**
     * Returns e (column index of entering variable) for next iteration.
     */
    private int nexte() {
      return bland();
    }

    /**
     * Bland's method returns e for smallest k for which c > 0.
     * If no such index exists, returns n.
     */
    private int bland() {
      for (int k=0; k<m+n; ++k) {
        int e = ij[k];
        if (e<n && c[e]>0.0)
          return e;
      }
      return n;
    }

    /**
     * Returns an auxiliary slack form with n+1 variables. Used when
     * initializing a slack form to have a feasible basic solution.
     * @param i0 row index at which to put the auxiliary variable x0.
     */
    private SlackForm makeAuxiliary(int i0) {
      double[][] aa = new double[m][1+n];
      double[] bb = new double[m];
      double[] cc = new double[1+n];
      cc[0] = -1.0;
      for (int i=0; i<m; ++i) {
        bb[i] = b[i];
        aa[i][0] = -1.0;
        for (int j=0; j<n; ++j) {
          aa[i][j+1] = a[i][j];
        }
      }
      SlackForm sf = new SlackForm(aa,bb,cc);
      //sf.dumpState();
      //trace("i0="+i0);
      sf.pivot(i0,0);
      //sf.dumpState();
      return sf;
    }

    /**
     * Solves this auxiliary slack form. If the auxiliary variable x0 is
     * zero, then a feasible solution exists and this method removes x0
     * from this slack form, so that x1 becomes x0, x2 becomes x1, ....
     * @return true, if feasible; false, otherwise.
     */
    private boolean solveAuxiliary() {
      double[] x = new double[n];
      solve(x);
      //dumpState();

      // If the auxiliary variable x0 is non-zero, then infeasible.
      if (x[0]!=0.0)
        return false;

      // If x0 is basic, use one degenerate pivot to make it non-basic.
      int l = ij[0]-n;
      int e = ij[0];
      if (l>=0) {
        e = 0;
        while (a[l][e]==0.0)
          ++e;
        pivot(l,e);
      }
      //dumpState();

      // Remove the column containing x0 from this auxiliary short form.
      // After removal, some arrays will have an extra element that will be
      // ignored, because here we decrement the number of columns n.
      // Note that indices k of all variables are decremented by one; after
      // removal, the indices k are consistent with those in the short form
      // from which this auxiliary short form was constructed.
      --n;
      for (int i=0; i<m; ++i) {
        for (int j=0,jj=0; j<n; ++j,++jj) {
          if (jj==e) ++jj;
          a[i][j] = a[i][jj];
        }
      }
      for (int j=0,jj=0; j<n; ++j,++jj) {
        if (jj==e) ++jj;
        c[j] = c[jj];
      }
      for (int i=0; i<m; ++i) {
        km[i] -= 1;
        ij[km[i]] = n+i;
      }
      for (int j=0,jj=0; j<n; ++j,++jj) {
        if (jj==e) ++jj;
        kn[j] = kn[jj]-1;
        ij[kn[j]] = j;
      }
      //dumpState();
      return true;
    }

    /**
     * Initialize this slack form so that its basic solution is feasible.
     * @return true, if feasible; false, otherwise.
     */
    private boolean init() {

      // Find index of smallest b.
      int imin = 0;
      double bmin = b[0];
      for (int i=1; i<m; ++i) {
        if (b[i]<bmin) {
          imin = i;
          bmin = b[i];
        }
      }

      // If smallest b is non-negative, then basic solution is feasible.
      if (bmin>=0.0)
        return true;

      // Otherwise, if the basic solution is not feasible, make an auxiliary
      // slack form with n+1 non-basic variables. If this slack form has a
      // feasible solution, then the auxiliary slack form will have a 
      // feasible basic solution.
      SlackForm sf = makeAuxiliary(imin);

      // Solve the auxiliary slack form and remove the auxiliary x0.
      boolean feasible = sf.solveAuxiliary();
      if (!feasible)
        return false;

      // Update the objective function in the auxiliary short form.
      // This part is tricky, because we use coefficients and mappings
      // from both this slack form and the auxiliary slack form.
      sf.v = 0.0;
      for (int i=0; i<m; ++i) { // for all rows in auxiliary slack form, ...
        int k = sf.km[i]; // variable k is basic in auxiliary slack form
        int j = ij[k]; // column index j in this slack form
        if (j<n) { // if non-basic in this slack form, ...
          sf.v += c[j]*sf.b[i];
          for (int jj=0; jj<n; ++jj)
            sf.c[jj] -= c[j]*sf.a[i][jj];
        }
      }
      for (int j=0; j<n; ++j) { // for all columns of this slack form, ...
        int k = kn[j]; // variable k is non-basic in this slack form
        int jj = sf.ij[k]; // column index j in auxiliary slack form
        if (jj<n) // if non-basic in auxiliary slack form, ...
          sf.c[jj] += c[j];
      }

      // Copy everything from the auxiliary slack form to this one.
      for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j)
          a[i][j] = sf.a[i][j];
        b[i] = sf.b[i];
        km[i] = sf.km[i];
        ij[km[i]] = n+i;
      }
      for (int j=0; j<n; ++j) {
        c[j] = sf.c[j];
        kn[j] = sf.kn[j];
        ij[kn[j]] = j;
      }
      v = sf.v;

      //dumpState();
      return true;
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static class LinearProgram {
    public double[][] a;
    public double[] b;
    public double[] c;
    public double[] x;
    public double v;
  }
  private static class LinearProgram1 extends LinearProgram {
    LinearProgram1() {
      a = new double[][]{{1.0, 1.0, 3.0},
                         {2.0, 2.0, 5.0},
                         {4.0, 1.0, 2.0}};
      b = new double[]{30.0, 24.0, 36.0};
      c = new double[]{3.0,1.0,2.0};
      x = new double[]{8.0,4.0,0.0};
      v = 28.0;
    }
  }
  private static class LinearProgram2 extends LinearProgram {
    LinearProgram2() {
      a = new double[][]{{2.0,-1.0},
                         {1.0,-5.0}};
      b = new double[]{2.0,-4.0};
      c = new double[]{2.0,-1.0};
      x = new double[]{1.55556,1.11111};
      v = 2.0;
    }
  }
  private static class LinearProgram3 extends LinearProgram {
    LinearProgram3() {
      a = new double[][]{{ 2.0,-1.0, 2.0},
                         { 2.0,-3.0, 1.0},
                         {-1.0, 1.0,-2.0}};
      b = new double[]{ 4.0,-5.0,-1.0};
      c = new double[]{ 1.0,-1.0, 1.0};
      x = new double[]{ 0.0, 2.8, 3.4};
      v = 0.6;
    }
  }
  private static class LinearProgram4 extends LinearProgram {
    LinearProgram4() {
      a = new double[][]{{2.0,1.0}};
      b = new double[]{2.0};
      c = new double[]{1.0,1.0};
      x = new double[]{0.0,2.0};
      v = 2.0;
    }
  }
  private static class LinearProgram5 extends LinearProgram {
    LinearProgram5() {
      a = new double[][]{{-1.0,-2.0},
                         { 2.0, 4.0}};
      b = new double[]{-1.0, 4.0};
      c = new double[]{1.0,1.0};
      x = new double[]{2.0,0.0};
      v = 2.0;
    }
  }

  public static void main(String[] args) {
    //test(new LinearProgram1());
    //test(new LinearProgram2());
    //test(new LinearProgram3());
    test(new LinearProgram4());
    //test(new LinearProgram5());
  }
  private static void test(LinearProgram lp) {
    int n = lp.c.length;
    double[] x = new double[n];
    SimplexSolver ss = new SimplexSolver();
    double v = ss.solve(lp.a,lp.b,lp.c,x);
    trace("solution v: "+v);
    trace("expected v: "+lp.v); 
    trace("solution x:"); dump(x);
    trace("expected x:"); dump(lp.x);
  }
}
