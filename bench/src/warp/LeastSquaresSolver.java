/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import dnp.*;

/**
 * A least-squares solver for linear inverse problems.
 * <p>
 * Let d denote a data vector that is related to a model vector m
 * by a linear operator G, such that d ~ Gm. Let Cd denote the
 * data covariance of d-Gm, and let Cm denote the model covariance
 * of m-m0, where m is the model and m0 is the a-priori model. Then
 * the least-squares solution m to the inverse problem is
 * <pre>
 * m = m0 + Cm G' inv(G Cm G' + Cd) (d - G m0). 
 * </pre>
 * Multiplication by inv(...) is performed by solving a system of
 * linear equations using the iterative method of conjugate-gradients.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.02.14
 */
public class LeastSquaresSolver {

  public interface CovarianceOperator {
    public Vec apply(Vec vx);
  }

  public interface LinearOperator {
    public Vec apply(Vec vx);
    public Vec applyTranspose(Vec vx);
  }

  public LeastSquaresSolver() {
    _g = new IdentityOperator();
    _cm = new IdentityOperator();
    _cd = null;
    _m0 = null;
  }

  public void setLinearOperator(LinearOperator g) {
    _g = g;
  }

  public void setModelCovariance(CovarianceOperator cm) {
    _cm = cm;
  }

  public void setDataCovariance(CovarianceOperator cd) {
    _cd = cd;
  }

  public void setPriorModel(Vec m0) {
    _m0 = m0;
  }

  public Vec solve(Vec d) {
    if (_m0!=null) {
      d = d.clone();
      d.add(1.0,_g.apply(_m0),-1.0);
    }
    OperatorA a = new OperatorA(_g,_cm,_cd);
    Vec x = d.clone();
    CgSolver cgs = new CgSolver(0.01,5);
    CgSolver.Info info = cgs.solve(a,d,x);
    x = _cm.apply(_g.applyTranspose(x));
    if (_m0!=null)
      x.add(1.0,_m0,1.0);
    return x;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private LinearOperator _g;
  private CovarianceOperator _cm;
  private CovarianceOperator _cd;
  private Vec _m0;

  private static class OperatorA implements CgSolver.A {
    OperatorA(
      LinearOperator g, 
      CovarianceOperator cm, 
      CovarianceOperator cd)
    {
      _g = g;
      _cm = cm;
      _cd = cd;
    }
    public void apply(Vec x, Vec y) {
      Vec gx = _g.applyTranspose(x); trace("norm2(gx)="+gx.norm2());
      Vec cmgx = _cm.apply(gx); trace("norm2(cmgx)="+cmgx.norm2());
      Vec gcmgx = _g.apply(cmgx); trace("norm2(gcmgx)="+gcmgx.norm2());
      y.add(0.0,gcmgx,1.0);
      trace("norm2(y)="+y.norm2());
      //y.add(0.0,_g.apply(_cm.apply(_g.applyTranspose(x))),1.0);
      if (_cd!=null)
        y.add(1.0,_cd.apply(x),1.0);
    }
    private LinearOperator _g;
    private CovarianceOperator _cm;
    private CovarianceOperator _cd;
  }
  private static void trace(String s) {
    System.out.println(s);
  }
  private static class IdentityOperator 
    implements LinearOperator,CovarianceOperator
  {
    public Vec apply(Vec x) {
      return x.clone();
    }
    public Vec applyTranspose(Vec x) {
      return x.clone();
    }
  }
}
