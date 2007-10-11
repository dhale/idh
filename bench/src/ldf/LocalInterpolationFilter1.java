/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local interpolation filter via conjugate gradient iterations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.10.11
 */
public class LocalInterpolationFilter1 {

  /**
   * Constructs a local interpolation filter.
   * @param sigma half-width of smoothing filter used to interpolate.
   * @param small stop when L2 norm of residuals decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalInterpolationFilter1(double sigma, double small, int niter) {
    _sigma = sigma;
    _small = (float)small;
    _niter = niter;
  }

  public void apply(float[] ds, byte[] f, float[] x) {
    int n1 = x.length;
    LocalSmoothingFilter1 lsf = new LocalSmoothingFilter1(_sigma,n1,ds);
    float[] b = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      b[i1] = (f[i1]!=0)?x[i1]:0.0f;
    lsf.applyTranspose(b,b);
    solve(lsf,f,b,x);
    lsf.apply(x,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _sigma; // smoothing filter half-width
  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations

  private void solve(
    LocalSmoothingFilter1 lsf, byte[] f, float[] b, float[] x) 
  {
    Operator1 a = new Smoother1(lsf,f);
    solve(a,b,x);
  }

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator1 {
    public void apply(float[] x, float[] y);
  }

  private static class Smoother1 implements Operator1 {
    Smoother1(LocalSmoothingFilter1 lsf, byte[] f) {
      _lsf = lsf;
      _f = f;
    }
    public void apply(float[] x, float[] y) {
      int n1 = x.length;
      float[] t = new float[n1];
      _lsf.apply(x,y);
      for (int i1=0; i1<n1; ++i1) {
        if (_f[i1]!=0)
          y[i1] = 0.0f;
      }
      _lsf.applyTranspose(y,y);
      Array.sub(x,y,y);
    }
    private LocalSmoothingFilter1 _lsf;
    private byte[] _f;
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   */
  private void solve(Operator1 a, float[] b, float[] x) {
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

  private static void scopy(float[] x, float[] y) {
    Array.copy(x,y);
  }

  // Returns the dot product x'y.
  private static float sdot(float[] x, float[] y) {
    int n1 = x.length;
    float d = 0.0f;
    for (int i1=0; i1<n1; ++i1)
      d += x[i1]*y[i1];
    return d;
  }

  // Computes y = y + ax.
  private static void saxpy(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += a*x[i1];
  }

  // Computes y = x + ay.
  private static void sxpay(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = a*y[i1]+x[i1];
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }
} 
