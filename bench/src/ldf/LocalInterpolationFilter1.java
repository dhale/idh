/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.ArrayMath;

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
    _sigma = (float)sigma;
    _small = (float)small;
    _niter = niter;
  }

  public void apply(float[] ds, byte[] f, float[] x) {
    int n1 = x.length;

    // Sub- and super-diagonal of tridiagonal matrix A in arrays a and c.
    float[] a = new float[n1];
    float[] c = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      a[i1] = c[i1] = -1.0f;
    if (ds!=null) {
      float dsa,dsc;
      dsc = 0.50f*(ds[1]+ds[0]);
      c[0] *= dsc*dsc;
      for (int i1=1; i1<n1-1; ++i1) {
        dsa = 0.50f*(ds[i1-1]+ds[i1]);
        dsc = 0.50f*(ds[i1+1]+ds[i1]);
        a[i1] *= dsa*dsa;
        c[i1] *= dsc*dsc;
      }
      dsa = 0.50f*(ds[n1-1]+ds[n1-2]);
      a[n1-1] *= dsa*dsa;
    }

    // Diagonal of tridiagonal matrix A in array b.
    float[] b = new float[n1];
    b[0] = -c[0];
    for (int i1=1; i1<n1-1; ++i1) {
      b[i1] = -(a[i1]+c[i1]);
    }
    b[n1-1] = -a[n1-1];

    // Solve A x = b, where b is x with zeros for missing samples.
    float ot = 1.0f/b[0];
    x[0] = x[0]*ot;
    for (int i1=1; i1<n1; ++i1) {
      float ai = (f[i1  ]==0)?a[i1  ]:0.0f;
      float bi = (f[i1  ]==0)?b[i1  ]:1.0f;
      float ci = (f[i1-1]==0)?c[i1-1]:0.0f;
      b[i1] = ci*ot;
      ot = 1.0f/(bi-ai*b[i1]);
      x[i1] = (x[i1]-ai*x[i1-1])*ot;
    }
    for (int i1=n1-1; i1>0; --i1)
      x[i1-1] -= b[i1]*x[i1];
  }

  public void applySmoother(float[] ds, byte[] f, float[] x) {
    int n1 = x.length;
    LocalSmoothingFilter1 lsf = new LocalSmoothingFilter1(_sigma,n1,ds);
    float[] b = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      b[i1] = (f[i1]!=0)?x[i1]:0.0f;
    lsf.applyTranspose(b,b);
    solve(lsf,f,b,x);
    lsf.apply(x,x);
  }

  public void applyCg(float[] ds, byte[] f, float[] x) {
    int n1 = x.length;
    float[] t = ArrayMath.neg(x);
    float[] b = new float[n1];
    mask(0,f,t,t);
    laplacian(ds,t,b);
    mask(1,f,b,b);
    solve(ds,f,b,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // smoothing filter half-width
  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations

  private void solve(
    LocalSmoothingFilter1 lsf, byte[] f, float[] b, float[] x) 
  {
    Operator1 a = new Smoother1(lsf,f);
    solve(a,b,x);
  }

  private void solve(
    float[] ds, byte[] f, float[] b, float[] x) 
  {
    Operator1 a = new Laplacian1(f,ds);
    solve(a,b,x);
  }

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator1 {
    public void apply(float[] x, float[] y);
  }

  private static class Laplacian1 implements Operator1 {
    Laplacian1(byte[] f, float[] ds) {
      _f = f;
      _ds = ds;
      _t = new float[f.length];
    }
    public void apply(float[] x, float[] y) {
      mask(1,_f,x,_t);
      laplacian(_ds,_t,y);
      mask(1,_f,y,y);
    }
    private byte[] _f;
    private float[] _ds,_t;
  }
  private static void laplacian(float[] ds, float[] x, float[] y) {
    int n1 = x.length;
    if (ds!=null) {
      for (int i1=1; i1<n1; ++i1) {
        float di = 0.5f*(ds[i1]+ds[i1-1]);
        float xy = di*(x[i1]-x[i1-1]);
        y[i1  ]  = xy;
        y[i1-1] -= xy;
      }
    } else {
      for (int i1=1; i1<n1; ++i1) {
        float xy = x[i1]-x[i1-1];
        y[i1  ]  = xy;
        y[i1-1] -= xy;
      }
    }
  }
  private static void mask(int zero, byte[] f, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = (f[i1]==zero)?0.0f:x[i1];
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
      ArrayMath.sub(x,y,y);
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
    ArrayMath.copy(x,y);
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
