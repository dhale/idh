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
 * Local anisotropic diffusion filter via conjugate gradient iterations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.05
 */
public class LocalDiffusionFilterCg extends LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   * @param small stop when sum of residuals squared decreases by this factor.
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

  protected void solveInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveInlineSsor(ds,v1,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  public static void main(String[] args) {
    testOperators();
  }

  private static void testOperators() {
    int n1 = 100;
    int n2 = 101;
    DirectionalLaplacianFilter dlf = new DirectionalLaplacianFilter(1.0);
    float[][] ds = Array.randfloat(n1,n2);
    float[][] v1 = Array.randfloat(n1,n2);
    Operator a = new InlineOperator(dlf,ds,v1);
    Operator m = new InlineSsorOperator(dlf,ds,v1);
    testSpd(n1,n2,a);
    testSpd(n1,n2,m);
  }

  private static void testSpd(int n1, int n2, Operator a) {
    float[][] x = Array.sub(Array.randfloat(n1,n2),0.5f);
    float[][] y = Array.sub(Array.randfloat(n1,n2),0.5f);
    float[][] ax = Array.zerofloat(n1,n2);
    float[][] ay = Array.zerofloat(n1,n2);
    a.apply(x,ax);
    a.apply(y,ay);
    float xax = sdot(x,ax);
    float yay = sdot(y,ay);
    float yax = sdot(y,ax);
    float xay = sdot(x,ay);
    System.out.println("xax="+xax+" yay="+yay);
    System.out.println("yax="+yax+" xay="+xay);
  }

  private float _sigma; // maximum filter half-width
  private float _small; // stop iterations when rr decreases by this factor
  private int _niter; // number of iterations
  private DirectionalLaplacianFilter _dlf;

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator {
    public void apply(float[][] x, float[][] y);
  }

  private static class InlineOperator implements Operator {
    InlineOperator(
      DirectionalLaplacianFilter dlf, float[][] ds, float[][] v1) 
    {
      _dlf = dlf;
      _ds = ds;
      _v1 = v1;
    }
    public void apply(float[][] x, float[][] y) {
      Array.copy(x,y);
      _dlf.applyInline(_ds,_v1,x,y);
    }
    private float[][] _ds,_v1;
    private DirectionalLaplacianFilter _dlf;
  }

  private static class InlineSsorOperator implements Operator {
    InlineSsorOperator(
      DirectionalLaplacianFilter dlf, float[][] ds, float[][] v1) 
    {
      _s33 = dlf.makeInlineStencil33(ds,v1);
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

  private void solveInlineSimple(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    Operator a = new InlineOperator(_dlf,ds,v1);
    solve(a,x,y);
  }

  private void solveInlineSsor(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    Operator a = new InlineOperator(_dlf,ds,v1);
    Operator m = new InlineSsorOperator(_dlf,ds,v1);
    solve(a,m,x,y);
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   */
  private void solve(Operator a, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] d = new float[n2][n1];
    float[][] q = new float[n2][n1];
    float[][] r = new float[n2][n1];
    Array.copy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,r,q); // q = b-Ax
    Array.copy(r,d);
    float delta = sdot(r,r);
    float deltaBegin = delta;
    float deltaSmall = delta*_small;
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
  private void solve(Operator a, Operator m, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] d = new float[n2][n1];
    float[][] q = new float[n2][n1];
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    Array.copy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,r,q); // q = b-Ax
    m.apply(r,d);
    float delta = sdot(r,d);
    float deltaBegin = delta;
    float deltaSmall = delta*_small;
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

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }
} 
