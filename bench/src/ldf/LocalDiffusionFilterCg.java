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
  //
  protected void solveInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveInlineSsor(ds,v1,x,y);
  }

  private void solveInlineSimple(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    final float[][] dsf = ds;
    final float[][] v1f = v1;
    Operator a = new Operator() {
      public void apply(float[][] x, float[][] y) {
        Array.copy(x,y);
        _dlf.applyInline(dsf,v1f,x,y);
      }
    };
    solveCg(a,x,y);
  }

  private void solveInlineSsor(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;

    // The operator I+G'DG
    final float[][] dsf = ds;
    final float[][] v1f = v1;
    Operator a = new Operator() {
      public void apply(float[][] x, float[][] y) {
        Array.copy(x,y);
        _dlf.applyInline(dsf,v1f,x,y);
      }
    };

    // SSOR preconditioner.
    final DirectionalLaplacianFilter.Stencil33 s33 =
      _dlf.makeInlineStencil33(ds,v1);
    Operator m = new Operator() {
      public void apply(float[][] x, float[][] y) {
        /*
        Array.copy(x,y);
        */
        int n1 = x[0].length;
        int n2 = x.length;
        int n1m = n1-1;
        int n2m = n2-1;
        float[][] a = new float[n1][9];
        for (int i2=0; i2<n2; ++i2) {
          int i2m = max(i2-1,0);
          int i2p = min(i2+1,n2m);
          float[] x0 = x[i2 ];
          float[] ym = y[i2m];
          float[] y0 = y[i2 ];
          float[] yp = y[i2p];
          s33.get(i2,a);
          for (int i1=0; i1<n1; ++i1) {
            int i1m = max(i1-1,0);
            int i1p = min(i1+1,n1m);
            float[] ai = a[i1];
            float d = 1.0f+ai[4];
            y0[i1] = (x0[i1] - 
                      ai[0]*ym[i1m] - ai[3]*y0[i1m] - ai[6]*yp[i1m] -
                      ai[1]*ym[i1 ]                 - ai[7]*yp[i1 ] -
                      ai[2]*ym[i1p] - ai[5]*y0[i1p] - ai[8]*yp[i1p])/d;
          }
        }
        for (int i2=n2-1; i2>=0; --i2) {
          int i2m = max(i2-1,0);
          int i2p = min(i2+1,n2m);
          float[] x0 = x[i2 ];
          float[] ym = y[i2m];
          float[] y0 = y[i2 ];
          float[] yp = y[i2p];
          s33.get(i2,a);
          for (int i1=n1-1; i1>=0; --i1) {
            int i1m = max(i1-1,0);
            int i1p = min(i1+1,n1m);
            float[] ai = a[i1];
            float d = 1.0f+ai[4];
            y0[i1] = (x0[i1] - 
                      ai[0]*ym[i1m] - ai[3]*y0[i1m] - ai[6]*yp[i1m] -
                      ai[1]*ym[i1 ]                 - ai[7]*yp[i1 ] -
                      ai[2]*ym[i1p] - ai[5]*y0[i1p] - ai[8]*yp[i1p])/d;
          }
        }
      }
    };
    solveCg(a,m,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // maximum filter half-width
  private float _small; // stop iterations when rr decreases by this factor
  private int _niter; // number of iterations
  private DirectionalLaplacianFilter _dlf;

  /**
   * Symmetric positive-definite operator to be inverted.
   */
  private static interface Operator {
    public void apply(float[][] x, float[][] y);
  }

  /**
   * Solves Ay = x via conjugate gradient iterations.
   */
  private void solveCg(Operator a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    float[][] t = new float[n2][n1];
    a.apply(y,t);
    double rr = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], r2 = r[i2], s2 = s[i2], t2 = t[i2];
      for (int i1=0; i1<n1; ++i1) {
        float ri = x2[i1]-t2[i1];
        r2[i1] = ri;
        s2[i1] = ri;
        rr += ri*ri;
      }
    }
    trace("solveCg: r="+rr);
    int miter;
    double rrsmall = rr*_small;
    for (miter=0; miter<_niter && rr>rrsmall; ++miter) {
      a.apply(s,t);
      double st = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] s2 = s[i2], t2 = t[i2];
        for (int i1=0; i1<n1; ++i1)
          st += s2[i1]*t2[i1];
      }
      float alpha = (float)(rr/st);
      double rrold = rr;
      rr = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] y2 = y[i2], r2 = r[i2], s2 = s[i2], t2 = t[i2];
        for (int i1=0; i1<n1; ++i1) {
          y2[i1] += alpha*s2[i1];
          r2[i1] -= alpha*t2[i1];
          rr += r2[i1]*r2[i1];
        }
      }
      if (rr<=rrsmall)
        break;
      float beta = (float)(rr/rrold);
      for (int i2=0; i2<n2; ++i2) {
        float[] r2 = r[i2], s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = r2[i1]+beta*s2[i1];
      }
    }
    trace("  miter="+miter+" rr="+rr);
  }

  /**
   * Solves Ay = x via conjugate gradient iterations with preconditioner M.
   */
  private void solveCg(Operator a, Operator m, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    float[][] t = new float[n2][n1];
    float[][] w = new float[n2][n1];
    a.apply(y,t);
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], r2 = r[i2], t2 = t[i2];
      for (int i1=0; i1<n1; ++i1) {
        r2[i1] = x2[i1]-t2[i1];
      }
    }
    m.apply(r,s);
    double rr = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] r2 = r[i2], s2 = s[i2];
      for (int i1=0; i1<n1; ++i1) {
        rr += r2[i1]*s2[i1];
      }
    }
    trace("solveCgPc: r="+rr);
    int miter;
    double rrsmall = rr*_small;
    for (miter=0; miter<_niter && rr>rrsmall; ++miter) {
      a.apply(s,t);
      double st = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] s2 = s[i2], t2 = t[i2];
        for (int i1=0; i1<n1; ++i1)
          st += s2[i1]*t2[i1];
      }
      float alpha = (float)(rr/st);
      for (int i2=0; i2<n2; ++i2) {
        float[] y2 = y[i2], r2 = r[i2], s2 = s[i2], t2 = t[i2];
        for (int i1=0; i1<n1; ++i1) {
          y2[i1] += alpha*s2[i1];
          r2[i1] -= alpha*t2[i1];
        }
      }
      m.apply(r,w);
      double rrold = rr;
      rr = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] r2 = r[i2], w2 = w[i2];
        for (int i1=0; i1<n1; ++i1) {
          rr += r2[i1]*w2[i1];
        }
      }
      if (rr<=rrsmall)
        break;
      float beta = (float)(rr/rrold);
      for (int i2=0; i2<n2; ++i2) {
        float[] w2 = w[i2], s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = w2[i1]+beta*s2[i1];
      }
    }
    trace("  miter="+miter+" rr="+rr);
  }

  /*
  private void applyInverseSpdSsor(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    float[][] w = new float[n2][n1]; // s
    Array.zero(y);
    Array.copy(x,r);
    applyInverseSsor(u1,u2,r,s);
    float rr = dot(r,s);
    float small = rr*0.00001f;
    trace("small="+small);
    int niter;
    for (niter=0; niter<100 && rr>small; ++niter) {
      applyForwardSpd(u1,u2,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      applyInverseSsor(u1,u2,r,w);
      float rrold = rr;
      rr = dot(r,w);
      float beta = rr/rrold;
      for (int i2=0; i2<n2; ++i2) {
        float[] w2 = w[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = w2[i1]+beta*s2[i1];
      }
      trace("niter="+niter+" rr="+rr);
    }
  }
  */

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }
} 
