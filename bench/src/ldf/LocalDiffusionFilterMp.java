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
 * Local anisotropic diffusion filter via minimum-phase factorization.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.12
 */
public class LocalDiffusionFilterMp extends LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   * @param small stop when L2 norm of residuals decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalDiffusionFilterMp(double sigma) {
    super(sigma);
    _sigma = (float)sigma;
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void solveInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // A local inline filter approximated with minimum-phase factors.
  // Factors are tabulated as a function of sigma and theta.
  private static class FactoredFilter2 {
    FactoredFilter2() {
      int ntheta = 21;
      double f
      int maxlag = 4;
      /*
      int nlag = maxlag+2+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
        lag2[ilag] = (ilag<=maxlag)?0:1;
      }
      */
      int nlag = 4+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=1)?ilag:ilag-2-maxlag;
        lag2[ilag] = (ilag<=1)?0:1;
      }
      float[][] w = new float[3][3];
      float[][] v = new float[3][3];
      float[][] t = new float[3][3];
      float[][] r = new float[3][3];
      t[1][1] = 1.0f;
      LocalDipFilter ldf = new LocalDipFilter(Factor.NOT);
      CausalFilter cf = new CausalFilter(lag1,lag2);
      for (int itheta=0; itheta<NTHETA; ++itheta) {
        float theta = FTHETA+itheta*DTHETA;
        for (int iwidth=0; iwidth<NSIGMA; ++iwidth) {
          float width = FSIGMA+iwidth*DSIGMA;
          Array.fill(width,w);
          Array.fill(-sin(theta),v);
          ldf.applyForwardNot(w,v,t,r);
          cf.factorWilsonBurg(100,0.000001f,r);
          _atable[itheta][iwidth] = cf.getA();
        }
      }
      _lcf = new LocalCausalFilter(lag1,lag2);
    }
    void applyForward(float[][] vw, float[][] v1, float[][] x, float[][] y) {
      A2 a2 = new A2(_atable,vw,v1);
      _lcf.apply(a2,x,y);
      _lcf.applyTranspose(a2,y,y);
    }
    void applyInverse(float[][] vw, float[][] v1, float[][] x, float[][] y) {
      A2 a2 = new A2(_atable,vw,v1);
      _lcf.applyInverseTranspose(a2,x,y);
      _lcf.applyInverse(a2,y,y);
    }

    private static final float SIGMA_MIN =  1.0f;
    private static final float SIGMA_MAX = 20.0f;
    private static final int NSIGMA = 20;
    private static final float FSIGMA = SIGMA_MIN;
    private static final float DSIGMA = (SIGMA_MAX-SIGMA_MIN)/(float)(NSIGMA-1);
    private static final float SSIGMA = 0.9999f/DSIGMA;

    private static final float THETA_MIN = -0.5f*FLT_PI;
    private static final float THETA_MAX = 0.5f*FLT_PI
    private static final int NTHETA = 21;
    private static final float FTHETA = THETA_MIN;
    private static final float DTHETA = (THETA_MAX_THETA_MIN)/(float)(NTHETA-1);
    private static final float STHETA = 0.9999f/DTHETA;

    private static class A2 implements LocalCausalFilter.A2 {
      A2(float[][][] atable, float[][] ds, float[][] v1) {
        _at = atable;
        _ds = ds;
        _v1 = v1;
      }
      public void get(int i1, int i2, float[] a) {
        float theta = -asin(_v1[i2][i1]);
        float sigma = min(max(_ds[i2][i1],SIGMA_MIN),SIGMA_MAX);
        float s = (sigma-FSIGMA)*SSIGMA;
        float t = (theta-FTHETA)*STHETA;
        int is = (int)s;
        int it = (int)t;
        float s1 = s-(float)is;
        float t1 = t-(float)it;
        float s0 = 1.0f-s1;
        float t0 = 1.0f-t1;
        float[] a00 = _at[it  ][is  ];
        float[] a01 = _at[it  ][is+1];
        float[] a10 = _at[it+1][is  ];
        float[] a11 = _at[it+1][is+1];
        int n = a00.length;
        for (int j=0; j<n; ++j)
          a[j] = t0*(s0*a00[j]+s1*a01[j])+t1*(s0*a10[j]+s1*a11[j]);
      }
      private float[][][] _at;
      private float[][] _vw;
      private float[][] _v1;
    }
    private LocalCausalFilter _lcf;
    private float[][][] _atable = new float[NTHETA][NSIGMA][];
  }

  private float _sigma; // maximum filter half-width

  // Inverse factored filter.
  private void applyInverseFac(
    float[][] vw, float[][] v1, float[][] x, float[][] y) 
  {
    if (_ff2==null)
      _ff2 = new FactoredFilter2();
    _ff2.applyInverse(vw,v1,x,y);
  }

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
