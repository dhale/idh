/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic interpolation filter via conjugate gradient iterations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.03
 */
public class LocalInterpolationFilter {

  /**
   * Constructs a local interpolation filter.
   * @param aniso anisotropic-to-isotropic diffusion ratio.
   * @param small stop when L2 norm of residuals decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalInterpolationFilter(double aniso, double small, int niter) {
    // aniso = sigmad/sigmae
    double sigmad = (aniso>1.0)?1.0:aniso;
    double sigmae = (aniso<1.0)?1.0:1.0/aniso;
    _dlf = new DirectionalLaplacianFilter(sigmad,sigmae);
    _small = (float)small;
    _niter = niter;
  }

  public void apply(
    float[][] ds, float[][] v1, byte[][] f, float[][] x) 
  {
    int n1 = x[0].length;
    int n2 = x.length;

    // Sub- and super-diagonal of tridiagonal matrix A in arrays a and c.
    float[][] a = new float[n2][n1];
    float[][] c = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2 ) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      for (int i1=0; i1<n1; ++i1) {
        float v1i = v1[i2 ][i1], v2i = 1.0f-v1i*v1i;
        float v1m = v1[i2m][i1], v2m = 1.0f-v1m*v1m;
        float v1p = v1[i2p][i1], v2p = 1.0f-v1p*v1p;
        a[i2][i1] = -0.5f*(v2m+v2i);
        c[i2][i1] = -0.5f*(v2p+v2i);
        if (ds!=null) {
          float dsa = 0.5f*(ds[i2m][i1]+ds[i2][i1]);
          float dsc = 0.5f*(ds[i2p][i1]+ds[i2][i1]);
          a[i2][i1] *= dsa*dsa;
          c[i2][i1] *= dsc*dsc;
        }
      }
    }

    // Diagonal of tridiagonal matrix A in array b.
    float[][] b = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      if (i2<n2-1) {
        for (int i1=0; i1<n1; ++i1)
          b[i2][i1] -= a[i2][i1];
      }
      if (i2>0) {
        for (int i1=0; i1<n1; ++i1)
          b[i2][i1] -= c[i2][i1];
      }
    }

    // Sinc interpolator.
    SincInterpolator si = SincInterpolator.fromErrorAndLength(0.01f,8);
    si.setExtrapolation(SincInterpolator.Extrapolation.ZERO);
    si.setUniformSampling(n1,1.0f,0.0f);
    float[] u1 = Array.rampfloat(0.0f,1.0f,n1);
    float[] t1 = new float[n1];
    float[] x1 = new float[n1];

    // Solve A x = b, where b is x with zeros where samples are missing.
    float[] ot = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      ot[i1] = 1.0f/b[0][i1];
      x[0][i1] = x[0][i1]*ot[i1];
    }
    for (int i2=1; i2<n2; ++i2) {
      maket(v1[i2],u1,t1);
      si.setUniformSamples(x[i2-1]);
      si.interpolate(n1,t1,x1);
      for (int i1=0; i1<n1; ++i1) {
        float ai = (f[i2  ][i1]==0)?a[i2  ][i1]:0.0f;
        float bi = (f[i2  ][i1]==0)?b[i2  ][i1]:1.0f;
        float ci = (f[i2-1][i1]==0)?c[i2-1][i1]:0.0f;
        b[i2][i1] = ci*ot[i1];
        ot[i1] = 1.0f/(bi-ai*b[i2][i1]);
        x[i2][i1] = (x[i2][i1]-ai*x1[i1])*ot[i1];
      }
    }
    for (int i2=n2-1; i2>0; --i2) {
      /*
      maket2(v1[i2],u1,t1);
      si.setUniformSamples(x[i2]);
      si.interpolate(n1,t1,x1);
      for (int i1=0; i1<n1; ++i1)
        x[i2-1][i1] -= b[i2][i1]*x1[i1];
      */
      for (int i1=0; i1<n1; ++i1)
        x1[i1] = -b[i2][i1]*x[i2][i1];
      maket(v1[i2],u1,t1);
      si.setUniformSamples(x[i2-1]);
      si.accumulate(n1,t1,x1);
    }
  }

  public void applyLinear(
    float[][] ds, float[][] es, float[][] v1, byte[][] f, float[][] x) 
  {
    trace("x min="+Array.min(x)+" max="+Array.max(x));
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] t = new float[n2][n1];
    float[][] b = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        t[i2][i1] = (f[i2][i1]!=0)?-x[i2][i1]:0.0f;
      }
    }
    trace("t min="+Array.min(t)+" max="+Array.max(t));
    _dlf.applyLinear(ds,es,v1,t,b);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        b[i2][i1] = (f[i2][i1]!=0)?0.0f:b[i2][i1];
      }
    }
    t = null;
    trace("b min="+Array.min(b)+" max="+Array.max(b));
    solveLinear(ds,es,v1,f,b,x);
  }

  public void applyLinear(
    float sigma, float aniso,
    float[][] ds, float[][] es, float[][] v1, 
    byte[][] f, float[][] x) 
  {
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(sigma,aniso);
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] t = new float[n2][n1];
    float[][] b = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        t[i2][i1] = (f[i2][i1]!=0)?x[i2][i1]:0.0f;
      }
    }
    lsf.applyPassTranspose(ds,es,v1,t,b);
    solve(lsf,ds,es,v1,f,b,t);
    lsf.applyPass(ds,es,v1,t,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean PARALLEL = true;

  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations
  private DirectionalLaplacianFilter _dlf;

  private static void maket(float[] v1, float[] u1, float[] t1) {
    int n1 = v1.length;
    float t1min = 0.0f;
    float t1max = (float)(n1-1);
    for (int i1=0; i1<n1; ++i1) {
      float v1i = v1[i1];
      float v2i = sqrt(1.0f-v1i*v1i);
      float vi = v1i/v2i;
      t1[i1] = max(t1min,min(t1max,u1[i1]-vi));
    }
  }

  private static void maket2(float[] v1, float[] u1, float[] t1) {
    int n1 = v1.length;
    float t1min = 0.0f;
    float t1max = (float)(n1-1);
    for (int i1=0; i1<n1; ++i1) {
      float v1i = v1[i1];
      float v2i = sqrt(1.0f-v1i*v1i);
      float vi = v1i/v2i;
      t1[i1] = max(t1min,min(t1max,u1[i1]+vi));
    }
  }

  private void solve(
    LocalSmoothingFilter lsf, 
    float[][] ds, float[][] es, float[][] v1, 
    byte[][] f, float[][] x, float[][] y) 
  {
    Operator2 a = new LinearSmoother2(lsf,ds,es,v1,f);
    solve(a,x,y);
  }

  private void solveLinear(
    float[][] ds, float[][] es, float[][] v1, 
    byte[][] f, float[][] x, float[][] y) 
  {
    solveLinearSimple(ds,es,v1,f,x,y);
  }

  private void solveLinear(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    solveLinearSimple(is,iw,x,y);
  }

  private void solvePlanar(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    solvePlanarSimple(is,iu,x,y);
  }

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator2 {
    public void apply(float[][] x, float[][] y);
  }
  private static interface Operator3 {
    public void apply(float[][][] x, float[][][] y);
  }

  private static class LinearSmoother2 implements Operator2 {
    LinearSmoother2(
      LocalSmoothingFilter lsf, 
      float[][] ds, float[][] es, float[][] v1, byte[][] f) 
    {
      _lsf = lsf;
      _ds = ds;
      _es = es;
      _v1 = v1;
      _f = f;
    }
    public void apply(float[][] x, float[][] y) {
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] t = new float[n2][n1];
      _lsf.applyPass(_ds,_es,_v1,x,t);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (_f[i2][i1]==0)
            t[i2][i1] = 0.0f;
        }
      }
      _lsf.applyPassTranspose(_ds,_es,_v1,t,y);
    }
    private float[][] _ds,_es,_v1;
    private byte[][] _f;
    private LocalSmoothingFilter _lsf;
  }

  private static class LinearOperator2 implements Operator2 {
    LinearOperator2(
      DirectionalLaplacianFilter dlf, 
      float[][] ds, float[][] es, float[][] v1, byte[][] f) 
    {
      _dlf = dlf;
      _ds = ds;
      _es = es;
      _v1 = v1;
      _f = f;
    }
    public void apply(float[][] x, float[][] y) {
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] t = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[i2][i1] = (_f[i2][i1]!=0)?0.0f:x[i2][i1];
        }
      }
      szero(y);
      _dlf.applyLinear(_ds,_es,_v1,t,y);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          y[i2][i1] = (_f[i2][i1]!=0)?0.0f:y[i2][i1];
        }
      }
    }
    private float[][] _ds,_es,_v1;
    private byte[][] _f;
    private DirectionalLaplacianFilter _dlf;
  }
  private static class LinearOperator3 implements Operator3 {
    LinearOperator3(
      DirectionalLaplacianFilter dlf, byte[][][] is, short[][][] iw) 
    {
      _dlf = dlf;
      _is = is;
      _iw = iw;
    }
    public void apply(float[][][] x, float[][][] y) {
      scopy(x,y);
      _dlf.applyLinear(_is,_iw,x,y);
    }
    private byte[][][] _is;
    private short[][][] _iw;
    private DirectionalLaplacianFilter _dlf;
  }
  private static class PlanarOperator3 implements Operator3 {
    PlanarOperator3(
      DirectionalLaplacianFilter dlf, byte[][][] is, short[][][] iu) 
    {
      _dlf = dlf;
      _is = is;
      _iu = iu;
    }
    public void apply(float[][][] x, float[][][] y) {
      scopy(x,y);
      _dlf.applyPlanar(_is,_iu,x,y);
    }
    private byte[][][] _is;
    private short[][][] _iu;
    private DirectionalLaplacianFilter _dlf;
  }

  private void solveLinearSimple(
    float[][] ds, float[][] es, float[][] v1, 
    byte[][] f, float[][] x, float[][] y) 
  {
    Operator2 a = new LinearOperator2(_dlf,ds,es,v1,f);
    solve(a,x,y);
  }
  private void solveLinearSimple(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    Operator3 a = new LinearOperator3(_dlf,is,iw);
    solve(a,x,y);
  }
  private void solvePlanarSimple(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    Operator3 a = new PlanarOperator3(_dlf,is,iu);
    solve(a,x,y);
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   */
  private void solve(Operator2 a, float[][] b, float[][] x) {
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
  private void solve(Operator3 a, float[][][] b, float[][][] x) {
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
  private void solve(Operator2 a, Operator2 m, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] d = new float[n2][n1];
    float[][] q = new float[n2][n1];
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
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
      trace("  iter="+iter+" delta="+delta);
    }
    trace("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }
  private void solve(Operator3 a, Operator3 m, float[][][] b, float[][][] x) {
    int n1 = b[0][0].length;
    int n2 = b[0].length;
    int n3 = b.length;
    float[][][] d = new float[n3][n2][n1];
    float[][][] q = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    float[][][] s = new float[n3][n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
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

  private static void scopy(float[][] x, float[][] y) {
    Array.copy(x,y);
  }
  private static void scopy(float[][][] x, float[][][] y) {
    if (PARALLEL) {
      scopyP(x,y);
    } else {
      scopyS(x,y);
    }
  }
  private static void scopyS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      scopy(x[i3],y[i3]);
  }
  private static void scopyP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
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

  private static void szero(float[][] x) {
    Array.zero(x);
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
  private static float sdot(float[][][] x, float[][][] y) {
    if (PARALLEL) {
      return sdotP(x,y);
    } else {
      return sdotS(x,y);
    }
  }
  private static float sdotS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    float d = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      d += sdot(x[i3],y[i3]);
    return d;
  }
  private static float sdotP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicFloat ad = new AtomicFloat(0.0f);
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
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
  private static void saxpy(float a, float[][][] x, float[][][] y) {
    if (PARALLEL) {
      saxpyP(a,x,y);
    } else {
      saxpyS(a,x,y);
    }
  }
  private static void saxpyS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpy(a,x[i3],y[i3]);
  }
  private static void saxpyP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
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
  private static void sxpay(float a, float[][][] x, float[][][] y) {
    if (PARALLEL) {
      sxpayP(a,x,y);
    } else {
      sxpayS(a,x,y);
    }
  }
  private static void sxpayS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sxpay(a,x[i3],y[i3]);
  }
  private static void sxpayP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
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

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    testOperators();
  }

  private static void testOperators() {
    int n1 = 100;
    int n2 = 101;
    DirectionalLaplacianFilter dlf = new DirectionalLaplacianFilter(1.0);
    float[][] ds = Array.randfloat(n1,n2);
    float[][] es = Array.randfloat(n1,n2);
    float[][] v1 = Array.randfloat(n1,n2);
    byte[][] f = Array.zerobyte(n1,n2);
    Operator2 a = new LinearOperator2(dlf,ds,es,v1,f);
    testSpd(n1,n2,a);
  }

  private static void testSpd(int n1, int n2, Operator2 a) {
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
    System.out.println("xax="+xax+" yay="+yay+" (should be positive)");
    System.out.println("yax="+yax+" xay="+xay+" (should be equal)");
  }
} 
