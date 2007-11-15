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
 * This version uses a preconditioner based on incomplete Cholesky factors.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.13
 */
public class LocalInterpolationFilterIc {

  /**
   * Constructs a local interpolation filter.
   * @param small stop when L2 norm of residuals decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalInterpolationFilterIc(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  public void apply(LocalDiffusionTensors2 ldt, byte[][] f, float[][] x) {
    Operator2[] op = makeOperators(ldt,f);
    Operator2 a = op[0];
    Operator2 m = op[1];
    float[][] b = makeB(ldt,f,x);
    solve(a,m,b,x);
    //solve(a,b,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations

  /**
   * Copies flagged samples in a specified array; zeros other samples.
   * @param flag the flag for elements to be copied.
   * @param f input array of flags.
   * @param x input array; may be same array as y.
   * @param y output array; may be same array as x.
   */
  private static void copy(int flag, byte[][] f, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        y[i2][i1] = (f[i2][i1]==flag)?x[i2][i1]:0.0f;
      }
    }
  }

  private static float[][] makeB(
    LocalDiffusionTensors2 ldt, byte[][] f, float[][] x) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] t = new float[n2][n1];
    float[][] b = new float[n2][n1];
    copy(1,f,x,t); // t = Kx
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    ldk.apply(ldt,t,b); // b = G'DGKx
    copy(0,f,b,b); // b = MG'DGKx
    Array.sub(t,b,b); // b = (K-MG'DGK)x
    return b;
  }

  private static Operator2[] makeOperators(
    LocalDiffusionTensors2 ldt, byte[][] f) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    float[][][] s = ldk.getCoefficients(ldt);
    float[][] s00 = s[0];
    float[][] s0p = s[1];
    float[][] spm = s[2];
    float[][] sp0 = s[3];
    float[][] spp = s[4];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (f[i2][i1]!=0) {
          s00[i2][i1] = 1.0f;
          s0p[i2][i1] = 0.0f;
          spm[i2][i1] = 0.0f;
          sp0[i2][i1] = 0.0f;
          spp[i2][i1] = 0.0f;
          if (0<i1)
            s0p[i2  ][i1-1] = 0.0f;
          if (0<i2) {
            sp0[i2-1][i1  ] = 0.0f;
            if (i1<n1-1) 
              spm[i2-1][i1+1] = 0.0f;
            if (0<i1) 
              spp[i2-1][i1-1] = 0.0f;
          }
        }
      }
    }
    LocalSpd9Filter a = new LocalSpd9Filter(s);
    LocalSpd9Filter m = new LocalSpd9Filter(s);
    m.factorIC0();
    return new Operator2[]{new A2(a), new M2(m)};
  }

  private static interface Operator2 {
    public void apply(float[][] x, float[][] y);
  }

  private static class A2 implements Operator2 {
    A2(LocalSpd9Filter lsf) {
      _lsf = lsf;
    }
    public void apply(float[][] x, float[][] y) {
      _lsf.apply(x,y);
    }
    private LocalSpd9Filter _lsf;
  }

  private static class M2 implements Operator2 {
    M2(LocalSpd9Filter lsf) {
      _lsf = lsf;
    }
    public void apply(float[][] x, float[][] y) {
      _lsf.applyInverse(x,y);
    }
    private LocalSpd9Filter _lsf;
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
      trace("  iter="+iter+" delta="+delta);
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

  private static void scopy(float[][] x, float[][] y) {
    Array.copy(x,y);
  }
  private static void szero(float[][] x) {
    Array.zero(x);
  }
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

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    //testOperators();
    testSolve();
  }

  private static void testSolve() {
    int n1 = 11;
    int n2 = 11;
    float s0 = 1.0f;
    float s1 = 0.0f;
    //float[][] d0 = Array.randfloat(n1,n2);
    //float[][] d1 = Array.randfloat(n1,n2);
    //float[][] v1 = Array.sub(Array.randfloat(n1,n2),0.5f);
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = new LocalDiffusionTensors2(s0,s1,d0,d1,v1);
    byte[][] f = Array.zerobyte(n1,n2);
    f[n2/2][n1/2] = 1;
    float small = 0.0001f;
    int niter = 100;
    LocalInterpolationFilterIc lif = 
      new LocalInterpolationFilterIc(small,niter);
    float[][] x = Array.zerofloat(n1,n2);
    x[n2/2][n1/2] = 1.0f;
    lif.apply(ldt,f,x);
    edu.mines.jtk.mosaic.SimplePlot.asPixels(x);
  }

  private static void testOperators() {
    int n1 = 9;
    int n2 = 11;
    float s0 = 1.0f;
    float s1 = 0.0f;
    //float[][] d0 = Array.randfloat(n1,n2);
    //float[][] d1 = Array.randfloat(n1,n2);
    //float[][] v1 = Array.sub(Array.randfloat(n1,n2),0.5f);
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = new LocalDiffusionTensors2(s0,s1,d0,d1,v1);
    byte[][] f = Array.zerobyte(n1,n2);
    for (int i1=0; i1<n1; ++i1)
      f[n2/2][i1] = 1;
    Operator2[] op = makeOperators(ldt,f);
    Operator2 a = op[0];
    Operator2 m = op[1];
    testSpd(n1,n2,a);
    testSpd(n1,n2,m);
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
