/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.World;
import edu.mines.jtk.sgl.test.TestFrame;
import edu.mines.jtk.util.ArrayMath;
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

  public void apply(DiffusionTensors3 ldt, byte[][][] f, float[][][] x) {
    Operator3[] op = makeOperators(ldt,f);
    Operator3 a = op[0];
    Operator3 m = op[1];
    float[][][] b = makeB(ldt,f,x);
    solve(a,m,b,x);
    //solve(a,b,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final LocalDiffusionKernel _ldk = 
    new LocalDiffusionKernel(1.0/12.0);

  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations

  ///////////////////////////////////////////////////////////////////////////
  // 2D solver

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
      _lsf.applyApproximateInverse(x,y);
    }
    private LocalSpd9Filter _lsf;
  }

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
    _ldk.apply(ldt,t,b); // b = G'DGKx
    copy(0,f,b,b); // b = MG'DGKx
    ArrayMath.sub(t,b,b); // b = (K-MG'DGK)x
    return b;
  }

  private static Operator2[] makeOperators(
    LocalDiffusionTensors2 ldt, byte[][] f) 
  {
    int n1 = f[0].length;
    int n2 = f.length;

    // First make A = G'DG, which is symmetric positive semidefinite.
    float[][][] s = _ldk.getCoefficients(ldt);

    // Then make A = K+MG'DGM, which should be symmetric positive definite.
    // (It will be SPD iff one or more of the flags in f are non-zero.)
    float[][] s00 = s[0];
    float[][] s0p = s[1];
    float[][] spm = s[2];
    float[][] sp0 = s[3];
    float[][] spp = s[4];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (f[i2][i1]!=0) {
          s00[i2][i1] = 1.0f; // A = (K+MG'DGM)
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
    LocalSpd9Filter lsf = new LocalSpd9Filter(s,0.001);
    //SimplePlot.asPixels(lsf.getMatrix());
    return new Operator2[]{new A2(lsf), new M2(lsf)};
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
    float deltaSmall = deltaBegin*_small*_small;
    //float deltaSmall = sdot(b,b)*_small*_small;
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
    //m.apply(b,s);
    float delta = sdot(r,d);
    float deltaBegin = delta;
    float deltaSmall = deltaBegin*_small*_small;
    //float deltaSmall = sdot(s,s)*_small*_small;
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
    ArrayMath.copy(x,y);
  }
  private static void szero(float[][] x) {
    ArrayMath.zero(x);
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

  ///////////////////////////////////////////////////////////////////////////
  // 3D solver

  private static interface Operator3 {
    public void apply(float[][][] x, float[][][] y);
  }
  private static class A3 implements Operator3 {
    A3(LocalSpd27Filter lsf) {
      _lsf = lsf;
    }
    public void apply(float[][][] x, float[][][] y) {
      _lsf.apply(x,y);
    }
    private LocalSpd27Filter _lsf;
  }
  private static class M3 implements Operator3 {
    M3(LocalSpd27Filter lsf) {
      _lsf = lsf;
    }
    public void apply(float[][][] x, float[][][] y) {
      _lsf.applyApproximateInverse(x,y);
    }
    private LocalSpd27Filter _lsf;
  }

  /**
   * Copies flagged samples in a specified array; zeros other samples.
   * @param flag the flag for elements to be copied.
   * @param f input array of flags.
   * @param x input array; may be same array as y.
   * @param y output array; may be same array as x.
   */
  private static void copy(
    int flag, byte[][][] f, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          y[i3][i2][i1] = (f[i3][i2][i1]==flag)?x[i3][i2][i1]:0.0f;
        }
      }
    }
  }

  private static float[][][] makeB(
    DiffusionTensors3 ldt, byte[][][] f, float[][][] x) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] t = new float[n3][n2][n1];
    float[][][] b = new float[n3][n2][n1];
    copy(1,f,x,t); // t = Kx
    _ldk.apply(ldt,t,b); // b = G'DGKx
    copy(0,f,b,b); // b = MG'DGKx
    ArrayMath.sub(t,b,b); // b = (K-MG'DGK)x
    return b;
  }
  // (K+MG'DGM)x = (K-MG'DGK)x
  // MG'DGMx + MG'DGKx = MG'DGx = 0


  private static Operator3[] makeOperators(
    DiffusionTensors3 ldt, byte[][][] f) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // First make A = G'DG, which is symmetric positive semidefinite.
    float[][][][] s = _ldk.getCoefficients(ldt);

    // Then make A = K+MG'DGM, which should be symmetric positive definite.
    // (It will be SPD iff one or more of the flags in f are non-zero.)
    float[][][] s000 = s[0];
    float[][][] s00p = s[1];
    float[][][] s0pm = s[2];
    float[][][] s0p0 = s[3];
    float[][][] s0pp = s[4];
    float[][][] spmm = s[5];
    float[][][] spm0 = s[6];
    float[][][] spmp = s[7];
    float[][][] sp0m = s[8];
    float[][][] sp00 = s[9];
    float[][][] sp0p = s[10];
    float[][][] sppm = s[11];
    float[][][] spp0 = s[12];
    float[][][] sppp = s[13];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (f[i3][i2][i1]!=0) {
            s000[i3][i2][i1] = 1.0f; // A = (K+MG'DGM)
            s00p[i3][i2][i1] = 0.0f;
            s0pm[i3][i2][i1] = 0.0f;
            s0p0[i3][i2][i1] = 0.0f;
            s0pp[i3][i2][i1] = 0.0f;
            spmm[i3][i2][i1] = 0.0f;
            spm0[i3][i2][i1] = 0.0f;
            spmp[i3][i2][i1] = 0.0f;
            sp0m[i3][i2][i1] = 0.0f;
            sp00[i3][i2][i1] = 0.0f;
            sp0p[i3][i2][i1] = 0.0f;
            sppm[i3][i2][i1] = 0.0f;
            spp0[i3][i2][i1] = 0.0f;
            sppp[i3][i2][i1] = 0.0f;
            if (0<i1)
              s00p[i3][i2][i1-1] = 0.0f;
            if (0<i2) {
              s0p0[i3][i2-1][i1  ] = 0.0f;
              if (i1<n1-1) 
                s0pm[i3][i2-1][i1+1] = 0.0f;
              if (0<i1) 
                s0pp[i3][i2-1][i1-1] = 0.0f;
            }
            if (0<i3) {
              sp00[i3-1][i2][i1] = 0.0f;
              if (i1<n1-1) 
                sp0m[i3-1][i2][i1+1] = 0.0f;
              if (0<i1)
                sp0p[i3-1][i2][i1-1] = 0.0f;
              if (i2<n2-1) {
                spm0[i3-1][i2+1][i1] = 0.0f;
                if (i1<n1-1) 
                  spmm[i3-1][i2+1][i1+1] = 0.0f;
                if (0<i1) 
                  spmp[i3-1][i2+1][i1-1] = 0.0f;
              }
              if (0<i2) {
                spp0[i3-1][i2-1][i1  ] = 0.0f;
                if (i1<n1-1) 
                  sppm[i3-1][i2-1][i1+1] = 0.0f;
                if (0<i1) 
                  sppp[i3-1][i2-1][i1-1] = 0.0f;
              }
            }
          }
        }
      }
    }
    LocalSpd27Filter lsf = new LocalSpd27Filter(s,0.001);
    //SimplePlot.asPixels(lsf.getMatrix());
    return new Operator3[]{new A3(lsf), new M3(lsf)};
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   */
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
    float deltaSmall = deltaBegin*_small*_small;
    //float deltaSmall = sdot(b,b)*_small*_small;
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
    //m.apply(b,s);
    float delta = sdot(r,d);
    float deltaBegin = delta;
    float deltaSmall = deltaBegin*_small*_small;
    //float deltaSmall = sdot(s,s)*_small*_small;
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

  private static void scopy(float[][][] x, float[][][] y) {
    ArrayMath.copy(x,y);
  }
  private static void szero(float[][][] x) {
    ArrayMath.zero(x);
  }
  private static float sdot(float[][][] x, float[][][] y) {
    int n3 = x.length;
    float d = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      d += sdot(x[i3],y[i3]);
    return d;
  }
  private static void saxpy(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpy(a,x[i3],y[i3]);
  }
  private static void sxpay(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sxpay(a,x[i3],y[i3]);
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void plotPixels(float[][] x) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setSize(650,600);
    PixelsView pv = sp.addPixels(x);
    //pv.setClips(-1.0f,1.0f);
    pv.setColorModel(ColorMap.JET);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    int n2 = x.length;
    float[] y = new float[n2];
    for (int i2=0; i2<n2; ++i2)
      y[i2] = x[i2][i2];
    SimplePlot.asPoints(y);
 
  }

  private static void plot3d(float[][][] x) {
    x = ArrayMath.copy(x);
    ImagePanelGroup ipg = new ImagePanelGroup(x);
    ipg.setColorModel(ColorMap.JET);
    World world = new World();
    world.addChild(ipg);
    TestFrame frame = new TestFrame(world);
    frame.setVisible(true);
  }

  private static void testOperators2() {
    int n1 = 9;
    int n2 = 11;
    float s0 = 1.0f;
    float s1 = 0.0f;
    //float[][] d0 = ArrayMath.randfloat(n1,n2);
    //float[][] d1 = ArrayMath.randfloat(n1,n2);
    //float[][] v1 = ArrayMath.sub(ArrayMath.randfloat(n1,n2),0.5f);
    float[][] d0 = ArrayMath.fillfloat(1.0f,n1,n2);
    float[][] d1 = ArrayMath.fillfloat(1.0f,n1,n2);
    float[][] v1 = ArrayMath.fillfloat(sqrt(0.5f),n1,n2);
    LocalDiffusionTensors2 ldt = new LocalDiffusionTensors2(s0,s1,d0,d1,v1);
    byte[][] f = ArrayMath.zerobyte(n1,n2);
    for (int i1=0; i1<n1; ++i1)
      f[n2/2][i1] = 1;
    Operator2[] op = makeOperators(ldt,f);
    Operator2 a = op[0];
    Operator2 m = op[1];
    testSpd(n1,n2,a);
    testSpd(n1,n2,m);
  }

  private static void testOperators3() {
    int n1 = 9;
    int n2 = 10;
    int n3 = 11;
    float theta = FLT_PI*2.0f/8.0f;
    float phi = FLT_PI*0.0f/8.0f;
    DiffusionTensors3 ldt = makeLinearDiffusionTensors3(n1,n2,n3,theta,phi);
    byte[][][] f = ArrayMath.zerobyte(n1,n2,n3);
    for (int i1=0; i1<n1; ++i1)
      f[n3/2][n2/2][i1] = 1;
    Operator3[] op = makeOperators(ldt,f);
    Operator3 a = op[0];
    Operator3 m = op[1];
    testSpd(n1,n2,n3,a);
    testSpd(n1,n2,n3,m);
  }

  private static void testSpd(int n1, int n2, Operator2 a) {
    float[][] x = ArrayMath.sub(ArrayMath.randfloat(n1,n2),0.5f);
    float[][] y = ArrayMath.sub(ArrayMath.randfloat(n1,n2),0.5f);
    float[][] ax = ArrayMath.zerofloat(n1,n2);
    float[][] ay = ArrayMath.zerofloat(n1,n2);
    a.apply(x,ax);
    a.apply(y,ay);
    float xax = sdot(x,ax);
    float yay = sdot(y,ay);
    float yax = sdot(y,ax);
    float xay = sdot(x,ay);
    System.out.println("xax="+xax+" yay="+yay+" (should be positive)");
    System.out.println("yax="+yax+" xay="+xay+" (should be equal)");
  }

  private static void testSpd(int n1, int n2, int n3, Operator3 a) {
    float[][][] x = ArrayMath.sub(ArrayMath.randfloat(n1,n2,n3),0.5f);
    float[][][] y = ArrayMath.sub(ArrayMath.randfloat(n1,n2,n3),0.5f);
    float[][][] ax = ArrayMath.zerofloat(n1,n2,n3);
    float[][][] ay = ArrayMath.zerofloat(n1,n2,n3);
    a.apply(x,ax);
    a.apply(y,ay);
    float xax = sdot(x,ax);
    float yay = sdot(y,ay);
    float yax = sdot(y,ax);
    float xay = sdot(x,ay);
    System.out.println("xax="+xax+" yay="+yay+" (should be positive)");
    System.out.println("yax="+yax+" xay="+xay+" (should be equal)");
  }

  private static void testSolve2() {
    int n1 = 101;
    int n2 = 101;
    float s0 = 0.01f;
    float s1 = 1.00f;
    //float[][] d0 = ArrayMath.randfloat(n1,n2);
    //float[][] d1 = ArrayMath.randfloat(n1,n2);
    //float[][] v1 = ArrayMath.sub(ArrayMath.randfloat(n1,n2),0.5f);
    float theta = FLT_PI*0.5f/8.0f;
    float ctheta = cos(theta);
    float stheta = sin(theta);
    float[][] d0 = ArrayMath.fillfloat(1.0f,n1,n2);
    float[][] d1 = ArrayMath.fillfloat(0.0f,n1,n2);
    float[][] v1 = ArrayMath.fillfloat(stheta,n1,n2);
    LocalDiffusionTensors2 ldt = new LocalDiffusionTensors2(s0,s1,d0,d1,v1);
    byte[][] f = ArrayMath.zerobyte(n1,n2);
    float[][] x = ArrayMath.zerofloat(n1,n2);
    /*
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        //if (i2==0 || i2==n2/2 || i2==n2-1) {
        //if (i2==n2/2 || i2==n2-1) {
        //if (i1==0 || i1==n1-1 || i2==0 || i2==n2-1) {
        //if (i1==1*n1/4 && i2==n2/2) {
        //  //x[i2][i1] = sin(0.1f*FLT_PI*(float)(i1*ctheta-i2*stheta));
        //  x[i2][i1] = 1.0f;
        //  f[i2][i1] = 1;
        //}
        //if (i2==1*n2/3 || i2==2*n2/3) {
        if (i2==1*n2/3) {
          x[i2][i1] = ((i1/10)%2==0)?0.0f:1.0f;
          f[i2][i1] = 1;
        }
      }
    }
    */
    x[0][0] = 0.0f; f[0][0] = 1;
    x[0][n1-1] = 0.0f; f[0][n1-1] = 1;
    x[n2-1][0] = 0.0f; f[n2-1][0] = 1;
    x[n2-1][n1-1] = 0.0f; f[n2-1][n1-1] = 1;
    x[n2/2][n1/2] = 1.0f; f[n2/2][n1/2] = 1;
    plotPixels(x);
    float small = 0.0001f;
    int niter = 1000;
    LocalInterpolationFilterIc lif = 
      new LocalInterpolationFilterIc(small,niter);
    lif.apply(ldt,f,x);
    trace("x min="+ ArrayMath.min(x)+" max="+ ArrayMath.max(x));
    plotPixels(x);
    //LocalDiffusionKernel ldk = new LocalDiffusionKernel(1.0/12.0);
    //float[][] y = ArrayMath.zerofloat(n1,n2);
    //ldk.apply(ldt,x,y);
    //plotPixels(y);
  }

  private static void testSolve3() {
    int n1 = 101;
    int n2 = 101;
    int n3 = 21;
    float theta = FLT_PI*2.0f/8.0f;
    float phi = FLT_PI*0.0f/8.0f;
    float ctheta = cos(theta);
    float stheta = sin(theta);
    DiffusionTensors3 ldt = makePlanarDiffusionTensors3(n1,n2,n3,theta,phi);
    //DiffusionTensors3 ldt = makeLinearDiffusionTensors3(n1,n2,n3,theta,phi);
    byte[][][] f = ArrayMath.zerobyte(n1,n2,n3);
    float[][][] x = ArrayMath.zerofloat(n1,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          //if (i2==n2/2) {
          if (i2==n2/2 && i3==n3/2) {
            //x[i3][i2][i1] = sin(0.1f*FLT_PI*(float)(i1*ctheta-i2*stheta));
            x[i3][i2][i1] = sin(0.1f*FLT_PI*(float)i1);
            f[i3][i2][i1] = 1;
          }
        }
      }
    }
    plot3d(x);
    float small = 0.0001f;
    int niter = 1000;
    LocalInterpolationFilterIc lif = 
      new LocalInterpolationFilterIc(small,niter);
    lif.apply(ldt,f,x);
    trace("x min="+ ArrayMath.min(x)+" max="+ ArrayMath.max(x));
    plot3d(x);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    float[][][] y = ArrayMath.zerofloat(n1,n2,n3);
    ldk.apply(ldt,x,y);
    plot3d(y);
  }

  private static DiffusionTensors3 makeLinearDiffusionTensors3(
    int n1, int n2, int n3, float theta, float phi) 
  {
    DiffusionTensors3 ldt = new DiffusionTensors3(n1,n2,n3,1.0f,1.0f,1.0f);
    float d1 = 1.000f;
    float d2 = 0.000f;
    float d3 = 0.000f;
    float w1 = cos(theta);
    float w2 = sin(theta)*cos(phi);
    float w3 = sin(theta)*sin(phi);
    float[] d = {d1,d2,d3};
    float[] w = {w1,w2,w3};
    float[] u = makeOrthogonalVector(w);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          ldt.setCoefficients(i1,i2,i3,d);
          ldt.setEigenvectorU(i1,i2,i3,u);
          ldt.setEigenvectorW(i1,i2,i3,w);
        }
      }
    }
    return ldt;
  }
  private static DiffusionTensors3 makePlanarDiffusionTensors3(
    int n1, int n2, int n3, float theta, float phi) 
  {
    DiffusionTensors3 ldt = new DiffusionTensors3(n1,n2,n3,1.0f,1.0f,1.0f);
    float d1 = 0.000f;
    float d2 = 1.000f;
    float d3 = 0.000f;
    float u1 = cos(theta);
    float u2 = sin(theta)*cos(phi);
    float u3 = sin(theta)*sin(phi);
    float[] d = {d1,d2,d3};
    float[] u = {u1,u2,u3};
    float[] w = makeOrthogonalVector(u);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          ldt.setCoefficients(i1,i2,i3,d);
          ldt.setEigenvectorU(i1,i2,i3,u);
          ldt.setEigenvectorW(i1,i2,i3,w);
        }
      }
    }
    return ldt;
  }
  private static java.util.Random r = new java.util.Random();
  private static float[] makeRandomCoefficients() {
    float d1 = r.nextFloat();
    float d2 = r.nextFloat();
    float d3 = r.nextFloat();
    float ds = 1.0f/(d1+d2+d3);
    return new float[]{d1*ds,d2*ds,d3*ds};
  }
  private static float[] makeRandomVector() {
    float a = r.nextFloat()-0.5f;
    float b = r.nextFloat()-0.5f;
    float c = r.nextFloat()-0.5f;
    float s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }
  private static float[] makeOrthogonalVector(float[] v1) {
    float a1 = v1[0];
    float b1 = v1[1];
    float c1 = v1[2];
    float a2 = r.nextFloat()-0.5f;
    float b2 = r.nextFloat()-0.5f;
    float c2 = r.nextFloat()-0.5f;
    float d11 = a1*a1+b1*b1+c1*c1;
    float d12 = a1*a2+b1*b2+c1*c2;
    float s = d12/d11;
    float a = a2-s*a1;
    float b = b2-s*b1;
    float c = c2-s*c1;
    s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }

  public static void main(String[] args) {
    //testOperators2();
    //testOperators3();
    testSolve2();
    //testSolve3();
  }
} 
