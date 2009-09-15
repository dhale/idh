/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp.test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;

/**
 * Tests {@link dnp.LsqrSolver}.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public class LsqrSolverTest extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(LsqrSolverTest.class);
    junit.textui.TestRunner.run(suite);
  }

  // A simple first test.
  public void test1() {
    float atol = 0.0f;
    float btol = 0.0f;
    float ctol = 0.0f;
    int maxi = 4;
    LsqrSolver ls = new LsqrSolver(atol,btol,ctol,maxi);
    float[] b = {4.0f,11.0f};
    float[] x = {1.0f,2.0f};
    VecArrayFloat1 vb = new VecArrayFloat1(b);
    VecArrayFloat1 vx = new VecArrayFloat1(x);
    float[] xtrue = vx.clone().getArray();
    float damp = 0.0f;
    LsqrSolver.A a = new Test1();
    LsqrSolver.Info info = ls.solve(a,damp,vb,vx);
    //println("x="); dump(x);
    //println("xtrue="); dump(xtrue);
    //printInfo(info);
    assertTrue(info.stop==LsqrSolver.Stop.RTOL_EPSILON);
    assertEquals(3,info.niter);
    assertTrue(info.rnorm<1.0e-7);
  }
  private static class Test1 implements LsqrSolver.A {
    public void apply(Vec vx, Vec vy) {
      float[] x = ((VecArrayFloat1)vx).getArray();
      float[] y = ((VecArrayFloat1)vy).getArray();
      y[0] += 2.0f*x[0]+1.0f*x[1];
      y[1] += 3.0f*x[0]+4.0f*x[1];
    }
    public void applyTranspose(Vec vy, Vec vx) {
      float[] y = ((VecArrayFloat1)vy).getArray();
      float[] x = ((VecArrayFloat1)vx).getArray();
      x[0] += 2.0f*y[0]+3.0f*y[1];
      x[1] += 1.0f*y[0]+4.0f*y[1];
    }
  }

  public void testPS() {
    //testPS(1,1,1,1,0.0);
    //testPS(2,1,1,1,0.0);
    //testPS(40,40,4,4,0.0);
    //testPS(40,40,4,4,0.01);
    //testPS(80,40,4,4,0.01);
    testPS(20,10,1,1,0.001,LsqrSolver.Stop.ATOL,10,1.0e-10);
  }

  // Test published by Paige and Saunders.
  private static class TestPS implements LsqrSolver.A {
    public int m,n;
    public int ndup,npow;
    public double damp,acond,rnorm;
    public double[] x,b,d,hy,hz,w;
    public TestPS(int m, int n, int ndup, int npow, double damp, double[] x) {
      this.m = m;
      this.n = n;
      this.ndup = ndup;
      this.npow = npow;
      this.damp = damp;
      this.x = x;
      double[] b = this.b = new double[m];
      double[] d = this.d = new double[n];
      double[] hy = this.hy = new double[m];
      double[] hz = this.hz = new double[n];
      double[] w = this.w = new double[m];

      // Make two vectors with norm 1 for Householder transformations.
      double dampsq = damp*damp;
      double fourpi = 4.0*PI;
      double alfa = fourpi/m;
      double beta = fourpi/n;
      for (int i=0; i<m; ++i)
        hy[i] = sin((i+1)*alfa);
      for (int i=0; i<n; ++i)
        hz[i] = cos((i+1)*beta);
      VecArrayDouble1 vhy = new VecArrayDouble1(hy);
      VecArrayDouble1 vhz = new VecArrayDouble1(hz);
      alfa = vhy.norm2();
      beta = vhz.norm2();
      vhy.scale(-1.0/alfa);
      vhz.scale(-1.0/beta);

      // Set the diagonal matrix D, which contains the singular values of A.
      for (int i=0; i<n; ++i) {
        int j = (i+ndup)/ndup;
        double t = j*ndup;
        t /= n;
        d[i] = pow(t,npow);
      }
      acond = sqrt((d[n-1]*d[n-1]+dampsq)/(d[0]*d[0]+dampsq));

      // Compute the residual vector, storing it in b.
      hprod(hz,x,b);
      for (int i=0; i<n; ++i)
        b[i] *= dampsq/d[i];
      double t = 1.0;
      for (int i=n; i<m; ++i) {
        int j = i-n;
        b[i] = (t*j)/m;
        t = -t;
      }
      hprod(hy,b,b);
      VecArrayDouble1 vb = new VecArrayDouble1(b);
      VecArrayDouble1 vx = new VecArrayDouble1(x);
      double bnorm = vb.norm2();
      double xnorm = vx.norm2();
      rnorm = sqrt(bnorm*bnorm+dampsq*xnorm*xnorm);

      // Now compute the true b = r + A*x.
      apply(x,b);
    }
    public void apply(Vec vx, Vec vy) {
      double[] x = ((VecArrayDouble1)vx).getArray();
      double[] y = ((VecArrayDouble1)vy).getArray();
      apply(x,y);
    }
    public void applyTranspose(Vec vy, Vec vx) {
      double[] y = ((VecArrayDouble1)vy).getArray();
      double[] x = ((VecArrayDouble1)vx).getArray();
      applyTranspose(y,x);
    }
    private void apply(double[] x, double[] y) {
      // y += A*x, where A = Hy*D*Hz
      hprod(hz,x,w);
      for (int i=0; i<n; ++i)
        w[i] *= d[i];
      for (int i=n; i<m; ++i)
        w[i] = 0.0;
      hprod(hy,w,w);
      for (int i=0; i<m; ++i)
        y[i] += w[i];
    }
    private void applyTranspose(double[] y, double[] x) {
      // x += A'*y, where A = Hy*D*Hz
      hprod(hy,y,w);
      for (int i=0; i<n; ++i)
        w[i] *= d[i];
      hprod(hz,w,w);
      for (int i=0; i<n; ++i)
        x[i] += w[i];
    }
    private void hprod(double[] hz, double[] x, double[] y) {
      // Householder transformation: y = (I-2*Hz*Hz')*x
      int n = hz.length;
      double s = 0.0;
      for (int i=0; i<n; ++i)
         s += hz[i]*x[i];
      s += s;
      for (int i=0; i<n; ++i)
        y[i] = x[i]-s*hz[i];
    }
  }
  private static void testPS(
    int m, int n, int ndup, int npow, double damp,
    LsqrSolver.Stop stop, int niter, double arnorm) {
    double[] xtrue = new double[n];
    for (int j=0; j<n; ++j)
      xtrue[j] = n-j-1;
    TestPS tps = new TestPS(m,n,ndup,npow,damp,xtrue);
    double acond = tps.acond;
    double rnorm = tps.rnorm;
    double[] b = tps.b;
    double atol = 1.0e-6;
    double btol = atol;
    double ctol = 0.1/acond;
    int maxi = m+n+50;
    LsqrSolver ls = new LsqrSolver(atol,btol,ctol,maxi);
    double[] x = new double[n];
    VecArrayDouble1 vb = new VecArrayDouble1(b);
    VecArrayDouble1 vx = new VecArrayDouble1(x);
    LsqrSolver.Info info = ls.solve(tps,damp,vb,vx);
    //println("x="); dump(x);
    //println("xtrue="); dump(xtrue);
    //printInfo(info);
    assertTrue(info.stop==stop);
    assertEquals(niter,info.niter);
    assertTrue(info.arnorm<arnorm);
  }

  private static void println(String s) {
    System.out.println(s);
  }
  private static void printInfo(LsqrSolver.Info info) {
    println("info:");
    println("    stop="+info.stop);
    println("   niter="+info.niter);
    println("   anorm="+info.anorm);
    println("   acond="+info.acond);
    println("   rnorm="+info.rnorm);
    println("  arnorm="+info.arnorm);
    println("   xnorm="+info.xnorm);
  }
}
