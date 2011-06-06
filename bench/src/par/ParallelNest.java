/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
*****************************************************************************
Benchmark nested parallel loops with edu.mines.jtk.util.Parallel.
@author Dave Hale, Colorado School of Mines
@version 2011.06.04
****************************************************************************/

import java.util.Random;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

public class ParallelNest {

  public static void main(String[] args) {
    int[][] ns = {
      {1000, 50000,     5},
      {1000,  5000,    50},
      {1000,   500,   500},
      {1000,    50,  5000},
      {1000,     5, 50000}};
    for (int i=0; i<5; i+=1) {
      int n1 = ns[i][0], n2 = ns[i][1], n3 = ns[i][2];
      benchSolr(n1,n2,n3);
    }
  }

  private static final int maxtest = 3;
  private static final double maxtime = 3.0;

  private static void solrS(
    float a1, float a2, float b0, float b1, float b2,
    float[] x, float[] y)
  {
    int n = x.length;
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim2 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = b0*xi+b1*xim1+b2*xim2-a1*yim1-a2*yim2;
      y[i] = yi;
      yim2 = yim1;
      yim1 = yi;
      xim2 = xim1;
      xim1 = xi;
    }
  }
  private static void solrS(
    float a1, float a2, float b0, float b1, float b2,
    float[][] x, float[][] y) 
  {
    int n = x.length;
    for (int i=0; i<n; ++i)
      solrS(a1,a2,b0,b1,b2,x[i],y[i]);
  }
  private static void solrS(
    float a1, float a2, float b0, float b1, float b2,
    float[][][] x, float[][][] y) 
  {
    int n = x.length;
    for (int i=0; i<n; ++i)
      solrS(a1,a2,b0,b1,b2,x[i],y[i]);
  }
  private static void solrP(
    final float a1, final float a2, 
    final float b0, final float b1, final float b2,
    final float[][] x, final float[][] y) 
  {
    int n = x.length;
    loop(n,new LoopInt() {
      public void compute(int i) {
        solrS(a1,a2,b0,b1,b2,x[i],y[i]);
      }
    });
  }
  private static void solrP(
    final float a1, final float a2, 
    final float b0, final float b1, final float b2,
    final float[][][] x, final float[][][] y) 
  {
    int n = x.length;
    loop(n,new LoopInt() {
      public void compute(int i) {
        solrP(a1,a2,b0,b1,b2,x[i],y[i]);
      }
    });
  }
  private static void benchSolr(int n1, int n2, int n3) {
    System.out.println("Array solr: n1="+n1+" n2="+n2+" n3="+n3);
    int niter,rate;
    double mflop2 = 9.0e-6*n1*n2;
    double mflop3 = 9.0e-6*n1*n2*n3;
    Stopwatch sw = new Stopwatch();
    float a1 = -1.8f;
    float a2 = 0.81f;
    float b0 = 2.0f;
    float b1 = -3.2f;
    float b2 = 1.28f;
    float[][][] x = fillfloat(1.0f,n1,n2,n3);
    float[][][] ys = fillfloat(1.0f,n1,n2,n3);
    float[][][] yp = fillfloat(1.0f,n1,n2,n3);
    for (int ntest=0; ntest<maxtest; ++ntest) {
      sw.restart();
      for (niter=0; sw.time()<maxtime; ++niter)
        solrS(a1,a2,b0,b1,b2,x[niter%n3],ys[niter%n3]);
      sw.stop();
      rate = (int)((niter*mflop2)/sw.time());
      System.out.println("2D S: rate = "+rate);
      sw.restart();
      for (niter=0; sw.time()<maxtime; ++niter)
        solrP(a1,a2,b0,b1,b2,x[niter%n3],yp[niter%n3]);
      sw.stop();
      rate = (int)((niter*mflop2)/sw.time());
      System.out.println("2D P: rate = "+rate);
      assertEquals(ys[0],yp[0]);
      sw.restart();
      for (niter=0; sw.time()<maxtime; ++niter)
        solrS(a1,a2,b0,b1,b2,x,ys);
      sw.stop();
      rate = (int)((niter*mflop3)/sw.time());
      System.out.println("3D S: rate = "+rate);
      sw.restart();
      for (niter=0; sw.time()<maxtime; ++niter)
        solrP(a1,a2,b0,b1,b2,x,yp);
      sw.stop();
      rate = (int)((niter*mflop3)/sw.time());
      System.out.println("3D P: rate = "+rate);
      assertEquals(ys,yp);
    }
  }

  private static void assertEquals(float[] e, float[] a) {
    int n = e.length;
    for (int i=0; i<n; ++i)
      assert e[i]==a[i]:"array elements are equal";
  }
  private static void assertEquals(float[][] e, float[][] a) {
    int n = e.length;
    for (int i=0; i<n; ++i)
      assertEquals(e[i],a[i]);
  }
  private static void assertEquals(float[][][] e, float[][][] a) {
    int n = e.length;
    for (int i=0; i<n; ++i)
      assertEquals(e[i],a[i]);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
