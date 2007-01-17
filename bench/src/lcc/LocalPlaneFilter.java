/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import java.util.Random;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local plane prediction error filtering of images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.12.12
 */
public class LocalPlaneFilter {

  public enum Type {
    QUAD,
    HALE0,
    HALE1,
    HALE2,
    HALE3,
    FOMEL1,
    FOMEL2,
  };

  public LocalPlaneFilter(double sigma) {
    this(sigma,Type.HALE1);
  }

  public LocalPlaneFilter(double sigma, Type type) {
    _rgfGradient = new RecursiveGaussianFilter(1.0);
    _rgfSmoother = new RecursiveGaussianFilter(sigma);
    if (type==Type.QUAD) {
      _filter = new QuadFilter();
    } else if (type==Type.HALE0) {
      _filter = new Hale0Filter();
    } else if (type==Type.HALE1) {
      _filter = new Hale1Filter();
    } else if (type==Type.HALE2) {
      _filter = new Hale2Filter();
    } else if (type==Type.FOMEL1) {
      _filter = new Fomel1Filter();
    } else if (type==Type.FOMEL2) {
      _filter = new Fomel2Filter();
    }
  }

  public float[][][] find(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;

    // Gradient.
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    _rgfGradient.apply1X(x,g1);
    _rgfGradient.applyX1(x,g2);

    // Gradient products.
    float[][] g11 = new float[n2][n1];
    float[][] g12 = g1;
    float[][] g22 = g2;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        g11[i2][i1] = g1i*g1i;
        g12[i2][i1] = g1i*g2i;
        g22[i2][i1] = g2i*g2i;
      }
    }
    
    // Smoothed gradient products comprise the structure tensor.
    float[][] gtt = new float[n2][n1];
    _rgfSmoother.apply0X(g11,gtt);
    _rgfSmoother.applyX0(gtt,g11);
    _rgfSmoother.apply0X(g12,gtt);
    _rgfSmoother.applyX0(gtt,g12);
    _rgfSmoother.apply0X(g22,gtt);
    _rgfSmoother.applyX0(gtt,g22);

    // For each sample, the eigenvector corresponding to the largest
    // eigenvalue is normal to the plane. The size of that eigenvalue
    // is a measure of planarity in the interval [0,1].
    float[][][] p = new float[3][][];
    p[0] = g11;
    p[1] = g12;
    p[2] = g22;
    float[][] a = new float[2][2];
    float[][] v = new float[2][2];
    float[] d = new float[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = g11[i2][i1];
        a[0][1] = g12[i2][i1];
        a[1][0] = g12[i2][i1];
        a[1][1] = g22[i2][i1];
        Eigen.solveSymmetric22(a,v,d);
        p[0][i2][i1] = d[0]/sqrt(d[0]*d[0]+d[1]*d[1]);
        float v1 = v[0][0];
        float v2 = v[0][1];
        if (v1<0.0f) {
          v1 = -v1;
          v2 = -v2;
        }
        p[1][i2][i1] = v1;
        p[2][i2][i1] = v2;
      }
    }
    return p;
  }

  public void applyForward(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float[] c = new float[6];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1];
    for (i1=1; i1<n1-1; ++i1) {
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1]+c[2]*x2[i1-1];
    }
    _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[1]*x2[i1]+c[2]*x2[i1-1];

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]
             + c[3]*x2m[i1+1]+c[4]*x2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
        y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]+c[2]*x2 [i1-1]
               + c[3]*x2m[i1+1]+c[4]*x2m[i1]+c[5]*x2m[i1-1];
      }
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[1]*x2 [i1]+c[2]*x2 [i1-1]
             + c[4]*x2m[i1]+c[5]*x2m[i1-1];
    }
  }

  public void applyInverse(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    TridiagonalFMatrix tm = new TridiagonalFMatrix(n1);
    float[] ta = tm.a();
    float[] tb = tm.b();
    float[] tc = tm.c();
    float[] r = new float[n1];
    float[] c = new float[6];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      r[i1] = x2[i1];
      ta[i1] = c[2];
      tb[i1] = c[1];
      tc[i1] = c[0];
    }
    tm.solve(r,y2);

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      int i1 = 0;
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      tb[i1] = c[1];
      tc[i1] = c[0];
      r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
        ta[i1] = c[2];
        tb[i1] = c[1];
        tc[i1] = c[0];
        r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      }
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      ta[i1] = c[2];
      tb[i1] = c[1];
      r[i1] = x2[i1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      tm.solve(r,y2);
    }
  }

  public void applyForwardF(float[][][] p, float[][] x, float[][] y) {
    Hale3Filter filter = new Hale3Filter();
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float[] c = new float[9];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2m = null;
    float[] x2p = x[i2+1];
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    /*
    filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[0]*x2p[i1+1]+c[1]*x2p[i1]
           + c[3]* x2[i1+1]+c[4]* x2[i1];
    for (i1=1; i1<n1-1; ++i1) {
      filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[0]*x2p[i1+1]+c[1]*x2p[i1]+c[2]*x2p[i1-1]
             + c[3]* x2[i1+1]+c[4]* x2[i1]+c[5]* x2[i1-1];
    }
    filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[1]*x2p[i1]+c[2]*x2p[i1-1]
           + c[4]* x2[i1]+c[5]* x2[i1-1];
    */
    for (i1=0; i1<n1; ++i1)
      y2[i1] = 0.0f;

    // For all 0<i2<n2-1, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2-1; ++i2) {
      x2m = x[i2-1];
      x2p = x[i2+1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      /*
      filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[0]*x2p[i1+1]+c[1]*x2p[i1]
             + c[3]* x2[i1+1]+c[4]* x2[i1]
             + c[6]*x2m[i1+1]+c[7]*x2m[i1];
      */
      y2[0] = 0.0f;
      for (i1=1; i1<n1-1; ++i1) {
        filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
        y2[i1] = c[0]*x2p[i1+1]+c[1]*x2p[i1]+c[2]*x2p[i1-1]
               + c[3]* x2[i1+1]+c[4]* x2[i1]+c[5]* x2[i1-1]
               + c[6]*x2m[i1+1]+c[7]*x2m[i1]+c[8]*x2m[i1-1];
      }
      /*
      filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[1]*x2p[i1]+c[2]*x2p[i1-1]
             + c[4]* x2[i1]+c[5]* x2[i1-1]
             + c[7]*x2m[i1]+c[8]*x2m[i1-1];
      */
      y2[n1-1] = 0.0f;
    }

    // For i2=n2-1, x[i2+1][i1] = y[i2+1][i1] = x[i2][-1] = x[i2][n1] = 0.
    x2m = x[i2-1];
    x2p = null;
    x2 = x[i2];
    y2 = y[i2];
    /*
    i1 = 0;
    filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[3]* x2[i1+1]+c[4]* x2[i1]
           + c[6]*x2m[i1+1]+c[7]*x2m[i1];
    for (i1=1; i1<n1-1; ++i1) {
      filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[3]* x2[i1+1]+c[4]* x2[i1]+c[5]* x2[i1-1]
             + c[6]*x2m[i1+1]+c[7]*x2m[i1]+c[8]*x2m[i1-1];
    }
    filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[4]* x2[i1]+c[5]* x2[i1-1]
           + c[7]*x2m[i1]+c[8]*x2m[i1-1];
    */
    for (i1=0; i1<n1; ++i1)
      y2[i1] = 0.0f;
  }

  public void applyForwardX(float[][][] p, float[][] x, float[][] y) {
    if (_lcf==null)
      makeLocalCausalFilter();
    A2 a2 = new A2(index(p));
    _lcf.apply(a2,x,y);
    _lcf.applyTranspose(a2,y,y);
  }
  public void applyInverseX(float[][][] p, float[][] x, float[][] y) {
    if (_lcf==null)
      makeLocalCausalFilter();
    A2 a2 = new A2(index(p));
    _lcf.applyInverseTranspose(a2,x,y);
    _lcf.applyInverse(a2,y,y);
  }
  private static int NTHETA = 65;
  private static float FTHETA = -0.5f*FLT_PI;
  private static float DTHETA = FLT_PI/(float)(NTHETA-1);;
  private static float STHETA = 0.9999f/DTHETA;
  private LocalCausalFilter _lcf;
  private static int[][] xindex(float[][][] p) {
    Random r = new Random(314159);
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    int[][] i = new int[n2][n1];
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float theta = asin(p2[i2][i1]);
        i[i2][i1] = (int)(0.5f+(theta-FTHETA)*STHETA);
        /*
        float thetai = (theta-FTHETA)*STHETA;
        int itheta = (int)(thetai);
        float dtheta = thetai-(float)itheta;
        if (r.nextFloat()<dtheta) ++itheta;
        i[i2][i1] = itheta;
        */
      }
    }
    return i;
  }
  private static int NP2 = 65;
  private static float FP2 = -1.0f;
  private static float DP2 = 2.0f/(float)(NP2-1);;
  private static float SP2 = 0.9999f/DP2;
  private static int[][] index(float[][][] p) {
    Random r = new Random(314159);
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    int[][] i = new int[n2][n1];
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float p2i = p2[i2][i1];
        i[i2][i1] = (int)(0.5f+(p2i-FP2)*SP2);
        /*
        float thetai = (theta-FTHETA)*STHETA;
        int itheta = (int)(thetai);
        float dtheta = thetai-(float)itheta;
        if (r.nextFloat()<dtheta) ++itheta;
        i[i2][i1] = itheta;
        */
      }
    }
    return i;
  }
  private void xmakeLocalCausalFilter() {
    int maxlag = 6;
    int nlag = maxlag+2+maxlag;
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    for (int ilag=0; ilag<nlag; ++ilag) {
      lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
      lag2[ilag] = (ilag<=maxlag)?0:1;
    }
    for (int itheta=0; itheta<NTHETA; ++itheta) {
      float theta = FTHETA+itheta*DTHETA;
      float c = cos(theta);
      float s = sin(theta);
      float m12 = 0.5f*(c-s);
      float p12 = 0.5f*(c+s);
      float[][] r = {
        {    -m12*m12,   -2.0f*m12*p12,      -p12*p12},
        {2.0f*m12*p12,          1.001f,  2.0f*m12*p12},
        {    -p12*p12,   -2.0f*m12*p12,      -m12*m12}
      };
      CausalFilter cf = new CausalFilter(lag1,lag2);
      cf.factorWilsonBurg(100,0.000001f,r);
      _aTable[itheta] = cf.getA();
    }
    _lcf = new LocalCausalFilter(lag1,lag2);
  }
  private void makeLocalCausalFilter() {
    int maxlag = 6;
    int nlag = maxlag+2+maxlag;
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    for (int ilag=0; ilag<nlag; ++ilag) {
      lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
      lag2[ilag] = (ilag<=maxlag)?0:1;
    }
    for (int ip2=0; ip2<NTHETA; ++ip2) {
      float p2 = FP2+ip2*DP2;
      float p1 = sqrt(1.0f-p2*p2);
      float m12 = 0.5f*(p1-p2);
      float p12 = 0.5f*(p1+p2);
      float[][] r = {
        {    -m12*m12,   -2.0f*m12*p12,      -p12*p12},
        {2.0f*m12*p12,          1.001f,  2.0f*m12*p12},
        {    -p12*p12,   -2.0f*m12*p12,      -m12*m12}
      };
      CausalFilter cf = new CausalFilter(lag1,lag2);
      cf.factorWilsonBurg(100,0.000001f,r);
      _aTable[ip2] = cf.getA();
    }
    _lcf = new LocalCausalFilter(lag1,lag2);
  }
  //private float[][] _aTable = new float[NTHETA][];
  private float[][] _aTable = new float[NP2][];
  private class A2 implements LocalCausalFilter.A2 {
    A2(int[][] index) {
      _index = index;
    }
    public void get(int i1, int i2, float[] a) {
      Array.copy(_aTable[_index[i2][i1]],a);
    }
    private int[][] _index;
  }

  public void xapplyForwardX(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float p1i,p2i,pm,pp,pmy,ppy;
    Array.copy(x,y);

    // Apply plane filter.
    for (int i2=n2-1; i2>0; --i2) {
      float[] y2 = y[i2];
      float[] y2m = y[i2-1];
      for (int i1=n1-1; i1>0; --i1) {
        p1i = p1[i2][i1]; p2i = p2[i2][i1]; 
        pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
        y2[i1] = pm*(y2[i1]-y2m[i1-1])+pp*(y2[i1-1]-y2m[i1]);
      }
      p1i = p1[i2][0]; p2i = p2[i2][0]; 
      pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
      y2[0] = pm*y2[0]-pp*y2m[0];
    }
    for (int i1=n1-1; i1>0; --i1) {
      p1i = p1[0][i1]; p2i = p2[0][i1]; 
      pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
      y[0][i1] = pm*y[0][i1]+pp*y[0][i1-1];
    }
    p1i = p1[0][0]; p2i = p2[0][0]; 
    pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
    y[0][0] = pm*y[0][0];

    // Apply transpose of plane filter (with care for variable coefficients).
    p1i = p1[0][0]; p2i = p2[0][0]; 
    pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
    pmy = pm*y[0][0]; ppy = pp*y[0][0];
    y[0][0] = pmy;
    for (int i1=1; i1<n1; ++i1) {
      p1i = p1[0][i1]; p2i = p2[0][i1]; 
      pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
      pmy = pm*y[0][i1]; ppy = pp*y[0][i1];
      y[0][i1] = pmy;
      y[0][i1-1] += ppy;
    }
    for (int i2=1; i2<n2; ++i2) {
      float[] y2 = y[i2];
      float[] y2m = y[i2-1];
      p1i = p1[i2][0]; p2i = p2[i2][0]; 
      pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
      pmy = pm*y2[0]; ppy = pp*y2[0];
      y2[0] = pmy;
      y2m[0] -= ppy;
      for (int i1=1; i1<n1; ++i1) {
        p1i = p1[i2][i1]; p2i = p2[i2][i1]; 
        pm = 0.5f*(p1i-p2i); pp = 0.5f*(p1i+p2i);
        pmy = pm*y2[i1]; ppy = pp*y2[i1];
        y2[i1] = pmy;
        y2[i1-1] += ppy;
        y2m[i1] -= ppy;
        y2m[i1-1] -= pmy;
      }
    }

    // Stability.
    float ax = 0.001f;
    float ay = 1.0f-ax;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        y2[i1] = ay*y2[i1]+ax*x2[i1];
    }
  }

  // Conjugate-gradients without pre-conditioning.
  public void xxapplyInverseX(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    Array.zero(y);
    Array.copy(x,r);
    Array.copy(r,s);
    float rr = dot(r,r);
    float small = rr*0.00001f;
    System.out.println("small="+small);
    int niter;
    for (niter=0; niter<200 && rr>small; ++niter) {
      xapplyForwardX(p,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      float rrold = rr;
      rr = dot(r,r);
      float beta = rr/rrold;
      for (int i2=0; i2<n2; ++i2) {
        float[] r2 = r[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = r2[i1]+beta*s2[i1];
      }
      System.out.println("niter="+niter+" rr="+rr);
    }
  }
  // Conjugate-gradients with pre-conditioning.
  public void xapplyInverseX(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    float[][] u = new float[n2][n1]; // s
    Array.zero(y);
    Array.copy(x,r);
    applyInverseX(p,r,s);
    float rr = dot(r,s);
    float small = rr*0.00001f;
    System.out.println("small="+small);
    int niter;
    for (niter=0; niter<200 && rr>small; ++niter) {
      xapplyForwardX(p,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      applyInverseX(p,r,u);
      float rrold = rr;
      rr = dot(r,u);
      float beta = rr/rrold;
      for (int i2=0; i2<n2; ++i2) {
        float[] u2 = u[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = u2[i1]+beta*s2[i1];
      }
      System.out.println("niter="+niter+" rr="+rr);
    }
  }
  private static float dot(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    double s = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        s += x2[i1]*y2[i1];
    }
    return (float)s;
  }
  private static void saxpy(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        y2[i1] += a*x2[i1];
    }
  }

  public void smooth(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int i2 = 0;
    for (i2=0; i2<n2-1; ++i2) {
      float[] yn = y[i2];
      float[] xn = x[i2];
      float[] xp = x[i2+1];
      int i1 = 0;
      float xnm = 0.0f;
      float xpm = 0.0f;
      float xnn = xn[i1];
      float xpn = xp[i1];
      float xnp = xn[i1+1];
      float xpp = xp[i1+1];
      yn[i1] = 0.25f*(xnn+xpn+xnp+xpp);
      for (i1=1; i1<n1-1; ++i1) {
        xnm = xnn;
        xpm = xpn;
        xnn = xnp;
        xpn = xpp;
        xnp = xn[i1+1];
        xpp = xp[i1+1];
        yn[i1] = 0.125f*(xnm+xpm+2.0f*(xnn+xpn)+xnp+xpp);
      }
      xnm = xnn;
      xpm = xpn;
      xnn = xnp;
      xpn = xpp;
      yn[i1] = 0.25f*(xnm+xpm+xnn+xpn);
    }
    float[] yn = y[i2];
    float[] xn = x[i2];
    int i1 = 0;
    float xnm = 0.0f;
    float xnn = xn[i1];
    float xnp = xn[i1+1];
    yn[i1] = 0.5f*(xnn+xnp);
    for (i1=1; i1<n1-1; ++i1) {
      xnm = xnn;
      xnn = xnp;
      xnp = xn[i1+1];
      yn[i1] = 0.25f*(xnm+2.0f*xnn+xnp);
    }
    xnm = xnn;
    xnn = xnp;
    yn[i1] = 0.5f*(xnm+xnn);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float P99 = 1.0000f;
  //private static final float P98 = P99*P99;
  private static final float P49 = P99/2.0f;
  //private static final float P48 = P98/2.0f;

  private RecursiveGaussianFilter _rgfGradient;
  private RecursiveGaussianFilter _rgfSmoother;
  private Filter _filter;

  /**
   * Different plane filters have different coefficients.
   * For forward filtering, the coefficients are used as follows:
   * <pre><code>
   * y[i2][i1] = c[0]*x[i2  ][i1+1]+c[1]*x[i2  ][i1]+c[2]*x[i2  ][i1-1]
   *           + c[3]*x[i2-1][i1+1]+c[4]*x[i2-1][i1]+c[5]*x[i2-1][i1-1]
   * </code></pre>
   * For inverse filtering, a corresponding tridiagonal system of 
   * equations is solved. Note that some filters may not be invertible. 
   * The corresponding finite-difference stencil is
   * <pre><code>
   *        i2-1  i2
   *   i1-1 c[5] c[2]
   *   i1   c[4] c[1]
   *   i1+1 c[3] c[0]
   * </code></pre>
   */
  private interface Filter {
    /**
     * Gets coefficients for specified normal unit-vector.
     * @param n1 component of normal for 1st dimension.
     * @param n2 component of normal for 2nd dimension.
     * @param c array of computed coefficients.
     */
    public void getCoefficients(float n1, float n2, float[] c);
  }
  private static class QuadFilter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float np = n1+n2;
      float nm = n1-n2;
      c[0] = 0.0f;
      c[1] = nm;
      c[2] = np;
      c[3] = 0.0f;
      c[4] = -np;
      c[5] = -nm;
    }
  }
  private static class Hale0Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      if (n2<0.0f) {
        c[0] = 0.0f;
        c[1] = n1-n2;
        c[2] = P99*n2;
      } else {
        c[0] = -P99*n2;
        c[1] = n1+n2;
        c[2] = 0.0f;
      }
      c[3] = 0.0f;
      c[4] = -P99*n1;
      c[5] = 0.0f;
    }
  }
  private static class Hale1Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      c[0] = -P49*n2*(n1+n2);
      c[1] =  1.0f;
      c[2] =  P49*n2*(n1-n2);
      c[3] =  0.0f;
      c[4] = -P99*n1*n1;
      c[5] =  0.0f;
    }
  }
  private static class Hale2Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float np = n1+n2;
      float nm = n1-n2;
      c[0] = 0.25f*nm*np;
      c[1] = 0.50f;
      c[2] = c[0];
      c[3] = -0.25f*np*np;
      c[4] = -0.50f*nm*np;
      c[5] = -0.25f*nm*nm;
    }
  }
  private static class Hale3Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float np = 0.5f*(n1+n2);
      float nm = 0.5f*(n1-n2);
      c[0] =      -nm*nm;
      c[1] = -2.0f*nm*np;
      c[2] =      -np*np;
      c[3] =  2.0f*nm*np;
      c[4] = 1.001f;
      c[5] = c[3];
      c[6] = c[2];
      c[7] = c[1];
      c[8] = c[0];
    }
  }
  private static class Fomel1Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      if (n1<N1MIN)
        n1 = N1MIN;
      float s = n2/n1;
      c[0] = (1.0f-s)*(2.0f-s)/12.0f;
      c[1] = (2.0f+s)*(2.0f-s)/6.0f;
      c[2] = (1.0f+s)*(2.0f+s)/12.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
    private static final float N1MIN = 0.001f;
  }
  private static class Fomel2Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float t1 = 2.0f*n1;
      c[0] = (n1-n2)*(t1-n2)/12.0f;
      c[1] = (t1+n2)*(t1-n2)/6.0f;
      c[2] = (n1+n2)*(t1+n2)/12.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
  }
}
