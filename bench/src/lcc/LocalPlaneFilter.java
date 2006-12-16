/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

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

  public LocalPlaneFilter(double sigma) {
    _rgfGradient = new RecursiveGaussianFilter(1.0);
    _rgfSmoother = new RecursiveGaussianFilter(sigma);
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

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    float q1 = p1[i2][i1];
    float q2 = p2[i2][i1];
    float qp = q1+q2;
    float qm = q1-q2;
    float qs = q1*q1;
    y2[i1] = x2[i1]-P999*(0.5f*q2*qp*x2[i1+1]);
    for (i1=1; i1<n1-1; ++i1) {
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      qp = q1+q2;
      qm = q1-q2;
      y2[i1] = x2[i1]-P999*(0.5f*q2*(qp*x2[i1+1]-qm*x2[i1-1]));
    }
    q1 = p1[i2][i1];
    q2 = p2[i2][i1];
    qm = q1-q2;
    y2[i1] = x2[i1]-P999*(-0.5f*q2*qm*x2[i1-1]);

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      qp = q1+q2;
      qs = q1*q1;
      y2[i1] = x2[i1]-P999*(0.5f*q2*qp*x2[i1+1]+qs*x2m[i1]);
      for (i1=1; i1<n1-1; ++i1) {
        q1 = p1[i2][i1];
        q2 = p2[i2][i1];
        qp = q1+q2;
        qm = q1-q2;
        qs = q1*q1;
        y2[i1] = x2[i1]-P999*(0.5f*q2*(qp*x2[i1+1]-qm*x2[i1-1])+qs*x2m[i1]);
      }
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      qm = q1-q2;
      qs = q1*q1;
      y2[i1] = x2[i1]-P999*(-0.5f*q2*qm*x2[i1-1]+qs*x2m[i1]);
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

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
      float q1 = p1[i2][i1];
      float q2 = p2[i2][i1];
      float qp = q1+q2;
      float qm = q1-q2;
      r[i1] = x2[i1];
      ta[i1] =  P999*0.5f*q2*qm;
      tb[i1] = 1.0f;
      tc[i1] = -P999*0.5f*q2*qp;
    }
    tm.solve(r,y2);

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        float q1 = p1[i2][i1];
        float q2 = p2[i2][i1];
        float qp = q1+q2;
        float qm = q1-q2;
        float qs = q1*q1;
        r[i1] = x2[i1]+P999*qs*y2m[i1];
        ta[i1] =  P999*0.5f*q2*qm;
        tb[i1] = 1.0f;
        tc[i1] = -P999*0.5f*q2*qp;
      }
      tm.solve(r,y2);
    }
  }

  public void applyForwardX(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    float q1 = p1[i2][i1];
    float q2 = p2[i2][i1];
    float qp = q1+q2;
    float qm = q1-q2;
    y2[i1] = x2[i1]+0.5f*P999*(qm*qp*(x2[i1+1]));
    for (i1=1; i1<n1-1; ++i1) {
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      qp = q1+q2;
      qm = q1-q2;
      y2[i1] = x2[i1]+0.5f*P999*(qm*qp*(x2[i1+1]+x2[i1-1]));
    }
    q1 = p1[i2][i1];
    q2 = p2[i2][i1];
    qp = q1+q2;
    qm = q1-q2;
    y2[i1] = x2[i1]+0.5f*P999*(qm*qp*(x2[i1-1]));

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      qp = q1+q2;
      qm = q1-q2;
      y2[i1] = x2[i1]+0.5f*P999*(qm*qp*(x2[i1+1]-2.0f*x2m[i1]))
                     -0.5f*P998*(qp*qp*x2m[i1+1]);
      for (i1=1; i1<n1-1; ++i1) {
        q1 = p1[i2][i1];
        q2 = p2[i2][i1];
        qp = q1+q2;
        qm = q1-q2;
        y2[i1] = x2[i1]+0.5f*P999*(qm*qp*(x2[i1+1]+x2[i1-1]-2.0f*x2m[i1]))
                       -0.5f*P998*(qp*qp*x2m[i1+1]+qm*qm*x2m[i1-1]);
      }
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      qp = q1+q2;
      qm = q1-q2;
      y2[i1] = x2[i1]+0.5f*P999*(qm*qp*(x2[i1-1]-2.0f*x2m[i1]))
                     -0.5f*P998*(qm*qm*x2m[i1-1]);
    }
  }

  public void applyInverseX(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    TridiagonalFMatrix tm = new TridiagonalFMatrix(n1);
    float[] ta = tm.a();
    float[] tb = tm.b();
    float[] tc = tm.c();
    float[] r = new float[n1];

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
      float q1 = p1[i2][i1];
      float q2 = p2[i2][i1];
      float qp = q1+q2;
      float qm = q1-q2;
      r[i1] = x2[i1];
      ta[i1] = tc[i1] = 0.5f*P999*qm*qp;
      tb[i1] = 1.0f;
    }
    tm.solve(r,y2);

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        float q1 = p1[i2][i1];
        float q2 = p2[i2][i1];
        float qp = q1+q2;
        float qm = q1-q2;
        r[i1] = x2[i1]+P999*qm*qp*y2m[i1];
        if (i1>0)
          r[i1] += 0.5f*P998*qm*qm*y2m[i1-1];
        if (i1<n1-1)
          r[i1] += 0.5f*P998*qp*qp*y2m[i1+1];
        ta[i1] = tc[i1] = 0.5f*P999*qm*qp;
        tb[i1] = 1.0f;
      }
      tm.solve(r,y2);
    }
  }

  public void applyForwardF(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float tiny = 1.0f-P999;

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    float q1 = max(p1[i2][i1],tiny);
    float q2 = p2[i2][i1];
    float s = -q2/q1;
    float bm = (1.0f-s)*(2.0f-s)/12.0f;
    float b0 = (2.0f-s)*(2.0f+s)/6.0f;
    float bp = (1.0f+s)*(2.0f+s)/12.0f;
    y2[i1] = b0*x2[i1]+bp*x2[i1+1];
    for (i1=1; i1<n1-1; ++i1) {
      q1 = max(p1[i2][i1],tiny);
      q2 = p2[i2][i1];
      s = -q2/q1;
      bm = (1.0f-s)*(2.0f-s)/12.0f;
      b0 = (2.0f-s)*(2.0f+s)/6.0f;
      bp = (1.0f+s)*(2.0f+s)/12.0f;
      y2[i1] = b0*x2[i1]+bm*x2[i1-1]+bp*x2[i1+1];
    }
    q1 = max(p1[i2][i1],tiny);
    q2 = p2[i2][i1];
    s = -q2/q1;
    bm = (1.0f-s)*(2.0f-s)/12.0f;
    b0 = (2.0f-s)*(2.0f+s)/6.0f;
    bp = (1.0f+s)*(2.0f+s)/12.0f;
    y2[i1] = b0*x2[i1]+bm*x2[i1-1];

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      q1 = max(p1[i2][i1],tiny);
      q2 = p2[i2][i1];
      s = -q2/q1;
      bm = (1.0f-s)*(2.0f-s)/12.0f;
      b0 = (2.0f-s)*(2.0f+s)/6.0f;
      bp = (1.0f+s)*(2.0f+s)/12.0f;
      y2[i1] = b0*(x2[i1]-x2m[i1])
              +bm*(        -x2m[i1+1])
              +bp*(x2[i1+1]          );
      for (i1=1; i1<n1-1; ++i1) {
        q1 = max(p1[i2][i1],tiny);
        q2 = p2[i2][i1];
        s = -q2/q1;
        bm = (1.0f-s)*(2.0f-s)/12.0f;
        b0 = (2.0f-s)*(2.0f+s)/6.0f;
        bp = (1.0f+s)*(2.0f+s)/12.0f;
        y2[i1] = b0*(x2[i1]-x2m[i1])
                +bm*(x2[i1-1]-x2m[i1+1])
                +bp*(x2[i1+1]-x2m[i1-1]);
      }
      q1 = max(p1[i2][i1],tiny);
      q2 = p2[i2][i1];
      s = -q2/q1;
      bm = (1.0f-s)*(2.0f-s)/12.0f;
      b0 = (2.0f-s)*(2.0f+s)/6.0f;
      bp = (1.0f+s)*(2.0f+s)/12.0f;
      y2[i1] = b0*(x2[i1]-x2m[i1])
              +bm*(x2[i1-1]          )
              +bp*(        -x2m[i1-1]);
    }
  }

  public void applyInverseF(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    TridiagonalFMatrix tm = new TridiagonalFMatrix(n1);
    float[] ta = tm.a();
    float[] tb = tm.b();
    float[] tc = tm.c();
    float[] r = new float[n1];

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
      float q1 = p1[i2][i1];
      float q2 = p2[i2][i1];
      float qp = q1+q2;
      float qm = q1-q2;
      r[i1] = x2[i1];
      ta[i1] = tc[i1] = 0.5f*P999*qm*qp;
      tb[i1] = 1.0f;
    }
    tm.solve(r,y2);

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        float q1 = p1[i2][i1];
        float q2 = p2[i2][i1];
        float qp = q1+q2;
        float qm = q1-q2;
        r[i1] = x2[i1]+P999*qm*qp*y2m[i1];
        if (i1>0)
          r[i1] += 0.5f*P998*qm*qm*y2m[i1-1];
        if (i1<n1-1)
          r[i1] += 0.5f*P998*qp*qp*y2m[i1+1];
        ta[i1] = tc[i1] = 0.5f*P999*qm*qp;
        tb[i1] = 1.0f;
      }
      tm.solve(r,y2);
    }
  }

  public void applyForwardS(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    float q1 = p1[i2][i1];
    float q2 = p2[i2][i1];
    y2[i1] = q1*x2[i1];
    if (q2<0.0f) {
      y2[i1] -= q2*x2[i1];
    } else {
      y2[i1] += q2*(x2[i1]-P999*x2[i1+1]);
    }
    for (i1=1; i1<n1-1; ++i1) {
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      y2[i1] = q1*x2[i1];
      if (q2<0.0f) {
        y2[i1] -= q2*(x2[i1]-P999*x2[i1-1]);
      } else {
        y2[i1] += q2*(x2[i1]-P999*x2[i1+1]);
      }
    }
    q1 = p1[i2][i1];
    q2 = p2[i2][i1];
    y2[i1] = q1*x2[i1];
    if (q2<0.0f) {
      y2[i1] -= q2*(x2[i1]-P999*x2[i1-1]);
    } else {
      y2[i1] += q2*x2[i1];
    }

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      y2[i1] = q1*(x2[i1]-P999*x2m[i1]);
      if (q2<0.0f) {
        y2[i1] -= q2*x2[i1];
      } else {
        y2[i1] += q2*(x2[i1]-P999*x2[i1+1]);
      }
      for (i1=1; i1<n1-1; ++i1) {
        q1 = p1[i2][i1];
        q2 = p2[i2][i1];
        y2[i1] = q1*(x2[i1]-P999*x2m[i1]);
        if (q2<0.0f) {
          y2[i1] -= q2*(x2[i1]-P999*x2[i1-1]);
        } else {
          y2[i1] += q2*(x2[i1]-P999*x2[i1+1]);
        }
      }
      q1 = p1[i2][i1];
      q2 = p2[i2][i1];
      y2[i1] = q1*(x2[i1]-P999*x2m[i1]);
      if (q2<0.0f) {
        y2[i1] -= q2*(x2[i1]-P999*x2[i1-1]);
      } else {
        y2[i1] += q2*x2[i1];
      }
    }
  }

  public void applyInverseS(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    TridiagonalFMatrix tm = new TridiagonalFMatrix(n1);
    float[] ta = tm.a();
    float[] tb = tm.b();
    float[] tc = tm.c();
    float[] r = new float[n1];

    // For i2=0, assume that x[i2-1][i1] = y[i2-1][i1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
      float q1 = p1[i2][i1];
      float q2 = p2[i2][i1];
      r[i1] = x2[i1];
      if (q2<0.0f) {
        ta[i1] = P999*q2;
        tb[i1] = q1-q2;
        tc[i1] = 0.0f;
      } else {
        ta[i1] = 0.0f;
        tb[i1] = q1+q2;
        tc[i1] = -P999*q2;
      }
    }
    tm.solve(r,y2);

    // For all i2>0, ...
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        float q1 = p1[i2][i1];
        float q2 = p2[i2][i1];
        r[i1] = x2[i1]+P999*q1*y2m[i1];
        if (q2<0.0f) {
          ta[i1] = P999*q2;
          tb[i1] = q1-q2;
          tc[i1] = 0.0f;
        } else {
          ta[i1] = 0.0f;
          tb[i1] = q1+q2;
          tc[i1] = -P999*q2;
        }
      }
      tm.solve(r,y2);
    }
  }
  /*
  q2 from - to + (unstable?)
  0  -1   1   0   0   0   0   0
  0   0  -1   1   0   0   0   0
  0   0   0   0   1  -1   0   0
  0   0   0   0   0   1  -1   0

  q2 from + to - (stable?)
  0   0   1  -1   0   0   0   0
  0   0   0   1  -1   0   0   0
  0   0   0  -1   1   0   0   0
  0   0   0   0  -1   1   0   0
  */

  ///////////////////////////////////////////////////////////////////////////
  // private

  private RecursiveGaussianFilter _rgfGradient;
  private RecursiveGaussianFilter _rgfSmoother;
  private static final float P999 = 0.99999f;
  private static final float P998 = P999*P999;
}
