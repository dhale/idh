/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.logging.Logger;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates and applies shifts to flatten features in 2D and 3D images.
 * In 2D, the shifts (in samples) are functions s(x1,x2) that flatten
 * image features. Specifically, for a 2D input image f(x1,x2), flattening 
 * is performed by computing g(x1,x2) = f(x1-s(x1,x2),x2). The shifts 
 * s(x1,x2) added to x1 in this expression are computed from slopes
 * measured in the input image f(x1,x2) so that features in the output
 * image g(x1,x2) vary as little as possible with x2.
 *
 * Flattens features in 3D images in the same way, but using shifts
 * s(x1,x2,x3) to compute g(x1,x2,x3) = f(x1-s(x1,x2,x3),x2,x3).
 * Features in the flattened image g(x1,x2,x3 vary little with x2 and x3.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.08
 */
public class FlattenerLsqr {

  public FlattenerLsqr(double sigma, double epsilon) {
    _sigma = (float)sigma;
    _epsilon = (float)epsilon;
  }

  public float[][] findShifts(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] p2 = new float[n2][n1]; // slopes p2
    float[][] el = new float[n2][n1]; // linearities
    LocalSlopeFinder lsf = new LocalSlopeFinder(_sigma,10.0);
    lsf.findSlopes(f,p2,el);
    return findShifts(p2,el);
  }

  public float[][] findShifts(float[][] p2, float[][] el) {
    int n1 = p2[0].length;
    int n2 = p2.length;
    float[][] r = new float[n2][2*n1]; // right-hand side
    float[][] s = new float[n2][n1]; // the shifts
    makeRhs(_epsilon,p2,el,r);
    double atol = 0.00;
    double btol = 0.01;
    double ctol = 0.1*sqrt(FLT_EPSILON);
    int maxi = _niter;
    double damp = 0.000;
    LsqrSolver ls = new LsqrSolver(atol,btol,ctol,maxi);
    A2 a2 = new A2(_epsilon,p2,el);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    VecArrayFloat2 vs = new VecArrayFloat2(s);
    LsqrSolver.Info info = ls.solve(a2,damp,vr,vs);
    invertShifts(s);
    return s;
  }

  public float[][][] findShifts(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] p2 = new float[n3][n2][n1]; // slopes p2
    float[][][] p3 = new float[n3][n2][n1]; // slopes p3
    float[][][] ep = new float[n3][n2][n1]; // planarities
    LocalSlopeFinder lsf = new LocalSlopeFinder(_sigma,10.0);
    lsf.findSlopes(f,p2,p3,ep);
    return findShifts(p2,p3,ep);
  }

  public float[][][] findShifts(
    float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    int n1 = p2[0][0].length;
    int n2 = p2[0].length;
    int n3 = p2.length;
    float[][][] r = new float[n3][n2][n1]; // right-hand side
    float[][][] s = new float[n3][n2][n1]; // the shifts
    makeRhs(_epsilon,p2,p3,ep,r);
    A3 a = new A3(_epsilon,p2,p3,ep);
    CgLinearSolver cls = new CgLinearSolver(1000,_small);
    cls.solve(a,r,s);
    //invertShifts(s);
    return s;
  }

  public float[][] applyShifts(float[][] f, float[][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      sub(r,s[i2],t);
      si.setUniformSamples(f[i2]);
      si.interpolate(n1,t,g[i2]);
    }
    return g;
  }

  public float[][][] applyShifts(float[][][] f, float[][][] s) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][][] g = zerofloat(n1,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        sub(r,s[i3][i2],t);
        si.setUniformSamples(f[i3][i2]);
        si.interpolate(n1,t,g[i3][i2]);
      }
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean PARALLEL = true; // false for single-threaded

  private float _sigma = 8.0f; // smoothing half-width to estimate slopes
  private float _epsilon = 0.1f; // penalty for ds/dx1 term in least-squares
  private float _small = 0.01f; // stop iterations when residuals are small
  private int _niter = 10000; // maximum number of iterations

  private static class A2 implements LsqrSolver.A {
    A2(float epsilon, float[][] p2, float[][] el) {
      _epsilon = epsilon;
      _p2 = p2;
      _el = el;
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      applyInternal(_epsilon,_p2,_el,x,y);
    }
    public void applyTranspose(Vec vy, Vec vx) {
      float[][] y = ((VecArrayFloat2)vy).getArray();
      float[][] x = ((VecArrayFloat2)vx).getArray();
      applyTransposeInternal(_epsilon,_p2,_el,y,x);
    }
    private float _epsilon;
    private float[][] _p2;
    private float[][] _el;
  }

  private static class A3 implements CgLinearSolver.A3 {
    A3(float epsilon, float[][][] p2, float[][][] p3, float[][][] ep) {
      _epsilon = epsilon;
      _p2 = p2;
      _p3 = p3;
      _ep = ep;
    }
    public void apply(float[][][] x, float[][][] y) {
      applyLhs(_epsilon,_p2,_p3,_ep,x,y);
    }
    private float _epsilon;
    private float[][][] _p2;
    private float[][][] _p3;
    private float[][][] _ep;
  }

  private static void cleanShifts(float[] s) {
    int n1 = s.length;
    for (int i1=1; i1<n1; ++i1) {
      if (s[i1]<=s[i1-1]-1.00f)
        s[i1] = s[i1-1]-0.99f;
    }
  }

  private static void invertShifts(
    InverseInterpolator ii, float[] u, float[] t, float[] s) 
  {
    cleanShifts(s);
    int n1 = s.length;
    for (int i1=0; i1<n1; ++i1)
      s[i1] += u[i1];
    ii.invert(s,t);
    float tmin = -5.0f;
    float tmax = n1-1+5.0f;
    for (int i1=0; i1<n1; ++i1) {
      if (t[i1]<tmin) t[i1] = tmin;
      if (t[i1]>tmax) t[i1] = tmax;
      s[i1] = u[i1]-t[i1];
    }
  }

  private static void invertShifts(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    float[] u = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    InverseInterpolator ii = new InverseInterpolator(n1,n1);
    for (int i2=0; i2<n2; ++i2)
      invertShifts(ii,u,t,s[i2]);
  }

  private static void invertShifts(float[][][] s) {
    int n1 = s[0][0].length;
    int n2 = s[0].length;
    int n3 = s.length;
    float[] u = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    InverseInterpolator ii = new InverseInterpolator(n1,n1);
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        invertShifts(ii,u,t,s[i3][i2]);
  }

  private static void makeRhs(
    float epsilon, float[][] p2, float[][] el, float[][] r)
  {
    int n1 = p2[0].length;
    int n2 = p2.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0,j1=0; i1<n1; ++i1,j1+=2) {
        float eli = el[i2][i1];
        float p2i = p2[i2][i1];
        float e2i = eli*p2i;
        r[i2][j1  ] = 0.0f;
        r[i2][j1+1] = -e2i;
      }
    }
  }
  private static void applyInternal(
    float eps, float[][] p2, float[][] el, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1,j1=2; i1<n1; ++i1,j1+=2) {
        float eli = el[i2][i1];
        float p2i = p2[i2][i1];
        float e2i = eli*p2i;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.5f*(xa-xb);
        float x2 = 0.5f*(xa+xb);
        float y1 = eps*x1;
        float y2 = e2i*x1+eli*x2;
        y[i2][j1  ] += y1;
        y[i2][j1+1] += y2;
      }
    }
  }
  private static void applyTransposeInternal(
    float eps, float[][] p2, float[][] el, float[][] y, float[][] x) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1,j1=2; i1<n1; ++i1,j1+=2) {
        float eli = el[i2][i1];
        float p2i = p2[i2][i1];
        float e2i = eli*p2i;
        float y1 = y[i2][j1  ];
        float y2 = y[i2][j1+1];
        float x1 = eps*y1+e2i*y2;
        float x2 =        eli*y2;
        float xa = 0.5f*(x1+x2);
        float xb = 0.5f*(x1-x2);
        x[i2  ][i1  ] += xa;
        x[i2  ][i1-1] -= xb;
        x[i2-1][i1  ] += xb;
        x[i2-1][i1-1] -= xa;
      }
    }
  }
  private static void testApplyInternal() {
    //int n1 = 251;
    //int n2 = 357;
    int n1 = 5;
    int n2 = 6;
    float[][] u = sub(randfloat(2*n1,n2),0.5f);
    float[][] w = sub(randfloat(n1,n2),0.5f);
    float[][] x = sub(randfloat(n1,n2),0.5f);
    float[][] p2 = sub(randfloat(n1,n2),0.5f);
    float[][] el = randfloat(n1,n2);
    float eps = 0.01f;
    float[][] y = new float[n2][2*n1];
    float[][] z = new float[n2][n1];
    zero(y); applyInternal(eps,p2,el,x,y);
    zero(z); applyTransposeInternal(eps,p2,el,y,z);
    float pd = sum(mul(x,z));
    System.out.println("positive-semidefinite xax="+sum(mul(x,z)));
    float wax = sum(mul(w,z));
    zero(y); applyInternal(eps,p2,el,w,y);
    zero(z); applyTransposeInternal(eps,p2,el,y,z);
    float xaw = sum(mul(x,z));
    System.out.println("symmetric w'A'Ax="+wax+" x'A'Aw="+xaw);
    zero(y); applyInternal(eps,p2,el,x,y);
    float uax = sum(mul(u,y));
    zero(w); applyTransposeInternal(eps,p2,el,u,w);
    float xau = sum(mul(x,w));
    System.out.println("adjoint u'(Ax)="+uax+" x'(A'u)="+xau);
  }
  private static void testSolver() {
    int n1 = 11;
    int n2 = 9;
    float eps = 0.01f;
    float[][] p2 = sub(randfloat(n1,n2),0.5f);
    float[][] el = fillfloat(1.0f,n1,n2);
    float[][] r = new float[n2][2*n1]; // right-hand side
    float[][] s = new float[n2][n1]; // the shifts
    float[][] t = sub(randfloat(n1,n2),0.5f);
    double atol = 0.00;
    double btol = 0.00;
    double ctol = 0.00;
    int maxi = 10000;
    double damp = 0.00;
    LsqrSolver ls = new LsqrSolver(atol,btol,ctol,maxi);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    VecArrayFloat2 vs = new VecArrayFloat2(s);
    VecArrayFloat2 vt = new VecArrayFloat2(t);
    A2 a2 = new A2(eps,p2,el);
    a2.apply(vt,vr);
    LsqrSolver.Info info = ls.solve(a2,damp,vr,vs);
    plot(s);
    plot(t);
  }
  private static void plot(float[][] x) {
    edu.mines.jtk.mosaic.SimplePlot sp =
      new edu.mines.jtk.mosaic.SimplePlot();
    edu.mines.jtk.mosaic.PixelsView pv = sp.addPixels(x);
    //pv.setClips(-0.5f,0.5f);
    sp.addColorBar();
  }
  public static void main(String[] args) {
    //testApplyInternal();
    testSolver();
  }

  private static void makeRhs(
    float epsilon, 
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    szero(y);
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float epi = ep[i3][i2][i1];
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float b12 = p2i*epi;
          float b13 = p3i*epi;
          float b22 = epi;
          float b33 = epi;
          // float x1 = 0.0f;
          float x2 = -0.25f*p2i;
          float x3 = -0.25f*p3i;
          float y1 = b12*x2+b13*x3;
          float y2 = b22*x2;
          float y3 = b33*x3;
          float ya = y1+y2+y3;
          float yb = y1-y2+y3;
          float yc = y1+y2-y3;
          float yd = y1-y2-y3;
          y[i3  ][i2  ][i1  ] += ya;
          y[i3  ][i2  ][i1-1] -= yd;
          y[i3  ][i2-1][i1  ] += yb;
          y[i3  ][i2-1][i1-1] -= yc;
          y[i3-1][i2  ][i1  ] += yc;
          y[i3-1][i2  ][i1-1] -= yb;
          y[i3-1][i2-1][i1  ] += yd;
          y[i3-1][i2-1][i1-1] -= ya;
        }
      }
    }
  }

  private static void applyLhs(
    float epsilon, 
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    szero(y);
    if (PARALLEL) {
      applyLhsParallel(epsilon,p2,p3,ep,x,y);
    } else {
      applyLhsSerial(epsilon,p2,p3,ep,x,y);
    }
  }

  private static void applyLhsSerial(
    float epsilon, 
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3)
      applyLhsSlice3(i3,epsilon,p2,p3,ep,x,y);
  }

  private static void applyLhsParallel(
    final float epsilon, 
    final float[][][] p2, final float[][][] p3, final float[][][] ep, 
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;

    // i3 = 1, 3, 5, ...
    final AtomicInteger a1 = new AtomicInteger(1);
    Thread[] thread1 = Threads.makeArray();
    for (int ithread=0; ithread<thread1.length; ++ithread) {
      thread1[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a1.getAndAdd(2); i3<n3; i3=a1.getAndAdd(2))
            applyLhsSlice3(i3,epsilon,p2,p3,ep,x,y);
        }
      });
    }
    Threads.startAndJoin(thread1);

    // i3 = 2, 4, 6, ...
    final AtomicInteger a2 = new AtomicInteger(2);
    Thread[] thread2 = Threads.makeArray();
    for (int ithread=0; ithread<thread2.length; ++ithread) {
      thread2[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a2.getAndAdd(2); i3<n3; i3=a2.getAndAdd(2))
            applyLhsSlice3(i3,epsilon,p2,p3,ep,x,y);
        }
      });
    }
    Threads.startAndJoin(thread2);
  }

  private static void applyLhsSlice3(
    int i3, float epsilon, 
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float epsilons = epsilon*epsilon;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float epi = ep[i3][i2][i1];
        float p2i = p2[i3][i2][i1];
        float p3i = p3[i3][i2][i1];
        float eps = epi*epi;
        float p2s = p2i*p2i;
        float p3s = p3i*p3i;
        float d11 = epsilons+(p2s+p3s)*eps;
        float d12 = p2i*eps;
        float d13 = p3i*eps;
        float d22 = eps;
        float d23 = 0.0f;
        float d33 = eps;
        float x000 = x00[i1 ];
        float x001 = x00[i1m];
        float x010 = x01[i1 ];
        float x011 = x01[i1m];
        float x100 = x10[i1 ];
        float x101 = x10[i1m];
        float x110 = x11[i1 ];
        float x111 = x11[i1m];
        float xa = x000-x111;
        float xb = x001-x110;
        float xc = x010-x101;
        float xd = x100-x011;
        float x1 = 0.0625f*(xa-xb+xc+xd);
        float x2 = 0.0625f*(xa+xb-xc+xd);
        float x3 = 0.0625f*(xa+xb+xc-xd);
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2+d23*x3;
        float y3 = d13*x1+d23*x2+d33*x3;
        float ya = y1+y2+y3;
        float yb = y1-y2+y3;
        float yc = y1+y2-y3;
        float yd = y1-y2-y3;
        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;
      }
    }
  }

  private static void szero(float[][] x) {
    zero(x);
  }
  private static void szero(float[][][] x) {
    if (PARALLEL) {
      szeroP(x);
    } else {
      szeroS(x);
    }
  }
  private static void szeroS(float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      szero(x[i3]);
  }
  private static void szeroP(final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            szero(x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }
}
