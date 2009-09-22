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
 * Requires estimates of slopes of locally linear or planar features in
 * in the 2D or 3D images, respectively.
 *
 * In 2D, the shifts (in samples) are functions s(x1,x2) that flatten
 * image features. Specifically, for a 2D image f(x1,x2), flattening 
 * is performed by computing g(x1,x2) = f(x1-s(x1,x2),x2). The shifts 
 * s(x1,x2) added to x1 in this expression are computed from local 
 * slopes dx1/dx2 and linearities provided for the 2D image f(x1,x2). 
 * After flattening, locally linear features with those slopes are
 * flat. In other words, corresponding features in the output image 
 * g(x1,x2) are orthogonal to the x1 axis.
 *
 * Flattening for 3D images works in the same way, but using shifts
 * s(x1,x2,x3) to compute g(x1,x2,x3) = f(x1-s(x1,x2,x3),x2,x3). The 
 * shifts are computed using local slopes dx1/dx2 and dx1/dx3 and 
 * planarities provided for the 3D image f(x1,x2,x3). After flattening, 
 * locally planar features with those slopes are flat; corresponding
 * features in the output image g(x1,x2,x3) are orthogonal to the x1 axis.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.19
 */
public class FlattenerCg {

  public FlattenerCg() {
  }

  public FlattenerCg(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public float[][] findShifts(float[][] p2) {
    return findShifts(p2,null);
  }

  public float[][] findShifts(float[][] p2, float[][] el) {
    int n1 = p2[0].length;
    int n2 = p2.length;
    float[][] r = new float[n2][n1]; // right-hand side
    float[][] s = new float[n2][n1]; // the shifts
    makeRhs(p2,el,r);
    //zeroReference(r);
    smoothTranspose(_sigma1,_sigma2,el,r);
    A2 a2 = new A2(_sigma1,_sigma2,_epsilon,p2,el);
    Vec2 vr = new Vec2(r);
    Vec2 vs = new Vec2(s);
    CgSolver cs = new CgSolver(_small,_niter);
    CgSolver.Info info = cs.solve(a2,vr,vs);
    System.out.println("vs.a0="+vs.a0);
    System.out.println("vs.a1="+vs.a1);
    smooth(_sigma1,_sigma2,el,s);
    checkShifts(s);
    invertShifts(s);
    return s;
  }
  private static void checkShifts(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    double s0 = 0.0;
    double s1 = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        s0 += s[i2][i1];
        s1 += i1*s[i2][i1];
      }
    }
    s0 = s0/n1/n2;
    s1 = s1*0.5/(n1-1)/n1/n2;
    System.out.println("s0="+s0+" s1="+s1);
  }

  public float[][][] findShifts(
    float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    int n1 = p2[0][0].length;
    int n2 = p2[0].length;
    int n3 = p2.length;
    float[][][] r = new float[n3][n2][n1]; // right-hand side
    float[][][] s = new float[n3][n2][n1]; // the shifts
    makeRhs(p2,p3,ep,r);
    //zeroReference(r);
    smoothTranspose(_sigma1,_sigma2,ep,r);
    A3 a3 = new A3(_sigma1,_sigma2,_epsilon,p2,p3,ep);
    Vec3 vr = new Vec3(r);
    Vec3 vs = new Vec3(s);
    CgSolver cs = new CgSolver(_small,_niter);
    CgSolver.Info info = cs.solve(a3,vr,vs);
    smooth(_sigma1,_sigma2,ep,s);
    checkShifts(s);
    invertShifts(s);
    return s;
  }
  private static void checkShifts(float[][][] s) {
    int n1 = s[0][0].length;
    int n2 = s[0].length;
    int n3 = s.length;
    double s0 = 0.0;
    double s1 = 0.0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          s0 += s[i3][i2][i1];
          s1 += i1*s[i3][i2][i1];
        }
      }
    }
    s0 = s0/n1/n2/n3;
    s1 = s1*0.5/(n1-1)/n1/n2/n3;
    System.out.println("s0="+s0+" s1="+s1);
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

  private static final boolean PARALLEL = true; // true if multi-threaded

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _epsilon = 0.000f; // damping for stability?
  private float _small = 0.01f; // stop CG iterations if residuals are small
  private int _niter = 1000; // maximum number of CG iterations

  // Vectors with extra elements for Lagrange multipliers.
  private static class Vec2 implements Vec {
    public VecArrayFloat2 v;
    public float a0,a1;
    public Vec2(float[][] a) {
      v = new VecArrayFloat2(a);
      a0 = 0.0f;
      a1 = 0.0f;
    }
    public float[][] getArray() {
      return v.getArray();
    }
    public double epsilon() {
      return v.epsilon();
    }
    public Vec2 clone() {
      Vec2 v2 = new Vec2();
      v2.v = v.clone();
      v2.a0 = a0;
      v2.a1 = a1;
      return v2;
    }
    public double dot(Vec vthat) {
      Vec2 v2 = (Vec2)vthat;
      return v.dot(v2.v)+a0*v2.a0+a1*v2.a1;
    }
    public double norm2() {
      double vnorm = v.norm2();
      return sqrt(a0*a0+a1*a1+vnorm*vnorm);
    }
    public void zero() {
      v.zero();
      a0 = 0.0f;
      a1 = 0.0f;
    }
    public void scale(double s) {
      v.scale(s);
      a0 *= s;
      a1 *= s;
    }
    public void add(double sthis, Vec vthat, double sthat) {
      Vec2 v2 = (Vec2)vthat;
      v.add(sthis,v2.v,sthat);
      float fthis = (float)sthis;
      float fthat = (float)sthat;
      a0 = a0*fthis+v2.a0*fthat;
      a1 = a1*fthis+v2.a1*fthat;
    }
    private Vec2() {
    }
  }
  private static class Vec3 implements Vec {
    public VecArrayFloat3 v;
    public float a0,a1;
    public Vec3(float[][][] a) {
      v = new VecArrayFloat3(a);
      a0 = 0.0f;
      a1 = 0.0f;
    }
    public float[][][] getArray() {
      return v.getArray();
    }
    public double epsilon() {
      return v.epsilon();
    }
    public Vec3 clone() {
      Vec3 v3 = new Vec3();
      v3.v = v.clone();
      v3.a0 = a0;
      v3.a1 = a1;
      return v3;
    }
    public double dot(Vec vthat) {
      Vec3 v3 = (Vec3)vthat;
      return v.dot(v3.v)+a0*v3.a0+a1*v3.a1;
    }
    public double norm2() {
      double vnorm = v.norm2();
      return sqrt(a0*a0+a1*a1+vnorm*vnorm);
    }
    public void zero() {
      v.zero();
      a0 = 0.0f;
      a1 = 0.0f;
    }
    public void scale(double s) {
      v.scale(s);
      a0 *= s;
      a1 *= s;
    }
    public void add(double sthis, Vec vthat, double sthat) {
      Vec3 v3 = (Vec3)vthat;
      v.add(sthis,v3.v,sthat);
      float fthis = (float)sthis;
      float fthat = (float)sthat;
      a0 = a0*fthis+v3.a0*fthat;
      a1 = a1*fthis+v3.a1*fthat;
    }
    private Vec3() {
    }
  }

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(float sigma1, float sigma2, float epsilon, 
       float[][] p2, float[][] el) 
    {
      int n1 = p2[0].length;
      int n2 = p2.length;
      _zr = new float[n1];
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _epsilon = epsilon;
      _p2 = p2;
      _el = el;
      double delta0 = 0.0/n1/n2; // turned off!
      double delta1 = 0.5*delta0/(n1-1);
      double first1 = -(n1/2)*delta1;
      _e0 = filldouble(delta0,n1);
      _e1 = rampdouble(first1,delta1,n1);
    }
    public void apply(Vec vx, Vec vy) {
      Vec2 v2x = (Vec2)vx;
      Vec2 v2y = (Vec2)vy;
      Vec2 v2z = v2x.clone();
      v2y.zero();
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = v2z.getArray();
      smooth(_sigma1,_sigma2,_el,z);
      //zeroReference(z,_zr);
      applyLhs(_p2,_el,z,y);
      double ya0 = 0.0;
      double ya1 = 0.0;
      double za0 = v2z.a0;
      double za1 = v2z.a1;
      int n1 = x[0].length;
      int n2 = x.length;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          double zi = z[i2][i1];
          double e0i = _e0[i1];
          double e1i = _e1[i1];
          y[i2][i1] += (float)(za0*e0i+za1*e1i);
          ya0 += zi*e0i;
          ya1 += zi*e1i;
        }
      }
      v2y.a0 = (float)ya0;
      v2y.a1 = (float)ya1;
      //setReference(_zr,y);
      smoothTranspose(_sigma1,_sigma2,_el,y);
      if (_epsilon>0.0f)
        v2y.v.add(1.0,v2x.v,_epsilon*_epsilon);
    }
    private float _sigma1;
    private float _sigma2;
    private float _epsilon;
    private float[][] _p2;
    private float[][] _el;
    private float[] _zr;
    private double[] _e0;
    private double[] _e1;
  }
  private static class A3 implements CgSolver.A {
    A3(float sigma1, float sigma23, float epsilon, 
       float[][][] p2, float[][][] p3, float[][][] ep) 
    {
      int n1 = p2[0][0].length;
      int n2 = p2[0].length;
      int n3 = p2.length;
      _zr = new float[n1];
      _sigma1 = sigma1;
      _sigma23 = sigma23;
      _epsilon = epsilon;
      _p2 = p2;
      _p3 = p3;
      _ep = ep;
      double delta0 = 0.0/n1/n2/n3; // turned off!
      double delta1 = 0.5*delta0/(n1-1);
      double first1 = -(n1/2)*delta1;
      _e0 = filldouble(delta0,n1);
      _e1 = rampdouble(first1,delta1,n1);
    }
    public void apply(Vec vx, Vec vy) {
      Vec3 v3x = (Vec3)vx;
      Vec3 v3y = (Vec3)vy;
      Vec3 v3z = v3x.clone();
      v3y.zero();
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = v3z.getArray();
      smooth(_sigma1,_sigma23,_ep,z);
      //zeroReference(z,_zr);
      applyLhs(_p2,_p3,_ep,z,y);
      double ya0 = 0.0;
      double ya1 = 0.0;
      double za0 = v3z.a0;
      double za1 = v3z.a1;
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            double zi = z[i3][i2][i1];
            double e0i = _e0[i1];
            double e1i = _e1[i1];
            y[i3][i2][i1] += (float)(za0*e0i+za1*e1i);
            ya0 += zi*e0i;
            ya1 += zi*e1i;
          }
        }
      }
      v3y.a0 = (float)ya0;
      v3y.a1 = (float)ya1;
      //setReference(_zr,y);
      smoothTranspose(_sigma1,_sigma23,_ep,y);
      if (_epsilon>0.0f)
        v3y.v.add(1.0,v3x.v,_epsilon*_epsilon);
    }
    private float _sigma1;
    private float _sigma23;
    private float _epsilon;
    private float[][][] _p2;
    private float[][][] _p3;
    private float[][][] _ep;
    private float[] _zr;
    private double[] _e0;
    private double[] _e1;
  }

  private static void setReference(float[] r, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    int i2 = n2/2;
    for (int i1=0; i1<n1; ++i1)
      x[i2][i1] = r[i1];
  }
  private static void setReference(float[] r, float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int i2 = n2/2;
    int i3 = n3/2;
    for (int i1=0; i1<n1; ++i1)
      x[i3][i2][i1] = r[i1];
  }
  private static void zeroReference(float[][] x) {
    zeroReference(x,null);
  }
  private static void zeroReference(float[][] x, float[] r) {
    int n1 = x[0].length;
    int n2 = x.length;
    int i2 = n2/2;
    for (int i1=0; i1<n1; ++i1) {
      if (r!=null) r[i1] = x[i2][i1];
      x[i2][i1] = 0.0f;
    }
  }
  private static void zeroReference(float[][][] x) {
    zeroReference(x,null);
  }
  private static void zeroReference(float[][][] x, float[] r) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int i2 = n2/2;
    int i3 = n3/2;
    for (int i1=0; i1<n1; ++i1) {
      if (r!=null) r[i1] = x[i3][i2][i1];
      x[i3][i2][i1] = 0.0f;
    }
  }

  // Smoothing operators (and transposes) for 2D and 3D images.
  private static void smooth(
    float sigma1, float sigma2, 
    float[][] s, float[][] x) 
  {
    smooth1(sigma1,x);
    smooth2(sigma2,s,x);
  }
  private static void smoothTranspose(
    float sigma1, float sigma2, 
    float[][] s, float[][] x) 
  {
    smooth2(sigma2,s,x);
    smooth1(sigma1,x);
  }
  private static void smooth(
    float sigma1, float sigma23,
    float[][][] s, float[][][] x) 
  {
    smooth1(sigma1,x);
    smooth2(sigma23,s,x);
    smooth3(sigma23,s,x);
  }
  private static void smoothTranspose(
    float sigma1, float sigma23, 
    float[][][] s, float[][][] x) 
  {
    smooth3(sigma23,s,x);
    smooth2(sigma23,s,x);
    smooth1(sigma1,x);
  }

  // Smoothers for dimension 1.
  private static void smooth1(float a, float[] x, float[] y) {
    int n1 = x.length;
    float s = (1.0f-a)/(1.0f+a);
    y[0] = s*x[0];
    for (int i1=1; i1<n1; ++i1)
      y[i1] = a*y[i1-1]+s*x[i1];
    float yip1 = 0.0f;
    for (int i1=n1-2; i1>=0; --i1) {
      float yi = a*(yip1+s*x[i1+1]);
      y[i1] += yi;
      yip1 = yi;
    }
  }
  private static void smooth1(float sigma, float[][] x) {
    float sigmaMin = sqrt(2.0f)/2.0f;
    if (sigma<sigmaMin)
      return;
    float a = (sigma-sigmaMin)/(sigma+sigmaMin);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] y = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      smooth1(a,x[i2],y);
      copy(y,x[i2]);
    }
  }
  private static void xsmooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = new float[n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      copy(x[i2],xt);
      lsf.apply(c,xt,x[i2]);
    }
  }
  private static void smooth1(float sigma, float[][][] x) {
    if (PARALLEL) {
      smooth1P(sigma,x);
    } else {
      smooth1S(sigma,x);
    }
  }
  private static void smooth1S(float sigma, float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      smooth1(sigma,x[i3]);
  }
  private static void smooth1P(final float sigma, final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger ai = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=ai.getAndIncrement(); i3<n3; i3=ai.getAndIncrement())
            smooth1(sigma,x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Smoothers for dimension 2.
  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }
  private static void smooth2(float sigma, float[][][] s, float[][][] x) {
    if (PARALLEL) {
      smooth2P(sigma,s,x);
    } else {
      smooth2S(sigma,s,x);
    }
  }
  private static void smooth2S(float sigma, float[][][] s, float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }
  }
  private static void smooth2P(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    final AtomicInteger ai = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=ai.getAndIncrement(); i3<n3; i3=ai.getAndIncrement()) {
            float[][] s3 = (s!=null)?s[i3]:null;
            float[][] x3 = x[i3];
            smooth2(sigma,s3,x3);
          }
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Smoothers for dimension 3.
  private static void smooth3(float sigma, float[][][] s, float[][][] x) {
    if (PARALLEL) {
      smooth3P(sigma,s,x);
    } else {
      smooth3S(sigma,s,x);
    }
  }
  private static void smooth3S(float sigma, float[][][] s, float[][][] x) {
    int n2 = x[0].length;
    int n3 = x.length;
    float[][] s2 = (s!=null)?new float[n3][]:null;
    float[][] x2 = new float[n3][];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }
  }
  private static void smooth3P(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    final AtomicInteger ai = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          float[][] s2 = (s!=null)?new float[n3][]:null;
          float[][] x2 = new float[n3][];
          for (int i2=ai.getAndIncrement(); i2<n2; i2=ai.getAndIncrement()) {
            for (int i3=0; i3<n3; ++i3) {
              if (s!=null)
                s2[i3] = s[i3][i2];
              x2[i3] = x[i3][i2];
            }
            smooth2(sigma,s2,x2);
          }
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  private static void makeRhs(float[][] p2, float[][] el, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = (el!=null)?el[i2][i1]:1.0f;
        float p2i = p2[i2][i1];
        float b12 = p2i*eli;
        float b22 = eli;
        float x2 = -0.5f*p2i*eli;
        float y1 = b12*x2;
        float y2 = b22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  private static void applyLhs(
    float[][] p2, float[][] el, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = (el!=null)?el[i2][i1]:1.0f;
        float p2i = p2[i2][i1];
        float els = eli*eli;
        float p2s = p2i*p2i;
        float d11 = p2s*els;
        float d12 = p2i*els;
        float d22 = els;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.25f*(xa-xb);
        float x2 = 0.25f*(xa+xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void makeRhs(
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float epi = (ep!=null)?ep[i3][i2][i1]:1.0f;
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float b12 = p2i*epi;
          float b13 = p3i*epi;
          float b22 = epi;
          float b33 = epi;
          float x2 = -0.25f*p2i*epi;
          float x3 = -0.25f*p3i*epi;
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
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    if (PARALLEL) {
      applyLhsP(p2,p3,ep,x,y);
    } else {
      applyLhsS(p2,p3,ep,x,y);
    }
  }

  private static void applyLhsS(
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3)
      applyLhsSlice3(i3,p2,p3,ep,x,y);
  }

  private static void applyLhsP(
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
            applyLhsSlice3(i3,p2,p3,ep,x,y);
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
            applyLhsSlice3(i3,p2,p3,ep,x,y);
        }
      });
    }
    Threads.startAndJoin(thread2);
  }

  private static void applyLhsSlice3(
    int i3,
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
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
        float epi = (ep!=null)?ep[i3][i2][i1]:1.0f;
        float p2i = p2[i3][i2][i1];
        float p3i = p3[i3][i2][i1];
        float eps = epi*epi;
        float p2s = p2i*p2i;
        float p3s = p3i*p3i;
        float d11 = (p2s+p3s)*eps;
        float d12 = p2i*eps;
        float d13 = p3i*eps;
        float d22 = eps;
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
        float y2 = d12*x1+d22*x2       ;
        float y3 = d13*x1       +d33*x3;
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

  // Post-processing and inversion of computed shifts.
  private static void cleanShifts(float[] s) {
    int n1 = s.length;
    for (int i1=1; i1<n1; ++i1) {
      if (s[i1]<=s[i1-1]-0.99f)
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
    float tmin = -n1;
    float tmax = n1+n1;
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
}
