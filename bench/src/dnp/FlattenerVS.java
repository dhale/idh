/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// FOR DEVELOPMENT ONLY
import java.awt.image.IndexColorModel;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * Flattening with vectors shifts of features in 2D and 3D images.
 * Requires estimates of slopes of locally linear or planar features
 * in the 2D or 3D images, respectively.
 * <p>
 * In 2D, the shifts (in samples) are functions s1(x1,x2) and s2(x1,x2)
 * that flatten image features. Specifically, for a 2D image f(u1,u2), 
 * flattening is performed by computing 
 * g(u1,u2) = f(x1-s1(x1,x2),x2-s2(x1,x2)). 
 * The shifts s1 and s2 subracted from x1 in this expression are computed 
 * from local slopes dx1/dx2 and linearities provided for the 2D image 
 * f(x1,x2). After flattening, locally linear features with those slopes 
 * are flat. In other words, corresponding features in the output image 
 * g(x1,x2) are orthogonal to the x1 axis.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.03.03
 */
public class FlattenerVS {

  public FlattenerVS() {
  }

  public FlattenerVS(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public float[][][] findShifts(double rotate, float[][] p2) {
    return findShifts(rotate,p2,null);
  }

  public float[][][] findShifts(double rotate, float[][] p2, float[][] el) {
    int n1 = p2[0].length;
    int n2 = p2.length;
    float[][][] p = new float[3][n2][n1];
    float[][] u1 = p[0], u2 = p[1], a = p[2];
    float ai = (float)rotate;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float dip = atan(p2[i2][i1]);
        u1[i2][i1] =  cos(dip);
        u2[i2][i1] = -sin(dip);
        a[i2][i1] = ai;
      }
    }
    float[][][] pinit = copy(p);
    float[][][] r = new float[2][n2][n1]; // right-hand side
    float[][][] s = new float[2][n2][n1]; // the shifts
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    VecArrayFloat3 vs = new VecArrayFloat3(s);
    Smoother2 smoother = new Smoother2(n1,n2,_sigma1,_sigma2,el);
    A2 a2 = new A2(smoother,p);
    CgSolver solver = new CgSolver(_small,_niter);
    for (int iter=0; iter<10; ++iter) {
      if (iter>0) {
        float[][][] t = copy(s);
        smoother.apply(t);
        copy(pinit,p);
        adjustParameters(t,p);
        vr.zero();
      }
      makeRhsS(p,r);
      smoother.applyTranspose(r);
      solver.solve(a2,vr,vs);
    }
    smoother.apply(s);
    return s;
  }
  private static class A2 implements CgSolver.A {
    A2(Smoother2 smoother, float[][][] p) {
      _smoother = smoother;
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      v3y.zero();
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = v3z.getArray();
      _smoother.apply(z);
      applyLhsS(_p,z,y);
      _smoother.applyTranspose(y);
    }
    private Smoother2 _smoother;
    private float[][][] _p;
  }

  public float[][][] findShiftsA(float[][] p2) {
    return findShiftsA(p2,null);
  }
  public float[][][] findShiftsA(float[][] p2, float[][] el) {
    int n1 = p2[0].length;
    int n2 = p2.length;
    float[][][] p = new float[2][n2][n1];
    float[][] u1 = p[0], u2 = p[1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float dip = atan(p2[i2][i1]);
        u1[i2][i1] =  cos(dip);
        u2[i2][i1] = -sin(dip);
      }
    }
    float[][][] pinit = copy(p);
    float[][][] r = new float[3][n2][n1]; // right-hand side
    float[][][] s = new float[3][n2][n1]; // shifts and a
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    VecArrayFloat3 vs = new VecArrayFloat3(s);
    Smoother2 smoother = new Smoother2(n1,n2,_sigma1,_sigma2,el);
    A2A a2a = new A2A(smoother,p);
    CgSolver solver = new CgSolver(_small,_niter);
    for (int iter=0; iter<10; ++iter) {
      if (iter>0) {
        float[][][] t = copy(s);
        smoother.apply(t);
        copy(pinit,p);
        adjustParameters(t,p);
        vr.zero();
      }
      makeRhsSA(p,r);
      smoother.applyTranspose(r);
      solver.solve(a2a,vr,vs);
    }
    smoother.apply(s);
    return s;
  }
  private static class A2A implements CgSolver.A {
    A2A(Smoother2 smoother, float[][][] p) {
      _smoother = smoother;
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      v3y.zero();
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = v3z.getArray();
      _smoother.apply(z);
      applyLhsSA(_p,z,y);
      _smoother.applyTranspose(y);
    }
    private Smoother2 _smoother;
    private float[][][] _p;
  }

  /**
   * Applies the specified shifts using sinc interpolation.
   * The returned array is a sampling of g(x) = f(x-s(x)).
   * @param s array {s1,s2} of shifts.
   * @param f input array to which shifts are to be applied.
   * @return array with shifts applied.
   */
  public static float[][] applyShifts(float[][][] s, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] s1 = s[0], s2 = s[1];
    SincInterpolator si = new SincInterpolator();
    si.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = si.interpolate(i1-s1[i2][i1],i2-s2[i2][i1]);
      }
    }
    return g;
  }

  /**
   * Applies the specified inverse shifts using linear interpolation.
   * The returned array is a sampling of g(x) = f(x+s(x)).
   * @param f array to which shifts are applied.
   * @param s array {s1,s2} of inverse shifts.
   * @return array of inverse shifts.
   */
  public static float[][] applyInverseShiftsL(float[][][] s, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] s1 = s[0], s2 = s[1];
    float[][] g = new float[n2][n1];
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = li.interpolate(i1+s1[i2][i1],i2+s2[i2][i1]);
      }
    }
    return g;
  }

  /**
   * Gets inverse of the specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return array {s1,s2} of inverse shifts.
   */
  public static float[][][] getInverseShifts(float[][][] s) {
    int n1 = s[0][0].length;
    int n2 = s[0].length;
    float[][] r1 = s[0], r2 = s[1];
    float[][] s1 = copy(r1);
    float[][] s2 = copy(r2);
    LinearInterpolator li1 = new LinearInterpolator();
    li1.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li1.setUniform(n1,1.0,0.0,n2,1.0,0.0,r1);
    LinearInterpolator li2 = new LinearInterpolator();
    li2.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li2.setUniform(n1,1.0,0.0,n2,1.0,0.0,r2);
    for (int i2=0; i2<n2; ++i2) {
      double x2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        double x1 = i1;
        double s1i = s1[i2][i1];
        double s2i = s2[i2][i1];
        double s1p,s2p,ds1,ds2;
        do {
          double u1 = x1+s1i;
          double u2 = x2+s2i;
          s1p = s1i;
          s2p = s2i;
          s1i = li1.interpolate(u1,u2);
          s2i = li2.interpolate(u1,u2);
          ds1 = s1i-s1p;
          ds2 = s2i-s2p;
        } while (ds1*ds1+ds2*ds2>0.0001);
        s1[i2][i1] = (float)s1i;
        s2[i2][i1] = (float)s2i;
      }
    }
    return new float[][][]{s1,s2};
  }

  /**
   * Gets partial derivatives of the specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return array {s11,s12,s21,s22} of partial derivatives.
   */
  public static float[][][] getDerivatives(float[][][] s) {
    int n1 = s[0][0].length;
    int n2 = s[0].length;
    float[][] s1 = s[0], s2 = s[1];
    float[][][] t = new float[4][n2][n1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float s100 = s1[i2  ][i1  ];
        float s101 = s1[i2  ][i1-1];
        float s110 = s1[i2-1][i1  ];
        float s111 = s1[i2-1][i1-1];
        float s200 = s2[i2  ][i1  ];
        float s201 = s2[i2  ][i1-1];
        float s210 = s2[i2-1][i1  ];
        float s211 = s2[i2-1][i1-1];
        float s1a = s100-s111;
        float s1b = s101-s110;
        float s2a = s200-s211;
        float s2b = s201-s210;
        float s11 = 0.5f*(s1a-s1b);
        float s12 = 0.5f*(s1a+s1b);
        float s21 = 0.5f*(s2a-s2b);
        float s22 = 0.5f*(s2a+s2b);
        t[0][i2][i1] = s11;
        t[1][i2][i1] = s12;
        t[2][i2][i1] = s21;
        t[3][i2][i1] = s22;
      }
    }
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    for (int i=0; i<4; ++i) {
      float[][] ti = copy(n1-1,n2-1,1,1,t[i]);
      li.setUniform(n1-1,1.0,0.5,n2-1,1.0,0.5,ti);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[i][i2][i1] = li.interpolate(i1,i2);
        }
      }
    }
    return t;
  }

  /**
   * Gets determinants of the Jacobian from specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return the determinants.
   */
  public static float[][] getDeterminantsFromShifts(float[][][] s) {
    float[][][] ds = getDerivatives(s);
    return getDeterminantsFromDerivatives(ds);
  }

  /**
   * Gets determinants of the Jacobian from specified derivatives of shifts.
   * @param ds array {s11,s12,s21,s22} of derivatives of shifts.
   * @return the determinants.
   */
  public static float[][] getDeterminantsFromDerivatives(float[][][] ds) {
    int n1 = ds[0][0].length;
    int n2 = ds[0].length;
    float[][] s11 = ds[0], s12 = ds[1], s21 = ds[2], s22 = ds[3];
    float[][] d = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float s11i = s11[i2][i1];
        float s12i = s12[i2][i1];
        float s21i = s21[i2][i1];
        float s22i = s22[i2][i1];
        d[i2][i1] = (1.0f-s11i)*(1.0f-s22i)-s12i*s21i;
      }
    }
    return d;
  }

  /**
   * Gets unit normal vectors from specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return array {u1,u2} of unit normal vectors.
   */
  public static float[][][] getNormalsFromShifts(float[][][] s) {
    float[][][] ds = getDerivatives(s);
    return getNormalsFromDerivatives(ds);
  }

  /**
   * Gets unit normal vectors from specified derivatives of shifts.
   * @param ds array {s11,s12,s21,s22} of derivatives of shifts.
   * @return array {u1,u2} of unit normal vectors.
   */
  public static float[][][] getNormalsFromDerivatives(float[][][] ds) {
    int n1 = ds[0][0].length;
    int n2 = ds[0].length;
    float[][] s12 = ds[1], s22 = ds[3];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = (1.0f-s22[i2][i1]);
        float u2i = s12[i2][i1];
        float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
        u1i *= usi;
        u2i *= usi;
        u1[i2][i1] = u1i;
        u2[i2][i1] = u2i;
      }
    }
    return new float[][][]{u1,u2};
  }

  /**
   * Gets the coefficient a from specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return the coefficients a
   */
  public static float[][] getAFromShifts(float[][][] s) {
    float[][][] ds = getDerivatives(s);
    return getAFromDerivatives(ds);
  }

  /**
   * Gets the coefficient a from derivatives of shifts.
   * @param ds array {s11,s12,s21,s22} of derivatives of shifts.
   * @return the coefficients a
   */
  public static float[][] getAFromDerivatives(float[][][] ds) {
    final float w1 = 1.001f;
    final float w2 = 1.001f;
    final float w3 = 1.001f;
    final float w4 = 1.001f;
    //final float w1 = W1;
    //final float w2 = W2;
    //final float w3 = W3;
    //final float w4 = W4;
    int n1 = ds[0][0].length;
    int n2 = ds[0].length;
    float[][][] u = getNormalsFromDerivatives(ds);
    float[][] u1 = u[0], u2 = u[1];
    float[][] s11 = ds[0], s12 = ds[1], s21 = ds[2], s22 = ds[3];
    float[][] anum = new float[n2][n1];
    float[][] aden = new float[n2][n1];
    float[][] a = anum;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float c1 =  (1.0f-u1i);
        float c2 = -(1.0f-u1i)*u2i/u1i;
        float c3 = -u2i;
        float c4 =  (1.0f-u1i);
        float b1 = s11[i2][i1];
        float b2 = s12[i2][i1]-u2i/u1i;
        float b3 = s21[i2][i1];
        float b4 = s22[i2][i1];
        anum[i2][i1] = w1*b1*c1+w2*b2*c2+w3*b3*c3+w4*b4*c4;
        aden[i2][i1] = w1*c1*c1+w2*c2*c2+w3*c3*c3+w4*c4*c4;
      }
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(40.0);
    rgf.apply00(anum,anum);
    rgf.apply00(aden,aden);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[i2][i1] = anum[i2][i1]/aden[i2][i1];
      }
    }
    return a;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _small = 0.005f; // stop CG iterations if residuals are small
  private int _niter = 20; // maximum number of CG iterations

  // Conjugate-gradient operators.

  // Smoothers used as preconditioners.
  private static class Smoother2 {
    public Smoother2(
      int n1, int n2, float sigma1, float sigma2, float[][] el) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _el = el;
    }
    public void apply(float[][][] x) {
      int nx = x.length;
      for (int i=0; i<nx; ++i) {
        float sigma1 = (i<2)?_sigma1:4.0f*_sigma1;
        float sigma2 = (i<2)?_sigma2:4.0f*_sigma2;
        smooth1(sigma1,x[i]);
        smooth2(sigma2,_el,x[i]);
        if (i<2)
          subtractMean(x[i]);
      }
    }
    public void applyTranspose(float[][][] x) {
      int nx = x.length;
      for (int i=0; i<nx; ++i) {
        if (i<2)
          subtractMean(x[i]);
        float sigma1 = (i<2)?_sigma1:4.0f*_sigma1;
        float sigma2 = (i<2)?_sigma2:4.0f*_sigma2;
        smooth2(sigma2,_el,x[i]);
        smooth1(sigma1,x[i]);
      }
    }
    private float _sigma1,_sigma2;
    private float[][] _el;
  }

  // Subtracts the mean from the specified array.
  private static void subtractMean(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float xavg = 0.0f;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        xavg += x[i2][i1];
    xavg /= (float)n1*(float)n2;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] -= xavg;
  }

  // Smoothing for dimension 1.
  private static void smooth1(float a, float[] x, float[] y) {
    int n1 = x.length;
    float s = (1.0f-a)/(1.0f+a);
    y[0] = x[0]*s/(1.0f-a); // zero-slope b.c.
    for (int i1=1; i1<n1; ++i1)
      y[i1] = a*y[i1-1]+s*x[i1];
    float yip1 = x[n1-1]*s*a/(1.0f-a); // zero-slope b.c.
    y[n1-1] += yip1;
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

  // Smoothing for dimension 2.
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

  private void adjustParameters(float[][][] s, float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    int np = p.length;
    float[][] s1 = s[0], s2 = s[1];
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    for (int ip=0; ip<np; ++ip) {
      float[][] q = new float[n2][n1];
      li.setUniform(n1,1.0,0.0,n2,1.0,0.0,p[ip]);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          q[i2][i1] = li.interpolate(i1-s1[i2][i1],i2-s2[i2][i1]);
        }
      }
      p[ip] = q;
    }
    float[][] u1 = p[0], u2 = p[1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float uss = 1.0f/sqrt(u1i*u1i+u2i*u2i);
        u1[i2][i1] *= uss;
        u2[i2][i1] *= uss;
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      //u1[i2][1] = u1[i2][2]; u1[i2][n1-2] = u1[i2][n1-3];
      //u2[i2][1] = u2[i2][2]; u2[i2][n1-2] = u2[i2][n1-3];
      u1[i2][0] = u1[i2][1]; u1[i2][n1-1] = u1[i2][n1-2];
      u2[i2][0] = u2[i2][1]; u2[i2][n1-1] = u2[i2][n1-2];
    }
    for (int i1=0; i1<n1; ++i1) {
      //u1[1][i1] = u1[2][i1]; u1[n2-2][i1] = u1[n2-3][i1];
      //u2[1][i1] = u2[2][i1]; u2[n2-2][i1] = u2[n2-3][i1];
      u1[0][i1] = u1[1][i1]; u1[n2-1][i1] = u1[n2-2][i1];
      u2[0][i1] = u2[1][i1]; u2[n2-1][i1] = u2[n2-2][i1];
    }
  }

  private static final float W1 = 0.001f; // r11 or u1*r11+u2*r21
  private static final float W2 = 1.001f; // r12 or u1*r12+u2*r22
  private static final float W3 = 1.001f; // r21 or u1*r21-u2*r11
  private static final float W4 = 1.001f; // r22 or u1*r22-u2*r12

  // Simple four equations.
  private static void makeRhsS(
    float[][][] p, float[][][] y) 
  {
    final float w1 = W1*0.5f; // weights include
    final float w2 = W2*0.5f; // scaling by 0.5
    final float w3 = W3*0.5f; // in the gradient
    final float w4 = W4*0.5f; // transpose
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] u1 = p[0], u2 = p[1], a = p[2];
    float[][] y1 = y[0], y2 = y[1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float ai = a[i2][i1];
        float x11 = ai*(1.0f-u1i);
        float x12 = (1.0f-ai*(1.0f-u1i))*u2i/u1i;
        float x21 = -ai*u2i;;
        float x22 = ai*(1.0f-u1i);
        float y11 = w1*x11;
        float y12 = w2*x12;
        float y21 = w3*x21;
        float y22 = w4*x22;
        float y1a = y11+y12;
        float y1b = y11-y12;
        float y2a = y21+y22;
        float y2b = y21-y22;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }
  private static void applyLhsS(
    float[][][] p, float[][][] x, float[][][] y) 
  {
    final float w1 = W1*0.5f*0.5f; // weights include
    final float w2 = W2*0.5f*0.5f; // scaling by 0.5 in
    final float w3 = W3*0.5f*0.5f; // both gradient and 
    final float w4 = W4*0.5f*0.5f; // gradient transpose
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float[][] u1 = p[0], u2 = p[1];
    float[][] x1 = x[0], x2 = x[1];
    float[][] y1 = y[0], y2 = y[1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float x100 = x1[i2  ][i1  ];
        float x101 = x1[i2  ][i1-1];
        float x110 = x1[i2-1][i1  ];
        float x111 = x1[i2-1][i1-1];
        float x200 = x2[i2  ][i1  ];
        float x201 = x2[i2  ][i1-1];
        float x210 = x2[i2-1][i1  ];
        float x211 = x2[i2-1][i1-1];
        float x1a = x100-x111;
        float x1b = x101-x110;
        float x2a = x200-x211;
        float x2b = x201-x210;
        float x11 = x1a-x1b;
        float x12 = x1a+x1b;
        float x21 = x2a-x2b;
        float x22 = x2a+x2b;
        float y11 = w1*x11;
        float y12 = w2*x12;
        float y21 = w3*x21;
        float y22 = w4*x22;
        float y1a = y11+y12;
        float y1b = y11-y12;
        float y2a = y21+y22;
        float y2b = y21-y22;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }


  // Transformed four equations.
  private static void makeRhsT(
    float[][][] p, float[][][] y) 
  {
    final float w1 = W1*0.5f; // weights include
    final float w2 = W2*0.5f; // scaling by 0.5
    final float w3 = W3*0.5f; // in the gradient
    final float w4 = W4*0.5f; // transpose
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] u1 = p[0], u2 = p[1], a = p[2];
    float[][] y1 = y[0], y2 = y[1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float ai = a[i2][i1];
        float t1 = ai*(u1i-1.0f);
        float t2 = u2i;
        float t3 = -ai*u2i;
        float t4 = (ai*(1.0f-u1i)-u2i*u2i)/u1i;
        t1 *= w1;
        t2 *= w2;
        t3 *= w3;
        t4 *= w4;
        float y11 = u1i*t1-u2i*t3;
        float y12 = u1i*t2-u2i*t4;
        float y21 = u1i*t3+u2i*t1;
        float y22 = u1i*t4+u2i*t2;
        float y1a = y11+y12;
        float y1b = y11-y12;
        float y2a = y21+y22;
        float y2b = y21-y22;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }
  private static void applyLhsT(
    float[][][] p, float[][][] x, float[][][] y) 
  {
    final float w1 = W1*0.5f*0.5f; // weights include
    final float w2 = W2*0.5f*0.5f; // scaling by 0.5 in
    final float w3 = W3*0.5f*0.5f; // both gradient and 
    final float w4 = W4*0.5f*0.5f; // gradient transpose
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float[][] u1 = p[0], u2 = p[1];
    float[][] x1 = x[0], x2 = x[1];
    float[][] y1 = y[0], y2 = y[1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float x100 = x1[i2  ][i1  ];
        float x101 = x1[i2  ][i1-1];
        float x110 = x1[i2-1][i1  ];
        float x111 = x1[i2-1][i1-1];
        float x200 = x2[i2  ][i1  ];
        float x201 = x2[i2  ][i1-1];
        float x210 = x2[i2-1][i1  ];
        float x211 = x2[i2-1][i1-1];
        float x1a = x100-x111;
        float x1b = x101-x110;
        float x2a = x200-x211;
        float x2b = x201-x210;
        float x11 = x1a-x1b;
        float x12 = x1a+x1b;
        float x21 = x2a-x2b;
        float x22 = x2a+x2b;
        float t1 = u1i*x11+u2i*x21;
        float t2 = u1i*x12+u2i*x22;
        float t3 = u1i*x21-u2i*x11;
        float t4 = u1i*x22-u2i*x12;
        t1 *= w1;
        t2 *= w2;
        t3 *= w3;
        t4 *= w4;
        float y11 = u1i*t1-u2i*t3;
        float y12 = u1i*t2-u2i*t4;
        float y21 = u1i*t3+u2i*t1;
        float y22 = u1i*t4+u2i*t2;
        float y1a = y11+y12;
        float y1b = y11-y12;
        float y2a = y21+y22;
        float y2b = y21-y22;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }

  // Simple four equations, solving for a.
  private static void makeRhsSA(
    float[][][] p, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] u1 = p[0], u2 = p[1];
    float[][] y1 = y[0], y2 = y[1], ya = y[2];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float s1 = -(1.0f-u1i);
        float s2 =  (1.0f-u1i)*u2i/u1i;
        float s3 =  u2i;
        float s4 = -(1.0f-u1i);
        float t1 = 0.0f;
        float t2 = u2i/u1i;
        float t3 = 0.0f;
        float t4 = 0.0f;
        t1 *= W1;
        t2 *= W2;
        t3 *= W3;
        t4 *= W4;
        float y11 = t1;
        float y12 = t2;
        float y21 = t3;
        float y22 = t4;
        float yai = s1*t1+s2*t2+s3*t3+s4*t4;
        float y1a = 0.5f*(y11+y12);
        float y1b = 0.5f*(y11-y12);
        float y2a = 0.5f*(y21+y22);
        float y2b = 0.5f*(y21-y22);
        ya[i2  ][i1  ]  = yai;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }
  private static void applyLhsSA(
    float[][][] p, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float[][] u1 = p[0], u2 = p[1];
    float[][] x1 = x[0], x2 = x[1], xa = x[2];
    float[][] y1 = y[0], y2 = y[1], ya = y[2];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float xai  = xa[i2  ][i1  ];
        float x100 = x1[i2  ][i1  ];
        float x101 = x1[i2  ][i1-1];
        float x110 = x1[i2-1][i1  ];
        float x111 = x1[i2-1][i1-1];
        float x200 = x2[i2  ][i1  ];
        float x201 = x2[i2  ][i1-1];
        float x210 = x2[i2-1][i1  ];
        float x211 = x2[i2-1][i1-1];
        float x1a = x100-x111;
        float x1b = x101-x110;
        float x2a = x200-x211;
        float x2b = x201-x210;
        float x11 = 0.5f*(x1a-x1b);
        float x12 = 0.5f*(x1a+x1b);
        float x21 = 0.5f*(x2a-x2b);
        float x22 = 0.5f*(x2a+x2b);
        float s1 = -(1.0f-u1i);
        float s2 =  (1.0f-u1i)*u2i/u1i;
        float s3 =  u2i;
        float s4 = -(1.0f-u1i);
        float t1 = x11+xai*s1;
        float t2 = x12+xai*s2;
        float t3 = x21+xai*s3;
        float t4 = x22+xai*s4;
        t1 *= W1;
        t2 *= W2;
        t3 *= W3;
        t4 *= W4;
        float y11 = t1;
        float y12 = t2;
        float y21 = t3;
        float y22 = t4;
        float yai = s1*t1+s2*t2+s3*t3+s4*t4;
        float y1a = 0.5f*(y11+y12);
        float y1b = 0.5f*(y11-y12);
        float y2a = 0.5f*(y21+y22);
        float y2b = 0.5f*(y21-y22);
        ya[i2  ][i1  ]  = yai;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }

  // FOR DEVELOPMENT ONLY
  private static final IndexColorModel gray = ColorMap.GRAY; 
  private static final IndexColorModel jet = ColorMap.JET; 
  private static void plot(float[][] x) {
    plot(x,null,null);
  }
  private static void plot(float[][] x, String title) {
    plot(x,title,null);
  }
  private static void plot(float[][] x, String title, IndexColorModel icm) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    if (title!=null)
      sp.setTitle(title);
    PixelsView pv = sp.addPixels(x);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    if (icm!=null)
      pv.setColorModel(icm);
    sp.addColorBar();
    sp.setSize(1400,800);
  }

  public static void main(String[] args) {
    int n1 = 251;
    int n2 = 501;
    //float[][] x = rampfloat(-n2/2.0f,0.0f,1.0f,n1,n2);
    float[][] x = fillfloat(1.0f,n1,n2);
    float[][] y = copy(x);
    smooth1(100.0f,y);
    //smooth2(100.0f,null,y);
    plot(x,"x",jet);
    plot(y,"y",jet);
  }
}
