/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Estimates displacement vectors for two images. For example, given two 
 * 2-D images f(x1,x2) and g(x1,x2), a shift finder estimates two vector 
 * components of displacement u1(x1,x2) and u2(x1,x2) such that 
 * f(x1,x2) ~ g(x1+u1(x1,x2),x2+u2(x1,x2)).
 * <p>
 * Like the images f and g, the components of displacement are sampled
 * functions of coordinates x1 and x2. That is, displacements may vary 
 * from sample to sample. The components u1 and u2 represent displacements 
 * in the x1 and x2 coordinate directions, respectively.
 * <p>
 * This shift finder estimates each component of displacement using local
 * cross-correlations. For each image sample, the estimated shift is that 
 * which yields the maximum correlation coefficient. This coefficient is
 * found by quadratic interpolation of correlation functions sampled at
 * integer lags.
 * <p>
 * Methods are provided to find and compensate for each component of shift 
 * sequentially. As each component is found, that component can be removed 
 * from the image g before estimating another component. For example, again 
 * for 2-D images f(x1,x2) and g(x1,x2), we might first estimate u1(x1,x2). 
 * If we then compute an image h(x1,x2) = g(x1+u1(x1,x2),x2), we can use
 * f(x1,x2) and h(x1,x2) to estimate u2(x1,x2). By repeating this process
 * sequentially, we obtain estimates for both u1(x1,x2) and u2(x1,x2) such
 * that f(x1,x2) ~ g(x1+u1(x1,x2),x2+u2(x1,x2)).
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.11.18
 */
public class ShiftFinder {

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has the
   * same correlation window half-width for all dimensions.
   * @param sigma the correlation window half-width; must not be less than 1.
   */
  public ShiftFinder(double sigma) {
    this(sigma,sigma,sigma);
  }

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has half-width 
   * sigma1 for the 1st dimension and half-width sigma2 for 2nd and higher 
   * dimensions.
   * @param sigma1 correlaton window half-width for 0st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd and higher 
   *  dimensions; must not be less than 1.
   */
  public ShiftFinder(double sigma1, double sigma2) {
    this(sigma1,sigma2,sigma2);
  }

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has half-width 
   * sigma1 for the 1st dimension, half-width sigma2 for the 2nd dimension, 
   * and half-width sigma3 for 3rd and higher dimensions.
   * @param sigma1 correlation window half-width for 1st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd dimension;
   *  must not be less than 1.
   * @param sigma3 correlation window half-width for 3rd and higher 
   *  dimensions; must not be less than 1.
   */
  public ShiftFinder(double sigma1, double sigma2, double sigma3) {
    Check.argument(sigma1>=1.0,"sigma1>=1.0");
    Check.argument(sigma2>=1.0,"sigma2>=1.0");
    Check.argument(sigma3>=1.0,"sigma3>=1.0");
    _lcfSimple = new LocalCorrelationFilter(
      LocalCorrelationFilter.Type.SIMPLE,
      LocalCorrelationFilter.Window.GAUSSIAN,
      sigma1,sigma2,sigma3);
    _lcfSymmetric = new LocalCorrelationFilter(
      LocalCorrelationFilter.Type.SYMMETRIC,
      LocalCorrelationFilter.Window.GAUSSIAN,
      sigma1,sigma2,sigma3);
    _rgfSmooth = new RecursiveGaussianFilter(1.0);
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    _sigma1 = sigma1;
    _sigma2 = sigma2;
    _sigma3 = sigma3;
  }

  /**
   * Finds shifts in the 1st dimension.
   * @param min1 the minimum shift.
   * @param max1 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find1(
    int min1, int max1, float[][] f, float[][] g, float[][] u) 
  {
    findShifts(1,min1,max1,f,g,u);
  }

  /**
   * Finds shifts in the 2nd dimension.
   * @param min2 the minimum shift.
   * @param max2 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find2(
    int min2, int max2, float[][] f, float[][] g, float[][] u) 
  {
    findShifts(2,min2,max2,f,g,u);
  }

  /**
   * Finds shifts in the 1st dimension.
   * @param min1 the minimum shift.
   * @param max1 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find1(
    int min1, int max1, float[][][] f, float[][][] g, float[][][] u) 
  {
    findShifts(1,min1,max1,f,g,u);
  }

  /**
   * Finds shifts in the 2nd dimension.
   * @param min2 the minimum shift.
   * @param max2 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find2(
    int min2, int max2, float[][][] f, float[][][] g, float[][][] u) 
  {
    findShifts(2,min2,max2,f,g,u);
  }

  /**
   * Finds shifts in the 3rd dimension.
   * @param min3 the minimum shift.
   * @param max3 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find3(
    int min3, int max3, float[][][] f, float[][][] g, float[][][] u) 
  {
    findShifts(3,min3,max3,f,g,u);
  }

  /**
   * Applies specified shift in the 1st dimension.
   * @param du array of changes to shifts in 1st dimension.
   * @param g input array.
   * @param h output array.
   * @param u1 updated array of shifts in 1st dimension.
   * @param u2 updated array of shifts in 2nd dimension.
   */
  public void shift1(
    float[][] du, float[][] g, float[][] h, float[][] u1, float[][] u2) 
  {
    int n1 = g[0].length;
    int n2 = g.length;
    float[] u = new float[n1];
    float[] x = new float[n1];
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      x[i1] = (float)(i1);
    _si.setUniformSampling(n1,1.0,0.0);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        u[i1] = u2[i2][i1];
        float d = du[i2][i1];
        y[i1] = x[i1]+d;
        u1[i2][i1] += d;
      }
      _si.setUniformSamples(g[i2]);
      _si.interpolate(n1,y,h[i2]);
      _si.setUniformSamples(u);
      _si.interpolate(n1,y,u2[i2]);
    }
  }

  /**
   * Applies specified shift in the 2nd dimension.
   * @param du array of changes to shifts in 2nd dimension.
   * @param g input array.
   * @param h output array.
   * @param u1 updated array of shifts in 1st dimension.
   * @param u2 updated array of shifts in 2nd dimension.
   */
  public void shift2(
    float[][] du, float[][] g, float[][] h, float[][] u1, float[][] u2) 
  {
    int n1 = g[0].length;
    int n2 = g.length;
    float[] u = new float[n2];
    float[] x = new float[n2];
    float[] y = new float[n2];
    float[] gt = new float[n2];
    float[] ht = new float[n2];
    float[] ut = new float[n2];
    for (int i2=0; i2<n2; ++i2)
      x[i2] = (float)(i2);
    _si.setUniformSampling(n2,1.0,0.0);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        gt[i2] = g[i2][i1];
        u[i2] = u1[i2][i1];
        float d = du[i2][i1];
        y[i2] = x[i2]+d;
        u2[i2][i1] += d;
      }
      _si.setUniformSamples(gt);
      _si.interpolate(n2,y,ht);
      _si.setUniformSamples(u);
      _si.interpolate(n2,y,ut);
      for (int i2=0; i2<n2; ++i2) {
        h[i2][i1] = ht[i2];
        u1[i2][i1] = ut[i2];
      }
    }
  }

  /**
   * Applies specified shift in the 1st dimension.
   * @param du array of changes to shifts in 1st dimension.
   * @param g input array.
   * @param h output array.
   * @param u1 updated array of shifts in 1st dimension.
   * @param u2 updated array of shifts in 2nd dimension.
   * @param u3 updated array of shifts in 3rd dimension.
   */
  public void shift1(
    float[][][] du, float[][][] g, float[][][] h, 
    float[][][] u1, float [][][] u2, float[][][] u3) 
  {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    float[] u = new float[n1];
    float[] v = new float[n1];
    float[] x = new float[n1];
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      x[i1] = (float)(i1);
    _si.setUniformSampling(n1,1.0,0.0);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          u[i1] = u2[i3][i2][i1];
          v[i1] = u3[i3][i2][i1];
          float d = du[i3][i2][i1];
          y[i1] = x[i1]+d;
          u1[i3][i2][i1] += d;
        }
        _si.setUniformSamples(g[i3][i2]);
        _si.interpolate(n1,y,h[i3][i2]);
        _si.setUniformSamples(u);
        _si.interpolate(n1,y,u2[i3][i2]);
        _si.setUniformSamples(v);
        _si.interpolate(n1,y,u3[i3][i2]);
      }
    }
  }

  /**
   * Applies specified shift in the 2nd dimension.
   * @param du array of changes to shifts in 2nd dimension.
   * @param g input array.
   * @param h output array.
   * @param u1 updated array of shifts in 1st dimension.
   * @param u2 updated array of shifts in 2nd dimension.
   * @param u3 updated array of shifts in 3rd dimension.
   */
  public void shift2(
    float[][][] du, float[][][] g, float[][][] h, 
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    float[] u = new float[n2];
    float[] v = new float[n2];
    float[] x = new float[n2];
    float[] y = new float[n2];
    float[] gt = new float[n2];
    float[] ht = new float[n2];
    float[] ut = new float[n2];
    float[] vt = new float[n2];
    for (int i2=0; i2<n2; ++i2)
      x[i2] = (float)(i2);
    _si.setUniformSampling(n2,1.0,0.0);
    for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        for (int i2=0; i2<n2; ++i2) {
          gt[i2] = g[i3][i2][i1];
          u[i2] = u1[i3][i2][i1];
          v[i2] = u3[i3][i2][i1];
          float d = du[i3][i2][i1];
          y[i2] = x[i2]+d;
          u2[i3][i2][i1] += d;
        }
        _si.setUniformSamples(gt);
        _si.interpolate(n2,y,ht);
        _si.setUniformSamples(u);
        _si.interpolate(n2,y,ut);
        _si.setUniformSamples(v);
        _si.interpolate(n2,y,vt);
        for (int i2=0; i2<n2; ++i2) {
          h[i3][i2][i1] = ht[i2];
          u1[i3][i2][i1] = ut[i2];
          u3[i3][i2][i1] = vt[i2];
        }
      }
    }
  }

  /**
   * Applies specified shift in the 3rd dimension.
   * @param du array of changes to shifts in 3rd dimension.
   * @param g input array.
   * @param h output array.
   * @param u1 updated array of shifts in 1st dimension.
   * @param u2 updated array of shifts in 2nd dimension.
   * @param u3 updated array of shifts in 3rd dimension.
   */
  public void shift3(
    float[][][] du, float[][][] g, float[][][] h, 
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    float[] u = new float[n3];
    float[] v = new float[n3];
    float[] x = new float[n3];
    float[] y = new float[n3];
    float[] gt = new float[n3];
    float[] ht = new float[n3];
    float[] ut = new float[n3];
    float[] vt = new float[n3];
    for (int i3=0; i3<n3; ++i3)
      x[i3] = (float)(i3);
    _si.setUniformSampling(n3,1.0,0.0);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int i3=0; i3<n3; ++i3) {
          gt[i3] = g[i3][i2][i1];
          u[i3] = u1[i3][i2][i1];
          v[i3] = u2[i3][i2][i1];
          float d = du[i3][i2][i1];
          y[i3] = x[i3]+d;
          u3[i3][i2][i1] += d;
        }
        _si.setUniformSamples(gt);
        _si.interpolate(n3,y,ht);
        _si.setUniformSamples(u);
        _si.interpolate(n3,y,ut);
        _si.setUniformSamples(v);
        _si.interpolate(n3,y,vt);
        for (int i3=0; i3<n3; ++i3) {
          h[i3][i2][i1] = ht[i3];
          u1[i3][i2][i1] = ut[i3];
          u2[i3][i2][i1] = vt[i3];
        }
      }
    }
  }

  /**
   * Applies local prediction-error (spectal whitening) filters.
   * The input and output arrays f and g can be the same array.
   * @param f the input array.
   * @param g the output array.
   */
  public void whiten(float[][] f, float[][] g) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] r00 = new float[n2][n1];
    float[][] rpm = new float[n2][n1];
    float[][] rp0 = new float[n2][n1];
    float[][] r0p = new float[n2][n1];
    _lcfSymmetric.setInputs(f,f);
    _lcfSymmetric.correlate( 0, 0,r00);
    _lcfSymmetric.correlate( 1,-1,rpm);
    _lcfSymmetric.correlate( 1, 0,rp0);
    _lcfSymmetric.correlate( 0, 1,r0p);
    float[][] s = rp0;
    float[][] t = r0p;
    for (int i2=0; i2<n2; ++i2)
      s[i2][0] = 0.0f;
    for (int i1=0; i1<n1; ++i1)
      s[0][i1] = 0.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        double b1 = rp0[i2][i1];
        double b2 = r0p[i2][i1];
        double a11 = r00[i2][i1];
        double a21 = rpm[i2][i1];
        double a22 = a11;
        double l11 = sqrt(a11);
        double l21 = a21/l11;
        double d22 = a22-l21*l21;
        double x1 = 0.0;
        double x2 = 0.0;
        if (d22>0.0) {
          double l22 = sqrt(d22);
          double v1 = b1/l11;
          double v2 = (b2-l21*v1)/l22;
          x2 = v2/l22;
          x1 = (v1-l21*x2)/l11;
        }
        float ap0 = (float)x1;
        float a0p = (float)x2;
        s[i2][i1] = f[i2][i1]
                  - ap0*f[i2][i1-1]
                  - a0p*f[i2-1][i1];
      }
    }
    _rgfSmooth.apply0X(s,t);
    _rgfSmooth.applyX0(t,g);
  }

  /**
   * Applies local prediction-error (spectal whitening) filters.
   * The input and output arrays f and g can be the same array.
   * @param f the input array.
   * @param g the output array.
   */
  public void whiten(float[][][] f, float[][][] g) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] r000 = new float[n3][n2][n1];
    float[][][] rpm0 = new float[n3][n2][n1];
    float[][][] rp0m = new float[n3][n2][n1];
    float[][][] r0pm = new float[n3][n2][n1];
    float[][][] rp00 = new float[n3][n2][n1];
    float[][][] r0p0 = new float[n3][n2][n1];
    float[][][] r00p = new float[n3][n2][n1];
    float[][][] s = rp00;
    float[][][] t = r0p0;
    _lcfSymmetric.setInputs(f,f);
    _lcfSymmetric.correlate( 0, 0, 0,r000);
    _lcfSymmetric.correlate( 1,-1, 0,rpm0);
    _lcfSymmetric.correlate( 1, 0,-1,rp0m);
    _lcfSymmetric.correlate( 0, 1,-1,r0pm);
    _lcfSymmetric.correlate( 1, 0, 0,rp00);
    _lcfSymmetric.correlate( 0, 1, 0,r0p0);
    _lcfSymmetric.correlate( 0, 0, 1,r00p);
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        s[i3][i2][0] = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      for (int i1=0; i1<n1; ++i1)
        s[i3][0][i1] = 0.0f;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        s[0][i2][i1] = 0.0f;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          double b1 = rp00[i3][i2][i1];
          double b2 = r0p0[i3][i2][i1];
          double b3 = r00p[i3][i2][i1];
          double a11 = r000[i3][i2][i1];
          double a21 = rpm0[i3][i2][i1];
          double a31 = rp0m[i3][i2][i1];
          double a32 = r0pm[i3][i2][i1];
          double a22 = a11;
          double a33 = a11;
          double x1 = 0.0;
          double x2 = 0.0;
          double x3 = 0.0;
          double l11 = sqrt(a11);
          double l21 = a21/l11;
          double l31 = a31/l11;
          double d22 = a22-l21*l21;
          if (d22>0.0) {
            double l22 = sqrt(d22);
            double l32 = (a32-l31*l21)/l22;
            double d33 = a33-l31*l31-l32*l32;
            if (d33>0.0) {
              double l33 = sqrt(d33);
              double v1 = b1/l11;
              double v2 = (b2-l21*v1)/l22;
              double v3 = (b3-l31*v1-l32*v2)/l33;
              x3 = v3/l33;
              x2 = (v2-l32*x3)/l22;
              x1 = (v1-l21*x2-l31*x3)/l11;
            }
          }
          float ap00 = (float)x1;
          float a0p0 = (float)x2;
          float a00p = (float)x3;
          s[i3][i2][i1] = f[i3][i2][i1]
                        - ap00*f[i3][i2][i1-1]
                        - a0p0*f[i3][i2-1][i1]
                        - a00p*f[i3-1][i2][i1];
        }
      }
    }
    _rgfSmooth.apply0XX(s,t);
    _rgfSmooth.applyX0X(t,s);
    _rgfSmooth.applyXX0(s,g);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private LocalCorrelationFilter _lcfSimple;
  private LocalCorrelationFilter _lcfSymmetric;
  private RecursiveGaussianFilter _rgfSmooth;
  private SincInterpolator _si;
  private double _sigma1,_sigma2,_sigma3;

  private void findShifts(
    int dim, int min, int max, float[][] f, float[][] g, float[][] u) 
  {
    int n1 = f[0].length;
    int n2 = f.length;

    // Default shifts are zero.
    Array.zero(u);

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][][] c = new float[3][n2][n1];

    // Array for current correlation maximum values.
    float[][] cmax = new float[n2][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = (dim==1)?min:0;
    int lag2 = (dim==2)?min:0;
    lcf.correlate(lag1,lag2,c[1]);
    lcf.normalize(lag1,lag2,c[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[][] ca = (lag>min)?c[(i+0)%3]:c[(i+2)%3];
      float[][] cb =           c[(i+1)%3];
      float[][] cc = (lag<max)?c[(i+2)%3]:c[(i+0)%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = (dim==1)?lag+1:0;
        lag2 = (dim==2)?lag+1:0;
        lcf.correlate(lag1,lag2,cc);
        lcf.normalize(lag1,lag2,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, update the correlation maximum value and displacement
      // using quadratic interpolation of three correlation values.
      for (int i2=0; i2<n2; ++i2) {
        float[] ca2 = ca[i2];
        float[] cb2 = cb[i2];
        float[] cc2 = cc[i2];
        for (int i1=0; i1<n1; ++i1) {
          float ai = ca2[i1];
          float bi = cb2[i1];
          float ci = cc2[i1];
          if (bi>=ai && bi>=ci) {
            double c0 = bi;
            double c1 = 0.5*(ci-ai);
            double c2 = 0.5*(ci+ai)-bi;
            double up = (c2<0.0)?-0.5*c1/c2:0.0;
            double cp = c0+up*(c1+up*c2);
            if (cp>cmax[i2][i1]) {
              cmax[i2][i1] = (float)cp;
              u[i2][i1] = (float)(lag+up);
            }
          }
        }
      }
    }
  }

  private void findShifts(
    int dim, int min, int max, float[][][] f, float[][][] g, float[][][] u) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // Default shifts are zero.
    Array.zero(u);

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][][][] c = new float[3][n3][n2][n1];

    // Array for current correlation maximum values.
    float[][][] cmax = new float[n3][n2][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = (dim==1)?min:0;
    int lag2 = (dim==2)?min:0;
    int lag3 = (dim==3)?min:0;
    lcf.correlate(lag1,lag2,lag3,c[1]);
    lcf.normalize(lag1,lag2,lag3,c[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[][][] ca = (lag>min)?c[(i+0)%3]:c[(i+2)%3];
      float[][][] cb =           c[(i+1)%3];
      float[][][] cc = (lag<max)?c[(i+2)%3]:c[(i+0)%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = (dim==1)?lag+1:0;
        lag2 = (dim==2)?lag+1:0;
        lag3 = (dim==3)?lag+1:0;
        lcf.correlate(lag1,lag2,lag3,cc);
        lcf.normalize(lag1,lag2,lag3,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, update the correlation maximum value and displacement
      // using quadratic interpolation of three correlation values.
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] ca32 = ca[i3][i2];
          float[] cb32 = cb[i3][i2];
          float[] cc32 = cc[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float ai = ca32[i1];
            float bi = cb32[i1];
            float ci = cc32[i1];
            if (bi>=ai && bi>=ci) {
              double c0 = bi;
              double c1 = 0.5*(ci-ai);
              double c2 = 0.5*(ci+ai)-bi;
              double up = (c2<0.0)?-0.5*c1/c2:0.0;
              double cp = c0+up*(c1+up*c2);
              if (cp>cmax[i3][i2][i1]) {
                cmax[i3][i2][i1] = (float)cp;
                u[i3][i2][i1] = (float)(lag+up);
              }
            }
          }
        }
      }
    }
  }
}
