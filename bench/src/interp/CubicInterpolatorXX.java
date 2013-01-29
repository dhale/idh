/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

import util.SimplexSolver;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING ONLY!
import java.awt.*;
import javax.swing.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.dsp.*;

/**
 * Piecewise cubic interpolation of a function y(x).
 * <p>
 * Piecewise cubic interpolators differ in the method they use to compute
 * slopes y'(x) at specified x (knots). The classic cubic spline computes the
 * slopes to obtain a continuous second derivative at the knots. These splines
 * often yield unacceptable wiggliness (overshoot) between the knots. A linear
 * spline yields no overshoot, but has discontinuous first (and higher)
 * derivatives. A monotonic spline has continuous first derivatives and yields
 * monotonic interpolation (with no overshoot) where function values at the
 * knots are monotonically increasing or decreasing.
 * <p>
 * For x outside the range of values specified when an interpolator was
 * constructed, the interpolator <em>extrapolates</em> using the cubic
 * polynomial corresponding to the knot nearest to x.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.20
 */
public class CubicInterpolatorXX {

  /**
   * The method used to compute 1st derivatives y'(x).
   */
  public enum Method {
    /**
     * The interpolated y(x) is continuous, but has discontinuous 1st and
     * higher derivatives. This method is equivalent to (though less efficient
     * than) simple piecewise linear interpolation.
     */
    LINEAR,
    /**
     * The interpolated y(x) is continuous with continuous 1st derivative, but
     * may have discontinuous 2nd and higher-order derivatives. This method
     * preserves monotonicity. In intervals where specified y(x) are
     * monotonic, the interpolated values y(x) are also monotonic.
     */
    MONOTONIC,
    /**
     * The interpolated y(x) is continuous with continuous 1st and 2nd
     * derivatives, but may have discontinuous 3rd and higher-order
     * derivatives.
     */
    SPLINE,
  }

  /**
   * Constructs an interpolator with default method monotonic.
   * @param x array of values at which y(x) are specified.
   *  These values must be monotonically increasing or decreasing, 
   *  with no equal values. (In other words, the array must be 
   *  monotonic-definite.)
   * @param y array of function values y(x).
   */
  public CubicInterpolatorXX(float[] x, float[] y) {
    this(Method.MONOTONIC,x,y);
  }

  /**
   * Constructs an interpolator.
   * @param method interpolation method: LINEAR, MONOTONIC, or SPLINE.
   * @param x array of values at which y(x) are specified.
   *  These values must be monotonically increasing or decreasing, 
   *  with no equal values. (In other words, the array must be 
   *  monotonic-definite.)
   * @param y array of function values y(x).
   */
  public CubicInterpolatorXX(Method method, float[] x, float[] y) {
    this(method,x.length,x,y);
  }

  /**
   * Constructs an interpolator.
   * @param method interpolation method: LINEAR, MONOTONIC, or SPLINE.
   * @param n number of x and y(x) values specified.
   * @param x array[n] of values at which y(x) are specified.
   *  These values must be monotonically increasing or decreasing, 
   *  with no equal values. (In other words, the array must be 
   *  monotonic-definite.)
   * @param y array[n] of function values y(x).
   */
  public CubicInterpolatorXX(Method method, int n, float[] x, float[] y) {
    Check.argument(isMonotonic(x), "array x is monotonic");
    _x = copy(n,x);
    _y = copy(n,y);
    _y1 = new float[n];
    _y2 = new float[n];
    _y3 = new float[n];
    if (n==1) { // if only one knot, then constant
      // all derivatives are zero
    } else if (n==2) { // else if only two knots, then linear
      _y1[0] = _y1[1] = (_y[1]-_y[0])/(_x[1]-_x[0]);
    } else {
      if (method==Method.LINEAR) {
        initLinear(_x,_y,_y1);
      } else if (method==Method.MONOTONIC) {
        initMonotonic(_x,_y,_y1);
        compute2ndAnd3rdDerivatives(_x,_y,_y1,_y2,_y3);
      } else if (method==Method.SPLINE) {
        initSpline(_x,_y,_y1);
        compute2ndAnd3rdDerivatives(_x,_y,_y1,_y2,_y3);
      } else {
        throw new IllegalArgumentException("unknown method");
      }
    }
  }

  /**
   * Constructs an interpolator with specified 1st derivatives y'(x).
   * @param x array of values at which y(x) are specified.
   *  These values must be monotonically increasing or decreasing, 
   *  with no equal values. (In other words, the array must be 
   *  monotonic-definite.)
   * @param y array of function values y(x).
   * @param y1 array of 1st derivatives y'(x).
   */
  public CubicInterpolatorXX(float[] x, float[] y, float[] y1) {
    Check.argument(isMonotonic(x), "array x is monotonic");
    int n = x.length;
    _x = copy(n,x);
    _y = copy(n,y);
    _y1 = copy(n,y1);
    _y2 = new float[n];
    _y3 = new float[n];
    compute2ndAnd3rdDerivatives(_x,_y,_y1,_y2,_y3);
  }

  /**
   * Interpolates a function value y(x).
   * Same as {@link #interpolate0(float)}.
   * @param x value at which to interpolate.
   * @return interpolated function value y(x).
   */
  public float interpolate(float x) {
    return interpolate0(x);
  }

  /**
   * Interpolates a function value y(x).
   * @param x value at which to interpolate.
   * @return interpolated function value y(x).
   */
  public float interpolate0(float x) {
    int i = index(x);
    float dx = x-_x[i];
    return _y[i]+dx*(_y1[i]+dx*(_y2[i]*FLT_O2+dx*(_y3[i]*FLT_O6)));
  }

  /**
   * Interpolates the first derivative y'(x).
   * @param x value at which to interpolate.
   * @return interpolated first derivative y'(x).
   */
  public float interpolate1(float x) {
    int i = index(x);
    float dx = x-_x[i];
    return _y1[i]+dx*(_y2[i]+dx*(_y3[i]*FLT_O2));
  }

  /**
   * Interpolates the second derivative y''(x).
   * @param x value at which to interpolate.
   * @return interpolated second derivative y''(x).
   */
  public float interpolate2(float x) {
    int i = index(x);
    float dx = x-_x[i];
    return _y2[i]+dx*_y3[i];
  }

  /**
   * Interpolates the third derivative y'''(x).
   * @param x value at which to interpolate.
   * @return interpolated third derivative y'''(x).
   */
  public float interpolate3(float x) {
    int i = index(x);
    return _y3[i];
  }


  /**
   * Returns an array of interpolated function values y(x).
   * Same as {@link #interpolate0(float[])}.
   * @param x array of values at which to interpolate.
   * @return array of interpolated function values.
   */
  public float[] interpolate(float[] x) {
    return interpolate0(x);
  }

  /**
   * Returns an array of interpolated function values y(x).
   * @param x array of values at which to interpolate.
   * @return array of interpolated function values.
   */
  public float[] interpolate0(float[] x) {
    float[] y = new float[x.length];
    interpolate0(x.length,x,y);
    return y;
  }


  /**
   * Interpolates an array of function values y(x).
   * Same as {@link #interpolate0(float[],float[])}.
   * @param x array of values at which to interpolate.
   * @param y array of interpolated function values.
   */
  public void interpolate(float[] x, float[] y) {
    interpolate0(x,y);
  }

  /**
   * Interpolates an array of function values y(x).
   * @param x array of values at which to interpolate.
   * @param y array of interpolated function values.
   */
  public void interpolate0(float[] x, float[] y) {
    interpolate0(x.length,x,y);
  }

  /**
   * Returns an array of interpolated first derivatives y'(x).
   * @param x array of values at which to interpolate.
   * @return array of interpolated first derivatives y'(x).
   */
  public float[] interpolate1(float[] x) {
    float[] y = new float[x.length];
    interpolate1(x.length,x,y);
    return y;
  }

  /**
   * Interpolates an array of first derivatives y'(x).
   * @param x array of values at which to interpolate.
   * @param y array of interpolated first derivatives y'(x).
   */
  public void interpolate1(float[] x, float[] y) {
    interpolate1(x.length,x,y);
  }

  /**
   * Returns an array of interpolated second derivatives y''(x).
   * @param x array of values at which to interpolate.
   * @return array of interpolated second derivatives y''(x).
   */
  public float[] interpolate2(float[] x) {
    float[] y = new float[x.length];
    interpolate2(x.length,x,y);
    return y;
  }

  /**
   * Interpolates an array of second derivatives y''(x).
   * @param x array of values at which to interpolate.
   * @param y array of interpolated second derivatives y''(x).
   */
  public void interpolate2(float[] x, float[] y) {
    interpolate2(x.length,x,y);
  }

  /**
   * Returns an array of interpolated third derivatives y'''(x).
   * @param x array of values at which to interpolate.
   * @return array of interpolated third derivatives y'''(x).
   */
  public float[] interpolate3(float[] x) {
    float[] y = new float[x.length];
    interpolate3(x.length,x,y);
    return y;
  }

  /**
   * Interpolates an array of third derivatives y'''(x).
   * @param x array of values at which to interpolate.
   * @param y array of interpolated third derivatives y'''(x).
   */
  public void interpolate3(float[] x, float[] y) {
    interpolate3(x.length,x,y);
  }

  /**
   * Interpolates an array of function values y(x).
   * Same as {@link #interpolate0(int,float[],float[])}.
   * @param n number of values to interpolate.
   * @param x array of values at which to interpolate.
   * @param y array of interpolated function values.
   */
  public void interpolate(int n, float[] x, float[] y) {
    interpolate0(n,x,y);
  }

  /**
   * Interpolates an array of function values y(x).
   * @param n number of values to interpolate.
   * @param x array of values at which to interpolate.
   * @param y array of interpolated function values.
   */
  public void interpolate0(int n, float[] x, float[] y) {
    int[] js = {0};
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      int j = index(xi,js);
      float dx = xi-_x[j];
      y[i] = _y[j]+dx*(_y1[j]+dx*(_y2[j]*FLT_O2+dx*(_y3[j]*FLT_O6)));
    }
  }

  /**
   * Interpolates an array of first derivatives y'(x).
   * @param n number of derivatives to interpolate.
   * @param x array of values at which to interpolate.
   * @param y array of interpolated first derivatives y'(x).
   */
  public void interpolate1(int n, float[] x, float[] y) {
    int[] js = {0};
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      int j = index(xi,js);
      float dx = xi-_x[j];
      y[i] = _y1[j]+dx*(_y2[j]+dx*(_y3[j]*FLT_O2));
    }
  }

  /**
   * Interpolates an array of second derivatives y''(x).
   * @param n number of derivatives to interpolate.
   * @param x array of values at which to interpolate.
   * @param y array of interpolated second derivatives y''(x).
   */
  public void interpolate2(int n, float[] x, float[] y) {
    int[] js = {0};
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      int j = index(xi,js);
      float dx = xi-_x[j];
      y[i] = _y2[j]+dx*_y3[j];
    }
  }

  /**
   * Interpolates an array of third derivatives y'''(x).
   * @param n number of derivatives to interpolate.
   * @param x array of values at which to interpolate.
   * @param y array of interpolated third derivatives y'''(x).
   */
  public void interpolate3(int n, float[] x, float[] y) {
    int[] js = {0};
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      int j = index(xi,js);
      y[i] = _y3[j];
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float FLT_O2 = 1.0f/2.0f;
  private static final float FLT_O6 = 1.0f/6.0f;

  private float[] _x,_y,_y1,_y2,_y3;
  private int _index; // index from most recent interpolation

  private int index(float x) {
    int index = binarySearch(_x,x,_index);
    if (index<0) 
      index = (index<-1)?-2-index:0;
    _index = index;
    return index;
  }

  private int index(float x, int[] i) {
    int index = binarySearch(_x,x,i[0]);
    if (index<0) 
      index = (index<-1)?-2-index:0;
    i[0] = index;
    return index;
  }

  /**
   * Computes y2 and y3 from x, y0, and y1.
   */
  private static void compute2ndAnd3rdDerivatives(
    float[] x, float[] y, float[] y1, float[] y2, float[] y3) 
  {
    int n = x.length;
    for (int i=0; i<n-1; ++i) {
      float h2 = x[i+1]-x[i];
      float del2 = (y[i+1]-y[i])/h2;
      float divdf3 = y1[i]+y1[i+1]-2.0f*del2;
      y2[i] = 2.0f*(del2-y1[i]-divdf3)/h2;
      y3[i] = (divdf3/h2)*(6.0f/h2);
    }
    y2[n-1] = y2[n-2]+(x[n-1]-x[n-2])*y3[n-2];
    y3[n-1] = y3[n-2];
  }

  /**
   * Computes cubic interpolation coefficients for linear interpolation.
   */
  private static void initLinear(float[] x, float[] y, float[] y1) {
    int n = x.length;
    for (int i=0; i<n-1; ++i)
      y1[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
    y1[n-1] = y1[n-2];
  }

  /**
   * Computes 1st derivatives via the Fritsch-Carlson method, 
   * which preserves monotonicity.
   * <p>
   * The Fritsch-Carlson method yields continuous 1st derivatives, but 2nd
   * and 3rd derivatives are discontinuous.  The method will yield a
   * monotonic interpolant for monotonic data.  1st derivatives are set to
   * zero wherever first divided differences change sign.
   * <p>
   * For more information, see Fritsch, F. N., and Carlson, R. E., 1980, 
   * Monotone piecewise cubic interpolation:  SIAM J. Numer. Anal., v. 17,
   * n. 2, p. 238-246.
   * <p>
   * Also, see the book by Kahaner, D., Moler, C., and Nash, S., 1989, 
   * Numerical Methods and Software, Prentice Hall.  This function was 
   * derived from SUBROUTINE PCHEZ contained on the diskette that comes 
   * with the book.
   */
  private static void initMonotonic(float[] x, float[] y, float[] y1) {
    int n = x.length;
    
    // Derivatives at left and right ends.
    y1[0] = leftShapePreserving(x,y);
    y1[n-1] = rightShapePreserving(x,y);

    // For all interior knots, ...
    for (int i=1; i<n-1; ++i) {

      // Compute intervals and slopes.
      float h1 = x[i]-x[i-1];
      float h2 = x[i+1]-x[i];
      float hsum = h1+h2;
      float del1 = (y[i]-y[i-1])/h1;
      float del2 = (y[i+1]-y[i])/h2;

      // If a local extremum at this knot, zero derivative.
      if (del1*del2<=0.0f) {
        y1[i] = 0.0f;
      }
      
      // Otherwise, use Butland's formula:
      //      3*(h1+h2)*del1*del2 
      // -------------------------------
      // ((2*h1+h2)*del1+(h1+2*h2)*del2)
      // computed as follows to reduce rounding errors.
      else {
        float dmin = min(abs(del1),abs(del2));
        float dmax = max(abs(del1),abs(del2));
        float drat1 = del1/dmax;
        float drat2 = del2/dmax;
        float hsum3 = hsum+hsum+hsum;
        float w1 = (hsum+h1)/hsum3;
        float w2 = (hsum+h2)/hsum3;
        y1[i] = dmin/(w1*drat1+w2*drat2);
      }
    }
  }

  /**
   * Computes cubic spline interpolation coefficients for interpolation 
   * with continuous second derivatives.
   */
  private static void initSpline(float[] x, float[] y, float[] y1) {
    int n = x.length;
    
    // Derivatives at left and right ends.
    float y1l = leftShapePreserving(x,y);
    float y1r = rightShapePreserving(x,y);
    
    // Compute tridiagonal system coefficients and right-hand-side.
    // a = lower-diagonal, b = diagonal = 2, c = upper-diagonal = 1-a
    // r = right-hand side, w = work array
    float[] a = new float[n];
    float[] r = new float[n];
    float[] w = new float[n];
    a[0] = 1.0f; // = 1-c[0], because c[0] = 0
    r[0] = 2.0f*y1l;
    for (int i=1; i<n-1; ++i) {
      float h1 = x[i]-x[i-1];
      float h2 = x[i+1]-x[i];
      float del1 = (y[i]-y[i-1])/h1;
      float del2 = (y[i+1]-y[i])/h2;
      float alpha = h2/(h1+h2);
      a[i] = alpha;
      r[i] = 3.0f*(alpha*del1+(1.0f-alpha)*del2);
    }
    a[n-1] = 0.0f;
    r[n-1] = 2.0f*y1r;
    
    // Solve tridiagonal system for slopes.
    float t = 2.0f; // = b[0]
    y1[0] = r[0]/t;
    for (int i=1; i<n; ++i) {
      w[i] = (1.0f-a[i-1])/t; // = c[i-1]/t
      t = 2.0f-a[i]*w[i]; // = b[i]-a[i]*w[i]
      y1[i] = (r[i]-a[i]*y1[i-1])/t;
    }
    for (int i=n-2; i>=0; --i) 
      y1[i] -= w[i+1]*y1[i+1];
  }

  /**
   * Left-end slope via shape-preserving 3-point formula.
   */
  private static float leftShapePreserving(float[] x, float[] y) {
    float h1 = x[1]-x[0];
    float h2 = x[2]-x[1];
    float hsum = h1+h2;
    float del1 = (y[1]-y[0])/h1;
    float del2 = (y[2]-y[1])/h2;
    float w1 = (h1+hsum)/hsum;
    float w2 = -h1/hsum;
    float sleft = w1*del1+w2*del2;
    if (sleft*del1<=0.0f) {
      sleft = 0.0f;
    } else if (del1*del2<0.0f) {
      float dmax = 3.0f*del1;
      if (abs(sleft)>abs(dmax)) 
        sleft = dmax;
    }
    return sleft;
  }

  /**
   * Right-end slope via shape-preserving 3-point formula.
   */
  private static float rightShapePreserving(float[] x, float[] y) {
    int n = x.length;
    float h1 = x[n-2]-x[n-3];
    float h2 = x[n-1]-x[n-2];
    float hsum = h1+h2;
    float del1 = (y[n-2]-y[n-3])/h1;
    float del2 = (y[n-1]-y[n-2])/h2;
    float w1 = -h2/hsum;
    float w2 = (h2+hsum)/hsum;
    float sright = w1*del1+w2*del2;
    if (sright*del2<=0.0f) {
      sright = 0.0f;
    } else if (del1*del2<0.0f) {
      float dmax = 3.0f*del2;
      if (abs(sright)>abs(dmax)) 
        sright = dmax;
    }
    return sright;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
    public void run() {
     test1();
     test2();
     test3();
     test4();
     //testMono();
    }});
  }
  private static void test1() {
    float[] x = {0.000f, 1.000f, 2.000f, 3.000f, 4.000f};
    float[] y = {0.000f, 0.000f, 1.000f, 1.000f, 2.000f};
    testMethods(x,y);
  }
  private static void test2() {
    float[] x = {0.0f, 2.0f, 3.0f, 5.0f, 6.0f, 8.0f,
                 9.0f, 11.0f, 12.0f, 14.0f, 15.0f};
    float[] y = {10.0f, 10.0f, 10.0f, 10.0f, 10.0f, 10.0f, 
                 10.5f, 15.0f, 50.0f, 60.0f, 85.0f};
    testMethods(x,y);
  }
  private static void test3() {
    float[] x = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 4.5f,
                 6.0f, 7.0f, 7.3f, 9.0f, 10.0f, 11.0f};
    float[] y = {0.0f, 1.0f, 4.8f, 6.0f, 8.0f, 13.0f,
                 14.0f, 15.5f, 18.0f, 19.0f, 23.0f, 24.1f};
    testMethods(x,y);
  }
  private static void test4() {
    float[] x = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 4.5f,
                 6.0f, 7.0f, 7.3f, 9.0f, 10.0f, 11.0f};
    float[] y = {0.0f, 1.0f, 4.8f, 6.0f, 5.0f, 13.0f,
                 14.0f, 15.5f, 18.0f, 17.0f, 23.0f, 24.1f};
    testMethods(x,y);
  }
  private static void testMethods(float[] x, float[] y) {
    CubicInterpolatorXX.Method[] methods = {
      CubicInterpolatorXX.Method.LINEAR,
      CubicInterpolatorXX.Method.MONOTONIC,
      CubicInterpolatorXX.Method.SPLINE
    };
    Color[] colors = {
      Color.RED,
      Color.GREEN,
      Color.BLUE,
    };
    int nm = methods.length;
    float xmin = min(x);
    float xmax = max(x);
    int nx = 1001;
    float dx = (xmax-xmin)/(nx-1);
    float fx = xmin;
    float[] xi = rampfloat(fx,dx,nx);
    float[] yi = new float[nx];
    SimplePlot sp = new SimplePlot();
    for (int im=0; im<nm; ++im) {
      CubicInterpolatorXX.Method method = methods[im];
      CubicInterpolatorXX ci = new CubicInterpolatorXX(method,x,y);
      ci.interpolate1(xi,yi);
      PointsView pv = sp.addPoints(xi,yi);
      pv.setLineColor(colors[im]);
    }
  }
  private static void testMono() {
    float[] x = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 4.5f,
                 6.0f, 7.0f, 7.3f, 9.0f, 10.0f, 11.0f};
    float[] y = {0.0f, 1.0f, 4.8f, 6.0f, 8.0f, 13.0f,
                 14.0f, 15.5f, 18.0f, 19.0f, 23.0f, 24.1f};
    Color[] colors = {
      Color.RED,
      Color.BLUE,
    };
    float xmin = min(x);
    float xmax = max(x);
    int nx = 1001;
    float dx = (xmax-xmin)/(nx-1);
    float fx = xmin;
    float[] xi = rampfloat(fx,dx,nx);
    float[] yi = new float[nx];
    float[] zi = new float[nx];
    SimplePlot sp = new SimplePlot();
    CubicInterpolatorXX ci = new CubicInterpolatorXX(x,y);
    ci.interpolate0(xi,yi);
    ci.interpolate1(xi,zi);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(20.0);
    ref.apply(zi,zi);
    zi[0] = yi[0];
    for (int ix=1; ix<nx; ++ix)
      zi[ix] = zi[ix-1]+dx*zi[ix];
    PointsView pv = sp.addPoints(xi,yi);
    pv.setLineColor(Color.RED);
    pv = sp.addPoints(xi,zi);
    pv.setLineColor(Color.BLUE);
    pv = sp.addPoints(x,y);
    pv.setLineStyle(PointsView.Line.NONE);
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
  }
}
