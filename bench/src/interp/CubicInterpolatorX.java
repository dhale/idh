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
 * @version 2012.12.27
 */
public class CubicInterpolatorX {

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

    /**
     * New method.
     */
    WOLBERG
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
  public CubicInterpolatorX(float[] x, float[] y, float[] y1) {
    Check.argument(isMonotonic(x), "array x is monotonic");
    int n = x.length;
    _xd = new float[n];
    _yd = new float[n][4];
    for (int i=0; i<n; ++i) {
      _xd[i] = x[i];
      _yd[i][0] = y[i];
      _yd[i][1] = y1[i];
    }
    compute2ndAnd3rdDerivatives(_xd,_yd);
  }

  /**
   * Constructs an interpolator with default method monotonic.
   * @param x array of values at which y(x) are specified.
   *  These values must be monotonically increasing or decreasing, 
   *  with no equal values. (In other words, the array must be 
   *  monotonic-definite.)
   * @param y array of function values y(x).
   */
  public CubicInterpolatorX(float[] x, float[] y) {
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
  public CubicInterpolatorX(Method method, float[] x, float[] y) {
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
  public CubicInterpolatorX(Method method, int n, float[] x, float[] y) {
    Check.argument(isMonotonic(x), "array x is monotonic");
    _xd = new float[n];
    _yd = new float[n][4];
    for (int i=0; i<n; ++i) {
      _xd[i] = x[i];
      _yd[i][0] = y[i];
    }
    if (method==Method.LINEAR) {
      initLinear(n,_xd,_yd);
    } else if (method==Method.MONOTONIC) {
      initMonotonic(n,_xd,_yd);
    } else if (method==Method.SPLINE) {
      initSpline(n,_xd,_yd);
    } else if (method==Method.WOLBERG) {
      initWolberg(n,_xd,_yd);
    } else {
      assert false;
    }
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
    return interpolate0(x-_xd[i],_yd[i]);
  }

  /**
   * Interpolates the first derivative y'(x).
   * @param x value at which to interpolate.
   * @return interpolated first derivative y'(x).
   */
  public float interpolate1(float x) {
    int i = index(x);
    return interpolate1(x-_xd[i],_yd[i]);
  }

  /**
   * Interpolates the second derivative y''(x).
   * @param x value at which to interpolate.
   * @return interpolated second derivative y''(x).
   */
  public float interpolate2(float x) {
    int i = index(x);
    return interpolate2(x-_xd[i],_yd[i]);
  }

  /**
   * Interpolates the third derivative y'''(x).
   * @param x value at which to interpolate.
   * @return interpolated third derivative y'''(x).
   */
  public float interpolate3(float x) {
    int i = index(x);
    return interpolate3(x-_xd[i],_yd[i]);
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
      int j = index(x[i],_xd,js);
      y[i] = interpolate0(x[i]-_xd[j],_yd[j]);
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
      int j = index(x[i],_xd,js);
      y[i] = interpolate1(x[i]-_xd[j],_yd[j]);
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
      int j = index(x[i],_xd,js);
      y[i] = interpolate2(x[i]-_xd[j],_yd[j]);
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
      int j = index(x[i],_xd,js);
      y[i] = interpolate3(x[i]-_xd[j],_yd[j]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private.

  private float[] _xd; // array[n] of x.
  private float[][] _yd; // array[n][4] of y, y', y'', and y'''.
  private int _index; // index from most recent interpolation

  private int index(float x) {
    int index = binarySearch(_xd,x,_index);
    if (index<0) 
      index = (index<-1)?-2-index:0;
    _index = index;
    return index;
  }

  private static int index(float x, float[] xs, int[] i) {
    int index = binarySearch(xs,x,i[0]);
    if (index<0) 
      index = (index<-1)?-2-index:0;
    i[0] = index;
    return index;
  }

  private static final float FLT_O2 = 1.0f/2.0f;
  private static final float FLT_O6 = 1.0f/6.0f;
  private static float interpolate0(float dx, float[] yd) {
    return yd[0]+dx*(yd[1]+dx*(yd[2]*FLT_O2+dx*(yd[3]*FLT_O6)));
  }
  private static float interpolate1(float dx, float[] yd) {
    return yd[1]+dx*(yd[2]+dx*(yd[3]*FLT_O2));
  }
  private static float interpolate2(float dx, float[] yd) {
    return yd[2]+dx*yd[3];
  }
  private static float interpolate3(float dx, float[] yd) {
    return yd[3];
  }

  /**
   * Computes cubic interpolation coefficients for linear interpolation.
   */
  private static void initLinear(int n, float[] x, float[][] y) {

    // If n=1, then use constant interpolation.
    if (n==1) {
      y[0][1] = 0.0f;
      y[0][2] = 0.0f;
      y[0][3] = 0.0f;
      return;
    }

    // Compute slopes.
    for (int i=0; i<n-1; ++i) {
      y[i][1] = (y[i+1][0]-y[i][0])/(x[i+1]-x[i]);
      y[i][2] = y[i][3] = 0.0f;
    }
    y[n-1][1] = y[n-2][1];
    y[n-1][2] = y[n-1][3] = 0.0f;
  }

  /**
   * Computes cubic interpolation coefficients via the Fritsch-Carlson 
   * method, which preserves monotonicity.
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
  private static void initMonotonic(int n, float[] x, float[][] y) {

    // If n=1, then use constant interpolation.
    if (n==1) {
      y[0][1] = 0.0f;
      y[0][2] = 0.0f;
      y[0][3] = 0.0f;
      return;

    // Else, if n=2, then use linear interpolation.
    } else if (n==2) {
      y[0][1] = y[1][1] = (y[1][0]-y[0][0])/(x[1]-x[0]);
      y[0][2] = y[1][2] = 0.0f;
      y[0][3] = y[1][3] = 0.0f;
      return;
    }

    // Set left end derivative via shape-preserving 3-point formula.
    float h1 = x[1]-x[0];
    float h2 = x[2]-x[1];
    float hsum = h1+h2;
    float del1 = (y[1][0]-y[0][0])/h1;
    float del2 = (y[2][0]-y[1][0])/h2;
    float w1 = (h1+hsum)/hsum;
    float w2 = -h1/hsum;
    y[0][1] = w1*del1+w2*del2;
    if (y[0][1]*del1<=0.0f)
      y[0][1] = 0.0f;
    else if (del1*del2<0.0f) {
      float dmax = 3.0f*del1;
      if (abs(y[0][1])>abs(dmax)) y[0][1] = dmax;
    }

    // Loop over interior points.
    for (int i=1; i<n-1; ++i) {

      // Compute intervals and slopes.
      h1 = x[i]-x[i-1];
      h2 = x[i+1]-x[i];
      hsum = h1+h2;
      del1 = (y[i][0]-y[i-1][0])/h1;
      del2 = (y[i+1][0]-y[i][0])/h2;

      // If not strictly monotonic, zero derivative.
      if (del1*del2<=0.0f) {
        y[i][1] = 0.0f;
      
      // Else, if strictly monotonic, use Butland's formula:
      //      3*(h1+h2)*del1*del2 
      // -------------------------------
      // ((2*h1+h2)*del1+(h1+2*h2)*del2)
      // computed as follows to avoid roundoff error
      } else {
        float dmin = min(abs(del1),abs(del2));
        float dmax = max(abs(del1),abs(del2));
        float drat1 = del1/dmax;
        float drat2 = del2/dmax;
        float hsum3 = hsum+hsum+hsum;
        w1 = (hsum+h1)/hsum3;
        w2 = (hsum+h2)/hsum3;
        y[i][1] = dmin/(w1*drat1+w2*drat2);
      }
    }

    // Set right end derivative via shape-preserving 3-point formula.
    w1 = -h2/hsum;
    w2 = (h2+hsum)/hsum;
    y[n-1][1] = w1*del1+w2*del2;
    if (y[n-1][1]*del2<=0.0f)
      y[n-1][1] = 0.0f;
    else if (del1*del2<0.0f) {
      float dmax = 3.0f*del2;
      if (abs(y[n-1][1])>abs(dmax)) y[n-1][1] = dmax;
    }

    compute2ndAnd3rdDerivatives(x,y);
  }

  /**
   * Computes cubic spline interpolation coefficients for interpolation 
   * with continuous second derivatives.
   */
  private static void initSpline(int n, float[] x, float[][] y) {

    // If n=1, then use constant interpolation.
    if (n==1) {
      y[0][1] = 0.0f;
      y[0][2] = 0.0f;
      y[0][3] = 0.0f;
      return;

    // Else, if n=2, then use linear interpolation.
    } else if (n==2) {
      y[0][1] = y[1][1] = (y[1][0]-y[0][0])/(x[1]-x[0]);
      y[0][2] = y[1][2] = 0.0f;
      y[0][3] = y[1][3] = 0.0f;
      return;
    }
    
    // Set left end derivative via shape-preserving 3-point formula.
    float h1 = x[1]-x[0];
    float h2 = x[2]-x[1];
    float hsum = h1+h2;
    float del1 = (y[1][0]-y[0][0])/h1;
    float del2 = (y[2][0]-y[1][0])/h2;
    float w1 = (h1+hsum)/hsum;
    float w2 = -h1/hsum;
    float sleft = w1*del1+w2*del2;
    if (sleft*del1<=0.0f)
      sleft = 0.0f;
    else if (del1*del2<0.0f) {
      float dmax = 3.0f*del1;
      if (abs(sleft)>abs(dmax)) sleft = dmax;
    }

    // Set right end derivative via shape-preserving 3-point formula.
    h1 = x[n-2]-x[n-3];
    h2 = x[n-1]-x[n-2];
    hsum = h1+h2;
    del1 = (y[n-2][0]-y[n-3][0])/h1;
    del2 = (y[n-1][0]-y[n-2][0])/h2;
    w1 = -h2/hsum;
    w2 = (h2+hsum)/hsum;
    float sright = w1*del1+w2*del2;
    if (sright*del2<=0.0f)
      sright = 0.0f;
    else if (del1*del2<0.0f) {
      float dmax = 3.0f*del2;
      if (abs(sright)>abs(dmax)) sright = dmax;
    }
    
    // Compute tridiagonal system coefficients and right-hand-side.
    float[] work = new float[n];
    work[0] = 1.0f;
    y[0][2] = 2.0f*sleft;
    for (int i=1; i<n-1; ++i) {
      h1 = x[i]-x[i-1];
      h2 = x[i+1]-x[i];
      del1 = (y[i][0]-y[i-1][0])/h1;
      del2 = (y[i+1][0]-y[i][0])/h2;
      float alpha = h2/(h1+h2);
      work[i] = alpha;
      y[i][2] = 3.0f*(alpha*del1+(1.0f-alpha)*del2);
    }
    work[n-1] = 0.0f;
    y[n-1][2] = 2.0f*sright;
    
    // Solve tridiagonal system for slopes.
    float t = 2.0f;
    y[0][1] = y[0][2]/t;
    for (int i=1; i<n; ++i) {
      y[i][3] = (1.0f-work[i-1])/t;
      t = 2.0f-work[i]*y[i][3];
      y[i][1] = (y[i][2]-work[i]*y[i-1][1])/t;
    }
    for (int i=n-2; i>=0; --i) 
      y[i][1] -= y[i+1][3]*y[i+1][1];

    compute2ndAnd3rdDerivatives(x,y);
  }

  private static void initWolberg(int nxy, float[] x, float[][] y) {

    // If nxy=1, then use constant interpolation.
    if (nxy==1) {
      y[0][1] = 0.0f;
      y[0][2] = 0.0f;
      y[0][3] = 0.0f;
      return;

    // Else, if nxy=2, then use linear interpolation.
    } else if (nxy==2) {
      y[0][1] = y[1][1] = (y[1][0]-y[0][0])/(x[1]-x[0]);
      y[0][2] = y[1][2] = 0.0f;
      y[0][3] = y[1][3] = 0.0f;
      return;
    }

    int n = nxy+(nxy-2); // n variables include nxy-2 slack variables
    int m = 4*(nxy-1)+2*(nxy-2); // m constraint equations
    double[][] a = new double[m][n]; // indexed by i,j (row,column)
    double[] b = new double[m]; // indexed by i (row)
    double[] c = new double[n]; // indexed by j (column)

    // Slopes in the nxy-1 intervals between the nxy knots. The derivatives at
    // the knots will be products of non-negative scale factors and these
    // slopes. The scale factors will be the solution to a linear program that
    // we construct below. Because we need nxy derivatives, but have only
    // nxy-1 intervals, we set the last slope equal to the next-to-last slope. 
    double[] p = new double[nxy];
    for (int k=0; k<nxy-1; ++k)
      p[k] = (y[k+1][0]-y[k][0])/(x[k+1]-x[k]);
    p[nxy-1] = p[nxy-2];

    // Monotonicity constraints on the nxy derivatives. We have four such
    // constraints for each of the nxy-1 intervals between the nxy knots.
    for (int k=0,i=0,j=0; k<nxy-1; ++k,i+=4,++j) {
      if (k>=0) 
        continue;
      double pm = (k>0)?p[k-1]:0.0;
      double p0 = p[k];
      double p1 = p[k+1];

      // If interval is on either side of a local extremum, then
      // impose no constraints; both A and b are already zeroed.
      if (pm*p0<0.0 || p0*p1<0.0)
        continue;

      // If slope is zero, then derivatives at each end must be zero.
      // For each end we use two inequality constraints to force the
      // corresponding derivative to be zero.
      if (p0==0.0) {
        a[i+0][j] =  1.0;                       b[i+0] = 0.0;
        a[i+1][j] = -1.0;                       b[i+1] = 0.0;
                           a[i+2][j+1] =  1.0;  b[i+2] = 0.0;
                           a[i+3][j+1] = -1.0;  b[i+3] = 0.0;
      } 
      
      // Otherwise, neither end is an extremum and slope is non-zero, so
      // ensure monotonically increasing or decreasing in this interval.
      else {
        /*
        a[i+0][j] =     p0;  a[i+0][j+1] =    -p1;  b[i+0] = 3.0*p0;
        a[i+1][j] = 2.0*p0;  a[i+1][j+1] =     p1;  b[i+1] = 9.0*p0;
        a[i+2][j] =     p0;  a[i+2][j+1] = 2.0*p1;  b[i+2] = 9.0*p0;
        a[i+3][j] =    -p0;  a[i+3][j+1] =     p1;  b[i+3] = 3.0*p0;
        */
        a[i+0][j] =  1.0;  a[i+0][j+1] = -1.0;  b[i+0] = 3.0*p0;
        a[i+1][j] =  2.0;  a[i+1][j+1] =  1.0;  b[i+1] = 9.0*p0;
        a[i+2][j] =  1.0;  a[i+2][j+1] =  2.0;  b[i+2] = 9.0*p0;
        a[i+3][j] = -1.0;  a[i+3][j+1] =  1.0;  b[i+3] = 3.0*p0;
      }
    }

    // Absolute-value constraints associated with nxy-2 slack variables sk.
    // The slack variables sk equal absolute values of discontinuities in 2nd
    // derivatives at the nxy-2 interior knots. We represent each of these
    // equalities with two inequalities.
    for (int k=0,i=4*(nxy-1),j=0,jsk=nxy; k<nxy-2; ++k,i+=2,++j,++jsk) {
      double d0 = x[k+1]-x[k  ];
      double d1 = x[k+2]-x[k+1];
      double p0 = p[k  ];
      double p1 = p[k+1];
      /*
      double p2 = p[k+2];
      double a0 = p0/d0;
      double a1 = 2.0*(p1/d0+p1/d1);
      double a2 = p2/d1;
      */
      double a0 = 1.0/d0;
      double a1 = 2.0/d0+2.0/d1;
      double a2 = 1.0/d1;
      double bi = 3.0*(p0/d0+p1/d1);
      a[i  ][j+0] =  a0; 
      a[i  ][j+1] =  a1; 
      a[i  ][j+2] =  a2;
      a[i  ][jsk] = -1.0;
      b[i  ] =  bi;
      a[i+1][j+0] = -a0; 
      a[i+1][j+1] = -a1; 
      a[i+1][j+2] = -a2;
      a[i+1][jsk] = -1.0;
      b[i+1] = -bi;
    }

    // Objective function = negated sum of all slack variables sk. The sum is
    // negated because the simple solver will maximize (not minimize) it.
    for (int j=nxy; j<n; ++j)
      c[j] = -1.0;
    dump(a);
    dump(b);
    dump(c);
    SimplexSolver ss = new SimplexSolver();
    double[] z = ss.solve(a,b,c);
    dump(z);
    for (int k=0; k<nxy; ++k)
      //y[k][1] = (float)(z[k]*p[k]);
      y[k][1] = (float)(z[k]);
    dump(y);

    compute2ndAnd3rdDerivatives(x,y);
  }

  private static void compute2ndAnd3rdDerivatives(float[] x, float[][] y) {
    int n = x.length;
    for (int i=0; i<n-1; ++i) {
      float h2 = x[i+1]-x[i];
      float del2 = (y[i+1][0]-y[i][0])/h2;
      float divdf3 = y[i][1]+y[i+1][1]-2.0f*del2;
      y[i][2] = 2.0f*(del2-y[i][1]-divdf3)/h2;
      y[i][3] = (divdf3/h2)*(6.0f/h2);
    }
    y[n-1][2] = y[n-2][2]+(x[n-1]-x[n-2])*y[n-2][3];
    y[n-1][3] = y[n-2][3];
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
    public void run() {
      test1();
      //test2();
      //test3();
    }});
  }
  private static void test1() {
    float[] x = {0.000f, 1.000f, 2.000f, 3.000f, 4.000f};
    float[] y = {0.000f, 0.001f, 1.000f, 1.001f, 2.000f};
    testMethods(x,y);
  }
  private static void test2() {
    float[] x = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 4.5f,
                 6.0f, 7.0f, 7.3f, 9.0f, 10.0f, 11.0f};
    float[] y = {0.0f, 1.0f, 4.8f, 6.0f, 8.0f, 13.0f,
                 14.0f, 15.5f, 18.0f, 19.0f, 23.0f, 24.1f};
    testMethods(x,y);
  }
  private static void test3() {
    float[] x = {0.0f, 2.0f, 3.0f, 5.0f, 6.0f, 8.0f,
                 9.0f, 11.0f, 12.0f, 14.0f, 15.0f};
    //float[] y = {10.0f, 10.001f, 10.002f, 10.003f, 10.004f, 10.005f, 
    //             10.5f, 15.0f, 50.0f, 60.0f, 85.0f};
    float[] y = {10.0f, 10.0f, 10.0f, 10.0f, 10.0f, 10.0f, 
                 10.5f, 15.0f, 50.0f, 60.0f, 85.0f};
    testMethods(x,y);
  }
  private static void testMethods(float[] x, float[] y) {
    CubicInterpolatorX.Method[] methods = {
      CubicInterpolatorX.Method.LINEAR,
      CubicInterpolatorX.Method.MONOTONIC,
      CubicInterpolatorX.Method.WOLBERG,
      CubicInterpolatorX.Method.SPLINE
    };
    Color[] colors = {
      Color.RED,
      Color.GREEN,
      Color.BLUE,
      Color.MAGENTA
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
      CubicInterpolatorX.Method method = methods[im];
      CubicInterpolatorX ci = new CubicInterpolatorX(method,x,y);
      ci.interpolate0(xi,yi);
      PointsView pv = sp.addPoints(xi,yi);
      pv.setLineColor(colors[im]);
    }
  }
}
