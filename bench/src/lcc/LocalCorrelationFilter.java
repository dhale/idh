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
 * Local cross-correlation of two arrays with seamless overlapping windows.
 * Given two input arrays f and g and a specified lag, this filter computes 
 * an output array c of local cross-correlation coefficients, one for each
 * sample in the input arrays f and g.
 * <p>
 * Two types of cross-correlation are implemented. Both types can be 
 * normalized to obtain cross-correlation coefficients with magnitudes 
 * that do not exceed one. The normalization varies, depending on the 
 * type of cross-correlation.
 * <p>
 * <em>Simple</em> cross-correlation computes an array of products 
 * h[j] = f[j]*g[j+lag] and then filters this array of products with a 
 * window. The resulting correlation cfg[k,lag] is not symmetric with 
 * respect to lag; cfg[k,-lag] = cgf[k-lag,lag] != cgf[k,lag]. For
 * simple cross-correlation, normalization scale factors vary with lag
 * and should be applied before picking correlation peaks.
 * <p>
 * <em>Symmetric</em> cross-correlation computes an array of products
 * h[j] = f[j-lag/2]*g[j+lag/2] and therefore requires interpolation
 * between samples for odd lags. (For efficiency, we interpolate the 
 * products h, not the inputs f and g.) The resulting correlation is 
 * symmetric with respect to lag; cfg[k,lag] = cgf[k,-lag]. Moreover,
 * when inputs f and g are the same, each local auto-correlation has a
 * Fourier transform (a power spectrum) that is positive-semidefinite.
 * This property is important in applications such as local prediction 
 * filtering. For symmetric cross-correlation, normalization scale 
 * factors do not vary with lag, and therefore need not be applied 
 * before picking correlation peaks.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.08.11
 */
public class LocalCorrelationFilter {

  /**
   * Cross-correlations windows.
   * The default window is GAUSSIAN.
   */
  enum Window {
    GAUSSIAN,
    RECTANGLE
  };

  /**
   * Cross-correlations types.
   * The default type is SYMMETRIC.
   */
  enum Type {
    SIMPLE,
    SYMMETRIC
  };

  /**
   * Construct a correlation filter with specified parameters.
   * When applied to multi-dimensional arrays, the filter has the
   * same half-width for all dimensions.
   * @param type the correlation type.
   * @param window the correlation window.
   * @param sigma the correlation window half-width; must not be less than 1.
   */
  public LocalCorrelationFilter(Type type, Window window, double sigma) {
    this(type,window,sigma,sigma,sigma);
  }

  /**
   * Construct a correlation filter with specified parameters.
   * When applied to multi-dimensional arrays, the filter has half-width 
   * sigma1 for the 1st dimension and half-width sigma2 for 2nd and higher 
   * dimensions.
   * @param type the correlation type.
   * @param window the correlation window.
   * @param sigma1 correlation window half-width for 1st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd and higher 
   *  dimensions; must not be less than 1.
   */
  public LocalCorrelationFilter(
    Type type, Window window, double sigma1, double sigma2) 
  {
    this(type,window,sigma1,sigma2,sigma2);
  }

  /**
   * Construct a correlation filter with specified parameters.
   * When applied to multi-dimensional arrays, the filter has half-width 
   * sigma1 for the 1st dimension, half-width sigma2 for the 2nd dimension,
   * and half-width sigma3 for 3rd and higher dimensions.
   * @param type the correlation type.
   * @param window the correlation window.
   * @param sigma1 correlation window half-width for 1st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd dimension;
   *  must not be less than 1.
   * @param sigma3 correlation window half-width for 3rd and higher 
   * dimensions; must not be less than 1.
   */
  public LocalCorrelationFilter(
    Type type, Window window, double sigma1, double sigma2, double sigma3) 
  {
    Check.argument(sigma1>=1.0,"sigma1>=1.0");
    Check.argument(sigma2>=1.0,"sigma2>=1.0");
    Check.argument(sigma3>=1.0,"sigma3>=1.0");
    _type = type;
    _window = window;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
    _sigma3 = sigma3;
    if (window==Window.GAUSSIAN) {
      _f1 = new GaussianFilter(sigma1);
      _f2 = new GaussianFilter(sigma2);
      _f3 = new GaussianFilter(sigma3);
    } else {
      _f1 = new RectangleFilter(sigma1);
      _f2 = new RectangleFilter(sigma2);
      _f3 = new RectangleFilter(sigma3);
    }
  }

  /**
   * Sets the input arrays to be cross-correlated.
   * The input arrays f and g can be the same array.
   * @param f the input array f; by reference, not copied.
   * @param g the input array g; by reference, not copied.
   */
  public void setInputs(float[] f, float[] g) {
    if (f==null || g==null) {
      _dimension = _n1 = _n2 = _n3 = 0;
      _f = null;
      _g = null;
    } else {
      Check.argument(f.length==g.length,"f.length==g.length");
      _dimension = 1;
      _n1 = f.length;
      _n2 = _n3 = 0;
      _f = new float[1][1][];
      _g = new float[1][1][];
      _f[0][0] = f;
      _g[0][0] = g;
    }
    _s = null;
  }

  /**
   * Sets the input arrays to be cross-correlated.
   * The input arrays f and g can be the same array.
   * @param f the input array f; by reference, not copied.
   * @param g the input array g; by reference, not copied.
   */
  public void setInputs(float[][] f, float[][] g) {
    if (f==null || g==null) {
      _dimension = _n1 = _n2 = _n3 = 0;
      _f = null;
      _g = null;
    } else {
      Check.argument(f[0].length==g[0].length,"f[0].length==g[0].length");
      Check.argument(f.length==g.length,"f.length==g.length");
      Check.argument(Array.isRegular(f),"f is regular");
      Check.argument(Array.isRegular(g),"g is regular");
      _dimension = 2;
      _n1 = f[0].length;
      _n2 = f.length;
      _n3 = 0;
      _f = new float[1][][];
      _g = new float[1][][];
      _f[0] = f;
      _g[0] = g;
    }
    _s = null;
  }

  /**
   * Sets the input arrays to be cross-correlated.
   * The input arrays f and g can be the same array.
   * @param f the input array f; by reference, not copied.
   * @param g the input array g; by reference, not copied.
   */
  public void setInputs(float[][][] f, float[][][] g) {
    if (f==null || g==null) {
      _dimension = _n1 = _n2 = _n3 = 0;
      _f = null;
      _g = null;
    } else {
      Check.argument(
        f[0][0].length==g[0][0].length,"f[0][0].length==g[0][0].length");
      Check.argument(f[0].length==g[0].length,"f[0].length==g[0].length");
      Check.argument(f.length==g.length,"f.length==g.length");
      Check.argument(Array.isRegular(f),"f is regular");
      Check.argument(Array.isRegular(g),"g is regular");
      _dimension = 3;
      _n1 = f[0][0].length;
      _n2 = f[0].length;
      _n3 = f.length;
      _f = f;
      _g = g;
    }
    _s = null;
  }

  /**
   * Correlates the current inputs for the specified lag.
   * @param lag the correlation lag.
   * @param c the output array; cannot be the same as inputs f or g.
   */
  public void correlate(int lag, float[] c) {
    checkDimensions(c);
    correlate(lag,_f[0][0],_g[0][0],c);
  }

  /**
   * Correlates the current inputs for the specified lag.
   * @param lag1 the lag in the 1st dimension.
   * @param lag2 the lag in the 2nd dimension.
   * @param c the output array; cannot be the same as inputs f or g.
   */
  public void correlate(int lag1, int lag2, float[][] c) {
    checkDimensions(c);
    correlate(lag1,lag2,_f[0],_g[0],c);
  }

  /**
   * Correlates the current inputs for the specified lag.
   * @param lag1 the lag in the 1st dimension.
   * @param lag2 the lag in the 2nd dimension.
   * @param lag3 the lag in the 3rd dimension.
   * @param c the output array; cannot be the same as inputs f or g.
   */
  public void correlate(int lag1, int lag2, int lag3, float[][][] c) {
    checkDimensions(c);
    correlate(lag1,lag2,lag3,_f,_g,c);
  }

  /**
   * Normalizes the cross-correlation for a specified lag.
   * @param lag the lag.
   * @param c the cross-correlation to be modified.
   */
  public void normalize(int lag, float[] c) {
    checkDimensions(c);
    if (_s==null)
      updateNormalize();
    int n1 = _n1;
    int l1 = lag;
    if (_type==Type.SIMPLE) {
      float[] sf = _s[0][0][0];
      float[] sg = _s[1][0][0];
      int i1min = max(0,-l1);
      int i1max = min(n1,n1-l1);
      for (int i1=0; i1<i1min; ++i1) {
        c[i1] *= sf[i1]*sg[0];
      }
      for (int i1=i1min; i1<i1max; ++i1) {
        c[i1] *= sf[i1]*sg[i1+l1];
      }
      for (int i1=i1max; i1<n1; ++i1) {
        c[i1] *= sf[i1]*sg[n1-1];
      }
    } else if (_type==Type.SYMMETRIC) {
      float[] s = _s[0][0][0];
      for (int i1=0; i1<n1; ++i1) {
        c[i1] *= s[i1];
      }
    }
  }

  /**
   * Normalizes the cross-correlation for a specified lag.
   * @param lag1 the lag.
   * @param c the cross-correlation to be modified.
   */
  public void normalize(int lag1, int lag2, float[][] c) {
    checkDimensions(c);
    if (_s==null)
      updateNormalize();
    int n1 = _n1;
    int n2 = _n2;
    int l1 = lag1;
    int l2 = lag2;
    if (_type==Type.SIMPLE) {
      float[][] sf = _s[0][0];
      float[][] sg = _s[1][0];
      int i1min = max(0,-l1);
      int i1max = min(n1,n1-l1);
      for (int i2=0; i2<n2; ++i2) {
        float[] c2 = c[i2];
        float[] sf2 = sf[i2];
        float[] sg2 = sg[max(0,min(n2-1,i2+l2))];
        for (int i1=0; i1<i1min; ++i1) {
          c2[i1] *= sf2[i1]*sg2[0];
        }
        for (int i1=i1min; i1<i1max; ++i1) {
          c2[i1] *= sf2[i1]*sg2[i1+l1];
        }
        for (int i1=i1max; i1<n1; ++i1) {
          c2[i1] *= sf2[i1]*sg2[n1-1];
        }
      }
    } else if (_type==Type.SYMMETRIC) {
      float[][] s = _s[0][0];
      for (int i2=0; i2<n2; ++i2) {
        float[] c2 = c[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1) {
          c2[i1] *= s2[i1];
        }
      }
    }
  }

  /**
   * Normalizes the cross-correlation for a specified lag.
   * @param lag the lag.
   * @param c the cross-correlation to be modified.
   */
  public void normalize(int lag1, int lag2, int lag3, float[][][] c) {
    checkDimensions(c);
    if (_s==null)
      updateNormalize();
    int n1 = _n1;
    int n2 = _n2;
    int n3 = _n3;
    int l1 = lag1;
    int l2 = lag2;
    int l3 = lag3;
    if (_type==Type.SIMPLE) {
      float[][][] sf = _s[0];
      float[][][] sg = _s[1];
      int i1min = max(0,-l1);
      int i1max = min(n1,n1-l1);
      for (int i3=0; i3<n3; ++i3) {
        float[][] c3 = c[i3];
        float[][] sf3 = sf[i3];
        float[][] sg3 = sg[max(0,min(n3-1,i3+l3))];
        for (int i2=0; i2<n2; ++i2) {
          float[] c32 = c3[i2];
          float[] sf32 = sf3[i2];
          float[] sg32 = sg3[max(0,min(n2-1,i2+l2))];
          for (int i1=0; i1<i1min; ++i1) {
            c32[i1] *= sf32[i1]*sg32[0];
          }
          for (int i1=i1min; i1<i1max; ++i1) {
            c32[i1] *= sf32[i1]*sg32[i1+l1];
          }
          for (int i1=i1max; i1<n1; ++i1) {
            c32[i1] *= sf32[i1]*sg32[n1-1];
          }
        }
      }
    } else if (_type==Type.SYMMETRIC) {
      float[][][] s = _s[0];
      for (int i3=0; i3<n3; ++i3) {
        float[][] c3 = c[i3];
        float[][] s3 = s[i3];
        for (int i2=0; i2<n2; ++i2) {
          float[] c32 = c3[i2];
          float[] s32 = s3[i2];
          for (int i1=0; i1<n1; ++i1) {
            c32[i1] *= s32[i1];
          }
        }
      }
    }
  }

  /** 
   * Removes bias by subtracting local means from the specified array.
   * @param f the input array.
   * @return the output array, with bias subtracted.
   */
  public float[] unbias(float[] f) {
    int n1 = f.length;
    float[] t = new float[n1];
    _f1.apply(f,t);
    Array.sub(f,t,t);
    return t;
  }

  /** 
   * Removes bias by subtracting local means from the specified array.
   * @param f the input array.
   * @return the output array, with bias subtracted.
   */
  public float[][] unbias(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] t = new float[n2][n1];
    _f1.apply1(f,t);
    _f2.apply2(t,t);
    Array.sub(f,t,t);
    return t;
  }

  /** 
   * Removes bias by subtracting local means from the specified array.
   * @param f the input array.
   * @return the output array, with bias subtracted.
   */
  public float[][][] unbias(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] t = new float[n3][n2][n1];
    _f1.apply1(f,t);
    _f2.apply2(t,t);
    _f3.apply3(t,t);
    Array.sub(f,t,t);
    return t;
  }

  /**
   * Searches for lags for which cross-correlations are maximized.
   * @param min minimum lag
   * @param max maximum lag
   * @param lag output array of lags
   */
  public void findMaxLags(int min, int max, byte[] lag) {
    checkDimension(1);

    // Initialize arrays of lags.
    int n = _n1;
    for (int i=0; i<n; ++i)
      lag[i] = 0;

    // Array for cross-correlations.
    float[] c = new float[n];
    float[] cmax = new float[n];
    for (int i=0; i<n; ++i)
      cmax[i] = -FLT_MAX;

    // Search begins in the middle of the specified range of lags.
    Lags lags = new Lags(min,max);
    int l = (min+max)/2;

    // While lags remain to be processed, ...
    boolean done = false;
    while (!done) {

      // Apply correlation filter for this lag.
      correlate(l,c);
      if (_type==Type.SIMPLE)
        normalize(l,c);

      // Correlations have been computed for this lag, 
      // but no maxima have yet been found.
      lags.markLag(l);

      // Look for maxima; if found, mark this lag accordingly.
      boolean foundMax = false;
      for (int i=0; i<n; ++i) { 
        float ci = c[i];
        if (ci>cmax[i]) {
          cmax[i] = ci;
          lag[i] = (byte)l;
          foundMax = true;
        }
      }
      if (foundMax)
        lags.markMax(l);

      // Which lag to process next?
      int[] ls = lags.nextLag();
      if (ls==null) {
        done = true;
      } else {
        l = ls[0];
      }
    }
  }

  /**
   * Searches for lags for which cross-correlations are maximized.
   * @param min1 minimum lag in 1st dimension
   * @param max1 maximum lag in 1st dimension
   * @param min2 minimum lag in 2nd dimension
   * @param max2 maximum lag in 2nd dimension
   * @param lag1 output array of lags in the 1st dimension.
   * @param lag2 output array of lags in the 2nd dimension.
   */
  public void findMaxLags(
    int min1, int max1, int min2, int max2,
    byte[][] lag1, byte[][] lag2) 
  {
    checkDimension(2);

    // Initialize arrays of lags.
    int n1 = _n1;
    int n2 = _n2;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        lag1[i2][i1] = 0;
        lag2[i2][i1] = 0;
      }
    }

    // Array for cross-correlations.
    float[][] c = new float[n2][n1];
    float[][] cmax = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        cmax[i2][i1] = -FLT_MAX;

    // Search begins in the middle of the specified range of lags.
    Lags lags = new Lags(min1,max1,min2,max2);
    int l1 = (min1+max1)/2;
    int l2 = (min2+max2)/2;

    // While lags remain to be processed, ...
    boolean done = false;
    while (!done) {

      // Apply correlation filter for this lag.
      correlate(l1,l2,c);
      if (_type==Type.SIMPLE)
        normalize(l1,l2,c);

      // Correlations have been computed for this lag, 
      // but no maxima have yet been found.
      lags.markLag(l1,l2);

      // Look for maxima; if found, mark this lag accordingly.
      boolean foundMax = false;
      for (int i2=0; i2<n2; ++i2) { 
        float[] c2 = c[i2];
        float[] cmax2 = cmax[i2];
        for (int i1=0; i1<n1; ++i1) { 
          float ci = c2[i1];
          if (ci>cmax2[i1]) {
            cmax[i2][i1] = ci;
            lag1[i2][i1] = (byte)l1;
            lag2[i2][i1] = (byte)l2;
            foundMax = true;
          }
        }
      }
      if (foundMax)
        lags.markMax(l1,l2);

      // Which lag to process next?
      int[] ls = lags.nextLag();
      if (ls==null) {
        done = true;
      } else {
        l1 = ls[0];
        l2 = ls[1];
      }
    }
  }

  public void findMaxLags(
    int min1, int max1, int min2, int max2, int min3, int max3,
    byte[][][] lag1, byte[][][] lag2, byte[][][] lag3) 
  {
    checkDimension(3);

    // Initialize arrays of lags.
    int n1 = _n1;
    int n2 = _n2;
    int n3 = _n3;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          lag1[i3][i2][i1] = 0;
          lag2[i3][i2][i1] = 0;
          lag3[i3][i2][i1] = 0;
        }
      }
    }

    // Array for cross-correlations.
    float[][][] c = new float[n3][n2][n1];
    float[][][] cmax = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          cmax[i3][i2][i1] = -FLT_MAX;

    // Search begins in the middle of the specified range of lags.
    Lags lags = new Lags(min1,max1,min2,max2,min3,max3);
    int l1 = (min1+max1)/2;
    int l2 = (min2+max2)/2;
    int l3 = (min3+max3)/2;

    // While lags remain to be processed, ...
    boolean done = false;
    while (!done) {

      // Apply correlation filter for this lag.
      System.out.println("findMaxLags: l1="+l1+" l2="+l2+" l3="+l3);
      correlate(l1,l2,l3,c);
      if (_type==Type.SIMPLE)
        normalize(l1,l2,l3,c);

      // Correlations have been computed for this lag, 
      // but no maxima have yet been found.
      lags.markLag(l1,l2,l3);

      // Look for maxima; if found, mark this lag accordingly.
      boolean foundMax = false;
      for (int i3=0; i3<n3; ++i3) { 
        for (int i2=0; i2<n2; ++i2) { 
          float[] c32 = c[i3][i2];
          float[] cmax32 = cmax[i3][i2];
          byte[] lag132 = lag1[i3][i2];
          byte[] lag232 = lag2[i3][i2];
          byte[] lag332 = lag3[i3][i2];
          for (int i1=0; i1<n1; ++i1) { 
            float ci = c32[i1];
            if (ci>cmax32[i1]) {
              cmax32[i1] = ci;
              lag132[i1] = (byte)l1;
              lag232[i1] = (byte)l2;
              lag332[i1] = (byte)l3;
              foundMax = true;
            }
          }
        }
      }
      if (foundMax)
        lags.markMax(l1,l2,l3);

      // Which lag to process next?
      int[] ls = lags.nextLag();
      if (ls==null) {
        done = true;
      } else {
        l1 = ls[0];
        l2 = ls[1];
        l3 = ls[2];
      }
    }
  }

  public void refineLags(byte[] l, float[] u) {
    int n = _n1;

    // Minimum and maximum lags.
    int min = l[0];
    int max = l[0];
    for (int i=1; i<n; ++i) {
      int lag = l[i];
      if (lag<min) min = lag;
      if (lag>max) max = lag;
    }
    System.out.println("refineLags:");
    System.out.println("  min="+min+" max="+max);

    // Coefficients for quadratic fit.
    float[] c = new float[n];
    float[] a1 = new float[n];
    float[] a2 = new float[n];
    for (int lag=min-1; lag<=max+1; ++lag) {
      correlate(lag,c);
      if (_type==Type.SIMPLE)
        normalize(lag,c);
      for (int i=0; i<n; ++i) {
        int k = lag-l[i];
        if (-1<=k && k<=1) {
          k += 1;
          float[] ck = C1[k];
          float ci = c[i];
          a1[i] += ck[1]*ci;
          a2[i] += ck[2]*ci;
        }
      }
    }

    // Refined lags.
    for (int i=0; i<n; ++i) {
      float a1i = a1[i];
      float a2i = a2[i];
      float w = 0.0f;
      if (a2i<0.0)
        w = -0.5f*a1i/a2i;
      if (w<-1.0f) {
        w = -1.0f;
      } else if (w>1.0f) {
        w = 1.0f;
      }
      u[i] = (float)(w+l[i]);
    }
  }

  public void refineLags(
    byte[][] l1, byte[][] l2,
    float[][] u1, float[][] u2) 
  {
    refineLags(l1,l2,u1,u2,null);
  }

  public void refineLags(
    byte[][] l1, byte[][] l2,
    float[][] u1, float[][] u2,
    float[][][] q)
  {
    int n1 = _n1;
    int n2 = _n2;

    // Minimum and maximum lags.
    int min1 = l1[0][0];
    int max1 = l1[0][0];
    int min2 = l2[0][0];
    int max2 = l2[0][0];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        int lag1 = l1[i2][i1];
        int lag2 = l2[i2][i1];
        if (lag1<min1) min1 = lag1;
        if (lag1>max1) max1 = lag1;
        if (lag2<min2) min2 = lag2;
        if (lag2>max2) max2 = lag2;
      }
    }
    System.out.println("refineLags:");
    System.out.println("  min1="+min1+" max1="+max1);
    System.out.println("  min2="+min2+" max2="+max2);

    // Coefficients for quadratic fit.
    int nbad = 0;
    float[][][][] ca = new float[n2][n1][3][3];
    float[][] c = new float[n2][n1];
    float[][] a0 = new float[n2][n1];
    float[][] a1 = new float[n2][n1];
    float[][] a2 = new float[n2][n1];
    float[][] a3 = new float[n2][n1];
    float[][] a4 = new float[n2][n1];
    float[][] a5 = new float[n2][n1];
    for (int lag2=min2-1; lag2<=max2+1; ++lag2) {
      for (int lag1=min1-1; lag1<=max1+1; ++lag1) {
        correlate(lag1,lag2,c);
        if (_type==Type.SIMPLE || q!=null)
          normalize(lag1,lag2,c);
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            int k1 = lag1-l1[i2][i1];
            int k2 = lag2-l2[i2][i1];
            if (-1<=k1 && k1<=1 && -1<=k2 && k2<=1) {
              int k = (k1+1)+3*(k2+1);
              float[] ck = C2[k];
              float ci = c[i2][i1];
              ca[i2][i1][k2+1][k1+1] = ci;
              a0[i2][i1] += ck[0]*ci;
              a1[i2][i1] += ck[1]*ci;
              a2[i2][i1] += ck[2]*ci;
              a3[i2][i1] += ck[3]*ci;
              a4[i2][i1] += ck[4]*ci;
              a5[i2][i1] += ck[5]*ci;
            }
          }
        }
      }
    }

    // Cholesky decomposition solves 2x2 system for refined lags.
    int i1min = -1;
    int i2min = -1;
    float d0min = FLT_MAX;
    int i1max = -1;
    int i2max = -1;
    float w1max = 0.0f;
    float w2max = 0.0f;
    float[][] a = new float[2][2];
    float[][] v = new float[2][2];
    float[] d = new float[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        double aa0 = a0[i2][i1];
        double aa1 = a1[i2][i1];
        double aa2 = a2[i2][i1];
        double aa3 = a3[i2][i1];
        double aa4 = a4[i2][i1];
        double aa5 = a5[i2][i1];
        double w1 = 0.0;
        double w2 = 0.0;
        double b1 = aa1;
        double b2 = aa2;
        double a21 = -aa3;
        double a11 = -2.0*aa4;
        double a22 = -2.0*aa5;
        a[0][0] = (float)a11;  a[0][1] = (float)a21;
        a[1][0] = (float)a21;  a[1][1] = (float)a22;
        Eigen.solveSymmetric22(a,v,d);
        float e = 0.001f*max(d[0],d[1]);
        if (d[0]<d0min) {
          d0min = d[0];
          i1min = i1;
          i2min = i2;
        }
        d[0] = max(e,d[0]);
        d[1] = max(e,d[1]);
        a11 = d[0]*v[0][0]*v[0][0]+d[1]*v[1][0]*v[1][0];
        a21 = d[0]*v[0][0]*v[0][1]+d[1]*v[1][0]*v[1][1];
        a22 = d[0]*v[0][1]*v[0][1]+d[1]*v[1][1]*v[1][1];
        boolean pd = false;
        double d11 = a11;
        if (d11>0.0) {
          double l11 = sqrt(d11);
          double l21 = a21/l11;
          double d22 = a22-l21*l21;
          if (d22>0.0) {
            double l22 = sqrt(d22);
            double v1 = b1/l11;
            double v2 = (b2-l21*v1)/l22;
            w2 = v2/l22;
            w1 = (v1-l21*w2)/l11;
            if (w1>w1max || w2>w2max) {
              w1max = (float)w1;
              w2max = (float)w2;
              i1max = i1;
              i2max = i2;
            }
            if (w1<-0.5) {
              w1 = -0.5;
            } else if (w1>0.5) {
              w1 = 0.5;
            }
            if (w2<-0.5) {
              w2 = -0.5;
            } else if (w2>0.5) {
              w2 = 0.5;
            }
            pd = true;
          }
        }
        if (!pd)
          System.out.println("!pd i1="+i1+" i2="+i2);
        if (abs(w1)==1.0 || abs(w2)==1.0) {
          //System.out.println("i1="+i1+" i2="+i2+" w1="+w1+" w2="+w2);
          //Array.dump(ca[i2][i1]);
          ++nbad;
        }

        // Refined lags.
        u1[i2][i1] = (float)(w1+l1[i2][i1]);
        u2[i2][i1] = (float)(w2+l2[i2][i1]);

        // Optional fitting coefficients.
        if (q!=null) {
          q[0][i2][i1] = 0.0f;
          q[1][i2][i1] = 0.0f;
          q[2][i2][i1] = 0.0f;
          q[3][i2][i1] = 0.0f;
          double cp = aa0+aa1*w1+aa2*w2+aa3*w1*w2+aa4*w1*w1+aa5*w2*w2;
          if (cp>0.0) {
            q[0][i2][i1] = (float)cp;
            q[1][i2][i1] = (float)(0.5*a11);
            q[2][i2][i1] = (float)(0.5*a21);
            q[3][i2][i1] = (float)(0.5*a22);
          }
        }
      }
    }
    System.out.println("refineLags: nbad="+nbad);
    System.out.println(
      "refineLags: d0min="+d0min + " i1min="+i1min + " i2min="+i2min);
    Array.dump(ca[i2min][i1min]);
    System.out.println("lags: l1="+l1[i2min][i1min]+" l2="+l2[i2min][i1min]);
    System.out.println("refineLags: w1max="+w1max + " w2max="+w2max);
    System.out.println("refineLags: i1max="+i1max + " i2max="+i2max);
    Array.dump(ca[i2max][i1max]);
    System.out.println("lags: l1="+l1[i2max][i1max]+" l2="+l2[i2max][i1max]);
  }

  public void refineLags(
    byte[][][] l1, byte[][][] l2, byte[][][] l3,
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n1 = _n1;
    int n2 = _n2;
    int n3 = _n3;

    // Minimum and maximum lags.
    int min1 = l1[0][0][0];
    int max1 = l1[0][0][0];
    int min2 = l2[0][0][0];
    int max2 = l2[0][0][0];
    int min3 = l3[0][0][0];
    int max3 = l3[0][0][0];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          int lag1 = l1[i3][i2][i1];
          int lag2 = l2[i3][i2][i1];
          int lag3 = l3[i3][i2][i1];
          if (lag1<min1) min1 = lag1;
          if (lag1>max1) max1 = lag1;
          if (lag2<min2) min2 = lag2;
          if (lag2>max2) max2 = lag2;
          if (lag3<min3) min3 = lag3;
          if (lag3>max3) max3 = lag3;
        }
      }
    }
    System.out.println("refineLags:");
    System.out.println("  min1="+min1+" max1="+max1);
    System.out.println("  min2="+min2+" max2="+max2);
    System.out.println("  min3="+min3+" max3="+max3);

    // Coefficients for quadratic fit.
    float[][][] c = new float[n3][n2][n1];
    float[][][] a1 = new float[n3][n2][n1];
    float[][][] a2 = new float[n3][n2][n1];
    float[][][] a3 = new float[n3][n2][n1];
    float[][][] a4 = new float[n3][n2][n1];
    float[][][] a5 = new float[n3][n2][n1];
    float[][][] a6 = new float[n3][n2][n1];
    float[][][] a7 = new float[n3][n2][n1];
    float[][][] a8 = new float[n3][n2][n1];
    float[][][] a9 = new float[n3][n2][n1];
    for (int lag3=min3-1; lag3<=max3+1; ++lag3) {
      for (int lag2=min2-1; lag2<=max2+1; ++lag2) {
        for (int lag1=min1-1; lag1<=max1+1; ++lag1) {
          System.out.println("("+lag1+","+lag2+","+lag3+")");
          correlate(lag1,lag2,lag3,c);
          if (_type==Type.SIMPLE)
            normalize(lag1,lag2,lag3,c);
          for (int i3=0; i3<n3; ++i3) {
            for (int i2=0; i2<n2; ++i2) {
              for (int i1=0; i1<n1; ++i1) {
                int k1 = lag1-l1[i3][i2][i1];
                int k2 = lag2-l2[i3][i2][i1];
                int k3 = lag3-l3[i3][i2][i1];
                if (-1<=k1 && k1<=1 && -1<=k2 && k2<=1 && -1<=k3 &&  k3<=1) {
                  int k = (k1+1)+3*(k2+1)+9*(k3+1);
                  float[] ck = C3[k];
                  float ci = c[i3][i2][i1];
                  a1[i3][i2][i1] += ck[1]*ci;
                  a2[i3][i2][i1] += ck[2]*ci;
                  a3[i3][i2][i1] += ck[3]*ci;
                  a4[i3][i2][i1] += ck[4]*ci;
                  a5[i3][i2][i1] += ck[5]*ci;
                  a6[i3][i2][i1] += ck[6]*ci;
                  a7[i3][i2][i1] += ck[7]*ci;
                  a8[i3][i2][i1] += ck[8]*ci;
                  a9[i3][i2][i1] += ck[9]*ci;
                }
              }
            }
          }
        }
      }
    }

    // Cholesky decomposition solves 3x3 system for refined lags.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          double w1 = 0.0;
          double w2 = 0.0;
          double w3 = 0.0;
          double b1 = a1[i3][i2][i1];
          double b2 = a2[i3][i2][i1];
          double b3 = a3[i3][i2][i1];
          double a21 = -a4[i3][i2][i1];
          double a31 = -a5[i3][i2][i1];
          double a32 = -a6[i3][i2][i1];
          double a11 = -2.0*a7[i3][i2][i1];
          double a22 = -2.0*a8[i3][i2][i1];
          double a33 = -2.0*a9[i3][i2][i1];
          double d11 = a11;
          if (d11>0.0) {
            double l11 = sqrt(d11);
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
                w3 = v3/l33;
                w2 = (v2-l32*w3)/l22;
                w1 = (v1-l21*w2-l31*w3)/l11;
                if (w1<-1.0) {
                  w1 = -1.0;
                } else if (w1>1.0) {
                  w1 = 1.0;
                }
                if (w2<-1.0) {
                  w2 = -1.0;
                } else if (w2>1.0) {
                  w2 = 1.0;
                }
                if (w3<-1.0) {
                  w3 = -1.0;
                } else if (w3>1.0) {
                  w3 = 1.0;
                }
              }
            }
          }
          u1[i3][i2][i1] = (float)(w1+l1[i3][i2][i1]);
          u2[i3][i2][i1] = (float)(w2+l2[i3][i2][i1]);
          u3[i3][i2][i1] = (float)(w3+l3[i3][i2][i1]);
        }
      }
    }
  }

  public void applyWindow(
    int jf, float[] f, 
    int jg, float[] g) 
  {
    float[] w1 = makeGaussianWindow(_sigma1);
    applyWindow(w1,jf,f,jg,g);
  }

  public void applyWindow(
    int j1f, int j2f, float[][] f, 
    int j1g, int j2g, float[][] g) 
  {
    float[] w1 = makeGaussianWindow(_sigma1);
    float[] w2 = makeGaussianWindow(_sigma2);
    applyWindow(w1,w2,j1f,j2f,f,j1g,j2g,g);
  }

  public void applyWindow(
    int j1f, int j2f, int j3f, float[][][] f, 
    int j1g, int j2g, int j3g, float[][][] g) 
  {
    float[] w1 = makeGaussianWindow(_sigma1);
    float[] w2 = makeGaussianWindow(_sigma2);
    float[] w3 = makeGaussianWindow(_sigma3);
    applyWindow(w1,w2,w3,j1f,j2f,j3f,f,j1g,j2g,j3g,g);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  // Information about lags used when searching for correlation maxima.
  // Lags for which correlations have been computed are marked. Lags for 
  // which correlation maxima have been found are marked differently.
  // The next lag to process is one that is not marked, but that is 
  // adjacent to a lag for which a maximum has been found. If no such 
  // next lag exists, then the next lag is null.
  private static class Lags {
    Lags(int min1, int max1) {
      this(min1,max1,0,0,0,0);
    }
    Lags(int min1, int max1, int min2, int max2) {
      this(min1,max1,min2,max2,0,0);
    }
    Lags(int min1, int max1, int min2, int max2, int min3, int max3) {
      int nl1 = 1+max1-min1;
      int nl2 = 1+max2-min2;
      int nl3 = 1+max3-min3;
      _min1 = min1;
      _max1 = max1;
      _min2 = min2;
      _max2 = max2;
      _min3 = min3;
      _max3 = max3;
      _mark = new byte[nl3][nl2][nl1];
    }
    void markLag(int l1) {
      markLag(l1,0,0);
    }
    void markLag(int l1, int l2) {
      markLag(l1,l2,0);
    }
    void markLag(int l1, int l2, int l3) {
      _mark[l3-_min3][l2-_min2][l1-_min1] = -1;
    }
    void markMax(int l1) {
      markMax(l1,0,0);
    }
    void markMax(int l1, int l2) {
      markMax(l1,l2,0);
    }
    void markMax(int l1, int l2, int l3) {
      _mark[l3-_min3][l2-_min2][l1-_min1] = 1;
    }
    boolean isMarkedLag(int l1) {
      return isMarkedLag(l1,0,0);
    }
    boolean isMarkedLag(int l1, int l2) {
      return isMarkedLag(l1,l2,0);
    }
    boolean isMarkedLag(int l1, int l2, int l3) {
      return !inBounds(l1,l2,l3) || _mark[l3-_min3][l2-_min2][l1-_min1]!=0;
    }
    boolean isMarkedMax(int l1) {
      return isMarkedMax(l1,0,0);
    }
    boolean isMarkedMax(int l1, int l2) {
      return isMarkedMax(l1,l2,0);
    }
    boolean isMarkedMax(int l1, int l2, int l3) {
      return inBounds(l1,l2,l3) && _mark[l3-_min3][l2-_min2][l1-_min1]==1;
    }
    boolean inBounds(int l1, int l2, int l3) {
      return l1>=_min1 && l1<=_max1 &&
             l2>=_min2 && l2<=_max2 &&
             l3>=_min3 && l3<=_max3;
    }
    int[] nextLag() {
      for (int l3=_min3; l3<=_max3; ++l3) {
        for (int l2=_min2; l2<=_max2; ++l2) {
          for (int l1=_min1; l1<=_max1; ++l1) {
            if (isMarkedMax(l1,l2,l3)) {
              for (int k3=l3-1; k3<=l3+1; ++k3) {
                for (int k2=l2-1; k2<=l2+1; ++k2) {
                  for (int k1=l1-1; k1<=l1+1; ++k1) {
                    if (!isMarkedLag(k1,k2,k3)) {
                      return new int[]{k1,k2,k3};
                    }
                  }
                }
              }
            }
          }
        }
      }
      return null;
    }
    int _min1,_max1,_min2,_max2,_min3,_max3;
    byte[][][] _mark;
  }

  // Fractions used in tables of coefficients below.
  private static final float C00 = 0.0f;
  private static final float C11 = 1.0f/1.0f;
  private static final float C12 = 1.0f/2.0f;
  private static final float C13 = 1.0f/3.0f;
  private static final float C14 = 1.0f/4.0f;
  private static final float C16 = 1.0f/6.0f;
  private static final float C19 = 1.0f/9.0f;
  private static final float C29 = 2.0f/9.0f;
  private static final float C59 = 5.0f/9.0f;
  private static final float C000 = 0.0f;
  private static final float C109 = 1.0f/ 9.0f;
  private static final float C112 = 1.0f/12.0f;
  private static final float C118 = 1.0f/18.0f;
  private static final float C127 = 1.0f/27.0f;
  private static final float C227 = 2.0f/27.0f;
  private static final float C427 = 4.0f/27.0f;
  private static final float C727 = 7.0f/27.0f;

  // Coefficients for 1-D lag refinement. Assume a sampled correlation
  // maximum at integer lag l. Let u be the fractional component of a
  // refined lag l+u. Near its maximum, we approximate the correlation
  // function c(l+u) by a quadratic function fit that interpolates the
  // three sampled correlation values surrounding the sample c(l).
  // The quadratic function is c(l+u) = a0 + a1*u + a2*u*u.
  // The three sampled correlation values c(l-1), c(l), and c(l+1)
  // yield three equations for the three coefficients a0, a1, and a2.
  // The rows of this array contain the weights that we apply to each
  // of the three correlation values in the computation of the three
  // quadratic coefficients. For example, the first (top) row contains
  // the weights applied to the sampled correlation value c(l-1).
  // When refining correlation lags, we use this table of weights to
  // accumulate the contributions of the three correlation values nearest
  // to each sampled correlation maximum. For lag refinement, we need only
  // the two coefficients a1 and a2, for the peak of the quadratic is at
  // u = -0.5*a1/a2.
  private static final float[][] C1 = {
  //  a0    a1    a2
    { C00, -C12,  C12}, // -1
    { C11,  C00, -C11}, //  0
    { C00,  C12,  C12}, //  1
  };


  // Coefficients for 2-D lag refinement. Assume a sampled correlation 
  // maximum at integer lag (l1,l2). Let (u1,u2) be fractional components 
  // of a refined lag (l1+u1,l2+u2). Near its maximum, we approximate the 
  // correlation function c(l1+u1,l2+u2) by a quadratic function least-
  // squares fit to the nine sampled correlation values surrounding the 
  // sampled c(l1,l2).
  // The quadratic function is
  // c(l1+u1,l2+u2) = a0 + a1*u1 + a2*u2 + a3*u1*u2 + a4*u1*u1 + a5*u2*u2
  // The nine sampled correlation values yield nine equations for the six
  // coefficients a0, a1, a2, a3, a4, and a5.
  // By QR decomposition of this overdetermined system of equations, we 
  // obtained the following array of constants. The rows of this array 
  // contain the weights that we apply to each of the nine correlation 
  // values in the computation of the six quadratic coefficients. For
  // example, the first (top) row contains the weights applied to the 
  // sampled correlation value c(l1-1,l2-1).
  // When refining correlation lags, we use this table of weights to 
  // accumulate the contributions of the nine correlation values nearest
  // to each sampled correlation maximum. For lag refinement, we need only 
  // the five coefficients a1, a2, a3, a4, and a5.
  private static final float[][] C2_LEAST_SQUARES = {
  //  a0    a1    a2    a3    a4    a5
    {-C19, -C16, -C16,  C14,  C16,  C16}, // (-1,-1)
    { C29,  C00, -C16,  C00, -C13,  C16}, // ( 0,-1)
    {-C19,  C16, -C16, -C14,  C16,  C16}, // ( 1,-1)
    { C29, -C16,  C00,  C00,  C16, -C13}, // (-1, 0)
    { C59,  C00,  C00,  C00, -C13, -C13}, // ( 0, 0)
    { C29,  C16,  C00,  C00,  C16, -C13}, // ( 1, 0)
    {-C19, -C16,  C16, -C14,  C16,  C16}, // (-1, 1)
    { C29,  C00,  C16,  C00, -C13,  C16}, // ( 0, 1)
    {-C19,  C16,  C16,  C14,  C16,  C16}, // ( 1, 1)
  };
  private static final float[][] C2_FINITE_DIFFERENCE = {
    { C00,  C00,  C00,  C14,  C00,  C00}, // (-1,-1)
    { C00,  C00, -C12,  C00,  C00,  C12}, // ( 0,-1)
    { C00,  C00,  C00, -C14,  C00,  C00}, // ( 1,-1)
    { C00, -C12,  C00,  C00,  C12,  C00}, // (-1, 0)
    { C11,  C00,  C00,  C00, -C11, -C11}, // ( 0, 0)
    { C00,  C12,  C00,  C00,  C12,  C00}, // ( 1, 0)
    { C00,  C00,  C00, -C14,  C00,  C00}, // (-1, 1)
    { C00,  C00,  C12,  C00,  C00,  C12}, // ( 0, 1)
    { C00,  C00,  C00,  C14,  C00,  C00}, // ( 1, 1)
  };
  private static final float[][] C2_FD_WITHOUT_CROSS_DERIVATIVES = {
    { C00,  C00,  C00,  C00,  C00,  C00}, // (-1,-1)
    { C00,  C00, -C12,  C00,  C00,  C12}, // ( 0,-1)
    { C00,  C00,  C00,  C00,  C00,  C00}, // ( 1,-1)
    { C00, -C12,  C00,  C00,  C12,  C00}, // (-1, 0)
    { C11,  C00,  C00,  C00, -C11, -C11}, // ( 0, 0)
    { C00,  C12,  C00,  C00,  C12,  C00}, // ( 1, 0)
    { C00,  C00,  C00,  C00,  C00,  C00}, // (-1, 1)
    { C00,  C00,  C12,  C00,  C00,  C12}, // ( 0, 1)
    { C00,  C00,  C00,  C00,  C00,  C00}, // ( 1, 1)
  };
  private static final float[][] C2 = C2_LEAST_SQUARES;
  //private static final float[][] C2 = C2_FINITE_DIFFERENCE;
  //private static final float[][] C2 = C2_FD_WITHOUT_CROSS_DERIVATIVES;

  // Coefficients for 3-D lag refinement. Here we fit 27 sampled correlation 
  // values with 10 coefficients of a 3-D quadratic function
  // c(l1+u1,l2+u2,l3+u3) = a0 + a1*u1    + a2*u2    + a3*u3    +
  //                             a4*u1*u2 + a5*u1*u3 + a6*u2*u3 +
  //                             a7*u1*u1 + a8*u2*u2 + a9*u3*u3
  // For lag refinement, we need only the last nine coefficients.
  private static final float[][] C3 = {
  //   a0     a1     a2     a3     a4     a5     a6     a7     a8     a9 
    {-C227, -C118, -C118, -C118,  C112,  C112,  C112,  C118,  C118,  C118},
    { C127,  C000, -C118, -C118,  C000,  C000,  C112, -C109,  C118,  C118}, 
    {-C227,  C118, -C118, -C118, -C112, -C112,  C112,  C118,  C118,  C118}, 
    { C127, -C118,  C000, -C118,  C000,  C112,  C000,  C118, -C109,  C118},
    { C427,  C000,  C000, -C118,  C000,  C000,  C000, -C109, -C109,  C118},
    { C127,  C118,  C000, -C118,  C000, -C112,  C000,  C118, -C109,  C118},
    {-C227, -C118,  C118, -C118, -C112,  C112, -C112,  C118,  C118,  C118}, 
    { C127,  C000,  C118, -C118,  C000,  C000, -C112, -C109,  C118,  C118}, 
    {-C227,  C118,  C118, -C118,  C112, -C112, -C112,  C118,  C118,  C118}, 
    { C127, -C118, -C118,  C000,  C112,  C000,  C000,  C118,  C118, -C109}, 
    { C427,  C000, -C118,  C000,  C000,  C000,  C000, -C109,  C118, -C109}, 
    { C127,  C118, -C118,  C000, -C112,  C000,  C000,  C118,  C118, -C109}, 
    { C427, -C118,  C000,  C000,  C000,  C000,  C000,  C118, -C109, -C109}, 
    { C727,  C000,  C000,  C000,  C000,  C000,  C000, -C109, -C109, -C109}, 
    { C427,  C118,  C000,  C000,  C000,  C000,  C000,  C118, -C109, -C109}, 
    { C127, -C118,  C118,  C000, -C112,  C000,  C000,  C118,  C118, -C109}, 
    { C427,  C000,  C118,  C000,  C000,  C000,  C000, -C109,  C118, -C109}, 
    { C127,  C118,  C118,  C000,  C112,  C000,  C000,  C118,  C118, -C109}, 
    {-C227, -C118, -C118,  C118,  C112, -C112, -C112,  C118,  C118,  C118}, 
    { C127,  C000, -C118,  C118,  C000,  C000, -C112, -C109,  C118,  C118}, 
    {-C227,  C118, -C118,  C118, -C112,  C112, -C112,  C118,  C118,  C118}, 
    { C127, -C118,  C000,  C118,  C000, -C112,  C000,  C118, -C109,  C118}, 
    { C427,  C000,  C000,  C118,  C000,  C000,  C000, -C109, -C109,  C118}, 
    { C127,  C118,  C000,  C118,  C000,  C112,  C000,  C118, -C109,  C118}, 
    {-C227, -C118,  C118,  C118, -C112, -C112,  C112,  C118,  C118,  C118}, 
    { C127,  C000,  C118,  C118,  C000,  C000,  C112, -C109,  C118,  C118}, 
    {-C227,  C118,  C118,  C118,  C112,  C112,  C112,  C118,  C118,  C118},
  };

  private Window _window = Window.GAUSSIAN; // window for correlations
  private Type _type = Type.SYMMETRIC; // correlation type
  private double _sigma1,_sigma2,_sigma3; // window half-widths
  private Filter _f1,_f2,_f3; // filters used to implement windows
  private int _dimension; // dimension of input arrays; 0 if no inputs
  private int _n1,_n2,_n3; // array lengths
  private float[][][] _f,_g; // inputs f and g; by reference
  private float[][][][] _s; // normalization scale factors

  // Kaiser-windowed sinc interpolation coefficients for half-sample shifts.
  private static float S1 =  0.6157280f;
  private static float S2 = -0.1558022f;
  private static float S3 =  0.0509014f;
  private static float S4 = -0.0115417f;
  private static float[] S = {S4,S3,S2,S1,S1,S2,S3,S4};

  private void correlate(int lag, float[] f, float[] g, float[] c) {
    Check.argument(f!=c,"f!=c");
    Check.argument(g!=c,"g!=c");
    int n1 = f.length;
    int l1 = lag;

    // Conventional lags for f and g for simple correlation.
    int l1f = 0;
    int l1g = l1;

    // Shifted lags for symmetric correlation.
    if (_type==Type.SYMMETRIC) {
      // Examples of symmetric lags:
      // lag  ...  -2  -1   0   1   2  ...
      // l1f  ...  -1  -1   0   0   1  ...
      // l1g  ...  -1   0   0   1   1  ...
      l1f = (l1>=0)?(l1+0)/2:(l1-1)/2;
      l1g = (l1>=0)?(l1+1)/2:(l1+0)/2;
    }

    // Scale factor so that center of window = 1.
    double scale1 = 1.0;
    if (_window==Window.GAUSSIAN) {
      scale1 *= sqrt(2.0*PI)*_sigma1;
    } else {
      scale1 *= 1.0+2.0*_sigma1;
    }

    // If symmetric correlation, need extra lag-dependent scaling.
    // This scaling accounts for the separation (by lag samples) of 
    // the two windows implicitly applied to f and g. The filter we
    // apply below to the correlation product h is the product of 
    // those two windows.
    if (_type==Type.SYMMETRIC) {
      if (_window==Window.GAUSSIAN) {
        scale1 *= exp((-0.125*l1*l1)/(_sigma1*_sigma1));
      } else {
        scale1 *= max(0.0,1.0+2.0*_sigma1-abs(l1))/(1.0+2.0*_sigma1);
      }
    }
    float scale = (float)scale1;

    // Correlation product.
    float[] h = new float[n1];
    int i1min = max(0,l1f,-l1g);
    int i1max = min(n1,n1+l1f,n1-l1g);
    for (int i1=i1min; i1<i1max; ++i1) {
      h[i1] = scale*f[i1-l1f]*g[i1+l1g];
    }

    // If Gaussian and symmetric and odd lag, delay (shift) by 1/2 sample.
    if (_window==Window.GAUSSIAN && _type==Type.SYMMETRIC) {
      if (l1f!=l1g) {
        shift(h,c);
        Array.copy(c,h);
      }
    }

    // Filter correlation product with window. For symmetric correlations
    // with a rectangle window, the width of the product rectangle depends 
    // on the lag, so we construct a new rectangle filter for each lag.
    Filter f1 = _f1;
    if (_window==Window.RECTANGLE && _type==Type.SYMMETRIC)
      f1 = new RectangleFilter(_sigma1,l1);
    f1.apply(h,c);
  }

  private void correlate(
    int lag1, int lag2, float[][] f, float[][] g, float[][] c) 
  {
    Check.argument(f!=c,"f!=c");
    Check.argument(g!=c,"g!=c");
    int n1 = f[0].length;
    int n2 = f.length;
    int l1 = lag1;
    int l2 = lag2;

    // Conventional lags for f and g for simple correlation.
    int l1f = 0;
    int l1g = l1;
    int l2f = 0;
    int l2g = l2;

    // Shifted lags for symmetric correlation.
    if (_type==Type.SYMMETRIC) {
      l1f = (l1>=0)?(l1+0)/2:(l1-1)/2;
      l1g = (l1>=0)?(l1+1)/2:(l1+0)/2;
      l2f = (l2>=0)?(l2+0)/2:(l2-1)/2;
      l2g = (l2>=0)?(l2+1)/2:(l2+0)/2;
    }

    // Scale factor so that center of window = 1.
    double scale1 = 1.0;
    double scale2 = 1.0;
    if (_window==Window.GAUSSIAN) {
      scale1 *= sqrt(2.0*PI)*_sigma1;
      scale2 *= sqrt(2.0*PI)*_sigma2;
    } else {
      scale1 *= 1.0+2.0*_sigma1;
      scale2 *= 1.0+2.0*_sigma2;
    }

    // If symmetric correlation, need extra lag-dependent scaling.
    // This scaling accounts for the separation (by lag samples) of 
    // the two windows implicitly applied to f and g. The filter we
    // apply below to the correlation product h is the product of 
    // those two windows.
    if (_type==Type.SYMMETRIC) {
      if (_window==Window.GAUSSIAN) {
        scale1 *= exp((-0.125*l1*l1)/(_sigma1*_sigma1));
        scale2 *= exp((-0.125*l2*l2)/(_sigma2*_sigma2));
      } else {
        scale1 *= max(0.0,1.0+2.0*_sigma1-abs(l1))/(1.0+2.0*_sigma1);
        scale2 *= max(0.0,1.0+2.0*_sigma2-abs(l2))/(1.0+2.0*_sigma2);
      }
    }
    float scale = (float)(scale1*scale2);

    // Correlation product.
    float[][] h = new float[n2][n1];
    int i1min = max(0,l1f,-l1g);
    int i1max = min(n1,n1+l1f,n1-l1g);
    int i2min = max(0,l2f,-l2g);
    int i2max = min(n2,n2+l2f,n2-l2g);
    for (int i2=i2min; i2<i2max; ++i2) {
      float[] f2 = f[i2-l2f];
      float[] g2 = g[i2+l2g];
      float[] h2 = h[i2];
      for (int i1=i1min; i1<i1max; ++i1) {
        h2[i1] = scale*f2[i1-l1f]*g2[i1+l1g];
      }
    }

    // If Gaussian and symmetric and odd lag, delay (shift) by 1/2 sample.
    if (_window==Window.GAUSSIAN && _type==Type.SYMMETRIC) {
      if (l1f!=l1g) {
        shift1(h,c);
        Array.copy(c,h);
      }
      if (l2f!=l2g) {
        shift2(h,c);
        Array.copy(c,h);
      }
    }

    // Filter correlation product with window. For symmetric correlations
    // with a rectangle window, the width of the product rectangle depends 
    // on the lag, so we construct a new rectangle filter for each lag.
    Filter f1 = _f1;
    Filter f2 = _f2;
    if (_window==Window.RECTANGLE && _type==Type.SYMMETRIC) {
      f1 = new RectangleFilter(_sigma1,l1);
      f2 = new RectangleFilter(_sigma2,l2);
    }
    f1.apply1(h,c);
    Array.copy(c,h);
    f2.apply2(h,c);
  }

  private void correlate(
    int lag1, int lag2, int lag3, float[][][] f, float[][][] g, float[][][] c) 
  {
    Check.argument(f!=c,"f!=c");
    Check.argument(g!=c,"g!=c");
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    int l1 = lag1;
    int l2 = lag2;
    int l3 = lag3;

    // Conventional lags for f and g for simple correlation.
    int l1f = 0;
    int l1g = l1;
    int l2f = 0;
    int l2g = l2;
    int l3f = 0;
    int l3g = l3;

    // Shifted lags for symmetric correlation.
    if (_type==Type.SYMMETRIC) {
      l1f = (l1>=0)?(l1+0)/2:(l1-1)/2;
      l1g = (l1>=0)?(l1+1)/2:(l1+0)/2;
      l2f = (l2>=0)?(l2+0)/2:(l2-1)/2;
      l2g = (l2>=0)?(l2+1)/2:(l2+0)/2;
      l3f = (l3>=0)?(l3+0)/2:(l3-1)/2;
      l3g = (l3>=0)?(l3+1)/2:(l3+0)/2;
    }

    // Scale factor so that center of window = 1.
    double scale1 = 1.0;
    double scale2 = 1.0;
    double scale3 = 1.0;
    if (_window==Window.GAUSSIAN) {
      scale1 *= sqrt(2.0*PI)*_sigma1;
      scale2 *= sqrt(2.0*PI)*_sigma2;
      scale3 *= sqrt(2.0*PI)*_sigma3;
    } else {
      scale1 *= 1.0+2.0*_sigma1;
      scale2 *= 1.0+2.0*_sigma2;
      scale3 *= 1.0+2.0*_sigma3;
    }

    // If symmetric correlation, need extra lag-dependent scaling.
    // This scaling accounts for the separation (by lag samples) of 
    // the two windows implicitly applied to f and g. The filter we
    // apply below to the correlation product h is the product of 
    // those two windows.
    if (_type==Type.SYMMETRIC) {
      if (_window==Window.GAUSSIAN) {
        scale1 *= exp((-0.125*l1*l1)/(_sigma1*_sigma1));
        scale2 *= exp((-0.125*l2*l2)/(_sigma2*_sigma2));
        scale3 *= exp((-0.125*l3*l3)/(_sigma3*_sigma3));
      } else {
        scale1 *= max(0.0,1.0+2.0*_sigma1-abs(l1))/(1.0+2.0*_sigma1);
        scale2 *= max(0.0,1.0+2.0*_sigma2-abs(l2))/(1.0+2.0*_sigma2);
        scale3 *= max(0.0,1.0+2.0*_sigma3-abs(l3))/(1.0+2.0*_sigma3);
      }
    }
    float scale = (float)(scale1*scale2*scale3);

    // Correlation product.
    float[][][] h = new float[n3][n2][n1];
    int i1min = max(0,l1f,-l1g);
    int i1max = min(n1,n1+l1f,n1-l1g);
    int i2min = max(0,l2f,-l2g);
    int i2max = min(n2,n2+l2f,n2-l2g);
    int i3min = max(0,l3f,-l3g);
    int i3max = min(n3,n3+l3f,n3-l3g);
    for (int i3=i3min; i3<i3max; ++i3) {
      float[][] f3 = f[i3-l3f];
      float[][] g3 = g[i3+l3g];
      float[][] h3 = h[i3];
      for (int i2=i2min; i2<i2max; ++i2) {
        float[] f32 = f3[i2-l2f];
        float[] g32 = g3[i2+l2g];
        float[] h32 = h3[i2];
        for (int i1=i1min; i1<i1max; ++i1) {
          h32[i1] = scale*f32[i1-l1f]*g32[i1+l1g];
        }
      }
    }

    // If Gaussian and symmetric and odd lag, delay (shift) by 1/2 sample.
    if (_window==Window.GAUSSIAN && _type==Type.SYMMETRIC) {
      if (l1f!=l1g) {
        shift1(h,c);
        Array.copy(c,h);
      }
      if (l2f!=l2g) {
        shift2(h,c);
        Array.copy(c,h);
      }
      if (l3f!=l3g) {
        shift3(h,c);
        Array.copy(c,h);
      }
    }

    // Filter correlation product with window. For symmetric correlations
    // with a rectangle window, the width of the product rectangle depends 
    // on the lag, so we construct a new rectangle filter for each lag.
    Filter f1 = _f1;
    Filter f2 = _f2;
    Filter f3 = _f3;
    if (_window==Window.RECTANGLE && _type==Type.SYMMETRIC) {
      f1 = new RectangleFilter(_sigma1,l1);
      f2 = new RectangleFilter(_sigma2,l2);
      f3 = new RectangleFilter(_sigma3,l3);
    }
    f1.apply1(h,c);
    Array.copy(c,h);
    f2.apply2(h,c);
    Array.copy(c,h);
    f3.apply3(h,c);
  }

  private void updateNormalize() {
    if (_dimension==0)
      return;
    int ns = (_type==Type.SIMPLE)?2:1;
    int n1 = max(1,_n1);
    int n2 = max(1,_n2);
    int n3 = max(1,_n3);
    _s = new float[ns][n3][n2][n1];
    if (_type==Type.SIMPLE) {
      if (_dimension==1) {
        float[] f = _f[0][0];
        float[] g = _g[0][0];
        float[] sf = _s[0][0][0];
        float[] sg = _s[1][0][0];
        correlate(0,f,f,sf);
        correlate(0,g,g,sg);
        Array.sqrt(sf,sf);
        Array.sqrt(sg,sg);
        Array.div(1.0f,sf,sf);
        Array.div(1.0f,sg,sg);
      } else if (_dimension==2) {
        float[][] f = _f[0];
        float[][] g = _g[0];
        float[][] sf = _s[0][0];
        float[][] sg = _s[1][0];
        correlate(0,0,f,f,sf);
        correlate(0,0,g,g,sg);
        Array.sqrt(sf,sf);
        Array.sqrt(sg,sg);
        Array.div(1.0f,sf,sf);
        Array.div(1.0f,sg,sg);
      } else {
        float[][][] f = _f;
        float[][][] g = _g;
        float[][][] sf = _s[0];
        float[][][] sg = _s[1];
        correlate(0,0,0,f,f,sf);
        correlate(0,0,0,g,g,sg);
        Array.sqrt(sf,sf);
        Array.sqrt(sg,sg);
        Array.div(1.0f,sf,sf);
        Array.div(1.0f,sg,sg);
      }
    } else {
      if (_dimension==1) {
        float[] f = _f[0][0];
        float[] g = _g[0][0];
        float[] s = _s[0][0][0];
        float[] sf = s;
        float[] sg = new float[_n1];
        correlate(0,f,f,sf);
        correlate(0,g,g,sg);
        Array.mul(sf,sg,s);
        Array.sqrt(s,s);
        Array.div(1.0f,s,s);
      } else if (_dimension==2) {
        float[][] f = _f[0];
        float[][] g = _g[0];
        float[][] s = _s[0][0];
        float[][] sf = s;
        float[][] sg = new float[_n2][_n1];
        correlate(0,0,f,f,sf);
        correlate(0,0,g,g,sg);
        Array.mul(sf,sg,s);
        Array.sqrt(s,s);
        Array.div(1.0f,s,s);
      } else {
        float[][][] f = _f;
        float[][][] g = _g;
        float[][][] s = _s[0];
        float[][][] sf = s;
        float[][][] sg = new float[_n3][_n2][_n1];
        correlate(0,0,0,f,f,sf);
        correlate(0,0,0,g,g,sg);
        Array.mul(sf,sg,s);
        Array.sqrt(s,s);
        Array.div(1.0f,s,s);
      }
    }
  }

  private static void shift(float[] f, float[] g) {
    int n1 = f.length;
    int i1b,i1e;

    // Rolling on.
    i1b = 0;
    i1e = min(4,n1);
    for (int i1=i1b; i1<i1e; ++i1) {
      int ib = max(0,4-i1);
      int ie = min(8,4-i1+n1);
      g[i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        g[i1] += S[i]*f[i1+i-4];
    }

    // Middle.
    i1b = 4;
    i1e = n1-3;
    for (int i1=i1b; i1<i1e; ++i1) {
      g[i1] = S4*(f[i1-4]+f[i1+3]) +
              S3*(f[i1-3]+f[i1+2]) +
              S2*(f[i1-2]+f[i1+1]) +
              S1*(f[i1-1]+f[i1  ]);
    }

    // Rolling off.
    i1b = max(0,n1-3);
    i1e = n1;
    for (int i1=i1b; i1<i1e; ++i1) {
      int ib = max(0,4-i1);
      int ie = min(8,4-i1+n1);
      g[i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        g[i1] += S[i]*f[i1+i-4];
    }
  }

  private static void shift1(float[][] f, float[][] g) {
    int n2 = f.length;
    for (int i2=0; i2<n2; ++i2)
      shift(f[i2],g[i2]);
  }

  private static void shift2(float[][] f, float[][] g) {
    int n2 = f.length;
    int n1 = f[0].length;
    int i2b,i2e;

    // Rolling on.
    i2b = 0;
    i2e = min(4,n2);
    for (int i2=i2b; i2<i2e; ++i2) {
      int ib = max(0,4-i2);
      int ie = min(8,4-i2+n2);
      for (int i1=0; i1<n1; ++i1)
        g[i2][i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        for (int i1=0; i1<n1; ++i1)
          g[i2][i1] += S[i]*f[i2+i-4][i1];
    }

    // Middle.
    i2b = 4;
    i2e = n2-3;
    for (int i2=i2b; i2<i2e; ++i2) {
      float[] g2 = g[i2];
      float[] fm4 = f[i2-4];
      float[] fm3 = f[i2-3];
      float[] fm2 = f[i2-2];
      float[] fm1 = f[i2-1];
      float[] fp0 = f[i2  ];
      float[] fp1 = f[i2+1];
      float[] fp2 = f[i2+2];
      float[] fp3 = f[i2+3];
      for (int i1=0; i1<n1; ++i1)
        g2[i1] = S4*(fm4[i1]+fp3[i1]) +
                 S3*(fm3[i1]+fp2[i1]) +
                 S2*(fm2[i1]+fp1[i1]) +
                 S1*(fm1[i1]+fp0[i1]);
    }

    // Rolling off.
    i2b = max(0,n2-3);
    i2e = n2;
    for (int i2=i2b; i2<i2e; ++i2) {
      int ib = max(0,4-i2);
      int ie = min(8,4-i2+n2);
      for (int i1=0; i1<n1; ++i1)
        g[i2][i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        for (int i1=0; i1<n1; ++i1)
          g[i2][i1] += S[i]*f[i2+i-4][i1];
    }
  }

  private static void shift1(float[][][] f, float[][][] g) {
    int n3 = f.length;
    for (int i3=0; i3<n3; ++i3)
      shift1(f[i3],g[i3]);
  }

  private static void shift2(float[][][] f, float[][][] g) {
    int n3 = f.length;
    for (int i3=0; i3<n3; ++i3)
      shift2(f[i3],g[i3]);
  }

  private static void shift3(float[][][] f, float[][][] g) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][] f2 = new float[n3][];
    float[][] g2 = new float[n3][];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
        f2[i3] = f[i3][i2];
        g2[i3] = g[i3][i2];
      }
      shift2(f2,g2);
    }
  }

  private static float[] makeGaussianWindow(double sigma) {
    int m = 1+2*(int)(4.0*sigma);
    int j = (m-1)/2;
    float[] w = new float[m];
    double s = -0.5/(sigma*sigma);
    for (int i=0; i<m; ++i) {
      double x = i-j;
      w[i] = (float)exp(s*x*x);
    }
    return w;
  }

  /**
   * Multiplies a specified array by a specified window. The window 
   * must have odd length, so that its center is uniquely defined.
   * Multiplication is limited to within the bounds of the specified arrays.
   * @param w the window array.
   * @param jf the index of the input array f at which to center the window.
   * @param f the input array.
   * @param jg the index of the output array g corresponding to index jf.
   * @param g the output array.
   */
  private static void applyWindow(
    float[] w, 
    int jf, float[] f, 
    int jg, float[] g) 
  {
    int nf = f.length;
    int ng = g.length;
    int nw = w.length;
    jf -= (nw-1)/2;
    jg -= (nw-1)/2;
    int imin = max(0,-jf,-jg);
    int imax = min(nw,nf-jf,ng-jg);
    for (int i=imin; i<imax; ++i)
      g[jg+i] = w[i]*f[jf+i];
  }

  /**
   * Multiplies a specified array by a specified separable window. The
   * multi-dimensional window is the product of multiple 1-D windows. All 
   * windows must have odd length, to make their centers uniquely defined.
   * Multiplication is limited to within the bounds of the specified arrays.
   * @param w1 the window array for the 1st dimension.
   * @param w2 the window array for the 2nd dimension.
   * @param j1f the index in the 1st dimension of the input array f at 
   *  which to center the window.
   * @param j2f the index in the 2nd dimension of the input array f at 
   *  which to center the window.
   * @param f the input array.
   * @param j1g the index of the output array g corresponding to index j1f.
   * @param j2g the index of the output array g corresponding to index j2f.
   * @param g the output array.
   */
  private static void applyWindow(
    float[] w1, float[] w2, 
    int j1f, int j2f, float[][] f, 
    int j1g, int j2g, float[][] g) 
  {
    int n1f = f[0].length;
    int n1g = g[0].length;
    int n1w = w1.length;
    j1f -= (n1w-1)/2;
    j1g -= (n1w-1)/2;
    int i1min = max(0,-j1f,-j1g);
    int i1max = min(n1w,n1f-j1f,n1g-j1g);
    int n2f = f.length;
    int n2g = g.length;
    int n2w = w2.length;
    j2f -= (n2w-1)/2;
    j2g -= (n2w-1)/2;
    int i2min = max(0,-j2f,-j2g);
    int i2max = min(n2w,n2f-j2f,n2g-j2g);
    for (int i2=i2min; i2<i2max; ++i2) {
      float w2i = w2[i2];
      float[] f2 = f[j2f+i2];
      float[] g2 = g[j2g+i2];
      for (int i1=i1min; i1<i1max; ++i1) {
        g2[j1g+i1] = w2i*w1[i1]*f2[j1f+i1];
      }
    }
  }

  private static void applyWindow(
    float[] w1, float[] w2, float[] w3,
    int j1f, int j2f, int j3f, float[][][] f, 
    int j1g, int j2g, int j3g, float[][][] g) 
  {
    int n1f = f[0].length;
    int n1g = g[0].length;
    int n1w = w1.length;
    j1f -= (n1w-1)/2;
    j1g -= (n1w-1)/2;
    int i1min = max(0,-j1f,-j1g);
    int i1max = min(n1w,n1f-j1f,n1g-j1g);
    int n2f = f.length;
    int n2g = g.length;
    int n2w = w2.length;
    j2f -= (n2w-1)/2;
    j2g -= (n2w-1)/2;
    int i2min = max(0,-j2f,-j2g);
    int i2max = min(n2w,n2f-j2f,n2g-j2g);
    int n3f = f.length;
    int n3g = g.length;
    int n3w = w3.length;
    j3f -= (n3w-1)/2;
    j3g -= (n3w-1)/2;
    int i3min = max(0,-j3f,-j3g);
    int i3max = min(n3w,n3f-j3f,n3g-j3g);
    for (int i3=i3min; i3<i3max; ++i3) {
      float w3i = w3[i3];
      float[][] f3 = f[j3f+i3];
      float[][] g3 = g[j3g+i3];
      for (int i2=i2min; i2<i2max; ++i2) {
        float w32i = w3i*w2[i2];
        float[] f32 = f3[j2f+i2];
        float[] g32 = g3[j2g+i2];
        for (int i1=i1min; i1<i1max; ++i1) {
          g32[j1g+i1] = w32i*w1[i1]*f32[j1f+i1];
        }
      }
    }
  }

  private void checkDimension(int dimension) {
    Check.state(_dimension==dimension,"dimension is valid");
  }

  private void checkDimensions(float[] c) {
    Check.argument(_n1==c.length,"array length is valid");
    checkDimension(1);
  }

  private void checkDimensions(float[][] c) {
    Check.argument(Array.isRegular(c),"c is regular");
    Check.argument(_n1==c[0].length,"array dimension 1 is valid");
    Check.argument(_n2==c.length,"array dimension 2 is valid");
    checkDimension(2);
  }

  private void checkDimensions(float[][][] c) {
    Check.argument(Array.isRegular(c),"c is regular");
    Check.argument(_n1==c[0][0].length,"array dimension 1 is valid");
    Check.argument(_n2==c[0].length,"array dimension 2 is valid");
    Check.argument(_n3==c.length,"array dimension 3 is valid");
    checkDimension(3);
  }

  // This interface makes it easier to implement different windows.
  private interface Filter {
    public void apply(float[] x, float[] y);
    public void apply1(float[][] x, float[][] y);
    public void apply2(float[][] x, float[][] y);
    public void apply1(float[][][] x, float[][][] y);
    public void apply2(float[][][] x, float[][][] y);
    public void apply3(float[][][] x, float[][][] y);
  }
  private class RectangleFilter implements Filter {
    public RectangleFilter(double sigma) {
      this(sigma,0);
    }
    public RectangleFilter(double sigma, int lag) {
      int n = (int)round(1+2*sigma);
      int m = max(0,(n-1-abs(lag))/2);
      int l = (lag%2==0)?-m:-m-1;
      _rrf = new RecursiveRectangleFilter(l,m);
    }
    public void apply(float[] x, float[] y) {
      _rrf.apply(x,y);
    }
    public void apply1(float[][] x, float[][] y) {
      _rrf.apply1(x,y);
    }
    public void apply2(float[][] x, float[][] y) {
      _rrf.apply2(x,y);
    }
    public void apply1(float[][][] x, float[][][] y) {
      _rrf.apply1(x,y);
    }
    public void apply2(float[][][] x, float[][][] y) {
      _rrf.apply2(x,y);
    }
    public void apply3(float[][][] x, float[][][] y) {
      _rrf.apply3(x,y);
    }
    private RecursiveRectangleFilter _rrf;
  }
  private class GaussianFilter implements Filter {
    public GaussianFilter(double sigma) {
      _rgf = new RecursiveGaussianFilter(sigma);
    }
    public void apply(float[] x, float[] y) {
      _rgf.apply0(x,y);
    }
    public void apply1(float[][] x, float[][] y) {
      _rgf.apply0X(x,y);
    }
    public void apply2(float[][] x, float[][] y) {
      _rgf.applyX0(x,y);
    }
    public void apply1(float[][][] x, float[][][] y) {
      _rgf.apply0XX(x,y);
    }
    public void apply2(float[][][] x, float[][][] y) {
      _rgf.applyX0X(x,y);
    }
    public void apply3(float[][][] x, float[][][] y) {
      _rgf.applyXX0(x,y);
    }
    private RecursiveGaussianFilter _rgf;
  }

  /**
   * Like {@link #apply(int,int,int,int,float[],float[],float[][])}, but
   * uses conventional windowing and FFTs to perform the cross-correlations.
   * Best for small numbers of cross-correlation windows.
   * The number of lags is nl1 = l1max-l1min+1.
   * @param l1min the minimum lag in the 1st dimension.
   * @param l1max the maximum lag in the 1st dimension.
   * @param j1c the sample index of the first correlation.
   * @param k1c the sample stride between correlations.
   * @param f the 1st input array; can be the same as g.
   * @param g the 2nd input array; can be the same as f.
   * @param c the output array; cannot be the same as f or g.
   */
  // Not yet tested!
  private void applyFft(
    int l1min, int l1max, int j1c, int k1c,
    float[] f, float[] g, float[][] c)
  {
    float[] w1 = makeGaussianWindow(_sigma1);
    int n1f = f.length;
    int n1c = c.length;
    int n1w = w1.length;
    int n1h = (n1w-1)/2;
    int n1l = l1max-l1min+1;
    int n1p = n1w+max(-l1min,n1l-1+l1min);
    int n1fft = FftReal.nfftFast(n1p);
    int n1pad = n1fft+2;
    FftReal fft = new FftReal(n1fft);
    float[] fpad = new float[n1pad];
    float[] gpad = new float[n1pad];
    int j1f = max(0, l1min);
    int j1g = max(0,-l1min);
    for (int i1c=0; i1c<n1c; ++i1c) {
      int m1c = j1c+i1c*k1c;
      Array.zero(fpad);
      Array.zero(gpad);
      applyWindow(w1,m1c,f,n1h+j1f,fpad);
      applyWindow(w1,m1c,g,n1h+j1g,gpad);
      fft.realToComplex(-1,fpad,fpad);
      fft.realToComplex(-1,gpad,gpad);
      for (int i1=0; i1<n1pad; i1+=2) {
        float fr = fpad[i1  ];
        float fi = fpad[i1+1];
        float gr = gpad[i1  ];
        float gi = gpad[i1+1];
        gpad[i1  ] = fr*gr+fi*gi;
        gpad[i1+1] = fr*gi-fi*gr;
      }
      fft.complexToReal(1,gpad,gpad);
      float s = 1.0f/(float)n1fft;
      float[] cc = c[i1c];
      for (int i1=0; i1<n1l; ++i1)
        cc[i1] = s*gpad[i1];
    }
  }

  // Not yet tested!
  private void applyFft(
    int l1min, int l1max, int j1c, int k1c,
    int l2min, int l2max, int j2c, int k2c,
    float[][] f, float[][] g, float[][][] c)
  {
    float[] w1 = makeGaussianWindow(_sigma1);
    float[] w2 = makeGaussianWindow(_sigma2);
    int j1f = max(0, l1min);
    int j1g = max(0,-l1min);
    int n1f = f[0].length;
    int n1g = g[0].length;
    int n1c = c[0].length;
    int n1w = w1.length;
    int n1h = (n1w-1)/2;
    int n1l = l1max-l1min+1;
    int n1p = n1w+max(-l1min,n1l-1+l1min);
    int n1fft = FftReal.nfftFast(n1p);
    int n1pad = n1fft+2;
    int j2f = max(0, l2min);
    int j2g = max(0,-l2min);
    int n2f = f.length;
    int n2g = g.length;
    int n2c = c.length;
    int n2w = w2.length;
    int n2h = (n2w-1)/2;
    int n2l = l2max-l2min+1;
    int n2p = n2w+max(-l2min,n2l-1+l2min);
    int n2fft = FftComplex.nfftFast(n2p);
    int n2pad = n2fft*2;
    FftReal fft1 = new FftReal(n1fft);
    FftComplex fft2 = new FftComplex(n2fft);
    float[][] fpad = new float[n2pad][n1pad];
    float[][] gpad = new float[n2pad][n1pad];
    for (int i2c=0; i2c<n2c; ++i2c) {
      int m2c = j2c+i2c*k2c;
      for (int i1c=0; i1c<n1c; ++i1c) {
        int m1c = j1c+i1c*k1c;
        Array.zero(fpad);
        Array.zero(gpad);
        applyWindow(w1,w2,m1c,m2c,f,n1h+j1f,n2h+j2f,fpad);
        applyWindow(w1,w2,m1c,m2c,g,n1h+j1g,n2h+j2g,gpad);
        fft1.realToComplex1(-1,n2p,fpad,fpad);
        fft1.realToComplex1(-1,n2p,gpad,gpad);
        fft2.complexToComplex2(-1,n1fft/2+1,fpad,fpad);
        fft2.complexToComplex2(-1,n1fft/2+1,gpad,gpad);
        for (int i2=0; i2<n2fft; ++i2) {
          float[] fpad2 = fpad[i2];
          float[] gpad2 = gpad[i2];
          for (int i1=0; i1<n1pad; i1+=2) {
            float fr = fpad2[i1  ];
            float fi = fpad2[i1+1];
            float gr = gpad2[i1  ];
            float gi = gpad2[i1+1];
            gpad2[i1  ] = fr*gr+fi*gi;
            gpad2[i1+1] = fr*gi-fi*gr;
          }
        }
        fft2.complexToComplex2(1,n1fft/2+1,gpad,gpad);
        fft1.realToComplex1(1,n2l,gpad,gpad);
        float s = 1.0f/((float)n1fft*(float)n2fft);
        float[] cc = c[i2c][i1c];
        for (int i2=0,ic=0; i2<n2l; ++i2) {
          float[] gpad2 = gpad[i2];
          for (int i1=0; i1<n1l; ++i1,++ic) {
            cc[ic] = s*gpad2[i1];
          }
        }
      }
    }
  }

  // Not yet tested!
  private void applyFft(
    int l1min, int l1max, int j1c, int k1c,
    int l2min, int l2max, int j2c, int k2c,
    int l3min, int l3max, int j3c, int k3c,
    float[][][] f, float[][][] g, float[][][][] c)
  {
    float[] w1 = makeGaussianWindow(_sigma1);
    float[] w2 = makeGaussianWindow(_sigma2);
    float[] w3 = makeGaussianWindow(_sigma3);
    int j1f = max(0, l1min);
    int j1g = max(0,-l1min);
    int n1f = f[0][0].length;
    int n1g = g[0][0].length;
    int n1c = c[0][0].length;
    int n1w = w1.length;
    int n1h = (n1w-1)/2;
    int n1l = l1max-l1min+1;
    int n1p = n1w+max(-l1min,n1l-1+l1min);
    int n1fft = FftReal.nfftFast(n1p);
    int n1pad = n1fft+2;
    int j2f = max(0, l2min);
    int j2g = max(0,-l2min);
    int n2f = f[0].length;
    int n2g = g[0].length;
    int n2c = c[0].length;
    int n2w = w2.length;
    int n2h = (n2w-1)/2;
    int n2l = l2max-l2min+1;
    int n2p = n2w+max(-l2min,n2l-1+l2min);
    int n2fft = FftComplex.nfftFast(n2p);
    int n2pad = n2fft*2;
    int j3f = max(0, l3min);
    int j3g = max(0,-l3min);
    int n3f = f.length;
    int n3g = g.length;
    int n3c = c.length;
    int n3w = w3.length;
    int n3h = (n3w-1)/2;
    int n3l = l3max-l3min+1;
    int n3p = n3w+max(-l3min,n3l-1+l3min);
    int n3fft = FftComplex.nfftFast(n3p);
    int n3pad = n3fft*2;
    FftReal fft1 = new FftReal(n1fft);
    FftComplex fft2 = new FftComplex(n2fft);
    FftComplex fft3 = new FftComplex(n3fft);
    float[][][] fpad = new float[n3pad][n2pad][n1pad];
    float[][][] gpad = new float[n3pad][n2pad][n1pad];
    for (int i3c=0; i3c<n3c; ++i3c) {
      int m3c = j3c+i3c*k3c;
      for (int i2c=0; i2c<n2c; ++i2c) {
        int m2c = j2c+i2c*k2c;
        for (int i1c=0; i1c<n1c; ++i1c) {
          int m1c = j1c+i1c*k1c;
          Array.zero(fpad);
          Array.zero(gpad);
          applyWindow(w1,w2,w3,m1c,m2c,m3c,f,n1h+j1f,n2h+j2f,j3f+n3h,fpad);
          applyWindow(w1,w2,w3,m1c,m2c,m3c,g,n1h+j1g,n2h+j2g,j3g+n3h,gpad);
          fft1.realToComplex1(-1,n2p,n3p,fpad,fpad);
          fft1.realToComplex1(-1,n2p,n3p,gpad,gpad);
          fft2.complexToComplex2(-1,n1fft/2+1,n3p,fpad,fpad);
          fft2.complexToComplex2(-1,n1fft/2+1,n3p,gpad,gpad);
          fft3.complexToComplex3(-1,n1fft/2+1,n2fft,fpad,fpad);
          fft3.complexToComplex3(-1,n1fft/2+1,n2fft,gpad,gpad);
          for (int i3=0; i3<n3fft; ++i3) {
            float[][] fpad3 = fpad[i3];
            float[][] gpad3 = gpad[i3];
            for (int i2=0; i2<n2fft; ++i2) {
              float[] fpad32 = fpad3[i2];
              float[] gpad32 = gpad3[i2];
              for (int i1=0; i1<n1pad; i1+=2) {
                float fr = fpad32[i1  ];
                float fi = fpad32[i1+1];
                float gr = gpad32[i1  ];
                float gi = gpad32[i1+1];
                gpad32[i1  ] = fr*gr+fi*gi;
                gpad32[i1+1] = fr*gi-fi*gr;
              }
            }
          }
          fft3.complexToComplex3(1,n1fft/2+1,n2fft,gpad,gpad);
          fft2.complexToComplex2(1,n1fft/2+1,n3l,gpad,gpad);
          fft1.realToComplex1(1,n2l,n3l,gpad,gpad);
          float s = 1.0f/((float)n1fft*(float)n2fft*(float)n3fft);
          float[] cc = c[i3c][i2c][i1c];
          for (int i3=0,ic=0; i3<n3l; ++i3) {
            float[][] gpad3 = gpad[i3];
            for (int i2=0; i2<n2l; ++i2) {
              float[] gpad32 = gpad3[i2];
              for (int i1=0; i1<n1l; ++i1,++ic) {
                cc[ic] = s*gpad32[i1];
              }
            }
          }
        }
      }
    }
  }

  private static void testQ2() {
    float a0 = 0.0f;
    float a1 = 0.0f;
    float a2 = 0.0f;
    float a3 = 0.0f;
    float a4 = 0.0f;
    float a5 = 0.0f;
    //float[][] ca = new float[3][3];
    float[][] ca = {
      { 0.705192f,  0.791807f,  0.669320f},
      { 0.679065f,  0.792324f,  0.696673f},
      { 0.649260f,  0.788553f,  0.721605f}
    };
    for (int k2=-1; k2<=1; ++k2) {
      float u2 = (float)k2;
      for (int k1=-1; k1<=1; ++k1) {
        float u1 = (float)k1;
        //float c = q2(u1,u2);
        //ca[k2+1][k1+1] = c;
        float c = ca[k2+1][k1+1];
        int k = (k1+1)+3*(k2+1);
        float[] ck = C2[k];
        a0 += ck[0]*c;
        a1 += ck[1]*c;
        a2 += ck[2]*c;
        a3 += ck[3]*c;
        a4 += ck[4]*c;
        a5 += ck[5]*c;
      }
    }
    Array.dump(ca);
    System.out.println("a0="+a0);
    System.out.println("a1="+a1);
    System.out.println("a2="+a2);
    System.out.println("a3="+a3);
    System.out.println("a4="+a4);
    System.out.println("a5="+a5);
    ca = new float[21][5];
    for (int i2=0; i2<21; ++i2) {
      for (int i1=0; i1<5; ++i1) {
        float u1 = (float)i1-2.0f;
        float u2 = (float)i2-10.0f;
        float c = a0+a1*u1+a2*u2+a3*u1*u2+a4*u1*u1+a5*u2*u2;
        ca[i2][i1] = c;
      }
    }
    Array.dump(ca);

    // Cholesky decomposition solves 2x2 system for refined lags.
    boolean pd = false;
    double aa0 = a0;
    double aa1 = a1;
    double aa2 = a2;
    double aa3 = a3;
    double aa4 = a4;
    double aa5 = a5;
    double w1 = 0.0;
    double w2 = 0.0;
    double b1 = aa1;
    double b2 = aa2;
    double a21 = -aa3;
    double a11 = -2.0*aa4;
    double a22 = -2.0*aa5;
    float[][] aa = {{(float)a11,(float)a21},{(float)a21,(float)a22}};
    Array.dump(aa);
    double d11 = a11;
    if (d11>0.0) {
      double l11 = sqrt(d11);
      double l21 = a21/l11;
      double d22 = a22-l21*l21;
      if (d22>0.0) {
        double l22 = sqrt(d22);
        double v1 = b1/l11;
        double v2 = (b2-l21*v1)/l22;
        w2 = v2/l22;
        w1 = (v1-l21*w2)/l11;
        pd = true;
      }
    }
    float u1 = (float)w1;
    float u2 = (float)w2;
    System.out.println("pd="+pd+" u1="+u1+" u2="+u2);
  }
  private static float q2(float u1, float u2) {
    // c00 + (u1-u1p)*c11*(u1-u1p) +
    //   2.0*(u1-u1p)*c12*(u2-u2p) +
    //       (u2-u2p)*c22*(u2-u2p)
    float c00 =  1.00f;
    float u1p =  0.00f;
    float u2p = -0.50f;
    float c11 = -0.90f;
    float c12 = -0.80f;
    float c22 = -0.90f;
    float d1 = u1-u1p;
    float d2 = u2-u2p;
    return c00+d1*c11*d1+2.0f*d1*c12*d2+d2*c22*d2;
  }
  public static void main(String[] args) {
    testQ2();
  }
}
