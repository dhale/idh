/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Warps a sequence or image, with care taken to prevent aliasing.
 * <p>
 * Warping of a sequence x[i] is defined by y[i] = x(u[i]), where u[i] is a
 * warping sequence of values (with dimensionless units of samples) that need
 * not be integers. Warping may include any combination of shifting
 * (u[i]-u[i-1] = 1), stretching (u[i]-u[i-1]&lt;1), and squeezing
 * (u[i]-u[i-1]&gt;1).
 * <p>
 * If squeezing, anti-alias filtering is performed as follows. First, a factor
 * by which to temporarily oversample the output is set to the smaller of the
 * maximum value of u[i]-u[i-1] or a specified limit. Second, the warping
 * sequence u[i] is oversampled by this factor using piecewise-cubic and
 * monotonicity-preserving interpolation. Third, sinc interpolation is used to
 * compute an oversampled output y[i] = x(u[i]). Finally, an anti-alias filter
 * is applied before subsampling the output sequence y[i] so that the number
 * of output samples equals the number of input samples.
 * <p>
 * As part of the warping process, output values y[i] may be optionally scaled
 * by u[i]-u[i-1]. This amplitude scaling compensates for squeezing and
 * stretching.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.04.21
 */
public class WarpingFilter {

  /**
   * Sets an upper bound on oversampling used to prevent aliasing.
   * If less than or equal to 1.0, no oversampling will be performed.
   * The default bound is 10.0.
   * @param oversamplingLimit upper bound on the oversampling factor.
   */
  public void setOversamplingLimit(double oversamplingLimit) {
    _oversamplingLimit = (float)oversamplingLimit;
  }

  /**
   * Determines whether or not to scale amplitudes by du when warping.
   * The default is true.
   * @param amplitudeScaling true, for amplitude scaling; false, otherwise.
   */
  public void setAmplitudeScaling(boolean amplitudeScaling) {
    _amplitudeScaling = amplitudeScaling;
  }

  /**
   * Returns a warped sequence y[i] = x(u[i]).
   * The number of samples in the returned sequence y[i] equals the number of
   * values in the warping sequence u[i], which need not equal the number of
   * samples in the input sequence x[i].
   * @param u the warping sequence u[i].
   * @param x the input sequence x[i].
   * @return the warped sequence y[i].
   */
  public float[] apply(float[] u, float[] x) {
    int nu = u.length;
    int nx = x.length;
    int ny = nu;
    float r = 1.0f;
    if (_oversamplingLimit>1.0f)
      r = min(_oversamplingLimit,dumax(u));

    // If oversampling the output y to prevent aliasing, then
    // compute oversampled values u[i] at which to interpolate.
    if (r>1.0f) {
      nu = 1+(int)(r*(nu-1));
      u = monoResample(u,nu);
      //trace("r="+r+" nx="+nx+" ny="+ny+" nu="+nu);
    }

    // Sinc interpolation for y[i] = x(u[i]).
    float[] y = new float[nu];
    _si.interpolate(nx,1.0,0.0,x,nu,u,y);

    // Optional scaling of amplitudes by u'(t). If oversampling, then scaling
    // by r is necessary because the time sampling interval is now 1/r. So we
    // must multiply by r when approximating u'(t) with (u[i]-u[i-1])/(1/r).
    if (_amplitudeScaling) {
      y[0] *= r*(u[1]-u[0]);
      for (int i=1; i<nu; ++i)
        y[i] *= r*(u[i]-u[i-1]);
    }

    // If oversampled, apply anti-alias filtering before subsampling to obtain
    // the output sequence y[i].
    if (r>1.0f) {
      BandPassFilter aaf = new BandPassFilter(0.0,0.5/r,0.10/r,0.01);
      aaf.apply(y,y);
      y = sincResample(y,ny);
    }
    return y;
  }

  /**
   * Returns a warped image y[i2][i1] = x[i2](u[i2][i1]).
   * The dimensions of the output image y equal those of warping image u, and
   * need not equal the dimensions of the input image x.
   * @param u the warping image u[i2][i1].
   * @param x the input image x[i2][i1].
   * @return the warped image y[i2][i1].
   */
  public float[][] apply(float[][] u, float[][] x) {
    int n2 = u.length;
    float[][] y = new float[n2][];
    for (int i2=0; i2<n2; ++i2)
      y[i2] = apply(u[i2],x[i2]);
    return y;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final SincInterpolator _si =
    SincInterpolator.fromErrorAndFrequency(0.01,0.40);

  private float _oversamplingLimit = 10.0f;
  private boolean _amplitudeScaling = true;

  /**
   * Uniformly resamples a sequence x using piecewise cubic interpolation. The
   * interpolating function has continuous 1st derivatives and preserves any
   * monotonicity in the input sequence. The number of input and output 
   * samples equal the lengths of the specified input and output arrays x and
   * y, respectively. After resampling, y[0] = x[0], y[y.length-1] =
   * x[x.length-1], and all other samples will be uniformly spaced between
   * these end samples.
   */
  private void monoResample(float[] x, float[] y) {
    int nx = x.length;
    int ny = y.length;
    float dy = (nx-1.0f)/(ny-1.0f);
    float[] tx = rampfloat(0.0f,1.0f,nx);
    float[] ty = rampfloat(0.0f,dy,ny);
    CubicInterpolator ci = new CubicInterpolator(tx,x);
    ci.interpolate(ty,y);
    y[0] = x[0];
    y[ny-1] = x[nx-1];
  }
  private float[] monoResample(float[] x, int ny) {
    float[] y = new float[ny];
    monoResample(x,y);
    return y;
  }

  /**
   * Uniformly resamples a sequence x using sinc interpolation. The number of
   * input samples equals the length of the specified input array x. The
   * number of output samples equals the length of the specified output array
   * y. After resampling, y[0] = x[0], y[y.length-1] = x[x.length-1], and all
   * other samples will be uniformly spaced between these end samples.
   */
  private void sincResample(float[] x, float[] y) {
    int nx = x.length;
    int ny = y.length;
    double dy = (nx-1.0)/(ny-1.0);
    _si.interpolate(nx,1.0,0.0,x,ny,dy,0.0,y);
    y[0] = x[0];
    y[ny-1] = x[nx-1];
  }
  private float[] sincResample(float[] x, int ny) {
    float[] y = new float[ny];
    sincResample(x,y);
    return y;
  }

  /**
   * Returns the largest squeezing factor r = max (u[i]-u[i-1]). If less than
   * or equal to one, then no squeezing is implied by the warping sequence
   * u[i].
   */
  private float dumax(float[] u) {
    int n = u.length;
    float r = 0.0f;
    for (int i=1; i<n; ++i) {
      float du = u[i]-u[i-1];
      if (r<du)
        r = du;
    }
    return r;
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
