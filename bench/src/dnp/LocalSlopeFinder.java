/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;

/**
 * Estimates local slopes of features in 2D and 3D images.
 * 
 * For a 2D image f(x1,x2), slope p2(x1,x2) is the ratio dx1/dx2 of 
 * linear features nearest the image sample at point (x1,x2). An 
 * estimate of the linearity (a number in [0,1]) of features nearest 
 * that sample may also be computed.
 * 
 * Likewise, for a 3D image f(x1,x2,x3), slopes p2(x1,x2,x3) and
 * p3(x1,x2,x3) are the ratios dx1/dx2 and dx1/dx3, respectively, 
 * of planar features nearest the image sample at point (x1,x2,x3). 
 * An estimate of the planarity (a number in [0,1]) of features 
 * nearest that sample may also be computed.
 *
 * All slopes are measured in unitless samples (dx1) per sample (dx2).
 * Minimum and maximum slopes may be specified, and default min and 
 * max bounds are -100 and 100 samples per sample, respectively.
 * Estimated slopes will be within the default or specified bounds.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.09
 */
public class LocalSlopeFinder {

  /**
   * Constructs a local slope finder with specified smoothing.
   * @param sigma1 half-width of smoother in 1st dimension.
   */
  public LocalSlopeFinder(double sigma1) {
    this(sigma1,100.0);
  }

  /**
   * Constructs a local slope finder with specified smoothing and bounds.
   * @param sigma1 half-width of smoother in 1st dimension.
   * @param pmax maximum slope returned by this slope finder.
   *  The minimum slope will be the negative of this maximum.
   */
  public LocalSlopeFinder(double sigma1, double pmax) {
    _sigma1 = (float)sigma1;
    setBounds(pmax);
  }

  /**
   * Constructs a local slope finder with specified smoothing and bounds.
   * @param sigma1 half-width of smoother in 1st dimension.
   * @param p2min minimum slope p2 returned by this slope finder.
   * @param p2max maximum slope p2 returned by this slope finder.
   * @param p3min minimum slope p3 returned by this slope finder.
   * @param p3max maximum slope p3 returned by this slope finder.
   */
  public LocalSlopeFinder(double sigma1,
    double p2min, double p2max,
    double p3min, double p3max) 
  {
    _sigma1 = (float)sigma1;
    setBounds(p2min,p2max,p3min,p3max);
  }

  /**
   * Sets the half-width of smoothers in 2nd and higher dimensions.
   * @param sigma2 the half-width.
   */
  public void setSigma2(double sigma2) {
    _sigma2 = _sigma3 = (float)sigma2;
  }

  /**
   * Sets bounds on slopes returned by this slope finder.
   * @param pmax maximum slope returned by this slope finder.
   *  The minimum slope will be the negative of this maximum.
   */
  public void setBounds(double pmax) {
    setBounds(-pmax,pmax,-pmax,pmax);
  }

  /**
   * Sets bounds on slopes returned by this slope finder.
   * @param p2min minimum slope p2 returned by this slope finder.
   * @param p2max maximum slope p2 returned by this slope finder.
   * @param p3min minimum slope p3 returned by this slope finder.
   * @param p3max maximum slope p3 returned by this slope finder.
   */
  public void setBounds(
    double p2min, double p2max,
    double p3min, double p3max) 
  {
    _p2min = (float)p2min;
    _p2max = (float)p2max;
    _p3min = (float)p3min;
    _p3max = (float)p3max;
  }
 
  /**
   * Finds slopes of features in the specified 2D image.
   * @param f array[n2][n1] of input image samples.
   * @param p2 array[n2][n1] of output slopes.
   */
  public void findSlopes(float[][] f, float[][] p2) {
    findSlopes(f,p2,null);
  }
 
  /**
   * Finds slopes of features in the specified 2D image.
   * Optionally estimates the linearities of image features.
   * @param f array[n2][n1] of input image samples.
   * @param p2 array[n2][n1] of output slopes.
   * @param el if not null, array[n2][n1] of output linearities.
   */
  public void findSlopes(float[][] f, float[][] p2, float[][] el) {
    int n1 = f[0].length;
    int n2 = f.length;

    // Normal vectors and linearities.
    float[][] u1 = new float[n2][n1];
    float[][] u2 = p2;
    LocalOrientFilter lof = new LocalOrientFilter(_sigma1,_sigma2);
    lof.applyForNormalLinear(f,u1,u2,el);

    // Compute slopes from normal vectors.
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        if (-u2i<_p2min*u1i) u2i = -_p2min*u1i;
        if (-u2i>_p2max*u1i) u2i = -_p2max*u1i;
        if (u1i==0.0f) {
          p2[i2][i1] = (u2i<0.0f)?_p2max:_p2min;
        } else {
          p2[i2][i1] = -u2i/u1i;
        }
      }
    }
  }
 
  /**
   * Finds slopes of features in the specified 3D image.
   * Optionally estimates the planarities of image features.
   * @param f array[n3][n2][n1] of input image samples.
   * @param p2 array[n3][n2][n1] of output slopes p2.
   * @param p3 array[n3][n2][n1] of output slopes p3.
   * @param ep if not null, array[n3][n2][n1] of output planarities.
   */
  public void findSlopes(
    float[][][] f, float[][][] p2, float[][][] p3, float[][][] ep) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // Normal vectors and linearities.
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = p2;
    float[][][] u3 = p3;
    LocalOrientFilter lof = new LocalOrientFilter(_sigma1,_sigma2,_sigma3);
    lof.applyForNormalPlanar(f,u1,u2,u3,ep);

    // Compute slopes from normal vectors.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          if (-u2i<_p2min*u1i) u2i = -_p2min*u1i;
          if (-u2i>_p2max*u1i) u2i = -_p2max*u1i;
          if (-u3i<_p3min*u1i) u3i = -_p3min*u1i;
          if (-u3i>_p3max*u1i) u3i = -_p3max*u1i;
          if (u1i==0.0f) {
            p2[i3][i2][i1] = (u2i<0.0f)?_p2max:_p2min;
            p3[i3][i2][i1] = (u3i<0.0f)?_p3max:_p3min;
          } else {
            p2[i3][i2][i1] = -u2i/u1i;
            p3[i3][i2][i1] = -u3i/u1i;
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1; // half-width of smoother in 1st dimension
  private float _sigma2 = 1.0f; // smoothing half-width in 2nd dimension
  private float _sigma3 = 1.0f; // smoothing half-width in 3rd dimension
  private float _p2min,_p2max; // min and max slopes in 2nd dimension
  private float _p3min,_p3max; // min and max slopes in 3rd dimension
}
