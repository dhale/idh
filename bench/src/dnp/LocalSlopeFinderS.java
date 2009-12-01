/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

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
 * @version 2009.09.23
 */
public class LocalSlopeFinderS {

  /**
   * Constructs a local slope finder with specified smoothing.
   * @param sigma1 half-width of smoother in 1st dimension.
   * @param sigma2 half-width of smoother in 2nd dimension.
   */
  public LocalSlopeFinderS(double sigma1, double sigma2) {
    this(sigma1,sigma2,100.0);
  }

  /**
   * Constructs a local slope finder with specified smoothing and bounds.
   * @param sigma1 half-width of smoother in 1st dimension.
   * @param sigma2 half-width of smoother in 2nd dimension.
   * @param pmax maximum slope returned by this slope finder.
   *  The minimum slope will be the negative of this maximum.
   */
  public LocalSlopeFinderS(double sigma1, double sigma2, double pmax) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    setBounds(pmax);
  }

  /**
   * Constructs a local slope finder with specified smoothing and bounds.
   * @param sigma1 half-width of smoother in 1st dimension.
   * @param sigma2 half-width of smoother in 2nd dimension.
   * @param p2min minimum slope p2 returned by this slope finder.
   * @param p2max maximum slope p2 returned by this slope finder.
   * @param p3min minimum slope p3 returned by this slope finder.
   * @param p3max maximum slope p3 returned by this slope finder.
   */
  public LocalSlopeFinderS(
    double sigma1, double sigma2,
    double p2min, double p2max,
    double p3min, double p3max) 
  {
    _sigma1 = (float)sigma1;
    setBounds(p2min,p2max,p3min,p3max);
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

  public void filter(
    int dir, float[][] f, float[][] g, float[][] p2, float[][] el)
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float sigmaMin = sqrt(2.0f)/2.0f;
    LocalCorrelationFilter lcf = 
      new LocalCorrelationFilter(
        LocalCorrelationFilter.Type.SIMPLE,
        LocalCorrelationFilter.Window.GAUSSIAN,
        _sigma1);
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[] gt = zerofloat(n1);
    float[] u = zerofloat(n1);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float a = (_sigma2>sigmaMin)?(_sigma2-sigmaMin)/(_sigma2+sigmaMin):0.0f;
    float b = 1.0f-a;
    int lmin = (int)(_p2min-1.0f);
    int lmax = (int)(_p2max+1.0f);
    int i2b = (dir>0)?0:n2-1;
    int i2e = (dir>0)?n2:-1;
    int i2d = (dir>0)?1:-1;
    float pscale = (dir>0)?-1.0f:1.0f;
    copy(f[i2b],g[i2b]);
    for (int i2=i2b+i2d; i2!=i2e; i2+=i2d) {
      findShifts(lcf,-lmax,-lmin,f[i2],g[i2-i2d],u,el[i2]);
      //fill(1.0f,el[i2]); zero(u); zero(p2[i2]);
      mul(pscale,u,p2[i2]);
      add(r,u,u);
      si.setUniformSamples(g[i2-i2d]);
      si.interpolate(n1,u,gt);
      for (int i1=0; i1<n1; ++i1) {
        //float e = pow(el[i2][i1],2.0f);
        float e = el[i2][i1];
        float c = a*e;
        g[i2][i1] = c*gt[i1]+(1.0f-c)*f[i2][i1];
        el[i2][i1] = e;
      }
    }
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
    float sigmaMin = sqrt(2.0f)/2.0f;
    LocalCorrelationFilter lcf = 
      new LocalCorrelationFilter(
        LocalCorrelationFilter.Type.SIMPLE,
        LocalCorrelationFilter.Window.GAUSSIAN,
        _sigma1);
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[][] gm = new float[n2][n1];
    float[][] gp = new float[n2][n1];
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] gt = zerofloat(n1);
    float[] cmp = zerofloat(n1);
    float[] cpm = zerofloat(n1);
    float[] ump = zerofloat(n1);
    float[] upm = zerofloat(n1);
    float[] u = ump;
    float[] c = cmp;
    float[] t = cpm;

    // Parameters for left-right and right-left filters.
    float a = (_sigma2>sigmaMin)?(_sigma2-sigmaMin)/(_sigma2+sigmaMin):0.0f;
    float b = 1.0f-a;
    int lmin = (int)(_p2min-1.0f);
    int lmax = (int)(_p2max+1.0f);

    // Filter from left to right.
    copy(f[0],gm[0]);
    for (int i2=1; i2<n2; ++i2) {
      findShifts(lcf,-lmax,-lmin,f[i2],gm[i2-1],u,c);
      add(r,u,t);
      si.setUniformSamples(gm[i2-1]);
      si.interpolate(n1,t,gt);
      for (int i1=0; i1<n1; ++i1)
        gm[i2][i1] = a*gt[i1]+b*f[i2][i1];
    }

    // Filter from right to left.
    copy(f[n2-1],gp[n2-1]);
    for (int i2=n2-2; i2>=0; --i2) {
      findShifts(lcf,lmin,lmax,f[i2],gp[i2+1],u,c);
      add(r,u,t);
      si.setUniformSamples(gp[i2+1]);
      si.interpolate(n1,t,gt);
      for (int i1=0; i1<n1; ++i1)
        gp[i2][i1] = a*gt[i1]+b*f[i2][i1];
    }

    // Use the two filtered images to find slopes.
    lmin = (int)(2.0f*_p2min-1.0f);
    lmax = (int)(2.0f*_p2max+1.0f);
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(0,i2-1);
      int i2p = min(i2+1,n2-1);
      findShifts(lcf, lmin, lmax,gm[i2m],gp[i2p],ump,cmp);
      findShifts(lcf,-lmax,-lmin,gp[i2p],gm[i2m],upm,cpm);
      float scale = 0.5f/(i2p-i2m);
      for (int i1=0; i1<n1; ++i1)
        p2[i2][i1] = scale*(ump[i1]-upm[i1]);
        //p2[i2][i1] = 0.5f*ump[i1]; 
      if (el!=null) {
        for (int i1=0; i1<n1; ++i1)
          el[i2][i1] = 0.5f*(cmp[i1]+cpm[i1]);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1; // half-width of smoother in 1st dimension
  private float _sigma2; // half-width of smoother in 2nd and 3rd dimension
  private float _p2min,_p2max; // min and max slopes in 2nd dimension
  private float _p3min,_p3max; // min and max slopes in 3rd dimension

  private static void findShifts(
    LocalCorrelationFilter lcf, int min, int max,
    float[] f, float[] g, float[] u, float[] c)
  {
    int n1 = f.length;

    // Default shifts and correlation coefficients are zero.
    zero(u);
    zero(c);

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][] cbuf = new float[3][n1];

    // Correlate for min lag.
    lcf.setInputs(f,g);
    int lag1 = min;
    lcf.correlate(lag1,cbuf[1]);
    lcf.normalize(lag1,cbuf[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[] ca = (lag>min)?cbuf[(i  )%3]:cbuf[(i+2)%3];
      float[] cb =           cbuf[(i+1)%3];
      float[] cc = (lag<max)?cbuf[(i+2)%3]:cbuf[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = lag+1;
        lcf.correlate(lag1,cc);
        lcf.normalize(lag1,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, update the correlation maximum value and displacement
      // using quadratic interpolation of three correlation values.
      for (int i1=0; i1<n1; ++i1) {
        float ai = ca[i1];
        float bi = cb[i1];
        float ci = cc[i1];
        if (bi>=ai && bi>=ci) {
          double c0 = bi;
          double c1 = 0.5*(ci-ai);
          double c2 = 0.5*(ci+ai)-bi;
          double up = (c2<0.0)?-0.5*c1/c2:0.0;
          double cp = c0+up*(c1+up*c2);
          if (cp>c[i1]) {
            c[i1] = (float)cp;
            u[i1] = (float)(lag+up);
          }
        }
      }
    }
  }
}
