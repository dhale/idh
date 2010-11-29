/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Bilateral filtering.
 * @author Dave Hale, Colorado School of Mines
 * @version 2010.11.09
 */
public class BilateralFilter {

  /**
   * The range function of difference between input sample values x. 
   * Half-widths of all range functions are scaled so that they are 
   * comparable in their rejection of outliers.
   */
  public enum Type {

    /**
     * Gauss's function:
     * f(x) = exp(-(x^2)/(2.0*sigma^2)).
     * For this function, sigma = sigmaRange, where sigmaRange is the 
     * half-width of the range function for the bilateral filter.
     */
    GAUSS,

    /**
     * Huber's minmax (modified L1) function.
     * f(x) = 1/sigma, if |x|&lt;sigma; 1/|x|, otherwise.
     * For this function, sigma = sigmaRange, where sigmaRange is the 
     * half-width of the range function for the bilateral filter.
     */
    HUBER,

    /**
     * Tukey's biweight function.
     * f(x) = (1-(x^2)/(sigma^2))^2, if |x|&lt;sigma; 0, otherwise.
     * For this function, sigma = sigmaRange * sqrt(5), where sigmaRange 
     * is the half-width of the range function for the bilateral filter.
     */
    TUKEY
  }

  /**
   * Returns samples of the specified type of range function f(x).
   * @param type type of range function f(x).
   * @param sigma half-width of range function f(x).
   * @param sx sampling of x for which to compute f(x).
   * @return array of computed f(x).
   */
  public static float[] sampleRangeFunction(
    Type type, double sigma, Sampling sx) 
  {
    Fx fx = null;
    if (type==Type.GAUSS) {
      fx = new GaussFunction(sigma); 
    } else if (type==Type.HUBER) {
      fx = new HuberFunction(sigma); 
    } else if (type==Type.TUKEY) {
      fx = new TukeyFunction(sigma);
    }
    int nx = sx.getCount();
    float[] f = new float[nx];
    for (int ix=0; ix<nx; ++ix) {
      float xi = (float)sx.getValue(ix);
      f[ix] = fx.eval(xi);
    }
    return f;
  }

  /**
   * Constructs a bilateral filter with specified half-widths.
   * @param sigmaSpace half-width of spatial filter.
   * @param sigmaRange half-width of range function.
   */
  public BilateralFilter(double sigmaSpace, double sigmaRange) {
    setType(Type.TUKEY);
    setSigmaSpace(sigmaSpace);
    setSigmaRange(sigmaRange);
  }

  /**
   * Sets the half-width of the spatial filter.
   * @param sigmaSpace half-width of spatial filter.
   */
  public void setSigmaSpace(double sigmaSpace) {
    _sigmaS = sigmaSpace;
    _rgf = null;
  }

  /**
   * Sets the half-width of the range function.
   * @param sigmaSpace half-width of range function.
   */
  public void setSigmaRange(double sigmaRange) {
    _sigmaR = sigmaRange;
    _fx.setSigma(sigmaRange);
  }

  /**
   * Sets the type of range function.
   * The default type is Type.TUKEY.
   * @param type type of range function.
   */
  public void setType(Type type) {
    if (type==Type.GAUSS) {
      _fx = new GaussFunction(_sigmaR);
    } else if (type==Type.HUBER) {
      _fx = new HuberFunction(_sigmaR);
    } else if (type==Type.TUKEY) {
      _fx = new TukeyFunction(_sigmaR);
    }
  }

  /**
   * Applies this filter. The spatial part of the filter is shift-invariant.
   * @param x array of input samples.
   * @param y array of output samples.
   */
  public void apply(float[] x, float[] y) {
    _rgf = getRgf();
    apply(_rgf,_fx,x,y);
  }

  /**
   * Applies this filter. Scale factors modify the half-width sigmaSpace in 
   * the spatial part of the filter, so that it may not be shift-invariant.
   * @param s array of scale factors.
   * @param x array of input samples.
   * @param y array of output samples.
   */
  public void apply(float[] s, float[] x, float[] y) {
    _lsf = getLsf();
    _lsf.setFactors(s);
    apply(_lsf,_fx,x,y);
  }

  public void apply(float[][] x, float[][] y) {
    _rgf = getRgf();
    apply(_rgf,_fx,x,y);
  }

  public void apply(Tensors2 d, float[][] x, float[][] y) {
    apply(d,null,x,y);
  }

  public void apply(Tensors2 d, float[][] s, float[][] x, float[][] y) {
    _lsf = getLsf();
    _lsf.setTensors(d);
    _lsf.setFactors(s);
    apply(_lsf,_fx,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _sigmaS; // half-width of spatial filter
  private double _sigmaR; // half-width of range function
  private Fx _fx; // range function of input sample differences
  private Rgf _rgf; // recursive Gaussian filter
  private Lsf _lsf; // local smoothing filter
  
  // Range functions f(x) of differences between input sample values.
  private interface Fx {
    public float eval(float x);
    public float getSigma();
    public void setSigma(double sigma);
  }
  private static class GaussFunction implements Fx {
    public GaussFunction(double sigma) {
      setSigma(sigma);
    }
    public float eval(float x) {
      return (float)Math.exp(_scale*x*x);
    }
    public float getSigma() {
      return _sigma;
    }
    public void setSigma(double sigma) {
      _sigma = (float)sigma;
      _scale = -0.5f/(_sigma*_sigma);
    }
    private float _sigma;
    private float _scale; 
  }
  private static class HuberFunction implements Fx {
    public HuberFunction(double sigma) {
      setSigma(sigma);
    }
    public float eval(float x) {
      float ax = (x>=0.0f)?x:-x;
      return (ax<=_sigma)?1.0f:_sigma/ax;
    }
    public float getSigma() {
      return _sigma;
    }
    public void setSigma(double sigma) {
      _sigma = (float)(sigma);
    }
    private float _sigma;
  }
  private static class TukeyFunction implements Fx {
    public TukeyFunction(double sigma) {
      setSigma(sigma);
    }
    public float eval(float x) {
      float ax = (x>=0.0f)?x:-x;
      if (ax<_sigma) {
        float sx = ax*_scale;
        float tx = 1.0f-sx*sx;
        return tx*tx;
      } else {
        return 0.0f;
      }
    }
    public float getSigma() {
      return _sigma;
    }
    public void setSigma(double sigma) {
      _sigma = (float)(sigma*sqrt(5.0));
      _scale = 1.0f/_sigma;
    }
    private float _sigma;
    private float _scale; 
  }

  // Interfaces and abstract class for spatial filters.
  private interface F1 {
    public void apply(float[] x, float[] y);
  }
  private interface F2 {
    public void apply(float[][] x, float[][] y);
  }
  private interface F3 {
    public void apply(float[][][] x, float[][][] y);
  }
  private abstract class SpatialFilter implements F1,F2,F3 {
    public float sigma;
    public SpatialFilter(double sigma) {
      this.sigma = (float)sigma;
    }
    public abstract void apply(float[] x, float[] y);
    public abstract void apply(float[][] x, float[][] y);
    public abstract void apply(float[][][] x, float[][][] y);
  }

  // A recursive Gaussian spatial filter.
  private class Rgf extends SpatialFilter {
    public Rgf(double sigma) {
      super(sigma);
      _rgf = new RecursiveGaussianFilter(sigma);
    }
    public void apply(float[] x, float[] y) {
      _rgf.apply0(x,y);
    }
    public void apply(float[][] x, float[][] y) {
      _rgf.apply00(x,y);
    }
    public void apply(float[][][] x, float[][][] y) {
      _rgf.apply000(x,y);
    }
    private RecursiveGaussianFilter _rgf;
  }
  private Rgf getRgf() {
    if (_rgf==null || _sigmaS!=_rgf.sigma)
      _rgf = new Rgf(_sigmaS);
    return _rgf;
  }

  // A local smoothing spatial filter.
  private class Lsf extends SpatialFilter {
    public Lsf(double sigma) {
      super(sigma);
      _lsf = new LocalSmoothingFilter(0.01,200);
      _c = (float)(0.5*sigma*sigma);
    }
    public void setTensors(Tensors2 t2) {
      _t2 = t2;
    }
    public void setTensors(Tensors3 t3) {
      _t3 = t3;
    }
    public void setFactors(float[] s1) {
      _s1 = s1;
    }
    public void setFactors(float[][] s2) {
      _s2 = s2;
    }
    public void setFactors(float[][][] s3) {
      _s3 = s3;
    }
    public void apply(float[] x, float[] y) {
      _lsf.apply(_c,_s1,x,y);
    }
    public void apply(float[][] x, float[][] y) {
      _lsf.apply(_t2,_c,_s2,x,y);
    }
    public void apply(float[][][] x, float[][][] y) {
      _lsf.apply(_t3,_c,_s3,x,y);
    }
    private LocalSmoothingFilter _lsf;
    private Tensors2 _t2;
    private Tensors3 _t3;
    private float[] _s1;
    private float[][] _s2;
    private float[][][] _s3;
    private float _c;
  }
  private Lsf getLsf() {
    if (_lsf==null || _sigmaS!=_lsf.sigma)
      _lsf = new Lsf(_sigmaS);
    return _lsf;
  }

  // Computes sampling of x values. Bilateral filtering is performed by
  // linear interpolation of filter outputs for these sampled x values.
  private static Sampling samplingX(float[] x, Fx fx) {
    return samplingX(min(x),max(x),fx);
  }
  private static Sampling samplingX(float[][] x, Fx fx) {
    return samplingX(min(x),max(x),fx);
  }
  private static Sampling samplingX(float[][][] x, Fx fx) {
    return samplingX(min(x),max(x),fx);
  }
  private static Sampling samplingX(float xmin, float xmax, Fx fx) {
    float sigma = fx.getSigma();
    sigma = max(sigma,0.001f*(xmax-xmin));
    int nxs = 2+(int)((xmax-xmin)/sigma);
    float fxs = xmin;
    float dxs = (xmax-xmin)/(nxs-1);
    System.out.println("nxs="+nxs);
    return new Sampling(nxs,dxs,fxs);
  }

  // Apply bilateral filter for abstract spatial filters and range function.
  private static void apply(F1 f1, Fx fx, float[] x, float[] y) {
    int n1 = x.length;
    float[] yn = new float[n1];
    float[] yd = new float[n1];
    float[] tn = new float[n1];
    float[] td = y;
    Sampling sx = samplingX(x,fx);
    float dx = (float)sx.getDelta();
    int nx = sx.getCount();
    for (int kx=0; kx<nx; ++kx) {
      float xk = (float)sx.getValue(kx);
      scale(fx,xk,x,tn,td);
      f1.apply(copy(tn),tn);
      f1.apply(copy(td),td);
      accum(xk-dx,xk,xk+dx,dx,x,tn,td,yn,yd);
    }
    div(yn,yd,y);
  }
  private static void apply(F2 f2, Fx fx, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] yn = new float[n2][n1];
    float[][] yd = new float[n2][n1];
    float[][] tn = new float[n2][n1];
    float[][] td = y;
    Sampling sx = samplingX(x,fx);
    float dx = (float)sx.getDelta();
    int nx = sx.getCount();
    for (int kx=0; kx<nx; ++kx) {
      float xk = (float)sx.getValue(kx);
      scale(fx,xk,x,tn,td);
      f2.apply(copy(tn),tn);
      f2.apply(copy(td),td);
      accum(xk-dx,xk,xk+dx,dx,x,tn,td,yn,yd);
    }
    div(yn,yd,y);
  }
  private static void apply(F3 f3, Fx fx, float[][][] x, float[][][] y) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] yn = new float[n3][n2][n1];
    float[][][] yd = new float[n3][n2][n1];
    float[][][] tn = new float[n3][n2][n1];
    float[][][] td = y;
    Sampling sx = samplingX(x,fx);
    float dx = (float)sx.getDelta();
    int nx = sx.getCount();
    for (int kx=0; kx<nx; ++kx) {
      float xk = (float)sx.getValue(kx);
      scale(fx,xk,x,tn,td);
      f3.apply(copy(tn),tn);
      f3.apply(copy(td),td);
      accum(xk-dx,xk,xk+dx,dx,x,tn,td,yn,yd);
    }
    div(yn,yd,y);
  }
  private static void scale(Fx fx, float xk, 
    float[] x, float[] tn, float[] td) 
  {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1) {
      float xi = x[i1];
      float t = fx.eval(xi-xk);
      tn[i1] = t*xi;
      td[i1] = t;
    }
  }
  private static void scale(Fx fx, float xk, 
    float[][] x, float[][] tn, float[][] td)
  {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      scale(fx,xk,x[i2],tn[i2],td[i2]);
  }
  private static void scale(Fx fx, float xk, 
    float[][][] x, float[][][] tn, float[][][] td)
  {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      scale(fx,xk,x[i3],tn[i3],td[i3]);
  }

  // xa <= xb <= xc
  private static void accum(
    float xa, float xb, float xc, float dx, float[] x, 
    float[] tn, float[] td, float[] yn, float[] yd) 
  {
    int n1 = x.length;
    float sab = 1.0f/max(dx,xb-xa);
    float sbc = 1.0f/max(dx,xc-xb);
    for (int i1=0; i1<n1; ++i1) {
      float xi = x[i1];
      if (xa<xi && xi<xc) {
        float w = 1.0f-((xi<=xb)?(xb-xi)*sab:(xi-xb)*sbc);
        yn[i1] += w*tn[i1];
        yd[i1] += w*td[i1]; 
      }
    }
  }
  private static void accum(
    float xa, float xb, float xc, float dx, float[][] x, 
    float[][] tn, float[][] td, float[][] yn, float[][] yd) 
  {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      accum(xa,xb,xc,dx,x[i2],tn[i2],td[i2],yn[i2],yd[i2]);
  }
  private static void accum(
    float xa, float xb, float xc, float dx, float[][][] x, 
    float[][][] tn, float[][][] td, float[][][] yn, float[][][] yd) 
  {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      accum(xa,xb,xc,dx,x[i3],tn[i3],td[i3],yn[i3],yd[i3]);
  }
}
