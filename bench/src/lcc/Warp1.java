package lcc;

import java.awt.*;
import java.util.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

import util.*;

public class Warp1 {

  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Warp1(args);
      }
    });
  }

  private int _fontSize = 24;
  private int _width = 900;
  private int _height = 800;

  private Warp1(String[] args) {

    // Small sigma (< delay) yields lots of bad peaks; bias ~ 0.4 sample
    //testDelay(0.15f,0.25f,4.5f,4.0f);

    // Small sigma (> delay) yields one bad peak; bias ~ 0.1 sample
    //testDelay(0.15f,0.25f,4.5f,6.0f);

    // Slightly larger sigma eliminates bad peak; bias ~ 0.1 sample
    testDelay(0.15f,0.25f,4.5f,10.0f);
    //testDelay(0.15f,0.25f,4.5f,8.0f);

    // Larger sigma has no bad peaks; bias < 0.01 sample
    //testDelay(0.15f,0.25f,4.5f,32.0f);

    //testWarp(0.15f,0.25f,4.5f,10.0f);
  }

  /**
   * Tests warp of filtered random sequence.
   * @param fl frequency at which filter amplitude response ~ 0.99
   * @param fh frequency at which filter amplitude response ~ 0.01
   * @param delay the maximum delay such that g(t) = f(t-delay).
   * @param sigma half-width of Gaussian window.
   */
  private void testWarp(float fl, float fh, float delay, float sigma) {
    int n = 315;
    int lmax =  8;
    int lmin = -lmax;
    int nlag = 1+lmax-lmin;
    float[] f = makeRandom(n,fl,fh);
    float[][] xyd = computeWarpFunctions(n);
    float[] x = xyd[0];
    float[] y = xyd[1];
    float[] d = xyd[2];
    float[] g = warp(x,f);
    float[][] c = new float[nlag][n];
    float[] u = new float[n];
    lcc(sigma,lmin,lmax,f,g,c,u);
    //f = applyGaussianWindow(79,sigma,f);
    //g = applyGaussianWindow(79,sigma,g);
    plot(f,g,c,u,d);
  }

  /**
   * Tests constant delay of filtered random sequence.
   * @param fl frequency at which filter amplitude response ~ 0.99
   * @param fh frequency at which filter amplitude response ~ 0.01
   * @param delay the delay such that g(t) = f(t-delay).
   * @param sigma half-width of Gaussian window.
   */
  private void testDelay(float fl, float fh, float delay, float sigma) {
    int n = 315;
    int lmax =  8;
    int lmin = -lmax;
    int nlag = 1+lmax-lmin;
    //float[] f = makeRandom(n,fl,fh);
    float[] f = makeSweep(n,0.5f*fl,0.5f*fl);
    float[] g = delay(delay,f);
    float[][] c = new float[nlag][n];
    float[] u = new float[n];
    lcc(sigma,lmin,lmax,f,g,c,u);
    float[] d = Array.fillfloat(delay,n);
    //f = applyGaussianWindow(79,sigma,f);
    //g = applyGaussianWindow(79,sigma,g);
    plot(f,g,c,u,d);
  }

  private static float gauss(double s, double x) {
    double e = x/s;
    return (float)exp(-0.5f*e*e);
  }

  private static float[] applyGaussianWindow(double x, double s, float[] f) {
    int n = f.length;
    float[] g = new float[n];
    for (int i=0; i<n; ++i)
      g[i] = f[i]*gauss(s,i-x);
    return g;
  }

  private float[] makeRandom(int n, double fl, double fh) {

    // Random sequence in [-1,1].
    float[] f = Array.randfloat(new Random(314159),n);
    f = Array.sub(Array.mul(2.0f,f),1.0f);

    // Lowpass filter.
    ButterworthFilter lp = new ButterworthFilter(0.01,0.01,0.02,0.99);
    lp.applyForward(f,f);

    // Highpass filter.
    ButterworthFilter hp = new ButterworthFilter(fl,0.99,fh,0.01);
    hp.applyForward(f,f);

    return f;
  }

  private float[] makeSweep(int n, double fl, double fh) {
    float[] f = new float[n];
    for (int i=0; i<n; ++i) {
      double phase = 2.0*DBL_PI*(fl+0.5*i*(fh-fl)/(n-1))*i;
      f[i] = (float)cos(phase);
    }
    return f;
  }

  private float[] delay(double d, float[] f) {
    int n = f.length;
    float[] g = new float[n];
    SincInterpolator si = new SincInterpolator();
    si.setUniform(n,1.0,0.0,f);
    si.interpolate(n,1.0,-d,g);
    return g;
  }

  private void lcc(
    double sigma, int lmin, int lmax, float[] f, float[] g,
    float[][] c, float[] u)
  {
    int n = f.length;
    int nlag = 1+lmax-lmin;
    LocalCorrelationFilter.Type type = 
      LocalCorrelationFilter.Type.SIMPLE;
    LocalCorrelationFilter.Window window = 
      LocalCorrelationFilter.Window.RECTANGLE;
    LocalCorrelationFilter lcf = 
      new LocalCorrelationFilter(type,window,sigma);
    lcf.setInputs(f,g);
    for (int lag=lmin; lag<=lmax; ++lag) {
      int k = lag-lmin;
      lcf.correlate(lag,c[k]);
      lcf.normalize(lag,c[k]);
    }
    byte[] l = new byte[n];
    lcf.findMaxLags(lmin,lmax,l);
    lcf.refineLags(l,u);
  }

  private void plot(float[] f, float[] g, float[][] c, float[] u, float[] d) {
    int n = f.length;
    Sampling s = new Sampling(n,1.0,0.0);
    int nlag = c.length;
    Sampling slag = new Sampling(nlag,1.0,-(nlag-1)/2);
    PlotPanel panel = new PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP);
    SequenceView fv = panel.addSequence(0,0,s,f);
    SequenceView gv = panel.addSequence(1,0,s,g);
    PixelsView cv = panel.addPixels(2,0,s,slag,c);
    cv.setColorModel(ColorMap.JET);
    PointsView dv = panel.addPoints(2,0,s,d);
    dv.setLineColor(Color.WHITE);
    dv.setLineStyle(PointsView.Line.DOT);
    PointsView uv = panel.addPoints(2,0,s,u);
    uv.setLineColor(Color.WHITE);
    panel.setHLabel("sample");
    panel.setVLabel(0,"f");
    panel.setVLabel(1,"g");
    panel.setVLabel(2,"lag");
    frame(panel);
  }

  private PlotFrame frame(PlotPanel panel) {
    PlotFrame frame = new PlotFrame(panel);
    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    frame.setFontSize(_fontSize);
    frame.setSize(_width,_height);
    frame.setVisible(true);
    return frame;
  }

  /**
   * Computes functions that define warping, unwarping and displacement.
   * Warping is q(y) = p(x(y)). Unwarping is p(x) = q(y(x)). The displacement
   * implied by warping is e = y(x)-x. (The displacement implied by unwarping
   * is x(y)-y, but is not computed.) All functions are uniformly-sampled.
   */
  private static float[][] computeWarpFunctions(int n) {
    float a = 0.20f;
    float b = (float)((n-1)/2);
    float s = 50.0f;
    float[] x = warpGauss(a,b,s,n);
    float[] y = unwarpGauss(a,b,s,n);
    float[] z = Array.rampfloat(0.0f,1.0f,n);
    float[] d = Array.sub(y,z);
    return new float[][]{x,y,d};
  }

  /**
   * Adjusts estimated displacements for symmetric correlations.
   */
  private static void adjustu(float[] u) {
    int n = u.length;
    float[] x = Array.rampfloat(0.0f,1.0f,n);
    CubicInterpolator.Method method = CubicInterpolator.Method.LINEAR;
    CubicInterpolator ci = new CubicInterpolator(method,n,x,u);
    for (int i=0; i<n; ++i) {
      float xa = x[i];
      for (int iter=0; iter<5; ++iter) {
        xa = x[i]+0.5f*ci.interpolate(xa);
      }
      //if ((xa-x[i])*u[i]>0.0f)
        u[i] = ci.interpolate(xa);
    }
  }

  /**
   * Uses interpolation to compute q(i) = p(x[i]).
   */
  private static float[] warp(float[] x, float[] p) {
    int n = p.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniform(n,1.0,0.0,p);
    float[] q = new float[n];
    si.interpolate(n,x,q);
    return q;
  }

  /**
   * Computes a uniformly-sampled warping function based on a gaussian. 
   * The warping function is x = x(y) for uniformly-sampled y. The 
   * function x is defined so that x equals y when y equals b, 
   * x is less than y for y less than b, and x is greater than y for y 
   * greater than b. Thus the parameter b controls the location of 
   * displacment. The parameter a controls the maximum displacement, and 
   * the parameter s controls the spatial extent of displacement.
   */
  private static float[] warpGauss(float a, float b, float s, int n) {
    float[] x = new float[n];
    for (int i=0; i<n; ++i) {
      float y = (float)i-b;
      x[i] = b+y*(1.0f+a*gauss(s,y));
    }
    return x;
  }

  /**
   * Computes a uniformly-sampled unwarping function based on a gaussian. 
   * The unwarping function is y = y(x). In other words, the unwarping 
   * function y(x) is the inverse of the warping function x(y). The 
   * inverse y(x) is computed by iteration, because x(y) and y(x) are 
   * transcendental functions.
   */
  private static float[] unwarpGauss(float a, float b, float s, int n) {
    float[] y = new float[n];
    for (int i=0; i<n; ++i) {
      float x = (float)i-b;
      float yi = x;
      float yp;
      do {
        yp = yi;
        yi = x/(1.0f+a*gauss(s,yp));
      } while (abs(yi-yp)>0.001f);
      y[i] = yi+b;
    }
    return y;
  }
}
