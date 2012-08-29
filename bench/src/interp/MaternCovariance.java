/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

//import static java.lang.Math.*;

// Testing only.
import static edu.mines.jtk.util.ArrayMath.*;
import javax.swing.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * The Matern (Whittle-Matern) covariance function of distance.
 * <p>
 * This function is parameterized by three values: shape, sigma,
 * and range. When this function is evaluated, distances are
 * scaled by 2*sqrt(shape)/range, not simply 1/range. This scaling
 * enables variation of the shape parameter without changing the
 * range.
 *
 * @author Dave Hale, Colorado School of Mines.
 */
public class MaternCovariance implements Covariance {

  /**
   * Constructs the covariance for specified parameters.
   * @param shape the shape; must be an integer multiple of 1/2.
   * @param sigma the sqrt(variance) for distance r = 0.
   * @param range the effective range.
   */
  public MaternCovariance(double shape, double sigma, double range) {
    if (shape==(int)shape) {
      _n = (int)shape;
      _half = false;
    } else if ((shape-0.5)==(int)(shape-0.5)) {
      _n = (int)(shape-0.5);
      _half = true;
    } else {
      throw new IllegalArgumentException(
        "shape = integer multiple of 1/2 is required");
    }
    double gammav = _half?gamma2(_n):gamma(_n);
    _scalec = sigma*sigma/gammav/pow(2.0,shape-1.0);
    _scaler = 2.0*sqrt(shape)/range;
    _shape = shape;
    _sigma = sigma;
    _range = range;
  }

  /**
   * Gets the shape for this covariance function.
   * @return the shape.
   */
  public double getShape() {
    return _shape;
  }

  /**
   * Gets the sigma for this covariance function.
   * @return the sigma.
   */
  public double getSigma() {
    return _sigma;
  }

  /**
   * Gets the range for this covariance function.
   * @return the range.
   */
  public double getRange() {
    return _range;
  }

  public double evaluate(double r) {
    if (r==0.0)
      return _sigma*_sigma;
    double sr = _scaler*r;
    double kvr = _half?SpecialMath.kn2(_n,sr):SpecialMath.kn(_n,sr);
    return _scalec*pow(sr,_shape)*kvr;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _shape;
  private double _sigma;
  private double _range;
  private int _n;
  private boolean _half;
  private double _scalec;
  private double _scaler;

  private static double factorial(int n) {
    double g = 1.0;
    while (n>0)
      g *= n--;
    return g;
  }
  private static double gamma(int n) { // gamma(n)
    return factorial(n-1);
  }
  private static double gamma2(int n) { // gamma(n+1/2)
    return sqrt(PI)*factorial(2*n)/factorial(n)/pow(4.0,n);
  }

  /**/
  ///////////////////////////////////////////////////////////////////////////
  // test

  public static void trace(String s) {
    System.out.println(s);
  }

  public static void plot(Sampling sx, float[] y, 
    double ymin, double ymax, String title) 
  {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    PointsView pv = sp.addPoints(sx,y);
    sp.setVLimits(ymin,ymax);
    sp.setSize(790,700);
    sp.setTitle(title);
  }

  public static void plot(Sampling sx, float[][] y, 
    double ymin, double ymax, String title) 
  {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    for (int iy=0; iy<y.length; ++iy)
      sp.addPoints(sx,y[iy]);
    sp.setVLimits(ymin,ymax);
    sp.setSize(790,700);
    sp.setTitle(title);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
    public void run() {
      go();
    }});
  }

  public static void go() {
    double sigma = 2.0;
    double range = 2.0;
    double xmin = 0.0;
    double xmax = 2.0*sigma;
    double ymin = 0.0;
    double ymax = sigma*sigma+0.1;
    int nx = 1001;
    double dx = (xmax-xmin)/(nx-1);
    double fx = xmin;
    Sampling sx = new Sampling(nx,dx,fx);
    float[] y05 = new float[nx];
    float[] y10 = new float[nx];
    float[] y15 = new float[nx];
    float[] y20 = new float[nx];
    float[] ygg = new float[nx];
    MaternCovariance m05 = new MaternCovariance(0.5,sigma,range);
    MaternCovariance m10 = new MaternCovariance(1.0,sigma,range);
    MaternCovariance m15 = new MaternCovariance(1.5,sigma,range);
    MaternCovariance m20 = new MaternCovariance(2.0,sigma,range);
    for (int ix=0; ix<nx; ++ix) {
      float x = (float)sx.getValue(ix);
      y05[ix] = (float)m05.evaluate(x);
      y10[ix] = (float)m10.evaluate(x);
      y15[ix] = (float)m15.evaluate(x);
      y20[ix] = (float)m20.evaluate(x);
      ygg[ix] = (float)(sigma*sigma*exp(-x*x/(range*range)));
    }
    plot(sx,y05,ymin,ymax,"m05(x)");
    plot(sx,y10,ymin,ymax,"m10(x)");
    plot(sx,y15,ymin,ymax,"m15(x)");
    plot(sx,y20,ymin,ymax,"m20(x)");
    plot(sx,ygg,ymin,ymax,"Gaussian");
    plot(sx,new float[][]{y05,y10,y15,y20,ygg},ymin,ymax,"All");
  }
  /**/
}
