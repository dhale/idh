/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

//import static java.lang.Math.*;
import edu.mines.jtk.util.Check;
import static interp.SpecialMath.*;
import edu.mines.jtk.opt.BrentZeroFinder;

// Testing only.
import static edu.mines.jtk.util.ArrayMath.*;
import javax.swing.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * The Matern (Whittle-Matern) covariance function of distance.
 * <p>
 * This function is parameterized by three values: sigma, shape, and 
 * range. The range parameter is the so-called "effective range", so
 * that a change in shape does not change the extent of the function
 * significantly.
 *
 * @author Dave Hale, Colorado School of Mines.
 * @version 2012.09.15
 */
public class MaternCovariance implements Covariance {

  /**
   * Constructs the covariance for specified parameters.
   * @param sigma the sqrt(variance) for distance r = 0.
   * @param shape the shape; must be greater than zero.
   * @param range the effective range; must be greater than zero.
   */
  public MaternCovariance(double sigma, double shape, double range) {
    Check.argument(shape>0.0,"shape > 0");
    Check.argument(range>0.0,"range > 0");
    _shape = shape;
    _sigma = sigma;
    _range = range;
    init();
  }

  private void init() {
    double rmax =  16.0; // c(r)/c(0) < 0.001, for r > rmax
    double kmax = 256.0; // C(k)/C(0) < 0.001, for k > kmax
    int n = (int)(4.0*rmax*kmax); // oversample for linear interpolation
    int nfft = FftReal.nfftSmall(n);
    FftReal fft = new FftReal(nfft);
    int nr = nfft/2;
    int nk = nfft/2+1;
    double dr = rmax/(nr-1);
    double dk = 2.0*PI/(nfft*dr);
    double spower = -0.5-_shape;
    double kscale = 0.25/_shape;
    float[] c = new float[nfft+2];
    for (int ik=0; ik<nk; ++ik) {
      double k = ik*dk;
      c[2*ik] = (float)pow(1.0+k*k*kscale,spower);
    }
    fft.complexToReal(1,c,c);
    float[] cr = new float[nr];
    float cs = (float)(_sigma*_sigma/c[0]);
    for (int ir=0; ir<nr; ++ir)
      cr[ir] = cs*c[ir];
    _li = new LinearInterpolator();
    _li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    _li.setUniform(nr,dr*_range,0.0,cr);
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
    return _li.interpolate(r);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _shape;
  private double _sigma;
  private double _range;
  private LinearInterpolator _li;

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
    double sigma = 1.0;
    double range = 1.0;
    double rmin = 0.0;
    double rmax = 2*range;
    double cmin = 0.0;
    double cmax = sigma*sigma+0.1;
    int nr = 1001;
    double dr = (rmax-rmin)/(nr-1);
    double fr = rmin;
    Sampling sr = new Sampling(nr,dr,fr);
    double[] shapes = {0.2,0.5,1.0,1.5,2.0,8.0};
    int nshape = shapes.length;
    float[][] c = new float[nshape][nr];
    for (int ishape=0; ishape<nshape; ++ishape) {
      double shape = shapes[ishape];
      MaternCovariance mc = new MaternCovariance(sigma,shape,range);
      for (int ir=0; ir<nr; ++ir) {
        double r = sr.getValue(ir);
        c[ishape][ir] = (float)mc.evaluate(r);
      }
    }
    plot(sr,c,cmin,cmax,"All");
  }
  /**/
}
