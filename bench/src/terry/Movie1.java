/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package terry;

import static java.lang.Math.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * Displays a function f(x;t) as a movie for sampled times t.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.01.05
 */
public class Movie1 extends Movie {

  /**
   * A sampled function f(x;t).
   */
  public interface SampledFunction {

    /**
     * Computes initial function values f(x;t=0).
     * @param f output array of function values.
     */
    public void init(float[] f);

    /**
     * Given f(x;t), computes function values f(x;t+dt).
     * @param t the current time.
     * @param dt the time increment.
     * @param f input/output array of function values.
     */
    public void step(double t, double dt, float[] f);
  }

  /**
   * Constructs a movie for a specified sampling function.
   * @param nx the number x samples.
   * @param dx the x sampling interval.
   * @param dt the t sampling interval.
   * @param sf the sampled function.
   */
  public Movie1(int nx, double dx, double dt, SampledFunction sf) {
    super(dt);
    _nx = nx;
    _dx = dx;
    _sx = new Sampling(_nx,_dx,0.0);
    _f = new float[nx];
    _sf = sf;
    initFrame(950,400);
  }

  /**
   * Sets the range of function values displayed. By default, this range 
   * is computed automatically so that all function values are visible.
   * @param fmin the minimum function valued displayed.
   * @param fmax the maximum function valued displayed.
   */
  public void setRange(float fmin, float fmax) {
    _panel.setVLimits(fmin,fmax);
  }

  /**
   * Sets automatic calculation of the range of function values displayed.
   * The range will be computed dynamically to include the current minimum 
   * and maximum function values.
   */
  public void setAutomaticRange() {
    _panel.setVLimitsDefault();
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void initSampledFunction() {
    _sf.init(_f);
  }
  protected void stepSampledFunction(double t, double dt) {
    _sf.step(t,dt,_f);
  }
  protected PlotPanel initPanel() {
    _panel = new PlotPanel(1,1);
    _panel.setHLabel("x");
    _panel.setHFormat("%1.6f");
    _panel.setVLabel("f(x,t)");
    _view = _panel.addPoints(_sx,_f);
    return _panel;
  }
  protected void updateView() {
    _view.set(_sx,_f);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nx;
  private double _dx;
  private Sampling _sx;
  private SampledFunction _sf;
  private float[] _f;
  private PointsView _view;
  private PlotPanel _panel;
}
