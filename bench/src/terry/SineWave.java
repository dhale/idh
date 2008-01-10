/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package terry;

import static java.lang.Math.*;

/**
 * An example of using the class Movie1 to animate a function f(x;t).
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.01.05
 */
public class SineWave implements Movie1.SampledFunction {

  public static void main(String[] args) {
    int nx = 101;
    double dx = 0.01;
    double dt = 1.0;
    SineWave sw = new SineWave(nx,dx,3);
    Movie1 movie = new Movie1(nx,dx,dt,sw);
  }

  /**
   * Constructs a sine wave sin(k*x).
   * The wavenumber k depends on the specified number of cycles.
   * @param nx number of x samples.
   * @param dx x sampling interval.
   * @param cycles number of cycles.
   */
  public SineWave(int nx, double dx, double cycles) {
    _nx = nx;
    _dx = dx;
    _k = 2.0*PI*cycles/(nx*dx);
  }

  /**
   * Initializes the sampled function f(x;t=0).
   * This implementation computes f(x;t=0) = sin(k*x).
   * @param f output array of sampled function values.
   */
  public void init(float[] f) {
    for (int ix=0; ix<_nx; ++ix) {
      double x = ix*_dx;
      f[ix] = (float)sin(_k*x);
    }
  }

  /**
   * Updates the sampled function f(x;t) for one time step.
   * Given f(x;t), computes f(x;t+dt).
   * This implementation models a sine wave moving to the right.
   * @param t the current time.
   * @param dt the time increment.
   * @param f input/output array of sampled function values.
   */
  public void step(double t, double dt, float[] f) {
    float flast = f[_nx-1];
    for (int ix=_nx-1; ix>0; --ix)
      f[ix] = f[ix-1];
    f[0] = flast;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nx; // number of x samples
  private double _dx; // x sampling interval
  private double _k; // wavenumber k in sin(k*x)
}
