/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package terry;

import static java.lang.Math.*;

/**
 * An example of using the class Movie2 to animate a function f(x,y;t).
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.01.05
 */
public class PlaneWave implements Movie2.SampledFunction {

  public static void main(String[] args) {
    int nx = 101;
    int ny = 101;
    double dx = 0.01;
    double dy = 0.01;
    double dt = 1.0;
    double angle = 30.0; // degrees clockwise from horizontal
    double cycles = 6.0;
    PlaneWave pw = new PlaneWave(nx,dx,ny,dy,angle,cycles);
    Movie2 movie = new Movie2(nx,dx,ny,dy,dt,pw);
  }

  /**
   * Constructs a sine wave sin(k*x).
   * The wavenumber k depends on the specified number of cycles.
   * @param nx number of x samples.
   * @param dx x sampling interval.
   * @param ny number of y samples.
   * @param dy y sampling interval.
   * @param angle angle of wave propagation, in degrees.
   * @param cycles number of cycles.
   */
  public PlaneWave(
    int nx, double dx, int ny, double dy, 
    double angle, double cycles) {
    _nx = nx;
    _dx = dx;
    _ny = ny;
    _dy = dy;
    double a = angle*PI/180.0;
    double xs = nx*dx;
    double ys = ny*dy;
    double ds = sqrt(xs*xs+ys*ys);
    double velocity = min(dx,dy);
    double k = 2.0*PI*cycles/ds;
    _kx = k*cos(a);
    _ky = k*sin(a);
    _w = k*velocity;
  }

  /**
   * Initializes the sampled function f(x,y;t=0).
   * @param f output array of sampled function values.
   */
  public void init(float[][] f) {
    computeWaves(0.0,f);
  }

  /**
   * Updates the sampled function f(x;t) for one time step.
   * Given f(x,y;t), computes f(x,y;t+dt).
   * @param t the current time.
   * @param dt the time increment.
   * @param f input/output array of sampled function values.
   */
  public void step(double t, double dt, float[][] f) {
    computeWaves(t+dt,f);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nx,_ny; // number of x samples
  private double _dx,_dy; // x sampling interval
  private double _kx,_ky; // wavenumbers kx and ky
  private double _w; // frequency w

  private void computeWaves(double t, float[][] f) {
    for (int iy=0; iy<_ny; ++iy) {
      double y = iy*_dy;
      for (int ix=0; ix<_nx; ++ix) {
        double x = ix*_dx;
        f[iy][ix] = (float)sin(_kx*x+_ky*y-_w*t);
      }
    }
  }
}
