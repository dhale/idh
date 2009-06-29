/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fd;

import java.awt.*;
import javax.swing.*;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.ArrayMath;
import static edu.mines.jtk.util.MathPlus.exp;

/**
 * Demonstrates finite-difference approximations of the 1-D wave equation.
 * <p>
 * Demonstrates three different methods for approximating spatial
 * derivatives. The different methods correspond to different schemes 
 * for handling the factor 1/density in between the 1st and 2nd spatial
 * derivative. For all methods, discretization errors are 2nd-order in 
 * both space and time.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.24
 */
public class Wave1 {

  /**
   * Method used to approximate the term with d/dx 1/density d/dx.
   * The symmetric method uses a symmetric negative-semidefinite
   * approximation. This method is simplest and most accurate.
   * Two product methods are implemented for comparison. The first 
   * product method uses finite difference approximations to 
   * 1/density d^2/dx^2 + (d/dx 1/density)*d/dx.
   * The second product uses the same approximations for
   * 1/density d^2/dx^2 - 1/density^2 (d/dx density)*d/dx.
   * In both product methods, 1st derivatives such as d/dx are
   * approximated by centered 2-point finite-difference stencils, and 
   * the 2nd derivatives d^2/dx^2 are approximated by the usual 3-point
   * stencil.
   */
  public enum Method {
    SYMMETRIC,
    PRODUCT1,
    PRODUCT2
  }

  /**
   * Source function.
   */
  public interface Source {

    /**
     * Adds the source function for specified time t to the solution f.
     * @param t the current time.
     * @param f array containing solution.
     */
    public void add(float t, float[] f);
  }

  /**
   * Constructs a demo for the specified method and modeling parameters.
   * Uses a default source function at the center of model.
   * @param method the method used to approximate spatial derivatives.
   * @param dt the time sampling interval.
   * @param dx the spatial sampling interval.
   * @param d array of densities.
   * @param v array of velocities.
   */
  public Wave1(Method method, float dt, float dx, float[] d, float[] v) {
    this(new DefaultSource(dt,dx,d.length),method,dt,dx,d,v);
  }

  /**
   * Constructs a demo for the specified method and modeling parameters.
   * @param source the source function.
   * @param method the method used to approximate spatial derivatives.
   * @param dt the time sampling interval.
   * @param dx the spatial sampling interval.
   * @param d array of densities.
   * @param v array of velocities.
   */
  public Wave1(
    Source source, Method method, 
    float dt, float dx, float[] d, float[] v) 
  {
    _source = source;
    _method = method;
    _dt = dt;
    _dx = dx;
    _nx = d.length;
    _od = new float[_nx];
    _dvs = new float[_nx];
    for (int ix=0; ix<_nx; ++ix) {
      _od[ix] = 1.0f/d[ix];
      _dvs[ix] = d[ix]*(v[ix]*v[ix])*(dt*dt)/(dx*dx);
    }
    _fm = new float[_nx];
    _f0 = new float[_nx];
    _fp = new float[_nx];
  }

  /**
   * Gets the current time.
   */
  public float time() {
    return _t;
  }

  /**
   * Steps f(x,t) forward the specified number of time steps.
   * @param nt the number of times steps.
   * @return array containing the current approximation to f(x).
   */
  public float[] step(int nt) {
    for (int it=0; it<nt; ++it)
      step();
    return _f0;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Source _source; // source function
  private Method _method; // method used for spatial derivatives
  private float _t; // current time
  private float _dt; // time sampling interval
  private float _dx; // spatial sampling interval
  private int _nx; // number of samples in spatial dimension
  private float[] _od; // array[n] of 1/density
  private float[] _dvs; // array[n] of density * velocity * velocity
  private float[] _fm,_f0,_fp; // f(x;t-1), f(x;t), and f(x;t+1)

  private static class DefaultSource implements Source {
    public DefaultSource(float dt, float dx, int nx) {
      _dt = dt;
      _dx = dx;
      _nx = nx;
      _sigmat = 5.0f*dt;
      _sigmax = 5.0f*dx;
    }
    public void add(float t, float[] f) {
      float t0 = 5.0f*_sigmat;
      float tn = (t-t0)/_sigmat;
      float st = -tn*exp(-0.5f*tn*tn);
      float x0 = 0.5f*_dx*(_nx-1);
      for (int ix=0; ix<_nx; ++ix) {
        float xi = ix*_dx;
        float xn = (xi-x0)/_sigmax;
        float sx = exp(-0.5f*xn*xn);
        f[ix] += st*sx;
      }
    }
    private int _nx;
    private float _dt,_dx,_sigmat,_sigmax;
  }

  /**
   * Update solution for time step.
   */
  private float[] step() {

    // 2nd spatial derivative (including derivative of 1/density).
    if (_method==Method.SYMMETRIC) {
      dx2Symmetric(_od,_f0,_fp);
    } else if (_method==Method.PRODUCT1) {
      dx2Product1(_od,_f0,_fp);
    } else if (_method==Method.PRODUCT2) {
      dx2Product2(_od,_f0,_fp);
    }
    
    // Integrate the time-derivatives.
    for (int ix=0; ix<_nx; ++ix)
      _fp[ix] = 2.0f*_f0[ix]-_fm[ix]+_dvs[ix]*_fp[ix];

    // Accumulate source function s(x,t) for current time t.
    _source.add(_t,_fp);

    // Time has increased by the time sampling interval.
    float[] ft = _fm;
    _fm = _f0;
    _f0 = _fp;
    _fp = ft;
    _t += _dt;
    return _f0;
  }

  // 2nd spatial derivative: d/dx 1/density d/dx.
  private static void dx2Symmetric(float[] od, float[] f, float[] g) {
    int nx = od.length;
    float odi,df,dg;
    g[0] = 0.0f;
    for (int ix=1; ix<nx; ++ix) {
      //odi = 0.5f*(od[ix]+od[ix-1]);
      odi = od[ix];
      df  = f[ix  ];
      df -= f[ix-1];
      dg = -odi*df;
      g[ix-1] -= dg;
      g[ix  ]  = dg;
    }
  }

  // 2nd spatial derivative: a common finite-difference method
  private static void dx2Product1(float[] od, float[] f, float[] g) {
    int nx = od.length;
    float odi,dod,df1,df2;
    g[0] = g[nx-1] = 0.0f; // TODO: zero-slope boundary conditions
    for (int ix=1; ix<nx-1; ++ix) {
      odi = od[ix];
      dod = 0.5f*(od[ix+1]-od[ix-1]);
      df1 = 0.5f*( f[ix+1]- f[ix-1]);
      df2 = f[ix+1]-2.0f*f[ix]+f[ix-1];
      g[ix] = odi*df2+dod*df1;
    }
  }

  // 2nd spatial derivative: a second common finite-difference method
  private static void dx2Product2(float[] od, float[] f, float[] g) {
    int nx = od.length;
    float odi,ddi,df1,df2;
    g[0] = g[nx-1] = 0.0f; // TODO: zero-slope boundary conditions
    for (int ix=1; ix<nx-1; ++ix) {
      odi = od[ix];
      ddi = 0.5f*(1.0f/od[ix+1]-1.0f/od[ix-1]);
      df1 = 0.5f*( f[ix+1]- f[ix-1]);
      df2 = f[ix+1]-2.0f*f[ix]+f[ix-1];
      g[ix] = odi*df2-odi*odi*ddi*df1;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void plotWaves(float[][] w) {
    Color[] colors = {Color.RED,Color.GREEN,Color.BLUE};
    PlotPanel pp = new PlotPanel(1,1,
      PlotPanel.Orientation.X1RIGHT_X2UP,
      PlotPanel.AxesPlacement.LEFT_BOTTOM);
    pp.setVLimits(-30.0,30.0);
    int nx = w[0].length;
    double dx = 1.0;
    double fx = -0.5f*(nx-1)*dx;
    Sampling sx = new Sampling(nx,dx,fx);
    for (int iw=0; iw<w.length; ++iw) {
      PointsView pv = pp.addPoints(sx,w[iw]);
      pv.setLineColor(colors[iw%3]);
      pv.setLineWidth(3);
    }
    PlotFrame pf = new PlotFrame(pp);
    pf.setSize(1000,800);
    pf.setVisible(true);
  }

  private static void test1() {
    int nx = 801;
    float[] d = new float[nx];
    float[] v = new float[nx];
    float dl = 1.0f; // density on left
    float dr = 0.2f; // density on right
    float rc = (dr-dl)/(dr+dl);
    float tc = 1.0f+rc;
    System.out.println("rc = "+rc+" tc = "+tc);
    for (int ix=0; ix<nx; ++ix) {
      v[ix] = 1.0f;
      d[ix] = (ix<3*nx/4)?dl:dr;
    }
    float dt = 0.5f;
    float dx = 1.0f;

    Wave1.Method[] methods = {
      Wave1.Method.SYMMETRIC,
      Wave1.Method.PRODUCT1,
      Wave1.Method.PRODUCT2
    };
    int nmethod = methods.length;
    int nt = 5;
    int mt = 125;
    float[][][] f = new float[nt][nmethod][nx];
  
    for (int imethod=0; imethod<nmethod; ++imethod) {
      Wave1 wave1 = new Wave1(methods[imethod],dt,dx,d,v);
      for (int it=0; it<nt; ++it) {
        float[] fi = wave1.step(mt);
        float[] ft = ArrayMath.copy(nx/2,fi);
        float fmax = ArrayMath.max(ft);
        System.out.println("fmax="+fmax+" fr="+fmax*rc+" ft="+fmax*tc);
        ArrayMath.copy(fi,f[it][imethod]);
      }
    }
    for (int it=0; it<nt; ++it) {
      plotWaves(f[it]);
    }
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        test1();
      }
    });
  }
}
