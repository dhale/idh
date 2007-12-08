package sw;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Estimates horizontal components of shifts from vertical time shifts.
 * The estimates are based on Hatchell and Bournes' simple approximations 
 * that relate small vertical time shifts to small perturbations in velocity. 
 * The velocity perturbations vary laterally with the time shifts, and this 
 * lateral velocity variation causes horizontal shifts that we can estimate 
 * from the vertical time shifts.
 * @author Dave Hale
 * @version 2007.11.22
 */
public class ShiftsXY {

  /**
   * Constructs shifts for specified samplings.
   * @param sx sampling of horizontal x dimension (3rd dimension).
   * @param sy sampling of horizontal y dimension (2nd dimension).
   * @param st sampling of vertical time dimension (1st dimension).
   */
  public ShiftsXY(Sampling sx, Sampling sy, Sampling st) {
    _sx = sx;
    _sy = sy;
    _st = st;
  }

  /**
   * Sets the lower bound on time used in integrals over time.
   * @param tmin the minimum time; the default is zero.
   */
  public void setMinimumTime(float tmin) {
    _itmin = _st.indexOfNearest(tmin);
  }

  /**
   * Estimates horizontal shifts in x direction.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of horizontal shifts in x direction.
   */
  public float[][][] getDeltaX(float r, float[] v0, float[][][] deltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dx = _sx.getDelta();
    double dt = _st.getDelta();
    float[][][] deltax = new float[nx][ny][nt];
    float scale = (float)(0.125*r/(1.0+r)*dt/dx);
    for (int ix=1; ix<nx-1; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float d = 0.0f;
        for (int it=_itmin; it<nt; ++it) {
          float ddtdx = deltat[ix+1][iy][it]-deltat[ix-1][iy][it];
          float v0s = v0[it]*v0[it];
          deltax[ix][iy][it] = d = d-scale*v0s*ddtdx;
        }
      }
    }
    return deltax;
  }

  /**
   * Estimates horizontal shifts in y direction.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of horizontal shifts in y direction.
   */
  public float[][][] getDeltaY(float r, float[] v0, float[][][] deltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dy = _sy.getDelta();
    double dt = _st.getDelta();
    float[][][] deltay = new float[nx][ny][nt];
    float scale = (float)(0.125*r/(1.0+r)*dt/dy);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=1; iy<ny-1; ++iy) {
        float d = 0.0f;
        for (int it=_itmin; it<nt; ++it) {
          float ddtdy = deltat[ix][iy+1][it]-deltat[ix][iy-1][it];
          float v0s = v0[it]*v0[it];
          deltay[ix][iy][it] = d = d-scale*v0s*ddtdy;
        }
      }
    }
    return deltay;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  Sampling _sx,_sy,_st;
  int _itmin = 0;

  /**
   * Partial derivative g = df/dx.
   */
  private float[][][] ddx(float[][][] f) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dx = _sx.getDelta();
    float[][][] g = new float[nx][ny][nt];
    float scale = (float)(0.5/dx);
    for (int ix=1; ix<nx-1; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int it=0; it<nt; ++it) {
          g[ix][iy][it] = scale*(f[ix+1][iy][it]-f[ix-1][iy][it]);
        }
      }
    }
    return g;
  }

  /**
   * Partial derivative g = df/dy.
   */
  private float[][][] ddy(float[][][] f) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dy = _sy.getDelta();
    float[][][] g = new float[nx][ny][nt];
    float scale = (float)(0.5/dy);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=1; iy<ny-1; ++iy) {
        for (int it=0; it<nt; ++it) {
          g[ix][iy][it] = scale*(f[ix][iy+1][it]-f[ix][iy-1][it]);
        }
      }
    }
    return g;
  }

  /**
   * Partial derivative g = df/dt.
   */
  private float[][][] ddt(float[][][] f) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dt = _st.getDelta();
    float[][][] g = new float[nx][ny][nt];
    float scale = (float)(0.5/dt);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int it=1; it<nt-1; ++it) {
          g[ix][iy][it] = scale*(f[ix][iy][it+1]-f[ix][iy][it-1]);
        }
      }
    }
    return g;
  }

  /**
   * Partial derivative of velocity v with respect to either x or y.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param ddeltat partial derivative of deltat with respect to 
   *  vertical time t and either distance x or y.
   * @return partial derivative of v.
   */
  private float[][][] dv(float r, float[] v0, float[][][] ddeltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    float s = r/(1.0f+r);
    float[][][] dv = new float[nx][ny][nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int it=0; it<nt; ++it) {
          dv[ix][iy][it] = -s*v0[it]*ddeltat[ix][iy][it];
        }
      }
    }
    return dv;
  }

  /**
   * Angle theta measured from vertical in either x or y direction.
   * @param v0 velocity as a function of vertical time t.
   * @param dv partial derivative of velocity with respect to x or y.
   * @return angle theta
   */
  private float[][][] theta(float[] v0, float[][][] dv) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    float dt = (float)_st.getDelta();

    // Precompute dv0/dt * 1/v0 (same for all x and y).
    float[] dv0 = new float[nt];
    float hodt = (float)(0.5/dt);
    for (int it=1; it<nt-1; ++it)
      dv0[it] = hodt*(v0[it+1]-v0[it-1])/v0[it];

    // Compute theta by simple integration of non-linear ODE.
    float[][][] theta = new float[nx][ny][nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float t = theta[ix][iy][_itmin] = 0.0f;
        for (int it=_itmin+1; it<nt; ++it) {
          float s = dv0[it]*sin(t);
          float c = -0.5f*dv[ix][iy][it]*cos(t);
          theta[ix][iy][it] = t = t+(s+c)*dt;
        }
      }
    }
    return theta;
  }

  private float[][][] delta(float[] v0, float[][][] theta) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    float dt = (float)_st.getDelta();
    float[][][] delta = new float[nx][ny][nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float d = delta[ix][iy][_itmin] = 0.0f;
        for (int it=_itmin+1; it<nt; ++it) {
          delta[ix][iy][it] = d = d-0.5f*v0[it]*sin(theta[ix][iy][it])*dt;
        }
      }
    }
    return delta;
  }

  /**
   * Estimates horizontal shifts in x direction.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of horizontal shifts in x direction.
   */
  private float[][][] getDeltaXOld(float r, float[] v0, float[][][] deltat) {
    float[][][] ddt = ddt(deltat);
    deltat = null;
    float[][][] ddtdx = ddx(ddt);
    ddt = null;
    float[][][] dvdx = dv(r,v0,ddtdx);
    ddtdx = null;
    float[][][] theta = theta(v0,dvdx);
    dvdx = null;
    float[][][] delta = delta(v0,theta);
    theta = null;
    return delta;
  }

  /**
   * Estimates horizontal shifts in y direction.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of horizontal shifts in y direction.
   */
  private float[][][] getDeltaYOld(float r, float[] v0, float[][][] deltat) {
    float[][][] ddt = ddt(deltat);
    deltat = null;
    float[][][] ddtdy = ddy(ddt);
    ddt = null;
    float[][][] dvdy = dv(r,v0,ddtdy);
    ddtdy = null;
    float[][][] theta = theta(v0,dvdy);
    dvdy = null;
    float[][][] delta = delta(v0,theta);
    theta = null;
    return delta;
  }
}
