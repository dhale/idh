package sw;

import edu.mines.jtk.dsp.Sampling;

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
   * Estimates partial derivative of delta-x with respect to time t.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of partial derivatives d(delta-x)/dt.
   */
  public float[][][] ddxdt(float r, float[] v0, float[][][] deltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dx = _sx.getDelta();
    double dt = _st.getDelta();
    float[][][] ddxdt = new float[nx][ny][nt];
    float scale = (float)(0.125*r/(1.0+r)/dx);
    for (int ix=1; ix<nx-1; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int it=0; it<nt; ++it) {
          float ddt = deltat[ix+1][iy][it]-deltat[ix-1][iy][it];
          float v0s = v0[it]*v0[it];
          ddxdt[ix][iy][it] = -scale*v0s*ddt;
        }
      }
    }
    return ddxdt;
  }

  /**
   * Estimates partial derivative of delta-y with respect to time t.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of partial derivatives d(delta-y)/dt.
   */
  public float[][][] ddydt(float r, float[] v0, float[][][] deltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dy = _sy.getDelta();
    double dt = _st.getDelta();
    float[][][] ddydt = new float[nx][ny][nt];
    float scale = (float)(0.125*r/(1.0+r)/dy);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=1; iy<ny-1; ++iy) {
        for (int it=0; it<nt; ++it) {
          float ddt = deltat[ix][iy+1][it]-deltat[ix][iy-1][it];
          float v0s = v0[it]*v0[it];
          ddydt[ix][iy][it] = -scale*v0s*ddt;
        }
      }
    }
    return ddydt;
  }

  /**
   * Computes partial derivative of delta with respect to time t.
   * @param delta array of shifts delta.
   * @return array of partial derivatives d(delta)/dt.
   */
  public float[][][] dddt(float[][][] delta) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dt = _st.getDelta();
    float[][][] dddt = new float[nx][ny][nt];
    float scale = (float)(0.5/dt);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int it=1; it<nt-1; ++it) {
          float dd = delta[ix][iy][it+1]-delta[ix][iy][it-1];
          dddt[ix][iy][it] = scale*dd;
        }
      }
    }
    return dddt;
  }

  /**
   * Estimates horizontal shifts delta-x.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of horizontal shifts delta-x.
   */
  public float[][][] getDeltaX(float r, float[] v0, float[][][] deltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dt = _st.getDelta();
    float[][][] ddxdt = ddxdt(r,v0,deltat);
    float[][][] deltax = ddxdt;
    float scale = (float)(dt);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float d = deltax[ix][iy][0] = 0.0f;
        for (int it=1; it<nt; ++it) {
          deltax[ix][iy][it] = d = d+scale*ddxdt[ix][iy][it];
        }
      }
    }
    return deltax;
  }

  /**
   * Estimates horizontal shifts delta-y.
   * @param r Hatchell and Bournes' parameter R.
   * @param v0 velocity as a function of vertical time t.
   * @param deltat array of vertical time shifts.
   * @return array of horizontal shifts delta-y.
   */
  public float[][][] getDeltaY(float r, float[] v0, float[][][] deltat) {
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nt = _st.getCount();
    double dt = _st.getDelta();
    float[][][] ddydt = ddydt(r,v0,deltat);
    float[][][] deltay = ddydt;
    float scale = (float)(dt);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float d = deltay[ix][iy][0] = 0.0f;
        for (int it=1; it<nt; ++it) {
          deltay[ix][iy][it] = d = d+scale*ddydt[ix][iy][it];
        }
      }
    }
    return deltay;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  Sampling _sx,_sy,_st;
}
