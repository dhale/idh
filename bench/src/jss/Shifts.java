package jss;

import edu.mines.jtk.util.*;

/**
 * Estimates apparent and physical shifts from measured shifts.
 * @author Dave Hale
 * @version 2009.03.15
 */
public class Shifts {

  public static void shifts(
    float r, float dx, float dz,
    float[][] ux, float[][] uz,
    float[][] vx, float[][] vz,
    float[][] px, float[][] pz)
  {
    int nx = ux.length;
    int nz = ux[0].length;
    float[][] phi = Array.zerofloat(nz,nx);
    for (int ix=0; ix<nx; ++ix) {
      for (int iz=1; iz<nz; ++iz) {
        phi[ix][iz] = phi[ix][iz-1]+dz*uz[ix][iz];
      }
    }
    float rorp1 = r/(r+1.0f);
    for (int ix=1; ix<nx-1; ++ix) {
      for (int iz=0; iz<nz; ++iz) {
        vx[ix][iz] = -rorp1*(phi[ix+1][iz]-phi[ix-1][iz])/(2.0f*dx);
        vz[ix][iz] = rorp1*uz[ix][iz];
        px[ix][iz] = ux[ix][iz]-vx[ix][iz];
        pz[ix][iz] = uz[ix][iz]-vz[ix][iz];
      }
    }
  }
}
