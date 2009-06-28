package jss;

import edu.mines.jtk.util.*;

/**
 * Computations for apparent shifts. 
 * @author Dave Hale
 * @version 2009.03.15
 */
public class Shifts {

  public static void shifts(
    float r, float dx, float dz,
    float[][] uxa, float[][] uza,
    float[][] uxv, float[][] uzv,
    float[][] ux, float[][] uz)
  {
    int nx = ux.length;
    int nz = ux[0].length;
    float[][] phi = ArrayMath.zerofloat(nz,nx);
    for (int ix=0; ix<nx; ++ix) {
      for (int iz=1; iz<nz; ++iz) {
        phi[ix][iz] = phi[ix][iz-1]+dz*uza[ix][iz];
      }
    }
    float zscale = r/(r+1.0f);
    float xscale = -zscale*0.5f/dx;
    float[] xasum = new float[nz];
    float[] xvsum = new float[nz];
    for (int ix=1; ix<nx-1; ++ix) {
      for (int iz=0; iz<nz; ++iz) {
        uxv[ix][iz] = xscale*(phi[ix+1][iz]-phi[ix-1][iz]);
        uzv[ix][iz] = zscale*uza[ix][iz];
        ux[ix][iz] = uxa[ix][iz]-uxv[ix][iz];
        uz[ix][iz] = uza[ix][iz]-uzv[ix][iz];
        xasum[iz] += uxa[ix][iz];
        xvsum[iz] += uxv[ix][iz];
      }
    }
    //SimplePlot sp = SimplePlot.asPoints(xasum);
    //sp.addPoints(xvsum);
  }
}
