/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import static java.lang.Math.*;

/**
 * Quasi-uniform sampling of a unit-sphere.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.02
 */
public class UnitSphereSampling {

  public UnitSphereSampling(int nbits) {
    _nbits = nbits;

    // Number of samples along r and s axes, 
    // such that (2*m+1)^2 = n^2 <= 2^nbits.
    _m = (int)((pow(2.0,nbits/2.0)-1.0)/2.0);
    _n = 2*_m+1;

    // Sampling interval for r and s.
    _d = 1.0/_m;

    makeTables();
  }

  /**
   * Gets the sampled point for the specified index.
   * For efficiency, returns the array {x,y,z} of point coordinates 
   * by reference, not by copy. These coordinates must not be modified.
   * @param index the index of the sampled point.
   * @return array {x,y,z} of point coordinates; by reference, not by copy.
   */
  public float[] getPoint(int index) {
    float[][][] p = _pu;
    if (index<0) {
      index = -index;
      p = _pl;
    }
    int it = index/_n;
    int is = index-it*_n;
    return p[it][is];
  }

  /**
   * Gets the index of the sampled point nearest to the specified point.
   * Here, the nearest sampled point is that corresponding to the
   */
  public int getIndex(float x, float y, float z) {
    double ax = (x>=0.0f)?x:-x;
    double ay = (y>=0.0f)?y:-y;
    double az = (z>=0.0f)?z:-z;
    double scale = 1.0/(ax+ay+az);
    double r = x*scale;
    double s = y*scale;
    int ir = (int)((r+1.0)/_d);
    int is = (int)((s+1.0)/_d);
    int index = ir+is*_n;
    return (z>=0.0f)?index:-index;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nbits; // number of bits used in quantization
  private int _m; // number of samples for positive r and s, not including zero
  private int _n; // number of samples of r and s
  private double _d; // sampling interval for r and s
  private float[][][] _pu; // table of points in upper hemisphere (z>=0)
  private float[][][] _pl; // table of points in lower hemisphere (z<0)

  private void makeTables() {

    // Tables for points in upper and lower hemispheres.
    _pu = new float[_n][_n][3];
    _pl = new float[_n][_n][3];

    // For all sampled s on octahedron, ...
    for (int is=0,js=-_m,index=0; is<_n; ++is,++js) {

      // Planar coordinate s and |s|.
      double s = js*_d;
      double as = (s>=0.0)?s:-s;

      // For all sampled r on octahedron, ...
      for (int ir=0,jr=-_m; ir<_n; ++ir,++jr,++index) {

        // Skip samples outside the octahedral diamond.
        // Tables for these indices will be null.
        if (abs(ir)+abs(is)>_m) continue;

        // Planar coordinate r and |r|.
        double r = jr*_d;
        double ar = (r>=0.0)?r:-r;

        // Third coordinate t (t>=0) on octahedron.
        double t = max(0.0,1.0f-ar-as);

        // Coordinates of point in upper hemisphere (z>=0).
        double scale = 1.0/sqrt(s*s+r*r+t*t);
        float x = (float)(r*scale);
        float y = (float)(s*scale);
        float z = (float)(t*scale);

        // Store coordinates in tables.
        float[] pu = _pu[is][ir];
        float[] pl = _pl[is][ir];
        pu[0] = x;  pu[1] = y;  pu[2] =  z;
        pl[0] = x;  pl[1] = y;  pl[2] = -z;
      }
    }
  }
}
