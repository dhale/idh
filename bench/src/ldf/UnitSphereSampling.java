/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import static java.lang.Math.*;
import edu.mines.jtk.util.Check;

/**
 * Quasi-uniform sampling of a unit-sphere.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.02
 */
public class UnitSphereSampling {

  /**
   * Constructs a sampling for the specified number of bits.
   * Sample indices are signed integers with no more than this
   * number of bits, which includes the sign bit.
   * @param nbits the number of bits.
   */
  public UnitSphereSampling(int nbits) {
    initialize(nbits);
  }

  /**
   * Gets the sampled point for the specified index.
   * For efficiency, returns the array {x,y,z} of point coordinates 
   * by reference, not by copy. These coordinates must not be modified.
   * @param index the index of the sampled point.
   * @return array {x,y,z} of point coordinates; by reference, not by copy.
   */
  public float[] getPoint(int index) {
    float[][] p = _pu;
    if (index<0) {
      index = -index;
      p = _pl;
    }
    return p[index];
  }

  /**
   * Gets the index of the sampled point nearest to the specified point.
   * Here, the nearest sampled point is that nearest on the octahedron.
   * @param x the x-coordinate of the point.
   * @param y the y-coordinate of the point.
   * @param z the z-coordinate of the point.
   * @return the sample index.
   */
  public int getIndex(float x, float y, float z) {
    double ax = (x>=0.0f)?x:-x;
    double ay = (y>=0.0f)?y:-y;
    double az = (z>=0.0f)?z:-z;
    double scale = 1.0/(ax+ay+az);
    double r = x*scale;
    double s = y*scale;
    int ir = (int)(0.5+(r+1.0)/_d);
    int is = (int)(0.5+(s+1.0)/_d);
    int index = _ip[is][ir];
    return (z>=0.0f)?index:-index;
  }

  /**
   * Gets the index of the sampled point nearest to the specified point.
   * Here, the nearest sampled point is that nearest on the octahedron.
   * @param xyz the array {x,y,z} of point coordinates.
   * @return the sample index.
   */
  public int getIndex(float[] xyz) {
    return getIndex(xyz[0],xyz[1],xyz[2]);
  }

  /**
   * Gets the maximum sample index, a positive integer. The smallest
   * positive index is one. The smallest index is the negative of the 
   * maximum index, and the largest negative index is minus one.
   * @return the maximum index.
   */
  public int getMaxIndex() {
    return _nindex;
  }

  /**
   * Gets the number of points sampled on this unit sphere.
   * @return the number of points sampled.
   */
  public int countPointsSampled() {
    return _npoint;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // The unit sphere is projected onto an octahedron with corners that lie
  // on the sphere; that is, on the x, y, and z axes. The upper hemisphere
  // (z>=0) is projected onto the upper half of the octahedron. 
  //
  // Imagine that this upper half is flattened into the x and y plane, and
  // let r and s the projected coordinates in this flat plane. The upper 
  // hemisphere corresponds to a diamond in the r-s plane. The number of
  // points in this diamond is the number of points sampled for the upper
  // hemisphere, including the points on the equator for which z=0.
  //
  // The lower hemisphere (z<=0) is sampled in the same way, and its 
  // mapping also includes points on the equator. The number of samples
  // in these diamonds (the sampling resolution) is limited by the number 
  // of bits (nbits) in the signed integer sample indices.
  //
  // The r-s plane is sampled on an n by n grid, where n = 2*m+1 and m is 
  // the number of sampling intervals along each of the positive r and s 
  // axes. But not all samples in this grid are used.
  // 
  // In the example below, for nbits = 6, m = 3, and n = 7, only the 25
  // samples marked with X correspond to sampled points. Samples marked 
  // with 0 are unused, but are included in the grid for simplicity. The
  // indices for the upper hemisphere are in the range [1,25]. For the 
  // lower hemisphere, indices are in the range [-25,-1]. The index 0 is
  // unused. This range of indices [-25,25] fits in a signed 6-bit integer 
  // with range[-32:31].
  //
  //                  s
  //                  ^
  //                  |
  //     3  |0  0  0  X  0  0  0
  //     2  |0  0  X  X  X  0  0
  //     1  |0  X  X  X  X  X  0
  //     0  |X  X  X  X  X  X  X  ---> r
  //    -1  |0  X  X  X  X  X  0
  //    -2  |0  0  X  X  X  0  0
  //    -3  |0  0  0  X  0  0  0
  //         -------------------
  //        -3 -2 -1  0  1  2  3
  //
  // In a more practical example, nbits = 16, m = 127, and n = 255, with
  // sample indices in [-32513,-1] and [1,32513]. In this example, the
  // number of unique points sampled is 64008, which is less than the
  // maximum possible 65536 points that could be represented with 16 bits.
  
  private int _nbits; // number of bits used in quantization
  private int _m; // number of samples for positive r and s, not including zero
  private int _n; // number of samples of r and s
  private int _nindex; // number of positive/negative indices
  private int _npoint; // number of unique points
  private double _d; // sampling interval for r and s
  private double _od; // one over sampling interval = 1/d
  private float[][] _pu; // table of points in upper hemisphere (z>=0)
  private float[][] _pl; // table of points in lower hemisphere (z<=0)
  private int[][] _ip; // table[n][n] of point indices

  private void initialize(int nbits) {
    Check.argument(nbits>=3,"nbits>=3");
    Check.argument(nbits<=32,"nbits<=32");

    // Number of bits in sample indices, including the sign bit.
    _nbits = nbits;

    // Sampling of the r-s plane with an n by n grid.
    _m = 1;
    while (2*_m*(1+_m)<=(1<<(nbits-1)))
      ++_m;
    --_m;
    _n = 2*_m+1;

    // Sampling interval and its inverse for r and s.
    _d = 1.0/_m;
    _od = _m;

    // Number of positive/negative indices.
    _nindex = 1+2*_m*(1+_m);

    // Number of unique sampled points. The number of points on the equator 
    // is 4*n-2; these points for which z=0 appear in both tables for the
    // upper and lower hemispheres, for both positive and negative indices.
    // Here we do not count them twice.
    _npoint = 2*_nindex-4*_n+2;

    trace("m="+_m+" n="+_n+" nindex="+_nindex+" npoint="+_npoint);

    // Tables for points in upper and lower hemispheres.
    _pu = new float[1+_nindex][];
    _pl = new float[1+_nindex][];

    // Table of point indices.
    _ip = new int[_n][_n];

    // For all sampled s on flattened octahedron, ...
    for (int is=0,js=-_m,index=0; is<_n; ++is,++js) {

      // Planar coordinate s and |s|.
      double s = js*_d;
      double as = (s>=0.0)?s:-s;

      // For all sampled r on flattened octahedron, ...
      for (int ir=0,jr=-_m; ir<_n; ++ir,++jr) {

        // Process only samples the octahedral diamond corresponding
        // to the upper and lower hemispheres. Other points in the
        // table will be null.
        if (abs(jr)+abs(js)<=_m) {

          // Increment and store index in table.
          _ip[is][ir] = ++index;
          //trace("ir="+ir+" is="+is+" index="+index);

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
          float[] pu = _pu[index] = new float[3];
          float[] pl = _pl[index] = new float[3];
          pu[0] = x;  pu[1] = y;  pu[2] =  z;
          pl[0] = x;  pl[1] = y;  pl[2] = -z;
        }
      }
    }
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static java.util.Random _random = new java.util.Random();
  private static float[] randomPoint() {
    float x = -1.0f+2.0f*_random.nextFloat();
    float y = -1.0f+2.0f*_random.nextFloat();
    float z = -1.0f+2.0f*_random.nextFloat();
    float s = 1.0f/(float)sqrt(x*x+y*y+z*z);
    return new float[]{x*s,y*s,z*s};
  }

  private static void testByteIndex() {
    UnitSphereSampling uss = new UnitSphereSampling(16);
    int npoint = 10;
    for (int ipoint=0; ipoint<npoint; ++ipoint) {
      float[] p = randomPoint();
      int i = uss.getIndex(p);
      float[] q = uss.getPoint(i);
      trace("ipoint="+ipoint+" i="+i);
      edu.mines.jtk.util.Array.dump(p);
      edu.mines.jtk.util.Array.dump(q);
    }
  }
  public static void main(String[] args) {
    testByteIndex();
  }
}
