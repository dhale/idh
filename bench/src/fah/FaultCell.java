/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

import static fah.FaultGeometry.*;

/**
 * A fault cell is an oriented point located on a fault. Fault cells can
 * be linked to form fault skins, which may be used to analyze faults.
 * <p>
 * Fault cells are computed from images of fault likelihoods, strikes and
 * dips. Each fault cell is an oriented point located on a ridge in an image
 * of fault likelihood. Fault cells have indices (i1,i2,i3) that indicate
 * which image sample is nearest to this ridge. In this way an image sample
 * is associated with either no cell or one cell.
 * <p>
 * A fault cell has up to four neighbors ("nabors") that lie above, below,
 * left and right of the cell when viewed from above the fault, that is, when
 * looking from the hanging wall toward the footwall. Links to nabors enables
 * cells to form a skin of connected cells, which represents a fault.
 * <p>
 * Links to left and right cell nabors can be used to iterate over all cells
 * along a fault trace, a path of constant depth that is everywhere tangent to
 * fault strike. Likewise, links to cell nabors above and below a cell can be
 * used to iterate up or down a fault. However, this simple up or down
 * iteration typically does not coincide with a fault curve that is everywhere
 * tangent to fault dip.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */

public class FaultCell {

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, the specified cells are represented by quads with specified
   * size, and colors corresponding to fault likelihoods or shifts.
   * @param size the size (in samples) of the quads.
   * @param smax if zero, colors correspond to fault likelihoods with a
   * maximum value of 1. If positive, colors correspond to minus-plus shifts
   * with the specified maximum shift. If negative, colors correspond to
   * plus-minus shifts with minimum equal to the specified shift.
   * @param cmap the colormap used to compute rgb colors from floats.
   * @param cells the cells for which to compute quads.
   */
  public static float[][] getXyzUvwRgb(
      float size, float smax, ColorMap cmap, FaultCell[] cells) {
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList fcl = new FloatList();
    size *= 0.5f;
    float[] qa = {0.0f,-size,-size};
    float[] qb = {0.0f, size,-size};
    float[] qc = {0.0f, size, size};
    float[] qd = {0.0f,-size, size};
    for (FaultCell cell:cells) {
      float x1 = cell.x1;
      float x2 = cell.x2;
      float x3 = cell.x3;
      float w1 = cell.w1;
      float w2 = cell.w2;
      float w3 = cell.w3;
      float fp = toRadians(cell.fp);
      float ft = toRadians(cell.ft);
      float cp = cos(fp);
      float sp = sin(fp);
      float ct = cos(ft);
      float st = sin(ft);
      float[] ra = rotatePoint(cp,sp,ct,st,qa);
      float[] rb = rotatePoint(cp,sp,ct,st,qb);
      float[] rc = rotatePoint(cp,sp,ct,st,qc);
      float[] rd = rotatePoint(cp,sp,ct,st,qd);
      float a1 = x1+ra[0], a2 = x2+ra[1], a3 = x3+ra[2];
      float b1 = x1+rb[0], b2 = x2+rb[1], b3 = x3+rb[2];
      float c1 = x1+rc[0], c2 = x2+rc[1], c3 = x3+rc[2];
      float d1 = x1+rd[0], d2 = x2+rd[1], d3 = x3+rd[2];
      xyz.add(a3); xyz.add(a2); xyz.add(a1);
      xyz.add(b3); xyz.add(b2); xyz.add(b1);
      xyz.add(c3); xyz.add(c2); xyz.add(c1);
      xyz.add(d3); xyz.add(d2); xyz.add(d1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      float fc = cell.fl;
      if (smax<0.0f) {
        fc = cell.spm;
      } else if (smax>0.0f) {
        fc = cell.smp;
      }
      fcl.add(fc);
      fcl.add(fc);
      fcl.add(fc);
      fcl.add(fc);
    }
    float[] fc = fcl.trim();
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  /**
   * Returns an array of packed (x,y,z) coordinates for a fault curve.
   * The fault curve is everywhere tangent to fault dip, and contains the
   * point for this cell. Returned coordinates are in above-to-below order.
   * @return array of packed (x,y,z) coordinates.
   */
  public float[] getFaultCurveXyz() {
    FloatList xyz = new FloatList();
    float[] p = new float[3];

    // Find xyz for this cell and above this cell by walking up dip.
    FaultCell cell = this;
    p[0] = x1; p[1] = x2; p[2] = x3;
    for (int j1=cell.i1; j1==cell.i1; --j1) {
      xyz.add(p[2]); xyz.add(p[1]); xyz.add(p[0]);
      cell = cell.walkUpDipFrom(p);
    }

    // Number of xyz found above this cell.
    int na = xyz.n/3-1;

    // Find xyz for cells below this one by walking down dip.
    cell = this;
    p[0] = x1; p[1] = x2; p[2] = x3;
    cell = cell.walkDownDipFrom(p); // in this walk, skip this cell
    for (int j1=cell.i1; j1==cell.i1; ++j1) {
      xyz.add(p[2]); xyz.add(p[1]); xyz.add(p[0]);
      cell = cell.walkDownDipFrom(p);
    }

    // Flip the order of all xyz found above this cell while walking up.
    float[] xyzs = xyz.trim();
    for (int ia=0,i=ia*3,ja=na-1,j=ja*3; ia<ja; ++ia,i+=3,--ja,j-=3) {
      float xi = xyzs[i  ];
      float yi = xyzs[i+1];
      float zi = xyzs[i+2];
      xyzs[i  ] = xyzs[j  ];
      xyzs[i+1] = xyzs[j+1];
      xyzs[i+2] = xyzs[j+2];
      xyzs[j  ] = xi;
      xyzs[j+1] = yi;
      xyzs[j+2] = zi;
    }
    return xyzs;
  }

  /**
   * Returns an array of packed (x,y,z) coordinates for a fault trace.
   * The fault trace is everywhere tangent to fault strike, and contains 
   * the point for this cell. Returned coordinates are left-to-right order.
   * @return array of packed (x,y,z) coordinates.
   */
  public float[] getFaultTraceXyz() {
    FloatList xyz = new FloatList();

    // First add coordinates for this cell.
    xyz.add(x3); xyz.add(x2); xyz.add(x1);

    // Add coordinates for cells to the left of this cell. Take care to handle
    // the case in which this cell is found while walking left.
    FaultCell c;
    for (c=this.cl; c!=null && c!=this; c=c.cl) {
      xyz.add(c.x3); xyz.add(c.x2); xyz.add(c.x1);
    }

    // Number of xyz found left of this cell.
    int nl = xyz.n/3-1;

    // If we did not end at this cell, then add coordinates of cells to
    // the right of this cell. Again, watch for this cell.
    if (c!=this) { // if we did not end at this cell, ...
      for (c=this.cr; c!=null && c!=this; c=c.cr) {
        xyz.add(c.x3); xyz.add(c.x2); xyz.add(c.x1);
      }
    }

    // Flip the order of all xyz found left of this cell.
    float[] xyzs = xyz.trim();
    for (int il=0,i=il*3,jl=nl-1,j=jl*3; il<jl; ++il,i+=3,--jl,j-=3) {
      float xi = xyzs[i  ];
      float yi = xyzs[i+1];
      float zi = xyzs[i+2];
      xyzs[i  ] = xyzs[j  ];
      xyzs[i+1] = xyzs[j+1];
      xyzs[i+2] = xyzs[j+2];
      xyzs[j  ] = xi;
      xyzs[j+1] = yi;
      xyzs[j+2] = zi;
    }
    return xyzs;
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  int i1,i2,i3; // cell indices
  float x1,x2,x3; // cell coordinates
  float fl,fp,ft; // likelihood, strike (phi) and dip (theta)
  float u1,u2,u3,us; // dip vector and scale factor 1/u1 = 1/sin(theta)
  float v1,v2,v3; // strike vector
  float w1,w2,w3; // normal vector
  FaultCell ca,cb,cl,cr; // nabors above, below, left and right
  FaultSkin skin; // if not null, the skin to which this cell belongs
  int i2m,i2p; // sample indices i2 for minus and plus sides of cell
  int i3m,i3p; // sample indices i3 for minus and plus sides of cell
  float[] emp; // array of alignment errors minus-plus
  float[] epm; // array of alignment errors plus-minus
  float smp; // shift from minus side to plus side of cell
  float spm; // shift from plus side to minus side of cell
  float s1,s2,s3; // fault dip-slip vector

  FaultCell(float x1, float x2, float x3, float fl, float fp, float ft) {
    this.x1 = x1; 
    this.x2 = x2; 
    this.x3 = x3;
    this.fl = fl; 
    this.fp = fp; 
    this.ft = ft;
    i1 = round(x1);
    i2 = round(x2);
    i3 = round(x3);
    float[] u = faultDipVectorFromStrikeAndDip(fp,ft);
    float[] v = faultStrikeVectorFromStrikeAndDip(fp,ft);
    float[] w = faultNormalVectorFromStrikeAndDip(fp,ft);
    u1 = u[0]; u2 = u[1]; u3 = u[2]; us = 1.0f/u1;
    v1 = v[0]; v2 = v[1]; v3 = v[2];
    w1 = w[0]; w2 = w[1]; w3 = w[2];

    // Indices (i2m,i2p) and (i3m,i3p) for minus-plus pairs of samples.
    // Cell normal vector w points from the minus side to the plus side.
    i2m = i2p = i2;
    i3m = i3p = i3;
    if (x2>i2) {
      ++i2p;
    } else if (x2<i2) {
      --i2m;
    }
    if (x3>i3) {
      ++i3p;
    } else if (x3<i3) {
      --i3m;
    }
    if ((i2p-i2m)*w2<0.0f) {
      int i2t = i2m; 
      i2m = i2p; 
      i2p = i2t;
    }
    if ((i3p-i3m)*w3<0.0f) {
      int i3t = i3m; 
      i3m = i3p; 
      i3p = i3t;
    }
  }

  /**
   * Returns the distance squared from this cell to a point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return distance squared.
   */
  float distanceTo(float p1, float p2, float p3) {
    return sqrt(distanceSquaredTo(p1,p2,p3));
  }

  /**
   * Returns the distance squared from this cell to a point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return distance squared.
   */
  float distanceSquaredTo(float p1, float p2, float p3) {
    float d1 = p1-x1;
    float d2 = p2-x2;
    float d3 = p3-x3;
    return d1*d1+d2*d2+d3*d3;
  }

  /**
   * Returns the signed distance from the plane of this cell to a point. The
   * distance is positive if the point lies on the side of the plane toward
   * which the normal vector points, negative if on the other side.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return signed distance.
   */
  float distanceFromPlaneTo(float p1, float p2, float p3) {
    return w1*(p1-x1)+w2*(p2-x2)+w3*(p3-x3);
  }

  /**
   * Gets the cell above this cell that is nearest to the specified point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return the cell above; null, if none.
   */
  FaultCell getCellAboveNearestTo(float p1, float p2, float p3) {
    FaultCell cla = (cl!=null)?cl.ca:null;
    FaultCell cra = (cr!=null)?cr.ca:null;
    return nearestCell(ca,cla,cra,p1,p2,p3);
  }

  /**
   * Gets the cell below this cell that is nearest to the specified point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return the cell below; null, if none.
   */
  FaultCell getCellBelowNearestTo(float p1, float p2, float p3) {
    FaultCell clb = (cl!=null)?cl.cb:null;
    FaultCell crb = (cr!=null)?cr.cb:null;
    return nearestCell(cb,clb,crb,p1,p2,p3);
  }

  /**
   * Walks a point up the fault along a curve tangent to fault dip.
   * The input point is assumed to lie in the plane of this cell,
   * and the output point will lie in the plane of the returned cell.
   * That returned cell will be one immediately above this cell, if
   * sufficient nabors exist, or this cell, otherwise.
   * @param p input and output array {p1,p2,p3} of point coordinates.
   * @return the cell with a plane that contains the output point.
   */
  FaultCell walkUpDipFrom(float[] p) {
    FaultCell cell = this;
    float p1 = p[0];
    float p2 = p[1];
    float p3 = p[2];
    assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;

    // Use dip vector of specified cell to walk the point up its plane.
    p1 -= 1.0f;
    p2 -= cell.us*cell.u2;
    p3 -= cell.us*cell.u3;

    // If a cell below is found, project point horizontally onto its plane.
    FaultCell ca = cell.getCellAboveNearestTo(p1,p2,p3);
    if (ca!=null) {
      cell = ca;
      float us = cell.us;
      float ws = us*us*cell.distanceFromPlaneTo(p1,p2,p3);
      p2 -= ws*cell.w2;
      p3 -= ws*cell.w3;
    }

    // Return updated point and cell.
    assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;
    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    return cell;
  }

  /**
   * Walks a point down the fault along a curve tangent to fault dip.
   * The input point is assumed to lie in the plane of this cell,
   * and the output point will lie in the plane of the returned cell.
   * That returned cell will be one immediately below this cell, if
   * sufficient cell nabors exist, or this cell, otherwise.
   * @param p input and output array {p1,p2,p3} of point coordinates.
   * @return the cell with a plane that contains the output point.
   */
  FaultCell walkDownDipFrom(float[] p) {
    FaultCell cell = this;
    float p1 = p[0];
    float p2 = p[1];
    float p3 = p[2];
    assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;

    // Use dip vector of specified cell to walk the point down its plane.
    p1 += 1.0f;
    p2 += cell.us*cell.u2;
    p3 += cell.us*cell.u3;

    // If a cell below is found, project point horizontally onto its plane.
    FaultCell cb = cell.getCellBelowNearestTo(p1,p2,p3);
    if (cb!=null) {
      cell = cb;
      float us = cell.us;
      float ws = us*us*cell.distanceFromPlaneTo(p1,p2,p3);
      p2 -= ws*cell.w2;
      p3 -= ws*cell.w3;
    }

    // Return updated point and cell.
    assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;
    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    return cell;
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  private static FaultCell nearestCell(
      FaultCell c1, FaultCell c2, FaultCell c3, 
      float p1, float p2, float p3) 
  {
    float ds1 = distanceSquared(c1,p1,p2,p3);
    float ds2 = distanceSquared(c2,p1,p2,p3);
    float ds3 = distanceSquared(c3,p1,p2,p3);
    float dsm = min(ds1,ds2,ds3);
    if (dsm==ds1) { 
      return c1;
    } else if (dsm==ds2) {
      return c2;
    } else {
      return c3;
    }
  }

  private static float distanceSquared(
      FaultCell c, float p1, float p2, float p3) {
    return c!=null ? c.distanceSquaredTo(p1,p2,p3) : Float.MAX_VALUE;
  }

  /**
   * Rotates a specified point by strike (phi) and dip (theta) angles. Uses
   * specified cosines (cp and ct) and sines (sp and st) of those angles. The
   * order of transformation is
   * (1) rotate around axis x3 by dip angle
   * (2) rotate around axis x1 by strike angle
   * Returns the coordinates of the rotated point.
   */
  private static float[] rotatePoint(
      float cp, float sp, float ct, float st, float[] x) {
    float x1 = x[0], x2 = x[1], x3 = x[2];
    float y1 =     ct*x1+   st*x2;
    float y2 = -cp*st*x1+cp*ct*x2+sp*x3;
    float y3 =  sp*st*x1-sp*ct*x2+cp*x3;
    return new float[]{y1,y2,y3};
  }

  private static class FloatList {
    public int n = 0;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
}
