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
 * which image sample is nearest to the ridge. An image sample corresponds to
 * either no cell or one cell.
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
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors.
   * In these arrays, cells are represented by quads with specified size.
   * @param size the size (in samples) of the quads.
   * @param cells the cells for which to compute quads.
   */
  public static float[][] getXyzUvwRgb(float size, FaultCell[] cells) {
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
      float fl = cell.fl;
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
      xyz.add(a3); xyz.add(a2); xyz.add(a1); fcl.add(fl);
      xyz.add(b3); xyz.add(b2); xyz.add(b1); fcl.add(fl);
      xyz.add(c3); xyz.add(c2); xyz.add(c1); fcl.add(fl);
      xyz.add(d3); xyz.add(d2); xyz.add(d1); fcl.add(fl);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
    }
    float[] fc = fcl.trim();
    float fcmin = 0.0f;
    float fcmax = 1.0f;
    ColorMap cmap = new ColorMap(fcmin,fcmax,ColorMap.JET);
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  int i1,i2,i3; // cell indices
  float x1,x2,x3; // cell coordinates
  float fl,fp,ft; // likelihood, strike (phi) and dip (theta)
  float u1,u2,u3; // dip vector
  float v1,v2,v3; // strike vector
  float w1,w2,w3; // normal vector
  float ws23; // scale factor 1/sqrt(w2*w2+w3*w3)
  FaultCell ca,cb,cl,cr; // nabors above, below, left and right
  FaultSkin skin; // if not null, the skin to which this cell belongs
  int i2m,i2p; // sample indices i2 for minus and plus sides of cell
  int i3m,i3p; // sample indices i3 for minus and plus sides of cell
  float[] emp; // array of alignment errors minus-plus
  float[] epm; // array of alignment errors plus-minus
  float smp; // shift from minus side to plus side of cell
  float spm; // shift from plus side to minus side of cell
  float t1,t2,t3; // fault slip vector

  FaultCell(float x1, float x2, float x3, float fl, float fp, float ft) {
    this.i1 = round(x1);
    this.i2 = round(x2);
    this.i3 = round(x3);
    this.x1 = x1; 
    this.x2 = x2; 
    this.x3 = x3;
    this.fl = fl; 
    this.fp = fp; 
    this.ft = ft;
    float[] u = faultDipVectorFromStrikeAndDip(fp,ft);
    float[] v = faultStrikeVectorFromStrikeAndDip(fp,ft);
    float[] w = faultNormalVectorFromStrikeAndDip(fp,ft);
    this.u1 = u[0]; this.u2 = u[1]; this.u3 = u[2];
    this.v1 = v[0]; this.v2 = v[1]; this.v3 = v[2];
    this.w1 = w[0]; this.w2 = w[1]; this.w3 = w[2];

    // Scale factor used to walk efficiently up and down a skin.
    ws23 = 1.0f/sqrt(w2*w2+w3*w3);

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
   * Computes alignment errors and initializes shifts. Computes both
   * minus-plus (emp) and plus-minus (epm) errors. The minus-plus errors
   * correspond to differences between the sample value on the minus side of
   * this cell and those for the plus sides of cells up and down dip from
   * this cell. The plus-minus errors are defined similarly.
   * <p>
   * Uses specified slopes to initialize both minus-plus and plus-minus shifts
   * to compensate for the fact that shifts are estimated using image samples
   * located a horizontal distance d away from this cell.
   * <p>
   * For lags where image sample values are unavailable, say, near surface
   * boundaries, errors are extrapolated from other lags, but are negated, so
   * that extrapolated errors can be detected and modified later after errors
   * for all cells in the surface have been computed.
   */
  /*
  void computeErrorsAndInitShifts(
    int lmax, float d, float[][][] f, float[][][] p2, float[][][] p3) 
  {
    int n1 = f[0][0].length;
    int nlag;
    float y1,y2,y3,gm,gp;
    Cell ci;

    // Errors for lag zero.
    emp = new float[lmax+1+lmax];
    epm = new float[lmax+1+lmax];
    float us = q.us23, u1 = us*q.u1, u2 = us*q.u2, u3 = us*q.u3;
    float ut = (i1-q.c1)*u1;
    float x1 = i1, x2 = q.c2-ut*u2, x3 = q.c3-ut*u3;

    float x2m = x2-d*w2, x3m = x3-d*w3;
    float x2p = x2+d*w2, x3p = x3+d*w3;
    float fm = valueAt(x1,x2m,x3m,f);
    float fp = valueAt(x1,x2p,x3p,f);
    float p2m = valueAt(x1,x2m,x3m,p2);
    float p2p = valueAt(x1,x2p,x3p,p2);
    float p3m = valueAt(x1,x2m,x3m,p3);
    float p3p = valueAt(x1,x2p,x3p,p3);
    float empl = emp[lmax] = error(fm,fp);
    float epml = epm[lmax] = error(fp,fm);

    // Initial shifts compensate for horizontal distance d.
    float s23 = d*((p2m+p2p)*w2+(p3m+p3p)*w3);
    smp = -s23;
    spm =  s23;

    // Errors for samples south; any extrapolated errors are negative.
    us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
    y1 = x1+1.0f; y2 = x2-u1*u2; y3 = x3-u1*u3;
    qi = quadSouth(q,y1,y2,y3);
    nlag = min(lmax,n1-1-i1);
    for (int ilag=1; ilag<=lmax; ++ilag) {
      if (qi!=null && ilag<=nlag) {
        us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
        ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
        y2 -= ut*u2; y3 -= ut*u3;
        gm = valueAt(y1,y2-d*u2,y3-d*u3,f);
        gp = valueAt(y1,y2+d*u2,y3+d*u3,f);
        empl = emp[lmax+ilag] = error(fm,gp);
        epml = epm[lmax+ilag] = error(fp,gm);
        y1 += 1.0f; y2 -= u1*u2; y3 -= u1*u3;
        qi = quadSouth(qi,y1,y2,y3);
      } else {
        emp[lmax+ilag] = -empl;
        epm[lmax+ilag] = -epml;
      }
    }

    // Errors for samples north; any extrapolated errors are negative.
    us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
    y1 = x1-1.0f; y2 = x2+u1*u2; y3 = x3+u1*u3;
    qi = quadNorth(q,y1,y2,y3);
    nlag = min(lmax,i1);
    for (int ilag=1; ilag<=lmax; ++ilag) {
      if (qi!=null && ilag<=nlag) {
        us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
        ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
        y2 -= ut*u2; y3 -= ut*u3;
        gm = valueAt(y1,y2-d*u2,y3-d*u3,f);
        gp = valueAt(y1,y2+d*u2,y3+d*u3,f);
        empl = emp[lmax-ilag] = error(fm,gp);
        epml = epm[lmax-ilag] = error(fp,gm);
        y1 -= 1.0f; y2 += u1*u2; y3 += u1*u3;
        qi = quadNorth(qi,y1,y2,y3);
      } else {
        emp[lmax-ilag] = -empl;
        epm[lmax-ilag] = -epml;
      }
    }
  }
  */

  /////////////////////////////////////////////////////////////////////////
  // private

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
    public int n;
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
