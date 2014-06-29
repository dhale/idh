/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import static edu.mines.jtk.util.ArrayMath.*;
import static fah.FaultUtil.*;

/**
 * A quad in a fault surface.
 * A quad references exactly four nodes and up to four quad nabors. Quad nodes
 * (na,nb,nc,nd) are ordered counter-clockwise, when viewed from the front
 * side of a quad, that is, from a viewpoint located at the tip of the quad
 * normal vector.
 * <p>
 * Each quad intersects exactly one edge in a 3D image sampling grid.
 * The 3D grid edges have three possible orientations: one vertical
 * edge in the direction of axis 1, and two horizontal edges in the
 * directions of axes 2 and 3. If a quad intersects a horizontal edge,
 * we say that the quad is "vertical", even though its four nodes may
 * not lie within a vertical plane. Vertical quads are most important
 * because we seek to estimate fault throws for faults that we assume
 * are more vertical than horizontal.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.03
 */
public class FaultQuad {

  int edge; // axis (1, 2, or 3) of edge intersected
  int i1,i2,i3; // lower indices of edge intersected
  FaultNode na,nb,nc,nd; // four nodes referenced by this quad
  FaultQuad qa,qb,qc,qd; // four quad nabors correspond to nodes
  float c1,c2,c3; // quad center point
  float u1,u2,u3; // quad normal vector
  float us23; // scale factor 1/sqrt(u2*u2+u3*u3)
  Vert v; // auxiliary data, if quad is vertical

  /**
   * Gets auxiliary data associated with this quad if it is vertical.
   * This method should be used only after all quads in a surface
   * have been linked and oriented. Quad nabors within the vert
   * data are not the same as quad nabors qa, qb, qc and qd.
   */
  Vert getVert() {
    if (v==null && isVertical())
      v = new Vert(this);
    return v;
  }

  /**
   * Constructs a quad that references four specified nodes.
   */
  FaultQuad(
      int edge, int i1, int i2, int i3, 
      FaultNode na, FaultNode nb, FaultNode nc, FaultNode nd) {
    this.edge = edge;
    this.i1 = i1;
    this.i2 = i2;
    this.i3 = i3;
    this.na = na;
    this.nb = nb;
    this.nc = nc;
    this.nd = nd;
  }

  public String toString() {
    return "("+c1+","+c2+","+c3+")";
  }

  boolean isVertical() {
    return edge!=1;
  }

  /**
   * Computes the plane for this quad if not already computed.
   */
  void computePlane() {
    if (u1==0.0f && u2==0.0f && u3==0.0f) {
      float ax1 = na.x1, ax2 = na.x2, ax3 = na.x3;
      float bx1 = nb.x1, bx2 = nb.x2, bx3 = nb.x3;
      float cx1 = nc.x1, cx2 = nc.x2, cx3 = nc.x3;
      float dx1 = nd.x1, dx2 = nd.x2, dx3 = nd.x3;
      float ca1 = cx1-ax1;
      float ca2 = cx2-ax2;
      float ca3 = cx3-ax3;
      float db1 = dx1-bx1;
      float db2 = dx2-bx2;
      float db3 = dx3-bx3;
      c1 = 0.25f*(ax1+bx1+cx1+dx1);
      c2 = 0.25f*(ax2+bx2+cx2+dx2);
      c3 = 0.25f*(ax3+bx3+cx3+dx3);
      u1 = ca3*db2-ca2*db3; // remember that sample coordinates
      u2 = ca1*db3-ca3*db1; // (x3,x2,x1) correspond to (x,y,z)
      u3 = ca2*db1-ca1*db2; // in right-handed coordinate system
      float us = sqrt(u1*u1+u2*u2+u3*u3);
      if (us>0.0f) {
        us = 1.0f/us;
        u1 *= us;
        u2 *= us;
        u3 *= us;
        us23 = sqrt(u2*u2+u3*u3);
        if (us23>0.0f)
          us23 = 1.0f/us23;
      }
    }
  }

  /**
   * Returns true, if this quad is good; false, otherwise.
   * The current criterion for goodness is that the quad normal
   * vector is consistent with the fault normal vectors stored
   * in each of the four nodes referenced by this quad. Here,
   * consistency depends on only the cosine-squared of the angle
   * between normal vectors, so that this test can be used before
   * this quad has been oriented.
   */
  boolean isGood() {
    computePlane();
    return goodNormal(na) &&
           goodNormal(nb) &&
           goodNormal(nc) &&
           goodNormal(nd);
  }
  private boolean goodNormal(FaultNode n) {
    float uu = u1*n.u1+u2*n.u2+u3*n.u3;
    //return uu*uu>0.80f;
    //return uu*uu>0.75f; // angle less than 30 degrees // OK for F3D
    //return uu*uu>0.50f; // angle less than 45 degrees
    //return uu*uu>0.25f; // angle less than 60 degrees
    return uu*uu>0.0f; // angle less than 90 degrees
  }

  /**
   * Returns true iff orientation of nabor is consistent with this quad.
   */
  boolean isOrientedLikeNabor(FaultQuad qn) {
    if (qn==qa) {
      if (na==qn.na) return nb==qn.nd; else
      if (na==qn.nb) return nb==qn.na; else
      if (na==qn.nc) return nb==qn.nb; else
      if (na==qn.nd) return nb==qn.nc; else
                     return false;
    } else if (qn==qb) {
      if (nb==qn.na) return nc==qn.nd; else
      if (nb==qn.nb) return nc==qn.na; else
      if (nb==qn.nc) return nc==qn.nb; else
      if (nb==qn.nd) return nc==qn.nc; else
                     return false;
    } else if (qn==qc) {
      if (nc==qn.na) return nd==qn.nd; else
      if (nc==qn.nb) return nd==qn.na; else
      if (nc==qn.nc) return nd==qn.nb; else
      if (nc==qn.nd) return nd==qn.nc; else
                     return false;
    } else if (qn==qd) {
      if (nd==qn.na) return na==qn.nd; else
      if (nd==qn.nb) return na==qn.na; else
      if (nd==qn.nc) return na==qn.nb; else
      if (nd==qn.nd) return na==qn.nc; else
                     return false;
    } else {
      return false;
    }
  }

  /**
   * Orients the specified quad nabor to be consistent with this quad.
   * Does nothing if the specified quad is null or is not a nabor of
   * this quad.
   */
  void orientNabor(FaultQuad qn) {
    if (!isOrientedLikeNabor(qn))
      qn.flip();
  }

  /**
   * Flips the orientation of this quad. Does not flip orientations
   * of the four nodes referenced by this quad.
   */
   void flip() {
     FaultNode ne = na; na = nd; nd = ne; // swap na-nd
     FaultNode nf = nb; nb = nc; nc = nf; // swap nb-nc
     FaultQuad qe = qa; qa = qc; qc = qe; // swap qa-qc
     u1 = -u1;
     u2 = -u2;
     u3 = -u3;
   }

  /**
   * Orients the four nodes referenced by this quad. Makes node normal 
   * and strike vectors consistent with the orientation of the quad. 
   * Nodes (na,nb,nc,nd) are ordered counter-clockwise, when viewed 
   * from the front side of a quad, from a viewpoint located at the 
   * tip of the quad normal vector.
   * <p>
   * Because nodes are shared by quads and their quad nabors, this
   * method should be called only after quads and their nabors have
   * been oriented.
   */
  void orientNodes() {
    orient(na);
    orient(nb);
    orient(nc);
    orient(nd);
  }
  void orient(FaultNode n) {
    computePlane();
    if (u1*n.u1+u2*n.u2+u3*n.u3<0.0f) {
      n.u1 = -n.u1;
      n.u2 = -n.u2;
      n.u3 = -n.u3;
      n.v1 = -n.v1;
      n.v2 = -n.v2;
      n.v3 = -n.v3;
    }
  }

  /**
   * Unlinks this quad from all of its quad nabors.
   */
  void unlink() {
    unlink(qa);
    unlink(qb);
    unlink(qc);
    unlink(qd);
  }
  void unlink(FaultQuad qn) {
    if (qn!=null) {
      if (this==qn.qa) qn.qa = null; else
      if (this==qn.qb) qn.qb = null; else
      if (this==qn.qc) qn.qc = null; else
      if (this==qn.qd) qn.qd = null;
      if (qn==qa) qa = null; else
      if (qn==qb) qb = null; else
      if (qn==qc) qc = null; else
      if (qn==qd) qd = null;
    }
  }

  void blocky() {
    na.blocky();
    nb.blocky();
    nc.blocky();
    nd.blocky();
  }

  /**
   * Auxiliary data associated with vertical quads.
   * Vertical quads are used to estimate fault displacements.
   * Before constructing these data, quads must be linked and oriented
   * so that east, north, west and south directions can be defined.
   */
  static class Vert {
    FaultQuad q; // the quad augmented by these data
    FaultQuad qe,qn,qw,qs; // quads east, north, west, and south
    int i2m,i2p; // sample indices i2 on minus and plus sides of quad
    int i3m,i3p; // sample indices i3 on minus and plus sides of quad
    float[] emp; // array of errors minus-plus
    float[] epm; // array of errors plus-minus
    float smp; // shift minus-plus
    float spm; // shift plus-minus
    float t1,t2,t3; // fault throw vector
    boolean tp; // true if fault throw is displacement of plus side

    /**
     * Constructs data for the specified vertical quad.
     */
    Vert(FaultQuad q) {
      assert q.isVertical();
      this.q = q;

      // Sample indices of minus-plus pairs (i2m,i2p) and (i3m,i3p).
      i2m = i2p = q.i2;
      i3m = i3p = q.i3;
      if (q.edge==2) {
        if (q.u2>0.0f)
          ++i2p;
        else
          ++i2m;
      } else if (q.edge==3) {
        if (q.u3>0.0f)
          ++i3p;
        else
          ++i3m;
      }

      // Use x1 coordinates to find quads east, north, west, south.
      FaultNode na = q.na, nb = q.nb, nc = q.nc, nd = q.nd;
      FaultQuad qa = q.qa, qb = q.qb, qc = q.qc, qd = q.qd;
      float ax1 = na.x1, bx1 = nb.x1, cx1 = nc.x1, dx1 = nd.x1;
      assert ax1!=bx1 || bx1!=cx1 || cx1!=dx1;
      if (bx1+cx1<=ax1+dx1) {
        qe = qa; qn = qb; qw = qc; qs = qd;
      } else if (cx1+dx1<=ax1+bx1) {
        qe = qb; qn = qc; qw = qd; qs = qa;
      } else if (ax1+dx1<=bx1+cx1) {
        qe = qc; qn = qd; qw = qa; qs = qb;
      } else if (ax1+bx1<=cx1+dx1) {
        qe = qd; qn = qa; qw = qb; qs = qc;
      } else { // one of the above must be true
        assert false;
      }
      assert qe==null || qe.isVertical():qe+" is vertical";
      assert qw==null || qw.isVertical():qw+" is vertical";

      // If north or south nabors are horizontal quads, then replace 
      // them with the quad on the other side of this one, if that
      // quad is vertical, or null, if that quad is also horizontal.
      // In linking vertical quads, we choose to bridge across one 
      // horizontal quad, if necessary, but never more than one.
      if (qn!=null && !qn.isVertical()) {
        qn = quadOnOtherSideOfNabor(qn);
        if (qn!=null && (!qn.isVertical() || qn.i1>=q.i1))
          qn = null;
      }
      if (qs!=null && !qs.isVertical()) {
        qs = quadOnOtherSideOfNabor(qs);
        if (qs!=null && (!qs.isVertical() || qs.i1<=q.i1))
          qs = null;
      }
    }
    private FaultQuad quadOnOtherSideOfNabor(FaultQuad qn) {
      if      (q==qn.qa) return qn.qc;
      else if (q==qn.qb) return qn.qd;
      else if (q==qn.qc) return qn.qa;
      else if (q==qn.qd) return qn.qb;
      else {assert false; return null;}
    }

    /**
     * Computes fault throw vector from vertical shifts.
     */
    void computeThrowFromShifts() {

      // Start on fault plane with integer x1 = i1.
      float us = q.us23, u1 = us*q.u1, u2 = us*q.u2, u3 = us*q.u3;
      float ut = (q.i1-q.c1)*u1;
      float x1 = q.i1, x2 = q.c2-ut*u2, x3 = q.c3-ut*u3;
      float y1 = x1, y2 = x2, y3 = x3;

      // Choose the largest (presumably non-negative) shift.
      // Throw corresponds to either the minus or plus side, 
      // depending on which shift is chosen.
      float s1 = max(smp,spm);
      tp = s1==smp;

      // If shift is small enough, do nothing.
      if (abs(s1)<=0.5f) {
        // do nothing
      }

      // Else, if shift is positive, walk south (downward) on quads.
      // Stop walking when we get to a quad corresponding to the 
      // shift, or when we can find no more quads to walk on.
      else if (s1>0.0f) {
        us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
        y1 += 1.0f; y2 -= u1*u2; y3 -= u1*u3;
        FaultQuad qi = quadSouth(q,y1,y2,y3);
        while (qi!=null && abs(y1-x1-s1)>0.5f) {
          us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
          ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
          y2 -= ut*u2; y3 -= ut*u3;
          y1 += 1.0f; y2 -= u1*u2; y3 -= u1*u3;
          qi = quadSouth(qi,y1,y2,y3);
        }
      }

      // Else, if shift is negative, walk north (upward) on quads. 
      // Stop walking when we get to a quad corresponding to the 
      // shift, or when we can find no more quads to walk on.
      else {
        us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
        y1 -= 1.0f; y2 += u1*u2; y3 += u1*u3;
        FaultQuad qi = quadNorth(q,y1,y2,y3);
        while (qi!=null && abs(y1-x1-s1)>0.5f) {
          us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
          ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
          y2 -= ut*u2; y3 -= ut*u3;
          y1 -= 1.0f; y2 += u1*u2; y3 += u1*u3;
          qi = quadNorth(qi,y1,y2,y3);
        }
      }

      // Now as close as we can get for the shift; compute the throw.
      float d1 = x1+s1-y1;
      y1 += d1; y2 -= d1*u1*u2; y3 -= d1*u1*u3;
      t1 = y1-x1; t2 = y2-x2; t3 = y3-x3;
    }

    float[] sampleFaultDip() {
      FaultQuad qi;
      float y1,y2,y3;

      // List of sample coordinates (x,y,z) = (x3,x2,x1).
      FloatList xyz = new FloatList();

      // First sample lies within fault plane with integer x1 = i1.
      float us = q.us23, u1 = us*q.u1, u2 = us*q.u2, u3 = us*q.u3;
      float ut = (q.i1-q.c1)*u1;
      float x1 = q.i1, x2 = q.c2-ut*u2, x3 = q.c3-ut*u3;
      xyz.add(x3); xyz.add(x2); xyz.add(x1);

      // Collect samples south; increment x1 by one for each sample.
      us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
      y1 = x1+1.0f; y2 = x2-u1*u2; y3 = x3-u1*u3;
      qi = quadSouth(q,y1,y2,y3);
      while (qi!=null) {
        us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
        ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
        y2 -= ut*u2; y3 -= ut*u3;
        xyz.add(y3); xyz.add(y2); xyz.add(y1);
        y1 += 1.0f; y2 -= u1*u2; y3 -= u1*u3;
        qi = quadSouth(qi,y1,y2,y3);
      }

      // Collect samples north; decrement x1 by one for each sample.
      us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
      y1 = x1-1.0f; y2 = x2+u1*u2; y3 = x3+u1*u3;
      qi = quadNorth(q,y1,y2,y3);
      while (qi!=null) {
        us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
        ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
        y2 -= ut*u2; y3 -= ut*u3;
        xyz.add(y3); xyz.add(y2); xyz.add(y1);
        y1 -= 1.0f; y2 += u1*u2; y3 += u1*u3;
        qi = quadNorth(qi,y1,y2,y3);
      }

      return xyz.trim();
    }

    /**
     * Computes alignment errors and initializes shifts.
     * Computes both minus-plus (emp) and plus-minus (epm) errors.
     * The minus-plus errors correspond to differences between
     * the sample value on the minus side of this quad and those 
     * for the plus sides of quads up and down dip from this quad.
     * The plus-minus errors are defined similarly. 
     * <p>
     * Uses the specified slopes to initialize both minus-plus and 
     * plus-minus shifts to compensate for the fact that shifts are 
     * estimated using image samples located a horizontal distance 
     * d away from this quad.
     * <p>
     * For lags where image sample values are unavailable, say, 
     * near surface boundaries, errors are extrapolated from other
     * lags, but are negated, so that extrapolated errors can be
     * detected and modified later after errors for all quads in
     * the surface have been computed.
     */
    void computeErrorsAndInitShifts(
      int lmax, float d, float[][][] f, float[][][] p2, float[][][] p3) 
    {
      int n1 = f[0][0].length;
      int i1 = q.i1;
      FaultQuad qi;
      int nlag;
      float y1,y2,y3,gm,gp;

      // Errors for lag zero.
      emp = new float[lmax+1+lmax];
      epm = new float[lmax+1+lmax];
      float us = q.us23, u1 = us*q.u1, u2 = us*q.u2, u3 = us*q.u3;
      float ut = (i1-q.c1)*u1;
      float x1 = i1, x2 = q.c2-ut*u2, x3 = q.c3-ut*u3;
      float x2m = x2-d*u2, x3m = x3-d*u3;
      float x2p = x2+d*u2, x3p = x3+d*u3;
      float fm = valueAt(x1,x2m,x3m,f);
      float fp = valueAt(x1,x2p,x3p,f);
      float p2m = valueAt(x1,x2m,x3m,p2);
      float p2p = valueAt(x1,x2p,x3p,p2);
      float p3m = valueAt(x1,x2m,x3m,p3);
      float p3p = valueAt(x1,x2p,x3p,p3);
      float empl = emp[lmax] = error(fm,fp);
      float epml = epm[lmax] = error(fp,fm);

      // Initial shifts compensate for horizontal distance d.
      float s23 = d*((p2m+p2p)*u2+(p3m+p3p)*u3);
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

    private static FaultQuad quadNorth(
        FaultQuad qi, float x1, float x2, float x3) {
      assert qi!=null && qi.isVertical();
      Vert vi = qi.getVert();
      FaultQuad qn = vi.qn;
      FaultQuad qe = vi.qe;
      FaultQuad qw = vi.qw;
      FaultQuad qen = (qe!=null)?qe.getVert().qn:null;
      FaultQuad qwn = (qw!=null)?qw.getVert().qn:null;
      return quadNearest(qn,qen,qwn,x1,x2,x3);
    }
    private static FaultQuad quadSouth(
        FaultQuad qi, float x1, float x2, float x3) {
      assert qi!=null && qi.isVertical();
      Vert vi = qi.getVert();
      FaultQuad qs = vi.qs;
      FaultQuad qe = vi.qe;
      FaultQuad qw = vi.qw;
      FaultQuad qes = (qe!=null)?qe.getVert().qs:null;
      FaultQuad qws = (qw!=null)?qw.getVert().qs:null;
      return quadNearest(qs,qes,qws,x1,x2,x3);
    }
    private static FaultQuad quadNearest(
        FaultQuad q1, FaultQuad q2, FaultQuad q3, 
        float x1, float x2, float x3) {
      float ds1 = ds(q1,x1,x2,x3);
      float ds2 = ds(q2,x1,x2,x3);
      float ds3 = ds(q3,x1,x2,x3);
      float dsm = min(ds1,ds2,ds3);
      if (dsm==ds1) { 
        return q1;
      } else if (dsm==ds2) {
        return q2;
      } else {
        return q3;
      }
    }
    private static float ds(FaultQuad q, float x1, float x2, float x3) {
      if (q==null)
        return Float.MAX_VALUE;
      float d1 = x1-q.c1;
      float d2 = x2-q.c2;
      float d3 = x3-q.c3;
      return d1*d1+d2*d2+d3*d3;
    }

    private static float error(float f, float g) {
      float fmg = f-g;
      return fmg*fmg;
    }

    private static float valueAt(float x1, float x2, float x3, float[][][]f) {
      int n1 = f[0][0].length;
      int n2 = f[0].length;
      int n3 = f.length;
      int i1 = max(0,min(n1-1,nint(x1)));
      int i2 = max(0,min(n2-1,nint(x2)));
      int i3 = max(0,min(n3-1,nint(x3)));
      return f[i3][i2][i1];
    }
    private static int nint(float x) {
      return (int)(x+0.5f);
    }
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
