/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Finds fault surfaces using fault likelihoods and orientations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.03
 */
public class FaultSurfer3 {

  public static float[][] slice1(int i1, float[][][] f) {
    int n2 = f[0].length;
    int n3 = f.length;
    float[][] fs = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        fs[i3][i2] = f[i3][i2][i1];
    return fs;
  }

  public float[][][] findShifts(
    double smax, Surf[] surfs, 
    float[][][] f, float[][][] p2, float[][][] p3) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] s = new float[n3][n2][n1];
    float[][][] c = new float[n3][n2][n1];
    for (Surf surf:surfs) {
      surf.computeShifts(smax,f,p2,p3);
      for (Quad quad:surf) {
        if (quad.isVertical()) {
          Vert v = quad.getVert();
          int i1 = quad.i1;
          int i2m = v.i2m, i2p = v.i2p;
          int i3m = v.i3m, i3p = v.i3p;
          s[i3m][i2m][i1] += v.smp;
          s[i3p][i2p][i1] += v.spm;
          c[i3m][i2m][i1] += 1.0f;
          c[i3p][i2p][i1] += 1.0f;
        }
      }
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (c[i3][i2][i1]>0.0f)
            s[i3][i2][i1] /= c[i3][i2][i1];
        }
      }
    }
    return s;
  }

  /**
   * Returns fault slip vectors for an array of fault surfaces.
   * Slip vectors are computed for only one side of faults.
   * Slips for the other side of the fault are set to zero.
   * Slips for samples not adjacent to faults are marked using
   * a specified value.
   * @param tmark the mark for slips not adjacent to a fault.
   * @param surfs array of fault surfaces.
   * @return array {s1,s2,s3} of components of fault slips.
   */
  public float[][][][] findSlips(float tmark, Surf[] surfs) {
    int n1 = _n1, n2 = _n2, n3 = _n3;
    float[][][] t1 = new float[n3][n2][n1];
    float[][][] t2 = new float[n3][n2][n1];
    float[][][] t3 = new float[n3][n2][n1];
    float[][][] ts = new float[n3][n2][n1];

    // Initially set all throw vectors to the specified mark.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t1[i3][i2][i1] = tmark;
          t2[i3][i2][i1] = tmark;
          t3[i3][i2][i1] = tmark;
        }
      }
    }

    // For all surfaces, ...
    for (Surf surf:surfs) {

      // Compute fault throws for all vertical quads in the surface.
      surf.computeThrows();

      // For all vertical quads in the surface, ...
      for (Quad quad:surf) {
        if (quad.isVertical()) {
          Vert v = quad.getVert();

          // Get sample indices for the zero and throw sides of the quad.
          int i1 = quad.i1;
          int i2z,i3z,i2t,i3t;
          if (v.tp) {
            i2z = v.i2m; i3z = v.i3m;
            i2t = v.i2p; i3t = v.i3p;
          } else {
            i2z = v.i2p; i3z = v.i3p;
            i2t = v.i2m; i3t = v.i3m;
          }

          // If throw on the zero side has not been set, zero it.
          if (t1[i3z][i2z][i1]==tmark) {
            t1[i3z][i2z][i1] = 0.0f;
            t2[i3z][i2z][i1] = 0.0f;
            t3[i3z][i2z][i1] = 0.0f;
          }

          // Set or accumulate throws on the throw side.
          if (t1[i3t][i2t][i1]==tmark) {
            t1[i3t][i2t][i1]  = v.t1;
            t2[i3t][i2t][i1]  = v.t2;
            t3[i3t][i2t][i1]  = v.t3;
            ts[i3t][i2t][i1]  = 1.0f;
          } else {
            t1[i3t][i2t][i1] += v.t1;
            t2[i3t][i2t][i1] += v.t2;
            t3[i3t][i2t][i1] += v.t3;
            ts[i3t][i2t][i1] += 1.0f;
          }
        }
      }
    }

    // Where more than one throw was accumulated, compute the average.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (ts[i3][i2][i1]>1.0f) {
            t1[i3][i2][i1] /= ts[i3][i2][i1];
            t2[i3][i2][i1] /= ts[i3][i2][i1];
            t3[i3][i2][i1] /= ts[i3][i2][i1];
          }
        }
      }
    }
    return new float[][][][]{t1,t2,t3};
  }

  ///////////////////////////////////////////////////////////////////////////
  // Node
  ///////////////////////////////////////////////////////////////////////////

  /**
   * A node is a point on a fault surface, with likelihood and orientation.
   */
  public static class Node {
    float fl,fp,ft; // fault likelihood, strike and dip
    float x1,x2,x3; // point on fault
    float u1,u2,u3; // normal vector
    float v1,v2,v3; // strike vector
    float w1,w2,w3; // dip vector
    void computeVectors() {
      float p = toRadians(fp);
      float t = toRadians(ft);
      float cp = cos(p); 
      float sp = sin(p);
      float ct = cos(t);
      float st = sin(t);
      this.u1 = -st;
      this.u2 = -sp*ct;
      this.u3 =  cp*ct;
      this.v1 = 0.0f;
      this.v2 = cp;
      this.v3 = sp;
      this.w1 = ct;
      this.w2 = -sp*st;
      this.w3 =  cp*st;
    }
    public String toString() {
      return "("+x1+","+x2+","+x3+")";
      //return "("+x1+","+x2+","+x3+"):("+fl+","+fp+","+ft+")";
    }
    void blocky() {
      x1 = (int)x1+0.5f;
      x2 = (int)x2+0.5f;
      x3 = (int)x3+0.5f;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Quad
  ///////////////////////////////////////////////////////////////////////////

  /**
   * A quad references exactly four nodes and up to four quad nabors.
   * Quad nodes (na,nb,nc,nd) are ordered counter-clockwise, when 
   * viewed from the front side of a quad, that is, from a viewpoint 
   * located at the tip of the quad normal vector.
   * <p>
   * Each quad intersects exactly one edge in a 3D image sampling grid.
   * The 3D grid edges have three possible orientations: one vertical
   * edge in the direction of axis 1, and two horizontal edges in the
   * directions of axes 2 and 3. If a quad intersects a horizontal edge,
   * we say that the quad is "vertical", even though its four nodes may
   * not lie within a vertical plane. Vertical quads are most important
   * because we seek to estimate fault throws for faults that we assume
   * are more vertical than horizontal.
   */
  public static class Quad {
    int edge; // axis (1, 2, or 3) of edge intersected
    int i1,i2,i3; // lower indices of edge intersected
    Node na,nb,nc,nd; // four nodes referenced by this quad
    Quad qa,qb,qc,qd; // four quad nabors correspond to nodes
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
    Quad(int edge, int i1, int i2, int i3, 
         Node na, Node nb, Node nc, Node nd) 
    {
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
    private boolean goodNormal(Node n) {
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
    boolean isOrientedLikeNabor(Quad qn) {
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
    void orientNabor(Quad qn) {
      if (!isOrientedLikeNabor(qn))
        qn.flip();
    }

    /**
     * Flips the orientation of this quad. Does not flip orientations
     * of the four nodes referenced by this quad.
     */
     void flip() {
       Node ne = na; na = nd; nd = ne; // swap na-nd
       Node nf = nb; nb = nc; nc = nf; // swap nb-nc
       Quad qe = qa; qa = qc; qc = qe; // swap qa-qc
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
    private void orient(Node n) {
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
    private void unlink(Quad qn) {
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
  }

  ///////////////////////////////////////////////////////////////////////////
  // Vert
  ///////////////////////////////////////////////////////////////////////////

  /**
   * Auxiliary data associated with vertical quads.
   * Vertical quads are used to estimate fault displacements.
   * Before constructing these data, quads must be linked and oriented
   * so that east, north, west and south directions can be defined.
   */
  public static class Vert {
    Quad q; // the quad augmented by these data
    Quad qe,qn,qw,qs; // quads east, north, west, and south
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
    Vert(Quad q) {
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
      Node na = q.na, nb = q.nb, nc = q.nc, nd = q.nd;
      Quad qa = q.qa, qb = q.qb, qc = q.qc, qd = q.qd;
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
    private Quad quadOnOtherSideOfNabor(Quad qn) {
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
        Quad qi = quadSouth(q,y1,y2,y3);
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
        Quad qi = quadNorth(q,y1,y2,y3);
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
      Quad qi;
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
      Quad qi;
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

    private static Quad quadNorth(Quad qi, float x1, float x2, float x3) {
      assert qi!=null && qi.isVertical();
      Vert vi = qi.getVert();
      Quad qn = vi.qn;
      Quad qe = vi.qe;
      Quad qw = vi.qw;
      Quad qen = (qe!=null)?qe.getVert().qn:null;
      Quad qwn = (qw!=null)?qw.getVert().qn:null;
      return quadNearest(qn,qen,qwn,x1,x2,x3);
    }
    private static Quad quadSouth(Quad qi, float x1, float x2, float x3) {
      assert qi!=null && qi.isVertical();
      Vert vi = qi.getVert();
      Quad qs = vi.qs;
      Quad qe = vi.qe;
      Quad qw = vi.qw;
      Quad qes = (qe!=null)?qe.getVert().qs:null;
      Quad qws = (qw!=null)?qw.getVert().qs:null;
      return quadNearest(qs,qes,qws,x1,x2,x3);
    }
    private static Quad quadNearest(
      Quad q1, Quad q2, Quad q3, float x1, float x2, float x3) 
    {
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
    private static float ds(Quad q, float x1, float x2, float x3) {
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

  ///////////////////////////////////////////////////////////////////////////
  // Surf
  ///////////////////////////////////////////////////////////////////////////

  /**
   * A fault surface consists of a set of linked and oriented quads.
   */
  public class Surf implements Iterable<Quad> {
    Surf(Quad[] quads) {
      _quads = quads;
    }
    public int size() {
      return _quads.length;
    }
    public Iterator<Quad> iterator() {
      return Arrays.asList(_quads).iterator();
    }
    public float[][] getXyzUvwRgb() {
      return FaultSurfer3.getXyzUvwRgb(_quads,_flmin);
    }
    public float[][] getXyzUvwRgbShifts(float smax) {
      return FaultSurfer3.getXyzUvwRgbShifts(_quads,smax);
    }

    /**
     * Replaces all fault likelihoods with the average for the surface.
     */
    public void smooth() {
      float flsum = 0.0f;
      int nlsum = 0;
      for (Quad q:_quads) {
        flsum += q.na.fl;
        flsum += q.nb.fl;
        flsum += q.nc.fl;
        flsum += q.nd.fl;
        nlsum += 4;
      }
      float flavg = flsum/nlsum;
      for (Quad q:_quads) {
        q.na.fl = flavg;
        q.nb.fl = flavg;
        q.nc.fl = flavg;
        q.nd.fl = flavg;
      }
    }

    /**
     * Orients this surface. Orientation is computed such that the 
     * average of u2 components of quad normal vectors is non-negative.
     */
    public void orientU2() {
      float u2sum = 0.0f;
      for (Quad q:_quads)
        u2sum += q.u2;
      if (u2sum<0.0) {
        for (Quad q:_quads)
          q.flip();
        for (Quad q:_quads)
          q.orientNodes();
      }
    }

    /**
     * Orients this surface. Orientation is computed such that the 
     * average of u1 components of quad normal vectors is negative.
     */
    public void orient() {
      float u1sum = 0.0f;
      for (Quad q:_quads)
        u1sum += q.u1;
      if (u1sum>0.0) {
        for (Quad q:_quads)
          q.flip();
        for (Quad q:_quads)
          q.orientNodes();
      }
    }

    /**
     * Moves all nodes to make a topologically equivalent blocky surface.
     * A blocky surface is useful to see which adjacent image samples lie
     * on opposite sides of a fault.
     */
    public void blocky() {
      for (Quad q:_quads)
        q.blocky();
    }

    public float[] sampleFaultDip() {
      float c1 = 0.0f;
      float c2 = 0.0f;
      float c3 = 0.0f;
      int nq = 0;
      for (Quad quad:_quads) {
        c1 += quad.c1;
        c2 += quad.c2;
        c3 += quad.c3;
        nq += 1;
      }
      c1 /= nq;
      c2 /= nq;
      c3 /= nq;
      float dsmin = Float.MAX_VALUE;
      Quad qmin = null;
      for (Quad quad:_quads) {
        if (quad.isVertical()) {
          float d1 = quad.c1-c1;
          float d2 = quad.c2-c2;
          float d3 = quad.c3-c3;
          float ds = d1*d1+d2*d2+d3*d3;
          if (ds<dsmin) {
            dsmin = ds;
            qmin = quad;
          }
        }
      }
      FloatList xyz = new FloatList();
      Vert vert = qmin.getVert();
      xyz.add(vert.sampleFaultDip());
      while (vert.qe!=null) {
        vert = vert.qe.getVert();
        xyz.add(vert.sampleFaultDip());
      }
      vert = qmin.getVert();
      while (vert.qw!=null) {
        vert = vert.qw.getVert();
        xyz.add(vert.sampleFaultDip());
      }
      return xyz.trim();
    }

    /**
     * Computes fault shifts for all vertical quads in this surface.
     */
    public void computeShifts(
      double smax, float[][][] f, float[][][] p2, float[][][] p3) 
    {
      int lmax = (int)smax;
      final float d = 2.0f;
      computeErrorsAndInitShifts(lmax,d,f,p2,p3);
      Quad[][][] qewns = getQuadsEwNs();
      Quad[][] qew = qewns[0];
      Quad[][] qns = qewns[1];
      extrapolateErrors(qns);
      DynamicWarping dw = new DynamicWarping(-lmax,lmax);
      dw.setStretchMax(0.25,0.25);
      findShiftsMP(dw,qew,qns);
      findShiftsPM(dw,qew,qns);
      filterShifts();
      smoothShifts();
      smoothShifts();
    }
    private void computeErrorsAndInitShifts(
      int lmax, float d, float[][][] f, float[][][] p2, float[][][] p3) 
    {
      for (Quad quad:_quads) {
        if (quad.isVertical()) {
          Vert vert = quad.getVert();
          vert.computeErrorsAndInitShifts(lmax,d,f,p2,p3);
        }
      }
    }
    private void smoothShifts() {
      for (Quad q:this) {
        if (q.isVertical()) {
          Vert v = q.getVert();
          float smp = 0.0f;
          float spm = 0.0f;
          float css = 0.0f;
          for (Quad qn:new Quad[]{v.qs,v.qn,v.qw,v.qe}) {
            if (qn!=null) {
              Vert vn = qn.getVert();
              smp += v.smp+vn.smp;
              spm += v.spm+vn.spm;
              css += 2.0f;
            }
          }
          v.smp = smp/css;
          v.spm = spm/css;
        }
      }
    }
    private void filterShifts() {
      for (Quad q:this) {
        if (q.isVertical()) {
          Vert v = q.getVert();
          if (v.smp*v.spm>0.0f) {
            v.smp = 0.0f;
            v.spm = 0.0f;
          }
        }
      }
    }

    /**
     * Computes throw vectors for all vertical quads in this surface.
     * Uses shifts that must already have been computed for each 
     * vertical quad.
     */
    public void computeThrows() {
      for (Quad quad:_quads) {
        if (quad.isVertical()) {
          Vert vert = quad.getVert();
          vert.computeThrowFromShifts();
        }
      }
    }

    /**
     * Returns east-west and north-south arrays of quads.
     */
    private Quad[][][] getQuadsEwNs() {
      HashSet<Quad> qewSet = new HashSet<Quad>(size());
      HashSet<Quad> qnsSet = new HashSet<Quad>(size());
      ArrayList<Quad[]> qewList = new ArrayList<Quad[]>();
      ArrayList<Quad[]> qnsList = new ArrayList<Quad[]>();
      for (Quad q:this) {
        if (q.isVertical()) {
          if (!qewSet.contains(q)) {
            Quad qi = q;
            for (Quad qe=qi.getVert().qe; qe!=null; qe=qi.getVert().qe)
              qi = qe;
            ArrayList<Quad> ql = new ArrayList<Quad>();
            for (; qi!=null; qi=qi.getVert().qw) {
              ql.add(qi);
              qewSet.add(qi);
            }
            qewList.add(ql.toArray(new Quad[0]));
          }
          if (!qnsSet.contains(q)) {
            Quad qi = q;
            for (Quad qn=qi.getVert().qn; qn!=null; qn=qi.getVert().qn)
              qi = qn;
            ArrayList<Quad> ql = new ArrayList<Quad>();
            for (; qi!=null; qi=qi.getVert().qs) {
              ql.add(qi);
              qnsSet.add(qi);
            }
            qnsList.add(ql.toArray(new Quad[0]));
          }
        }
      }
      Quad[][] qew = qewList.toArray(new Quad[0][]);
      Quad[][] qns = qnsList.toArray(new Quad[0][]);
      return new Quad[][][]{qew,qns};
    }

    private Quad[] _quads;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  /**
   * Constructs a fault surfer for specified likelihoods and orientations.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   */
  public FaultSurfer3(float[][][][] flpt) {
    _fl = flpt[0];
    _fp = flpt[1];
    _ft = flpt[2];
    _n1 = _fl[0][0].length;
    _n2 = _fl[0].length;
    _n3 = _fl.length;
  }

  /**
   * Sets a lower bound on fault likelihoods.
   * The default lower bounds is 0.1.
   * @param flmin the lower bound
   */
  public void setThreshold(double flmin) {
    _flmin = (float)flmin;
  }

  /**
   * Returns array of quads for faults, ridge surfaces in fault likelihood.
   * Returned quads may share nodes, but are neither linked nor oriented.
   */
  public Quad[] findQuads() {
    int n1 = _n1, n2 = _n2, n3 = _n3;
    float[][][] fl = _fl, fp = _fp, ft = _ft;

    // Gaussian smoothing for robust derivatives of fault likelihoods.
    float[][][] fs = new float[n3][n2][n1];
    _rgf.apply000(fl,fs);

    // Buffers to hold vectors h for two constant-i3 slices.
    float[][][] h1 = new float[2][n2][n1];
    float[][][] h2 = new float[2][n2][n1];
    float[][][] h3 = new float[2][n2][n1];

    // Array of nodes; will be non-null where faults are found.
    Node[][][] nodes = new Node[n3][n2][n1];

    // List of quads within fault surfaces.
    ArrayList<Quad> quads = new ArrayList<Quad>();

    // Process first constant-i3 slice.
    int i3 = INSET23;
    processSlices(i3,i3,_flmin,fs,fl,fp,ft,h1,h2,h3,nodes,quads);

    // Process all remaining pairs of adjacent constant-i3 slices.
    for (int j3=i3+1; j3<n3-INSET23; ++i3,++j3)
      processSlices(i3,j3,_flmin,fs,fl,fp,ft,h1,h2,h3,nodes,quads);

    // Complete computation of values for all non-null nodes. 
    completeNodes(nodes);

    // Clean the list of quads by removing any bad quads.
    trace("findQuads: before cleaning, nquad = "+quads.size());
    quads = cleanQuads(quads);
    trace("            after cleaning, nquad = "+quads.size());
    
    return quads.toArray(new Quad[0]);
  }

  private static class TwoQuads {
    Quad q1,q2;
    TwoQuads(Quad q) {
      q1 = q;
    }
    void linkOrUnlink(Edge e, Quad q) {
      if (q2==null) {
        e.link(q2=q,q1);
      } else if (q1!=null) {
        q1.unlink(q2);
        q1 = null;
      }
    }
  }

  /**
   * Returns an array of linked and oriented quads.
   */
  public Quad[] linkQuads(Quad[] quads) {

    // Link two quads that share an edge. However, if an edge is shared
    // by more than two quads, do not link any of them. This latter case
    // should be rare, so we first link quads and then later unlink them
    // if necessary.
    HashMap<Edge,TwoQuads> map = new HashMap<Edge,TwoQuads>(4*quads.length);
    for (Quad q:quads) {
      Edge eab = new Edge(q.na,q.nb);
      TwoQuads qab = map.get(eab);
      if (qab!=null) qab.linkOrUnlink(eab,q);
      else map.put(eab,new TwoQuads(q));
      Edge ebc = new Edge(q.nb,q.nc);
      TwoQuads qbc = map.get(ebc);
      if (qbc!=null) qbc.linkOrUnlink(ebc,q);
      else map.put(ebc,new TwoQuads(q));
      Edge ecd = new Edge(q.nc,q.nd);
      TwoQuads qcd = map.get(ecd);
      if (qcd!=null) qcd.linkOrUnlink(ecd,q);
      else map.put(ecd,new TwoQuads(q));
      Edge eda = new Edge(q.nd,q.na);
      TwoQuads qda = map.get(eda);
      if (qda!=null) qda.linkOrUnlink(eda,q);
      else map.put(eda,new TwoQuads(q));
    }
    trace("linkQuads: after linking all quads"); printStats(quads);

    // Unlink any quads that are folded on top of one another.
    // Normal vectors of quads and quad nabors must be consistent.
    unlinkFoldedQuads(quads);
    trace("linkQuads:  after unlinking folded quads"); printStats(quads);

    // Unlink any quads with insufficient nabors. Quad nabors are
    // insufficent if they do not reference at least one common node.
    // This criterion eliminates fins and bridges with width one quad.
    unlinkSkinnyQuads(quads);
    trace("linkQuads:  after unlinking skinny quads"); printStats(quads);

    // Remove any single quads that have no nabors.
    quads = removeSingleQuads(quads);
    trace("linkQuads:  after removing single quads"); printStats(quads);

    return quads;
  }

  /**
   * Randomly permutes the specified array of quads.
   * @param quads array of quads to be shuffled, in place.
   */
  public void shuffleQuads(Quad[] quads) {
    Random r = new Random(3);
    int n = quads.length;
    for (int i=n-1; i>0; --i) {
      int j = r.nextInt(i+1);
      Quad quadi = quads[i];
      quads[i] = quads[j];
      quads[j] = quadi;
    }
  }

  /**
   * Returns surfaces, collections of linked and oriented quads.
   * Assumes that the quads have already been linked to their nabors, 
   * and that each surface comprised of quads is orientable, so that 
   * all quads and their quad nabors can have consistent orientations.
   */
  public Surf[] findSurfs(Quad[] quads) {

    // List of surfs found while orienting quads.
    ArrayList<Surf> surfs = new ArrayList<Surf>();

    // Set of quads that have been oriented.
    HashSet<Quad> set = new HashSet<Quad>(quads.length);

    // Stack of quads that have been oriented, but with nabors to visit.
    ArrayDeque<Quad> stack = new ArrayDeque<Quad>();

    // Count surfaces that are not orientable.
    int nsno = 0;

    // For all quads, ...
    for (Quad quad:quads) {

      // If quad has not yet been oriented, ...
      if (!set.contains(quad)) {

        // Begin collecting quads for a new surface.
        ArrayList<Quad> list = new ArrayList<Quad>();

        // This quad will determine the orientation of its nabors.
        // Add this quad to the set of oriented quads, to the list
        // of quads for the current surface and to the stack of quads 
        // with nabors to visit.
        set.add(quad);
        list.add(quad);
        stack.push(quad);

        // Count of quads that cannot be oriented consistently.
        int nqno = 0;

        // While quad nabors have not yet been oriented, ...
        while (!stack.isEmpty()) {

          // Quad q on the stack is oriented, but has nabors to visit.
          Quad q = stack.pop();

          // Orient the nodes for the oriented quad q.
          q.orientNodes();

          // For all quad nabors qn of the oriented quad q, ...
          Quad[] qns = new Quad[]{q.qa,q.qb,q.qc,q.qd};
          for (Quad qn:qns) {
            if (qn!=null) {

              // If nabor has not yet been oriented, orient it, add it to
              // the set of oriented quads, and to the stack of quads with 
              // nabors to visit.
              if (!set.contains(qn)) {
                q.orientNabor(qn);
                set.add(qn);
                list.add(qn);
                stack.push(qn);
              }

              // Otherwise, if nabor has been oriented, then determine 
              // whether or not its orientation is consistent with that
              // for this quad. If not, then surface is not orientable.
              else if (!q.isOrientedLikeNabor(qn)) {
                ++nqno;
                q.unlink(qn);
                qn.unlink(q);
              }
            }
          }
        }

        // Construct and orient a new surface, and add list of surfaces.
        if (nqno>0) {
          ++nsno;
          trace("findSurfs: surf "+surfs.size()+" not orientable!");
        }
        Surf surf = new Surf(list.toArray(new Quad[0]));
        surf.orient();
        surfs.add(surf);
      }
    }
    trace("findSurfs: nquad="+quads.length+
                    " nsurf="+surfs.size()+
                    " nsno="+nsno+
                    " nset="+set.size());
    return surfs.toArray(new Surf[0]);
  }

  /**
   * Returns an array of surfs with the specified minimum number of quads.
   */
  public static Surf[] getSurfsWithSize(Surf[] surfs, int minSize) {
    ArrayList<Surf> surfList = new ArrayList<Surf>();
    for (Surf surf:surfs) {
      if (surf.size()>minSize)
        surfList.add(surf);
    }
    trace("getSurfsWithSize: input = "+surfs.length);
    trace("                 output = "+surfList.size());
    return surfList.toArray(new Surf[0]);
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of quad coordinates, normals and colors.
   */
  public static float[][] getXyzUvwRgb(Quad[] quads, float flmin) {
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList fcl = new FloatList();
    for (Quad quad:quads) {
      Node na = quad.na;
      Node nb = quad.nb;
      Node nc = quad.nc;
      Node nd = quad.nd;
      xyz.add(na.x3); xyz.add(na.x2); xyz.add(na.x1); fcl.add(na.fl);
      xyz.add(nb.x3); xyz.add(nb.x2); xyz.add(nb.x1); fcl.add(nb.fl);
      xyz.add(nc.x3); xyz.add(nc.x2); xyz.add(nc.x1); fcl.add(nc.fl);
      xyz.add(nd.x3); xyz.add(nd.x2); xyz.add(nd.x1); fcl.add(nd.fl);
      uvw.add(na.u3); uvw.add(na.u2); uvw.add(na.u1);
      uvw.add(nb.u3); uvw.add(nb.u2); uvw.add(nb.u1);
      uvw.add(nc.u3); uvw.add(nc.u2); uvw.add(nc.u1);
      uvw.add(nd.u3); uvw.add(nd.u2); uvw.add(nd.u1);
    }
    float[] fc = fcl.trim();
    float fcmin = flmin;
    float fcmax = 1.0f;
    ColorMap cmap = new ColorMap(fcmin,fcmax,ColorMap.JET);
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of quad coordinates, normals and colors.
   * Uses color to represent fault shifts.
   */
  public static class Sum {
    float s;
    int n;
    void add(float v) {
      s += v;
      n += 1;
    }
    float avg() {
      return n>0?s/n:0.0f;
    }
  }
  private static void acc(HashMap<Node,Sum> map, Node node, float v) {
    Sum s = map.get(node);
    if (s==null)
      map.put(node,s=new Sum());
    s.add(v);
  }
  private static float avg(HashMap<Node,Sum> map, Node node) {
    Sum s = map.get(node);
    return (s!=null)?s.avg():0.0f;
  }
  public static float[][] getXyzUvwRgbShifts(Quad[] quads, float smax) {
    HashMap<Node,Sum> map = new HashMap<Node,Sum>(2*quads.length);
    for (Quad quad:quads) {
      if (quad.isVertical()) {
        Vert vert = quad.getVert();
        //float s = smax>0.0?vert.smp:vert.spm;
        float s;
        if (smax>0.0f) {
          s = vert.smp>0.0f?vert.smp:vert.spm;
        } else {
          s = vert.spm>0.0f?vert.spm:vert.smp;
        }
        acc(map,quad.na,s);
        acc(map,quad.nb,s);
        acc(map,quad.nc,s);
        acc(map,quad.nd,s);
      }
    }
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList scl = new FloatList();
    for (Quad quad:quads) {
      Node na = quad.na; float sa = avg(map,na);
      Node nb = quad.nb; float sb = avg(map,nb);
      Node nc = quad.nc; float sc = avg(map,nc);
      Node nd = quad.nd; float sd = avg(map,nd);
      xyz.add(na.x3); xyz.add(na.x2); xyz.add(na.x1); scl.add(sa);
      xyz.add(nb.x3); xyz.add(nb.x2); xyz.add(nb.x1); scl.add(sb);
      xyz.add(nc.x3); xyz.add(nc.x2); xyz.add(nc.x1); scl.add(sc);
      xyz.add(nd.x3); xyz.add(nd.x2); xyz.add(nd.x1); scl.add(sd);
      uvw.add(na.u3); uvw.add(na.u2); uvw.add(na.u1);
      uvw.add(nb.u3); uvw.add(nb.u2); uvw.add(nb.u1);
      uvw.add(nc.u3); uvw.add(nc.u2); uvw.add(nc.u1);
      uvw.add(nd.u3); uvw.add(nd.u2); uvw.add(nd.u1);
    }
    float[] sc = scl.trim();
    float csmin = 0.0f;
    float csmax = abs(smax);
    ColorMap cmap = new ColorMap(csmin,csmax,ColorMap.JET);
    float[] rgb = cmap.getRgbFloats(sc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _fl; // fault likelihoods in [0,1]
  private float[][][] _fp; // fault strikes phi (degrees)
  private float[][][] _ft; // fault dips theta (degrees)
  private int _n1,_n2,_n3; // numbers of image samples
  private float _flmin = 0.1f; // fault likelihood threshold

  // Used to compute 1st and 2nd derivatives of fault likelihoods.
  private RecursiveGaussianFilter _rgf = new RecursiveGaussianFilter(1.0);

  /**
   * An edge na---nb is used to link quads that share nodes its na and nb.
   * Note that edge na---nb equals the edge nb---na. This definition of 
   * the method equals is necessary because quads must be linked before 
   * they can be oriented. 
   * <p>
   * Only two quads may be linked by an edge. Attempts to link more quads
   * will cause the two quads linked via this edge to be unlinked.
   */
  private static class Edge {
    Node na,nb;
    Edge(Node na, Node nb) {
      this.na = na;
      this.nb = nb;
    }
    public boolean equals(Object edge) {
      Edge that = (Edge)edge;
      return this.na==that.na && this.nb==that.nb ||
             this.na==that.nb && this.nb==that.na;
    }
    public int hashCode() {
      return na.hashCode()^nb.hashCode();
    }
    void link(Quad q1, Quad q2) {
      // NOTE: the logic in this method depends on the order in which
      // nodes are passed to the constructors for the quads q1 and q2.
      if (na==q1.na) {
        assert nb==q1.nb;
        q1.qa = q2;
        if (na==q2.na && nb==q2.nb) {
          q2.qa = q1;
        } else if (na==q2.nb && nb==q2.nc) {
          q2.qb = q1;
        } else if (na==q2.nd && nb==q2.nc) {
          q2.qc = q1;
        } else if (na==q2.na && nb==q2.nd) {
          q2.qd = q1;
        } else { // one of the above must be true!
          assert false;
        }
      } else if (na==q1.nb) {
        assert nb==q1.nc;
        q1.qb = q2;
        if (na==q2.na && nb==q2.nb) {
          q2.qa = q1;
        } else if (na==q2.nb && nb==q2.nc) {
          q2.qb = q1;
        } else if (na==q2.nd && nb==q2.nc) {
          q2.qc = q1;
        } else if (na==q2.na && nb==q2.nd) {
          q2.qd = q1;
        } else { // one of the above must be true!
          assert false;
        }
      } else if (na==q1.nc) {
        assert nb==q1.nd;
        q1.qc = q2;
        if (na==q2.nb && nb==q2.na) {
          q2.qa = q1;
        } else if (na==q2.nc && nb==q2.nb) {
          q2.qb = q1;
        } else if (na==q2.nc && nb==q2.nd) {
          q2.qc = q1;
        } else if (na==q2.nd && nb==q2.na) {
          q2.qd = q1;
        } else { // one of the above must be true!
          assert false;
        }
      } else if (na==q1.nd) {
        assert nb==q1.na;
        q1.qd = q2;
        if (na==q2.nb && nb==q2.na) {
          q2.qa = q1;
        } else if (na==q2.nc && nb==q2.nb) {
          q2.qb = q1;
        } else if (na==q2.nc && nb==q2.nd) {
          q2.qc = q1;
        } else if (na==q2.nd && nb==q2.na) {
          q2.qd = q1;
        } else { // one of the above must be true!
          assert false;
        }
      } else { // one of the above must be true!
        assert false;
      }
    }
  }

  // Nodes are inside a box inset by this many samples from image bounds.
  // This number reduces artifacts caused by image boundaries. Must be
  // at least one sample, to avoid array index out of bounds exceptions.
  // TODO: clean up computation of fault likelihoods to reduce artifacts.
  // Should be able to get within one sample of image bounds, although
  // approximations to Gaussian derivatives will be less accurate there.
  private static final int INSET23 = 5;
  //private static final int INSET23 = 1;

  /**
   * Uses vectors h to process one pair of adjacent constant-i3 slices.
   * Arrays h1[i3%2], h2[i3%2] and h3[i3%2] are input to this method.
   * Arrays h1[j3%2], h2[j3%2] and h3[j3%2] are computed by this method.
   * These arrays correspond to image slices i3 and j3 = i3+1. Vectors 
   * h with components (h1,h2,h3) are used to detect intersections of 
   * faults with edges of the image sampling grid. For each such
   * intersection, this method constructs a quad referencing four 
   * nodes on the fault.
   */
  private static void processSlices(
    int i3, int j3, float flmin, float[][][] fs,
    float[][][] fl, float[][][] fp, float[][][] ft,
    float[][][] h1, float[][][] h2, float[][][] h3,
    Node[][][] nodes, ArrayList<Quad> quads)
  {
    int n1 = nodes[0][0].length;
    int n2 = nodes[0].length;
    int n3 = nodes.length;
   
    // Have vectors h for slice i3; must compute them for slice j3.
    computeVectorsH(j3,flmin,fs,fl,fp,ft,h1[j3%2],h2[j3%2],h3[j3%2]);

    // Construct nodes by looking for intersections with grid edges.
    // We assume that edges within slice i3 have already been checked.
    // Therefore, we need only check edges within the slice j3, and 
    // edges between the two slices i3 and j3.
    for (int i2=INSET23,j2=i2+1; i2<n2-INSET23; ++i2,++j2) {
      j2 = min(j2,n2-INSET23-1);
      for (int i1=0,j1=i1+1; i1<n1; ++i1,++j1) {
        j1 = min(j1,n1-1);

        // Edge i1---j1 in slice j3.
        if (fl[j3][i2][i1]>=flmin && fl[j3][i2][j1]>=flmin)
          processEdge(i1,i2,j3,j1,i2,j3,fl,fp,ft,h1,h2,h3,nodes,quads);

        // Edge i2---j2 in slice j3.
        if (fl[j3][i2][i1]>=flmin && fl[j3][j2][i1]>=flmin)
          processEdge(i1,i2,j3,i1,j2,j3,fl,fp,ft,h1,h2,h3,nodes,quads);

        // Edge i3---j3 between slices i3 and j3.
        if (fl[i3][i2][i1]>=flmin && fl[j3][i2][i1]>=flmin)
          processEdge(i1,i2,i3,i1,i2,j3,fl,fp,ft,h1,h2,h3,nodes,quads);
      }
    }
  }

  /**
   * Computes vectors h = uu'g for one constant-i3 slice. 
   * Unit vector u is the fault normal vector; g is the gradient vector.
   */
  private static void computeVectorsH(
    int i3, float flmin, float[][][] fs,
    float[][][] fl, float[][][] fp, float[][][] ft,
    float[][] h1, float[][] h2, float[][] h3)
  {
    int n1 = fs[0][0].length;
    int n2 = fs[0].length;

    // For all samples in this constant-i3 slice, ...
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=1; i1<n1-1; ++i1) {

        // Components of vector h are initially zero.
        float h1i = 0.0f;
        float h2i = 0.0f;
        float h3i = 0.0f;

        // If fault likelihood not less than threshold, ...
        if (fl[i3][i2][i1]>=flmin) {

          // Normal vector u from fault strike and dip.
          float pr = toRadians(fp[i3][i2][i1]);
          float tr = toRadians(ft[i3][i2][i1]);
          float cp = cos(pr);
          float sp = sin(pr);
          float ct = cos(tr);
          float st = sin(tr);
          float u1 = -st;
          float u2 = -sp*ct;
          float u3 =  cp*ct;

          // Gradient of fault likelihood.
          float g1 = 0.5f*(fs[i3][i2][i1+1]-fs[i3][i2][i1-1]);
          float g2 = 0.5f*(fs[i3][i2+1][i1]-fs[i3][i2-1][i1]);
          float g3 = 0.5f*(fs[i3+1][i2][i1]-fs[i3-1][i2][i1]);

          // Scaled dot product u'g.
          float ug = u1*g1+u2*g2+u3*g3;

          // Components of vector h = uu'g.
          h1i = u1*ug;
          h2i = u2*ug;
          h3i = u3*ug;
        }

        h1[i2][i1] = h1i;
        h2[i2][i1] = h2i;
        h3[i2][i1] = h3i;
      }
    }
  }

  /**
   * Computes vectors h = (1-eu)uu'g = g-Tg for one constant-i3 slice. 
   * The matrix T is defined by Schultz et al. (2010), except that the 
   * eigenvectors of that matrix are determined by fault strike and dip, 
   * not the eigen-decomposition of the Hessian. The eigenvalues of T 
   * are computed as in Schultz. Because two of those eigenvalues are 1, 
   * we have h = g-TG = (I-T)g = (1-eu)uu'g, where eu is the non-unit 
   * eigenvalue of T as in Schultz and u is the unit vector normal to 
   * the fault. As we search for ridges, we set h = 0 if (1) the 
   * smallest eigenvalue ew of the Hessian matrix H is non-negative or 
   * (2) the fault normal vector u is not aligned with the eigenvector 
   * w of H corresponding to the most negative eigenvalue ew.
   */
  private static void computeVectorsHSchultz(
    int i3, float flmin, float[][][] fs,
    float[][][] fl, float[][][] fp, float[][][] ft,
    float[][] h1, float[][] h2, float[][] h3)
  {
    int n1 = fs[0][0].length;
    int n2 = fs[0].length;

    // Arrays for eigen-decompositions of Hessians.
    double[][] a = new double[3][3]; // Hessian matrix
    double[][] z = new double[3][3]; // eigenvectors
    double[] e = new double[3]; // eigenvalues

    // Schultz's threshold theta for small differences ev-ew.
    float theta = 0.01f; // OK for F3D

    // For all samples in this constant-i3 slice, ...
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=1; i1<n1-1; ++i1) {

        // Components of vector h are initially zero.
        float h1i = h1[i2][i1] = 0.0f;
        float h2i = h2[i2][i1] = 0.0f;
        float h3i = h3[i2][i1] = 0.0f;

        // Skip samples with small fault likelihoods.
        if (fl[i3][i2][i1]<flmin) continue;

        // Hessian matrix H.
        float h12 = 0.25f*(fs[i3][i2+1][i1+1]-fs[i3][i2-1][i1+1] +
                           fs[i3][i2-1][i1-1]-fs[i3][i2+1][i1-1]);
        float h13 = 0.25f*(fs[i3+1][i2][i1+1]-fs[i3-1][i2][i1+1] +
                           fs[i3-1][i2][i1-1]-fs[i3+1][i2][i1-1]);
        float h23 = 0.25f*(fs[i3+1][i2+1][i1]-fs[i3-1][i2+1][i1] +
                           fs[i3-1][i2-1][i1]-fs[i3+1][i2-1][i1]);
        float h11 = fs[i3][i2][i1-1]-2.0f*fs[i3][i2][i1]+fs[i3][i2][i1+1];
        float h22 = fs[i3][i2-1][i1]-2.0f*fs[i3][i2][i1]+fs[i3][i2+1][i1];
        float h33 = fs[i3-1][i2][i1]-2.0f*fs[i3][i2][i1]+fs[i3+1][i2][i1];

        // Eigen-decomposition of H.
        a[0][0] = h11; a[0][1] = h12; a[0][2] = h13;
        a[1][0] = h12; a[1][1] = h22; a[1][2] = h23;
        a[2][0] = h13; a[2][1] = h23; a[2][2] = h33;
        Eigen.solveSymmetric33(a,z,e);
        float ev = (float)e[1];
        float ew = (float)e[2];

        // If we might have a ridge in fault likelihood, ...
        if (ew<0.0f) {

          // Normal vector w from fault likelihood.
          float w1 = (float)z[2][0];
          float w2 = (float)z[2][1];
          float w3 = (float)z[2][2];

          // Normal vector u from fault strike and dip.
          float pr = toRadians(fp[i3][i2][i1]);
          float tr = toRadians(ft[i3][i2][i1]);
          float cp = cos(pr);
          float sp = sin(pr);
          float ct = cos(tr);
          float st = sin(tr);
          float u1 = -st;
          float u2 = -sp*ct;
          float u3 =  cp*ct;

          // Ensure normal vectors are consistent. 
          // Bounds on u-w angle: 
          //   0.00f for < 90 degrees
          //   0.25f for < 60 degrees // OK for F3D
          //   0.50f for < 45 degrees
          //   0.75f for < 30 degrees
          float uw = u1*w1+u2*w2+u3*w3;
          if (uw*uw>0.00f) { 

            // The non-unit eigenvalue eu of T.
            float eu = 0.0f;
            if (ev-ew<theta) {
              eu = 1.0f-(ev-ew)/theta;
              eu = eu*eu;
            }

            // Gradient of fault likelihood.
            float g1 = 0.5f*(fs[i3][i2][i1+1]-fs[i3][i2][i1-1]);
            float g2 = 0.5f*(fs[i3][i2+1][i1]-fs[i3][i2-1][i1]);
            float g3 = 0.5f*(fs[i3+1][i2][i1]-fs[i3-1][i2][i1]);

            // Scaled dot product (1-eu)u'g.
            float ug = (1.0f-eu)*(u1*g1+u2*g2+u3*g3);

            // Components of vector h = (1-eu)uu'g.
            h1i = u1*ug;
            h2i = u2*ug;
            h3i = u3*ug;
          }
        }
        h1[i2][i1] = h1i;
        h2[i2][i1] = h2i;
        h3[i2][i1] = h3i;
      }
    }
  }

  /**
   * Processes one edge i---j of the image sampling grid.
   * If a fault intersects this edge, then this method constructs a 
   * new intersecting quad referencing four nodes within the fault.
   */
  private static void processEdge(
    int i1, int i2, int i3, int j1, int j2, int j3,
    float[][][] fl, float[][][] fp, float[][][] ft,
    float[][][] h1, float[][][] h2, float[][][] h3,
    Node[][][] nodes, ArrayList<Quad> quads)
  {
    // Avoid division by zero below.
    final float hsdtiny = 1.0e-6f;

    // If no edge i---j to process, simply return.
    if (i1==j1 && i2==j2 && i3==j3)
      return;

    // Vectors ei = xi-xj and ej = xj-xi.
    float e1i = i1-j1;
    float e2i = i2-j2;
    float e3i = i3-j3;
    float e1j = j1-i1;
    float e2j = j2-i2;
    float e3j = j3-i3;

    // Vectors h at endpoints of edge i---j.
    float h1i = h1[i3%2][i2][i1];
    float h2i = h2[i3%2][i2][i1];
    float h3i = h3[i3%2][i2][i1];
    float h1j = h1[j3%2][j2][j1];
    float h2j = h2[j3%2][j2][j1];
    float h3j = h3[j3%2][j2][j1];

    // If a fault intersects edge i---j, ...
    if (h1i*h1j+h2i*h2j+h3i*h3j<0.0f &&
        h1i*e1i+h2i*e2i+h3i*e3i<0.0f &&
        h1j*e1j+h2j*e2j+h3j*e3j<0.0f) {

      // Fault attributes for indices i and j.
      float fli = fl[i3][i2][i1];
      float fpi = fp[i3][i2][i1];
      float fti = ft[i3][i2][i1];
      float flj = fl[j3][j2][j1];
      float fpj = fp[j3][j2][j1];
      float ftj = ft[j3][j2][j1];
      fpi = toRadians(fpi); 
      fpj = toRadians(fpj); 
      float cpi = cos(fpi), spi = sin(fpi);
      float cpj = cos(fpj), spj = sin(fpj);
      float cci = cpi*cpi, ssi = spi*spi, csi = cpi*spi;
      float ccj = cpj*cpj, ssj = spj*spj, csj = cpj*spj;

      // Weights for interpolation of values for indices i and j.
      float h1d = h1j-h1i;
      float h2d = h2j-h2i;
      float h3d = h3j-h3i;
      float hsd = h1d*h1d+h2d*h2d+h3d*h3d;
      float wi = (hsd>hsdtiny)?(h1j*h1d+h2j*h2d+h3j*h3d)/hsd:0.5f;
      //if (wi>0.99f) wi = 0.99f;
      float wj = 1.0f-wi;

      // Compute values on edge of sampling grid via linear interpolation.
      float el = wi*fli+wj*flj;
      float et = wi*fti+wj*ftj;
      float cc = wi*cci+wj*ccj;
      float ss = wi*ssi+wj*ssj;
      float cs = wi*csi+wj*csj;
      float x1 = wi*i1+wj*j1;
      float x2 = wi*i2+wj*j2;
      float x3 = wi*i3+wj*j3;

      // The edge intersected and the four nodes for a new quad.
      int edge;
      Node na,nb,nc,nd;
      if (i1<j1) { // if edge i1---j1, ...
        edge = 1;
        na = nodeAt(i1,i2  ,i3  ,nodes);
        nb = nodeAt(i1,i2-1,i3  ,nodes);
        nc = nodeAt(i1,i2-1,i3-1,nodes);
        nd = nodeAt(i1,i2  ,i3-1,nodes);
      } else if (i2<j2) { // else if edge i2---j2, ...
        edge = 2;
        na = nodeAt(i1  ,i2,i3  ,nodes);
        nb = nodeAt(i1-1,i2,i3  ,nodes);
        nc = nodeAt(i1-1,i2,i3-1,nodes);
        nd = nodeAt(i1  ,i2,i3-1,nodes);
      } else { // else edge i3---j3, ...
        edge = 3;
        na = nodeAt(i1  ,i2  ,i3,nodes);
        nb = nodeAt(i1-1,i2  ,i3,nodes);
        nc = nodeAt(i1-1,i2-1,i3,nodes);
        nd = nodeAt(i1  ,i2-1,i3,nodes);
      }

      // Accumulate values in those four nodes.
      na.x1 += x1; nb.x1 += x1; nc.x1 += x1; nd.x1 += x1;
      na.x2 += x2; nb.x2 += x2; nc.x2 += x2; nd.x2 += x2;
      na.x3 += x3; nb.x3 += x3; nc.x3 += x3; nd.x3 += x3;
      na.fl += el; nb.fl += el; nc.fl += el; nd.fl += el;
      na.ft += et; nb.ft += et; nc.ft += et; nd.ft += et;
      na.u1 += cc; nb.u1 += cc; nc.u1 += cc; nd.u1 += cc;
      na.u2 += cs; nb.u2 += cs; nc.u2 += cs; nd.u2 += cs;
      na.u3 += ss; nb.u3 += ss; nc.u3 += ss; nd.u3 += ss;
      na.v1 +=  1; nb.v1 +=  1; nc.v1 +=  1; nd.v1 +=  1;

      // Construct a new quad that references the four nodes.
      Quad quad = new Quad(edge,i1,i2,i3,na,nb,nc,nd);
      quads.add(quad);
    }
  }

  /**
   * Processes one edge i---j of the image sampling grid.
   * If a fault intersects this edge, then this method constructs a 
   * new intersecting quad referencing four nodes within the fault.
   */
  private static void processEdgeSchultz(
    int i1, int i2, int i3, int j1, int j2, int j3,
    float[][][] fl, float[][][] fp, float[][][] ft,
    float[][][] h1, float[][][] h2, float[][][] h3,
    Node[][][] nodes, ArrayList<Quad> quads)
  {
    // Avoid division by zero below.
    final float hsdtiny = 1.0e-6f;

    // If no edge i---j to process, simply return.
    if (i1==j1 && i2==j2 && i3==j3)
      return;

    // Vectors h at endpoints of edge i---j.
    float h1i = h1[i3%2][i2][i1];
    float h2i = h2[i3%2][i2][i1];
    float h3i = h3[i3%2][i2][i1];
    float h1j = h1[j3%2][j2][j1];
    float h2j = h2[j3%2][j2][j1];
    float h3j = h3[j3%2][j2][j1];

    // If a fault intersects edge i---j, ...
    if (h1i*h1j+h2i*h2j+h3i*h3j<0.0f) {

      // Fault attributes for indices i and j.
      float fli = fl[i3][i2][i1];
      float fpi = fp[i3][i2][i1];
      float fti = ft[i3][i2][i1];
      float flj = fl[j3][j2][j1];
      float fpj = fp[j3][j2][j1];
      float ftj = ft[j3][j2][j1];
      fpi = toRadians(fpi); 
      fpj = toRadians(fpj); 
      float cpi = cos(fpi), spi = sin(fpi);
      float cpj = cos(fpj), spj = sin(fpj);
      float cci = cpi*cpi, ssi = spi*spi, csi = cpi*spi;
      float ccj = cpj*cpj, ssj = spj*spj, csj = cpj*spj;

      // Weights for interpolation of values for indices i and j.
      float h1d = h1j-h1i;
      float h2d = h2j-h2i;
      float h3d = h3j-h3i;
      float hsd = h1d*h1d+h2d*h2d+h3d*h3d;
      float wi = (hsd>hsdtiny)?(h1j*h1d+h2j*h2d+h3j*h3d)/hsd:0.5f;
      //if (wi>0.99f) wi = 0.99f;
      float wj = 1.0f-wi;

      // Compute values on edge of sampling grid via linear interpolation.
      float el = wi*fli+wj*flj;
      float et = wi*fti+wj*ftj;
      float cc = wi*cci+wj*ccj;
      float ss = wi*ssi+wj*ssj;
      float cs = wi*csi+wj*csj;
      float x1 = wi*i1+wj*j1;
      float x2 = wi*i2+wj*j2;
      float x3 = wi*i3+wj*j3;

      // The edge intersected and the four nodes for a new quad.
      int edge;
      Node na,nb,nc,nd;
      if (i1<j1) { // if edge i1---j1, ...
        edge = 1;
        na = nodeAt(i1,i2  ,i3  ,nodes);
        nb = nodeAt(i1,i2-1,i3  ,nodes);
        nc = nodeAt(i1,i2-1,i3-1,nodes);
        nd = nodeAt(i1,i2  ,i3-1,nodes);
      } else if (i2<j2) { // else if edge i2---j2, ...
        edge = 2;
        na = nodeAt(i1  ,i2,i3  ,nodes);
        nb = nodeAt(i1-1,i2,i3  ,nodes);
        nc = nodeAt(i1-1,i2,i3-1,nodes);
        nd = nodeAt(i1  ,i2,i3-1,nodes);
      } else { // else edge i3---j3, ...
        edge = 3;
        na = nodeAt(i1  ,i2  ,i3,nodes);
        nb = nodeAt(i1-1,i2  ,i3,nodes);
        nc = nodeAt(i1-1,i2-1,i3,nodes);
        nd = nodeAt(i1  ,i2-1,i3,nodes);
      }

      // Accumulate values in those four nodes.
      na.x1 += x1; nb.x1 += x1; nc.x1 += x1; nd.x1 += x1;
      na.x2 += x2; nb.x2 += x2; nc.x2 += x2; nd.x2 += x2;
      na.x3 += x3; nb.x3 += x3; nc.x3 += x3; nd.x3 += x3;
      na.fl += el; nb.fl += el; nc.fl += el; nd.fl += el;
      na.ft += et; nb.ft += et; nc.ft += et; nd.ft += et;
      na.u1 += cc; nb.u1 += cc; nc.u1 += cc; nd.u1 += cc;
      na.u2 += cs; nb.u2 += cs; nc.u2 += cs; nd.u2 += cs;
      na.u3 += ss; nb.u3 += ss; nc.u3 += ss; nd.u3 += ss;
      na.v1 +=  1; nb.v1 +=  1; nc.v1 +=  1; nd.v1 +=  1;

      // Construct a new quad that references the four nodes.
      Quad quad = new Quad(edge,i1,i2,i3,na,nb,nc,nd);
      quads.add(quad);
    }
  }
  private static Node nodeAt(int i1, int i2, int i3, Node[][][] nodes) {
    Node node = nodes[i3][i2][i1];
    if (node==null)
      nodes[i3][i2][i1] = node = new Node();
    return node;
  }

  /**
   * Completes the averaging of values for all non-null nodes.
   */
  private void completeNodes(Node[][][] nodes) {
    int n1 = nodes[0][0].length;
    int n2 = nodes[0].length;
    int n3 = nodes.length;

    // Averaging of fault strikes requires eigen-decomposition.
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];

    // Complete the averaging of values for all nodes.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          Node node = nodes[i3][i2][i1];
          if (node!=null) {
            a[0][0] =  node.u1; a[1][0] = -node.u2;
            a[0][1] = -node.u2; a[1][1] =  node.u3;
            Eigen.solveSymmetric22(a,z,e);
            float fp = toDegrees(atan2(z[1][0],z[1][1]));
            if (fp<-90.0f) fp += 180.0f;
            if (fp> 90.0f) fp -= 180.0f;
            node.fp = fp;
            float s = 1.0f/node.v1;
            node.fl *= s;
            node.ft *= s;
            node.x1 *= s;
            node.x2 *= s;
            node.x3 *= s;
            node.computeVectors();
          }
        }
      }
    }
  }

  /**
   * Returns a list of good quads, omitting any bad quads.
   * This cleaning must be performed before quads are linked.
   * All links to quad nabors are assumed to be null.
   */
  private ArrayList<Quad> cleanQuads(ArrayList<Quad> quads) {
    ArrayList<Quad> quadsGood = new ArrayList<Quad>(quads.size());
    int nqf = _quadFilters.size();
    for (Quad quad:quads) {
      boolean good = quad.isGood();
      for (int iqf=0; good && iqf<nqf; ++iqf) {
        QuadFilter qf = _quadFilters.get(iqf);
        if (!qf.good(quad))
            good = false;
      }
      if (good)
        quadsGood.add(quad);
    }
    return quadsGood;
  }
  public interface QuadFilter {
    boolean good(Quad quad);
  }
  public void addQuadFilter(QuadFilter qf) {
    _quadFilters.add(qf);
  }
  public void removeQuadFilter(QuadFilter qf) {
    _quadFilters.remove(qf);
  }
  private ArrayList<QuadFilter> _quadFilters = new ArrayList<QuadFilter>();

  /**
   * Removes any quads that have no nabors.
   */
  private static Quad[] removeSingleQuads(Quad[] quads) {
    ArrayList<Quad> quadList = new ArrayList<Quad>(quads.length);
    for (Quad q:quads) {
      if (q.qa!=null || q.qb!=null || q.qc!=null || q.qd!=null)
        quadList.add(q);
    }
    return quadList.toArray(new Quad[0]);
  }

  /**
   * Unlink any quads that are folded on top of one another.
   */
  private static void unlinkFoldedQuads(Quad[] quads) {
    for (Quad q:quads) {
      unlinkIfFolded(q,q.qa);
      unlinkIfFolded(q,q.qb);
      unlinkIfFolded(q,q.qc);
      unlinkIfFolded(q,q.qd);
    }
  }
  private static void unlinkIfFolded(Quad q1, Quad q2) {
    final float uusmall = 0.00f; // folded if angle > 90 degrees
    //final float uusmall = 0.25f; // folded if angle > 60 degrees
    //final float uusmall = 0.50f; // folded if angle > 45 degrees
    //final float uusmall = 0.75f; // folded if angle > 30 degrees
    if (q1!=null && q2!=null) {
      boolean qq = q1.isOrientedLikeNabor(q2);
      float uu = q1.u1*q2.u1+q1.u2*q2.u2+q1.u3*q2.u3;
      if (qq) {
        if (uu<uusmall) {
          q1.unlink();
          q2.unlink();
        }
      } else {
        if (uu>-uusmall) {
          q1.unlink();
          q2.unlink();
        }
      }
    }
  }

  /**
   * Recursively unlinks any quads with insufficient nabors.
   */
  private static void unlinkSkinnyQuads(Quad[] quads) {
    for (Quad q:quads)
      unlinkIfSkinny(q);
  }
  private static void unlinkIfSkinny(Quad q) {
    if (q!=null &&
        (q.qa==null && q.qc==null ||
         q.qb==null && q.qd==null)) {
      Quad qa = q.qa;
      Quad qb = q.qb;
      Quad qc = q.qc;
      Quad qd = q.qd;
      q.unlink();
      unlinkIfSkinny(qa);
      unlinkIfSkinny(qb);
      unlinkIfSkinny(qc);
      unlinkIfSkinny(qd);
    }
  }

  /**
   * Prints statistics for the specified array of quads.
   */
  private static void printStats(Quad[] quads) {
    int[] nq = new int[5];
    for (Quad q:quads) {
      int mq = 0;
      if (q.qa!=null) ++mq;
      if (q.qb!=null) ++mq;
      if (q.qc!=null) ++mq;
      if (q.qd!=null) ++mq;
      ++nq[mq];
    }
    trace("  quad stats:     number of quads = "+quads.length);
    trace("    number of quads with 0 nabors = "+nq[0]);
    trace("    number of quads with 1 nabor  = "+nq[1]);
    trace("    number of quads with 2 nabors = "+nq[2]);
    trace("    number of quads with 3 nabors = "+nq[3]);
    trace("    number of quads with 4 nabors = "+nq[4]);
  }

  private static void findShiftsMP(
    DynamicWarping dw, Quad[][] qew, Quad[][] qns) 
  {
    findShifts(dw,true,qew,qns);
  }

  private static void findShiftsPM(
    DynamicWarping dw, Quad[][] qew, Quad[][] qns) 
  {
    findShifts(dw,false,qew,qns);
  }

  /**
   * Extrapolates errors emp and epm where they could not be computed.
   * Errors that could not be computed are negative, corresponding to
   * errors for smaller lags for which errors could be computed. 
   * (Errors can always be computed for zero lag.) For each lag with 
   * a negative error, this method first attempts to extrapolate using 
   * other errors for the same lag stored in north or south quad nabors.
   * This first extrapolation works best when shifts vary slowly north
   * and south; it is ideal for constant shifts up or down fault dips.
   * <p> 
   * If this first extrapolation is impossible, because the number of
   * north-south quad nabors is too small for some lag, then errors are 
   * extrapolated using the errors already computed for other lags.
   * Those errors are already stored in the quads, but are negative, so
   * we simply change their sign.
   */
  private static void extrapolateErrors(Quad[][] qns) {
    int mns = qns.length;

    // For all arrays of quads linked north-south, ...
    for (int ins=0; ins<mns; ++ins) {
      int lns = qns[ins].length;
      float[][] emp = new float[lns][];
      float[][] epm = new float[lns][];

      // Get arrays of errors for all quads in the array.
      for (int jns=0; jns<lns; ++jns) {
        Vert v = qns[ins][jns].getVert();
        emp[jns] = v.emp;
        epm[jns] = v.epm;
      }
      int lmax = (emp[0].length-1)/2;

      // For each array of errors, ...
      for (int jns=0; jns<lns; ++jns) {

        // For all lags, negative and positive, ...
        for (int ilag=1; ilag<=lmax; ++ilag) {
          int ilagm = lmax-ilag;
          int ilagp = lmax+ilag;

          // Extrapolate for negative lags.
          float empim = emp[jns][ilagm];
          float epmim = epm[jns][ilagm];
          for (int kns=jns; kns<lns && empim<0.0f; ++kns)
            if (emp[kns][ilagm]>=0.0f)
              empim = emp[kns][ilagm];
          for (int kns=jns; kns<lns && epmim<0.0f; ++kns)
            if (epm[kns][ilagm]>=0.0f)
              epmim = epm[kns][ilagm];
          if (empim<0.0f) empim = -empim;
          if (epmim<0.0f) epmim = -epmim;
          emp[jns][ilagm] = empim;
          epm[jns][ilagm] = epmim;

          // Extrapolate for positive lags.
          float empip = emp[jns][ilagp];
          float epmip = epm[jns][ilagp];
          for (int kns=jns; kns>=0 && empip<0.0f; --kns)
            if (emp[kns][ilagp]>=0.0f)
              empip = emp[kns][ilagp];
          for (int kns=jns; kns>=0 && epmip<0.0f; --kns)
            if (epm[kns][ilagp]>=0.0f)
              epmip = epm[kns][ilagp];
          if (empip<0.0f) empip = -empip;
          if (epmip<0.0f) epmip = -epmip;
          emp[jns][ilagp] = empip;
          epm[jns][ilagp] = epmip;
        }
      }
    }
  }

  private static void findShifts(
    DynamicWarping dw, boolean mp, Quad[][] qew, Quad[][] qns) 
  {

    // Arrays of arrays of errors, linked east to west.
    int mew = qew.length;
    float[][][] eew = new float[mew][][];
    for (int iew=0; iew<mew; ++iew) {
      int lew = qew[iew].length;
      eew[iew] = new float[lew][];
      for (int jew=0; jew<lew; ++jew) {
        Vert v = qew[iew][jew].getVert();
        eew[iew][jew] = mp?v.emp:v.epm;
      }
    }

    // Arrays of arrays of errors, linked north to south.
    int mns = qns.length;
    float[][][] ens = new float[mns][][];
    for (int ins=0; ins<mns; ++ins) {
      int lns = qns[ins].length;
      ens[ins] = new float[lns][];
      for (int jns=0; jns<lns; ++jns) {
        Vert v = qns[ins][jns].getVert();
        ens[ins][jns] = mp?v.emp:v.epm;
      }
    }

    // Smooth the errors in north-south and east-west directions.
    dw.accumulate1(ens,ens);
    dw.accumulate1(eew,eew);
    dw.accumulate1(ens,ens);
    dw.accumulate1(eew,eew);

    // Find shifts by accumulating and backtracking.
    for (int ins=0; ins<mns; ++ins) {
      float[][] dns = dw.accumulateForward(ens[ins]);
      float[] s = dw.findShiftsReverse(dns,ens[ins]);
      int lns = s.length;
      for (int jns=0; jns<lns; ++jns) {
        Vert v = qns[ins][jns].getVert();
        if (mp) {
          v.smp += s[jns];
        } else {
          v.spm += s[jns];
        }
      }
    }
  }

  private static float[] like(float[] x) {
    return new float[x.length];
  }
  private static float[][] like(float[][] x) {
    int n = x.length; 
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = like(x[i]);
    return y;
  }
  private static float[][][] like(float[][][] x) {
    int n = x.length; 
    float[][][] y = new float[n][][];
    for (int i=0; i<n; ++i)
      y[i] = like(x[i]);
    return y;
  }

  private static void trace(String s) {
    System.out.println(s);
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
