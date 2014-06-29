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
import static fah.FaultUtil.*;

/**
 * Finds fault surfaces using fault likelihoods and orientations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.03
 */
public class FaultSurfer3X {

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
      for (FaultQuad quad:surf) {
        if (quad.isVertical()) {
          FaultQuad.Vert v = quad.getVert();
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
      for (FaultQuad quad:surf) {
        if (quad.isVertical()) {
          FaultQuad.Vert v = quad.getVert();

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
  // Surf
  ///////////////////////////////////////////////////////////////////////////

  /**
   * A fault surface consists of a set of linked and oriented quads.
   */
  public class Surf implements Iterable<FaultQuad> {
    Surf(FaultQuad[] quads) {
      _quads = quads;
    }
    public int size() {
      return _quads.length;
    }
    public Iterator<FaultQuad> iterator() {
      return Arrays.asList(_quads).iterator();
    }
    public float[][] getXyzUvwRgb() {
      return FaultSurfer3X.getXyzUvwRgb(_quads,_flmin);
    }
    public float[][] getXyzUvwRgbShifts(float smax) {
      return FaultSurfer3X.getXyzUvwRgbShifts(_quads,smax);
    }

    /**
     * Replaces all fault likelihoods with the average for the surface.
     */
    public void smooth() {
      float flsum = 0.0f;
      int nlsum = 0;
      for (FaultQuad q:_quads) {
        flsum += q.na.fl;
        flsum += q.nb.fl;
        flsum += q.nc.fl;
        flsum += q.nd.fl;
        nlsum += 4;
      }
      float flavg = flsum/nlsum;
      for (FaultQuad q:_quads) {
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
      for (FaultQuad q:_quads)
        u2sum += q.u2;
      if (u2sum<0.0) {
        for (FaultQuad q:_quads)
          q.flip();
        for (FaultQuad q:_quads)
          q.orientNodes();
      }
    }

    /**
     * Orients this surface. Orientation is computed such that the 
     * average of u1 components of quad normal vectors is negative.
     */
    public void orient() {
      float u1sum = 0.0f;
      for (FaultQuad q:_quads)
        u1sum += q.u1;
      if (u1sum>0.0) {
        for (FaultQuad q:_quads)
          q.flip();
        for (FaultQuad q:_quads)
          q.orientNodes();
      }
    }

    /**
     * Moves all nodes to make a topologically equivalent blocky surface.
     * A blocky surface is useful to see which adjacent image samples lie
     * on opposite sides of a fault.
     */
    public void blocky() {
      for (FaultQuad q:_quads)
        q.blocky();
    }

    public float[] sampleFaultDip() {
      float c1 = 0.0f;
      float c2 = 0.0f;
      float c3 = 0.0f;
      int nq = 0;
      for (FaultQuad quad:_quads) {
        c1 += quad.c1;
        c2 += quad.c2;
        c3 += quad.c3;
        nq += 1;
      }
      c1 /= nq;
      c2 /= nq;
      c3 /= nq;
      float dsmin = Float.MAX_VALUE;
      FaultQuad qmin = null;
      for (FaultQuad quad:_quads) {
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
      FaultQuad.Vert vert = qmin.getVert();
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
      FaultQuad[][][] qewns = getQuadsEwNs();
      FaultQuad[][] qew = qewns[0];
      FaultQuad[][] qns = qewns[1];
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
      for (FaultQuad quad:_quads) {
        if (quad.isVertical()) {
          FaultQuad.Vert vert = quad.getVert();
          vert.computeErrorsAndInitShifts(lmax,d,f,p2,p3);
        }
      }
    }
    private void smoothShifts() {
      for (FaultQuad q:this) {
        if (q.isVertical()) {
          FaultQuad.Vert v = q.getVert();
          float smp = 0.0f;
          float spm = 0.0f;
          float css = 0.0f;
          for (FaultQuad qn:new FaultQuad[]{v.qs,v.qn,v.qw,v.qe}) {
            if (qn!=null) {
              FaultQuad.Vert vn = qn.getVert();
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
      for (FaultQuad q:this) {
        if (q.isVertical()) {
          FaultQuad.Vert v = q.getVert();
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
      for (FaultQuad quad:_quads) {
        if (quad.isVertical()) {
          FaultQuad.Vert vert = quad.getVert();
          vert.computeThrowFromShifts();
        }
      }
    }

    /**
     * Returns east-west and north-south arrays of quads.
     */
    private FaultQuad[][][] getQuadsEwNs() {
      HashSet<FaultQuad> qewSet = new HashSet<FaultQuad>(size());
      HashSet<FaultQuad> qnsSet = new HashSet<FaultQuad>(size());
      ArrayList<FaultQuad[]> qewList = new ArrayList<FaultQuad[]>();
      ArrayList<FaultQuad[]> qnsList = new ArrayList<FaultQuad[]>();
      for (FaultQuad q:this) {
        if (q.isVertical()) {
          if (!qewSet.contains(q)) {
            FaultQuad qi = q;
            for (FaultQuad qe=qi.getVert().qe; qe!=null; qe=qi.getVert().qe)
              qi = qe;
            ArrayList<FaultQuad> ql = new ArrayList<FaultQuad>();
            for (; qi!=null; qi=qi.getVert().qw) {
              ql.add(qi);
              qewSet.add(qi);
            }
            qewList.add(ql.toArray(new FaultQuad[0]));
          }
          if (!qnsSet.contains(q)) {
            FaultQuad qi = q;
            for (FaultQuad qn=qi.getVert().qn; qn!=null; qn=qi.getVert().qn)
              qi = qn;
            ArrayList<FaultQuad> ql = new ArrayList<FaultQuad>();
            for (; qi!=null; qi=qi.getVert().qs) {
              ql.add(qi);
              qnsSet.add(qi);
            }
            qnsList.add(ql.toArray(new FaultQuad[0]));
          }
        }
      }
      FaultQuad[][] qew = qewList.toArray(new FaultQuad[0][]);
      FaultQuad[][] qns = qnsList.toArray(new FaultQuad[0][]);
      return new FaultQuad[][][]{qew,qns};
    }

    private FaultQuad[] _quads;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  /**
   * Constructs a fault surfer for specified likelihoods and orientations.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   */
  public FaultSurfer3X(float[][][][] flpt) {
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
  public FaultQuad[] findQuads() {
    int n1 = _n1, n2 = _n2, n3 = _n3;
    float[][][] fl = _fl, fp = _fp, ft = _ft;
    float flmin = _flmin;

    // Directional derivatives of fault likelihood. Quads correspond to edges
    // of the image sampling grid where the 1st derivative changes sign and
    // the 2nd derivative is negative.
    float[][][][] fd = fd(new float[][][][]{fl,fp,ft});
    float[][][] f1 = fd[0];
    float[][][] f2 = fd[1];

    // Array of nodes; will be non-null where faults are found.
    FaultNode[][][] nodes = new FaultNode[n3][n2][n1];

    // List of quads within fault surfaces.
    ArrayList<FaultQuad> quads = new ArrayList<FaultQuad>();

    // Construct nodes by looking for intersections of quads with grid edges.
    // We assume that edges within slice i3 have already been checked.
    // Therefore, we need only check edges within the slice j3, and edges
    // between the two slices i3 and j3.
    for (int i3=1,j3=i3+1; i3<n3; ++i3,++j3) {
      j3 = min(j3,n3-1);
      for (int i2=1,j2=i2+1; i2<n2; ++i2,++j2) {
        j2 = min(j2,n2-1);
        for (int i1=1,j1=i1+1; i1<n1; ++i1,++j1) {
          j1 = min(j1,n1-1);

          // Edge i1---j1.
          if (fl[i3][i2][i1]>=flmin && 
              fl[i3][i2][j1]>=flmin &&
              f1[i3][i2][i1]*f1[i3][i2][j1]<0.0f &&
              f2[i3][i2][i1]<0.0f && f2[i3][i2][j1]<0.0f)
            processEdge(i1,i2,i3,j1,i2,i3,fl,fp,ft,fd,nodes,quads);

          // Edge i2---j2.
          if (fl[i3][i2][i1]>=flmin && 
              fl[i3][j2][i1]>=flmin &&
              f1[i3][i2][i1]*f1[i3][j2][i1]<0.0f &&
              f2[i3][i2][i1]<0.0f && f2[i3][j2][i1]<0.0f)
            processEdge(i1,i2,i3,i1,j2,i3,fl,fp,ft,fd,nodes,quads);

          // Edge i3---j3.
          if (fl[i3][i2][i1]>=flmin && 
              fl[j3][i2][i1]>=flmin &&
              f1[i3][i2][i1]*f1[j3][i2][i1]<0.0f &&
              f2[i3][i2][i1]<0.0f && f2[j3][i2][i1]<0.0f)
            processEdge(i1,i2,i3,i1,i2,j3,fl,fp,ft,fd,nodes,quads);
        }
      }
    }

    // Complete computation of values for all non-null nodes. 
    completeNodes(nodes);

    return quads.toArray(new FaultQuad[0]);
  }
  private static void processEdge(
    int i1, int i2, int i3, int j1, int j2, int j3,
    float[][][] fl, float[][][] fp, float[][][] ft, float[][][][] fd,
    FaultNode[][][] nodes, ArrayList<FaultQuad> quads)
  {
    float[][][] f1 = fd[0]; // 1st derivative
    float[][][] f2 = fd[1]; // 2nd derivative

    // Fault attributes for grid indices i and j.
    float fli = fl[i3][i2][i1];
    float fpi = fp[i3][i2][i1];
    float fti = ft[i3][i2][i1];
    float f1i = f1[i3][i2][i1];
    float flj = fl[j3][j2][j1];
    float fpj = fp[j3][j2][j1];
    float ftj = ft[j3][j2][j1];
    float f1j = f1[j3][j2][j1];
    float[] ui = faultNormalFromStrikeAndDip(fpi,fti);
    float[] uj = faultNormalFromStrikeAndDip(fpj,ftj);
    float u1i = ui[0];
    float u2i = ui[1];
    float u3i = ui[2];
    float u1j = uj[0];
    float u2j = uj[1];
    float u3j = uj[2];

    // Weights for interpolation of values for indices i and j.
    float wi = -f1j/(f1i-f1j);
    float wj = 1.0f-wi;

    // Compute values on edge of sampling grid via linear interpolation.
    float el = wi*fli+wj*flj;
    float x1 = wi*i1+wj*j1;
    float x2 = wi*i2+wj*j2;
    float x3 = wi*i3+wj*j3;
    float u1 = wi*u1i+wj*u1j;
    float u2 = wi*u2i+wj*u2j;
    float u3 = wi*u3i+wj*u3j;

    // The edge intersected and the four nodes for a new quad.
    int edge;
    FaultNode na,nb,nc,nd;
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
    na.accumulate(el,x1,x2,x3,u1,u2,u3);
    nb.accumulate(el,x1,x2,x3,u1,u2,u3);
    nc.accumulate(el,x1,x2,x3,u1,u2,u3);
    nd.accumulate(el,x1,x2,x3,u1,u2,u3);

    // Construct a new quad that references the four nodes.
    FaultQuad quad = new FaultQuad(edge,i1,i2,i3,na,nb,nc,nd);
    quads.add(quad);
  }

  private static class TwoQuads {
    FaultQuad q1,q2;
    TwoQuads(FaultQuad q) {
      q1 = q;
    }
    void linkOrUnlink(Edge e, FaultQuad q) {
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
  public FaultQuad[] linkQuads(FaultQuad[] quads) {

    // Link two quads that share an edge. However, if an edge is shared
    // by more than two quads, do not link any of them. This latter case
    // should be rare, so we first link quads and then later unlink them
    // if necessary.
    HashMap<Edge,TwoQuads> map = new HashMap<Edge,TwoQuads>(4*quads.length);
    for (FaultQuad q:quads) {
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
  public void shuffleQuads(FaultQuad[] quads) {
    Random r = new Random(3);
    int n = quads.length;
    for (int i=n-1; i>0; --i) {
      int j = r.nextInt(i+1);
      FaultQuad quadi = quads[i];
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
  public Surf[] findSurfs(FaultQuad[] quads) {

    // List of surfs found while orienting quads.
    ArrayList<Surf> surfs = new ArrayList<Surf>();

    // Set of quads that have been oriented.
    HashSet<FaultQuad> set = new HashSet<FaultQuad>(quads.length);

    // Stack of quads that have been oriented, but with nabors to visit.
    ArrayDeque<FaultQuad> stack = new ArrayDeque<FaultQuad>();

    // Count surfaces that are not orientable.
    int nsno = 0;

    // For all quads, ...
    for (FaultQuad quad:quads) {

      // If quad has not yet been oriented, ...
      if (!set.contains(quad)) {

        // Begin collecting quads for a new surface.
        ArrayList<FaultQuad> list = new ArrayList<FaultQuad>();

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
          FaultQuad q = stack.pop();

          // Orient the nodes for the oriented quad q.
          q.orientNodes();

          // For all quad nabors qn of the oriented quad q, ...
          FaultQuad[] qns = new FaultQuad[]{q.qa,q.qb,q.qc,q.qd};
          for (FaultQuad qn:qns) {
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
        Surf surf = new Surf(list.toArray(new FaultQuad[0]));
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
  public static float[][] getXyzUvwRgb(FaultQuad[] quads, float flmin) {
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList fcl = new FloatList();
    for (FaultQuad quad:quads) {
      FaultNode na = quad.na;
      FaultNode nb = quad.nb;
      FaultNode nc = quad.nc;
      FaultNode nd = quad.nd;
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
  private static void acc(
      HashMap<FaultNode,Sum> map, FaultNode node, float v) {
    Sum s = map.get(node);
    if (s==null)
      map.put(node,s=new Sum());
    s.add(v);
  }
  private static float avg(HashMap<FaultNode,Sum> map, FaultNode node) {
    Sum s = map.get(node);
    return (s!=null)?s.avg():0.0f;
  }
  public static float[][] getXyzUvwRgbShifts(FaultQuad[] quads, float smax) {
    HashMap<FaultNode,Sum> map = new HashMap<FaultNode,Sum>(2*quads.length);
    for (FaultQuad quad:quads) {
      if (quad.isVertical()) {
        FaultQuad.Vert vert = quad.getVert();
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
    for (FaultQuad quad:quads) {
      FaultNode na = quad.na; float sa = avg(map,na);
      FaultNode nb = quad.nb; float sb = avg(map,nb);
      FaultNode nc = quad.nc; float sc = avg(map,nc);
      FaultNode nd = quad.nd; float sd = avg(map,nd);
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
    FaultNode na,nb;
    Edge(FaultNode na, FaultNode nb) {
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
    void link(FaultQuad q1, FaultQuad q2) {
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

  private static FaultNode nodeAt(
      int i1, int i2, int i3, FaultNode[][][] nodes) {
    FaultNode node = nodes[i3][i2][i1];
    if (node==null)
      nodes[i3][i2][i1] = node = new FaultNode();
    return node;
  }

  /**
   * Completes the averaging of values for all non-null nodes.
   */
  private void completeNodes(FaultNode[][][] nodes) {
    int n1 = nodes[0][0].length;
    int n2 = nodes[0].length;
    int n3 = nodes.length;

    // Complete the averaging of values for all nodes.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          FaultNode node = nodes[i3][i2][i1];
          if (node!=null)
            node.complete();
        }
      }
    }
  }

  /**
   * Removes any quads that have no nabors.
   */
  private static FaultQuad[] removeSingleQuads(FaultQuad[] quads) {
    ArrayList<FaultQuad> quadList = new ArrayList<FaultQuad>(quads.length);
    for (FaultQuad q:quads) {
      if (q.qa!=null || q.qb!=null || q.qc!=null || q.qd!=null)
        quadList.add(q);
    }
    return quadList.toArray(new FaultQuad[0]);
  }

  /**
   * Unlink any quads that are folded on top of one another.
   */
  private static void unlinkFoldedQuads(FaultQuad[] quads) {
    for (FaultQuad q:quads) {
      unlinkIfFolded(q,q.qa);
      unlinkIfFolded(q,q.qb);
      unlinkIfFolded(q,q.qc);
      unlinkIfFolded(q,q.qd);
    }
  }
  private static void unlinkIfFolded(FaultQuad q1, FaultQuad q2) {
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
  private static void unlinkSkinnyQuads(FaultQuad[] quads) {
    for (FaultQuad q:quads)
      unlinkIfSkinny(q);
  }
  private static void unlinkIfSkinny(FaultQuad q) {
    if (q!=null &&
        (q.qa==null && q.qc==null ||
         q.qb==null && q.qd==null)) {
      FaultQuad qa = q.qa;
      FaultQuad qb = q.qb;
      FaultQuad qc = q.qc;
      FaultQuad qd = q.qd;
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
  private static void printStats(FaultQuad[] quads) {
    int[] nq = new int[5];
    for (FaultQuad q:quads) {
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
    DynamicWarping dw, FaultQuad[][] qew, FaultQuad[][] qns) 
  {
    findShifts(dw,true,qew,qns);
  }

  private static void findShiftsPM(
    DynamicWarping dw, FaultQuad[][] qew, FaultQuad[][] qns) 
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
  private static void extrapolateErrors(FaultQuad[][] qns) {
    int mns = qns.length;

    // For all arrays of quads linked north-south, ...
    for (int ins=0; ins<mns; ++ins) {
      int lns = qns[ins].length;
      float[][] emp = new float[lns][];
      float[][] epm = new float[lns][];

      // Get arrays of errors for all quads in the array.
      for (int jns=0; jns<lns; ++jns) {
        FaultQuad.Vert v = qns[ins][jns].getVert();
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
    DynamicWarping dw, boolean mp, FaultQuad[][] qew, FaultQuad[][] qns) 
  {

    // Arrays of arrays of errors, linked east to west.
    int mew = qew.length;
    float[][][] eew = new float[mew][][];
    for (int iew=0; iew<mew; ++iew) {
      int lew = qew[iew].length;
      eew[iew] = new float[lew][];
      for (int jew=0; jew<lew; ++jew) {
        FaultQuad.Vert v = qew[iew][jew].getVert();
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
        FaultQuad.Vert v = qns[ins][jns].getVert();
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
        FaultQuad.Vert v = qns[ins][jns].getVert();
        if (mp) {
          v.smp += s[jns];
        } else {
          v.spm += s[jns];
        }
      }
    }
  }

  /**
   * Computes 1st and 2nd directional derivatives of fault likelihood.
   * The direction is normal to the fault, computed using fault strike and
   * dip. The normal vector points upwards, from the footwall into the hanging
   * wall.
   */
  public static float[][][][] fd(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    float[][][] fs = new float[n3][n2][n1];
    rgf.apply000(f,fs);
    float[][][] d1 = new float[n3][n2][n1];
    float[][][] d2 = new float[n3][n2][n1];
    int inset23 = 5;
    for (int i3=inset23; i3<n3-inset23; ++i3) {
      int i3m = i3-1, i3p = i3+1;
      for (int i2=inset23; i2<n2-inset23; ++i2) {
        int i2m = i2-1, i2p = i2+1;
        for (int i1=1; i1<n1-1; ++i1) {
          int i1m = i1-1, i1p = i1+1;

          // Need 19 samples.
          float fmmi = fs[i3m][i2m][i1 ];
          float fmim = fs[i3m][i2 ][i1m];
          float fmii = fs[i3m][i2 ][i1 ];
          float fmip = fs[i3m][i2 ][i1p];
          float fmpi = fs[i3m][i2p][i1 ];
          float fimm = fs[i3 ][i2m][i1m];
          float fimi = fs[i3 ][i2m][i1 ];
          float fimp = fs[i3 ][i2m][i1p];
          float fiim = fs[i3 ][i2 ][i1m];
          float fiii = fs[i3 ][i2 ][i1 ];
          float fiip = fs[i3 ][i2 ][i1p];
          float fipm = fs[i3 ][i2p][i1m];
          float fipi = fs[i3 ][i2p][i1 ];
          float fipp = fs[i3 ][i2p][i1p];
          float fpmi = fs[i3p][i2m][i1 ];
          float fpim = fs[i3p][i2 ][i1m];
          float fpii = fs[i3p][i2 ][i1 ];
          float fpip = fs[i3p][i2 ][i1p];
          float fppi = fs[i3p][i2p][i1 ];

          // Gradient vector.
          float g1 = 0.5f*(fiip-fiim);
          float g2 = 0.5f*(fipi-fimi);
          float g3 = 0.5f*(fpii-fmii);

          // Hessian matrix.
          float h11 = fiip-2.0f*fiii+fiim;
          float h12 = 0.25f*(fipp-fipm-fimp+fimm);
          float h13 = 0.25f*(fpip-fpim-fmip+fmim);
          float h21 = h12;
          float h22 = fipi-2.0f*fiii+fimi;
          float h23 = 0.25f*(fppi-fpmi-fmpi+fmmi);
          float h31 = h13;
          float h32 = h23;
          float h33 = fpii-2.0f*fiii+fmii;

          // Fault normal vector, pointing upward.
          float pi = toRadians(p[i3][i2][i1]);
          float ti = toRadians(t[i3][i2][i1]);
          float cp = cos(pi); 
          float sp = sin(pi);
          float ct = cos(ti);
          float st = sin(ti);
          float u1 = -st;
          float u2 = -sp*ct;
          float u3 =  cp*ct;
          if (u1>0.0) {
            u1 = -u1;
            u2 = -u2;
            u3 = -u3;
          }

          // 1st and 2nd directional derivatives.
          d1[i3][i2][i1] = u1*g1+u2*g2+u3*g3;
          d2[i3][i2][i1] = 
              u1*h11*u1+u1*h12*u2+u1*h13*u3+
              u2*h21*u1+u2*h22*u2+u2*h23*u3+
              u3*h31*u1+u3*h32*u2+u3*h33*u3;
        }
      }
    }
    return new float[][][][]{d1,d2};
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
