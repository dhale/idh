/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.util.*;

import static edu.mines.jtk.util.ArrayMath.*;

import static fah.FaultGeometry.*;

/**
 * A fault surface is a collection of linked and oriented quads.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.29
 */
public class FaultSurf implements Iterable<FaultQuad> {

  /**
   * Constructs a fault surface from the specified fault quads.
   * @param quads array of fault quads.
   */
  FaultSurf(FaultQuad[] quads) {
    _quads = quads;
  }

  /**
   * Returns the number of quads in this fault surface.
   * @return the number of quads.
   */
  public int size() {
    return _quads.length;
  }

  /**
   * Returns an iterator for the quads in this fault surface.
   */
  public Iterator<FaultQuad> iterator() {
    return Arrays.asList(_quads).iterator();
  }

  /**
   * Orients this surface. Orientation is computed such that the 
   * average of w1 components of quad normal vectors is negative.
   */
  public void orient() {
    float w1sum = 0.0f;
    for (FaultQuad q:_quads)
      w1sum += q.w1;
    if (w1sum>0.0) {
      for (FaultQuad q:_quads)
        q.flip();
      for (FaultQuad q:_quads)
        q.orientNodes();
    }
  }

  /*
  public float[][] getXyzUvwRgb() {
    return FaultSurfer.getXyzUvwRgb(_quads,_flmin);
  }

  public float[][] getXyzUvwRgbShifts(float smax) {
    return FaultSurfer.getXyzUvwRgbShifts(_quads,smax);
  }
  */

  /**
   * Replaces all fault likelihoods with the average for this surface.
   */
  public void useAverageFaultLikelihood() {
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
   * Moves all nodes to make a topologically equivalent blocky surface.
   * A blocky surface is useful to see which adjacent image samples lie
   * on opposite sides of a fault.
   */
  public void blocky() {
    for (FaultQuad q:_quads)
      q.blocky();
  }

  /**
   * Returns xyz coordinates for up- and down-dip sampling of fault surface.
   * @return array of packed xyz coordinates of locations sampled.
   */
  public float[] sampleFaultDip() {

    // Average center coordinates for all quads in this surface.
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

    // Find the vertical quad nearest to that average.
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
   * @param smax the maximum shift.
   * @param f seismic image.
   * @param p2 slopes in the 2nd dimension.
   * @param p3 slopes in the 3rd dimension.
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
    zeroInconsistentShifts();
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
  private void zeroInconsistentShifts() {
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
   * Computes slip vectors for all vertical quads in this surface.
   * Uses shifts that must already have been computed for each 
   * vertical quad.
   */
  public void computeSlips() {
    for (FaultQuad quad:_quads) {
      if (quad.isVertical()) {
        FaultQuad.Vert vert = quad.getVert();
        vert.computeSlipFromShifts();
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
