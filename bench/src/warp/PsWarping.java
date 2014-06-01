/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
//package edu.mines.jtk.dsp;
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping for PP and PS seismic images.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.23
 */
public class PsWarping {

  /**
   * Constructs a warping for a pair of PP and PS traces.
   * @param gamin lower bound on average gamma = average Vp/Vs.
   * @param gamax upper bound on average gamma = average Vp/Vs.
   * @param spp time sampling for PP trace.
   * @param sps time sampling for PS trace.
   */
  public PsWarping(double gamin, double gamax, Sampling spp, Sampling sps) {
    this(gamin,gamax,spp,sps,null,null);
  }

  /**
   * Constructs a warping for a pair of 2D PP and PS images.
   * @param gamin lower bound on average gamma = average Vp/Vs.
   * @param gamax upper bound on average gamma = average Vp/Vs.
   * @param spp time sampling for PP image.
   * @param sps time sampling for PS image.
   * @param s2 sampling for 2nd image dimension.
   */
  public PsWarping(
    double gamin, double gamax, 
    Sampling spp, Sampling sps, Sampling s2) 
  {
    this(gamin,gamax,spp,sps,s2,null);
  }

  /**
   * Constructs a warping for a pair of 3D PP and PS images.
   * @param gamin lower bound on average gamma = average Vp/Vs.
   * @param gamax upper bound on average gamma = average Vp/Vs.
   * @param spp time sampling for PP image.
   * @param sps time sampling for PS image.
   * @param s2 sampling for 2nd image dimension.
   * @param s3 sampling for 3rd image dimension.
   */
  public PsWarping(
    double gamin, double gamax, 
    Sampling spp, Sampling sps, Sampling s2, Sampling s3) 
  {
    double c = 0.5*(1.0+gamin);
    double df = spp.getDelta();
    double ff = spp.getFirst();
    double lf = spp.getLast();
    int nf = spp.getCount();
    double dg = sps.getDelta()/c;
    double fg = sps.getFirst()/c;
    double lg = sps.getLast()/c;
    int ng = sps.getCount();
    double d1 = df;
    double f1 = ff;
    double l1 = lg*(1.0+gamin)/(1.0+gamax);
    int n1 = 1+(int)((l1-f1)/d1);
    double ds = d1;
    double fs = 0.0;
    double ls = lg-l1;
    int ns = 1+(int)((ls-fs)/ds);
    _sf = new Sampling(nf,df,ff);
    _sg = new Sampling(ng,dg,fg);
    _ss = new Sampling(ns,ds,fs);
    _s1 = new Sampling(n1,d1,f1);
    _s2 = s2;
    _s3 = s3;
    //_dw = new DynamicWarpingR(_ss.getFirst(),_ss.getLast(),_s1,_s2,_s3);
  }

  /**
   * Gets the sampling of shifts s used to compute alignment errors.
   * @return the sampling.
   */
  public Sampling getSamplingS() {
    return _ss;
  }

  /**
   * Gets the sampling in the 1st dimension.
   * @return the sampling.
   */
  public Sampling getSampling1() {
    return _s1;
  }

  /**
   * Gets the sampling in the 2nd dimension.
   * @return the sampling.
   */
  public Sampling getSampling2() {
    return _s2;
  }

  /**
   * Gets the sampling in the 3rd dimension.
   * @return the sampling.
   */
  public Sampling getSampling3() {
    return _s3;
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in 1st dimension.
   * @param r1max upper bound on strain in 1st dimension.
   */
  public void setStrainLimits(double r1min, double r1max) {
    setStrainLimits(r1min,r1max,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in 1st dimension.
   * @param r1max upper bound on strain in 1st dimension.
   * @param r2min lower bound on strain in 2nd dimension.
   * @param r2max upper bound on strain in 2nd dimension.
   */
  public void setStrainLimits(
    double r1min, double r1max,
    double r2min, double r2max)
  {
    setStrainLimits(r1min,r1max,r2min,r2max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in 1st dimension.
   * @param r1max upper bound on strain in 1st dimension.
   * @param r2min lower bound on strain in 2nd dimension.
   * @param r2max upper bound on strain in 2nd dimension.
   * @param r3min lower bound on strain in 3rd dimension.
   * @param r3max upper bound on strain in 3rd dimension.
   */
  public void setStrainLimits(
    double r1min, double r1max,
    double r2min, double r2max,
    double r3min, double r3max)
  {
    _r1min = r1min; _r1max = r1max;
    _r2min = r2min; _r2max = r2max;
    _r3min = r3min; _r3max = r3max;
  }

  /**
   * Sets the smoothness of shifts computed by this dynamic warping.
   * <em>
   * Units of smoothness are the same as those used in the corresponding
   * sampling of the 1st dimension.
   * </em>
   * Default smoothness values are 10 times the sampling intervals for
   * corresponding dimensions.
   * @param d1min smoothness in 1st dimension.
   */
  public void setSmoothness(double d1min) {
    double d2 = (_s2!=null)?_s2.getDelta():1.0;
    setSmoothness(d1min,10.0*d2);
  }

  /**
   * Sets the smoothness of shifts computed by this dynamic warping.
   * <em>
   * Units of smoothness are the same as those used in corresponding
   * samplings of 1st and 2nd dimensions.
   * </em>
   * Default smoothness values are 10 times the sampling intervals for
   * corresponding dimensions.
   * @param d1min smoothness in 1st dimension.
   * @param d2min smoothness in 2nd dimension.
   */
  public void setSmoothness(double d1min, double d2min) {
    double d3 = (_s3!=null)?_s3.getDelta():1.0;
    setSmoothness(d1min,d2min,10.0*d3);
  }

  /**
   * Sets the smoothness of shifts computed by this dynamic warping.
   * <em>
   * Units of smoothness are the same as those used in corresponding
   * samplings of 1st, 2nd and 3rd dimensions.
   * </em>
   * Default smoothness values are 10 times the sampling intervals for
   * corresponding dimensions.
   * @param d1min smoothness in 1st dimension.
   * @param d2min smoothness in 2nd dimension.
   * @param d3min smoothness in 3rd dimension.
   */
  public void setSmoothness(double d1min, double d2min, double d3min) {
    double d1 = (_s1!=null)?_s1.getDelta():1.0;
    double d2 = (_s2!=null)?_s2.getDelta():1.0;
    double d3 = (_s3!=null)?_s3.getDelta():1.0;
    _k1min = max(1,(int)ceil(d1min/d1));
    _k2min = max(1,(int)ceil(d2min/d2));
    _k3min = max(1,(int)ceil(d3min/d3));
  }

  /**
   * Returns shifts computed for specified 1D sequences.
   * @param sf sampling of 1st dimension for the seqeunce f.
   * @param f array of values for sequence f.
   * @param sg sampling of 1st dimension for the seqeunce g.
   * @param g array of values for sequence g.
   * @return array of shifts.
   */
  public float[] findShifts(
    Sampling sf, float[] f,
    Sampling sg, float[] g)
  {
    float[][] e = computeErrors(sf,f,sg,g);
    return findShifts(e);
  }

  /**
   * Returns shifts computed for specified 2D images.
   * @param sf sampling of 1st dimension for the image f.
   * @param f array of values for image f.
   * @param sg sampling of 1st dimension for the image g.
   * @param g array of values for image g.
   * @return array of shifts.
   */
  public float[][] findShifts(
    final Sampling sf, final float[][] f,
    final Sampling sg, final float[][] g)
  {
    // Samplings.
    final Sampling ss = _ss;
    final Sampling s1 = _s1;
    final Sampling s2 = sampling2(f);
    final int ns = ss.getCount();
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();

    // Quasi-uniform subsamplings.
    final int[] k1s = subsample(n1,_k1min);
    final int[] k2s = subsample(n2,_k2min);
    final int nk1 = k1s.length;
    final int nk2 = k2s.length;
    trace("k1s:"); dump(k1s);
    trace("k2s:"); dump(k2s);

    // ek[n2][nk1][ns] = errors subsampled in 1st dimension
    trace("findShifts: loop over i2 ...");
    final float[][][] ek = new float[n2][][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      // e1[n1][ns] are errors for 1st dimension
      float[][] e1 = computeErrors(sf,f[i2],sg,g[i2]);
      ek[i2] = subsampleErrors(_r1min,_r1max,k1s,ss,s1,e1);
    }});

    trace("findShifts: loop over ik1 ...");
    // ekk[nk2][nk1][ns] = errors subsampled in both 1st and 2nd dimensions
    final float[][][] ekk = new float[nk2][nk1][ns];
    Parallel.loop(nk1,new Parallel.LoopInt() {
    public void compute(int ik1) {

      // e2[n2][ns] = errors for 2nd dimension
      float[][] e2 = new float[n2][ns];
      for (int i2=0; i2<n2; ++i2) {
        for (int is=0; is<ns; ++is) {
          e2[i2][is] = ek[i2][ik1][is];
        }
      }

      // ek2[nk2][ns] = errors subsampled in 2nd dimension
      float[][] ek2 = subsampleErrors(_r2min,_r2max,k2s,ss,s2,e2);

      // copy to ekk[nk2][nk1][ns]
      for (int ik2=0; ik2<nk2; ++ik2) {
        for (int is=0; is<ns; ++is) {
          ekk[ik2][ik1][is] = ek2[ik2][is];
        }
      }
    }});

    trace("findShifts: smoothing subsampled errors ...");
    smoothSubsampledErrors(_r1min,_r1max,k1s,_r2min,_r2max,k2s,ss,s1,s2,ekk);

    trace("findShifts: finding shifts ...");
    float[][] ukk = new float[nk2][];
    for (int ik2=0; ik2<nk2; ++ik2) {
      ukk[ik2] = findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,ss,s1,ekk[ik2]);
    }
    trace("findShifts: interpolating shifts ...");
    float[][] u = interpolateShifts(s1,s2,k1s,k2s,ukk);
    trace("findShifts: ... done");
    return u;
  }

  /**
   * Returns alignment errors computed for specified sequences.
   * @param sf sampling of 1st dimension for the seqeunce f.
   * @param f array of values for sequence f.
   * @param sg sampling of 1st dimension for the seqeunce g.
   * @param g array of values for sequence g.
   * @return array of alignment errors.
   */
  public float[][] computeErrors(
    Sampling sf, float[] f,
    Sampling sg, float[] g)
  {
    Sampling ss = _ss;
    Sampling se = _s1;
    int ns = ss.getCount();
    int ne = se.getCount();
    int nf = sf.getCount();
    int ng = sg.getCount();
    float[][] e = new float[ne][ns];
    float[] fi = new float[ne];
    float[] gi = new float[ne];
    _si.interpolate(sf,f,se,fi);
    for (int is=0; is<ns; ++is) {
      _si.interpolate(
        ng,sg.getDelta(),sg.getFirst(),g,
        ne,se.getDelta(),se.getFirst()+ss.getValue(is),gi);
      for (int ie=0; ie<ne; ++ie)
        e[ie][is] = error(fi[ie],gi[ie]);
    }
    return e;
  }

  /**
   * Normalizes alignment errors, in place.
   * After normalizing, minimum error = 0 and maximum error = 1.
   */
  public static void normalizeErrors(float[][] e) {
    int ns = e[0].length;
    int ne = e.length;
    float emin = e[0][0];
    float emax = e[0][0];
    for (int ie=0; ie<ne; ++ie) {
      for (int is=0; is<ns; ++is) {
        float ei = e[ie][is];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Returns shifts estimated for specified alignment errors.
   * @param e array[n1][ns] of alignment errors.
   * @return array[n1] of shifts.
   */
  public float[] findShifts(float[][] e) {
    int ns = _ss.getCount();
    int n1 = _s1.getCount();
    double ds = _ss.getDelta();
    double d1 = _s1.getDelta();
    int k1min = min(_k1min,n1-1);
    int[] i1k = subsample(n1,k1min);
    //dump(i1k);
    int n1k = i1k.length;
    return findShiftsFromErrors(_r1min,_r1max,i1k,_ss,_s1,e);
  }

  /**
   * Returns uniformly sampled warped sequence h(x1) = g(x1+u(x1)).
   * @param sg sampling of the sequence g to be warped.
   * @param g array for the sequence g to be warped.
   * @return array for the warped sequence h.
   */
  public float[] applyShifts(Sampling sg, float[] g, float[] u) {
    Sampling s1 = _s1;
    int ng = sg.getCount();
    int n1 = s1.getCount();
    float[] h = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      double x1 = s1.getValue(i1)+u[i1];
      h[i1] = _si.interpolate(sg,g,x1);
    }
    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _sf,_sg,_ss,_s1,_s2,_s3;
  private double _r1min,_r2min,_r3min;
  private double _r1max,_r2max,_r3max;
  private int _k1min,_k2min,_k3min;
  private SincInterpolator _si;
  private float _epow = 2.0f;

  private static CubicInterpolator makeCubicInterpolator(
    float[] x, float[] y) 
  {
    return new CubicInterpolator(CubicInterpolator.Method.LINEAR,x,y);
  }
  private static CubicInterpolator makeCubicInterpolator(
    float[] x, float[] y, float[] yd) 
  {
    return new CubicInterpolator(x,y,yd);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  /**
   * Returns an approximately uniformly-sampled subset of indices in [0,n).
   * Indices in the subset are chosen to be approximately uniform, with the
   * difference between consecutive indices not less than the specified
   * minimum increment kmin. Because the first and last indices 0 and n-1 are
   * included in the subset, n must be greater than the minimum increment
   * kmin.
   * @param n number of indices in the set {0,1,2,...,n-1}.
   * @param kmin minimum increment between indices in the subset.
   * @return array of indices in the subset.
   */
  private static int[] subsample(int n, int kmin) {
    if (kmin>=n)
      kmin = n-1;
    int m = 1+(n-1)/kmin;
    double d = (double)(n-1)/(double)(m-1);
    int[] j = new int[m];
    for (int i=0; i<m; ++i)
      j[i] = (int)(i*d+0.5);
    return j;
  }
  private static void subsampleTest() {
    int kmin = 3;
    for (int n=kmin+1; n<5*kmin+1; ++n) {
      int[] j = subsample(n,kmin);
      System.out.println("n="+n);
      dump(j);
    }
  }
  private static void main(String[] args) {
    subsampleTest();
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }
  }

  /**
   * Subsamples alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @return array[nke][ns] of subsampled errors.
   */
  private static float[][] subsampleErrors(
    double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e) 
  {
    int ns = ss.getCount();
    int ne = se.getCount();
    int nke = kes.length;
    float[][] df = new float[nke][ns];
    float[][] dr = new float[nke][ns];
    accumulate( 1,rmin,rmax,kes,ss,se,e,df,null);
    accumulate(-1,rmin,rmax,kes,ss,se,e,dr,null);
    float[][] d = df;
    float scale = 1.0f/ne;
    for (int ike=0; ike<nke; ++ike) {
      int ke = kes[ike];
      for (int is=0; is<ns; ++is) {
        d[ike][is] = scale*(df[ike][is]+dr[ike][is]-e[ke][is]);
      }
    }
    return d;
  }

  private static void smoothSubsampledErrors(
    double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e) 
  {
    int ns = ss.getCount();
    int nke = kes.length;
    float[][] ef = new float[nke][ns];
    float[][] er = new float[nke][ns];
    accumulateSubsampled( 1,rmin,rmax,kes,ss,se,e,ef,null);
    accumulateSubsampled(-1,rmin,rmax,kes,ss,se,e,er,null);
    float scale = 1.0f/nke;
    for (int ike=0; ike<nke; ++ike) {
      int ke = kes[ike];
      for (int is=0; is<ns; ++is) {
        e[ike][is] = scale*(ef[ike][is]+er[ike][is]-e[ike][is]);
      }
    }
  }

  private static void smoothSubsampledErrors(
    double r1min, double r1max, int[] k1s,
    double r2min, double r2max, int[] k2s,
    Sampling ss, Sampling s1, Sampling s2, float[][][] e) 
  {
    int ns = ss.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;
    float[][] e2 = new float[nk2][];
    for (int ik1=0; ik1<nk1; ++ik1) {
      for (int ik2=0; ik2<nk2; ++ik2)
        e2[ik2] = e[ik2][ik1];
      smoothSubsampledErrors(r2min,r2max,k2s,ss,s2,e2);
      for (int ik2=0; ik2<nk2; ++ik2)
        e[ik2][ik1] = e2[ik2];
    }
    for (int ik2=0; ik2<nk2; ++ik2) {
      smoothSubsampledErrors(r1min,r1max,k1s,ss,s1,e[ik2]);
    }
  }

  /**
   * Finds shifts from alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @return array[ne] of shifts.
   */
  private static float[] findShiftsFromErrors(
    double rmin, double rmax, int[] kes,
    Sampling ss, Sampling se, float[][] e) 
  {
    int nke = kes.length;
    int ns = ss.getCount();
    float[][] d = new float[nke][ns];
    int[][] m = new int[nke][ns];
    accumulate(1,rmin,rmax,kes,ss,se,e,d,m);
    float[] uke = backtrackForShifts(kes,ss,se,d[nke-1],m);
    return interpolateShifts(se,kes,uke);
  }

  /**
   * Finds shifts from subsampled alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e array[nke][ns] of subsampled alignment errors.
   * @return array[ne] of shifts.
   */
  private static float[] findShiftsFromSubsampledErrors(
    double rmin, double rmax, int[] kes,
    Sampling ss, Sampling se, float[][] e) 
  {
    int nke = kes.length;
    int ns = ss.getCount();
    float[][] d = new float[nke][ns];
    int[][] m = new int[nke][ns];
    accumulateSubsampled(1,rmin,rmax,kes,ss,se,e,d,m);
    return backtrackForShifts(kes,ss,se,d[nke-1],m);
  }

  /**
   * Accumulates alignment errors in forward or reverse direction.
   * @param dir direction, 1 for forward, -1 for reverse.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @param d output array[nke][ns] of accumulated errors.
   * @param m output array[nke][ns] of minimizing moves; or null.
   */
  private static void accumulate(
    int dir, double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e, float[][] d, int[][] m) 
  {
    int ns = ss.getCount();
    double ds = ss.getDelta();
    double de = se.getDelta();
    int nke = kes.length;
    int iked = dir>0?1:-1;
    int ikeb = dir>0?0:nke-1;
    int ikee = dir>0?nke:-1;
    for (int is=0; is<ns; ++is)
      d[ikeb][is] = e[kes[ikeb]][is];
    float[] dprev = new float[ns];
    for (int ike=ikeb+iked; ike!=ikee; ike+=iked) {
      int je = kes[ike-iked];
      int ie = kes[ike];
      int me = ie-je;
      int msmin = (int) ceil(-rmax*me*de/ds);
      int msmax = (int)floor(-rmin*me*de/ds);
      if (msmin>msmax) {
        int mstmp = msmin;
        msmin = msmax;
        msmax = mstmp;
      }
      fill(Float.MAX_VALUE,d[ike]);
      for (int ms=msmin; ms<=msmax; ++ms) {
        int islo = max(0,-ms);
        int ishi = min(ns,ns-ms);
        for (int is=islo; is<ishi; ++is)
          dprev[is] = d[ike-iked][is+ms];
        if (m!=null)
          updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],m[ike]);
        else
          updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],null);
      }
    }
  }

  /**
   * Updates sums of errors for one shift between two error sample indices.
   * Computes sums of errors along linear trajectories according to:
   * <pre>
   * d[is] += sum from ke=ie to ke!=je of e[ke][is+(ie-ke)*ms/(ie-je)]
   * </pre>
   * Where the complicated last subscript (typically) is not an integer, this
   * method uses linear interpolation of the alignment errors e. After the
   * sums of errors have been computed for all shift indices is, this method
   * updates the minimum sum of errors and the corresponding change in shift.
   * @param ie error sample index at which to begin sum.
   * @param je error sample index at which to end (not) sum.
   * @param ms change in shift at error sample index je, not in sum.
   * @param e[ne][ns] input array of alignment errors.
   * @param d[ns] input/output array in which to accumulate errors.
   * @param dmin[ns] input/output array of minimum accumulated errors.
   * @param mmin[ns] input/output array of minimizing moves; or null.
   */
  private static void updateSumsOfErrors(
    int ie, int je, int ms, float[][] e, float[] d, float[] dmin, int[] mmin) 
  {
    int ns = d.length;
    int islo = max(0,-ms); // update only for is >= islo
    int ishi = min(ns,ns-ms); // update only for is < ishi
    for (int is=islo; is<ishi; ++is)
      d[is] += e[ie][is];
    int me = ie-je;
    int de = me>0?-1:1;
    if (ms==0) { // if no shift, no interpolation required
      for (int ke=ie+de; ke!=je; ke+=de) {
        for (int is=islo; is<ishi; ++is) {
          d[is] += e[ke][is];
        }
      }
    } else { // else, use linear interpolation of errors
      float r = (float)ms/(float)me; // strain
      for (int ke=ie+de; ke!=je; ke+=de) {
        float sk = r*(ie-ke);
        int ks = (int)sk;
        if (sk<0.0f) --ks;
        int ksa = ks+1;
        int ksb = ks;
        float wsa = sk-ks;
        float wsb = 1.0f-wsa;
        for (int is=islo; is<ishi; ++is)
          d[is] += wsa*e[ke][is+ksa]+wsb*e[ke][is+ksb];
      }
    }
    for (int is=islo; is<ishi; ++is) {
      if (d[is]<dmin[is]) {
        dmin[is] = d[is];
        if (mmin!=null)
          mmin[is] = ms;
      }
    }
  }

  /**
   * Accumulates subsampled errors in forward or reverse direction.
   * @param dir direction, 1 for forward, -1 for reverse.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[nke][ns] of subsampled errors.
   * @param d output array[nke][ns] of accumulated errors.
   * @param m output array[nke][ns] of minimizing moves; or null.
   */
  private static void accumulateSubsampled(
    int dir, double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e, float[][] d, int[][] m) 
  {
    int ns = ss.getCount();
    double ds = ss.getDelta();
    double de = se.getDelta();
    int nke = kes.length;
    int iked = dir>0?1:-1;
    int ikeb = dir>0?0:nke-1;
    int ikee = dir>0?nke:-1;
    for (int is=0; is<ns; ++is)
      d[ikeb][is] = e[ikeb][is];
    for (int ike=ikeb+iked; ike!=ikee; ike+=iked) {
      float[] dprev = d[ike-iked];
      int me = kes[ike]-kes[ike-iked];
      int msmin = (int) ceil(-rmax*me*de/ds);
      int msmax = (int)floor(-rmin*me*de/ds);
      if (msmin>msmax) {
        int mstmp = msmin;
        msmin = msmax;
        msmax = mstmp;
      }
      for (int is=0; is<ns; ++is) {
        float dmin = Float.MAX_VALUE;
        int mmin = -1;
        for (int ms=msmin; ms<=msmax; ++ms) {
          int js = is+ms;
          if (0<=js && js<ns) {
            float dj = dprev[js];
            if (dj<dmin) {
              dmin = dj;
              mmin = ms;
            }
          }
          d[ike][is] = dmin+e[ike][is];
          if (m!=null)
            m[ike][is] = mmin;
        }
      }
    }
  }

  /**
   * Returns shifts found by backtracking with precomputed moves.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param d array[ns] of last forward accumulated errors.
   * @param m array[nke][ns] of moves, changes in shifts.
   * @return array[ne] of shifts.
   */
  private static float[] backtrackForShifts(
    int[] kes, Sampling ss, Sampling se, float[] d, int[][] m) 
  {
    int nke = kes.length;
    int ns = ss.getCount();
    int ne = se.getCount();
    int ike = nke-1;
    float dmin = Float.MAX_VALUE;
    int imin = -1;
    for (int is=0; is<ns; ++is) {
      if (d[is]<dmin) {
        dmin = d[is];
        imin = is;
      }
    }
    int is = imin;
    float[] uke = new float[nke];
    uke[ike] = (float)ss.getValue(is);
    for (--ike; ike>=0; --ike) {
      is += m[ike+1][is];
      uke[ike] = (float)ss.getValue(is);
    }
    return uke;
  }

  private static float[] interpolateShifts(
    Sampling s1, int[] k1s, float[] uk) 
  {
    int n1 = s1.getCount();
    int nk1 = k1s.length;
    float[] xk1 = new float[nk1];
    for (int jk1=0; jk1<nk1; ++jk1)
      xk1[jk1] = (float)s1.getValue(k1s[jk1]);
    CubicInterpolator ci = makeCubicInterpolator(xk1,uk);
    float[] u = new float[n1];
    for (int j1=0; j1<n1; ++j1) {
      float x1 = (float)s1.getValue(j1);
      u[j1] = ci.interpolate(x1);
    }
    return u;
  }
  private static float[][] interpolateShifts(
    Sampling s1, Sampling s2, int[] k1s, int[] k2s, float[][] ukk) 
  {
    //trace("ukk:"); dump(ukk);
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;

    // Coarse sampling of 1st and 2nd dimensions.
    float[] xk1 = new float[nk1];
    for (int jk1=0; jk1<nk1; ++jk1)
      xk1[jk1] = (float)s1.getValue(k1s[jk1]);
    float[] xk2 = new float[nk2];
    for (int jk2=0; jk2<nk2; ++jk2)
      xk2[jk2] = (float)s2.getValue(k2s[jk2]);

    // Compute 1st derivatives in 1st dimension.
    float[][] vkk = new float[nk2][nk1];
    for (int jk2=0; jk2<nk2; ++jk2) {
      CubicInterpolator ci = makeCubicInterpolator(xk1,ukk[jk2]);
      ci.interpolate1(xk1,vkk[jk2]);
    }
      
    // Interpolate in 2nd dimension.
    float[] uk2 = new float[nk2];
    float[] vk2 = new float[nk2];
    float[][] uk = new float[n2][nk1];
    float[][] vk = new float[n2][nk1];
    for (int jk1=0; jk1<nk1; ++jk1) {
      for (int jk2=0; jk2<nk2; ++jk2) {
        uk2[jk2] = ukk[jk2][jk1];
        vk2[jk2] = vkk[jk2][jk1];
      }
      CubicInterpolator ciu = makeCubicInterpolator(xk2,uk2);
      CubicInterpolator civ = makeCubicInterpolator(xk2,vk2);
      for (int j2=0; j2<n2; ++j2) {
        float x2 = (float)s2.getValue(j2);
        uk[j2][jk1] = ciu.interpolate(x2);
        vk[j2][jk1] = civ.interpolate(x2);
      }
    }

    // Interpolate 1st dimension.
    float[][] u = new float[n2][n1];
    for (int j2=0; j2<n2; ++j2) {
      CubicInterpolator ci = makeCubicInterpolator(xk1,uk[j2],vk[j2]);
      for (int j1=0; j1<n1; ++j1) {
        float x1 = (float)s1.getValue(j1);
        u[j2][j1] = ci.interpolate(x1);
      }
    }
    return u;
  }

  private Sampling sampling2(float[][] f) {
    if (_s2!=null) {
      Check.argument(f.length==_s2.getCount(),"valid sampling2");
      return _s2;
    } else {
      return new Sampling(f.length,1.0,0.0);
    }
  }
  private Sampling sampling2(float[][][] f) {
    if (_s2!=null) {
      Check.argument(f[0].length==_s2.getCount(),"valid sampling2");
      return _s2;
    } else {
      return new Sampling(f[0].length,1.0,0.0);
    }
  }
  private Sampling sampling3(float[][][] f) {
    if (_s3!=null) {
      Check.argument(f.length==_s3.getCount(),"valid sampling3");
      return _s3;
    } else {
      return new Sampling(f.length,1.0,0.0);
    }
  }
}
