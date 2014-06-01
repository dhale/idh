/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;
import interp.*;

// FOR DEVELOPMENT ONLY
import edu.mines.jtk.mosaic.*;

/**
 * Dynamic warping to find shifts between two sequences or images.
 * For 1D sequences f(x) and g(x), dynamic warping computes shifts u(x) such
 * that f(x) ~ g(x+u(x)). For 2D images f(x1,x2) and g(x1,x2), it finds shifts
 * u(x1,x2) such that f(x1,x2) ~ g(x1+u(x1,x2),x2). Note that the shifts
 * u(x1,x2,...) are computed for only the first dimension of multi-dimensional
 * images. For example, if the 1st dimension of an image is time, then only
 * time shifts are computed.
 * <p>
 * Constraints are placed on strains, the rates at which shifts change in any
 * dimension. For example, strain r1 = du/d1 is the derivative of shift
 * u(x1,x2,...) with respect to the x1 coordinate, and is constrained to lie
 * between lower and upper bounds r1min and r1max. Default bounds are 
 * r1min = -1.0 and r1max = 1.0.
 * <p>
 * In many applications of warping, strains derived from estimated shifts may
 * be an important by-product. However, when computing shifts, only a finite
 * number of strains are permitted, and this quantization may yield
 * unrealistic estimates of strain. To address this problem, this dynamic
 * warping provides control over the sampling of strains, in the form of a
 * smoothness value. The number of strains sampled is proportional to this
 * value.
 * <p>
 * Smoothness values represent approximate intervals for a quasi-uniform
 * subsampling grid on which shifts are computed. For example, for a time
 * sampling interval of 4 ms, one might specify a smoothness of 200 ms, so
 * that shifts are computed on a grid that is 50 times more coarse than the 4
 * ms sampling grid. However, strains on this coarse grid would be sampled 50
 * times more finely. After initially computing shifts on the coarse grid,
 * shifts are interpolated onto the finer grid.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.12.18
 */
public class DynamicWarpingR {

  /**
   * Constructs a dynamic warping.
   * If this warping is used for 2D or 3D images, then default unit samplings
   * are assumed for 2nd and 3rd dimensions.
   * @param smin lower bound on shift.
   * @param smax upper bound on shift.
   * @param s1 sampling of shifts for 1st dimension.
   */
  public DynamicWarpingR(double smin, double smax, Sampling s1) {
    this(smin,smax,s1,null,null);
  }

  /**
   * Constructs a dynamic warping.
   * If this warping is used for 3D images, then default unit samplings
   * are assumed for the 3rd dimensions.
   * @param smin lower bound on shift.
   * @param smax upper bound on shift.
   * @param s1 sampling of shifts for 1st dimension.
   * @param s2 sampling of shifts for 2nd dimension.
   */
  public DynamicWarpingR(double smin, double smax, Sampling s1, Sampling s2) {
    this(smin,smax,s1,s2,null);
  }

  /**
   * Constructs a dynamic warping.
   * @param smin lower bound on shift.
   * @param smax upper bound on shift.
   * @param s1 sampling of shifts for 1st dimension.
   * @param s2 sampling of shifts for 2nd dimension.
   * @param s3 sampling of shifts for 3rd dimension.
   */
  public DynamicWarpingR(
    double smin, double smax, 
    Sampling s1, Sampling s2, Sampling s3) 
  {
    double ds = s1.getDelta(); // shift sampling interval
    int ismin = (int) ceil(smin/ds);
    int ismax = (int)floor(smax/ds);
    _ss = new Sampling(1+ismax-ismin,ds,ismin*ds);
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _r1min = -1.0;
    _r2min = -1.0;
    _r3min = -1.0;
    _r1max =  1.0;
    _r2max =  1.0;
    _r3max =  1.0;
    _k1min = 10;
    _k2min = 10;
    _k3min = 10;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
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

    trace("findShifts: smoothing in 1st dimension ...");
    final float[][][] ek = new float[n2][][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e1 = computeErrors(sf,f[i2],sg,g[i2]);
      ek[i2] = subsampleErrors(_r1min,_r1max,k1s,ss,s1,e1);
    }});
    normalizeErrors(ek);

    trace("findShifts: smoothing in 2nd dimension ...");
    final float[][][] ekk = new float[nk2][nk1][ns];
    Parallel.loop(nk1,new Parallel.LoopInt() {
    public void compute(int ik1) {
      float[][] e2 = new float[n2][ns];
      for (int i2=0; i2<n2; ++i2)
        for (int is=0; is<ns; ++is)
          e2[i2][is] = ek[i2][ik1][is];
      float[][] ek2 = subsampleErrors(_r2min,_r2max,k2s,ss,s2,e2);
      for (int ik2=0; ik2<nk2; ++ik2)
        for (int is=0; is<ns; ++is)
          ekk[ik2][ik1][is] = ek2[ik2][is];
    }});
    normalizeErrors(ekk);

    //trace("findShifts: smoothing subsampled errors ...");
    //smoothSubsampledErrors(_r1min,_r1max,k1s,_r2min,_r2max,k2s,ss,s1,s2,ekk);
    //normalizeErrors(ekk);

    //SimpleFrame frame = new SimpleFrame();
    //frame.addImagePanels(pow(ekk,0.1f)).setPercentiles(2,98);

    trace("findShifts: finding shifts ...");
    float[][] ukk = new float[nk2][];
    for (int ik2=0; ik2<nk2; ++ik2) {
      ukk[ik2] = findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,ss,s1,ekk[ik2]);
    }
    trace("findShifts: interpolating shifts ...");
    //float[][] u = interpolateShifts(s1,s2,k1s,k2s,ukk);
    float[][] u = interpolateShiftsBl(s1,s2,k1s,k2s,ukk);
    //float[][] u = interpolateShiftsIg(s1,s2,f,k1s,k2s,ukk);
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
   * @param u array of shifts.
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

  /**
   * Returns uniformly sampled warped image h(x1,x2) = g(x1+u(x1,x2),x2).
   * @param sg sampling of the sequence g to be warped.
   * @param g array for the sequence g to be warped.
   * @param u array of shifts.
   * @return array for the warped sequence h.
   */
  public float[][] applyShifts(Sampling sg, float[][] g, float[][] u) {
    int n2 = g.length;
    float[][] h = new float[n2][];
    for (int i2=0; i2<n2; ++i2)
      h[i2] = applyShifts(sg,g[i2],u[i2]);
    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _ss,_s1,_s2,_s3;
  private double _r1min,_r2min,_r3min;
  private double _r1max,_r2max,_r3max;
  private int _k1min,_k2min,_k3min;
  private SincInterpolator _si;
  private float _epow = 2.0f;

  private static CubicInterpolator makeInterpolator1(
    float[] x, float[] y) 
  {
    //return new CubicInterpolator(CubicInterpolator.Method.LINEAR,x,y);
    return new CubicInterpolator(CubicInterpolator.Method.SPLINE,x,y);
  }
  private static CubicInterpolator makeInterpolator2(
    float[] x, float[] y) 
  {
    return new CubicInterpolator(CubicInterpolator.Method.SPLINE,x,y);
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
   * Smooths alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param me number of errors summed per smoothed error.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @return array[ne][ns] of subsampled errors.
   */
  private static float[][] smoothErrors(
    double rmin, double rmax, int me, 
    Sampling ss, Sampling se, float[][] e) 
  {
    int ns = ss.getCount();
    int ne = se.getCount();
    float[][] df = new float[ne][ns];
    float[][] dr = new float[ne][ns];
    accumulate( 1,rmin,rmax,me,ss,se,e,df);
    accumulate(-1,rmin,rmax,me,ss,se,e,dr);
    float[][] d = df;
    float scale = 1.0f/ne;
    for (int ie=0; ie<ne; ++ie) {
      for (int is=0; is<ns; ++is) {
        d[ie][is] = scale*(df[ie][is]+dr[ie][is]-e[ie][is]);
      }
    }
    return d;
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
      int msmin,msmax;
      if (me>0) {
        msmin = (int) ceil(-rmax*me*de/ds);
        msmax = (int)floor(-rmin*me*de/ds);
      } else {
        msmin = (int) ceil(-rmin*me*de/ds);
        msmax = (int)floor(-rmax*me*de/ds);
      }
      if (msmin>msmax) {
        trace("ie="+ie+" je="+je+" me="+me+" msmin="+msmin+" msmax="+msmax);
      }
      assert msmin<=msmax:"msmin<=msmax";
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
   * Accumulates alignment errors in forward or reverse direction.
   * Does not subsample the accumulated errors.
   * @param dir direction, 1 for forward, -1 for reverse.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param me number of errors e summed per accumulated error d
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @param d output array[ne][ns] of accumulated errors.
   */
  private static void accumulate(
    int dir, double rmin, double rmax, int me, 
    Sampling ss, Sampling se, float[][] e, float[][] d) 
  {
    int ns = ss.getCount();
    int ne = se.getCount();
    double ds = ss.getDelta();
    double de = se.getDelta();
    int ied = dir>0?1:-1;
    int ieb = dir>0?0:ne-1;
    int iee = dir>0?ne:-1;
    for (int is=0; is<ns; ++is)
      d[ieb][is] = e[ieb][is];
    float[] dprev = new float[ns];
    for (int ie=ieb+ied; ie!=iee; ie+=ied) {
      int je = max(0,min(ne-1,ie-dir*me));
      int ke = ie-je;
      int msmin,msmax;
      if (ke>0) {
        msmin = (int) ceil(-rmax*ke*de/ds);
        msmax = (int)floor(-rmin*ke*de/ds);
      } else {
        msmin = (int) ceil(-rmin*ke*de/ds);
        msmax = (int)floor(-rmax*ke*de/ds);
      }
      fill(Float.MAX_VALUE,d[ie]);
      for (int ms=msmin; ms<=msmax; ++ms) {
        int islo = max(0,-ms);
        int ishi = min(ns,ns-ms);
        for (int is=islo; is<ishi; ++is)
          dprev[is] = d[ie-ied][is+ms];
        updateSumsOfErrors(ie,je,ms,e,dprev,d[ie],null);
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
      int msmin,msmax;
      if (me>0) {
        msmin = (int) ceil(-rmax*me*de/ds);
        msmax = (int)floor(-rmin*me*de/ds);
      } else {
        msmin = (int) ceil(-rmin*me*de/ds);
        msmax = (int)floor(-rmax*me*de/ds);
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
    CubicInterpolator ci = makeInterpolator1(xk1,uk);
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
      CubicInterpolator ci = makeInterpolator1(xk1,ukk[jk2]);
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
      CubicInterpolator ciu = makeInterpolator2(xk2,uk2);
      CubicInterpolator civ = makeInterpolator2(xk2,vk2);
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
  private static float[][] interpolateShiftsBl(
    Sampling s1, Sampling s2, int[] k1s, int[] k2s, float[][] ukk) 
  {
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

    // Interpolate.
    BilinearInterpolator2 bl = new BilinearInterpolator2(xk1,xk2,ukk);
    //BicubicInterpolator2 bl = new BicubicInterpolator2(
    //  BicubicInterpolator2.Method.MONOTONIC,
    //  BicubicInterpolator2.Method.SPLINE,
    //  xk1,xk2,ukk);
    float[][] u = new float[n2][n1];
    for (int j2=0; j2<n2; ++j2) {
      float x2 = (float)s2.getValue(j2);
      for (int j1=0; j1<n1; ++j1) {
        float x1 = (float)s1.getValue(j1);
        u[j2][j1] = bl.interpolate(x1,x2);
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

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  private static void normalizeErrors(float[][][] e) {
    final int ns = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int is=0; is<ns; ++is) {
          float ei = ef[i2][i1][is];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }
  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
    trace("shiftAndScale: emin="+emin+" emax="+emax);
    final int ns = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int is=0; is<ns; ++is) {
          ef[i2][i1][is] = (ef[i2][i1][is]-eshift)*escale;
        }
      }
    }});
  }

  private static float[][] interpolateShiftsIg(
    Sampling s1, Sampling s2, float[][] f,
    int[] k1s, int[] k2s, float[][] ukk)
  {
    //trace("ukk:"); dump(ukk);
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;

    // Interpolate 1st derivatives of u at coarsely-sampled x2.
    float[] x1k = new float[nk1];
    for (int ik1=0; ik1<nk1; ++ik1)
      x1k[ik1] = (float)s1.getValue(k1s[ik1]);
    int nk = n1*nk2;
    float[] x1s = new float[nk];
    float[] x2s = new float[nk];
    float[] vs = new float[nk];
    for (int ik2=0,ik=0; ik2<nk2; ++ik2) {
      CubicInterpolator ci = 
        new CubicInterpolator(CubicInterpolator.Method.LINEAR,x1k,ukk[ik2]);
      for (int i1=0; i1<n1; ++i1,++ik) {
        x1s[ik] = (float)s1.getValue(i1);
        x2s[ik] = (float)s2.getValue(k2s[ik2]);
        vs[ik] = ci.interpolate1(x1s[ik]);
      }
    }

    // Tensor-guided interpolation of 1st derivatives.
    LocalOrientFilter lof = new LocalOrientFilter(8.0,2.0);
    EigenTensors2 et = lof.applyForTensors(f);
    et.invertStructure(0.0,8.0);
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.addPixels(f);
    TensorsView tv = new TensorsView(s1,s2,et);
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT);
    sp.getPlotPanel().addTiledView(tv);
    BlendedGridder2 bg = new BlendedGridder2(vs,x1s,x2s);
    bg.setSmoothness(0.5);
    float[][] v = bg.grid(s1,s2);

    // Integrate 1st derivatives to find shifts.
    float[][] u = new float[n2][n1];
    float d1 = (float)s1.getDelta();
    for (int i2=0; i2<n2; ++i2) {
      u[i2][0] = 0.0f;
      for (int i1=1; i1<n1; ++i1) {
        u[i2][i1] = u[i2][i1-1]+v[i2][i1-1]*d1;
      }
    }
    return u;
  }

  private static class Cm implements LeastSquaresSolver.CovarianceOperator {
    Cm(SmoothCovariance sc, Sampling s1, Sampling s2, Tensors2 t) {
      _sc = sc;
      _s1 = s1;
      _s2 = s2;
      _t = t;
    }
    public Vec apply(Vec x) {
      VecArrayFloat2 vx = (VecArrayFloat2)x;
      float[][] ay = copy(vx.getArray());
      _sc.apply(_s1,_s2,_t,ay);
      return new VecArrayFloat2(ay);
    }
    private SmoothCovariance _sc;
    private Sampling _s1,_s2;
    private Tensors2 _t;
  }
  private static class G implements LeastSquaresSolver.LinearOperator {
    G(Sampling s1, Sampling s2, int[] k1s, int[] k2s) {
      _s1 = s1;
      _s2 = s2;
      _k1s = k1s;
      _k2s = k2s;
      _nk1 = k1s.length;
      _nk2 = k2s.length;
      _n1 = s1.getCount();
      _n2 = s2.getCount();
      _nk = _nk1*_nk2;
    }
    public Vec apply(Vec x) {
      VecArrayFloat2 vx = (VecArrayFloat2)x;
      float[][] ax = vx.getArray();
      ax = integrateForward((float)_s1.getDelta(),ax);
      float[] ay = new float[_nk];
      for (int ik2=0,ik=0; ik2<_nk2; ++ik2) {
        int i2 = _k2s[ik2];
        for (int ik1=0; ik1<_nk1; ++ik1,++ik) {
          int i1 = _k1s[ik1];
          ay[ik] = ax[i2][i1];
        }
      }
      return new VecArrayFloat1(ay);
    }
    public Vec applyTranspose(Vec x) {
      VecArrayFloat1 vx = (VecArrayFloat1)x;
      float[] ax = vx.getArray();
      float[][] ay = new float[_n2][_n1];
      for (int ik2=0,ik=0; ik2<_nk2; ++ik2) {
        int i2 = _k2s[ik2];
        for (int ik1=0; ik1<_nk1; ++ik1,++ik) {
          int i1 = _k1s[ik1];
          ay[i2][i1] = ax[ik];
        }
      }
      ay = integrateReverse((float)_s1.getDelta(),ay);
      return new VecArrayFloat2(ay);
    }
    private float[][] integrateForward(final float d1, final float[][] x) {
      final float[][] y = new float[_n2][_n1];
      Parallel.loop(_n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[] x2 = x[i2];
        float[] y2 = y[i2];
        y2[0] = x2[0]*d1;
        for (int i1=1; i1<_n1; ++i1)
          y2[i1] = y2[i1-1]+d1*x2[i1];
      }});
      return y;
    }
    private float[][] integrateReverse(final float d1, final float[][] x) {
      final float[][] y = new float[_n2][_n1];
      Parallel.loop(_n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[] x2 = x[i2];
        float[] y2 = y[i2];
        y2[_n1-1] = x2[_n1-1]*d1;
        for (int i1=_n1-2; i1>=0; --i1)
          y2[i1] = y2[i1+1]+d1*x2[i1];
      }});
      return y;
    }
    private Sampling _s1,_s2;
    private int[] _k1s,_k2s;
    private int _nk1,_nk2,_n1,_n2,_nk;
  }

  private static float[][] interpolateShiftsLs(
    Sampling s1, Sampling s2, float[][] f,
    int[] k1s, int[] k2s, float[][] ukk)
  {
    //trace("ukk:"); dump(ukk);
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;

    LocalOrientFilter lof = new LocalOrientFilter(8.0);
    EigenTensors2 et = lof.applyForTensors(f);
    et.invertStructure(0.0,2.0);
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.addPixels(f);
    TensorsView tv = new TensorsView(s1,s2,et);
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT);
    sp.getPlotPanel().addTiledView(tv);
    double range = 100.0; // TODO: compute this!
    SmoothCovariance sc = new SmoothCovariance(1.0,1.0,range,2);
    Cm cm = new Cm(sc,s1,s2,et);
    G g = new G(s1,s2,k1s,k2s);
    LeastSquaresSolver lss = new LeastSquaresSolver();
    lss.setLinearOperator(g);
    lss.setModelCovariance(cm);
    int nk = nk1*nk2;
    float[] au = new float[nk];
    for (int ik2=0,ik=0; ik2<nk2; ++ik2) {
      for (int ik1=0; ik1<nk1; ++ik1,++ik) {
        au[ik] = ukk[ik2][ik1];
      }
    }
    VecArrayFloat1 vu = new VecArrayFloat1(au);
    VecArrayFloat2 vv = (VecArrayFloat2)lss.solve(vu);
    float[][] av = vv.getArray();
    float[][] u = new float[n2][n1];
    float d1 = (float)s1.getDelta();
    for (int i2=0; i2<n2; ++i2) {
      u[i2][0] = av[i2][0];
      for (int i1=1; i1<n1; ++i1) {
        u[i2][i1] = u[i2][i1-1]+av[i2][i1-1]*d1;
      }
    }
    return u;
  }
}
