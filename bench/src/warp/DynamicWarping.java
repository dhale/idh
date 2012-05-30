package warp;

import java.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping of sequences and images.
 * <p>
 * For sequences f and g, dynamic warping finds a sequence of 
 * integer shifts u such that f[i1] ~ g[i1+u[i1]], subject to a 
 * bound b1 on strain, the rate at which the shifts u[i1] vary 
 * with sample index i1.
 * <p>
 * An increasing u[i1] = u[i1-1] + 1 implies that g[i1] between 
 * indices i1-1 and i1 is a stretched version of f[i1] ~ g[i1+u[i1]].
 * Values in f for indices i1 and i1-1 are one sample apart, but
 * corresponding values in g are two samples apart, which implies 
 * stretching by 100%. Likewise, a decreasing u[i1] = u[i1-1] - 1 
 * implies squeezing by 100%.
 * <p>
 * In practice, 100% strain (stretching or squeezing) may be extreme.
 * Therefore, the upper bound on strain may be smaller than one. For 
 * example, if the bound b1 = 0.5, then the local average strain is 
 * bounded by 0.5. This constraint is complicated by the fact that 
 * the shifts u[i1] are integers. The actual constraint for the bound 
 * b1 = 0.5 is |u[i1-1]-u[i1-2]| + |u[i1]-u[i1-1]| &le; 1.
 * <p>
 * For 2D images f and g, dynamic warping finds a 2D array of integer 
 * shifts u[i2][i1] such that f[i2][i1] ~ g[i2][i1+u[i2][i1]], 
 * subject to bounds b1 and b2 on strains, the rates at which shifts 
 * u[i2][i1] vary with samples indices i1 and i2, respectively.
 * <p>
 * Estimated shifts u are integers, but can be smoothed to obtain
 * non-integer shifts. The extent of smoothing along each dimension 
 * is inversely proportional to the strain limit for that dimension, 
 * and these extents can be scaled by specified factors for more or 
 * less smoothing. The default scale factors are zero, for no 
 * smoothing.
 * <p>
 * Dynamic image warping may require a large amount of memory.
 * Specifically, a temporary array of nlag*nimage floats is required,
 * where the number of lags nlag = 1+shiftMax-shiftMin and nimage is 
 * the number of floats in the image. For 3D images, the product 
 * nlag*nimage is assumed to be too large for the array to fit in 
 * random-access memory (RAM), and this large array is instead stored 
 * in a temporary random-access file on disk.
 * <p>
 * This class provides numerous methods, but typical applications
 * require only several of these, including the methods that find
 * and apply shifts. The many other methods are provided only for 
 * research and atypical applications.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.05.25
 */
public class DynamicWarping {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public DynamicWarping(int shiftMin, int shiftMax) {
    Check.argument(shiftMax-shiftMin>1,"shiftMax-shiftMin>1");
    _lmin = shiftMin;
    _lmax = shiftMax;
    _nl = 1+_lmax-_lmin;
    _si = new SincInterpolator();
    _ref1 = _ref2 = _ref3 = new RecursiveExponentialFilter(0.0);
  }

  /**
   * Sets bound on strain for all dimensions. Must be in (0,1].
   * The actual bound on strain is 1.0/ceil(1.0/strainMax), which
   * is less than the specified strainMax when 1.0/strainMax is not
   * an integer. The default bound on strain is 1.0 (100%).
   * @param strainMax the bound, a value less than or equal to one.
   */
  public void setStrainMax(double strainMax) {
    Check.argument(strainMax<=1.0,"strainMax<=1.0");
    Check.argument(strainMax>0.0,"strainMax>0.0");
    setStrainMax(strainMax,strainMax);
  }

  /**
   * Sets bound on strains in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    setStrainMax(strainMax1,strainMax2,strainMax2);
  }

  /**
   * Sets bound on strains in 1st, 2nd and 3rd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   * @param strainMax3 bound on strain in the 3rd dimension.
   */
  public void setStrainMax(
    double strainMax1, double strainMax2, double strainMax3) 
  {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax3<=1.0,"strainMax3<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    Check.argument(strainMax3>0.0,"strainMax3>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
    _bstrain2 = (int)ceil(1.0/strainMax2);
    _bstrain3 = (int)ceil(1.0/strainMax3);
    _ref1 = new RecursiveExponentialFilter(_usmooth1*2.0*_bstrain1);
    _ref2 = new RecursiveExponentialFilter(_usmooth2*2.0*_bstrain2);
    _ref3 = new RecursiveExponentialFilter(_usmooth3*2.0*_bstrain3);
  }

  /**
   * Sets the exponent used to compute alignment errors |f-g|^e.
   * The default exponent is 2.
   * @param e the exponent.
   */
  public void setErrorExponent(double e) {
    _epow = (float)e;
  }

  /**
   * Sets the number of nonlinear smoothings of alignment errors.
   * In dynamic warping, alignment errors are smoothed the specified 
   * number of times, along all dimensions (in order 1, 2, ...), 
   * before estimating shifts by accumulating and backtracking along 
   * only the 1st dimension. The default number of smoothings is zero.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setErrorSmoothing(int esmooth) {
    _esmooth = esmooth;
  }

  /**
   * Sets extent of smoothing filters used to smooth integer shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factor. Default 
   * factor is zero, for no smoothing.
   * @param usmooth extent of smoothing filter in all dimensions.
   */
  public void setShiftSmoothing(double usmooth) {
    setShiftSmoothing(usmooth,usmooth);
  }

  /**
   * Sets extents of smoothing filters used to smooth integer shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setShiftSmoothing(double usmooth1, double usmooth2) {
    setShiftSmoothing(usmooth1,usmooth2,usmooth2);
  }

  /**
   * Sets extents of smoothing filters used to smooth integer shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   * @param usmooth3 extent of smoothing filter in 3rd dimension.
   */
  public void setShiftSmoothing(
    double usmooth1, double usmooth2, double usmooth3) 
  {
    _usmooth1 = usmooth1;
    _usmooth2 = usmooth2;
    _usmooth3 = usmooth3;
    _ref1 = new RecursiveExponentialFilter(usmooth1*2.0*_bstrain1);
    _ref2 = new RecursiveExponentialFilter(usmooth2*2.0*_bstrain2);
    _ref3 = new RecursiveExponentialFilter(usmooth3*2.0*_bstrain3);
  }

  /**
   * Sets the directory used for temporary random-access files.
   * Such files are used for arrays too large to fit in memory.
   * The specified directory should be within a partition with a
   * large amount of unused space available.
   * @param directory the directory used for temporary files.
   */
  public void setTempFileDirectory(String directory) {
    _edir = new File(directory);
  }

  /**
   * Sets the directory used for temporary random-access files.
   * Such files are used for arrays too large to fit in memory.
   * The specified directory should be within a partition with a
   * large amount of unused space available.
   * @param directory the directory used for temporary files.
   */
  public void setTempFileDirectory(File directory) {
    _edir = directory;
  }

  /**
   * Computes and returns shifts for specified sequences.
   * @param f array for the sequence f.
   * @param g array for the sequence g.
   * @return array of shifts u.
   */
  public float[] findShifts(float[] f, float[] g) {
    float[] u = like(f);
    findShifts(f,g,u);
    return u;
  }

  /**
   * Computes and returns shifts for specified images.
   * @param f array for the image f.
   * @param g array for the image g.
   * @return array of shifts u.
   */
  public float[][] findShifts(float[][] f, float[][] g) {
    float[][] u = like(f);
    findShifts(f,g,u);
    return u;
  }

  /**
   * Computes and returns shifts for specified images.
   * @param f array for the image f.
   * @param g array for the image g.
   * @return array of shifts u.
   */
  public float[][][] findShifts(float[][][] f, float[][][] g) {
    float[][][] u = like(f);
    findShifts(f,g,u);
    return u;
  }

  /**
   * Computes shifts for specified sequences.
   * @param f input array for the sequence f.
   * @param g input array for the sequence g.
   * @param u output array of shifts u.
   */
  public void findShifts(float[] f, float[] g, float[] u) {
    float[][] e = computeErrors(f,g);
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    float[][] d = accumulateForward(e);
    backtrackReverse(d,e,u);
    smoothShifts(u,u);
  }

  /**
   * Computes shifts for specified images.
   * @param f input array for the image f.
   * @param g input array for the image g.
   * @param u output array of shifts u.
   */
  public void findShifts(float[][] f, float[][] g, float[][] u) {
    final float[][][] e = computeErrors(f,g);
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      accumulateForward(e[i2],d);
      backtrackReverse(d,e[i2],uf[i2]);
    }});
    smoothShifts(u,u);
  }

  /**
   * Computes shifts for specified images.
   * @param f input array for the image f.
   * @param g input array for the image g.
   * @param u output array of shifts u.
   */
  public void findShifts(float[][][] f, float[][][] g, float[][][] u) {
    System.out.println("findShifts: begin");
    AlignmentErrors3 ae = computeErrors(f,g);
    System.out.println("findShifts: errors computed");
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(ae);
    System.out.println("findShifts: errors smoothed");
    computeShifts(ae,u);
    System.out.println("findShifts: shifts computed");
    smoothShifts(u);
    System.out.println("findShifts: shifts smoothed");
    ae.delete();
  }

  /**
   * Returns a sequence warped by applying specified shifts.
   * @param u array of shifts.
   * @param g array for the sequence to be warped.
   * @return array for the warped sequence.
   */
  public float[] applyShifts(float[] u, float[] g) {
    float[] h = like(g);
    applyShifts(u,g,h);
    return h;
  }

  /**
   * Returns an image warped by applying specified shifts.
   * @param u array of shifts.
   * @param g array for the image to be warped.
   * @return array for the warped image.
   */
  public float[][] applyShifts(float[][] u, float[][] g) {
    float[][] h = like(g);
    applyShifts(u,g,h);
    return h;
  }

  /**
   * Returns an image warped by applying specified shifts.
   * @param u array of shifts.
   * @param g array for the image to be warped.
   * @return array for the warped image.
   */
  public float[][][] applyShifts(float[][][] u, float[][][] g) {
    float[][][] h = like(g);
    applyShifts(u,g,h);
    return h;
  }

  /**
   * Computes a sequence warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the sequence to be warped.
   * @param h output array for the warped sequence.
   */
  public void applyShifts(float[] u, float[] g, float[] h) {
    int n1 = u.length;
    _si.setUniformSampling(n1,1.0,0.0);
    _si.setUniformSamples(g);
    for (int i1=0; i1<n1; ++i1) {
      h[i1] = _si.interpolate(i1+u[i1]);
    }
  }

  /**
   * Computes an image warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the image to be warped.
   * @param h output array for the warped image.
   */
  public void applyShifts(float[][] u, float[][] g, float[][] h) {
    final int n1 = u[0].length;
    final int n2 = u.length;
    final float[][] uf = u;
    final float[][] gf = g;
    final float[][] hf = h;
    final Parallel.Unsafe<SincInterpolator> siu =
      new Parallel.Unsafe<SincInterpolator>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      SincInterpolator si = siu.get();
      if (si==null) {
        si = new SincInterpolator();
        si.setUniformSampling(n1,1.0,0.0);
        siu.set(si);
      }
      si.setUniformSamples(gf[i2]);
      for (int i1=0; i1<n1; ++i1) {
        hf[i2][i1] = si.interpolate(i1+uf[i2][i1]);
      }
    }});
  }

  /**
   * Computes an image warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the image to be warped.
   * @param h output array for the warped image.
   */
  public void applyShifts(float[][][] u, float[][][] g, float[][][] h) {
    int n3 = u.length;
    final float[][][] uf = u;
    final float[][][] gf = g;
    final float[][][] hf = h;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyShifts(uf[i3],gf[i3],hf[i3]);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // for research and atypical applications

  /**
   * Returns normalized alignment errors for all samples and lags.
   * Alignment errors are normalized to be in range [0,1].
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax]. Alignment errors are a monotonically
   * increasing function of |f[i1]-g[i1+il+shiftMin]|.
   * @param f array[n1] for the sequence f[i1].
   * @param g array[n1] for the sequence g[i1].
   * @return array[n1][nl] of squared errors.
   */
  public float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    float[][] e = new float[n1][_nl];
    computeErrors(f,g,e);
    normalizeErrors(e);
    return e;
  }

  /**
   * Returns normalized alignment errors for all samples and lags.
   * Alignment errors are normalized to be in range [0,1].
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax]. Alignment errors are a monotonically
   * increasing function of |f[i2][i1]-g[i2][i1+il+shiftMin]|.
   * @param f array[n2][n1] for the image f[i2][i1].
   * @param g array[n2][n1] for the sequence g[i2][i1].
   * @return array[n2][n1][nl] of squared errors.
   */
  public float[][][] computeErrors(float[][] f, float[][] g) {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] ff = f;
    final float[][] gf = g;
    final float[][][] ef = new float[n2][n1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(ff[i2],gf[i2],ef[i2]);
    }});
    normalizeErrors(ef);
    return ef;
  }

  /**
   * Returns smoothed (and normalized) alignment errors.
   * @param e array[n1][nl] of alignment errors.
   * @return array[n1][nl] of smoothed errors.
   */
  public float[][] smoothErrors(float[][] e) {
    float[][] es = like(e);
    smoothErrors(e,es);
    return es;
  }

  /**
   * Returns smoothed (and normalized) alignment errors.
   * @param e array[n2][n1][nl] of alignment errors.
   * @return array[n2][n1][nl] of smoothed errors.
   */
  public float[][][] smoothErrors(float[][][] e) {
    float[][][] es = like(e);
    smoothErrors(e,es);
    return es;
  }

  /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n1][nl] of alignment errors.
   * @param es output array[n1][nl] of smoothed errors.
   */
  public void smoothErrors(float[][] e, float[][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
  }

  /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothErrors(float[][][] e, float[][][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
    smoothErrors2(_bstrain2,es,es);
    normalizeErrors(es);
  }

  /**
   * Returns smoothed shifts.
   * @param u array of shifts to be smoothed.
   * @return array of smoothed shifts
   */
  public float[] smoothShifts(float[] u) {
    float[] us = like(u);
    smoothShifts(u,us);
    return us;
  }

  /**
   * Returns smoothed shifts.
   * @param u array of shifts to be smoothed.
   * @return array of smoothed shifts
   */
  public float[][] smoothShifts(float[][] u) {
    float[][] us = like(u);
    smoothShifts(u,us);
    return us;
  }

  /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothShifts(float[] u, float[] us) {
    _ref1.apply(u,us);
  }

  /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothShifts(float[][] u, float[][] us) {
    _ref1.apply1(u,us);
    _ref2.apply2(us,us);
  }

  /**
   * Returns errors accumulated in forward direction.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    accumulateForward(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][] accumulateReverse(float[][] e) {
    float[][] d = like(e);
    accumulateReverse(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in forward direction in 1st dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateForward1(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward1(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction in 1st dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateReverse1(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse1(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in forward direction in 2nd dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateForward2(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward2(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction in 2nd dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateReverse2(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse2(e,d);
    return d;
  }

  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  /**
   * Accumulates alignment errors in reverse direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse(float[][] e, float[][] d) {
    accumulate(-1,_bstrain1,e,d);
  }

  /**
   * Accumulates alignment errors in forward direction in 1st dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward1(float[][][] e, float[][][] d) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateForward(e[i2],d[i2]);
  }

  /**
   * Accumulates alignment errors in reverse direction in 1st dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse1(float[][][] e, float[][][] d) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateReverse(e[i2],d[i2]);
  }

  /**
   * Accumulates alignment errors in forward direction in 2nd dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward2(float[][][] e, float[][][] d) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] di1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ei1[i2] = e[i2][i1];
        di1[i2] = d[i2][i1];
      }
      accumulate( 1,_bstrain2,ei1,di1);
    }
  }

  /**
   * Accumulates alignment errors in reverse direction in 2nd dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse2(float[][][] e, float[][][] d) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] di1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ei1[i2] = e[i2][i1];
        di1[i2] = d[i2][i1];
      }
      accumulate(-1,_bstrain2,ei1,di1);
    }
  }

  /**
   * Returns shifts found by backtracking in reverse.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackReverse(d,e,u);
    return u;
  }

  /**
   * Returns shifts found by backtracking in reverse in 1st dimension.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[][] backtrackReverse1(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    backtrackReverse1(d,e,u);
    return u;
  }

  /**
   * Returns shifts found by backtracking in reverse in 2nd dimension.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[][] backtrackReverse2(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    backtrackReverse2(d,e,u);
    return u;
  }

  /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,_lmin,d,e,u);
  }

  /**
   * Computes shifts by backtracking in reverse direction in 1st dimension.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse1(float[][][] d, float[][][] e, float[][] u) {
    int n2 = d.length;
    for (int i2=0; i2<n2; ++i2)
      backtrackReverse(d[i2],e[i2],u[i2]);
  }

  /**
   * Computes shifts by backtracking in reverse direction in 2nd dimension.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse2(float[][][] d, float[][][] e, float[][] u) {
    int n1 = d[0].length;
    int n2 = d.length;
    float[][] di1 = new float[n2][];
    float[][] ei1 = new float[n2][];
    float[] ui1 = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        di1[i2] = d[i2][i1];
        ei1[i2] = e[i2][i1];
      }
      backtrack(-1,_bstrain2,_lmin,di1,ei1,ui1);
      for (int i2=0; i2<n2; ++i2)
        u[i2][i1] = ui1[i2];
    }
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin = e[0][0];
    float emax = e[0][0];
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
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

  /**
   * Returns the sum of errors for specified shifts, rounded to integers.
   * @param e array[n1][nl] of errors.
   * @param u array[n1] of shifts.
   * @return the sum of errors.
   */
  public float sumErrors(float[][] e, float[] u) {
    int n1 = e.length;
    int nl = e[0].length;
    float ul = 0.5f-_lmin;
    double sum = 0.0;
    for (int i1=0; i1<n1; ++i1) {
      int il = (int)(u[i1]+ul);
      il = max(0,min(nl-1,il));
      sum += e[i1][il];
    }
    return (float)sum;
  }

  /**
   * Returns the sum of errors for specified shifts, rounded to integers.
   * @param e array[n2][n1][nl] of errors.
   * @param u array[n2][n1] of shifts.
   * @return the sum of errors.
   */
  public float sumErrors(float[][][] e, float[][] u) {
    int n2 = e.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2)
      sum += sumErrors(e[i2],u[i2]);
    return (float)sum;
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n1][nl] of errors.
   * @return transposed array[nl][n1] of errors.
   */
  public static float[][] transposeLag(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] t = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        t[il][i1] = e[i1][il];
      }
    }
    return t;
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n2][n1][nl] of errors.
   * @return transposed array[nl][n2][n1] of errors.
   */
  public static float[][][] transposeLag(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[il][i2][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nl; // number of lags
  private int _lmin,_lmax; // min,max lags
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private int _esmooth = 0; // number of nonlinear smoothings of errors
  private double _usmooth1 = 0.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 0.0; // extent of smoothing shifts in 2nd dim
  private double _usmooth3 = 0.0; // extent of smoothing shifts in 3rd dim
  private int _bstrain1 = 1; // inverse of bound on strain in 1st dimension
  private int _bstrain2 = 1; // inverse of bound on strain in 2nd dimension
  private int _bstrain3 = 1; // inverse of bound on strain in 3rd dimension
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private RecursiveExponentialFilter _ref2; // for smoothing shifts
  private RecursiveExponentialFilter _ref3; // for smoothing shifts
  private SincInterpolator _si; // for warping with non-integer shifts
  private File _edir; // directory used for temporary files; null for default

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  /**
   * Computes alignment errors.
   * @param f input array[ni] for sequence f.
   * @param g input array[ni] for sequence g.
   * @param e output array[ni][nl] of alignment errors.
   */
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1 = f.length;
    int nl = _nl;
    int n1m = n1-1;
    // 0 <= il < nl, where il is index for lag
    // 0 <= i1 < n1, where i1 is index for sequence f
    // 0 <= j1 < n1, where j1 index for sequence g
    // j1 = i1+il+lmin, where il+lmin = lag
    // 0 <= i1+il+lmin < n1, so that j1 is in bounds
    // max(0,-lmin-i1) <= il < min(nl,n1-lmin-i1)
    // max(0,-lmin-il) <= i1 < min(n1,n1-lmin-il)
    // j1 = 0    => i1 = -lmin-il
    // j1 = n1-1 => i1 = n1-1-lmin-il
    for (int i1=0; i1<n1; ++i1) {
      int illo = min(nl-1,max(0,-_lmin-i1)); // see notes
      int ilhi = max(0,min(nl,n1-_lmin-i1)); // above
      for (int il=0,j1=i1+il+_lmin; il<illo; ++il,++j1)
        e[i1][il] = (j1>=0) ?
          error(f[i1],g[j1]) : 
          error(f[-_lmin-il],g[0]);
      for (int il=illo,j1=i1+il+_lmin; il<ilhi; ++il,++j1)
        e[i1][il] = error(f[i1],g[j1]);
      for (int il=ilhi,j1=i1+il+_lmin; il<nl; ++il,++j1)
        e[i1][il] = (j1<n1) ?
          error(f[i1],g[j1]) : 
          error(f[n1-1-_lmin-il],g[n1-1]);
    }
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param lmin minimum lag corresponding to lag index zero.
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private static void backtrack(
    int dir, int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?nim1:0;
    int is = (dir>0)?1:-1;
    int ii = ib;
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il+lmin;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=is) { // ii-1
        dm += e[kb][ilm1];
        dp += e[kb][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      ii += is;
      u[ii] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        for (int kb=ji; kb!=jb; kb+=is) {
          ii += is;
          u[ii] = il+lmin;
        }
      }
    }
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
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          ef[i2][i1][il] = (ef[i2][i1][il]-eshift)*escale;
        }
      }
    }});
  }

  /**
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors1(int b, float[][] e, float[][] es) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulate( 1,b,e,ef);
    accumulate(-1,b,e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }

  /**
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors1(int b, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      smoothErrors1(bf,ef[i2],esf[i2]);
    }});
  }

  /**
   * Smooths alignment errors in 2nd dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 2nd dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors2(int b, float[][][] e, float[][][] es) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int bf = b;
    final float[][][]  ef = e;
    final float[][][] esf = es;
    final Parallel.Unsafe<float[][][]> eeu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = eeu.get();
      if (ee==null) eeu.set(ee=new float[4][n2][nl]);
      float[][]  e1 = ee[0];
      float[][] es1 = ee[1];
      float[][] ef1 = ee[2];
      float[][] er1 = ee[3];
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  ef[i2][i1];
        es1[i2] = esf[i2][i1];
        for (int il=0; il<nl; ++il) {
          ef1[i2][il] = 0.0f;
          er1[i2][il] = 0.0f;
        }
      }
      accumulate( 1,bf,e1,ef1);
      accumulate(-1,bf,e1,er1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
        }
      }
    }});
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private static float[] like(float[] a) {
    return new float[a.length];
  }
  private static float[][] like(float[][] a) {
    return new float[a.length][a[0].length];
  }
  private static float[][][] like(float[][][] a) {
    return new float[a.length][a[0].length][a[0][0].length];
  }

  ///////////////////////////////////////////////////////////////////////////
  // for 3D image warping

  /**
   * File-based array of alignment errors for 3D image warping.
   * Alignment errors for 3D image warping are likely to require more
   * memory than is available. This class provides efficient access 
   * to slices of these arrays stored in a random-access file.
   */
  static class AlignmentErrors3 {
    AlignmentErrors3(File dir, int nl, int n1, int n2, int n3) {
      _nl = nl;
      _n1 = n1;
      _n2 = n2;
      _n3 = n3;
      _file = null;
      try {
        if (dir!=null) {
          _file = File.createTempFile("tmpae",".tmp",dir);
        } else {
          _file = File.createTempFile("tmpae",".tmp");
        }
        _af = new ArrayFile(_file,"rw");
      } catch (IOException exception) {
        throw new RuntimeException(exception);
      }
    }
    int getNL() { return (int)_nl; }
    int getN1() { return (int)_n1; }
    int getN2() { return (int)_n2; }
    int getN3() { return (int)_n3; }
    void close() {
      try {
        _af.close();
      } catch (IOException exception) {
        throw new RuntimeException(exception);
      }
    }
    void delete() {
      _file.delete();
    }
    void get2(int i2, float[][][] e) {
      try {
        for (int i3=0; i3<_n3; ++i3) {
          seek(i2,i3);
          _af.readFloats(e[i3]);
        }
      } catch (IOException exception) {
        throw new RuntimeException(exception);
      }
    }
    void set2(int i2, float[][][] e) {
      try {
        for (int i3=0; i3<_n3; ++i3) {
          seek(i2,i3);
          _af.writeFloats(e[i3]);
        }
      } catch (IOException exception) {
        throw new RuntimeException(exception);
      }
    }
    void get3(int i3, float[][][] e) {
      try {
        seek(0,i3);
        _af.readFloats(e);
      } catch (IOException exception) {
        throw new RuntimeException(exception);
      }
    }
    void set3(int i3, float[][][] e) {
      try {
        seek(0,i3);
        _af.writeFloats(e);
      } catch (IOException exception) {
        throw new RuntimeException(exception);
      }
    }
    private long _nl,_n1,_n2,_n3; // number of lags, dimensions of 3D array
    private ArrayFile _af; // random-access file for alignment errors
    private File _file; // file name for random-access file
    private void seek(int i2, int i3) throws IOException {
      _af.seek(4L*_nl*_n1*(i2+_n2*i3));
    }
  }
  // In-memory version, for debugging only.
  static class xAlignmentErrors3 {
    xAlignmentErrors3(File dir, int nl, int n1, int n2, int n3) {
      _nl = nl;
      _n1 = n1;
      _n2 = n2;
      _n3 = n3;
      _ae = new float[n3][n2][n1][nl];
    }
    int getNL() { return (int)_nl; }
    int getN1() { return (int)_n1; }
    int getN2() { return (int)_n2; }
    int getN3() { return (int)_n3; }
    void close() {
    }
    void delete() {
    }
    void get2(int i2, float[][][] e) {
      for (int i3=0; i3<_n3; ++i3)
        copy(_ae[i3][i2],e[i3]);
    }
    void set2(int i2, float[][][] e) {
      for (int i3=0; i3<_n3; ++i3)
        copy(e[i3],_ae[i3][i2]);
    }
    void get3(int i3, float[][][] e) {
      copy(_ae[i3],e);
    }
    void set3(int i3, float[][][] e) {
      copy(e,_ae[i3]);
    }
    private long _nl,_n1,_n2,_n3; // number of lags, dimensions of 3D array
    float[][][][] _ae;
  }

  private AlignmentErrors3 computeErrors(float[][][] f, float[][][] g) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] e = new float[n2][n1][_nl];
    AlignmentErrors3 ae = new AlignmentErrors3(_edir,_nl,n1,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
      final float[][] f3 = f[i3];
      final float[][] g3 = g[i3];
      Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        computeErrors(f3[i2],g3[i2],e[i2]);
      }});
      ae.set3(i3,e);
    }
    normalizeErrors(ae);
    return ae;
  }

  private static void normalizeErrors(AlignmentErrors3 ae) {
    final int nl = ae.getNL();
    final int n1 = ae.getN1();
    final int n2 = ae.getN2();
    final int n3 = ae.getN3();
    final float[][][] e = new float[n2][n1][nl];
    float emin =  Float.MAX_VALUE;
    float emax = -Float.MAX_VALUE;
    for (int i3=0; i3<n3; ++i3) {
      ae.get3(i3,e);
      MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
      public MinMax compute(int i2) {
        float emin =  Float.MAX_VALUE;
        float emax = -Float.MAX_VALUE;
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            float ei = e[i2][i1][il];
            if (ei<emin) emin = ei;
            if (ei>emax) emax = ei;
          }
        }
        return new MinMax(emin,emax);
      }
      public MinMax combine(MinMax mm1, MinMax mm2) {
        return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
      }});
      emin = min(emin,mm.emin);
      emax = max(emax,mm.emax);
      ae.set3(i3,e);
    }
    shiftAndScale(emin,emax,ae);
  }
  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }

  private static void shiftAndScale(
    float emin, float emax, AlignmentErrors3 ae) 
  {
    System.out.println("sAS: emin="+emin+" emax="+emax);
    final int nl = ae.getNL();
    final int n1 = ae.getN1();
    final int n2 = ae.getN2();
    final int n3 = ae.getN3();
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] e = new float[n2][n1][nl];
    for (int i3=0; i3<n3; ++i3) {
      ae.get3(i3,e);
      Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            e[i2][i1][il] = (e[i2][i1][il]-eshift)*escale;
          }
        }
      }});
      ae.set3(i3,e);
    }
  }

  private void smoothErrors(AlignmentErrors3 ae) {
    int nl = ae.getNL();
    int n1 = ae.getN1();
    int n2 = ae.getN2();
    int n3 = ae.getN3();
    float[][][] e = new float[n2][n1][nl];
    for (int i3=0; i3<n3; ++i3) {
      ae.get3(i3,e);
      smoothErrors1(_bstrain1,e,e);
      ae.set3(i3,e);
    }
    normalizeErrors(ae);
    for (int i3=0; i3<n3; ++i3) {
      ae.get3(i3,e);
      smoothErrors2(_bstrain2,e,e);
      ae.set3(i3,e);
    }
    normalizeErrors(ae);
    e = new float[n3][n1][nl];
    for (int i2=0; i2<n2; ++i2) {
      ae.get2(i2,e);
      smoothErrors2(_bstrain3,e,e);
      ae.set2(i2,e);
    }
    normalizeErrors(ae);
  }

  private void computeShifts(AlignmentErrors3 ae, float[][][] u) {
    final int nl = ae.getNL();
    final int n1 = ae.getN1();
    final int n2 = ae.getN2();
    final int n3 = ae.getN3();
    final float[][][] e = new float[n2][n1][nl];
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    for (int i3=0; i3<n3; ++i3) {
      ae.get3(i3,e);
      final float[][] ui3 = u[i3];
      Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[][] d = du.get();
        if (d==null) du.set(d=new float[n1][nl]);
        accumulateForward(e[i2],d);
        backtrackReverse(d,e[i2],ui3[i2]);
      }});
    }
  }

  private void smoothShifts(float[][][] u) {
    _ref1.apply1(u,u);
    _ref2.apply2(u,u);
    _ref3.apply3(u,u);
  }
}
