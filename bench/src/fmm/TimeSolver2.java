/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

// for testing only
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;

/**
 * A solver for 2D anisotropic eikonal equations. The non-linear equations 
 * are sqrt(grad(t) dot D*grad(t)) = 1, where t is the solution time field, 
 * and D denotes a positive-definite (velocity-squared) tensor field.
 * <p>
 * This solver uses an iterative method to compute the solution times t.
 * Iterations are similar to those described by Jeong and Whitaker (2007).
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.22
 */
public class TimeSolver2 {

  /**
   * Type of concurrency used when solving for times.
   */
  public enum Concurrency {
    PARALLEL,
    SERIAL
  };

  /**
   * An interface for classes of velocity-squared tensors. Each tensor is a
   * symmetric positive-definite 2-by-2 matrix {{d11,d12},{d12,d22}}.
   */
  public interface Tensors {

    /**
     * Gets tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param d array {d11,d12,d22} of tensor elements.
     */
    public void getTensor(int i1, int i2, float[] d);
  }

  /**
   * Constructs a solver with constant identity tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public TimeSolver2(int n1, int n2) {
    this(n1,n2,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors velocity-squared tensors.
   */
  public TimeSolver2(int n1, int n2, Tensors tensors) {
    init(n1,n2,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors velocity-squared tensors.
   */
  public TimeSolver2(float[][] t, Tensors tensors) {
    init(t[0].length,t.length,t,tensors);
  }

  /**
   * Sets the type of concurrency used to solve for times.
   * The default concurrency is parallel.
   * @param concurrency the type of concurrency.
   */
  public void setConcurrency(Concurrency concurrency) {
    _concurrency = concurrency;
  }

  /**
   * Zeros the time at the specified sample and updates times elsewhere.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @return the updated array of times; by reference, not by copy.
   */
  public float[][] zeroAt(int i1, int i2) {
    solveFrom(i1,i2);
    return _t;
  }

  /**
   * Gets the array of times computed by this solver.
   * @return array of times; by reference, not by copy.
   */
  public float[][] getTimes() {
    return _t;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Default time for samples not yet computed.
  private static final float INFINITY = Float.MAX_VALUE;

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.001f;

  private int _n1,_n2;
  private int _n1m,_n2m;
  private Tensors _tensors;
  private float[][] _t;
  private Sample[][] _s;
  private Concurrency _concurrency = Concurrency.PARALLEL;

  private void init(int n1, int n2, float[][] t, Tensors tensors) {
    _n1 = n1;
    _n2 = n2;
    _n1m = n1-1;
    _n2m = n2-1;
    _tensors = tensors;
    _t = (t!=null)?t:Array.fillfloat(INFINITY,n1,n2);
    _s = new Sample[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        _s[i2][i1] = new Sample(i1,i2);
  }

  // Diffusion tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, float[] d) {
      d[0] = 1.00f; // d11
      d[1] = 0.00f; // d12
      d[2] = 1.00f; // d22
    }
  }

  // Sample index offsets for four neighbor samples.
  // Must be consistent with the neighbor sets below.
  private static final int[] K1 = {-1, 1, 0, 0};
  private static final int[] K2 = { 0, 0,-1, 1};

  // Sets of neighbor sample offsets used to compute times. These must
  // be consistent with the offsets above. For example, when updating the 
  // neighbor with offsets {K1[1],K2[1]} = {1,0}, only the sets K1S[1] 
  // and K2S[1] are used. The sets K1S[4] and K2S[4] are special offsets 
  // for all four neighbors. Indices in each set are ordered so that tris
  // are first and edges last.
  private static final int[][] K1S = {
    { 1, 1, 1},
    {-1,-1,-1},
    {-1, 1, 0},
    {-1, 1, 0},
    {-1, 1,-1, 1,-1, 1, 0, 0}};
  private static final int[][] K2S = {
    {-1, 1, 0},
    {-1, 1, 0},
    { 1, 1, 1},
    {-1,-1,-1},
    {-1,-1, 1, 1, 0, 0,-1, 1}};

  // A sample has indices and a flag used to build the active list.
  private static class Sample {
    int i1,i2; // sample indices
    boolean absent; // used to build active lists
    Sample(int i1, int i2) {
      this.i1 = i1;
      this.i2 = i2;
    }
  }

  // List of active samples.
  private class ActiveList {
    void append(Sample s) {
      if (_n==_a.length)
        growTo(2*_n);
      _a[_n++] = s;
    }
    boolean isEmpty() {
      return _n==0;
    }
    int size() {
      return _n;
    }
    Sample get(int i) {
      return _a[i];
    }
    void clear() {
      _n = 0;
    }
    void markAllAbsent() {
      for (int i=0; i<_n; ++i)
        _a[i].absent = true;
    }
    void appendIfAbsent(ActiveList al) {
      if (_n+al._n>_a.length)
        growTo(2*(_n+al._n));
      int n = al._n;
      for (int i=0; i<n; ++i) {
        Sample s = al.get(i);
        if (s.absent) {
          _a[_n++] = s;
          s.absent = false;
        }
      }
    }
    void shuffle() { // experiment: randomizes order of samples in this list
      Random r = new Random();
      for (int i=0; i<_n; ++i) {
        int j = r.nextInt(_n);
        int k = r.nextInt(_n);
        Sample aj = _a[j];
        _a[j] = _a[k];
        _a[k] = aj;
      }
    }
    void dump() { // debugging: prints this list
      trace("ActiveList.dump: n="+_n);
      for (int i=0; i<_n; ++i) {
        Sample s = _a[i];
        trace(" s["+i+"] = ("+s.i1+","+s.i2+")");
      }
    }
    private int _n;
    private Sample[] _a = new Sample[1024];
    private void growTo(int capacity) {
      Sample[] a = new Sample[capacity];
      System.arraycopy(_a,0,a,0,_n);
      _a = a;
    }
  }

  /**
   * Zeros the time for the specified sample and recursively updates times 
   * for neighbor samples until all times have converged.
   */
  private void solveFrom(int i1, int i2) {

    // Zero the time for the specified sample.
    _t[i2][i1] = 0.0f;

    // Put the sample with zero time into the active list.
    ActiveList al = new ActiveList();
    Sample si = _s[i2][i1];
    al.append(si);

    // Complete the solve by processing the active list until it is empty.
    if (_concurrency==Concurrency.PARALLEL) {
      solveParallel(al);
    } else {
      solveSerial(al);
    }
  }

  /**
   * Solves for times by sequentially processing each sample in active list.
   */
  private void solveSerial(ActiveList al) {
    float[] d = new float[3];
    ActiveList bl = new ActiveList();
    int ntotal = 0;
    while (!al.isEmpty()) {
      //al.shuffle(); // demonstrate that solution depends on order
      int n = al.size();
      ntotal += n;
      for (int i=0; i<n; ++i) {
        Sample s = al.get(i);
        solveOne(s,bl,d);
      }
      bl.markAllAbsent();
      al.clear();
      al.appendIfAbsent(bl);
      bl.clear();
    }
    trace("solveSerial: ntotal="+ntotal);
    trace("             nratio="+(float)ntotal/(float)(_n1*_n2));
  }
  
  /**
   * Solves for times by processing samples in the active list in parallel.
   */
  private void solveParallel(final ActiveList al) {
    int nthread = Runtime.getRuntime().availableProcessors();
    /////////////////////////////////////////////////////////////////////////
    // Benchmarks: 07/24/2008
    // Anisotropic (angle=110.0,su=1.00,sv=0.01) constant tensor with 
    // zero time at center sample of 2D array.
    // Intel 2.4 GHz Core 2 Duo for size 2001^2
    // serial         5.4 s
    //nthread = 1; // 5.6 s
    //nthread = 2; // 3.0 s
    /////////////////////////////////////////////////////////////////////////
    ExecutorService es = Executors.newFixedThreadPool(nthread);
    CompletionService<Void> cs = new ExecutorCompletionService<Void>(es);
    ActiveList[] bl = new ActiveList[nthread];
    float[][] d = new float[nthread][];
    for (int ithread=0; ithread<nthread; ++ithread) {
      bl[ithread] = new ActiveList();
      d[ithread] = new float[3];
    }
    final AtomicInteger ai = new AtomicInteger();
    int ntotal = 0;
    while (!al.isEmpty()) {
      ai.set(0); // initialize the shared block index to zero
      final int n = al.size(); // number of samples in active (A) list
      ntotal += n;
      final int mb = 32; // size of blocks of samples
      final int nb = 1+(n-1)/mb; // number of blocks of samples
      int ntask = min(nb,nthread); // number of tasks (threads to be used)
      for (int itask=0; itask<ntask; ++itask) { // for each task, ...
        final ActiveList bltask = bl[itask]; // task-specific B list 
        final float[] dtask = d[itask]; // task-specific work array
        cs.submit(new Callable<Void>() { // submit new task
          public Void call() {
            for (int ib=ai.getAndIncrement(); ib<nb; ib=ai.getAndIncrement()) {
              int i = ib*mb; // beginning of block
              int j = min(i+mb,n); // beginning of next block (or end)
              for (int k=i; k<j; ++k) { // for each sample in block, ...
                Sample s = al.get(k); // get k'th sample from A list
                solveOne(s,bltask,dtask); // process the sample
              }
            }
            bltask.markAllAbsent(); // needed when merging B lists below
            return null;
          }
        });
      }
      try {
        for (int itask=0; itask<ntask; ++itask)
          cs.take();
      } catch (InterruptedException e) {
        throw new RuntimeException(e);
      }

      // Merge samples from all B lists to a new A list. As samples
      // are appended, their absent flags are set to false, so that 
      // each sample is appended no more than once to the new A list.
      al.clear();
      for (int itask=0; itask<ntask; ++itask) {
        al.appendIfAbsent(bl[itask]);
        bl[itask].clear();
      }
    }
    es.shutdown();
    trace("solveParallel: ntotal="+ntotal);
    trace("               nratio="+(float)ntotal/(float)(_n1*_n2));
  }

  /**
   * Processes one sample from the A list.
   * Appends samples not yet converged to the B list.
   */
  private void solveOne(Sample s, ActiveList bl, float[] d) {

    // Sample indices.
    int i1 = s.i1;
    int i2 = s.i2;

    // Current time and new time computed from all neighbors.
    float ti = _t[i2][i1];
    float gi = computeTime(i1,i2,K1S[4],K2S[4],d);
    _t[i2][i1] = gi;

    // If new and current times are close enough (converged), then ...
    if (ti-gi<=ti*EPSILON) {

      // For all four neighbors, ...
      for (int k=0; k<4; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;

        // Compute time for neighbor.
        float tj = _t[j2][j1];
        float gj = computeTime(j1,j2,K1S[k],K2S[k],d);

        // If computed time significantly less than neighbor's current time, ...
        if (tj-gj>tj*EPSILON) {

          // Replace the current time.
          _t[j2][j1] = gj;
          
          // Append neighbor to the B list.
          bl.append(_s[j2][j1]);
        }
      }
    }

    // Else, if not converged, append this sample to the B list.
    else {
      bl.append(s);
    }
  }

  /**
   * Solves a quadratic equation for a positive time t0.
   * The equation is:
   *   d11*s1*s1*(t1-t0)*(t1-t0) + 
   * 2*d12*s1*s2*(t1-t0)*(t2-t0) + 
   *   d22*s2*s2*(t2-t0)*(t2-t0) = 1
   * To reduce rounding errors, this method actually solves for u = t0-t1,
   * via the following equation:
   *   ds11*(u    )*(u    ) + 
   *   ds22*(u+t12)*(u+t12) +
   * 2*ds12*(u    )*(u+t12) = 1
   * It then returns t0 = t1+u. If no solution exists, because the 
   * discriminant is negative, this method returns INFINITY.
   */
  private static float solveQuadratic(
    float d11, float d12, float d22,
    float s1, float s2, float t1, float t2) 
  {
    double ds11 = d11*s1*s1;
    double ds12 = d12*s1*s2;
    double ds22 = d22*s2*s2;
    double t12 = t1-t2;
    double a = ds11+2.0*ds12+ds22;
    double b = 2.0*(ds12+ds22)*t12;
    double c = ds22*t12*t12-1.0;
    double d = b*b-4.0*a*c;
    if (d<0.0) 
      return INFINITY;
    double u = (-b+sqrt(d))/(2.0*a);
    return t1+(float)u; // t0 = t1+u
  }

  /**
   * Jeong's fast test for a valid solution time t0 to H(p1,p2,p3) = 1.
   * Parameters tm and tp are times for samples backward and forward of 
   * the sample with time t0. The parameter k is the index k1, k2, or k3
   * that was used to compute the time, and the parameter p is a critical
   * point of H(p1,p2,p3) for fixed p1, p2, or p3.
   */
  private static boolean isValid(
    int k, float p, float tm, float tp, float t0) 
  {
    float pm = t0-tm;
    float pp = tp-t0;
    int j = -1;
    if (pm<p && p<pp) { // (pm-p) < 0 and (pp-p) > 0
      j = 0;
    } else if (0.5f*(pm+pp)<p) { // (pm-p) < -(pp-p)
      j = 1;
    }
    return j==k;
  }

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor samples with specified indices.
   */
  private float g(int i1, int i2, int[] k1s, int[] k2s, float[] d) {

    // Current time i'th sample and its four neighbors. We must cache all of
    // these now because times may be changed concurrently in other threads.
    // That is, when this method returns, some of these times may be less 
    // than the values cached here, but the logic within this method call 
    // in the current thread will be consistent with these cached values.
    float tc = _t[i2][i1];
    float t1m = (i1>0   )?_t[i2][i1-1]:INFINITY;
    float t1p = (i1<_n1m)?_t[i2][i1+1]:INFINITY;
    float t2m = (i2>0   )?_t[i2-1][i1]:INFINITY;
    float t2p = (i2<_n2m)?_t[i2+1][i1]:INFINITY;

    // Tensor coefficients.
    _tensors.getTensor(i1,i2,d);
    float d11 = d[0];
    float d12 = d[1];
    float d22 = d[2];

    // For all relevant neighbor samples, ...
    for (int k=0; k<k1s.length; ++k) {
      int k1 = k1s[k];
      int k2 = k2s[k];

      // If (p1s,p2-) or (p1s,p2+), ...
      if (k1==0) {
        float t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        float t1 = t2;
        float s2 = k2;
        float s1 = -s2*d12/d11;
        float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
        if (t0<tc && t0>=t2) {
          float p2 = -s1*(t1-t0)*d12/d22;
          if (isValid(k2,p2,t2m,t2p,t0))
            return t0;
        }
      } 
      
      // else, if (p1-,p2s) or (p1+,p2s), ...
      else if (k2==0) {
        float t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        float t2 = t1;
        float s1 = k1;
        float s2 = -s1*d12/d22;
        float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
        if (t0<tc && t0>=t1) {
          float p1 = -s2*(t2-t0)*d12/d11;
          if (isValid(k1,p1,t1m,t1p,t0))
            return t0;
        }
      }
      
      // else, if (p1-,p2-), (p1+,p2-), (p1-,p2+) or (p1+,p2+), ...
      else {
        float t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        float t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        float s1 = k1;
        float s2 = k2;
        float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
        if (t0<tc && t0>=min(t1,t2)) {
          float p1 = -s2*(t2-t0)*d12/d11;
          float p2 = -s1*(t1-t0)*d12/d22;
          if (isValid(k1,p1,t1m,t1p,t0) &&
              isValid(k2,p2,t2m,t2p,t0))
            return t0;
        }
      }
    }

    return tc;
  }

  private float computeTime(int i1, int i2, int[] k1s, int[] k2s, float[] d) {

    // Current time i'th sample and its four neighbors. We must cache all of
    // these now because times may be changed concurrently in other threads.
    // That is, when this method returns, some of these times may be less 
    // than the values cached here, but the logic within this method call 
    // in the current thread will be consistent with these cached values.
    float tc = _t[i2][i1];
    float t1m = (i1>0   )?_t[i2][i1-1]:INFINITY;
    float t1p = (i1<_n1m)?_t[i2][i1+1]:INFINITY;
    float t2m = (i2>0   )?_t[i2-1][i1]:INFINITY;
    float t2p = (i2<_n2m)?_t[i2+1][i1]:INFINITY;

    // Tensor coefficients.
    _tensors.getTensor(i1,i2,d);
    float d11 = d[0];
    float d12 = d[1];
    float d22 = d[2];

    // For all relevant neighbor samples, ...
    for (int k=0; k<k1s.length; ++k) {
      int k1 = k1s[k];
      int k2 = k2s[k];
      float s1,s2,t1,t2;

      // ( 0,-1), ( 0,+1)
      if (k1==0) {
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t1 = t2;
        s2 = k2;
        s1 = -s2*d12/d11;
      } 

      // (-1, 0), (+1, 0)
      else if (k2==0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t2 = t1;
        s1 = k1;
        s2 = -s1*d12/d22;
      }

      // (-1,-1), (-1,+1), (+1,-1), (+1,+1)
      else {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        s1 = k1;
        s2 = k2;
      }

      // Compute time t0 from times t1 and t2.
      float t0 = computeTime(d11,d12,d22,s1,s2,t1,t2);
      if (t0<tc)
        return t0;
    }

    return tc;
  }
  private static float computeTime(
    float d11, float d12, float d22,
    float s1, float s2, float t1, float t2) 
  {
    double ds11 = d11*s1*s1;
    double ds12 = d12*s1*s2;
    double ds22 = d22*s2*s2;
    double t12 = t1-t2;
    double a = ds11+2.0*ds12+ds22;
    double b = 2.0*(ds12+ds22)*t12;
    double c = ds22*t12*t12-1.0;
    double d = b*b-4.0*a*c;
    if (d<0.0)
      return INFINITY;
    double u1 = (-b+sqrt(d))/(2.0*a);
    double u2 = u1+t12;
    if (ds11*u1+ds12*u2 < 0.0 ||
        ds12*u1+ds22*u2 < 0.0)
      return INFINITY;
    float t0 = t1+(float)u1;
    //if (t0<min(t1,t2))
    //  return INFINITY;
    return t0;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void plot(float[][] x) {
    float[][] y = Array.copy(x);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (y[i2][i1]==INFINITY) {
          y[i2][i1] = 0.0f;
        }
      }
    }
    SimplePlot sp = new SimplePlot();
    sp.setSize(800,790);
    PixelsView pv = sp.addPixels(y);
    pv.setColorModel(ColorMap.PRISM);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    //pv.setInterpolation(PixelsView.Interpolation.LINEAR);
  }

  private static class ConstantTensors implements TimeSolver2.Tensors {
    ConstantTensors(float d11, float d12, float d22) {
      _d11 = d11;
      _d12 = d12;
      _d22 = d22;
    }
    public void getTensor(int i1, int i2, float[] d) {
      d[0] = _d11;
      d[1] = _d12;
      d[2] = _d22;
    }
    private float _d11,_d12,_d22;
  }

  private static ConstantTensors makeTensors(float a, float su, float sv) {
    a *= FLT_PI/180.0f;
    float cosa = cos(a);
    float sina = sin(a);
    float d11 = su*cosa*cosa+sv*sina*sina;
    float d12 = (su-sv)*sina*cosa;
    float d22 = sv*cosa*cosa+su*sina*sina;
    return new ConstantTensors(d11,d12,d22);
  }

  private static float[][] computeSerial(
    int n1, int n2, int i1, int i2,
    TimeSolver2.Tensors tensors)
  {
    trace("computeSerial:");
    return computeTimes(
      n1,n2,i1,i2,tensors,TimeSolver2.Concurrency.SERIAL);
  }

  private static float[][] computeParallel(
    int n1, int n2, int i1, int i2,
    TimeSolver2.Tensors tensors)
  {
    trace("computeParallel:");
    return computeTimes(
      n1,n2,i1,i2,tensors,TimeSolver2.Concurrency.PARALLEL);
  }

  private static float[][] computeTimes(
    int n1, int n2, int i1, int i2,
    TimeSolver2.Tensors tensors, 
    TimeSolver2.Concurrency concurrency) 
  {
    TimeSolver2 ts = new TimeSolver2(n1,n2,tensors);
    ts.setConcurrency(concurrency);
    Stopwatch sw = new Stopwatch();
    sw.start();
    ts.zeroAt(i1,i2);
    sw.stop();
    float[][] t = ts.getTimes();
    trace("  time="+sw.time()+" sum="+Array.sum(t));
    return t;
  }

  private static void testConstant() {
    int n1 = 2001, n2 = 2001;
    ConstantTensors tensors = makeTensors(110.0f,1.000f,0.010f);
    int i1 = 2*(n1-1)/4, i2 = 2*(n2-1)/4;
    float[][] ts = computeSerial(n1,n2,i1,i2,tensors);
    float[][] tp = computeParallel(n1,n2,i1,i2,tensors);
    float[][] te = Array.div(Array.abs(Array.sub(tp,ts)),ts);
    te[i2][i1] = 0.0f;
    float temax = Array.max(te);
    trace("temax="+temax);
    trace("********************************************************");
    plot(ts);
    plot(tp);
    //plot(te);
    if (temax>0.1f)
      System.exit(-1);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
  private static void trace(int i1, int i2, String s) {
    //if (i1==2 && i2==2)
    //  trace(s);
  }

  public static void main(String[] args) {
    //SwingUtilities.invokeLater(new Runnable() {
    //  public void run() {
    //  for (;;)
        testConstant();
    //  }
    //});
  }
}
