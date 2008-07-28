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
   * A listener for time changes.
   */
  public interface Listener {

    /**
     * Called when time for the specified sample has decreased.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param t the decreased time.
     */
    public void timeDecreased(int i1, int i2, float t);
  }
  
  /**
   * Constructs a solver for the specified tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors velocity-squared tensors.
   */
  public TimeSolver2(int n1, int n2, Tensors2 tensors) {
    init(n1,n2,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors velocity-squared tensors.
   */
  public TimeSolver2(float[][] t, Tensors2 tensors) {
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
   * Sets the tensors used by this solver.
   * @param tensors the tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = tensors;
  }

  /**
   * Zeros the time at the specified sample and computes times for neighbors.
   * Times of neighbor samples are computed recursively while computed times 
   * are less than current times. In other words, this method only decreases 
   * times in the array of times referenced by this class. Finally, this 
   * method notifies any listeners of all samples with times decreased.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @return the modified array of times; by reference, not by copy.
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

  /**
   * Adds the specified listener.
   * @param listener the listener.
   */
  public void addListener(Listener listener) {
    _listeners.add(listener);
  }

  /**
   * Removes the specified listener.
   * @param listener the listener.
   */
  public void removeListener(Listener listener) {
    _listeners.remove(listener);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Default time for samples not yet computed.
  private static final float INFINITY = Float.MAX_VALUE;

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.001f;

  private int _n1,_n2;
  private int _n1m,_n2m;
  private Tensors2 _tensors;
  private float[][] _t;
  private Sample[][] _s;
  private Concurrency _concurrency = Concurrency.PARALLEL;
  private ArrayList<Listener> _listeners = new ArrayList<Listener>();
  private ArrayList<Sample> _stack = new ArrayList<Sample>(1024);

  private void init(int n1, int n2, float[][] t, Tensors2 tensors) {
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

  // Sample index offsets for four neighbor samples.
  // Must be consistent with the neighbor sets below.
  private static final int[] K1 = {-1, 1, 0, 0};
  private static final int[] K2 = { 0, 0,-1, 1};

  // Sets of neighbor sample offsets used to compute times. These must
  // be consistent with the offsets above. For example, when updating the 
  // neighbor with offsets {K1[1],K2[1]} = {1,0}, only the sets K1S[1] 
  // and K2S[1] are used. The sets K1S[4] and K2S[4] are special offsets 
  // for all four neighbors. Indices in each set are ordered so that tris
  // are first and edges last. Tris are defined by two non-zero offsets,
  // and edges are defined by one.
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
    int marked; // used to mark samples with decreased times
    boolean absent; // used to build active lists
    Sample(int i1, int i2) {
      this.i1 = i1;
      this.i2 = i2;
    }
  }

  // List of active samples.
  private class ActiveList {
    void append(Sample s) {
      s.marked = _marked;
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
    void setAllAbsent() {
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

  // Marks set during computation of times. For efficiency, do not
  // loop over all the marks to clear them before computing times.
  // Instead, modify the value that represents marked samples.
  private int _marked = 1;
  private void clearMarked() {
    if (_marked==Integer.MAX_VALUE) { // rarely!
      _marked = 1;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1].marked = 0;
        }
      }
    } else { // typically
      ++_marked;
    }
  }
  private void mark(Sample s) {
    s.marked = _marked;
  }
  private void unmark(Sample s) {
    s.marked -= 1;
  }
  private boolean isMarked(Sample s) {
    return s.marked==_marked;
  }

  private void fireTimesDecreasedFrom(int i1, int i2) {
    Sample si = _s[i2][i1];
    if (!isMarked(si))
      return;
    int nlistener = _listeners.size();
    if (nlistener==0)
      return;
    _stack.clear();
    _stack.add(si);
    while (!_stack.isEmpty()) {
      si = _stack.remove(_stack.size()-1);
      if (isMarked(si)) {
        unmark(si);
        i1 = si.i1;
        i2 = si.i2;
        float ti = _t[i2][i1];
        for (int i=0; i<nlistener; ++i)
          _listeners.get(i).timeDecreased(i1,i2,ti);
        for (int k=0; k<4; ++k) {
          int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
          int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
          Sample sj = _s[j2][j1];
          if (isMarked(sj))
            _stack.add(sj);
        }
      }
    }
  }

  /**
   * Zeros the time for the specified sample. Recursively updates times 
   * for neighbor samples until all times have converged, and then notifies
   * any listeners of all times decreased.
   */
  private void solveFrom(int i1, int i2) {

    // Clear all marked samples so we can tell which ones have times set.
    clearMarked();

    // Zero the time for the specified sample.
    _t[i2][i1] = 0.0f;

    // Put the sample with zero time into the active list.
    ActiveList al = new ActiveList();
    al.append(_s[i2][i1]);

    // Complete the solve by processing the active list until it is empty.
    if (_concurrency==Concurrency.PARALLEL) {
      solveParallel(al);
    } else {
      solveSerial(al);
    }

    // Notify any listeners of all times decreased.
    fireTimesDecreasedFrom(i1,i2);
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
      bl.setAllAbsent();
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
    // Benchmarks: 07/26/2008
    // Anisotropic (angle=110.0,su=1.00,sv=0.01) constant tensor with 
    // zero time at center sample of 2D array.
    // Intel 2.4 GHz Core 2 Duo for size 2001^2
    // serial         3.8 s
    //nthread = 1; // 4.0 s
    //nthread = 2; // 2.3 s
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
            bltask.setAllAbsent(); // needed when merging B lists below
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
    float ci = computeTime(i1,i2,K1S[4],K2S[4],d);
    _t[i2][i1] = ci;

    // If new and current times are close enough (converged), then ...
    if (ti-ci<=ti*EPSILON) {

      // For all four neighbors, ...
      for (int k=0; k<4; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;

        // Compute time for neighbor.
        float tj = _t[j2][j1];
        float cj = computeTime(j1,j2,K1S[k],K2S[k],d);

        // If computed time significantly less than neighbor's current time, ...
        if (tj-cj>tj*EPSILON) {

          // Replace the current time.
          _t[j2][j1] = cj;
          
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
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor samples with specified indices.
   */
  private float computeTime(int i1, int i2, int[] k1s, int[] k2s, float[] d) {
    _tensors.getTensor(i1,i2,d);
    float d11 = d[0];
    float d12 = d[1];
    float d22 = d[2];
    float e12 = 1.0f/(d11*d22-d12*d12);
    float tc = _t[i2][i1];
    float t1m = (i1>0   )?_t[i2][i1-1]:INFINITY;
    float t1p = (i1<_n1m)?_t[i2][i1+1]:INFINITY;
    float t2m = (i2>0   )?_t[i2-1][i1]:INFINITY;
    float t2p = (i2<_n2m)?_t[i2+1][i1]:INFINITY;
    for (int k=0; k<k1s.length; ++k) {
      int k1 = k1s[k];
      int k2 = k2s[k];
      float t0,t1,t2;
      if (k1!=0 && k2!=0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t0 = computeTime(d11,d12,d22,k1,k2,t1,t2);
      } else if (k1!=0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t0 = t1+sqrt(d22*e12);
      } else { // k2!=0
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t0 = t2+sqrt(d11*e12);
      }
      if (t0<tc)
        return t0;
    }
    return tc;
  }

  /**
   * Solves a 2D anisotropic eikonal equation for a positive time t0.
   * The equation is:
   *   d11*s1*s1*(t1-t0)*(t1-t0) + 
   * 2*d12*s1*s2*(t1-t0)*(t2-t0) + 
   *   d22*s2*s2*(t2-t0)*(t2-t0) = 1
   * To reduce rounding errors, this method actually solves for u = t0-t1,
   * via the following equation:
   *   ds11*(u    )*(u    ) + 
   *   ds22*(u+t12)*(u+t12) +
   * 2*ds12*(u    )*(u+t12) = 1
   * If a valid u can be computed, then the time returned is t0 = t1+u.
   * Otherwise, this method returns INFINITY.
   */
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
    return t1+(float)u1;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void plot(float[][] x, IndexColorModel icm) {
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
    pv.setColorModel(icm);
    //pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setInterpolation(PixelsView.Interpolation.LINEAR);
  }

  private static class ConstantTensors implements Tensors2 {
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

  private static ConstantTensors makeConstantTensors(
    float a, float su, float sv) 
  {
    a *= FLT_PI/180.0f;
    float cosa = cos(a);
    float sina = sin(a);
    float d11 = su*cosa*cosa+sv*sina*sina;
    float d12 = (su-sv)*sina*cosa;
    float d22 = sv*cosa*cosa+su*sina*sina;
    return new ConstantTensors(d11,d12,d22);
  }

  private static class TsaiTensors extends EigenTensors2 {
    TsaiTensors(int n1, int n2) {
      super(n1,n2);
      float a1 = 2.0f*FLT_PI;
      float a2 = 2.0f*FLT_PI;
      float d1 = 2.0f/(float)(n1-1);
      float d2 = 2.0f/(float)(n2-1);
      float f1 = -1.0f;
      float f2 = -1.0f;
      float[][] f = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        float x2 = f2+i2*d2;
        for (int i1=0; i1<n1; ++i1) {
          float x1 = f1+i1*d1;
          f[i2][i1] = cos(a1*x1)*sin(a2*x2);
          float e1 = -a1*sin(a1*x1)*sin(a2*x2);
          float e2 =  a2*cos(a1*x1)*cos(a2*x2);
          float den =  1.0f+e1*e1+e2*e2;
          float d11 =  1.0f-e1*e1/den;
          float d22 =  1.0f-e2*e2/den;
          float d12 = -2.0f*e1*e2/den;
          float[] d = {d11,d12,d22};
          setTensor(i1,i2,d);
        }
      }
      plot(f,ColorMap.JET);
    }
  }

  private static float[][] computeSerial(
    int n1, int n2, int i1, int i2, Tensors2 tensors)
  {
    trace("computeSerial:");
    return computeTimes(
      n1,n2,i1,i2,tensors,TimeSolver2.Concurrency.SERIAL);
  }

  private static float[][] computeParallel(
    int n1, int n2, int i1, int i2, Tensors2 tensors)
  {
    trace("computeParallel:");
    return computeTimes(
      n1,n2,i1,i2,tensors,TimeSolver2.Concurrency.PARALLEL);
  }

  private static float[][] computeTimes(
    int n1, int n2, int i1, int i2,
    Tensors2 tensors, TimeSolver2.Concurrency concurrency) 
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

  private static void benchConstant() {
    int n1 = 2001, n2 = 2001;
    ConstantTensors tensors = makeConstantTensors(110.0f,1.000f,0.010f);
    int i1 = 2*(n1-1)/4, i2 = 2*(n2-1)/4;
    for (;;) {
      computeSerial(n1,n2,i1,i2,tensors);
      computeParallel(n1,n2,i1,i2,tensors);
      trace("********************************************************");
    }
  }

  private static void testConstant() {
    int n1 = 1001, n2 = 1001;
    ConstantTensors tensors = makeConstantTensors(110.0f,1.000f,0.010f);
    int i1 = 2*(n1-1)/4, i2 = 2*(n2-1)/4;
    float[][] ts = computeSerial(n1,n2,i1,i2,tensors);
    float[][] tp = computeParallel(n1,n2,i1,i2,tensors);
    float[][] te = Array.div(Array.abs(Array.sub(tp,ts)),ts);
    te[i2][i1] = 0.0f;
    float temax = Array.max(te);
    trace("temax="+temax);
    plot(ts,ColorMap.PRISM);
    plot(tp,ColorMap.PRISM);
    plot(te,ColorMap.JET);
    if (temax>0.1f)
      System.exit(-1);
  }

  private static void testTsai() {
    int n1 = 101, n2 = 101;
    TsaiTensors tensors = new TsaiTensors(n1,n2);
    int i1 = 2*(n1-1)/4, i2 = 2*(n2-1)/4;
    TimeSolver2 ts = new TimeSolver2(n1,n2,tensors);
    ts.zeroAt(50,50);
    ts.zeroAt(10,25);
    float[][] t = ts.getTimes();
    //plot(t,ColorMap.PRISM);
    plot(t,ColorMap.JET);
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
        //benchConstant();
        //testConstant();
        testTsai();
    //  }
    //});
  }
}
