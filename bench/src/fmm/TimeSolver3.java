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
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.sgl.test.*;

/**
 * A solver for 3D anisotropic eikonal equations. The non-linear equations 
 * are sqrt(grad(t) dot D*grad(t)) = 1, where t is the solution time field, 
 * and D denotes a positive-definite (velocity-squared) tensor field.
 * <p>
 * This solver uses an iterative method to compute the solution times t.
 * Iterations are similar to those described by Jeong and Whitaker (2007).
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.22
 */
public class TimeSolver3 {

  /**
   * Type of concurrency used when solving for times.
   */
  public enum Concurrency {
    PARALLEL,
    SERIAL
  };

  /**
   * An interface for classes of velocity-squared tensors. Each tensor is a
   * symmetric positive-definite 3-by-3 matrix 
   * {{d11,d12,d13},{d12,d22,d23},{d13,d23,d33}}.
   */
  public interface Tensors {

    /**
     * Gets tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param i3 index for 3rd dimension.
     * @param d array {d11,d12,d13,d22,d23,d33} of tensor elements.
     */
    public void getTensor(int i1, int i2, int i3, float[] d);
  }

  /**
   * A listener for time changes.
   */
  public interface Listener {

    /**
     * Called when time for the specified sample has decreased.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param i3 index for 3rd dimension.
     * @param t the decreased time.
     */
    public void timeDecreased(int i1, int i2, int i3, float t);
  }

  /**
   * Constructs a solver with constant identity tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public TimeSolver3(int n1, int n2, int n3) {
    this(n1,n2,n3,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param tensors velocity-squared tensors.
   */
  public TimeSolver3(int n1, int n2, int n3, Tensors tensors) {
    init(n1,n2,n3,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors velocity-squared tensors.
   */
  public TimeSolver3(float[][][] t, Tensors tensors) {
    init(t[0][0].length,t[0].length,t.length,t,tensors);
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
   * Zeros the time at the specified sample and computes times for neighbors.
   * Times of neighbor samples are computed recursively while computed times 
   * are less than current times. In other words, this method only decreases 
   * times in the array of times referenced by this class. Finally, this 
   * method notifies any listeners of all samples with times decreased.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @param i3 index in 3rd dimension of time to zero.
   * @return the modified array of times; by reference, not by copy.
   */
  public float[][][] zeroAt(int i1, int i2, int i3) {
    solveFrom(i1,i2,i3);
    return _t;
  }

  /**
   * Gets the array of times computed by this solver.
   * @return array of times; by reference, not by copy.
   */
  public float[][][] getTimes() {
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

  // Times converged when fractional change is less than this value.
  private static final float EPSILON = 0.001f;

  private int _n1,_n2,_n3;
  private int _n1m,_n2m,_n3m;
  private Tensors _tensors;
  private float[][][] _t;
  private Sample[][][] _s;
  private Concurrency _concurrency = Concurrency.PARALLEL;
  private ArrayList<Listener> _listeners = new ArrayList<Listener>();
  private ArrayList<Sample> _stack = new ArrayList<Sample>(1024);

  private void init(int n1, int n2, int n3, float[][][] t, Tensors tensors) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _n1m = n1-1;
    _n2m = n2-1;
    _n3m = n3-1;
    _tensors = tensors;
    _t = (t!=null)?t:Array.fillfloat(INFINITY,n1,n2,n3);
    _s = new Sample[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          _s[i3][i2][i1] = new Sample(i1,i2,i3);
  }

  // Diffusion tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, int i3, float[] d) {
      d[0] = 1.00f; // d11
      d[1] = 0.00f; // d12
      d[2] = 0.00f; // d13
      d[3] = 1.00f; // d22
      d[4] = 0.00f; // d23
      d[5] = 1.00f; // d33
    }
  }

  // Sample index offsets for six neighbor samples.
  // Must be consistent with the neighbor sets below.
  private static final int[] K1 = {-1, 1, 0, 0, 0, 0};
  private static final int[] K2 = { 0, 0,-1, 1, 0, 0};
  private static final int[] K3 = { 0, 0, 0, 0,-1, 1};

  // Sets of neighbor sample offsets used to compute times. These must
  // be consistent with the offsets above. For example, when updating the 
  // neighbor with offsets {K1[1],K2[1],K3[1]} = {1,0,0}, only the sets 
  // K1S[1], K2S[1], and K3S[1] are used. The sets K1S[6], K2S[6], and 
  // K3S[6] are special offsets for all six neighbors. Indices in each
  // set are ordered so that tets are first, tris next, and edges last.
  private static final int[][] K1S = {
    { 1, 1, 1, 1, 1, 1, 1, 1, 1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1}, // A
    {-1, 1,-1, 1,-1, 1, 0, 0, 0}, // B
    {-1, 1,-1, 1,-1, 1, 0, 0, 0}, // B
    {-1,-1, 1, 1, 0, 0,-1, 1, 0}, // C
    {-1,-1, 1, 1, 0, 0,-1, 1, 0}, // C
    {-1, 1,-1, 1,-1, 1,-1, 1,             //    8 tets 
     -1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0, // + 12 tris
     -1, 1, 0, 0, 0, 0}};                 // +  6 edges = 26 cases
  private static final int[][] K2S = {
    {-1,-1, 1, 1, 0, 0,-1, 1, 0}, // C
    {-1,-1, 1, 1, 0, 0,-1, 1, 0}, // C
    { 1, 1, 1, 1, 1, 1, 1, 1, 1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1}, // A
    {-1, 1,-1, 1,-1, 1, 0, 0, 0}, // B
    {-1, 1,-1, 1,-1, 1, 0, 0, 0}, // B
    {-1,-1, 1, 1,-1,-1, 1, 1,
     -1,-1, 1, 1, 0, 0, 0, 0,-1, 1,-1, 1,
      0, 0,-1, 1, 0, 0}};
  private static final int[][] K3S = {
    {-1, 1,-1, 1,-1, 1, 0, 0, 0}, // B
    {-1, 1,-1, 1,-1, 1, 0, 0, 0}, // B
    {-1,-1, 1, 1, 0, 0,-1, 1, 0}, // C
    {-1,-1, 1, 1, 0, 0,-1, 1, 0}, // C
    { 1, 1, 1, 1, 1, 1, 1, 1, 1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1}, // A
    {-1,-1,-1,-1, 1, 1, 1, 1,
      0, 0, 0, 0,-1,-1, 1, 1,-1,-1, 1, 1,
      0, 0, 0, 0,-1, 1}};

  // A sample has indices and a flag used to build the active list.
  private static class Sample {
    int i1,i2,i3; // sample indices
    int marked; // used to mark samples when computing times
    boolean absent; // used to build active lists
    Sample(int i1, int i2, int i3) {
      this.i1 = i1;
      this.i2 = i2;
      this.i3 = i3;
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
        trace(" s["+i+"] = ("+s.i1+","+s.i2+","+s.i3+")");
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
      for (int i3=0; i3<_n3; ++i3) {
        for (int i2=0; i2<_n2; ++i2) {
          for (int i1=0; i1<_n1; ++i1) {
            _s[i3][i2][i1].marked = 0;
          }
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

  private void fireTimesDecreasedFrom(int i1, int i2, int i3) {
    Sample si = _s[i3][i2][i1];
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
        i3 = si.i3;
        float ti = _t[i3][i2][i1];
        for (int i=0; i<nlistener; ++i)
          _listeners.get(i).timeDecreased(i1,i2,i3,ti);
        for (int k=0; k<6; ++k) {
          int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
          int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
          int j3 = i3+K3[k];  if (j3<0 || j3>=_n3) continue;
          Sample sj = _s[j3][j2][j1];
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
  private void solveFrom(int i1, int i2, int i3) {

    // Zero the time for the specified sample.
    _t[i3][i2][i1] = 0.0f;

    // Put the sample with zero time into the active list.
    ActiveList al = new ActiveList();
    Sample si = _s[i3][i2][i1];
    al.append(si);

    // Complete the solve by processing the active list until it is empty.
    if (_concurrency==Concurrency.PARALLEL) {
      solveParallel(al);
    } else {
      solveSerial(al);
    }

    // Notify any listeners of all times decreased.
    fireTimesDecreasedFrom(i1,i2,i3);
  }

  /**
   * Solves for times by sequentially processing each sample in active list.
   */
  private void solveSerial(ActiveList al) {
    float[] d = new float[6];
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
    trace("             nratio="+(float)ntotal/(float)(_n1*_n2*_n3));
  }
  
  /**
   * Solves for times by processing samples in the active list in parallel.
   */
  private void solveParallel(final ActiveList al) {
    int nthread = Runtime.getRuntime().availableProcessors();
    /////////////////////////////////////////////////////////////////////////
    // Benchmarks: 07/26/2008
    // Anisotropic constant tensor with zero time at center sample. Tensor
    // coefficients are d11 = d22 = d33 = 1.0, d12 = d13 = d23 = 0.9, as
    // in Jeong and Whitaker's benchmark.
    // Intel 2.4 GHz Core 2 Duo for size 64^3
    // serial         1.2 s
    //nthread = 1; // 1.2 s
    //nthread = 2; // 0.6 s
    /////////////////////////////////////////////////////////////////////////
    ExecutorService es = Executors.newFixedThreadPool(nthread);
    CompletionService<Void> cs = new ExecutorCompletionService<Void>(es);
    ActiveList[] bl = new ActiveList[nthread];
    float[][] d = new float[nthread][];
    for (int ithread=0; ithread<nthread; ++ithread) {
      bl[ithread] = new ActiveList();
      d[ithread] = new float[6];
    }
    final AtomicInteger ai = new AtomicInteger();
    int ntotal = 0;
    while (!al.isEmpty()) {
      ai.set(0); // initialize the shared block index to zero
      final int n = al.size(); // number of samples in active (A) list
      ntotal += n;
      final int mb = 16; // size of blocks of samples
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
    trace("               nratio="+(float)ntotal/(float)(_n1*_n2*_n3));
  }

  /**
   * Processes one sample from the A list.
   * Appends samples not yet converged to the B list.
   */
  private void solveOne(Sample s, ActiveList bl, float[] d) {

    // Sample indices.
    int i1 = s.i1;
    int i2 = s.i2;
    int i3 = s.i3;

    // Current time and new time computed from all neighbors.
    float ti = _t[i3][i2][i1];
    float ci = computeTime(i1,i2,i3,K1S[6],K2S[6],K3S[6],d);
    _t[i3][i2][i1] = ci;

    // If new and current times are close enough (converged), then ...
    if (ti-ci<=ti*EPSILON) {

      // For all six neighbors, ...
      for (int k=0; k<6; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
        int j3 = i3+K3[k];  if (j3<0 || j3>=_n3) continue;

        // Compute time for neighbor.
        float tj = _t[j3][j2][j1];
        float cj = computeTime(j1,j2,j3,K1S[k],K2S[k],K3S[k],d);

        // If computed time significantly less than neighbor's current time, ...
        if (tj-cj>tj*EPSILON) {

          // Replace the current time.
          _t[j3][j2][j1] = cj;
          
          // Append neighbor to the B list.
          bl.append(_s[j3][j2][j1]);
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
   * Computations are limited to neighbor samples with specified offsets.
   */
  private float computeTime(
    int i1, int i2, int i3, int[] k1s, int[] k2s, int[] k3s, float[] d) 
  {
    _tensors.getTensor(i1,i2,i3,d);
    float d11 = d[0];
    float d12 = d[1];
    float d13 = d[2];
    float d22 = d[3];
    float d23 = d[4];
    float d33 = d[5];
    float o11 = 1.0f/d11;
    float o22 = 1.0f/d22;
    float o33 = 1.0f/d33;
    float d1212 = d12*d12;
    float d1213 = d12*d13;
    float d1223 = d12*d23;
    float d1313 = d13*d13;
    float d1323 = d13*d23;
    float d2323 = d23*d23;
    float a11 = d11-d1313*o33;
    float a12 = d12-d1323*o33;
    float a22 = d22-d2323*o33;
    float b11 = d11-d1212*o22;
    float b13 = d13-d1223*o22;
    float b33 = d33-d2323*o22;
    float c22 = d22-d1212*o11;
    float c23 = d23-d1213*o11;
    float c33 = d33-d1313*o11;
    float e12 = 1.0f/(a11*a22-a12*a12);
    float e13 = 1.0f/(b11*b33-b13*b13);
    float tc = _t[i3][i2][i1];
    float t1m = (i1>0   )?_t[i3][i2][i1-1]:INFINITY;
    float t1p = (i1<_n1m)?_t[i3][i2][i1+1]:INFINITY;
    float t2m = (i2>0   )?_t[i3][i2-1][i1]:INFINITY;
    float t2p = (i2<_n2m)?_t[i3][i2+1][i1]:INFINITY;
    float t3m = (i3>0   )?_t[i3-1][i2][i1]:INFINITY;
    float t3p = (i3<_n3m)?_t[i3+1][i2][i1]:INFINITY;
    for (int k=0; k<k1s.length; ++k) {
      int k1 = k1s[k];
      int k2 = k2s[k];
      int k3 = k3s[k];
      float t0,t1,t2,t3;
      if (k1!=0 && k2!=0 && k3!=0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t3 = (k3<0)?t3m:t3p;  if (t3==INFINITY) continue;
        t0 = computeTime(d11,d12,d13,d22,d23,d33,k1,k2,k3,t1,t2,t3);
      } else if (k1!=0 && k2!=0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t0 = computeTime(a11,a12,a22,k1,k2,t1,t2);
      } else if (k1!=0 && k3!=0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t3 = (k3<0)?t3m:t3p;  if (t3==INFINITY) continue;
        t0 = computeTime(b11,b13,b33,k1,k3,t1,t3);
      } else if (k2!=0 && k3!=0) {
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t3 = (k3<0)?t3m:t3p;  if (t3==INFINITY) continue;
        t0 = computeTime(c22,c23,c33,k2,k3,t2,t3);
      } else if (k1!=0) {
        t1 = (k1<0)?t1m:t1p;  if (t1==INFINITY) continue;
        t0 = t1+sqrt(a22*e12);
      } else if (k2!=0) {
        t2 = (k2<0)?t2m:t2p;  if (t2==INFINITY) continue;
        t0 = t2+sqrt(a11*e12);
      } else { // k3!=0
        t3 = (k3<0)?t3m:t3p;  if (t3==INFINITY) continue;
        t0 = t3+sqrt(b11*e13);
      }
      if (t0<tc)
        return t0;
    }
    return tc;
  }

  /**
   * Solves a 3D anisotropic eikonal equation for a positive time t0.
   * The equation is:
   *   d11*s1*s1*(t0-t1)*(t0-t1) + 
   *   d22*s2*s2*(t0-t2)*(t0-t2) +
   *   d33*s3*s3*(t0-t3)*(t0-t3) +
   * 2*d12*s1*s2*(t0-t1)*(t0-t2) + 
   * 2*d13*s1*s3*(t0-t1)*(t0-t3) + 
   * 2*d23*s2*s3*(t0-t2)*(t0-t3) = 1
   * To reduce rounding errors, this method actually solves for u = t0-t1,
   * via the following equation:
   *   ds11*(u    )*(u    ) + 
   *   ds22*(u+t12)*(u+t12) +
   *   ds33*(u+t13)*(u+t13) +
   * 2*ds12*(u    )*(u+t12) + 
   * 2*ds13*(u    )*(u+t13) + 
   * 2*ds23*(u+t12)*(u+t13) = 1
   * If a valid u can be computed, then the time returned is t0 = t1+u.
   * Otherwise, this method returns INFINITY.
   */
  private static float computeTime(
    float d11, float d12, float d13, float d22, float d23, float d33,
    float s1, float s2, float s3, float t1, float t2, float t3) 
  {
    double ds11 = d11*s1*s1;
    double ds22 = d22*s2*s2;
    double ds33 = d33*s3*s3;
    double ds12 = d12*s1*s2;
    double ds13 = d13*s1*s3;
    double ds23 = d23*s2*s3;
    double t12 = t1-t2;
    double t13 = t1-t3;
    double a = ds11+ds22+ds33+2.0*(ds12+ds13+ds23);
    double b = 2.0*((ds22+ds12+ds23)*t12+(ds33+ds13+ds23)*t13);
    double c = ds22*t12*t12+ds33*t13*t13+2.0*ds23*t12*t13-1.0;
    double d = b*b-4.0*a*c;
    if (d<0.0) 
      return INFINITY;
    double u1 = (-b+sqrt(d))/(2.0*a);
    double u2 = u1+t12;
    double u3 = u1+t13;
    if (ds11*u1+ds12*u2+ds13*u3 < 0.0 ||
        ds12*u1+ds22*u2+ds23*u3 < 0.0 ||
        ds13*u1+ds23*u2+ds33*u3 < 0.0)
      return INFINITY;
    return t1+(float)u1;
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

  private static void plot(float[][][] t, IndexColorModel cm) {
    int n1 = t[0][0].length;
    int n2 = t[0].length;
    int n3 = t.length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    ImagePanelGroup ipg = new ImagePanelGroup(s1,s2,s3,new SimpleFloat3(t));
    ipg.setColorModel(cm);
    World world = new World();
    world.addChild(ipg);
    TestFrame frame = new TestFrame(world);
    frame.setVisible(true);
  }

  private static class ConstantTensors implements TimeSolver3.Tensors {
    ConstantTensors(
      float d11, float d12, float d13, float d22, float d23, float d33) 
    {
      _d11 = d11;
      _d12 = d12;
      _d13 = d13;
      _d22 = d22;
      _d23 = d23;
      _d33 = d33;
    }
    public void getTensor(int i1, int i2, int i3, float[] d) {
      d[0] = _d11;
      d[1] = _d12;
      d[2] = _d13;
      d[3] = _d22;
      d[4] = _d23;
      d[5] = _d33;
    }
    private float _d11,_d12,_d13,_d22,_d23,_d33;
  }

  private static float[][][] computeSerial(
    int n1, int n2, int n3,
    int i1, int i2, int i3, 
    TimeSolver3.Tensors tensors)
  {
    trace("computeSerial:");
    return computeTimes(
      n1,n2,n3,i1,i2,i3,tensors,TimeSolver3.Concurrency.SERIAL);
  }

  private static float[][][] computeParallel(
    int n1, int n2, int n3,
    int i1, int i2, int i3, 
    TimeSolver3.Tensors tensors)
  {
    trace("computeParallel:");
    return computeTimes(
      n1,n2,n3,i1,i2,i3,tensors,TimeSolver3.Concurrency.PARALLEL);
  }

  private static float[][][] computeTimes(
    int n1, int n2, int n3,
    int i1, int i2, int i3, 
    TimeSolver3.Tensors tensors, 
    TimeSolver3.Concurrency concurrency) 
  {
    TimeSolver3 ts = new TimeSolver3(n1,n2,n3,tensors);
    ts.setConcurrency(concurrency);
    Stopwatch sw = new Stopwatch();
    sw.start();
    ts.zeroAt(i1,i2,i3);
    sw.stop();
    float[][][] t = ts.getTimes();
    trace("  time="+sw.time()+" sum="+Array.sum(t));
    return t;
  }

  private static void testConstant() {
    //int n1 = 101, n2 = 101, n3 = 101;
    int n1 = 64, n2 = 64, n3 = 64;
    float d11 = 1.000f, d12 = 0.900f, d13 = 0.900f,
                        d22 = 1.000f, d23 = 0.900f,
                                      d33 = 1.000f;
    ConstantTensors tensors = new ConstantTensors(d11,d12,d13,d22,d23,d33);
    int i1 = 2*(n1-1)/4, i2 = 2*(n2-1)/4, i3 = 2*(n3-1)/4;
    float[][][] ts = computeSerial(n1,n2,n3,i1,i2,i3,tensors);
    float[][][] tp = computeParallel(n1,n2,n3,i1,i2,i3,tensors);
    float[][][] te = Array.div(Array.abs(Array.sub(tp,ts)),ts);
    te[i3][i2][i1] = 0.0f;
    float temax = Array.max(te);
    trace("temax="+temax);
    trace("********************************************************");
    //plot(ts,ColorMap.PRISM);
    //plot(tp,ColorMap.PRISM);
    //plot(te,ColorMap.JET);
    if (temax>0.1f)
      System.exit(-1);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
  private static void trace(int i1, int i2, int i3, String s) {
    if (i1==0 && i2==2 && i3==1 ||
        i1==2 && i2==0 && i3==1)
      trace("i1="+i1+" i2="+i2+" i3="+i3+": "+s);
  }

  public static void main(String[] args) {
    //SwingUtilities.invokeLater(new Runnable() {
    //  public void run() {
        for (;;)
          testConstant();
    //  }
    //});
  }
}
