/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.Tensors2;
import static edu.mines.jtk.util.MathPlus.min;
import static edu.mines.jtk.util.MathPlus.sqrt;

/**
 * A time transform for a 2D anisotropic eikonal equation. The non-linear 
 * equation is grad(t) dot W*grad(t) = 1, where t is the solution time 
 * map and W denotes a positive-definite (velocity-squared) metric tensor 
 * field.
 * <p>
 * Times in a time map are flagged as either known or unknown. Times for 
 * known samples are typically zero (or near zero), and are never modified. 
 * Times for unknown samples are computed to be solutions to the eikonal 
 * equation. Each of these computed times represents the minimum time
 * to travel from the unknown sample to one of the known samples.
 * <p>
 * A separate map of marks may be computed as times are computed. The 
 * mark of each unknown sample equals the mark of that known sample for 
 * which the time between the two samples is minimized. For an isotropic 
 * homogeneous tensor field, time equals Euclidean distance, and a map of 
 * marks that represents a discrete Voronoi diagram can be computed by 
 * assigning each known sample a unique mark.
 * <p>
 * This transform uses an iterative sweeping method to compute the time map.
 * Iterations are similar to those described by Tsai et al., 2002.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.12.24
 */
public class TimeMapper2 {

  /**
   * Type of concurrency used by this transform.
   */
  public enum Concurrency {
    PARALLEL,
    SERIAL
  };
  
  /**
   * Constructs a time mapper for the specified tensor field.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors velocity-squared tensors.
   */
  public TimeMapper2(int n1, int n2, Tensors2 tensors) {
    init(n1,n2,tensors);
  }

  /**
   * Sets the tensors used by this time mapper.
   * @param tensors the tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = tensors;
  }

  /**
   * Applys this transform to the specified array of times and marks.
   * @param known array of flags; zero for unknown samples
   * @param times array of times in which unknown times are to be computed.
   * @param marks array of marks in which unknown marks are to be computed.
   */
  public void apply(byte[][] known, float[][] times, int[][] marks) {

    // Initialize unknown times to infinity.
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        //if (known[i2][i1]==0) {
          times[i2][i1] = INFINITY;
        //}
      }
    }

    // Indices of known samples in random order.
    short[][] kk = indexKnownSamples(known);
    short[] k1 = kk[0];
    short[] k2 = kk[1];
    shuffle(k1,k2);
    int nk = k1.length;

    // Active list of samples used to compute times.
    ActiveList al = new ActiveList();

    // For all known samples, ...
    for (int ik=0; ik<nk; ++ik) {
      int i1 = k1[ik];
      int i2 = k2[ik];
      //trace("processing known sample at i1="+i1+" i2="+i2+
      //      " mark="+marks[i2][i1]);
      times[i2][i1] = 0.0f;

      // Clear activated flags so we can tell which samples become activated.
      clearActivated();

      // Put the known sample into the active list.
      al.append(_s[i2][i1]);

      // Process the active list until empty.
      solve(times,al);

      // Mark samples that were activated while solving for times.
      markActivated(i1,i2,known,marks);
      //plot(marks);
    }
  }

  private void solve(float[][] times, ActiveList al) {
    if (_concurrency==Concurrency.PARALLEL) {
      solveParallel(times,al);
    } else {
      solveSerial(times,al);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Default time for samples not yet computed.
  private static final float INFINITY = Float.MAX_VALUE;

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.01f;
  private static final float ONE_MINUS_EPSILON = 1.0f-EPSILON;

  private int _n1,_n2;
  private Tensors2 _tensors;
  private Sample[][] _s;
  private Concurrency _concurrency = Concurrency.PARALLEL;
  private ArrayList<Sample> _stack = new ArrayList<Sample>(1024);

  private void init(int n1, int n2, Tensors2 tensors) {
    _n1 = n1;
    _n2 = n2;
    _tensors = tensors;
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
    int activated; // used to flag activated samples
    boolean absent; // used to build active lists
    Sample(int i1, int i2) {
      this.i1 = i1;
      this.i2 = i2;
    }
  }

  // List of active samples.
  private class ActiveList {
    void append(Sample s) {
      s.activated = _activated;
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

  // Flags set during computation of times. For efficiency, do not
  // loop over all the flags to clear them before computing times.
  // Instead, modify the value that represents activated samples.
  private int _activated = 1;
  private void clearActivated() {
    if (_activated==Integer.MAX_VALUE) { // rarely!
      _activated = 1;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1].activated = 0;
        }
      }
    } else { // typically
      ++_activated;
    }
  }
  private void setActivated(Sample s) {
    s.activated = _activated;
  }
  private void clearActivated(Sample s) {
    s.activated = 0;
  }
  private boolean wasActivated(Sample s) {
    return s.activated==_activated;
  }

  /**
   * Returns arrays of indices of known samples.
   * Includes only known samples adjacent to at least one unknown sample.
   * (Does not include known samples surrounded by other known samples.)
   */
  private short[][] indexKnownSamples(byte[][] known) {
    ShortStack ss1 = new ShortStack();
    ShortStack ss2 = new ShortStack();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (known[i2][i1]!=0) {
          for (int k=0; k<4; ++k) {
            int k1 = K1[k];
            int k2 = K2[k];
            int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
            int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
            if (known[j2][j1]==0) {
              ss1.push(i1);
              ss2.push(i2);
              break;
            }
          }
        }
      }
    }
    short[] i1 = ss1.array();
    short[] i2 = ss2.array();
    return new short[][]{i1,i2};
  }

  /**
   * Randomly (but consistently) shuffles the specified arrays of indices.
   */
  private static void shuffle(short[] i1, short[] i2) {
    int n = i1.length;
    Random r = new Random(314159);
    short ii;
    for (int i=n-1; i>0; --i) {
      int j = r.nextInt(i+1);
      ii = i1[i]; i1[i] = i1[j]; i1[j] = ii;
      ii = i2[i]; i2[i] = i2[j]; i2[j] = ii;
    }
  }

  // More efficient than ArrayStack<Short>.
  private static class ShortStack {
    void push(int k) {
      if (_n==_a.length) {
        short[] a = new short[2*_n];
        System.arraycopy(_a,0,a,0,_n);
        _a = a;
      }
      _a[_n++] = (short)k;
    }
    short pop() {
      return _a[--_n];
    }
    int size() {
      return _n;
    }
    void clear() {
      _n = 0;
    }
    boolean isEmpty() {
      return _n==0;
    }
    short[] array() {
      short[] a = new short[_n];
      System.arraycopy(_a,0,a,0,_n);
      return a;
    }
    private int _n = 0;
    private short[] _a = new short[2048];
  }

  /**
   * Initializes all unknown times to infinity.
   */
  private void init(byte[][] known, float[][] times) {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (known[i2][i1]==0) {
          times[i2][i1] = INFINITY;
        }
      }
    }
  }

  /**
   * Marks samples that were activated during computation of times from
   * a known sample with specified indices. Samples to mark are flagged
   * as activated. Activation flags are cleared as samples are marked.
   */
  private void markActivated(int i1, int i2, byte[][] known, int[][] marks) {
    ShortStack ss1 = new ShortStack();
    ShortStack ss2 = new ShortStack();
    ss1.push(i1);
    ss2.push(i2);
    Sample si = _s[i2][i1];
    int mark = marks[i2][i1];
    clearActivated(si);
    int nmark = 0;
    while (!ss1.isEmpty()) {
      i1 = ss1.pop();
      i2 = ss2.pop();
      if (known[i2][i1]==0)
        marks[i2][i1] = mark;
      ++nmark;
      for (int k=0; k<4; ++k) {
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
        Sample sj = _s[j2][j1];
        if (wasActivated(sj)) {
          ss1.push(j1);
          ss2.push(j2);
          clearActivated(sj);
        }
      }
    }
  }

  /**
   * Solves for times by sequentially processing each sample in active list.
   */
  private void solveSerial(float[][] times, ActiveList al) {
    float[] d = new float[3];
    ActiveList bl = new ActiveList();
    int ntotal = 0;
    while (!al.isEmpty()) {
      //al.shuffle(); // demonstrate that solution depends on order
      int n = al.size();
      ntotal += n;
      for (int i=0; i<n; ++i) {
        Sample s = al.get(i);
        solveOne(times,s,bl,d);
      }
      bl.setAllAbsent();
      al.clear();
      al.appendIfAbsent(bl);
      bl.clear();
    }
    //trace("solveSerial: ntotal="+ntotal);
    //trace("             nratio="+(float)ntotal/(float)(_n1*_n2));
  }
  
  /**
   * Solves for times by processing samples in the active list in parallel.
   */
  private void solveParallel(final float[][] times, final ActiveList al) {
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
    int niter = 0;
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
                solveOne(times,s,bltask,dtask); // process the sample
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
      ++niter;
      //if (niter%100==1)
      //  plot(_t,ColorMap.JET);
    }
    es.shutdown();
    //trace("solveParallel: ntotal="+ntotal);
    //trace("               nratio="+(float)ntotal/(float)(_n1*_n2));
  }

  /**
   * Processes one sample from the A list.
   * Appends samples not yet converged to the B list.
   */
  private void solveOne(float[][] times, Sample s, ActiveList bl, float[] d) {

    // Sample indices.
    int i1 = s.i1;
    int i2 = s.i2;

    // Current time and new time computed from all neighbors.
    float ti = times[i2][i1];
    float ci = computeTime(times,i1,i2,K1S[4],K2S[4],d);
    times[i2][i1] = ci;

    // If new and current times are close enough (converged), then ...
    if (ci>=ti*ONE_MINUS_EPSILON) {

      // For all four neighbors, ...
      for (int k=0; k<4; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;

        // Compute time for neighbor.
        float tj = times[j2][j1];
        float cj = computeTime(times,j1,j2,K1S[k],K2S[k],d);

        // If computed time significantly less than neighbor's current time, ...
        if (cj<tj*ONE_MINUS_EPSILON) {

          // Replace the current time.
          times[j2][j1] = cj;
          
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

  // Methods to get times for neighbors.
  /*
  private float t1m(float[][] t, int i1, int i2) {
    return (--i1>=0 && wasActivated(_s[i2][i1]))?t[i2][i1]:INFINITY;
  }
  private float t1p(float[][] t, int i1, int i2) {
    return (++i1<_n1 && wasActivated(_s[i2][i1]))?t[i2][i1]:INFINITY;
  }
  private float t2m(float[][] t, int i1, int i2) {
    return (--i2>=0 && wasActivated(_s[i2][i1]))?t[i2][i1]:INFINITY;
  }
  private float t2p(float[][] t, int i1, int i2) {
    return (++i2<_n2 && wasActivated(_s[i2][i1]))?t[i2][i1]:INFINITY;
  }
  */
  private float t1m(float[][] t, int i1, int i2) {
    return (--i1>=0)?t[i2][i1]:INFINITY;
  }
  private float t1p(float[][] t, int i1, int i2) {
    return (++i1<_n1)?t[i2][i1]:INFINITY;
  }
  private float t2m(float[][] t, int i1, int i2) {
    return (--i2>=0)?t[i2][i1]:INFINITY;
  }
  private float t2p(float[][] t, int i1, int i2) {
    return (++i2<_n2)?t[i2][i1]:INFINITY;
  }

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor samples with specified indices.
   */
  private float computeTime(
    float[][] t, int i1, int i2, int[] k1s, int[] k2s, float[] d) 
  {
    _tensors.getTensor(i1,i2,d);
    float d11 = d[0];
    float d12 = d[1];
    float d22 = d[2];
    float e12 = 1.0f/(d11*d22-d12*d12);
    float tc = t[i2][i1];
    float t1m = t1m(t,i1,i2);
    float t1p = t1p(t,i1,i2);
    float t2m = t2m(t,i1,i2);
    float t2p = t2p(t,i1,i2);
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

  private static void trace(String s) {
    System.out.println(s);
  }

  private static float[][] toFloat(byte[][] b) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        f[i2][i1] = (float)b[i2][i1];
    return f;
  }
  private static float[][] toFloat(int[][] i) {
    int n1 = i[0].length;
    int n2 = i.length;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        f[i2][i1] = (float)i[i2][i1];
    return f;
  }
  private static void plot(byte[][] b) {
    plot(toFloat(b));
  }
  private static void plot(int[][] i) {
    plot(toFloat(i));
  }
  private static void plot(float[][] f) {
    edu.mines.jtk.mosaic.SimplePlot sp =
      new edu.mines.jtk.mosaic.SimplePlot(
        edu.mines.jtk.mosaic.SimplePlot.Origin.UPPER_LEFT);
    sp.setSize(920,900);
    edu.mines.jtk.mosaic.PixelsView pv = sp.addPixels(f);
    pv.setColorModel(edu.mines.jtk.awt.ColorMap.JET);
    pv.setInterpolation(edu.mines.jtk.mosaic.PixelsView.Interpolation.NEAREST);
  }
}
