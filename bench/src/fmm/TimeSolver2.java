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
 * A solver for 2D anisotropic eikonal equations.
 * These non-linear equations are sqrt(grad(t) dot D*grad(t)) = 1, where t
 * denotes the solution time, and D denotes a diffusion (velocity-squared)
 * tensor field.
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
   * An interface for classes of diffusion tensors. Each tensor is a
   * symmetric positive-definite 2-by-2 matrix {{d11,d12},{d12,d22}}.
   */
  public interface Tensors {

    /**
     * Gets diffusion tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param d array {d11,d12,d22} of tensor elements.
     */
    public void getTensor(int i1, int i2, float[] d);
  }

  /**
   * Constructs a solver with constant identity diffusion tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public TimeSolver2(int n1, int n2) {
    this(n1,n2,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified diffusion tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors diffusion tensors.
   */
  public TimeSolver2(int n1, int n2, Tensors tensors) {
    init(n1,n2,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors diffusion tensors.
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
  private static final float EPSILON = 0.0001f;

  private int _n1,_n2;
  private Tensors _tensors;
  private float[][] _t;
  private Sample[][] _s;
  private int _active = 0;
  private Concurrency _concurrency = Concurrency.PARALLEL;

  private void init(int n1, int n2, float[][] t, Tensors tensors) {
    _n1 = n1;
    _n2 = n2;
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

  // A sample has indices and is either active or inactive.
  // For efficiency the active flag is an integer and not a boolean, 
  // so that we need not loop over all samples when initializing them 
  // to be inactive. See the comments for the method clearActive below.
  private static class Sample {
    int i1,i2; // sample indices
    int active; // determines whether this sample is active
    boolean absent; // determines whether this sample is in a list
    Sample(int i1, int i2) {
      this.i1 = i1;
      this.i2 = i2;
    }
  }

  // Returns true if specified sample is active; false, otherwise.
  private boolean isActive(int i1, int i2) {
    return _s[i2][i1].active==_active;
  }

  // Marks all samples inactive. For efficiency, we typically do not loop 
  // over all the samples to clear their active flags. Usually we simply
  // increment the active value with which the flags are compared.
  private void clearActive() {
    if (_active==Integer.MAX_VALUE) { // rarely!
      _active = 1;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1].active = 0;
        }
      }
    } else { // typically
      ++_active;
    }
  }

  // List of active samples.
  private class ActiveList {
    void append(Sample s) {
      if (_n==_a.length)
        growTo(2*_n);
      _a[_n++] = s;
      s.active = _active;
    }
    boolean isEmpty() {
      return _n==0;
    }
    int size() {
      return _n;
    }
    Sample get(int i) {
      Sample s = _a[i];
      assert s.active==_active:"sample in list is active";
      return s;
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
          s.active = _active;
        }
      }
    }
    private int _n;
    private Sample[] _a = new Sample[1024];
    private void growTo(int capacity) {
      //trace("SampleList: growing to capacity="+capacity);
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

    // All samples initially inactive.
    clearActive();

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
    while (!al.isEmpty()) {
      int n = al.size();
      for (int i=0; i<n; ++i)
        solveOne(i,al,bl,d);
      ActiveList tl = al; 
      al = bl; 
      bl = tl; 
      bl.clear();
    }
  }
  
  /**
   * Solves for times by processing samples in the active list in parallel.
   */
  private void solveParallel(final ActiveList al) {
    int nthread = Runtime.getRuntime().availableProcessors();
    /////////////////////////////////////////////////////////////////////////
    // Benchmarks: 07/22/2008
    // Intel 2.4 GHz Core 2 Duo for size 2001*2001
    // serial         5.1 s
    //nthread = 1; // 5.0 s
    //nthread = 2; // 2.8 s
    // Intel 2.4 GHz Core 2 Duo for size 4001*4001
    // serial         23.5 s (peak %CPU = 100)
    //nthread = 1; // 25.5 s (peak %CPU = 100)
    //nthread = 2; // 14.8 s (peak %CPU = 180)
    // Intel 3.0 GHz 2 * Quad Core Xeon ??? for size 2001*2001
    //serial          2.9 s
    //nthread = 1; // 3.5 s
    //nthread = 4; // 1.7 s
    //nthread = 8; // 1.2 s
    // Intel 3.0 GHz 2 * Quad Core Xeon ??? for size 4001*4001
    //serial          12.4 s (peak %CPU = 100)
    //nthread = 1; // 17.5 s (peak %CPU = 100)
    //nthread = 4; //  6.9 s (peak %CPU = 350)
    //nthread = 8; //  4.5 s (peak %CPU = 610)
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
    while (!al.isEmpty()) {
      ai.set(0); // initialize the shared block index to zero
      final int n = al.size(); // number of samples in active (A) list
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
              for (int k=i; k<j; ++k) // for each sample in block
                solveOne(k,al,bltask,dtask); // process sample in active list 
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

      // Merge samples from all B lists to a new A list. Ensure that 
      // a sample is appended no more than once to the new A list.
      al.clear();
      for (int itask=0; itask<ntask; ++itask) {
        al.appendIfAbsent(bl[itask]);
        bl[itask].clear();
      }
    }
    es.shutdown();
  }

  /**
   * Processes one sample in the active list al.
   * May append samples to the active list bl.
   */
  private void solveOne(int i, ActiveList al, ActiveList bl, float[] d) {

    // Get one sample from active list A.
    Sample si = al.get(i);
    int i1 = si.i1;
    int i2 = si.i2;

    // Current time and new time computed from all neighbors.
    float ti = _t[i2][i1];
    float gi = g(i1,i2,K1S[4],K2S[4],d);
    _t[i2][i1] = gi;

    // If new and current times are close enough (converged), then ...
    if (ti-gi<=ti*EPSILON) {

      // This sample is now inactive (but may be reactivated).
      si.active -= 1;

      // For all four neighbor samples, ...
      for (int k=0; k<4; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;

        // If neighbor is not active (not in lists A or B), ...
        if (!isActive(j1,j2)) {

          // Compute time for the neighbor.
          float gj = g(j1,j2,K1S[k],K2S[k],d);

          // If computed time less than the neighbor's current time, ...
          if (gj<_t[j2][j1]) {

            // Replace the current time.
            _t[j2][j1] = gj;
            
            // Append neighbor sample to the active list B.
            Sample sj = _s[j2][j1];
            bl.append(sj);
          }
        }
      }
    }

    // Else, if not converged, append this sample to the active list B.
    else {
      bl.append(si);
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
   * Jeong's fast tests for a valid solution time t0 to H(p1,p2) = 1.
   * Parameters tm and tp are times for samples backward and forward of 
   * the sample with time t0. The parameter k is the index k1 or k2 that 
   * was used to compute the time, and the parameter p is a critical
   * point of H(p1,p2) for fixed p1 or p2.
   */
  private static boolean isValid(
    float tm, float tp, float t0, int k, float p) 
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
  private boolean isValid1(int i1, int i2, int k1, float p1, float t0) {
    float tm = (i1>0    )?_t[i2][i1-1]:INFINITY;
    float tp = (i1<_n1-1)?_t[i2][i1+1]:INFINITY;
    return isValid(tm,tp,t0,k1,p1);
  }
  private boolean isValid2(int i1, int i2, int k2, float p2, float t0) {
    float tm = (i2>0    )?_t[i2-1][i1]:INFINITY;
    float tp = (i2<_n2-1)?_t[i2+1][i1]:INFINITY;
    return isValid(tm,tp,t0,k2,p2);
  }

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor samples with specified indices.
   */
  private float g(int i1, int i2, int[] k1s, int[] k2s, float[] d) {
    float tc = _t[i2][i1];

    // Get tensor coefficients.
    //float[] d = new float[3];
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
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        float t2 = _t[j2][i1];
        if (t2!=INFINITY) {
          float t1 = t2;
          float s2 = k2;
          float s1 = -s2*d12/d11;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tc && t0>=t2) {
            float p2 = -s1*(t1-t0)*d12/d22;
            if (isValid2(i1,i2,k2,p2,t0)) {
              return t0;
            }
          }
        }
      } 
      
      // else, if (p1-,p2s) or (p1+,p2s), ...
      else if (k2==0) {
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        float t1 = _t[i2][j1];
        if (t1!=INFINITY) {
          float t2 = t1;
          float s1 = k1;
          float s2 = -s1*d12/d22;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tc && t0>=t1) {
            float p1 = -s2*(t2-t0)*d12/d11;
            if (isValid1(i1,i2,k1,p1,t0)) {
              return t0;
            }
          }
        }
      } 
      
      // else, if (p1-,p2-), (p1+,p2-), (p1-,p2+) or (p1+,p2+), ...
      else {
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        float t1 = _t[i2][j1];
        float t2 = _t[j2][i1];
        if (t1!=INFINITY && t2!=INFINITY) {
          float s1 = k1;
          float s2 = k2;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tc && t0>=min(t1,t2)) {
            float p1 = -s2*(t2-t0)*d12/d11;
            float p2 = -s1*(t1-t0)*d12/d22;
            if (isValid1(i1,i2,k1,p1,t0) &&
                isValid2(i1,i2,k2,p2,t0)) {
              return t0;
            }
          }
        }
      }
    }

    return tc;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

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

  private static void testConstant() {
    int n1 = 2001;
    int n2 = 2001;
    float angle = FLT_PI*110.0f/180.0f;
    float su = 1.000f;
    float sv = 0.010f;
    //float sv = 1.000f;
    float cosa = cos(angle);
    float sina = sin(angle);
    float d11 = su*cosa*cosa+sv*sina*sina;
    float d12 = (su-sv)*sina*cosa;
    float d22 = sv*cosa*cosa+su*sina*sina;
    ConstantTensors dt = new ConstantTensors(d11,d12,d22);
    TimeSolver2 ts = new TimeSolver2(n1,n2,dt);
    ts.setConcurrency(TimeSolver2.Concurrency.PARALLEL);
    //ts.setConcurrency(TimeSolver2.Concurrency.SERIAL);
    Stopwatch sw = new Stopwatch();
    sw.start();
    ts.zeroAt(2*n1/4,2*n2/4);
    //ts.zeroAt(1*n1/4,1*n2/4);
    //ts.zeroAt(3*n1/4,3*n2/4);
    sw.stop();
    float[][] t = ts.getTimes();
    trace("time="+sw.time()+" sum="+Array.sum(t));
    //Array.dump(t);
    //plot(t);
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
      for (;;)
        testConstant();
    //  }
    //});
  }
}
