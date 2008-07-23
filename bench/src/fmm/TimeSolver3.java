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
 * A solver for 3D anisotropic eikonal equations.
 * These non-linear equations are sqrt(grad(t) dot D*grad(t)) = 1, where t
 * denotes the solution time, and D denotes a diffusion (velocity-squared)
 * tensor field.
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
   * An interface for classes of diffusion tensors. Each tensor is a
   * symmetric positive-definite 3-by-3 matrix 
   * {{d11,d12,d13},{d12,d22,d23},{d13,d23,d33}}.
   */
  public interface Tensors {

    /**
     * Gets diffusion tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param i3 index for 3rd dimension.
     * @param d array {d11,d12,d13,d22,d23,d33} of tensor elements.
     */
    public void getTensor(int i1, int i2, int i3, float[] d);
  }

  /**
   * Constructs a solver with constant identity diffusion tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public TimeSolver3(int n1, int n2, int n3) {
    this(n1,n2,n3,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified diffusion tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param tensors diffusion tensors.
   */
  public TimeSolver3(int n1, int n2, int n3, Tensors tensors) {
    init(n1,n2,n3,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors diffusion tensors.
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
   * Zeros the time at the specified sample and updates times elsewhere.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @param i3 index in 3rd dimension of time to zero.
   * @return the updated array of times; by reference, not by copy.
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

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Default time for samples not yet computed.
  private static final float INFINITY = Float.MAX_VALUE;

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.001f;

  private int _n1,_n2,_n3;
  private Tensors _tensors;
  private float[][][] _t;
  private Sample[][][] _s;
  private Concurrency _concurrency = Concurrency.PARALLEL;

  private void init(int n1, int n2, int n3, float[][][] t, Tensors tensors) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
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
      Sample s = _a[i];
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
        }
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
  }

  /**
   * Solves for times by sequentially processing each sample in active list.
   */
  private void solveSerial(ActiveList al) {
    float[] d = new float[6];
    ActiveList bl = new ActiveList();
    int ntotal = 0;
    while (!al.isEmpty()) {
      int n = al.size();
      ntotal += n;
      for (int i=0; i<n; ++i)
        solveOne(i,al,bl,d);
      ActiveList tl = al; 
      al = bl; 
      bl = tl; 
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
    // Benchmarks: 07/22/2008
    // Intel 2.4 GHz Core 2 Duo for size 101 * 101 * 101
    // serial         5.7 s
    //nthread = 1; // 5.8 s
    //nthread = 2; // 3.2 s
    // Intel 3.0 GHz 2 * Quad Core Xeon XXX for size 101 * 101 * 101
    //serial          4.1 s
    //nthread = 1; // 4.3 s
    //nthread = 4; // 1.3 s
    //nthread = 8; // 0.7 s
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
    trace("solveParallel: ntotal="+ntotal);
    trace("               nratio="+(float)ntotal/(float)(_n1*_n2*_n3));
  }

  /**
   * Processes one sample in the A list.
   * Appends samples not yet converged to the B list.
   */
  private void solveOne(int i, ActiveList al, ActiveList bl, float[] d) {

    // Get one sample from the A list.
    Sample si = al.get(i);
    int i1 = si.i1;
    int i2 = si.i2;
    int i3 = si.i3;

    // Current time and new time computed from all neighbors.
    float ti = _t[i3][i2][i1];
    float gi = g(i1,i2,i3,K1S[6],K2S[6],K3S[6],d);
    _t[i3][i2][i1] = gi;

    // If new and current times are close enough (converged), then ...
    if (ti-gi<=ti*EPSILON) {

      // For all six neighbor samples, ...
      for (int k=0; k<6; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
        int j3 = i3+K3[k];  if (j3<0 || j3>=_n3) continue;

        // Compute time for the neighbor.
        float tj = _t[j3][j2][j1];
        float gj = g(j1,j2,j3,K1S[k],K2S[k],K3S[k],d);

        // If computed time less than the neighbor's current time, ...
        if (tj-gj>tj*EPSILON) {

          // Replace the current time.
          _t[j3][j2][j1] = gj;
          
          // Append neighbor sample to the B list.
          Sample sj = _s[j3][j2][j1];
          bl.append(sj);
        }
      }
    }

    // Else, if not converged, append this sample to the B list.
    else {
      bl.append(si);
    }
  }

  /**
   * Solves a quadratic equation for a positive time t0.
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
   * It then returns t0 = t1+u. If no solution exists, because the 
   * discriminant is negative, this method returns INFINITY.
   */
  private static float solveQuadratic(
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
    double u = (-b+sqrt(d))/(2.0*a);
    return t1+(float)u; // t0 = t1+u
  }

  /**
   * Jeong's fast tests for a valid solution time t0 to H(p1,p2,p3) = 1.
   * Parameters tm and tp are times for samples backward and forward of 
   * the sample with time t0. The parameter k is the index k1, k2, or k3
   * that was used to compute the time, and the parameter p is a critical
   * point of H(p1,p2,p3) for fixed p1, p2, or p3.
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
  private boolean isValid1(int i1, int i2, int i3, int k1, float p1, float t0) {
    float tm = (i1>0    )?_t[i3][i2][i1-1]:INFINITY;
    float tp = (i1<_n1-1)?_t[i3][i2][i1+1]:INFINITY;
    return isValid(tm,tp,t0,k1,p1);
  }
  private boolean isValid2(int i1, int i2, int i3, int k2, float p2, float t0) {
    float tm = (i2>0    )?_t[i3][i2-1][i1]:INFINITY;
    float tp = (i2<_n2-1)?_t[i3][i2+1][i1]:INFINITY;
    return isValid(tm,tp,t0,k2,p2);
  }
  private boolean isValid3(int i1, int i2, int i3, int k3, float p3, float t0) {
    float tm = (i3>0    )?_t[i3-1][i2][i1]:INFINITY;
    float tp = (i3<_n3-1)?_t[i3+1][i2][i1]:INFINITY;
    return isValid(tm,tp,t0,k3,p3);
  }

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor samples with specified indices.
   */
  private float g(
    int i1, int i2, int i3, int[] k1s, int[] k2s, int[] k3s, float[] d) 
  {
    float tc = _t[i3][i2][i1];

    // Get tensor coefficients.
    _tensors.getTensor(i1,i2,i3,d);
    float d11 = d[0];
    float d12 = d[1];
    float d13 = d[2];
    float d22 = d[3];
    float d23 = d[4];
    float d33 = d[5];

    // For all relevant neighbor samples, ...
    for (int k=0; k<k1s.length; ++k) {
      int k1 = k1s[k];
      int k2 = k2s[k];
      int k3 = k3s[k];

      // (p1-,p2s,p3s), (p1+,p2s,p3s)
      if (k1!=0 && k2==0 && k3==0) {
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        float t1 = _t[i3][i2][j1];  if (t1==INFINITY) continue;
        float t2 = t1;
        float t3 = t1;
        float s1 = k1;
        float ddet = d22*d33-d23*d23;
        float s2 = (d23*d13-d12*d33)*s1/ddet;
        float s3 = (d23*d12-d13*d22)*s1/ddet;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=t1) {
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          if (isValid1(i1,i2,i3,k1,p1,t0)) {
            return t0;
          }
        }
      }

      // (p1s,p2-,p3s), (p1s,p2-,p3s)
      else if (k1==0 && k2!=0 && k3==0) {
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        float t2 = _t[i3][j2][i1];  if (t2==INFINITY) continue;
        float t1 = t2;
        float t3 = t2;
        float s2 = k2;
        float ddet = d11*d33-d13*d13;
        float s1 = (d13*d23-d12*d33)*s2/ddet;
        float s3 = (d13*d12-d23*d11)*s2/ddet;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=t2) {
          float t01 = t0-t1;
          float t03 = t0-t3;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          if (isValid2(i1,i2,i3,k2,p2,t0)) {
            return t0;
          }
        }
      }

      // (p1s,p2s,p3-), (p1s,p2s,p3-)
      else if (k1==0 && k2==0 && k3!=0) {
        int j3 = i3+k3;  if (j3<0 || j3>=_n3) continue;
        float t3 = _t[j3][i2][i1];  if (t3==INFINITY) continue;
        float t1 = t3;
        float t2 = t3;
        float s3 = k3;
        float ddet = d11*d22-d12*d12;
        float s1 = (d12*d23-d13*d22)*s3/ddet;
        float s2 = (d12*d13-d23*d11)*s3/ddet;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=t3) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid3(i1,i2,i3,k3,p3,t0)) {
            return t0;
          }
        }
      }

      // (p1s,p2-,p3-), (p1s,p2+,p3-), (p1s,p2+,p3-), (p1s,p2+,p3+)
      else if (k1==0 && k2!=0 && k3!=0) {
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        int j3 = i3+k3;  if (j3<0 || j3>=_n3) continue;
        float t2 = _t[i3][j2][i1];  if (t2==INFINITY) continue;
        float t3 = _t[j3][i2][i1];  if (t3==INFINITY) continue;
        float s2 = k2;
        float s3 = k3;
        float ds12 = d12*s2;
        float ds13 = d13*s3;
        float dnum = ds12*t2+ds13*t3;
        float dden = ds12+ds13;
        float t1;
        if (dden==0.0f) {
          if (dnum==0.0f)
            t1 = 0.5f*(t2+t3);
          else
            continue;
        } else {
          t1 = dnum/dden;
        }
        float s1 = -dden/d11;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=min(t2,t3)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid2(i1,i2,i3,k2,p2,t0) &&
              isValid3(i1,i2,i3,k3,p3,t0)) {
            return t0;
          }
        }
      }

      // (p1-,p2s,p3-), (p1+,p2s,p3-), (p1-,p2s,p3-), (p1+,p2s,p3+)
      else if (k1!=0 && k2==0 && k3!=0) {
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        int j3 = i3+k3;  if (j3<0 || j3>=_n3) continue;
        float t1 = _t[i3][i2][j1];  if (t1==INFINITY) continue;
        float t3 = _t[j3][i2][i1];  if (t3==INFINITY) continue;
        float s1 = k1;
        float s3 = k3;
        float ds12 = d12*s1;
        float ds23 = d23*s3;
        float dnum = ds12*t1+ds23*t3;
        float dden = ds12+ds23;
        float t2;
        if (dden==0.0f) {
          if (dnum==0.0f)
            t2 = 0.5f*(t1+t3);
          else
            continue;
        } else {
          t2 = dnum/dden;
        }
        float s2 = -dden/d22;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=min(t1,t3)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid1(i1,i2,i3,k1,p1,t0) &&
              isValid3(i1,i2,i3,k3,p3,t0)) {
            return t0;
          }
        }
      }

      // (p1-,p2-,p3s), (p1+,p2-,p3s), (p1-,p2+,p3s), (p1+,p2+,p3s)
      else if (k1!=0 && k2!=0 && k3==0) {
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        float t1 = _t[i3][i2][j1];  if (t1==INFINITY) continue;
        float t2 = _t[i3][j2][i1];  if (t2==INFINITY) continue;
        float s1 = k1;
        float s2 = k2;
        float ds13 = d13*s1;
        float ds23 = d23*s2;
        float dnum = ds13*t1+ds23*t2;
        float dden = ds13+ds23;
        float t3;
        if (dden==0.0f) {
          if (dnum==0.0f)
            t3 = 0.5f*(t1+t2);
          else
            continue;
        } else {
          t3 = dnum/dden;
        }
        float s3 = -dden/d33;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=min(t1,t2)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          if (isValid1(i1,i2,i3,k1,p1,t0) &&
              isValid2(i1,i2,i3,k2,p2,t0)) {
            return t0;
          }
        }
      }
      
      // (p1-,p2-,p3-), (p1+,p2-,p3-), (p1-,p2+,p3-), (p1+,p2+,p3-), 
      // (p1-,p2-,p3+), (p1+,p2-,p3+), (p1-,p2+,p3+), (p1+,p2+,p3+)
      else {
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        int j3 = i3+k3;  if (j3<0 || j3>=_n3) continue;
        float t1 = _t[i3][i2][j1];  if (t1==INFINITY) continue;
        float t2 = _t[i3][j2][i1];  if (t2==INFINITY) continue;
        float t3 = _t[j3][i2][i1];  if (t3==INFINITY) continue;
        float s1 = k1;
        float s2 = k2;
        float s3 = k3;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tc && t0>=min(t1,t2,t3)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid1(i1,i2,i3,k1,p1,t0) &&
              isValid2(i1,i2,i3,k2,p2,t0) &&
              isValid3(i1,i2,i3,k3,p3,t0)) {
            return t0;
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
    trace("********************************************************");
    int n1 = 101;
    int n2 = 101;
    int n3 = 101;
    //float s11 = 1.000f, s12 = 0.000f, s13 = 0.000f,
    //                    s22 = 1.000f, s23 = 0.000f,
    //                                  s33 = 1.000f;
    float s11 = 1.000f, s12 = 0.900f, s13 = 0.900f,
                        s22 = 1.000f, s23 = 0.900f,
                                      s33 = 1.000f;
    ConstantTensors tensors = new ConstantTensors(s11,s12,s13,s22,s23,s33);
    int i1 = 2*(n1-1)/4;
    int i2 = 2*(n2-1)/4;
    int i3 = 2*(n3-1)/4;
    float[][][] ts = computeSerial(n1,n2,n3,i1,i2,i3,tensors);
    float[][][] tp = computeParallel(n1,n2,n3,i1,i2,i3,tensors);
    float[][][] te = Array.div(Array.abs(Array.sub(tp,ts)),ts);
    te[i3][i2][i1] = 0.0f;
    float temax = Array.max(te);
    trace("temax="+temax);
    if (temax>0.1f)
      System.exit(-1);
    trace("********************************************************");
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
