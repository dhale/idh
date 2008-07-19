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
 * 3D implementation of Jeong and Whitakers' fast iterative method.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.16
 */
public class FimSolver3 {

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
  public FimSolver3(int n1, int n2, int n3) {
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
  public FimSolver3(int n1, int n2, int n3, Tensors tensors) {
    init(n1,n2,n3,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors diffusion tensors.
   */
  public FimSolver3(float[][][] t, Tensors tensors) {
    init(t[0][0].length,t[0].length,t.length,t,tensors);
  }

  /**
   * Sets the type of concurrency used by this solver.
   * The default concurrency is parallel.
   * @param parallel true, for parallel; false, for serial.
   */
  public void setParallel(boolean parallel) {
    _parallel = parallel;
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
  private static final float EPSILON = 0.01f;

  private int _n1,_n2,_n3;
  private Tensors _tensors;
  private float[][][] _t;
  private Sample[][][] _s;
  private int _active = 0;
  private boolean _parallel = true;

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
  // neighbor with offsets {K1[0],K2[0],K3[0]} = {-1,0,0}, only the
  // sets K1S[0], K2S[0], and K3S[0] are used. The sets K1S[6], K2S[6],
  // and K3S[6] are special offsets for all six neighbors.
  private static final int[][] K1S = {
    { 1, 1, 1, 1, 1, 1, 1, 1, 1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1}, // A
    { 0,-1, 1,-1, 1, 0, 0,-1, 1}, // B
    { 0,-1, 1,-1, 1, 0, 0,-1, 1}, // B
    { 0,-1,-1, 1, 1,-1, 1, 0, 0}, // C
    { 0,-1,-1, 1, 1,-1, 1, 0, 0}, // C
    {-1, 0, 1,-1, 0, 1,-1, 0, 1,
     -1, 0, 1,-1,    1,-1, 0, 1,
     -1, 0, 1,-1, 0, 1,-1, 0, 1}};
  private static final int[][] K2S = {
    { 0,-1,-1, 1, 1,-1, 1, 0, 0}, // C
    { 0,-1,-1, 1, 1,-1, 1, 0, 0}, // C
    { 1, 1, 1, 1, 1, 1, 1, 1, 1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1}, // A
    { 0,-1, 1,-1, 1, 0, 0,-1, 1}, // B
    { 0,-1, 1,-1, 1, 0, 0,-1, 1}, // B
    {-1,-1,-1, 0, 0, 0, 1, 1, 1,
     -1,-1,-1, 0,    0, 1, 1, 1,
     -1,-1,-1, 0, 0, 0, 1, 1, 1}};
  private static final int[][] K3S = {
    { 0,-1, 1,-1, 1, 0, 0,-1, 1}, // B
    { 0,-1, 1,-1, 1, 0, 0,-1, 1}, // B
    { 0,-1,-1, 1, 1,-1, 1, 0, 0}, // C
    { 0,-1,-1, 1, 1,-1, 1, 0, 0}, // C
    { 1, 1, 1, 1, 1, 1, 1, 1, 1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1}, // A
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,
      0, 0, 0, 0,    0, 0, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1}};

  // A sample has indices and is either active or inactive.
  private static class Sample {
    int i1,i2,i3; // sample indices
    int ia; // determines whether this sample is active
    Sample(int i1, int i2, int i3) {
      this.i1 = i1;
      this.i2 = i2;
      this.i3 = i3;
    }
  }

  // Returns true if specified sample is active; false, otherwise.
  private boolean isActive(int i1, int i2, int i3) {
    return _s[i3][i2][i1].ia==_active;
  }

  // Marks all samples inactive. For efficiency, we typically do not loop 
  // over all the samples to clear their active flags. Usually we simply
  // increment the active value with which the flags are compared.
  private void clearActive() {
    if (_active==Integer.MAX_VALUE) { // rarely!
      _active = 1;
      for (int i3=0; i3<_n3; ++i3) {
        for (int i2=0; i2<_n2; ++i2) {
          for (int i1=0; i1<_n1; ++i1) {
            _s[i3][i2][i1].ia = 0;
          }
        }
      }
    } else { // typically
      ++_active;
    }
  }

  // Queue of active samples.
  private class ActiveQueue {
    synchronized Sample get() {
      Sample s = _q.remove();
      s.ia -= 1;
      --_n;
      return s;
    }
    synchronized void put(int i1, int i2, int i3) {
      Sample s = _s[i3][i2][i1];
      _q.add(s);
      s.ia = _active;
      ++_n;
    }
    boolean isEmpty() {
      return _n==0;
    }
    int size() {
      return _n;
    }
    int _n;
    private ArrayQueue<Sample> _q = new ArrayQueue<Sample>(1024);
  }

  /**
   * Zeros the time for the specified sample and recursively updates times 
   * for neighbor samples until all times have converged.
   */
  private void solveFrom(int i1, int i2, int i3) {

    // All samples initially inactive.
    clearActive();

    // Zero the time for the specified sample.
    _t[i3][i2][i1] = 0.0f;

    // Put six neighbor samples into the active queue.
    ActiveQueue q = new ActiveQueue();
    for (int k=0; k<6; ++k) {
      int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
      int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
      int j3 = i3+K3[k];  if (j3<0 || j3>=_n3) continue;
      q.put(j1,j2,j3);
    }

    // Complete the solve by processing the active queue until empty.
    if (_parallel) {
      solveParallel(q);
    } else {
      solveSerial(q);
    }
  }

  /**
   * Solves for times by sequentially processing each sample in the queue.
   */
  private void solveSerial(ActiveQueue q) {
    while (!q.isEmpty()) {
      solveOne(q);
    }
  }
  
  /**
   * Solves for times by processing samples in the queue in parallel.
   */
  private void solveParallel(final ActiveQueue aq) {
    int ntask = Runtime.getRuntime().availableProcessors();
    ExecutorService es = Executors.newFixedThreadPool(ntask);
    CompletionService<Void> cs = new ExecutorCompletionService<Void>(es);
    final AtomicInteger ai = new AtomicInteger();
    while (!aq.isEmpty()) {
      final int nq = aq.size();
      ai.set(0);
      for (int itask=0; itask<ntask; ++itask) {
        cs.submit(new Callable<Void>() {
          public Void call() {
            for (int iq=ai.getAndIncrement(); iq<nq; iq=ai.getAndIncrement())
              solveOne(aq);
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
    }
    es.shutdown();
  }

  /**
   * Processes one sample from the active queue.
   */
  private void solveOne(ActiveQueue aq) {

    // Get one sample from active queue.
    Sample i = aq.get();
    int i1 = i.i1;
    int i2 = i.i2;
    int i3 = i.i3;

    // Current time and new time.
    float ti = _t[i3][i2][i1];
    float gi = g(i1,i2,i3,K1S[6],K2S[6],K3S[6]);
    _t[i3][i2][i1] = gi;

    // If new and current times are close (converged), then ...
    if (ti-gi<ti*EPSILON) {

      // For all six neighbor samples, ...
      for (int k=0; k<6; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
        int j3 = i3+K3[k];  if (j3<0 || j3>=_n3) continue;

        // If neighbor is not in the active queue, ...
        if (!isActive(j1,j2,j3)) {

          // Compute time for the neighbor.
          float gj = g(j1,j2,j3,K1S[k],K2S[k],K3S[k]);

          // If computed time less than the neighbor's current time, ...
          if (gj<_t[j3][j2][j1]) {

            // Replace the current time.
            _t[j3][j2][j1] = gj;
            
            // Put the neighbor sample into the active queue.
            aq.put(j1,j2,j3);
          }
        }
      }
    }

    // Else, if not converged, put this sample back into the active queue.
    else {
      aq.put(i1,i2,i3);
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
   * Returns a valid time t computed via the Godunov-Hamiltonian.
   * Computations are limited to neighbor indices in arrays k1s, k2s, and k3s.
   */
  private float g(int i1, int i2, int i3, int[] k1s, int[] k2s, int[] k3s) {
    float tmin = _t[i3][i2][i1];

    // Get tensor coefficients.
    float[] d = new float[6];
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
        float ddet = d22*d33-d23*d23;  if (ddet==0.0f) continue;
        float s2 = (d23*d13-d12*d33)*s1/ddet;
        float s3 = (d23*d12-d13*d22)*s1/ddet;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tmin && t0>=t1) {
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          if (isValid1(i1,i2,i3,k1,p1,t0)) {
            tmin = t0;
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
        float ddet = d11*d33-d13*d13;  if (ddet==0.0f) continue;
        float s1 = (d13*d23-d12*d33)*s2/ddet;
        float s3 = (d13*d12-d23*d11)*s2/ddet;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tmin && t0>=t2) {
          float t01 = t0-t1;
          float t03 = t0-t3;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          if (isValid2(i1,i2,i3,k2,p2,t0)) {
            tmin = t0;
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
        float ddet = d11*d22-d12*d12;  if (ddet==0.0f) continue;
        float s1 = (d12*d23-d13*d22)*s3/ddet;
        float s2 = (d12*d13-d23*d11)*s3/ddet;
        float t0 = solveQuadratic(d11,d12,d13,d22,d23,d33,s1,s2,s3,t1,t2,t3);
        if (t0<tmin && t0>=t3) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid3(i1,i2,i3,k3,p3,t0)) {
            tmin = t0;
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
        if (t0<tmin && t0>=min(t2,t3)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid2(i1,i2,i3,k2,p2,t0) &&
              isValid3(i1,i2,i3,k3,p3,t0)) {
            tmin = t0;
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
        if (t0<tmin && t0>=min(t1,t3)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid1(i1,i2,i3,k1,p1,t0) &&
              isValid3(i1,i2,i3,k3,p3,t0)) {
            tmin = t0;
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
        if (t0<tmin && t0>=min(t1,t2)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          if (isValid1(i1,i2,i3,k1,p1,t0) &&
              isValid2(i1,i2,i3,k2,p2,t0)) {
            tmin = t0;
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
        if (t0<tmin && t0>=min(t1,t2,t3)) {
          float t01 = t0-t1;
          float t02 = t0-t2;
          float t03 = t0-t3;
          float p1 = (d12*s2*t02+d13*s3*t03)/d11;
          float p2 = (d12*s1*t01+d23*s3*t03)/d22;
          float p3 = (d13*s1*t01+d23*s2*t02)/d33;
          if (isValid1(i1,i2,i3,k1,p1,t0) &&
              isValid2(i1,i2,i3,k2,p2,t0) &&
              isValid3(i1,i2,i3,k3,p3,t0)) {
            tmin = t0;
          }
        }
      }
    }

    return tmin;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  /*
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
  */

  private static class ConstantTensors implements FimSolver3.Tensors {
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

  private static void testConstant() {
    int n1 = 51;
    int n2 = 51;
    int n3 = 51;
    float d11 = 1.000f, d12 = 0.000f, d13 = 0.000f,
                        d22 = 1.000f, d23 = 0.000f,
                                      d33 = 1.000f;
    ConstantTensors dt = new ConstantTensors(d11,d12,d13,d22,d23,d33);
    FimSolver3 fs = new FimSolver3(n1,n2,n3,dt);
    fs.setParallel(true);
    Stopwatch sw = new Stopwatch();
    sw.start();
    fs.zeroAt(2*n1/4,2*n2/4,2*n3/4);
    //fs.zeroAt(1*n1/4,1*n2/4,1*n3/4);
    //fs.zeroAt(3*n1/4,3*n2/4,3*n3/4);
    sw.stop();
    trace("time="+sw.time());
    float[][][] t = fs.getTimes();
    //Array.dump(t);
    //plot(t);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
  private static void trace(int i1, int i2, int i3, String s) {
    if (i1==2 && i2==2 && i3==2)
      trace(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testConstant();
      }
    });
  }
}
