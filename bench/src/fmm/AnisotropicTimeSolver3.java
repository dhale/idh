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
 * This solver uses an iterative method to compute the solution times t. The 
 * iterations are similar to those described by Jeong and Whitaker (2007),
 * but the local method used here for computing each sampled time from the 
 * times of its neighbor samples is both simpler and faster. Note that the 
 * structure tensors S are inverses of the tensors used by Jeong and 
 * Whitaker.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.19
 */
public class AnisotropicTimeSolver3 {

  /**
   * Type of concurrency used when solving for times.
   */
  public enum Concurrency {
    PARALLEL,
    SERIAL
  };

  /**
   * An interface for classes of structure tensors. Each tensor is a
   * symmetric positive-definite 3-by-3 matrix 
   * {{s11,s12,s13},{s12,s22,s23},{s13,s23,s33}}.
   */
  public interface Tensors {

    /**
     * Gets structure tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param i3 index for 3rd dimension.
     * @param s array {s11,s12,s13,s22,s23,s33} of tensor elements.
     */
    public void getTensor(int i1, int i2, int i3, float[] s);
  }

  /**
   * Constructs a solver with constant identity structure tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   */
  public AnisotropicTimeSolver3(int n1, int n2, int n3) {
    this(n1,n2,n3,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified structure tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param tensors structure tensors.
   */
  public AnisotropicTimeSolver3(int n1, int n2, int n3, Tensors tensors) {
    init(n1,n2,n3,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors structure tensors.
   */
  public AnisotropicTimeSolver3(float[][][] t, Tensors tensors) {
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
  private static final float EPSILON = 0.0001f;

  private int _n1,_n2,_n3;
  private Tensors _tensors;
  private float[][][] _t;
  private Sample[][][] _s;
  private int _active = 0;
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

  // Structure tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, int i3, float[] s) {
      s[0] = 1.00f; // s11
      s[1] = 0.00f; // s12
      s[2] = 0.00f; // s13
      s[3] = 1.00f; // s22
      s[4] = 0.00f; // s23
      s[5] = 1.00f; // s33
    }
  }

  // Sample index offsets for the 6 neighbor samples.
  private static final int[] K1 = {-1, 1, 0, 0, 0, 0};
  private static final int[] K2 = { 0, 0,-1, 1, 0, 0};
  private static final int[] K3 = { 0, 0, 0, 0,-1, 1};

  // Indices of 8 neighbor tetrahedra. These indices limit the number of 
  // tetrahedra to consider when updating a neighbor sample. For example, 
  // when updating the neighbor with offsets {K1[3],K2[3],K3[3]} = {0,1,0}, 
  // we consider only tetrahedra with indices KT[3] = {2,3,6,7}. The last 
  // array contains the indices of all eight tetrahedra, and is used when 
  // computing the time at X0 from all six neighbor samples.
  private static final int[][] KT  = {
    {0,3,4,7},{1,2,5,6},{0,1,4,5},{2,3,6,7},{0,1,2,3},{4,5,6,7},
    {0,1,2,3,4,5,6,7}
  };

  // Sample index offsets for vertices X1 of 8 neighbor tetrahedra.
  private static final int[] K11 = { 1, 0,-1, 0, 1, 0,-1, 0};
  private static final int[] K12 = { 0, 1, 0,-1, 0, 1, 0,-1};
  private static final int[] K13 = { 0, 0, 0, 0, 0, 0, 0, 0};

  // Sample index offsets for vertices X2 of 8 neighbor tetrahedra.
  private static final int[] K21 = { 0,-1, 0, 1, 0,-1, 0, 1};
  private static final int[] K22 = { 1, 0,-1, 0, 1, 0,-1, 0};
  private static final int[] K23 = { 0, 0, 0, 0, 0, 0, 0, 0};

  // Sample index offsets for vertices X3 of 8 neighbor tetrahedra.
  private static final int[] K31 = { 0, 0, 0, 0, 0, 0, 0, 0};
  private static final int[] K32 = { 0, 0, 0, 0, 0, 0, 0, 0};
  private static final int[] K33 = { 1, 1, 1, 1,-1,-1,-1,-1};

  // A sample has indices and is either active or inactive.
  // For efficiency the active flag is an integer and not a boolean, 
  // so that we need not loop over all samples when initializing them 
  // to be inactive. See the comments for the method clearActive below.
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

    // Solve by processing the active queue until empty.
    if (_concurrency==Concurrency.PARALLEL) {
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
   * TODO: try using the fork-join framework in JSR-166.
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

    // Current time and new time computed from all six neighbors.
    float ti = _t[i3][i2][i1];
    float gi = g(i1,i2,i3,KT[6]);
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
          float gj = g(j1,j2,j3,KT[k]);

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

  // Returns the parameter a (alpha) that minimizes
  // t(a) = a*u1+sqrt((y2-a*y1)'*S*(y2-a*y1))
  // subject to the constraint 0 <= a <= 1.
  private static float computeAlpha(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float u1,
    float y11, float y12, float y13,
    float y21, float y22, float y23)
  {
    float z11 = s11*y11+s12*y12+s13*y13;
    float z12 = s12*y11+s22*y12+s23*y13;
    float z13 = s13*y11+s23*y12+s33*y13;
    float z21 = s11*y21+s12*y22+s13*y23;
    float z22 = s12*y21+s22*y22+s23*y23;
    float z23 = s13*y21+s23*y22+s33*y23;
    float d11 = y11*z11+y12*z12+y13*z13;
    float d12 = y11*z21+y12*z22+y13*z23;
    float d22 = y21*z21+y22*z22+y23*z23;
    return computeAlpha(u1,d11,d12,d22);
  }

  // Returns the parameter a (alpha) that minimizes
  // t(a) = a*u1+sqrt(d22-2.0f*a*d12+a*a*d11)
  // subject to the constraint 0 <= a <= 1.
  private static float computeAlpha(
    float u1, float d11, float d12, float d22) 
  {
    float alpha;
    float dd = d11*d22-d12*d12;
    if (dd<0.0f)
      dd = 0.0f;
    float du = d11-u1*u1;
    if (du<=0.0f) {
      alpha = (u1>=0.0f)?0.0f:1.0f;
    } else {
      alpha = (d12-u1*sqrt(dd/du))/d11;
      if (alpha<=0.0f) {
        alpha = 0.0f;
      } else if (alpha>=1.0f) {
        alpha = 1.0f;
      }
    }
    return alpha;
  }

  // Returns the time t+sqrt(y'*S*y).
  private static float computeTime(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float t, float y1, float y2, float y3)
  {
    float z1 = s11*y1+s12*y2+s13*y3;
    float z2 = s12*y1+s22*y2+s23*y3;
    float z3 = s13*y1+s23*y2+s33*y3;
    float dd = y1*z1+y2*z2+y3*z3;
    return t+sqrt(dd);
  }

  // Returns the minimum time 
  // t(a) = u2+a*u1+sqrt((y2-a*y1)'*S*(y2-a*y1))
  // subject to the constraint 0 <= a <= 1.
  private static float computeTime(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float u1, float u2,
    float y11, float y12, float y13,
    float y21, float y22, float y23)
  {
    float z11 = s11*y11+s12*y12+s13*y13;
    float z12 = s12*y11+s22*y12+s23*y13;
    float z13 = s13*y11+s23*y12+s33*y13;
    float z21 = s11*y21+s12*y22+s13*y23;
    float z22 = s12*y21+s22*y22+s23*y23;
    float z23 = s13*y21+s23*y22+s33*y23;
    float d11 = y11*z11+y12*z12+y13*z13;
    float d12 = y11*z21+y12*z22+y13*z23;
    float d22 = y21*z21+y22*z22+y23*z23;
    float alpha = computeAlpha(u1,d11,d12,d22);
    return u2+alpha*u1+sqrt(d22-2.0f*alpha*d12+alpha*alpha*d11);
  }

  // Returns the minimum time
  // t(a1,a2) = u3+a1*u1+a2*u2+sqrt((y3-a1*y1-a2*y2)'*S*(y3-a1*y1-a2*y2))
  // subject to the constraints 0 <= a1, 0 <= a2, and 0 <= a3 = 1-a1-a2.
  private static float computeTime(
    float s11, float s12, float s13, float s22, float s23, float s33,
    float u1, float u2, float u3,
    float y11, float y12, float y13,
    float y21, float y22, float y23,
    float y31, float y32, float y33)
  {
    // Inner products with respect to metric tensor S.
    float z11 = s11*y11+s12*y12+s13*y13;
    float z12 = s12*y11+s22*y12+s23*y13;
    float z13 = s13*y11+s23*y12+s33*y13;
    float z21 = s11*y21+s12*y22+s13*y23;
    float z22 = s12*y21+s22*y22+s23*y23;
    float z23 = s13*y21+s23*y22+s33*y23;
    float z31 = s11*y31+s12*y32+s13*y33;
    float z32 = s12*y31+s22*y32+s23*y33;
    float z33 = s13*y31+s23*y32+s33*y33;
    float d11 = y11*z11+y12*z12+y13*z13;
    float d12 = y11*z21+y12*z22+y13*z23;
    float d13 = y11*z31+y12*z32+y13*z33;
    float d22 = y21*z21+y22*z22+y23*z23;
    float d23 = y21*z31+y22*z32+y23*z33;
    float d33 = y31*z31+y32*z32+y33*z33;

    // The minimum lies on the line a1*alpha1 + a2*alpha2 + bb = 0.
    float a1 = u1*d12-u2*d11;
    float a2 = u1*d22-u2*d12;
    float alpha1,alpha2;

    // If a1 == a2 (== b == 0), solve easily for both alpha1 and alpha2.
    if (a1==0.0f && a2==0.0f) {
      float dd = d11*d22-d12*d12;
      alpha1 = (d13*d22-d12*d23)/dd;
      alpha2 = (d23*d11-d13*d12)/dd;
    }
    
    // Else, if a1 != a2, 
    else {
      float bb = u2*d13-u1*d23;
      float aa1 = (a1>=0.0f)?a1:-a1;
      float aa2 = (a2>=0.0f)?a2:-a2;

      // if abs(a1) <= abs(a2), solve for alpha1 first, then alpha2.
      if (aa1<=aa2) {
        float aoa = a1/a2;
        float boa = bb/a2;
        float v1 = u1-aoa*u2;
        float w21 = y31+boa*y21;
        float w22 = y32+boa*y22;
        float w23 = y33+boa*y23;
        float w11 = y11-aoa*y21;
        float w12 = y12-aoa*y22;
        float w13 = y13-aoa*y23;
        alpha1 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              v1,w11,w12,w13,w21,w22,w23);
        alpha2 = -boa-aoa*alpha1;
      }

      // Else if abs(a1) > abs(a2), solve for alpha2 first, then alpha1.
      else {
        float aoa = a2/a1;
        float boa = bb/a1;
        float v1 = u2-aoa*u1;
        float w21 = y31+boa*y11;
        float w22 = y32+boa*y12;
        float w23 = y33+boa*y13;
        float w11 = y21-aoa*y11;
        float w12 = y22-aoa*y12;
        float w13 = y23-aoa*y13;
        alpha2 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              v1,w11,w12,w13,w21,w22,w23);
        alpha1 = -boa-aoa*alpha2;
      } 
    }

    // The sum alpha1 + alpha2 + alpha3 = 1.
    float alpha3 = 1.0f-alpha1-alpha2;

    // Initial time is huge.
    float t = INFINITY;

    // If minimum is strictly inside the triangle 123, ...
    if (alpha1>0.0f && alpha2>0.0f && alpha3>0.0f) {
      t = u3+alpha1*u1+alpha2*u2 + 
        sqrt(d33+alpha1*alpha1*d11+alpha2*alpha2*d22 +
             2.0f*alpha1*alpha2*d12-2.0f*alpha1*d13-2.0f*alpha2*d23);

    // Else if the minimum is not strictly inside the triangle 123,
    // search for the minimum on the edges of the triangle.
    } else {

      // If minimum could be on the edge 23 (alpha1=0), ...
      if (alpha1<=0.0f) {
        alpha2 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              u2,y21,y22,y23,y31,y32,y33);
        float t23 = u3+alpha2*u2+sqrt(d33+alpha2*alpha2*d22-2.0f*alpha2*d23);
        if (t23<t) t = t23;
      }

      // If minimum could be on the edge 13 (alpha2=0), ...
      if (alpha2<=0.0f) {
        alpha1 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              u1,y11,y12,y13,y31,y32,y33);
        float t13 = u3+alpha1*u1+sqrt(d33+alpha1*alpha1*d11-2.0f*alpha1*d13);
        if (t13<t) t = t13;
      }

      // If minimum could be on the edge 12 (alpha3=0), ...
      if (alpha3<=0.0f) {
        alpha1 = computeAlpha(s11,s12,s13,s22,s23,s33,
                              u1-u2,
                              y11-y21,y12-y22,y13-y23,
                              y31-y21,y32-y22,y33-y23);
        alpha2 = 1.0f-alpha1;
        float t12 = u3+alpha1*u1+alpha2*u2 + 
          sqrt(d33+alpha1*alpha1*d11+alpha2*alpha2*d22 +
               2.0f*alpha1*alpha2*d12-2.0f*alpha1*d13-2.0f*alpha2*d23);
        if (t12<t) t = t12;
      }
    }
    return t;
  }

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor tets with specified indices.
   */
  private float g(int i1, int i2, int i3, int[] kt) {
    float tmin = _t[i3][i2][i1];

    // Get tensor coefficients.
    float[] s = new float[6];
    _tensors.getTensor(i1,i2,i3,s);
    float s11 = s[0];
    float s12 = s[1];
    float s13 = s[2];
    float s22 = s[3];
    float s23 = s[4];
    float s33 = s[5];

    // For all relevant neighbor tets, ...
    for (int it=0; it<kt.length; ++it) {
      int jt = kt[it];

      // Sample indices of vertices X0, X1, X2 and X3 of neighbor tet.
      int j01 = i1;
      int j02 = i2;
      int j03 = i3;
      int j11 = i1+K11[jt];  if (j11<0 || j11>=_n1) continue;
      int j12 = i2+K12[jt];  if (j12<0 || j12>=_n2) continue;
      int j13 = i3+K13[jt];  if (j13<0 || j13>=_n3) continue;
      int j21 = i1+K21[jt];  if (j21<0 || j21>=_n1) continue;
      int j22 = i2+K22[jt];  if (j22<0 || j22>=_n2) continue;
      int j23 = i3+K23[jt];  if (j23<0 || j23>=_n3) continue;
      int j31 = i1+K31[jt];  if (j31<0 || j31>=_n1) continue;
      int j32 = i2+K32[jt];  if (j32<0 || j32>=_n2) continue;
      int j33 = i3+K33[jt];  if (j33<0 || j33>=_n3) continue;

      // Times T0, T1, T2 and T3 for vertices X0, X1, X2 and X3 of tet.
      float t0; // to be computed below
      float t1 = _t[j13][j12][j11];
      float t2 = _t[j23][j22][j21];
      float t3 = _t[j33][j32][j31];

      // Use times < INFINITY in {T1,T2,T3} to compute candidate time T0.
      if (t1<INFINITY) {
        if (t2<INFINITY) {
          if (t3<INFINITY) { // use {T1,T2,T3}
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1-t3,t2-t3,t3,
                             j11-j31,j12-j32,j13-j33,
                             j21-j31,j22-j32,j23-j33,
                             j01-j31,j02-j32,j03-j33);
          } else {  // use {T1,T2}
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1-t2,t2,
                             j11-j21,j12-j22,j13-j23,
                             j01-j21,j02-j22,j03-j23);
          }
        } else {
          if (t3<INFINITY) { // use {T1,T3}
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1-t3,t3,
                             j11-j31,j12-j32,j13-j33,
                             j01-j31,j02-j32,j03-j33);
          } else { // use {T1}
            t0 = computeTime(s11,s12,s13,s22,s23,s33,
                             t1,j01-j11,j02-j12,j03-j13);
          }
        }
      } else if (t2<INFINITY) {
        if (t3<INFINITY) { // use {T2,T3}
          t0 = computeTime(s11,s12,s13,s22,s23,s33,
                           t2-t3,t3,
                           j21-j31,j22-j32,j23-j33,
                           j01-j31,j02-j32,j03-j33);
        } else { // use {T2}
          t0 = computeTime(s11,s12,s13,s22,s23,s33,
                           t2,j01-j21,j02-j22,j03-j23);
        }
      } else { // use {T3}
        t0 = computeTime(s11,s12,s13,s22,s23,s33,
                         t3,j01-j31,j02-j32,j03-j33);
      }

      // If candidate T0 smaller than the min time, update the min time.
      if (t0<tmin)
        tmin = t0;
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

  private static class ConstantTensors 
    implements AnisotropicTimeSolver3.Tensors 
  {
    ConstantTensors(
      float s11, float s12, float s13, float s22, float s23, float s33) 
    {
      _s11 = s11;
      _s12 = s12;
      _s13 = s13;
      _s22 = s22;
      _s23 = s23;
      _s33 = s33;
    }
    public void getTensor(int i1, int i2, int i3, float[] s) {
      s[0] = _s11;
      s[1] = _s12;
      s[2] = _s13;
      s[3] = _s22;
      s[4] = _s23;
      s[5] = _s33;
    }
    private float _s11,_s12,_s13,_s22,_s23,_s33;
  }

  private static void testConstant() {
    int n1 = 101;
    int n2 = 101;
    int n3 = 101;
    float s11 = 1.000f, s12 = 0.000f, s13 = 0.000f,
                        s22 = 1.000f, s23 = 0.000f,
                                      s33 = 1.000f;
    ConstantTensors st = new ConstantTensors(s11,s12,s13,s22,s23,s33);
    AnisotropicTimeSolver3 ats = new AnisotropicTimeSolver3(n1,n2,n3,st);
    ats.setConcurrency(AnisotropicTimeSolver3.Concurrency.PARALLEL);
    Stopwatch sw = new Stopwatch();
    sw.start();
    //ats.zeroAt(0*n1/4,0*n2/4,0*n3/4);
    ats.zeroAt(2*n1/4,2*n2/4,2*n3/4);
    //ats.zeroAt(1*n1/4,1*n2/4,1*n3/4);
    //ats.zeroAt(3*n1/4,3*n2/4,3*n3/4);
    sw.stop();
    float[][][] t = ats.getTimes();
    trace("time="+sw.time()+" sum="+Array.sum(t));
    //Array.dump(t);
    //plot(t);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
  private static void trace(int i1, int i2, int i3, String s) {
    if (i1==0 && i2==0 && i3==0)
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
