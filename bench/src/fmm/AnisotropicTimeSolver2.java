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
 * This solver uses an iterative method to compute the solution times t. The 
 * iterations are similar to those described by Jeong and Whitaker (2007),
 * but the local method used here for computing each sampled time from the 
 * times of its neighbor samples is both simpler and faster. Note that the 
 * structure tensors S are inverses of the tensors used by Jeong and 
 * Whitaker.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.19
 */
public class AnisotropicTimeSolver2 {

  /**
   * Type of concurrency used when solving for times.
   */
  public enum Concurrency {
    PARALLEL,
    SERIAL
  };

  /**
   * An interface for classes of structure tensors. Each tensor is a
   * symmetric positive-definite 2-by-2 matrix {{s11,s12},{s12,s22}}.
   */
  public interface Tensors {

    /**
     * Gets structure tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param s array {s11,s12,s22} of tensor elements.
     */
    public void getTensor(int i1, int i2, float[] s);
  }

  /**
   * Constructs a solver with constant identity structure tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public AnisotropicTimeSolver2(int n1, int n2) {
    this(n1,n2,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified structure tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors structure tensors.
   */
  public AnisotropicTimeSolver2(int n1, int n2, Tensors tensors) {
    init(n1,n2,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors structure tensors.
   */
  public AnisotropicTimeSolver2(float[][] t, Tensors tensors) {
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
  private static final float EPSILON = 0.01f;

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

  // Structure tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, float[] s) {
      s[0] = 1.00f; // s11
      s[1] = 0.00f; // s12
      s[2] = 1.00f; // s22
    }
  }

  // Times for each sample are computed using four neighbor samples that
  // form four triangles. The triangles are indexed as follows:
  //       2 ^
  //         *
  //       / | \
  //     / 1 | 0 \
  //   * - - X - - * >
  //     \ 2 | 3 /   1
  //       \ | /
  //         *
  // The symbol X represents the vertex X0 shared by all neighbor triangles,
  // and the symbols * represent the other two triangle vertices X1 and X2.

  // Sample index offsets for the 4 neighbor samples.
  private static final int[] K1 = {-1, 1, 0, 0};
  private static final int[] K2 = { 0, 0,-1, 1};

  // Indices of 4 neighbor triangles. These indices limit the number of 
  // triangles that must be considered when updating a neighbor sample. 
  // For example, when updating the neighbor with offsets {K1[3],K2[3]} 
  // = {0,1}, we use only triangles with indices KT[3] = {2,3}. The last 
  // array contains the indices of all four triangles, and is used when 
  // computing the time at X0 from all four neighbor samples.
  private static final int[][] KT  = {
    {3,0},{1,2},{0,1},{2,3},
    {0,1,2,3}
  };

  // Sample index offsets for vertices X1 of 4 neighbor triangles.
  private static final int[] K11 = { 1, 0,-1, 0};
  private static final int[] K12 = { 0, 1, 0,-1};

  // Sample index offsets for vertices X2 of 4 neighbor triangles.
  private static final int[] K21 = { 0,-1, 0, 1};
  private static final int[] K22 = { 1, 0,-1, 0};

  // A sample has indices and is either active or inactive.
  // For efficiency the active flag is an integer and not a boolean, 
  // so that we need not loop over all samples when initializing them 
  // to be inactive. See the comments for the method clearActive below.
  private static class Sample {
    int i1,i2; // sample indices
    int ia; // determines whether this sample is active
    Sample(int i1, int i2) {
      this.i1 = i1;
      this.i2 = i2;
    }
  }

  // Returns true if specified sample is active; false, otherwise.
  private boolean isActive(int i1, int i2) {
    return _s[i2][i1].ia==_active;
  }

  // Marks all samples inactive. For efficiency, we typically do not loop 
  // over all the samples to clear their active flags. Usually we simply
  // increment the active value with which the flags are compared.
  private void clearActive() {
    if (_active==Integer.MAX_VALUE) { // rarely!
      _active = 1;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1].ia = 0;
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
    synchronized void put(int i1, int i2) {
      Sample s = _s[i2][i1];
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
  private void solveFrom(int i1, int i2) {

    // All samples initially inactive.
    clearActive();

    // Zero the time for the specified sample.
    _t[i2][i1] = 0.0f;

    // Put four neighbor samples into the active queue.
    ActiveQueue q = new ActiveQueue();
    for (int k=0; k<4; ++k) {
      int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
      int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
      q.put(j1,j2);
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
    //ntask = 1; //  3.7 s
    //ntask = 2; //  5.1 s
    //ntask = 3; // 10.6 s
    //ntask = 4; // 11.6 s
    ExecutorService es = Executors.newFixedThreadPool(ntask);
    CompletionService<Void> cs = new ExecutorCompletionService<Void>(es);
    final AtomicInteger ai = new AtomicInteger();
    int nqtotal = 0;
    while (!aq.isEmpty()) {
      final int nq = aq.size();
      nqtotal += nq;
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
    trace("solveParallel: nqtotal="+nqtotal);
  }

  /**
   * Processes one sample from the active queue.
   */
  private void solveOne(ActiveQueue aq) {

    // Get one sample from active queue.
    Sample i = aq.get();
    int i1 = i.i1;
    int i2 = i.i2;

    // Current time and new time computed from all four neighbors.
    float ti = _t[i2][i1];
    float gi = g(i1,i2,KT[4]);
    _t[i2][i1] = gi;

    // If new and current times are close (converged), then ...
    if (ti-gi<ti*EPSILON) {

      // For all four neighbor samples, ...
      for (int k=0; k<4; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;

        // If neighbor is not in the active queue, ...
        if (!isActive(j1,j2)) {

          // Compute time for the neighbor.
          float gj = g(j1,j2,KT[k]);

          // If computed time less than the neighbor's current time, ...
          if (gj<_t[j2][j1]) {

            // Replace the current time.
            _t[j2][j1] = gj;
            
            // Put the neighbor sample into the active queue.
            aq.put(j1,j2);
          }
        }
      }
    }

    // Else, if not converged, put this sample back into the active queue.
    else {
      aq.put(i1,i2);
    }
  }

  // Returns the time t+sqrt(y'*S*y).
  private static float computeTime(
    float s11, float s12, float s22,
    float t, float y1, float y2)
  {
    float z1 = s11*y1+s12*y2;
    float z2 = s12*y1+s22*y2;
    float dd = y1*z1+y2*z2;
    return t+sqrt(dd);
  }

  // Returns the minimum time 
  // t(a) = u2+a*u1+sqrt((y2-a*y1)'*S*(y2-a*y1))
  // subject to the constraint 0 <= a <= 1.
  private static float computeTime(
    float s11, float s12, float s22, float u1, float u2,
    float y11, float y12, float y21, float y22)
  {
    float z11 = s11*y11+s12*y12;
    float z12 = s12*y11+s22*y12;
    float z21 = s11*y21+s12*y22;
    float z22 = s12*y21+s22*y22;
    float d11 = y11*z11+y12*z12;
    float d12 = y11*z21+y12*z22;
    float d22 = y21*z21+y22*z22;
    float alpha = computeAlpha(u1,d11,d12,d22);
    return u2+alpha*u1+sqrt(d22-2.0f*alpha*d12+alpha*alpha*d11);
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

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor tris with specified indices.
   */
  private float g(int i1, int i2, int[] kt) {
    float tmin = _t[i2][i1];

    // Get tensor coefficients.
    float[] s = new float[3];
    _tensors.getTensor(i1,i2,s);
    float s11 = s[0];
    float s12 = s[1];
    float s22 = s[2];

    // For all relevant neighbor tris, ...
    for (int it=0; it<kt.length; ++it) {
      int jt = kt[it];

      // Sample indices of vertices X0, X1, and X2 of neighbor tri.
      int j01 = i1;
      int j02 = i2;
      int j11 = i1+K11[jt];  if (j11<0 || j11>=_n1) continue;
      int j12 = i2+K12[jt];  if (j12<0 || j12>=_n2) continue;
      int j21 = i1+K21[jt];  if (j21<0 || j21>=_n1) continue;
      int j22 = i2+K22[jt];  if (j22<0 || j22>=_n2) continue;

      // Times T0, T1, and T2 for vertices X0, X1, and X2 of tri.
      float t0; // to be computed below
      float t1 = _t[j12][j11];
      float t2 = _t[j22][j21];

      // Use times < INFINITY in {T1,T2} to compute candidate time T0.
      if (t1<INFINITY) {
        if (t2<INFINITY) {
          t0 = computeTime(s11,s12,s22,t1-t2,t2,
                           j11-j21,j12-j22,
                           j01-j21,j02-j22);
        } else {
          t0 = computeTime(s11,s12,s22,t1,j01-j11,j02-j12);
        }
      } else {
        t0 = computeTime(s11,s12,s22,t2,j01-j21,j02-j22);
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

  private static class ConstantTensors 
    implements AnisotropicTimeSolver2.Tensors 
  {
    ConstantTensors(float s11, float s12, float s22) {
      _s11 = s11;
      _s12 = s12;
      _s22 = s22;
    }
    public void getTensor(int i1, int i2, float[] s) {
      s[0] = _s11;
      s[1] = _s12;
      s[2] = _s22;
    }
    private float _s11,_s12,_s22;
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
    //float su = 1.000f;
    float su = 0.010f;
    float sv = 1.000f;
    float cosa = cos(angle);
    float sina = sin(angle);
    float s11 = su*cosa*cosa+sv*sina*sina;
    float s12 = (su-sv)*sina*cosa;
    float s22 = sv*cosa*cosa+su*sina*sina;
    ConstantTensors st = new ConstantTensors(s11,s12,s22);
    AnisotropicTimeSolver2 ats = new AnisotropicTimeSolver2(n1,n2,st);
    ats.setConcurrency(AnisotropicTimeSolver2.Concurrency.PARALLEL);
    //ats.setConcurrency(AnisotropicTimeSolver2.Concurrency.SERIAL);
    Stopwatch sw = new Stopwatch();
    sw.start();
    ats.zeroAt(2*n1/4,2*n2/4);
    //ats.zeroAt(1*n1/4,1*n2/4);
    //ats.zeroAt(3*n1/4,3*n2/4);
    sw.stop();
    trace("time="+sw.time());
    float[][] t = ats.getTimes();
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
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testConstant();
      }
    });
  }
}
