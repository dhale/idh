/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

// for testing
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;

/**
 * 2D implementation of Jeong and Whitakers' fast iterative method.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.14
 */
public class FimSolver2 {

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
  public FimSolver2(int n1, int n2) {
    this(n1,n2,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified diffusion tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors diffusion tensors.
   */
  public FimSolver2(int n1, int n2, Tensors tensors) {
    init(n1,n2,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors diffusion tensors.
   */
  public FimSolver2(float[][] t, Tensors tensors) {
    init(t[0].length,t.length,t,tensors);
  }

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

  /**
   * Zeros the time at the specified sample and updates times elsewhere.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @return the updated array of times; by reference, not by copy.
   */
  public float[][] zeroAt(int i1, int i2) {
    updateFrom(i1,i2);
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

  private int _n1,_n2;
  private Tensors _tensors;
  private float[][] _t;
  private Sample[][] _s;
  private int _active;

  // Diffusion tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, float[] d) {
      d[0] = 1.00f; // d11
      d[1] = 0.00f; // d12
      d[2] = 1.00f; // d22
    }
  }

  // Default time for samples not yet computed.
  private static final float INFINITY = Float.MAX_VALUE;

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.01f;

  // Nominal number of threads to use in this solver.
  private static int NTHREAD = Runtime.getRuntime().availableProcessors();

  private static final int[] K1 = { 1, 0,-1, 0};
  private static final int[] K2 = { 0, 1, 0,-1};
  private static final int[][] K1S = {
    {-1,-1,-1},{-1, 0, 1},{ 1, 1, 1},{-1, 0, 1},{ 1, 0,-1, 0, 1,-1,-1, 1}};
  private static final int[][] K2S = {
    {-1, 0, 1},{-1,-1,-1},{-1, 0, 1},{ 1, 1, 1},{ 0, 1, 0,-1, 1, 1,-1,-1}};

  // A sample has indices and is either active or inactive.
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
      Sample s = _q.poll();
      s.ia -= 1;
      --_n;
      return s;
    }
    synchronized void put(int i1, int i2) {
      Sample s = _s[i2][i1];
      _q.offer(s);
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
    private LinkedList<Sample> _q = new LinkedList<Sample>();
  }
  /*
    if (_naborEs==null) {
      _naborR = new NaborRunnable[26];
      for (int k=0; k<26; ++k)
        _naborR[k] = new NaborRunnable();
      int nthread = Runtime.getRuntime().availableProcessors();
      nthread *= 4;
      _naborEs = Executors.newFixedThreadPool(nthread);
      _naborCs = new ExecutorCompletionService<Void>(_naborEs);
    }
  */

  private void update(ActiveQueue q) {

    // While the active queue is not empty, ...
    while (!q.isEmpty()) {

      // For all samples *currently* in the active queue, ...
      int nq = q.size();
      for (int iq=0; iq<nq; ++iq) {

        // Get sample from the active queue.
        Sample i = q.get();
        int i1 = i.i1;
        int i2 = i.i2;

        // Update for this sample.
        update(i1,i2,q);
      }
    }
  }

  private void update(int i1, int i2, ActiveQueue q) {

    // Current time and new time.
    float ti = _t[i2][i1];
    float gi = g(i1,i2,K1S[4],K2S[4]);
    _t[i2][i1] = gi;

    // If new and current times are close (converged), then ...
    if (ti-gi<ti*EPSILON) {

      // For all four neighbor samples, ...
      for (int k=0; k<4; ++k) {

        // Neighbor sample indices; skip if out of bounds.
        int j1 = i1+K1[k];
        int j2 = i2+K2[k];
        if (j1<0 || j1>=_n1) continue;
        if (j2<0 || j2>=_n2) continue;

        // If neighbor is not in the active queue, ...
        if (!isActive(j1,j2)) {

          // Compute time for the neighbor.
          float gj = g(j1,j2,K1S[k],K2S[k]);

          // If computed time less than the neighbor's current time, ...
          if (gj<_t[j2][j1]) {

            // Replace the current time.
            _t[j2][j1] = gj;
            
            // Put the neighbor sample into the active queue.
            q.put(j1,j2);
          }
        }
      }
    }

    // Else, if not converged, put this sample back into the active queue.
    else {
      q.put(i1,i2);
    }
  }

  private void updateFrom(int i1, int i2) {

    // All samples initially inactive.
    clearActive();

    // Zero the time for the specified sample.
    _t[i2][i1] = 0.0f;

    // Put four neighbor samples into the active queue.
    ActiveQueue q = new ActiveQueue();
    for (int k=0; k<4; ++k) {
      int j1 = i1+K1[k];
      int j2 = i2+K2[k];
      if (0<=j1 && j1<_n1 && 0<=j2 && j2<_n2)
        q.put(j1,j2);
    }

    // Complete the update by processing the active queue.
    update(q);
  }

  /**
   * Solves the quadratic equation
   *   d11*s1*s1*(t1-t0)*(t1-t0) + 
   * 2*d12*s1*s2*(t1-t0)*(t2-t0) + 
   *   d22*s2*s2*(t2-t0)*(t2-t0) = 1
   * for a positive time t0. If no solution exists, because the 
   * discriminant is negative, this method returns INFINITY.
   */
  private static float solveQuadratic(
    float d11, float d12, float d22,
    float s1, float s2, float t1, float t2) 
  {
    double ds11 = d11*s1*s1;
    double ds12 = d12*s1*s2;
    double ds22 = d22*s2*s2;
    double t12 = t1-t2; // reduce rounding errors by solving for u = t0-t1
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
   * Returns a valid time t computed via the Godunov-Hamiltonian.
   * Computations are limited to neighbor indices in arrays k1s and k2s.
   */
  private float g(int i1, int i2, int[] k1s, int[] k2s) {
    float tmin = _t[i2][i1];

    // Get tensor coefficients.
    float[] d = new float[3];
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
        int j2 = i2+k2;
        if (j2<0 || j2>=_n2) continue;
        float t2 = _t[j2][i1];
        if (t2!=INFINITY) {
          float t1 = t2;
          float s2 = k2;
          float s1 = -s2*d12/d11;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tmin && t0>=t2) {
            float p2 = -s1*(t1-t0)*d12/d22;
            if (isValid2(i1,i2,k2,p2,t0)) {
              tmin = t0;
            }
          }
        }
      } 
      
      // else, if (p1-,p2s) or (p1+,p2s), ...
      else if (k2==0) {
        int j1 = i1+k1;
        if (j1<0 || j1>=_n1) continue;
        float t1 = _t[i2][j1];
        if (t1!=INFINITY) {
          float t2 = t1;
          float s1 = k1;
          float s2 = -s1*d12/d22;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tmin && t0>=t1) {
            float p1 = -s2*(t2-t0)*d12/d11;
            if (isValid1(i1,i2,k1,p1,t0)) {
              tmin = t0;
            }
          }
        }
      } 
      
      // else, if (p1-,p2-), (p1+,p2-), (p1-,p2+) or (p1+,p2+), ...
      else {
        int j2 = i2+k2;
        int j1 = i1+k1;
        if (j1<0 || j1>=_n1) continue;
        if (j2<0 || j2>=_n2) continue;
        float t1 = _t[i2][j1];
        float t2 = _t[j2][i1];
        if (t1!=INFINITY && t2!=INFINITY) {
          float s1 = k1;
          float s2 = k2;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tmin && t0>=min(t1,t2)) {
            float p1 = -s2*(t2-t0)*d12/d11;
            float p2 = -s1*(t1-t0)*d12/d22;
            if (isValid1(i1,i2,k1,p1,t0) &&
                isValid2(i1,i2,k2,p2,t0)) {
              tmin = t0;
            }
          }
        }
      }
    }

    return tmin;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // unused

  // Tsai's tests for valid solutions. For high anisotropy, these
  // seem to be less robust than Jeong's; the trial queue tends to
  // become larger in some simple tests, so also more costly.
  private static final float TSAI_THRESHOLD = 0.0001f;
  private boolean isValid1(
    int i1, int i2, float d11, float d12, float d22, 
    float s1, float t1, float t0) 
  {
    float p1 = s1*(t1-t0);
    float t2m = (i2>0    )?_t[i2-1][i1]:INFINITY;
    float t2p = (i2<_n2-1)?_t[i2+1][i1]:INFINITY;
    float p2s = -p1*d12/d22;
    float p2m = t0-t2m;
    float p2p = t2p-t0;
    float p2 = p2m; // p2 = sgn max ( (p2m-p2s)+ , (p2p-p2s)- ) + p2s
    if (p2m<p2s && p2s<p2p) {
      p2 = p2s;
    } else if (0.5f*(p2m+p2p)<p2s) {
      p2 = p2p;
    }
    float h = d11*p1*p1+2.0f*d12*p1*p2+d22*p2*p2-1.0f;
    return abs(h)<TSAI_THRESHOLD;
  }
  private boolean isValid2(
    int i1, int i2, float d11, float d12, float d22, 
    float s2, float t2, float t0) 
  {
    float p2 = s2*(t2-t0);
    float t1m = (i1>0    )?_t[i2][i1-1]:INFINITY;
    float t1p = (i1<_n1-1)?_t[i2][i1+1]:INFINITY;
    float p1s = -p2*d12/d11;
    float p1m = t0-t1m;
    float p1p = t1p-t0;
    float p1 = p1m; // p1 = sgn max ( (p1m-p1s)+ , (p1p-p1s)- ) + p1s
    if (p1m<p1s && p1s<p1p) {
      p1 = p1s;
    } else if (0.5f*(p1m+p1p)<p1s) {
      p1 = p1p;
    }
    float h = d11*p1*p1+2.0f*d12*p1*p2+d22*p2*p2-1.0f;
    return abs(h)<TSAI_THRESHOLD;
  }

  /**
   * Returns a valid time t computed via the Godunov-Hamiltonian.
   * This version does not limit computations to only relevant neighbors.
   */
  private float g(int i1, int i2) {
    float tmin = _t[i2][i1];

    // Get tensor coefficients.
    float[] d = new float[3];
    _tensors.getTensor(i1,i2,d);
    float d11 = d[0];
    float d12 = d[1];
    float d22 = d[2];

    // For (p1-,p2-), (p1+,p2-), (p1-,p2+), (p1+,p2+)
    for (int k2=-1; k2<=1; k2+=2) {
      int j2 = i2+k2;
      if (j2<0 || j2>=_n2) continue;
      float t2 = _t[j2][i1];
      for (int k1=-1; k1<=1; k1+=2) {
        int j1 = i1+k1;
        if (j1<0 || j1>=_n1) continue;
        float t1 = _t[i2][j1];
        if (t1!=INFINITY && t2!=INFINITY) {
          float s1 = k1;
          float s2 = k2;
          float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
          if (t0<tmin && t0>=min(t1,t2)) {
            float p1 = -s2*(t2-t0)*d12/d11;
            float p2 = -s1*(t1-t0)*d12/d22;
            if (isValid1(i1,i2,k1,p1,t0) &&
                isValid2(i1,i2,k2,p2,t0)) {
              tmin = t0;
            }
          }
        }
      }
    }

    // For (p1s,p2-), (p1s,p2+)
    for (int k2=-1; k2<=1; k2+=2) {
      int j2 = i2+k2;
      if (j2<0 || j2>=_n2) continue;
      float t2 = _t[j2][i1];
      if (t2!=INFINITY) {
        float t1 = t2;
        float s2 = k2;
        float s1 = -s2*d12/d11;
        float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
        if (t0<tmin && t0>=t2) {
          float p2 = -s1*(t1-t0)*d12/d22;
          if (isValid2(i1,i2,k2,p2,t0)) {
            tmin = t0;
          }
        }
      }
    }
  
    // For (p1-,p2s), (p1+,p2s)
    for (int k1=-1; k1<=1; k1+=2) {
      int j1 = i1+k1;
      if (j1<0 || j1>=_n1) continue;
      float t1 = _t[i2][j1];
      if (t1!=INFINITY) {
        float t2 = t1;
        float s1 = k1;
        float s2 = -s1*d12/d22;
        float t0 = solveQuadratic(d11,d12,d22,s1,s2,t1,t2);
        if (t0<tmin && t0>=t1) {
          float p1 = -s2*(t2-t0)*d12/d11;
          if (isValid1(i1,i2,k1,p1,t0)) {
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

  private static class ConstantTensors implements FimSolver2.Tensors {
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
    int n1 = 1001;
    int n2 = 1001;
    float angle = FLT_PI*110.0f/180.0f;
    float su = 1.000f;
    float sv = 0.010f;
    float cosa = cos(angle);
    float sina = sin(angle);
    float d11 = su*cosa*cosa+sv*sina*sina;
    float d12 = (su-sv)*sina*cosa;
    float d22 = sv*cosa*cosa+su*sina*sina;
    trace("d11="+d11+" d12="+d12+" d22="+d22+" d="+(d11*d22-d12*d12));
    ConstantTensors dt = new ConstantTensors(d11,d12,d22);
    FimSolver2 fs = new FimSolver2(n1,n2,dt);
    Stopwatch sw = new Stopwatch();
    sw.start();
    fs.zeroAt(2*n1/4,2*n2/4);
    //fs.zeroAt(1*n1/4,1*n2/4);
    //fs.zeroAt(3*n1/4,3*n2/4);
    sw.stop();
    trace("time="+sw.time());
    float[][] t = fs.getTimes();
    //Array.dump(t);
    plot(t);
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
