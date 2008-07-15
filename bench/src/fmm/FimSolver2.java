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
   * @param dt diffusion tensors.
   */
  public FimSolver2(int n1, int n2, Tensors dt) {
    _n1 = n1;
    _n2 = n2;
    _dt = dt;
    _t = Array.fillfloat(INFINITY,n1,n2);
    _ia = new Index[_n1*_n2];
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param dt diffusion tensors.
   */
  public FimSolver2(float[][] t, Tensors dt) {
    _n1 = t[0].length;
    _n2 = t.length;
    _dt = dt;
    _t = t;
    _ia = new Index[_n1*_n2];
  }

  /**
   * Zeros the time at the specified sample and updates times elsewhere.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @return the updated array of times; by reference, not by copy.
   */
  public float[][] zeroAt(int i1, int i2) {
    _t[i2][i1] = 0.0f;
    for (int k=0; k<4; ++k)
      activate(i1+K1[k],i2+K2[k]);
    solve();
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

  // Diffusion tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, float[] d) {
      d[0] = 1.00f; // d11
      d[1] = 0.00f; // d12
      d[2] = 1.00f; // d22
    }
  }

  private static final float INFINITY = Float.MAX_VALUE;

  private static final int[] K1 = { 1, 0,-1, 0};
  private static final int[] K2 = { 0, 1, 0,-1};
  private static final int[][] K1S = {
    {-1,-1,-1},{-1, 0, 1},{ 1, 1, 1},{-1, 0, 1},{ 1, 0,-1, 0, 1,-1,-1, 1}};
  private static final int[][] K2S = {
    {-1, 0, 1},{-1,-1,-1},{-1, 0, 1},{ 1, 1, 1},{ 0, 1, 0,-1, 1, 1,-1,-1}};

  private static class Index {
    int i1,i2;
    Index(int i1, int i2) {
      this.i1 = i1;
      this.i2 = i2;
    }
    public boolean equals(Index o) {
      Index i = (Index)o;
      return i.i1==i1 && i.i2==i2;
    }
    public int hashCode() {
      return i1^i2;
    }
  }

  private int _n1,_n2;
  private Tensors _dt;
  private float[][] _t;
  private Index[] _ia;
  private HashMap<Index,Index> _as = new HashMap<Index,Index>(1024);
  private int[][] _mark;

  // Marks used during computation of times. For efficiency, we do not
  // loop over all the marks to clear them before beginning a fast 
  // marching loop. Instead, we simply modify the mark values.
  private int _far = 0; // samples with no time
  private int _trial = 1; // samples in active list
  private int _known = 2; // samples with a known time
  private void clearMarks() {
    if (_known+2>Integer.MAX_VALUE) { // if we must loop over all marks, ...
      _far = 0;
      _trial = 1;
      _known = 2;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _mark[i2][i1] = _far;
        }
      }
    } else {
      _far += 2; // all known samples instantly become far samples
      _trial +=2; // no samples are trial
      _known +=2; // no samples are known
    }
  }

  private void activate(int i1, int i2) {
    if (0<=i1 && i1<_n1 &&
        0<=i2 && i2<_n2) {
      Index i = new Index(i1,i2);
      _as.put(i,i);
    }
  }

  private void deactivate(Index i) {
    _as.remove(i);
  }

  /**
   * Jeong's fast iterative marching algorithm. While the active
   * set is not empty, this method computes times until converged.
   */
  private void solve() {

    // Index object used below to check for neighbors in active set.
    Index jj = new Index(-1,-1);

    // Tolerance for convergence.
    float epsilon = 0.01f; // tolerance for convergence

    // While the active set is not empty, ...
    for (int niter=0; !_as.isEmpty(); ++niter) {

      // DEBUG
      //if (niter%10==1)
      //  plot(_t);

      // Copy sample indices from the active set to a simple array.
      // The active set can now change as we loop over samples.
      int n = _as.size();
      _as.keySet().toArray(_ia);
      //trace("niter="+niter+" n="+n); // DEBUG

      // For all samples in the active set, ...
      for (int i=0; i<n; ++i) {

        // Sample indices.
        Index ii = _ia[i];
        int i1 = ii.i1;
        int i2 = ii.i2;

        // Current time and new time.
        float ti = _t[i2][i1];
        float gi = g(i1,i2,K1S[4],K2S[4]);
        _t[i2][i1] = gi;

        // If the new and old times are close (converged), then ...
        if (ti-gi<ti*epsilon) {

          // For all neighbor samples, ...
          for (int k=0; k<4; ++k) {

            // Neighbor sample indices; skip if out of bounds.
            int j1 = i1+K1[k];
            int j2 = i2+K2[k];
            if (j1<0 || j1>=_n1) continue;
            if (j2<0 || j2>=_n2) continue;

            // If the neighbor is not already in the active set, ...
            jj.i1 = j1;
            jj.i2 = j2;
            if (!_as.containsKey(jj)) {

              // Compute time for the neighbor.
              float gj = g(j1,j2,K1S[k],K2S[k]);

              // If computed time less than the neighbor's current time, ...
              if (gj<_t[j2][j1]) {

                // Replace the current time and activate the neighbor.
                _t[j2][j1] = gj;
                activate(j1,j2);
              }
            }
          }

          // Deactivate this sample (because it's time has converged).
          // It will be reactivated later if time for any of it's
          // neighbors changes.
          deactivate(ii);
        }
      }
    }
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

  private static int N_COMPUTE_TIME = 0;

  /**
   * Returns a valid time t computed via the Godunov-Hamiltonian.
   * Computations are limited to neighbor indices in arrays k1s and k2s.
   */
  private float g(int i1, int i2, int[] k1s, int[] k2s) {
    ++N_COMPUTE_TIME;
    float tmin = _t[i2][i1];

    // Get tensor coefficients.
    float[] d = new float[3];
    _dt.getTensor(i1,i2,d);
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
  // seem to be less robust than Jeong's; the active set tends to
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
    ++N_COMPUTE_TIME;
    float tmin = _t[i2][i1];

    // Get tensor coefficients.
    float[] d = new float[3];
    _dt.getTensor(i1,i2,d);
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
    pv.setColorModel(ColorMap.JET);
    //pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setInterpolation(PixelsView.Interpolation.LINEAR);
  }

  private static void testConstant() {
    int n1 = 501;
    int n2 = 501;
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
    int nct = N_COMPUTE_TIME;
    trace("time="+sw.time()+" nct="+nct);
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
