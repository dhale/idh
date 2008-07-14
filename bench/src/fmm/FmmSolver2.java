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
 * 2D implementation of Sethian's fast marching method.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.14
 */
public class FmmSolver2 {

  /**
   * An interface for classes of structure tensors. Each tensor is a
   * symmetric positive-definite 2-by-2 matrix {{d11,d12},{d12,d22}}.
   */
  public interface Tensors {

    /**
     * Gets structure tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param d array {d11,d12,d22} of tensor elements.
     */
    public void getTensor(int i1, int i2, float[] d);
  }

  /**
   * Constructs a solver with constant identity structure tensors.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public FmmSolver2(int n1, int n2) {
    this(n1,n2,new IdentityTensors());
  }
  
  /**
   * Constructs a solver for the specified structure tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param st structure tensors.
   */
  public FmmSolver2(int n1, int n2, Tensors st) {
    _n1 = n1;
    _n2 = n2;
    _st = st;
    _t = Array.fillfloat(INFINITY,n1,n2);
    _heap = new TimeHeap2(TimeHeap2.Type.MIN,_n1,_n2);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param st structure tensors.
   */
  public FmmSolver2(float[][] t, Tensors st) {
    _n1 = t[0].length;
    _n2 = t.length;
    _st = st;
    _t = t;
    _heap = new TimeHeap2(TimeHeap2.Type.MIN,_n1,_n2);
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

  // Sample index offsets for neighbor samples.
  private static final int[] K1 = {-1, 1, 0, 0};
  private static final int[] K2 = { 0, 0,-1, 1};

  // Sample index offsets for vertices X1 of the four neighbor tris.
  private static final int[] K11 = { 1, 0,-1, 0};
  private static final int[] K12 = { 0, 1, 0,-1};

  // Sample index offsets for vertices X2 of the four neighbor tris.
  private static final int[] K21 = { 0,-1, 0, 1};
  private static final int[] K22 = { 1, 0,-1, 0};

  // Indices of neighbor triangles.
  private static final int[][] KT  = {{1,2},{2,3},{0,3},{0,1}};

  private int _n1,_n2;
  private Tensors _st;
  private float[][] _t;
  private int[][] _mark; // samples are marked far, trial, or known
  private TimeHeap2 _heap;

  // Marks used during computation of times. For efficiency, we do not
  // loop over all the marks to clear them before beginning a fast 
  // marching loop. Instead, we simply modify the mark values.
  private int _far = 0; // samples with no time
  private int _trial = 1; // samples in min-heap with a proposed time
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

  /**
   * Updates the time for one sample using known times in neighbor tris.
   * The array kt contains indices in [0,3] of relevant neighbor tris.
   */
  private void updateTime(int j1, int j2, int[] kt) {

    // Elements of structure tensor.
    float[] s = new float[3];
    _st.getTensor(j1,j2,s);
    float s11 = s[0];
    float s12 = s[1];
    float s22 = s[2];

    // The current minimum time.
    float tmin = _t[j2][j1];

    // Initally assume that no computed time will be less than current min.
    boolean smallerTimeFound = false;

    // For all relevant neighbor tris, ...
    for (int it=0; it<kt.length; ++it) {
      int jt = kt[it];

      // Sample indices of vertices X0, X1, and X2 of neighbor tet.
      int j01 = j1;
      int j02 = j2;
      int j11 = j1-K11[jt];
      int j12 = j2-K12[jt];
      int j21 = j1-K21[jt];
      int j22 = j2-K22[jt];

      // All indices must be in bounds.
      if (j11<0 || j11>=_n1) continue;
      if (j12<0 || j12>=_n2) continue;
      if (j21<0 || j21>=_n1) continue;
      if (j22<0 || j22>=_n2) continue;

      // Which neighbor vertices are known? (At least one must be!)
      int m1 = _mark[j12][j11];
      int m2 = _mark[j22][j21];

      // Times T0, T1, and T2 at vertices X0, X1, and X2 of tri.
      float t0 = INFINITY;
      float t1 = _tk[j13][j12][j11];
      float t2 = _tk[j23][j22][j21];
      float t3 = _tk[j33][j32][j31];

      // Use only known times in {T1,T2} to compute candidate time T0.
      if (m1==_known) {
        if (m2==_known) { // use edge 12
          t0 = computeTime(s11,s12,s22,t1-t2,t2,
                           j11-j21,j12-j22,
                           j01-j21,j02-j22);
        } else { // use vertex 1
          t0 = computeTime(s11,s12,s22,t1,j01-j11,j02-j12);
        }
      } else if (m2==_known) { // use vertex 2
        t0 = computeTime(s11,s12,s22,t2,j01-j21,j02-j22);
      }

      // If computed time T0 is smaller than the min time, update the min time.
      if (t0<tmin) {
        tmin = t0;
        smallerTimeFound = true;
      }
    }

    // If a smaller time has been found, ...
    if (smallerTimeFound) {

      // If not been here before, so this sample not already in the min-heap, 
      // then insert this sample into the min-heap, and if the time list is
      // not null, save the current stored time so it can be restored later.
      if (_mark[j2][j1]!=_trial) {
        _mark[j2][j1] = _trial;
        _heap.insert(j1,j2,tmin);
      }

      // Else, simply reduce the time already stored in the min-heap.
      else {
        _heap.reduce(j1,j2,tmin);
      }

      // Store the smaller time.
      _t[j2][j1] = tmin;
    }
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

  private void solve() {

    // Index object used below to check for neighbors in active set.
    Index jj = new Index(-1,-1);

    // Tolerance for convergence.
    float epsilon = 0.01f; // tolerance for convergence

    // While the active set is not empty, ...
    for (int niter=0; !_as.isEmpty(); ++niter) {

      // DEBUG
      if (niter%10==1)
        plot(_t);

      // Copy sample indices from the active set to a simple array.
      // The active set can now change as we loop over samples.
      int n = _as.size();
      _as.keySet().toArray(_ia);
      trace("niter="+niter+" n="+n); // DEBUG

      // For all samples in the active set, ...
      for (int i=0; i<n; ++i) {

        // Sample indices.
        Index ii = _ia[i];
        int i1 = ii.i1;
        int i2 = ii.i2;

        // Current time and new time.
        float ti = _t[i2][i1];
        float gi = g(i1,i2);
        _t[i2][i1] = gi;

        // If the new and old times are close (converged), then ...
        if (ti-gi<ti*epsilon) {

          // For all neighbor samples, ...
          for (int k=0; k<4; ++k) {

            // Neighbor sample indices; skip if out of bounds.
            int j1 = i1+K1[k];
            int j2 = i2+K2[k];
            if (j1<0 || j1>=_n1 || j2<0 || j2>=_n2) 
              continue;

            // If the neighbor is not already in the active set, ...
            jj.i1 = j1;
            jj.i2 = j2;
            if (!_as.containsKey(jj)) {

              // Compute time for the neighbor.
              float gj = g(j1,j2);

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
   * Jeong's fast tests for a valid solution time t0.
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
   */
  private float g(int i1, int i2) {
    boolean jeongTest = true;
    //boolean jeongTest = false;
    float tmin = _t[i2][i1];

    // Get tensor coefficients.
    float[] d = new float[3];
    _st.getTensor(i1,i2,d);
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

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static class ConstantTensors implements FmmSolver2.Tensors {
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
    int n1 = 101;
    int n2 = 101;
    float angle = FLT_PI*110.0f/180.0f;
    float su = 1.000f;
    float sv = 0.001f;
    float cosa = cos(angle);
    float sina = sin(angle);
    float d11 = su*cosa*cosa+sv*sina*sina;
    float d12 = (su-sv)*sina*cosa;
    float d22 = sv*cosa*cosa+su*sina*sina;
    trace("d11="+d11+" d12="+d12+" d22="+d22+" d="+(d11*d22-d12*d12));
    ConstantTensors st = new ConstantTensors(d11,d12,d22);
    FmmSolver2 fs = new FmmSolver2(n1,n2,st);
    fs.zeroAt(1*n1/4,1*n2/4);
    //fs.zeroAt(2*n1/4,2*n2/4);
    fs.zeroAt(3*n1/4,3*n2/4);
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
