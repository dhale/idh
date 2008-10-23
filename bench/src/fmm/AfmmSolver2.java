/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * A solver for 2D anisotropic eikonal equations. The non-linear equations 
 * are sqrt(grad(t) dot D*grad(t)) = 1, where t is the solution time field, 
 * and D denotes a positive-definite (velocity-squared) tensor field.
 * <p>
 * This solver uses an anisotropic fast marching method like that described 
 * by Konukoglu et al. (2007), but with a different (faster?) local solver.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.30
 */
public class AfmmSolver2 {

  /**
   * A listener for time changes.
   */
  public interface Listener {

    /**
     * Called when time for the specified sample has decreased.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param t the decreased time.
     */
    public void timeDecreased(int i1, int i2, float t);
  }
  
  /**
   * Constructs a solver for the specified tensor field.
   * All times are initially infinite (very large).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors velocity-squared tensors.
   */
  public AfmmSolver2(int n1, int n2, Tensors2 tensors) {
    init(n1,n2,null,tensors);
  }
  
  /**
   * Constructs a solver for a specified array of times.
   * The array is referenced (not copied) by this solver.
   * @param t array of times to be updated by this solver; 
   * @param tensors velocity-squared tensors.
   */
  public AfmmSolver2(float[][] t, Tensors2 tensors) {
    init(t[0].length,t.length,t,tensors);
  }

  /**
   * Sets the tensors used by this solver.
   * @param tensors the tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = tensors;
  }

  /**
   * Zeros the time at the specified sample and computes times for neighbors.
   * Times of neighbor samples are computed recursively while computed times 
   * are less than current times. In other words, this method only decreases 
   * times in the array of times referenced by this class. Finally, this 
   * method notifies any listeners of all samples with times decreased.
   * @param i1 index in 1st dimension of time to zero.
   * @param i2 index in 2nd dimension of time to zero.
   * @return the modified array of times; by reference, not by copy.
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

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.0001f;
  private static final float ONE_MINUS_EPSILON = 1.0f-EPSILON;

  private int _n1,_n2; // dimensions for array of times
  private Tensors2 _tensors; // the tensor field
  private Sample[][] _s; // array of samples
  private float[][] _t; // times for all samples (same as times *in* samples)
  private TimeHeapX<Sample> _htrial; // min-heap of trial samples
  private TimeHeapX<Sample> _hchanged; // min-heap of changed samples
  private ArrayList<Listener> _listeners = new ArrayList<Listener>();
  private ArrayList<Sample> _stack = new ArrayList<Sample>(1024);

  // A sample has indices, a time, and a status.
  private static class Sample implements TimeHeapX.Entry {
    int i1,i2; // sample indices
    float t; // time for this sample (same as _t[i2][i1])
    int status; // far, trial, known or changed
    Sample(int i1, int i2, float t) {
      this.i1 = i1;
      this.i2 = i2;
      this.t = t;
    }
    public int i1() { return i1; }
    public int i2() { return i2; }
    public float t() { return t; }
  }
  private boolean isTrial(Sample s) {
    return s.status==_trial;
  }
  private boolean isKnown(Sample s) {
    return s.status==_known || s.status==_changed;
  }
  private boolean isChanged(Sample s) {
    return s.status==_changed;
  }
  private void setTime(Sample s, float t) {
    _t[s.i2][s.i1] = s.t = t;
  }
  private boolean haveChanged() {
    return !_hchanged.isEmpty();
  }
  private boolean haveTrial() {
    return !_htrial.isEmpty();
  }
  private Sample removeChanged() {
    Sample s = _hchanged.remove();
    s.status = _known;
    return s;
  }
  private Sample removeTrial() {
    Sample s = _htrial.remove();
    s.status = _known;
    return s;
  }
  private void addChanged(Sample s) {
    _hchanged.insert(s);
    s.status = _changed;
  }
  private void addTrial(Sample s) {
    _htrial.insert(s);
    s.status = _trial;
  }
  private void reduceChanged(Sample s) {
    _hchanged.reduce(s);
  }
  private void reduceTrial(Sample s) {
    _htrial.reduce(s);
  }
  private void clearHeaps() {
    _htrial.clear();
    _hchanged.clear();
  }

  // Status of samples during computation of times. For efficiency, we do not
  // loop over all samples to set their status to far before beginning a fast 
  // marching loop. Instead, we simply modify the four possible status values.
  private int _trial = 1; // samples in trial heap
  private int _known = 2; // samples with known time
  private int _changed = 3; // known samples in changed heap
  private void clearStatus() {
    if (_changed+3>Integer.MAX_VALUE) { // rarely!
      _trial = 1;
      _known = 2;
      _changed = 3;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1].status = 0;
        }
      }
    } else { // typically, ...
      _trial +=3; // no trial samples
      _known +=3; // no known samples
      _changed += 3; // no changed samples
    }
  }

  private void init(int n1, int n2, float[][] t, Tensors2 tensors) {
    _n1 = n1;
    _n2 = n2;
    _tensors = tensors;
    _t = (t!=null)?t:Array.fillfloat(INFINITY,n1,n2);
    _s = new Sample[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        _s[i2][i1] = new Sample(i1,i2,INFINITY);
        _t[i2][i1] = INFINITY;
      }
    }
    _htrial = new TimeHeapX<Sample>(TimeHeapX.Type.MIN,n1,n2);
    _hchanged = new TimeHeapX<Sample>(TimeHeapX.Type.MIN,n1,n2);
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

  /*
  private void fireTimesDecreasedFrom(int i1, int i2) {
    Sample si = _s[i2][i1];
    if (!isKnown(si))
      return;
    int nlistener = _listeners.size();
    if (nlistener==0)
      return;
    _stack.clear();
    _stack.add(si);
    setTrial(si);
    while (!_stack.isEmpty()) {
      si = _stack.remove(_stack.size()-1);
      i1 = si.i1;
      i2 = si.i2;
      float ti = _t[i2][i1];
      for (int i=0; i<nlistener; ++i)
        _listeners.get(i).timeDecreased(i1,i2,ti);
      for (int k=0; k<4; ++k) {
        int j1 = i1+K1[k];  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k];  if (j2<0 || j2>=_n2) continue;
        Sample sj = _s[j2][j1];
        if (isKnown(sj)) {
          _stack.add(sj);
          setTrial(sj);
        }
      }
    }
  }
  */


  /**
   * Zeros the time for the specified sample. Recursively updates times 
   * for neighbor samples until all times have converged, and then notifies
   * any listeners of all times decreased.
   */
  private void solveFrom(int i1, int i2) {

    // Clear status of all samples so we can tell which ones become known.
    clearStatus();

    // Clear both the trial and changed heaps.
    clearHeaps();

    // Zero the time for the specified sample and it to the trial heap.
    Sample si = _s[i2][i1];
    setTime(si,0.0f);
    addTrial(si);

    // Work array for tensor elements.
    float[] d = new float[3];

    // Statistics.
    int mtrial = 0;
    int ntrial = 0;
    int mchanged = 0;
    int nchanged = 0;
    float tplot = 500.0f;

    // Complete the solve by processing the heaps until empty.
    while (haveTrial() || haveChanged()) {
      mtrial = max(mtrial,_htrial.size());
      mchanged = max(mchanged,_hchanged.size());

      // Get sample from changed heap if not empty, else trial heap.
      si = haveChanged()?removeChanged():removeTrial();
      i1 = si.i1;
      i2 = si.i2;

      // For all neighbor samples, compute a new time. Update the changed
      // heap for any known samples with reduced time. Update the trial
      // heap for any unknown samples with reduced time. We can handle
      // both known and unknown neighbors with one loop, because times
      // for any one neighbor do not affect times for other neighbors.
      for (int k=0; k<4; ++k) {
        int j1 = i1+K1[k]; if (j1<0 || j1>=_n1) continue;
        int j2 = i2+K2[k]; if (j2<0 || j2>=_n2) continue;
        Sample sj = _s[j2][j1];
        float tj = sj.t;
        float cj = computeTime(j1,j2,K1S[k],K2S[k],d);
        if (isKnown(sj)) {
          if (cj<tj*ONE_MINUS_EPSILON) {
            setTime(sj,cj);
            if (isChanged(sj)) {
              reduceChanged(sj);
            } else {
              addChanged(sj);
              ++nchanged;
            }
          }
        } else {
          if (cj<tj) {
            setTime(sj,cj);
            if (isTrial(sj)) {
              reduceTrial(sj);
            } else {
              addTrial(sj);
              ++ntrial;
            }
          }
        }
      }
    }
    int ntotal = ntrial+nchanged;
    trace("mtrial="+mtrial+" mchanged="+mchanged);
    trace("ntrial="+ntrial+" nchanged="+nchanged);
    trace("ntotal="+ntotal+" nratio="+((double)ntotal/(double)(_n1*_n2)));

    // Notify any listeners of all times decreased.
    //fireTimesDecreasedFrom(i1,i2);
  }

  // Methods to get times for neighbors.
  private float t1m(int i1, int i2) {
    return (--i1>=0 && isKnown(_s[i2][i1]))?_t[i2][i1]:INFINITY;
  }
  private float t1p(int i1, int i2) {
    return (++i1<_n1 && isKnown(_s[i2][i1]))?_t[i2][i1]:INFINITY;
  }
  private float t2m(int i1, int i2) {
    return (--i2>=0 && isKnown(_s[i2][i1]))?_t[i2][i1]:INFINITY;
  }
  private float t2p(int i1, int i2) {
    return (++i2<_n2 && isKnown(_s[i2][i1]))?_t[i2][i1]:INFINITY;
  }

  /**
   * Returns a time t not greater than the current time for one sample.
   * Computations are limited to neighbor samples with specified indices.
   */
  private float computeTime(int i1, int i2, int[] k1s, int[] k2s, float[] d) {
    _tensors.getTensor(i1,i2,d);
    float d11 = d[0];
    float d12 = d[1];
    float d22 = d[2];
    float e12 = 1.0f/(d11*d22-d12*d12);
    float tc = _t[i2][i1];
    float t1m = t1m(i1,i2);
    float t1p = t1p(i1,i2);
    float t2m = t2m(i1,i2);
    float t2p = t2p(i1,i2);
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
        //tc = t0;
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
  private static void trace(int i1, int i2, String s) {
    //if (i1==2 && i2==2)
    //  trace(s);
  }
}
