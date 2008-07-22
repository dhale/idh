/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import static edu.mines.jtk.util.MathPlus.*;

// for testing
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

/**
 * A time map for 2D image processing. A time map is like a distance map
 * in which the concept of "nearest" is replaced by "least time". A time
 * map differs from a distance map when velocity varies as a function of
 * position and direction. 
 * <p>
 * More precisely, a time map is a generalized distance map in which 
 * distances are defined by a metric tensor field. (Distances defined 
 * this way are also called "geodesic distances.") For a constant field
 * defined by a single tensor S, time t = sqrt(x'Sx), where x' denotes
 * the transpose of a vector x between a point with coordinates x and a
 * reference point x = 0 at which time t = 0. It is easy to verify that 
 * t = t(x) satisfies the eikonal equation grad(t)'*inv(S)*grad(t) = 1.
 * <p>
 * Assuming a more general field in which tensors may vary from sample 
 * to sample, a time-map solves this non-linear differential equation 
 * numerically.
 * <p>
 * The metric tensor S has units of time-squared / length-squared, so
 * that S is a "slowness-squared tensor" or "sloth tensor".
 * <p>
 * A time map maintains for each sample in a 2-D array the minimum
 * traveltime to one or more specified features. When those features are 
 * themselves a subset of image samples, then times stored for samples in 
 * that subset are zero. In this case, a time map may also maintain the 
 * sample indices of the nearest (least-time) features.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.09
 */
public class TimeMap2 {

  /**
   * An interface for classes of sloth tensors. Each tensor is a
   * symmetric positive-definite 2-by-2 matrix {{s11,s12},{s12,s22}}.
   */
  public interface Tensors {

    /**
     * Gets sloth tensor elements for specified indices.
     * @param i1 index for 1st dimension.
     * @param i2 index for 2nd dimension.
     * @param s array {s11,s12,s22} of tensor elements.
     */
    public void getTensor(int i1, int i2, float[] s);
  }

  /**
   * Constructs a time map with constant identity sloth tensors.
   * In this case, time = distance, which is useful for testing.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public TimeMap2(int n1, int n2) {
    this(n1,n2,new IdentityTensors());
  }
  
  /**
   * Constructs a time map for the specified sloth tensor field.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param st the sloth tensors.
   */
  public TimeMap2(int n1, int n2, Tensors st) {
    _n1 = n1;
    _n2 = n2;
    _st = st;
    _tk = new float[n2][n1];
    _k1 = new int[n2][n1];
    _k2 = new int[n2][n1];
    _mark = new int[n2][n1];
    _imin = new int[n2][n1];
    _imax = new int[n2][n1];
    _hmin = new MinTimeHeap(this);
    _hmax = new MaxTimeHeap(this);
    initialize();
  }

  public void insertZerosAt(boolean[][] z) {

    // Build the max-heap of samples to be zeroed.
    _hmax.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        _imax[i2][i1] = -1;
        if (z[i2][i1]) {
          _hmax.insert(i1,i2,_tk[i2][i1]);
        }
      }
    }

    // Zero samples in the max-heap, from farthest to nearest.
    while (!_hmax.isEmpty()) {
      Entry e = _hmax.remove();
      int i1 = e.i1;
      int i2 = e.i2;
      trace("hmax removed: i1="+i1+" i2="+i2+" t="+e.t);
      insertZeroAt(i1,i2);
      plot(_tk,ColorMap.JET);
    }
  }

  /**
   * Updates this time map from the specified sample with time zero.
   * Specifically, updates all samples in this time map that are nearest to
   * the sample with specified indices, when compared with other samples
   * already known. Afterwards, the specified sample will have zero time.
   * time map that are nearest to the specified sample.
   * @param k1 index in 1st dimension of known sample.
   * @param k2 index in 2nd dimension of known sample.
   */
  public void insertZeroAt(int k1, int k2) {
    clearMarks();
    _hmin.clear();
    _mark[k2][k1] = _known;
    _tk[k2][k1] = 0.0f;
    _k1[k2][k1] = k1;
    _k2[k2][k1] = k2;
    updateNabors(k1,k2);
    if (_monitor!=null)
      _monitor.timeSet(k1,k2,k1,k2,_tk[k2][k1]);
    while (!_hmin.isEmpty()) {
      Entry e = _hmin.remove();
      int i1 = e.i1;
      int i2 = e.i2;
      float t = e.t;
      _mark[i2][i1] = _known;
      _hmax.reduce(i1,i2,t);
      if (_monitor!=null)
        _monitor.timeSet(i1,i2,_k1[i2][i1],_k2[i2][i1],_tk[i2][i1]);
      updateNabors(i1,i2);
    }
  }

  // Initialize all times known but invalid.
  public void initialize() {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        _mark[i2][i1] = _known;
        _tk[i2][i1] = TIME_INVALID;
        _k1[i2][i1] = -1;
        _k2[i2][i1] = -1;
      }
    }
  }

  /**
   * Initializes this time map to time zero for specified known samples.
   * Sets all other samples to time unknown.
   * @param k array of flags; true, if time is zero; false, otherwise.
   */
  public void initialize(boolean[][] k) {

    // Initialize times for all known samples and mark the rest unknown.
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (k[i2][i1]) {
          _mark[i2][i1] = _known;
          _tk[i2][i1] = 0.0f;
          _k1[i2][i1] = i1;
          _k2[i2][i1] = i2;
          if (_monitor!=null)
            _monitor.timeSet(i1,i2,_k1[i2][i1],_k2[i2][i1],_tk[i2][i1]);
        } else {
          _mark[i2][i1] = _far;
          _tk[i2][i1] = TIME_INVALID;
        }
      }
    }

    // Initialize the heap with nabors of all known samples.
    _hmin.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_mark[i2][i1]==_known) {
          updateNabors(i1,i2);
        }
      }
    }
  }

  /**
   * Computes times for unknown samples by extrapolating from known times.
   * This method does nothing if all times are known.
   */
  public void extrapolate() {
    while (!_hmin.isEmpty()) {
      Entry e = _hmin.remove();
      int i1 = e.i1;
      int i2 = e.i2;
      _mark[i2][i1] = _known;
      if (_monitor!=null)
        _monitor.timeSet(i1,i2,_k1[i2][i1],_k2[i2][i1],_tk[i2][i1]);
      updateNabors(i1,i2);
    }
  }

  /**
   * Sets a monitor for changes to this time map.
   * @param monitor the monitor; null, if none.
   */
  public void setMonitor(Monitor monitor) {
    _monitor = monitor;
  }

  /**
   * Gets the array of times maintained by this time map.
   * @return the array of times; by reference, not by copy.
   */
  public float[][] getTimes() {
    return _tk;
  }

  public interface Monitor {
    public void timeSet(int i1, int i2, int k1, int k2, float ts);
  }

  // The value for times not yet computed. Also the value returned by
  // methods that compute times from nabor times when a valid time
  // cannot be computed. We use the maximum possible float so that
  // it will be larger than any valid times we compute.
  private static final float TIME_INVALID = Float.MAX_VALUE;

  private int _far = 0;
  private int _trial = 1;
  private int _known = 2;
  private void clearMarks() {
    if (_known+2>Integer.MAX_VALUE) {
      _far = 0;
      _trial = 1;
      _known = 2;
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _mark[i2][i1] = _far;
        }
      }
    } else {
      _far += 2; // all known samples are now far
      _trial +=2; // no trial samples
      _known +=2; // no known samples
    }
  }

  private Tensors _st; // the sloth tensor field
  private int _n1,_n2; // map dimensions
  private float[][] _tk; // time to nearest painted sample
  private int[][] _k1,_k2; // indices of nearest painted sample
  private int[][] _mark; // samples are marked far, trial, or known
  private byte[][] _type; // fixed, extra, 
  private int[][] _imin,_imax; // indices for samples in min/max heaps
  private MinTimeHeap _hmin; // the min heap
  private MaxTimeHeap _hmax; // the max heap
  private Monitor _monitor; // not null if monitoring changes

  // Times for each sample are computed from one of eight nabor triangles.
  // These triangles are indexed as follows:
  //       2 ^
  //   * - - * - - *
  //   | \ 2 | 1 / | 
  //   | 3 \ | / 0 |
  //   * - - X - - * >
  //   | 4 / | \ 7 | 1
  //   | / 5 | 6 \ | 
  //   * - - * - - *
  // The symbol X represents the vertex X0 shared by all eight triangles. 
  // The symbol * represents the other two triangle vertices X1 and X2, 
  // which are indexed in counter-clockwise order around X0.

  // Sample index offsets for vertices X1 of the eight nabor triangles.
  private static final int[] K11 = { 1, 1, 0,-1,-1,-1, 0, 1};
  private static final int[] K12 = { 0, 1, 1, 1, 0,-1,-1,-1};

  // Sample index offsets for vertices X2 of the eight nabor triangles.
  private static final int[] K21 = { 1, 0,-1,-1,-1, 0, 1, 1};
  private static final int[] K22 = { 1, 1, 1, 0,-1,-1,-1, 0};

  // Components of vectors Y1 = X1-X2 for the eight nabor triangles.
  private static final float[] Y11 =
    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,-1.0f,-1.0f, 0.0f};
  private static final float[] Y12 =
    {-1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,-1.0f};

  // Components of vectors Y2 = X0-X2 for the eight nabor triangles.
  private static final float[] Y21 =
    {-1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f,-1.0f,-1.0f};
  private static final float[] Y22 =
    {-1.0f,-1.0f,-1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f};

  // Sloth tensors.
  private static class IdentityTensors implements Tensors {
    public void getTensor(int i1, int i2, float[] s) {
      s[0] = 1.0f; // s11
      s[1] = 0.0f; // s12
      s[2] = 1.0f; // s22
    }
  }

  private void updateNabors(int i1, int i2) {

    // For all eight nabors of (i1,i2) ...
    for (int k=0; k<8; ++k) {
      int k1 = K11[k];
      int k2 = K12[k];

      // Sample indices for this nabor; skip if out of bounds.
      int j1 = i1+k1;
      int j2 = i2+k2;
      if (j1<0 || j1>=_n1) continue;
      if (j2<0 || j2>=_n2) continue;

      // Skip this nabor if time is already known.
      if (_mark[j2][j1]==_known) continue;

      // Current time for this nabor. (May be unknown.)
      float tj = _tk[j2][j1];

      // If nabor not already in the trial heap, insert it.
      if (_mark[j2][j1]!=_trial) {
        _mark[j2][j1] = _trial;
        _hmin.insert(j1,j2,tj);
      }

      // Update time for this nabor.
      updateTime(j1,j2);
    }
  }

  /**
   * Updates the time for one sample using times at eight nabors.
   * @param i1 sample index in 1st dimension at which to compute the time.
   * @param i2 sample index in 2nd dimension at which to compute the time.
   */
  private void updateTime(int i1, int i2) {

    // Elements of sloth tensor.
    float[] s = new float[3];
    _st.getTensor(i1,i2,s);
    float s11 = s[0];
    float s12 = s[1];
    float s22 = s[2];

    // Current time and indices for the specified sample.
    float ti = _tk[i2][i1];
    int ki1 = _k1[i2][i1];
    int ki2 = _k2[i2][i1];

    // For all eight nabor triangles, ...
    for (int it=0; it<8; ++it) {
      int ik1 = -1;
      int ik2 = -1;

      // Sample indices of vertices X1 and X2 of nabor triangle.
      int i11 = i1+K11[it];
      int i12 = i2+K12[it];
      int i21 = i1+K21[it];
      int i22 = i2+K22[it];
      if (i11<0 || i11>=_n1) continue;
      if (i12<0 || i12>=_n2) continue;
      if (i21<0 || i21>=_n1) continue;
      if (i22<0 || i22>=_n2) continue;

      // Need at least one nabor with known time.
      int m1 = _mark[i12][i11];
      int m2 = _mark[i22][i21];
      if (m1!=_known && m2!=_known) continue;

      // Times T0, T1 and T2 at vertices X0, X1 and X2 of nabor triangle.
      float t0 = TIME_INVALID;
      float t1 = _tk[i12][i11];
      float t2 = _tk[i22][i21];

      // Components of vectors Y1 = X1-X2 and Y2 = X0-X2.
      float y11 = Y11[it];
      float y12 = Y12[it];
      float y21 = Y21[it];
      float y22 = Y22[it];

      // Inner products with respect to metric tensor S.
      float d11 = y11*s11*y11+y11*s12*y12+y12*s12*y11+y12*s22*y12;
      float d12 = y11*s11*y21+y11*s12*y22+y12*s12*y21+y12*s22*y22;
      float d22 = y21*s11*y21+y21*s12*y22+y22*s12*y21+y22*s22*y22;

      // Time T0 computed for one nabor triangle.
      if (m1!=_known) {
        t0 = t2+sqrt(d22); // a = 0
        ik1 = i21;
        ik2 = i22;
      } else if (m2!=_known) {
        t0 = t1+sqrt(d22-2.0f*d12+d11); // a = 1
        ik1 = i11;
        ik2 = i12;
      } else {
        float u1 = t1-t2;
        float u2 = t2;
        float dd = d11*d22-d12*d12;
        if (dd<0.0f) dd = 0.0f;
        float du = d11-u1*u1;
        if (du>0.0f) {
          float a = (d12-u1*sqrt(dd/du))/d11;
          if (a<=0.0f) { // a <= 0
            t0 = t2+sqrt(d22);
            ik1 = i21;
            ik2 = i22;
          } else if (a>=1.0f) { // a >= 1
            t0 = t1+sqrt(d22-2.0f*d12+d11);
            ik1 = i11;
            ik2 = i12;
          } else { // 0 < a < 1
            float da = d22-a*(2.0f*d12-a*d11);
            if (da<0.0f) da = 0.0f;
            t0 = u2+a*u1+sqrt(d22-2.0f*a*d12+a*a*d11);
            if (t1<t2) {
              ik1 = i11;
              ik2 = i12;
            } else {
              ik1 = i21;
              ik2 = i22;
            }
          }
        }
      }

      // If computed time T0 is smaller, update the current time.
      if (t0<ti) {
        ti = t0;
        ki1 = _k1[ik2][ik1];
        ki2 = _k2[ik2][ik1];
      }
    }

    // If computed time is smaller, reduce the current time.
    if (ti<_tk[i2][i1]) {
      _tk[i2][i1] = ti;
      _k1[i2][i1] = ki1;
      _k2[i2][i1] = ki2;
      _hmin.reduce(i1,i2,ti);
    }
  }

  // Used by the min-heap to maintain indices of samples in the heap.
  // Each heap entry has array indices (i1,i2). These methods access 
  // and update the corresponding index in the heap. They are necessary 
  // so that times for samples in the heap can be reduced without 
  // searching through the entire heap looking for the sample with 
  // indices (i1,i2).
  private int getMinTimeHeapIndex(int i1, int i2) {
    return _imin[i2][i1];
  }
  private void setMinTimeHeapIndex(Entry e, int i) {
    _imin[e.i2][e.i1] = i;
  }

  // Same for max-heap.
  private int getMaxTimeHeapIndex(int i1, int i2) {
    return _imax[i2][i1];
  }
  private void setMaxTimeHeapIndex(Entry e, int i) {
    _imax[e.i2][e.i1] = i;
  }

  // An entry in a heap has sample indices (i1,i2) and time t.
  private static class Entry {
    int i1,i2;
    float t;
  }

  // A min- or max- heap of times. This heap is special in that it 
  // maintains indices in a corresponding time map. For specified 
  // sample indices (i1,i2), those indices enable O(1) access to heap 
  // entries. Such fast access is important when reducing times in the 
  // time map.
  private static class TimeHeap {

    enum Type {MIN,MAX};

    // Constructs a heap with a corresponding time map.
    TimeHeap(TimeMap2 tmap, Type type) {
      _tmap = tmap;
      _type = type;
    }

    // Inserts a new entry with specified indices and time.
    void insert(int i1, int i2, float t) {
      int i = _n; // index at which to insert the entry
      if (_n==_e.length) // if necessary, ...
        grow(_n+1); // increase the capacity of this heap
      Entry ei = _e[i];
      if (ei==null) // if an entry does not already exist, ...
        ei = new Entry(); // make a new entry
      ei.i1 = i1;
      ei.i2 = i2;
      ei.t = t;
      set(i,ei);
      siftUp(i);
      ++_n;
    }

    // Reduces the time of the entry with specified indices.
    void reduce(int i1, int i2, float t) {
      int i = (_type==Type.MIN) ?
        _tmap.getMinTimeHeapIndex(i1,i2) :
        _tmap.getMaxTimeHeapIndex(i1,i2);
      if (0<=i && i<_n) {
        Entry ei = _e[i];
        ei.t = t;
        set(i,ei);
        siftUp(i);
      }
    }

    // Removes and returns the entry with smallest time.
    Entry remove() {
      Entry e0 = _e[0];
      --_n;
      set(0,_e[_n]);
      set(_n,e0);
      siftDown(0);
      return e0;
    }

    // Removes all entries from this heap.
    void clear() {
      _n = 0;
    }

    // Returns number of entries in this heap.
    int size() {
      return _n;
    }

    // Returns true if this heap is empty; false, otherwise.
    boolean isEmpty() {
      return _n==0;
    }

    private int _n; // number of entries in this heap
    private Entry[] _e = new Entry[1024]; // array of entries in this heap
    private TimeMap2 _tmap; // time map kept in sync with this heap
    private Type _type; // the type of this heap

    // If necessary, moves entry e[i] down so not greater/less than children.
    private void siftDown(int i) {
      Entry ei = _e[i]; // entry ei that may move down
      float eit = ei.t; // cached time for entry ei
      int m = _n>>>1; // number of entries with at least one child
      while (i<m) { // while not childless, ...
        int c = (i<<1)+1; // index of left child
        int r = c+1; // index of right child
        Entry ec = _e[c]; // initially assume left child smallest/largest
        if (_type==Type.MIN) { // if min-heap
          if (r<_n && _e[r].t<ec.t) // if right child smallest, ...
            ec = _e[c=r]; // the smaller of left and right children
          if (eit<=ec.t) // break if ei not greater than smaller child
            break;
        } else { // if max-heap
          if (r<_n && _e[r].t>ec.t) // if right child largest, ...
            ec = _e[c=r]; // the larger of left and right children
          if (eit>=ec.t) // break if ei not less than larger child
            break;
        }
        set(i,ec); // move smaller/larger child up
        i = c;
      }
      if (ei!=_e[i]) // if necessary, ...
        set(i,ei); // set ei where it belongs
    }

    // If necessary, moves entry e[i] up so not less/greater than parent.
    private void siftUp(int i) {
      Entry ei = _e[i]; // entry ei that may move up
      float eit = ei.t; // cached time for entry ei
      while (i>0) { // while a parent (not the root entry), ...
        int p = (i-1)>>>1; // index of parent
        Entry ep = _e[p]; // the parent
        if (_type==Type.MIN) { // if min-heap
          if (eit>=ep.t) // break if ei not less than parent
            break;
        } else {
          if (eit<=ep.t) // break if ei not greater than parent
            break;
        }
        set(i,ep); // ei less/greater than parent, so move parent down
        i = p;
      }
      if (ei!=_e[i]) // if necessary, ...
        set(i,ei); // set ei where it belongs
    }

    // Sets the i'th entry and updates the time map.
    private void set(int i, Entry ei) {
      _e[i] = ei;
      if (_type==Type.MIN)
        _tmap.setMinTimeHeapIndex(ei,i);
      else
        _tmap.setMaxTimeHeapIndex(ei,i);
    }

    // Grows this heap to have at least the specified capacity.
    private void grow(int minCapacity) {
      if (minCapacity<0) // overflow
        throw new OutOfMemoryError();
      int oldCapacity = _e.length;
      int newCapacity = oldCapacity*2;
      if (newCapacity<0) // overflow
        newCapacity = Integer.MAX_VALUE;
      if (newCapacity<minCapacity)
        newCapacity = minCapacity;
      Entry[] e = new Entry[newCapacity];
      System.arraycopy(_e,0,e,0,oldCapacity);
      _e = e;
    }

    // Dumps this heap to the console; leading spaces show level in tree.
    void dump() {
      dump("",0);
    }
    private void dump(String s, int i) {
      if (i<_n) {
        s = s+"  ";
        Entry e = _e[i];
        trace(s+""+e.i1+" "+e.i2+" "+e.t);
        dump(s,2*i+1);
        dump(s,2*i+2);
      }
    }
  }
  private static class MinTimeHeap extends TimeHeap {
    MinTimeHeap(TimeMap2 tmap) {
      super(tmap,Type.MIN);
    }
  }
  private static class MaxTimeHeap extends TimeHeap {
    MaxTimeHeap(TimeMap2 tmap) {
      super(tmap,Type.MAX);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void testMinTimeHeap() {
    int n1 = 5;
    int n2 = 6;
    int n = n1*n2;
    TimeMap2 tmap = new TimeMap2(n1,n2);
    MinTimeHeap hmin = new MinTimeHeap(tmap);
    float[] s = Array.randfloat(n);
    float[][] t = Array.reshape(n1,n2,s);
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        float ti = t[i2][i1];
        hmin.insert(i1,i2,ti);
        s[i] = ti;
      }
    }
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        s[i] -= 0.5f;
        t[i2][i1] -= 0.5f;
        hmin.reduce(i1,i2,t[i2][i1]);
      }
    }
    assert !hmin.isEmpty();
    assert hmin.size()==n;
    Array.quickSort(s);
    for (int i=0; i<n; ++i) {
      Entry e = hmin.remove();
      float ti = e.t;
      assert ti==s[i];
    }
    assert hmin.isEmpty();
    assert hmin.size()==0;
  }

  private static class PlotMonitor implements TimeMap2.Monitor {
    PlotMonitor(TimeMap2 tmap, float tmax) {
      _tmap = tmap;
      _tk = _tmap.getTimes();
      _n1 = _tk[0].length;
      _n2 = _tk.length;
      _tk = Array.fillfloat(tmax,_n1,_n2);
      _sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
      _sp.setSize(650,600);
      _pv = _sp.addPixels(tmap.getTimes());
      _pv.setColorModel(ColorMap.JET);
      //_pv.setColorModel(ColorMap.PRISM);
      _pv.setInterpolation(PixelsView.Interpolation.NEAREST);
      //_pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    }
    public void timeSet(int i1, int i2, int k1, int k2, float tk) {
      _tk[i2][i1] = tk;
      SwingUtilities.invokeLater(new Runnable() {
        public void run() {
          _pv.set(_tk);
        }
      });
      try {
        Thread.currentThread().sleep(100);
      } catch (InterruptedException e) {
      }
    }
    TimeMap2 _tmap;
    private int _n1,_n2;
    private SimplePlot _sp;
    private PixelsView _pv;
    private float[][] _tk;
  }

  private static void plot(float[][] f) {
    plot(f,null);
  }

  private static void plot(float[][] f, IndexColorModel cm) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    //sp.setSize(650,600);
    sp.setSize(1550,1500);
    PixelsView pv = sp.addPixels(f);
    if (cm==null) cm = ColorMap.JET;
    pv.setColorModel(cm);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  }

  private static class EigenTensors implements Tensors {
    EigenTensors(int n1, int n2, double s1, double s2, double v1) {
      float u2 = -(float)v1;
      float u1 = sqrt(1.0f-u2*u2);
      float a1 = (float)s1;
      float a2 = (float)s2;
      _et = new EigenTensors2(n1,n2);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          _et.setEigenvectorU(i1,i2,u1,u2);
          _et.setEigenvalues(i1,i2,a1,a2);
        }
      }
    }
    public void getTensor(int i1, int i2, float[] s) {
      _et.getTensor(i1,i2,s);
    }
    private EigenTensors2 _et;
  }

  private static float[][] timeMapExact(
    double s1, double s2, double v1, boolean[][] k)
  {
    int n1 = k[0].length;
    int n2 = k.length;
    float u1 = (float)sqrt(1.0-v1*v1);
    float u2 = (float)(-v1);
    float a1 = (float)s1;
    float a2 = (float)s2;
    float s11 = a1*u1*u1+a2;
    float s12 = a1*u1*u2   ;
    float s22 = a1*u2*u2+a2;
    int nk = 0;
    int[] k1 = new int[n1*n2];
    int[] k2 = new int[n1*n2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (k[i2][i1]) {
          k1[nk] = i1;
          k2[nk] = i2;
          ++nk;
        }
      }
    }
    float[][] t = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float tmin = Float.MAX_VALUE;
        for (int ik=0; ik<nk; ++ik) {
          float d1 = (float)(i1-k1[ik]);
          float d2 = (float)(i2-k2[ik]);
          float ti = sqrt(d1*s11*d1+2.0f*d1*s12*d2+d2*s22*d2);
          if (ti<tmin) {
            tmin = ti;
          }
        }
        t[i2][i1] = tmin;
      }
    }
    return t;
  }

  private static void testTimeMaps() {
    //double[] angles = {0.0,1*PI/12,2*PI/12,3*PI/12,4*PI/12};
    double su = 100;
    double sv = 1;
    double s1 = su-sv;
    double s2 = sv;
    double[] angles = {PI/8};
    for (double angle:angles)
      testTimeMap(s1,s2,sin(angle));
  }
  private static void testTimeMap(double s1, double s2, double v1) {
    int n1 = 128;
    int n2 = 128;
    boolean[][] known = new boolean[n2][n1];
    //known[n2/2][n1/2] = true;
    //known[   0][   0] = true;
    //known[   0][n1-1] = true;
    //known[n2-1][   0] = true;
    //known[n2-1][n1-1] = true;
    known[1*n2/4][1*n1/4] = true;
    known[3*n2/4][3*n1/4] = true;
    float[][] te = timeMapExact(s1,s2,v1,known);
    float tmax = Array.max(te);
    EigenTensors et = new EigenTensors(n1,n2,s1,s2,v1);
    TimeMap2 tmap = new TimeMap2(n1,n2,et);
    //PlotMonitor pm = new PlotMonitor(tmap,tmax);
    //tmap.setMonitor(pm);
    tmap.initialize(known);
    tmap.extrapolate();
    float[][] tk = tmap.getTimes();
    int[][] k1 = tmap._k1;
    int k1a = 1*n1/4;
    int k1b = 3*n1/4;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (k1[i2][i1]==k1a) {
          tk[i2][i1] = 1.0f;
        } else {
          tk[i2][i1] = 2.0f;
        }
      }
    }
    //diffStats(tk,te);
    plot(tk);
    plot(te);
  }

  private static float[][] getK2(int n1, int n2, TimeMap2 tmap) {
    int[][] k2 = tmap._k2;
    float[][] f2 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f2[i2][i1] = (float)k2[i2][i1];
      }
    }
    return f2;
  }

  private static float[][] readImage(int n1, int n2, String fileName) {
    String dataDir = "/data/seis/joe/";
    try {
      ArrayInputStream ais = new ArrayInputStream(dataDir+fileName);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      return null;
    }
  }

  private static TimeMap2.Tensors getStructureTensors(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    final float[][] u1 = new float[n2][n1];
    final float[][] u2 = new float[n2][n1];
    final float[][] eu = new float[n2][n1];
    final float[][] ev = new float[n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(4);
    lof.apply(x,null,u1,u2,null,null,eu,ev,null);
    final float[][] s1 = Array.div(Array.sub(eu,ev),eu);
    final float[][] s2 = Array.div(ev,eu);
    Array.mul(100.0f,s1,s1);
    return new TimeMap2.Tensors() {
      public void getTensor(int i1, int i2, float[] a) {
        _et.getTensor(i1,i2,a);
      }
      private EigenTensors2 _et = new EigenTensors2(u1,u2,s1,s2);
    };
  }

  private static class LensEigenTensors implements Tensors {
    LensEigenTensors(int n1, int n2, double s1, double s2, double v1) {
      float u2 = -(float)v1;
      float u1 = sqrt(1.0f-u2*u2);
      float a1 = (float)s1;
      float a2 = (float)s2;
      _et = new EigenTensors2(n1,n2);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float d1 = (float)(i1-n1/2);
          float d2 = (float)(i2-n2/2);
          float as = 1.0f-exp(-0.0001f*(d1*d1+d2*d2));
          _et.setEigenvectorU(i1,i2,u1,u2);
          _et.setEigenvalues(i1,i2,a1*as,a2*as);
        }
      }
    }
    public void getTensor(int i1, int i2, float[] s) {
      _et.getTensor(i1,i2,s);
    }
    private EigenTensors2 _et;
  }

  private static void testChannels() {
    int n1 = 200;
    int n2 = 200;
    float[][] x = readImage(n1,n2,"x174.dat");
    plot(x,ColorMap.GRAY);
    TimeMap2.Tensors st = getStructureTensors(x);
    int[] k1 = { 34, 92,172, 27, 25, 12, 81};
    int[] k2 = { 81,109,109,111,124,138,146};
    int nk = k1.length;
    TimeMap2 tmap;

    boolean[][] z = new boolean[n2][n1];
    for (int ik=0; ik<nk; ++ik)
      z[k2[ik]][k1[ik]] = true;
    tmap = new TimeMap2(n1,n2,st);
    tmap.insertZerosAt(z);
    
    tmap = new TimeMap2(n1,n2,st);
    for (int ik=0; ik<nk; ++ik) {
      tmap.insertZeroAt(k1[ik],k2[ik]);
      //plot(tmap.getTimes(),ColorMap.JET);
    }
    plot(tmap.getTimes(),ColorMap.JET);
  }

  private static float[][][] slowWay(int n1, int n2, TimeMap2.Tensors st) {
    float[][] tk = new float[n2][n1];
    Array.fill(Float.MAX_VALUE,tk);
    float[][] k2 = new float[n2][n1];
    boolean[][] known = new boolean[n2][n1];
    for (int j2=0; j2<n2; j2+=1) {
      for (int i2=0; i2<n2; ++i2)
        known[i2][n1-1] = (i2==j2)?true:false;
      TimeMap2 tmap = new TimeMap2(n1,n2,st);
      tmap.initialize(known);
      tmap.extrapolate();
      float[][] tkj = tmap.getTimes();
      int[][] k2j = tmap._k2;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (tkj[i2][i1]<tk[i2][i1]) {
            tk[i2][i1] = tkj[i2][i1];
            k2[i2][i1] = (float)k2j[i2][i1];
          }
        }
      }
    }
    return new float[][][]{tk,k2};
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    testChannels();
    //testTimeMaps();
    /*
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testTimeMaps();
      }
    });
    */
  }
}
