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
 * A 2D array of painted values.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.12
 */
public class Painting2 {

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
   * Constructs a painting with constant identity sloth tensors.
   * In this case, time = distance, which is useful for testing.
   * Painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param nv number of values painted for each sample
   */
  public Painting2(int n1, int n2, int nv) {
    this(n1,n2,nv,new IdentityTensors());
  }
  
  /**
   * Constructs a painting for the specified sloth tensor field.
   * All painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param nv number of values painted for each sample
   * @param st the sloth tensors.
   */
  public Painting2(int n1, int n2, int nv, Tensors st) {
    _n1 = n1;
    _n2 = n2;
    _nv = nv;
    _st = st;
    _k1 = new int[n2][n1];
    _k2 = new int[n2][n1];
    _tk = new float[n2][n1];
    _vk = new float[n2][n1][];
    _type = new byte[n2][n1];
    _mark = new int[n2][n1];
    _imin = new int[n2][n1];
    _imax = new int[n2][n1];
    _hmin = new MinTimeHeap(this);
    _hmax = new MaxTimeHeap(this);
    clear();
  }

  /**
   * Clears all painted values, including all fixed values.
   */
  public void clear() {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        _type[i2][i1] = CLEAR;
        _mark[i2][i1] = _known;
        _imin[i2][i1] = -1;
        _imax[i2][i1] = -1;
        _k1[i2][i1] = -1;
        _k2[i2][i1] = -1;
        _tk[i2][i1] = TIME_INVALID;
        _vk[i2][i1] = null;
      }
    }
  }

  /**
   * Paints the specified sample with one specified value at index zero.
   * Paints zero values for indices greater than zero.
   * @param k1 index in 1st dimension of painted sample.
   * @param k2 index in 2nd dimension of painted sample.
   * @param vk value at index zero for the painted sample.
   */
  public void paintAt(int k1, int k2, float vk) {
    _type[k2][k1] = FIXED;
    _k1[k2][k1] = k1;
    _k2[k2][k1] = k2;
    _tk[k2][k1] = TIME_INVALID;
    _vk[k2][k1] = new float[_nv];
    _vk[k2][k1][0] = vk;
  }

  /**
   * Paints the specified sample with specified values.
   * After painting, the specified sample is fixed; its values will
   * not change in any subsequent extrapolation or interpolation.
   * @param k1 index in 1st dimension of painted sample.
   * @param k2 index in 2nd dimension of painted sample.
   * @param vk array of values for painted sample.
   */
  public void paintAt(int k1, int k2, float[] vk) {
    _type[k2][k1] = FIXED;
    _k1[k2][k1] = k1;
    _k2[k2][k1] = k2;
    _tk[k2][k1] = TIME_INVALID;
    _vk[k2][k1] = new float[_nv];
    float[] vs = _vk[k2][k1];
    for (int iv=0; iv<_nv; ++iv)
      vs[iv] = vk[iv];
  }

  /**
   * Extrapolates values from all fixed (explicitly painted) samples.
   * After extrapolation, all samples are either fixed or extrapolated.
   */
  public void extrapolate() {

    // Clear all samples that are not fixed, and insert all fixed 
    // samples into the max-heap with huge (invalid) times that can
    // only get smaller. After the max-heap is built, any one of the
    // fixed samples could be at the top of the heap.
    _hmax.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        _tk[i2][i1] = TIME_INVALID;
        _imax[i2][i1] = -1;
        if (_type[i2][i1]==FIXED) {
          _hmax.insert(i1,i2,TIME_INVALID);
        } else {
          _type[i2][i1] = CLEAR;
        }
      }
    }

    // Extrapolate from all fixed samples, one at a time, in order of
    // decreasing time. The fixed sample with the largest time is at the
    // top of the max-heap. As we extrapolate from this fixed sample, times 
    // for other fixed samples that remain in the heap may be reduced, so 
    // that their order may change. We choose a decreasing order to reduce 
    // the number of samples that must be modified during extrapolation 
    // from each fixed sample. 
    while (!_hmax.isEmpty()) {

      // Remove from the max-heap the fixed sample with largest time.
      Entry ef = _hmax.remove();
      int k1 = ef.i1;
      int k2 = ef.i2;

      // The values to be extrapolated.
      float[] vk = _vk[k2][k1];

      // Mark all samples as far, mark the fixed sample as known with 
      // time zero, and update its neighbors.
      clearMarks();
      _hmin.clear();
      _mark[k2][k1] = _known;
      _k1[k2][k1] = k1;
      _k2[k2][k1] = k2;
      _tk[k2][k1] = 0.0f;
      updateNabors(k1,k2);

      // Extrapolate from the fixed sample to all samples that are
      // nearer to the fixed sample than to any other fixed sample.
      while (!_hmin.isEmpty()) {
        Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        float t = e.t;
        _mark[i2][i1] = _known;
        if (_type[i2][i1]==FIXED) {
          _hmax.reduce(i1,i2,t);
        } else {
          _type[i2][i1] = EXTRA;
          _k1[i2][i1] = k1;
          _k2[i2][i1] = k2;
          _vk[i2][i1] = vk;
        }
        updateNabors(i1,i2);
      }
      //plot(_tk); // DEBUG
    }
  }

  /**
   * Gets the array of times used in this painting.
   * These times are modified by extrapolation and interpolation.
   * @return array of times; by reference, not by copy.
   */
  public float[][] getTimes() {
    return _tk;
  }

  /**
   * Gets a copy of the painted values with index zero.
   * @return array of values; by copy, not by reference.
   */
  public float[][] getValues() {
    return getValues(0);
  }

  /**
   * Gets a copy of the painted values with specified index.
   * Zero values are returned for any samples not yet painted.
   * @param iv index of values to get.
   * @return array of values; by copy, not by reference.
   */
  public float[][] getValues(int iv) {
    float[][] v = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        float[] vk = _vk[i2][i1];
        v[i2][i1] = (vk!=null)?vk[iv]:0.0f;
      }
    }
    return v;
  }

  // The value for times not yet computed. Also the value returned by
  // methods that compute times from nabor times when a valid time
  // cannot be computed. We use the maximum possible float so that
  // it will be larger than any valid times we compute.
  private static final float TIME_INVALID = Float.MAX_VALUE;

  // Type of paint.
  private static final byte CLEAR = 0; // values null (not painted)
  private static final byte FIXED = 1; // values painted explicitly
  private static final byte EXTRA = 2; // values painted by extrapolation
  private static final byte INTER = 3; // values painted by interpolation

  // Marks used during computation of times.
  private int _far = 0; // samples with no time
  private int _trial = 1; // samples (in min-heap) with a proposed time
  private int _known = 2; // samples with a known time
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
      _far += 2; // known samples become far samples
      _trial +=2; // no trial samples
      _known +=2; // no known samples
    }
  }

  private Tensors _st; // the sloth tensor field
  private int _n1,_n2; // painting dimensions
  private int _nv; // number of painted values
  private float[][][] _vk; // painted values
  private float[][] _tk; // time to nearest painted sample
  private int[][] _k1,_k2; // indices of nearest painted sample
  private int[][] _mark; // samples are marked far, trial, or known
  private byte[][] _type; // fixed, extra, 
  private int[][] _imin,_imax; // indices for samples in min/max heaps
  private MinTimeHeap _hmin; // the min heap
  private MaxTimeHeap _hmax; // the max heap

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

    // Time for the specified sample.
    float ti = _tk[i2][i1];

    // For all eight nabors of specified sample at (i1,i2) ...
    for (int k=0; k<8; ++k) {
      int k1 = K11[k];
      int k2 = K12[k];

      // Sample indices (j1,j2) for this nabor; skip if out of bounds.
      int j1 = i1+k1;
      int j2 = i2+k2;
      if (j1<0 || j1>=_n1) continue;
      if (j2<0 || j2>=_n2) continue;

      // If time for nabor not already known, update it.
      if (_mark[j2][j1]!=_known)
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

    // The current minimum time.
    float tmin = _tk[i2][i1];

    // Initally assume that no computed time will be less than current min.
    boolean reduced = false;

    // For all eight nabor triangles, ...
    for (int it=0; it<8; ++it) {

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
      } else if (m2!=_known) {
        t0 = t1+sqrt(d22-2.0f*d12+d11); // a = 1
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
          } else if (a>=1.0f) { // a >= 1
            t0 = t1+sqrt(d22-2.0f*d12+d11);
          } else { // 0 < a < 1
            float da = d22-a*(2.0f*d12-a*d11);
            if (da<0.0f) da = 0.0f;
            t0 = u2+a*u1+sqrt(d22-2.0f*a*d12+a*a*d11);
          }
        }
      }

      // If computed time T0 is smaller than the min time, update the min time.
      if (t0<tmin) {
        tmin = t0;
        reduced = true;
      }
    }

    // If the minimum time has been reduced, ...
    if (reduced) {

      // Remember the minimum time.
      _tk[i2][i1] = tmin;

      // If this sample not already in the min-heap, insert it.
      if (_mark[i2][i1]!=_trial) {
        _mark[i2][i1] = _trial;
        _hmin.insert(i1,i2,tmin);
      } 

      // else, reduce the time already stored in the min-heap.
      else {
        _hmin.reduce(i1,i2,tmin);
      }
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
  // maintains indices in a corresponding painting. For specified 
  // sample indices (i1,i2), those indices enable O(1) access to heap 
  // entries. Such fast access is important when reducing times in the 
  // painting.
  private static class TimeHeap {

    enum Type {MIN,MAX};

    // Constructs a heap with a corresponding painting.
    TimeHeap(Painting2 painting, Type type) {
      _painting = painting;
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
        _painting.getMinTimeHeapIndex(i1,i2) :
        _painting.getMaxTimeHeapIndex(i1,i2);
      if (0<=i && i<_n) {
        Entry ei = _e[i];
        ei.t = t;
        set(i,ei);
        if (_type==Type.MIN) {
          siftUp(i);
        } else {
          siftDown(i);
        }
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
    private Painting2 _painting; // painting kept in sync with this heap
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

    // Sets the i'th entry and updates the painting.
    private void set(int i, Entry ei) {
      _e[i] = ei;
      if (_type==Type.MIN)
        _painting.setMinTimeHeapIndex(ei,i);
      else
        _painting.setMaxTimeHeapIndex(ei,i);
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
    MinTimeHeap(Painting2 painting) {
      super(painting,Type.MIN);
    }
  }
  private static class MaxTimeHeap extends TimeHeap {
    MaxTimeHeap(Painting2 painting) {
      super(painting,Type.MAX);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void testMinTimeHeap() {
    int n1 = 5;
    int n2 = 6;
    int n = n1*n2;
    Painting2 painting = new Painting2(n1,n2,0);
    MinTimeHeap hmin = new MinTimeHeap(painting);
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

  private static void plot(float[][] f) {
    plot(f,null);
  }

  private static void plot(float[][] f, IndexColorModel cm) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    //sp.setSize(650,600);
    sp.setSize(1050,1000);
    PixelsView pv = sp.addPixels(f);
    if (cm==null) cm = ColorMap.JET;
    pv.setColorModel(cm);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
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
          float as = exp(-0.0001f*(d1*d1+d2*d2));
          _et.setEigenvectorU(i1,i2,u1,u2);
          _et.setCoefficients(i1,i2,a1*as,a2*as);
        }
      }
    }
    public void getTensor(int i1, int i2, float[] s) {
      _et.getTensor(i1,i2,s);
    }
    private EigenTensors2 _et;
  }

  private static Painting2.Tensors getStructureTensors(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    final float[][] u1 = new float[n2][n1];
    final float[][] u2 = new float[n2][n1];
    final float[][] eu = new float[n2][n1];
    final float[][] ev = new float[n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(4);
    lof.apply(x,null,u1,u2,null,null,eu,ev,null);
    //final float[][] s1 = Array.div(Array.sub(eu,ev),eu);
    //final float[][] s2 = Array.div(ev,eu);
    final float[][] s1 = Array.sub(eu,ev);
    final float[][] s2 = Array.copy(ev);
    Array.mul(10.0f,s1,s1);
    return new Painting2.Tensors() {
      public void getTensor(int i1, int i2, float[] a) {
        _et.getTensor(i1,i2,a);
      }
      private EigenTensors2 _et = new EigenTensors2(u1,u2,s1,s2);
    };
  }

  private static void testChannels() {
    int n1 = 200;
    int n2 = 200;
    int nv = 1;
    float[][] x = readImage(n1,n2,"x174.dat");
    plot(x,ColorMap.GRAY);
    Painting2.Tensors st = getStructureTensors(x);
    //Painting2.Tensors st = new LensEigenTensors(n1,n2,0.0,1.0,1.0);

    int[] k1 =    { 34,  92, 172,  27,  25,  12,  81, 117,  94,  14,  44};
    int[] k2 =    { 81, 109, 109, 111, 124, 138, 146,  82, 122,  99, 162};
    float[] vk = {1.0f,2.0f,2.0f,2.0f,2.0f,3.0f,3.0f,0.0f,0.0f,0.0f,0.0f};
    int nk = vk.length;
    /*
    int nk = n2;
    int[] k1 = new int[nk];
    int[] k2 = new int[nk];
    float[] vk = new float[nk];
    for (int i2=0,ik=0; i2<n2; i2+=1,++ik) {
      k1[ik] = n1-1;
      k2[ik] = i2;
      vk[ik] = (float)i2;
    }
    int[] k1 =   {  n1-1,  n1-1};
    int[] k2 =   {1*n2/4,3*n2/4};
    float[] vk = {  1.0f,  2.0f};
    int nk = vk.length;
    */
    
    Painting2 p = new Painting2(n1,n2,nv,st);
    for (int ik=0; ik<nk; ++ik) {
      p.paintAt(k1[ik],k2[ik],vk[ik]);
    }
    p.extrapolate();
    plot(p.getTimes(),ColorMap.JET);
    plot(p.getValues(),ColorMap.JET);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testChannels();
      }
    });
  }
}
