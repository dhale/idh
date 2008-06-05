/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * 2D image painting using fast marching methods.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.05.27
 */
public class ImagePainter2 {

  public ImagePainter2(float[][] image) {
    _n1 = image[0].length;
    _n2 = image.length;
    _image = image;
    _tmap = new TimeMap(_n1,_n2);
  }

  private void initialize(float[][] paint, byte[][] flags) {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (flags[i2][i1]==TimeMap.FIXED) {
          
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2;
  private float[][] _image;
  private float[][] _paint;
  private byte[][] _flags;
  private TimeMap _tmap;

  // A time map is a distance map for which "nearest" means least time.
  private static class TimeMap {

    int _n1,_n2; // map dimensions
    float[][] _tk; // time to nearest painted (known) sample
    int[][] _k1,_k2; // indices of nearest painted (known) sample
    byte[][] _mark; // samples are marked FAR, TRIAL or KNOWN
    int[][] _imin,_imax; // indices for samples in min/max heaps
    MinTimeHeap _hmin; // the min heap

    TimeMap(int n1, int n2) {
      _n1 = n1;
      _n2 = n2;
      _tk = new float[n2][n1];
      _k1 = new int[n2][n1];
      _k2 = new int[n2][n1];
      _mark = new byte[n2][n1];
      _imin = new int[n2][n1];
      _imax = new int[n2][n1];
      _hmin = new MinTimeHeap(this);
    }
    void initialize(byte[][] flags) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (flags[i2][i1]==FIXED) {
            _mark[i2][i1] = KNOWN;
            _tk[i2][i1] = 0.0f;
            _k1[i2][i1] = i1;
            _k2[i2][i1] = i2;
          } else {
            _mark[i2][i1] = FAR;
            _tk[i2][i1] = TIME_UNKNOWN;
          }
        }
      }
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (_mark[i2][i1]==KNOWN) {
            initializeNabors(i1,i2);
          }
        }
      }
    }
    void extrapolate() {
      trace("extrapolate: heap size="+_hmin.size());
      while (!_hmin.isEmpty()) {
        Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        _mark[i2][i1] = KNOWN;
        updateNabors(i1,i2);
      }
    }
    private static final int[] K1 = { 1,-1,-1, 1, 1, 0,-1, 0};
    private static final int[] K2 = { 1, 1,-1,-1, 0, 1, 0,-1};
    void initializeNabors(int i1, int i2) {
      for (int k=0; k<8; ++k) {
        int k1 = K1[k];
        int k2 = K2[k];

        // Sample indices for this nabor; skip if out of bounds.
        int j1 = i1+k1;
        int j2 = i2+k2;
        if (j1<0 || j1>=_n1) continue;
        if (j2<0 || j2>=_n2) continue;

        // Skip this nabor if time is already known.
        if (_mark[j2][j1]==KNOWN) continue;

        // Compute time for this nabor.
        float e11 = 1.0f; // TODO: tensor coefficients
        float e12 = 0.0f;
        float e22 = 1.0f;
        float y1 = (float)(j1-i1);
        float y2 = (float)(j2-i2);
        float tj = _tk[j2][j1];
        float tc = sqrt(y1*e11*y1+y1*e12*y2+y2*e12*y1+y2*e22*y2);
        trace("  j1="+j1+" j2="+j2+" tj="+tj+" tc="+tc);
        _tk[j2][j1] = tc;
        if (_mark[j2][j1]!=TRIAL) {
          trace("  inserting tc="+tc);
          _hmin.insert(j1,j2,tc);
        } else {
          trace("  reducing tj="+tj+" to tc="+tc);
          _hmin.reduce(j1,j2,tc);
        }
      }
    }
    void updateNabors(int i1, int i2) {
      trace("updateNabors: i1="+i1+" i2="+i2);

      // For all eight nabors of (i1,i2) ...
      for (int k=0; k<8; ++k) {
        int k1 = K1[k];
        int k2 = K2[k];

        // Sample indices for this nabor; skip if out of bounds.
        int j1 = i1+k1;
        int j2 = i2+k2;
        if (j1<0 || j1>=_n1) continue;
        if (j2<0 || j2>=_n2) continue;

        // Skip this nabor if time is already known.
        if (_mark[j2][j1]==KNOWN) continue;

        // Current time for this nabor.
        float tj = _tk[j2][j1];

        // If nabor not already in the trial heap, insert it.
        if (_mark[j2][j1]==FAR) {
          _mark[j2][j1] = TRIAL;
          _hmin.insert(j1,j2,tj);
        }

        // Compute time for this nabor.
        float e11 = 1.0f; // TODO: tensor coefficients
        float e12 = 0.0f;
        float e22 = 1.0f;
        float tc = computeTime(j1,j2,e11,e12,e22,_tk);
        trace("  j1="+j1+" j2="+j2+" tc="+tc);

        // If computed time is smaller, reduce the current time.
        if (tc<tj) {
          trace("  reducing tj="+tj+" to tc="+tc);
          _tk[j2][j1] = tc;
          _hmin.reduce(j1,j2,tc);
        }
      }
    }
    int getMinTimeHeapIndex(int i1, int i2) {
      return _imin[i2][i1];
    }
    void setMinTimeHeapIndex(Entry e, int i) {
      _imin[e.i2][e.i1] = i;
    }

    // The value for times not yet computed. Also the value returned by
    // methods that compute times from nabor times when a valid time
    // cannot be computed. We use the maximum possible float so that
    // it will be larger than any valid times we compute.
    private static final float TIME_UNKNOWN = Float.MAX_VALUE;

    private static final int CLEAR = 0;
    private static final int FIXED = 1;
    private static final int EXTRA = 2;
    private static final int INTER = 3;

    private static final int FAR = 0;
    private static final int TRIAL = 1;
    private static final int KNOWN = 2;

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

    /**
     * Computes the time for one sample using times at eight nabors.
     * @param i1 sample index in 1st dimension at which to compute the time.
     * @param i2 sample index in 2nd dimension at which to compute the time.
     * @param e11 sloth tensor coefficient (1,1) at sample (i1,i2).
     * @param e12 sloth tensor coefficient (1,2) at sample (i1,i2).
     * @param e22 sloth tensor coefficient (2,2) at sample (i1,i2).
     * @param t array of times; referenced but not modified
     * @return the computed time.
     */
    private static float computeTime(
      int i1, int i2, float e11, float e12, float e22, float[][] t) 
    {
      trace("computeTime: i1="+i1+" i2="+i2);
      int n1 = t[0].length;
      int n2 = t.length;

      // Current time for the specified sample.
      float ti = t[i2][i1];

      // For all eight nabor triangles, ...
      for (int it=0; it<8; ++it) {
        trace("  it="+it);

        // Sample indices of vertices X0, X1 and X2 of nabor triangle.
        int i01 = i1;
        int i02 = i2;
        int i11 = i01+K11[it];
        int i12 = i02+K12[it];
        int i21 = i01+K21[it];
        int i22 = i02+K22[it];
        if (i11<0 || i11>=n1) continue;
        if (i12<0 || i12>=n2) continue;
        if (i21<0 || i21>=n1) continue;
        if (i22<0 || i22>=n2) continue;

        // Times T1 and T2 at vertices X1 and X2 of nabor triangle.
        float t0 = TIME_UNKNOWN;
        float t1 = t[i12][i11];
        float t2 = t[i22][i21];
        if (t1==TIME_UNKNOWN && t2==TIME_UNKNOWN) continue;

        // Components of vectors Y1 = X1-X2 and Y2 = X0-X2.
        float y11 = Y11[it];
        float y12 = Y12[it];
        float y21 = Y21[it];
        float y22 = Y22[it];

        // Dot products with respect to tensor E.
        float s11 = y11*e11*y11+y11*e12*y12+y12*e12*y11+y12*e22*y12;
        float s12 = y11*e11*y21+y11*e12*y22+y12*e12*y21+y12*e22*y22;
        float s22 = y21*e11*y21+y21*e12*y22+y22*e12*y21+y22*e22*y22;

        // Time T0 computed for one nabor triangle.
        if (t1==TIME_UNKNOWN) {
          t0 = t2+sqrt(s22); // a = 0
          trace("  t1 unknown: t0="+t0);
        } else if (t2==TIME_UNKNOWN) {
          t0 = t1+sqrt(s22-2.0f*s12+s11); // a = 1
          trace("  t2 unknown: t0="+t0);
        } else {
          trace("  t1 and t2 known");
          float u1 = t1-t2;
          float u2 = t2;
          float ss = s11*s22-s12*s12;
          if (ss<0.0f) ss = 0.0f;
          float su = s11-u1*u1;
          if (su>0.0f) {
            float a = (s12-u1*sqrt(ss/su))/s11;
            if (a<=0.0f) { // a <= 0
              t0 = t2+sqrt(s22);
              trace("    a <= 0: t0="+t0);
            } else if (a>=1.0f) { // a >= 1
              t0 = t1+sqrt(s22-2.0f*s12+s11);
              trace("    a >= 1: t0="+t0);
            } else { // 0 < a < 1
              float sa = s22-a*(2.0f*s12-a*s11);
              if (sa<0.0f) sa = 0.0f;
              t0 = u2+a*u1+sqrt(s22-2.0f*a*s12+a*a*s11);
              trace("    0 < a < 1: t0="+t0);
            }
          } else {
            trace("    su="+su);
          }
        }

        // If computed time T0 is smaller, update the current time.
        if (t0<ti) ti = t0;
      }

      return ti;
    }
  }

  // An entry in a min-heap or max-heap.
  private static class Entry {
    int i1,i2;
    float t;
  }

  // A min-heap of times. This heap is special in that it maintains
  // indices in a corresponding time map. For specified sample indices 
  // (i1,i2), those indices enable O(1) access to heap entries. Such
  // fast access is important when reducing times in the time map.
  private static class MinTimeHeap {

    private int _n; // number of entries in this heap
    private Entry[] _e = new Entry[8]; // array of entries in this heap
    private TimeMap _tmap; // time map kept in sync with this heap

    // Constructs a heap with a corresponding time map.
    MinTimeHeap(TimeMap tmap) {
      _tmap = tmap;
    }

    // Inserts a new entry with specified indices and time.
    void insert(int i1, int i2, float t) {
      int i = _n; // index at which to insert the entry
      if (_n>=_e.length) // if necessary, ...
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
      int i = _tmap.getMinTimeHeapIndex(i1,i2);
      Entry ei = _e[i];
      ei.t = t;
      set(i,ei);
      siftUp(i);
    }

    // Removes and returns the entry with smallest time.
    Entry remove() {
      Entry e0 = _e[0];
      --_n;
      set(0,_e[_n]);
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

    // If necessary, moves entry e[i] down so not greater than children.
    private void siftDown(int i) {
      Entry ei = _e[i]; // entry ei that may move down
      float eit = ei.t; // cached time for entry ei
      int m = _n>>>1; // number of entries with at least one child
      while (i<m) { // while not childless, ...
        int c = (i<<1)+1; // index of left child
        int r = c+1; // index of right child
        Entry ec = _e[c]; // initially assume left child smallest
        if (r<_n && _e[r].t<ec.t) // but if right child smallest, ...
          ec = _e[c=r]; // the smaller of left and right children
        if (eit<=ec.t) // break if ei not greater than smaller child
          break;
        set(i,ec); // ei greater than smaller child, so move smaller child up
        i = c;
      }
      if (ei!=_e[i]) // if necessary, ...
        set(i,ei); // set ei where it belongs
    }

    // If necessary, moves entry e[i] up so not less than parent.
    private void siftUp(int i) {
      Entry ei = _e[i]; // entry ei that may move up
      float eit = ei.t; // cached time for entry ei
      while (i>0) { // while a parent (not the root entry), ...
        int p = (i-1)>>>1; // index of parent
        Entry ep = _e[p]; // the parent
        if (eit>=ep.t) // break if ei not less than parent
          break;
        set(i,ep); // ei less than parent, so move parent down
        i = p;
      }
      if (ei!=_e[i]) // if necessary, ...
        set(i,ei); // set ei where it belongs
    }

    // Sets the i'th entry and updates the time map.
    private void set(int i, Entry ei) {
      _e[i] = ei;
      _tmap.setMinTimeHeapIndex(ei,i);
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
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void plot(float[][] f) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setSize(650,600);
    PixelsView pv = sp.addPixels(f);
    pv.setColorModel(ColorMap.JET);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  }

  private static void testMinTimeHeap() {
    int n1 = 5;
    int n2 = 6;
    int n = n1*n2;
    TimeMap tmap = new TimeMap(n1,n2);
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
      //trace("ti="+ti+" si="+s[i]);
      assert ti==s[i];
    }
    assert hmin.isEmpty();
    assert hmin.size()==0;
  }

  private static void testTimeMap() {
    int n1 = 11;
    int n2 = 11;
    TimeMap tmap = new TimeMap(n1,n2);
    byte[][] flags = new byte[n2][n1];
    flags[n2/2][n1/2] = TimeMap.FIXED;
    tmap.initialize(flags);
    //tmap.extrapolate();
    plot(tmap._tk);
  }

  private static void trace(String s) {
    //System.out.println(s);
  }

  public static void main(String[] args) {
    //testMinTimeHeap();
    testTimeMap();
  }
}
