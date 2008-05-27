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

/**
 * 2D image painting using fast marching methods.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.05.27
 */
public class ImagePainter2 {

  // A time map is a distance map for which "nearest" means least time.
  private class TimeMap {
    float[][] tk; // time to nearest painted (known) sample
    int[][] k1,k2; // indices of nearest painted (known) sample
    int[][] imin,imax; // indices for samples in min/max heaps
    int getMinTimeHeapIndex(int i1, int i2) {
      return imin[i2][i1];
    }
    void setMinTimeHeapIndex(Entry e, int i) {
      imin[e.i2][e.i1] = i;
    }
  }

  private float[][] _image;
  private float[][] _paint;
  private byte[][] _flags;
  private TimeMap _tmap;

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

    int n; // number of entries in this heap
    Entry[] e = new Entry[64]; // array of entries in this heap
    TimeMap tmap; // time map kept in sync with this heap

    // Constructs a heap with a corresponding time map.
    MinTimeHeap(TimeMap timeMap) {
      tmap = timeMap;
    }

    // Inserts a new entry with specified indices and time.
    void insert(int i1, int i2, float t) {
      int i = n; // index at which to insert the entry
      if (n>=e.length) // if necessary, ...
        grow(n+1); // increase the capacity of this heap
      Entry ei = e[i];
      if (ei==null) // if an entry does not already exist, ...
        ei = new Entry(); // make a new entry
      ei.i1 = i1;
      ei.i2 = i2;
      ei.t = t;
      set(i,e[i]);
      siftUp(i);
      ++n;
    }

    // Reduces the time of the entry with specified indices.
    void reduce(int i1, int i2, float t) {
      int i = tmap.getMinTimeHeapIndex(i1,i2);
      Entry ei = e[i];
      ei.t = t;
      set(i,ei);
      siftUp(i);
    }

    // Removes and returns the entry with smallest time.
    Entry remove() {
      Entry e0 = e[0];
      --n;
      set(0,e[n]);
      siftDown(0);
      return e0;
    }

    // Removes all entries from this heap.
    void clear() {
      n = 0;
    }

    // Returns true if this heap is empty; false, otherwise.
    boolean isEmpty() {
      return n==0;
    }

    // If necessary, moves entry e[i] down so not greater than children.
    private void siftDown(int i) {
      Entry ei = e[i]; // entry ei that may move down
      float eit = ei.t; // cached time for entry ei
      int m = n>>>1; // number of entries with at least one child
      while (i<m) { // while not childless, ...
        int c = (i<<1)+1; // index of left child
        int r = c+1; // index of right child
        Entry ec = e[c]; // initially assume left child smallest
        if (r<n && e[r].t<ec.t) // but if right child smallest, ...
          ec = e[c=r]; // the smaller of left and right children
        if (eit<=ec.t) // break if ei not greater than smaller child
          break;
        set(i,ec); // ei greater than smaller child, so move smaller child up
        i = c;
      }
      if (ei!=e[i]) // if necessary, ...
        set(i,ei); // set ei where it belongs
    }

    // If necessary, moves entry e[i] up so not less than parent.
    private void siftUp(int i) {
      Entry ei = e[i]; // entry ei that may move up
      float eit = ei.t; // cached time for entry ei
      while (i>0) { // while a parent (not the root entry), ...
        int p = (i-1)>>>1; // index of parent
        Entry ep = e[p]; // the parent
        if (eit>=ep.t) // break if ei not less than parent
          break;
        set(i,ep); // ei less than parent, so move parent down
        i = p;
      }
      if (ei!=e[i]) // if necessary, ...
        set(i,ei); // set ei where it belongs
    }

    // Sets the i'th entry and updates the time map.
    private void set(int i, Entry ei) {
      e[i] = ei;
      tmap.setMinTimeHeapIndex(ei,i);
    }

    // Grows this heap to have at least the specified capacity.
    private void grow(int minCapacity) {
      if (minCapacity<0) // overflow
        throw new OutOfMemoryError();
      int oldCapacity = e.length;
      int newCapacity = oldCapacity*2;
      if (newCapacity<0) // overflow
        newCapacity = Integer.MAX_VALUE;
      if (newCapacity<minCapacity)
        newCapacity = minCapacity;
      Entry[] enew = new Entry[newCapacity];
      System.arraycopy(e,0,enew,0,oldCapacity);
      e = enew;
    }
  }
}
