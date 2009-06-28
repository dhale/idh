/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.ArrayMath;

/**
 * A min- or max-heap of times sampled in a 1D array. 
 * Such a heap is typically used in fast marching methods. It enhances
 * a conventional heap by maintaining a map of heap entries in a 1D array 
 * of indices. Given an array index i, this index map enables O(1) 
 * access to heap entries. Such fast access is required as times in the
 * heap are reduced while marching.
 * <p>
 * Depending on the type of heap, min or max, the entry with either the 
 * smallest or largest time is stored at the top (root) of the heap. This
 * entry can be accessed with complexity O(1) and removed with complexity 
 * O(log N), where N is the number of entries in the heap. Complexity for
 * inserting new entries or reducing the times for existing entries is
 * O(log N).
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.09.30
 */
public class TimeHeap1 {

  /**
   * The heap type determines whether the entry at the top has the minimum
   * or maximum time. That top entry can be accessed with O(1) complexity
   * and removed in O(log N) complexity, where N is the number times in the
   * heap.
   */
  public enum Type {MIN,MAX};

  /**
   * An entry in the heap has sample indices i1 and a time t.
   */
  public static class Entry {
    public int i1;
    public float t;
  }

  /**
   * Constructs a heap with specified type and array dimensions.
   * @param type the type of heap.
   * @param n1 number of samples.
   */
  public TimeHeap1(Type type, int n1) {
    _type = type;
    _n1 = n1;
    _imap = new int[n1];
  }

  /**
   * Gets the type of this heap.
   * @return the type.
   */
  public Type getType() {
    return _type;
  }

  /**
   * Gets the number of samples.
   * @return the number of samples.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Inserts a new entry into this heap with specified index and time.
   * The heap must not already contain an entry with that index.
   * @param i1 the sample index.
   * @param t the time.
   */
  public void insert(int i1, float t) {
    int i = indexOf(i1); // index of entry with time to reduce
    Check.argument(i<0,"entry with index i1 is not in the heap");
    i = _n; // index at which to insert the entry
    if (_n==_e.length) // if necessary, ...
      grow(_n+1); // increase the capacity of this heap
    Entry ei = _e[i];
    if (ei==null) // if an unused entry does not already exist, ...
      ei = new Entry(); // construct a new entry
    ei.i1 = i1;
    ei.t = t;
    set(i,ei);
    siftUp(i);
    ++_n;
  }

  /**
   * Reduces the time of the entry in this heap with specified index.
   * This heap must already contain an entry with that index, and
   * the specified time must be less than the time for that entry.
   * @param i1 the sample index in 1st dimension.
   * @param t the reduced time.
   */
  public void reduce(int i1, float t) {
    int i = indexOf(i1); // index of entry with time to reduce
    Check.argument(i>=0,"entry with index i1 is in the heap");
    Check.argument(t<_e[i].t,"specified time less than time in heap");
    _e[i].t = t; // reduce the time
    if (_type==Type.MIN) { // for a min-heap, ...
      siftUp(i); // the entry may need to move up
    } else { // but for a max heap, ...
      siftDown(i); // the entry may need to move down
    }
  }

  /**
   * Removes and returns the heap entry with smallest/largest time.
   * The heap must not be empty.
   */
  public Entry remove() {
    Check.state(_n>0,"heap is not empty");
    Entry e0 = _e[0];
    --_n;
    if (_n>0) {
      set(0,_e[_n]);
      set(_n,e0);
      siftDown(0);
    }
    return e0;
  }

  /**
   * Removes all entries from this heap.
   */
  public void clear() {
    _n = 0;
  }

  /**
   * Returns the number of entries in this heap.
   */
  public int size() {
    return _n;
  }

  /**
   * Returns true if this heap is empty; false, otherwise.
   */
  public boolean isEmpty() {
    return _n==0;
  }

  /**
   * Dumps this heap to standard output; leading spaces show level in tree.
   */
  public void dump() {
    dump("",0);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Type _type; // heap type, either MIN or MAX
  private int _n1; // array dimension
  private int _n; // number of entries in this heap
  private int[] _imap; // maps array index i1 to heap index i
  private Entry[] _e = new Entry[1024]; // array of entries in this heap

  /**
   * Returns the heap index of the entry with array indices (i1,i2).
   * If such an entry is not in the heap, this method returns -1, but
   * throws an exception if the indices i1 or i2 are out of bounds.
   */
  private int indexOf(int i1) {
    int i = _imap[i1];
    if (i<0 || i>=_n) 
      return -1;
    Entry ei = _e[i];
    if (ei.i1!=i1)
      return -1;
    return i;
  }

  /**
   * Sets the i'th entry, and updates the index map accordingly.
   */
  private void set(int i, Entry ei) {
    _e[i] = ei;
    _imap[ei.i1] = i;
  }

  /**
   * If necessary, moves entry e[i] down so not greater/less than children.
   */
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

  /**
   * If necessary, moves entry e[i] up so not less/greater than parent.
   */
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

  /**
   * Grows this heap to have at least the specified capacity.
   */
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

  /**
   * Recursively dumps heap entries with leading spaces.
   */
  private void dump(String s, int i) {
    if (i<_n) {
      s = s+"  ";
      Entry e = _e[i];
      System.out.println(s+e.i1+" "+e.t);
      dump(s,2*i+1);
      dump(s,2*i+2);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    testHeap(new TimeHeap1(TimeHeap1.Type.MIN,9));
    testHeap(new TimeHeap1(TimeHeap1.Type.MAX,9));
  }

  private static void testHeap(TimeHeap1 heap) {
    int n1 = heap.getN1();
    int n = n1;
    float[] s = ArrayMath.randfloat(n);
    float[] t = ArrayMath.copy(s);
    for (int i1=0; i1<n1; ++i1) {
      float ti = t[i1];
      heap.insert(i1,ti);
      s[i1] = ti;
    }
    for (int i1=0; i1<n1; ++i1) {
      s[i1] -= 0.5f;
      t[i1] -= 0.5f;
      heap.reduce(i1,t[i1]);
    }
    assert !heap.isEmpty();
    assert heap.size()==n;
    ArrayMath.quickSort(s); // increasing order
    if (heap.getType()==TimeHeap1.Type.MAX)
      s = ArrayMath.reverse(s); // decreasing order
    for (int i=0; i<n; ++i) {
      Entry e = heap.remove();
      float ti = e.t;
      assert ti==s[i];
    }
    assert heap.isEmpty();
    assert heap.size()==0;
  }
}
