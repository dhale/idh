/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm.test;

import junit.framework.*;

import fmm.TimeHeap2;

/**
 * Tests {@link fmm.TimeHeap2}.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.13
 */
public class TimeHeap2Test extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(TimeHeap2Test.class);
    junit.textui.TestRunner.run(suite);
  }

  public void testMinHeap() {
    testHeap(new TimeHeap2(TimeHeap2.Type.MIN,9,11));
  }

  public void testMaxHeap() {
    testHeap(new TimeHeap2(TimeHeap2.Type.MAX,11,9));
  }

  private static void testHeap(TimeHeap2 heap) {
    TimeHeap2.Type type = heap.getType();
    int n1 = heap.getN1();
    int n2 = heap.getN2();
    int n = n1*n2;
    float[] s = edu.mines.jtk.util.Array.randfloat(n);
    float[][] t = edu.mines.jtk.util.Array.reshape(n1,n2,s);
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        float ti = t[i2][i1];
        heap.insert(i1,i2,ti);
        s[i] = ti;
      }
    }
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        s[i] -= 0.5f;
        t[i2][i1] -= 0.5f;
        heap.reduce(i1,i2,t[i2][i1]);
      }
    }
    assert !heap.isEmpty();
    assert heap.size()==n;
    edu.mines.jtk.util.Array.quickSort(s); // increasing order
    if (type==TimeHeap2.Type.MAX)
      s = edu.mines.jtk.util.Array.reverse(s); // decreasing order
    for (int i=0; i<n; ++i) {
      TimeHeap2.Entry e = heap.remove();
      float ti = e.t;
      assert ti==s[i];
    }
    assert heap.isEmpty();
    assert heap.size()==0;
  }
}
