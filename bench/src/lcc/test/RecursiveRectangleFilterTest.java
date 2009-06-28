/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc.test;

import junit.framework.*;
import edu.mines.jtk.dsp.Conv;
import edu.mines.jtk.util.ArrayMath;
import static edu.mines.jtk.util.MathPlus.*;

import lcc.RecursiveRectangleFilter;

/**
 * Tests {@link lcc.RecursiveRectangleFilter}.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.08.13
 */
public class RecursiveRectangleFilterTest extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(RecursiveRectangleFilterTest.class);
    junit.textui.TestRunner.run(suite);
  }

  public void test1() {
    int[] ns = {4,11};
    int[] ls = {-11,-4,4,11};
    int[] ms = {-11,-4,4,11};
    for (int i=0; i<ns.length; ++i) {
      int n = ns[i];
      for (int j=0; j<ls.length; ++j) {
        int l = ls[j];
        for (int k=0; k<ms.length; ++k) {
          int m = ms[k];
          if (m<l) 
            continue;
          RFilter rf = new RFilter(l,m);
          SFilter sf = new SFilter(l,m);
          float[] x = ArrayMath.randfloat(n);
          float[] ys = ArrayMath.randfloat(n);
          float[] yr = ArrayMath.randfloat(n);
          sf.apply(x,ys);
          rf.apply(x,yr);
          assertEqual(ys,yr);
        }
      }
    }
  }

  public void test2() {
    int[] ns = {4,11};
    int[] ls = {-11,-4,4,11};
    int[] ms = {-11,-4,4,11};
    for (int i2=0; i2<ns.length; ++i2) {
      for (int i1=0; i1<ns.length; ++i1) {
        int n1 = ns[i1];
        int n2 = ns[i2];
        for (int j=0; j<ls.length; ++j) {
          int l = ls[j];
          for (int k=0; k<ms.length; ++k) {
            int m = ms[k];
            if (m<l) 
              continue;
            RFilter rf = new RFilter(l,m);
            SFilter sf = new SFilter(l,m);
            float[][] x = ArrayMath.randfloat(n1,n2);
            float[][] ys = ArrayMath.randfloat(n1,n2);
            float[][] yr = ArrayMath.randfloat(n1,n2);
            sf.apply1(x,ys);
            rf.apply1(x,yr);
            assertEqual(ys,yr);
            sf.apply2(x,ys);
            rf.apply2(x,yr);
            assertEqual(ys,yr);
          }
        }
      }
    }
  }

  public void test3() {
    int[] ns = {4,11};
    int[] ls = {-11,-4,4,11};
    int[] ms = {-11,-4,4,11};
    for (int i3=0; i3<ns.length; ++i3) {
      for (int i2=0; i2<ns.length; ++i2) {
        for (int i1=0; i1<ns.length; ++i1) {
          int n1 = ns[i1];
          int n2 = ns[i2];
          int n3 = ns[i3];
          for (int j=0; j<ls.length; ++j) {
            int l = ls[j];
            for (int k=0; k<ms.length; ++k) {
              int m = ms[k];
              if (m<l) 
                continue;
              RFilter rf = new RFilter(l,m);
              SFilter sf = new SFilter(l,m);
              float[][][] x = ArrayMath.randfloat(n1,n2,n3);
              float[][][] ys = ArrayMath.randfloat(n1,n2,n3);
              float[][][] yr = ArrayMath.randfloat(n1,n2,n3);
              sf.apply1(x,ys);
              rf.apply1(x,yr);
              assertEqual(ys,yr);
              sf.apply2(x,ys);
              rf.apply2(x,yr);
              assertEqual(ys,yr);
              sf.apply3(x,ys);
              rf.apply3(x,yr);
              assertEqual(ys,yr);
            }
          }
        }
      }
    }
  }

  private void assertEqual(float[] re, float[] ra) {
    int n = re.length;
    float tolerance = (float)(n)*FLT_EPSILON;
    for (int i=0; i<n; ++i) {
      if (abs(re[i]-ra[i])>tolerance)
        System.out.println("i="+i);
      assertEquals(re[i],ra[i],tolerance);
    }
  }

  private void assertEqual(float[][] re, float[][] ra) {
    int n2 = re.length;
    int n1 = re[0].length;
    float tolerance = (float)(n1*n2)*FLT_EPSILON;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        assertEquals(re[i2][i1],ra[i2][i1],tolerance);
  }

  private void assertEqual(float[][][] re, float[][][] ra) {
    int n3 = re.length;
    int n2 = re[0].length;
    int n1 = re[0][0].length;
    float tolerance = (float)(n1*n2*n3)*FLT_EPSILON;
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          assertEquals(re[i3][i2][i1],ra[i3][i2][i1],tolerance);
  }

  private interface Filter {
    public void apply(float[] x, float[] y);
    public void apply1(float[][] x, float[][] y);
    public void apply2(float[][] x, float[][] y);
    public void apply1(float[][][] x, float[][][] y);
    public void apply2(float[][][] x, float[][][] y);
    public void apply3(float[][][] x, float[][][] y);
  }

  // Recursive filter.
  private class RFilter implements Filter {
    public RFilter(int l, int m) {
      _rrf = new RecursiveRectangleFilter(l,m);
    }
    public void apply(float[] x, float[] y) {
      _rrf.apply(x,y);
    }
    public void apply1(float[][] x, float[][] y) {
      _rrf.apply1(x,y);
    }
    public void apply2(float[][] x, float[][] y) {
      _rrf.apply2(x,y);
    }
    public void apply1(float[][][] x, float[][][] y) {
      _rrf.apply1(x,y);
    }
    public void apply2(float[][][] x, float[][][] y) {
      _rrf.apply2(x,y);
    }
    public void apply3(float[][][] x, float[][][] y) {
      _rrf.apply3(x,y);
    }
    private RecursiveRectangleFilter _rrf;
  }

  // Slow equivalent.
  private class SFilter implements Filter {
    public SFilter(int l, int m) {
      _l = l;
      _n = 1+m-l;
      float s = 1.0f/(float)(_n);
      _f = ArrayMath.fillfloat(s,_n);
    }
    public void apply(float[] x, float[] y) {
      int n = x.length;
      Conv.xcor(_n,_l,_f,n,0,x,n,0,y);
    }
    public void apply1(float[][] x, float[][] y) {
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] f = ArrayMath.reshape(_n,1,_f);
      Conv.xcor(_n,1,_l,0,f,n1,n2,0,0,x,n1,n2,0,0,y);
    }
    public void apply2(float[][] x, float[][] y) {
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] f = ArrayMath.reshape(1,_n,_f);
      Conv.xcor(1,_n,0,_l,f,n1,n2,0,0,x,n1,n2,0,0,y);
    }
    public void apply1(float[][][] x, float[][][] y) {
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      float[][][] f = ArrayMath.reshape(_n,1,1,_f);
      Conv.xcor(_n,1,1,_l,0,0,f,n1,n2,n3,0,0,0,x,n1,n2,n3,0,0,0,y);
    }
    public void apply2(float[][][] x, float[][][] y) {
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      float[][][] f = ArrayMath.reshape(1,_n,1,_f);
      Conv.xcor(1,_n,1,0,_l,0,f,n1,n2,n3,0,0,0,x,n1,n2,n3,0,0,0,y);
    }
    public void apply3(float[][][] x, float[][][] y) {
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      float[][][] f = ArrayMath.reshape(1,1,_n,_f);
      Conv.xcor(1,1,_n,0,0,_l,f,n1,n2,n3,0,0,0,x,n1,n2,n3,0,0,0,y);
    }
    private int _l,_n;
    private float[] _f;
  }
}
