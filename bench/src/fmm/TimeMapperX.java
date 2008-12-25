/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * A time transform for a 2D anisotropic eikonal equation. The non-linear 
 * equation is grad(t) dot W*grad(t) = 1, where t is the solution time 
 * map and W denotes a positive-definite (velocity-squared) metric tensor 
 * field.
 * <p>
 * Times in a time map are flagged as either known or unknown. Times for 
 * known samples are typically zero (or near zero), and are never modified. 
 * Times for unknown samples are computed to be solutions to the eikonal 
 * equation. Each of these computed times represents the minimum time
 * to travel from the unknown sample to one of the known samples.
 * <p>
 * A separate map of marks may be computed as times are computed. The 
 * mark of each unknown sample equals the mark of that known sample for 
 * which the time between the two samples is minimized. For an isotropic 
 * homogeneous tensor field, time equals Euclidean distance, and a map of 
 * marks that represents a discrete Voronoi diagram can be computed by 
 * assigning each known sample a unique mark.
 * <p>
 * This transform uses an iterative sweeping method to compute the time map.
 * Iterations are similar to those described by Tsai et al., 2002.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.12.07
 */
public class TimeMapperX {
  
  /**
   * Constructs a time mapper for the specified tensor field.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param tensors velocity-squared tensors.
   */
  public TimeMapperX(int n1, int n2, Tensors2 tensors) {
    _n1 = n1;
    _n2 = n2;
    _tensors = tensors;
  }

  /**
   * Sets the tensors used by this time mapper.
   * @param tensors the tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = tensors;
  }

  /**
   * Applys this transform to the specified array of times and marks.
   * @param known array of flags; zero for unknown samples
   * @param times array of times in which unknown times are to be computed.
   * @param marks array of marks in which unknown marks are to be computed.
   */
  public void apply(byte[][] known, float[][] times, int[][] marks) {
    init(known,times);
    int nsweep = 0;
    for (boolean again=true; nsweep<1000 && again; ++nsweep) {
      again = false;
      for (int ksweep=0; ksweep<NSWEEP; ++ksweep) {
        int k1 = K1SWEEP[ksweep];
        int k2 = K2SWEEP[ksweep];
        again = sweep(k1,k2,known,times) || again;
      }
    }
    trace("TimeMapperX.apply: nsweep="+nsweep);
    //perturb(times);
    domarks(known,times,marks);
  }
  private static void perturb(float[][] t) {
    Random random = new Random();
    int n1 = t[0].length;
    int n2 = t.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        t[i2][i1] *= 1.0f+0.01f*random.nextFloat();
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Default time for samples not yet computed.
  private static final float INFINITY = Float.MAX_VALUE;

  // Times are converged when the fractional change is less than this value.
  private static final float EPSILON = 0.001f;
  private static final float ONE_MINUS_EPSILON = 1.0f-EPSILON;

  // Indices (k1,k2) used in sweeps to compute times.
  private static final int[] K1SWEEP = {-1, 1,-1, 1};
  private static final int[] K2SWEEP = {-1,-1, 1, 1};
  private static final int NSWEEP = K1SWEEP.length;

  // Indices (k1,k2) of neighbors used to compute marks.
  //private static final int[] K1NABOR = {-1, 1, 0, 0,-1, 1,-1, 1};
  //private static final int[] K2NABOR = { 0, 0,-1, 1,-1,-1, 1, 1};
  private static final int[] K1NABOR = {-1, 1, 0, 0};
  private static final int[] K2NABOR = { 0, 0,-1, 1};
  private static final int NNABOR = K1NABOR.length;

  // More efficient than ArrayStack<Integer>.
  private static class IntStack {
    void push(int k) {
      if (_n==_a.length) {
        int[] a = new int[2*_n];
        for (int i=0; i<_n; ++i)
          a[i] = _a[i];
        _a = a;
      }
      _a[_n++] = k;
    }
    int pop() {
      return _a[--_n];
    }
    int size() {
      return _n;
    }
    void clear() {
      _n = 0;
    }
    boolean isEmpty() {
      return _n==0;
    }
    private int _n = 0;
    private int[] _a = new int[2048];
  }

  private int _n1,_n2;
  private Tensors2 _tensors;

  /**
   * Initializes times and returns corresponding array for directions.
   */
  private void init(byte[][] known, float[][] times) {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (known[i2][i1]==0)
          times[i2][i1] = INFINITY;
      }
    }
  }

  /**
   * Performs one sweep of updates of specified times.
   * The indices k1 and k2 are sample offsets that determine which
   * samples are used to compute the update, as well as the sweep
   * direction. Valid (k1,k2) are (-1,-1), (-1,1), (1,-1), and (1,1).
   * Returns true iff change in time for at least one sample is significant.
   */ 
  private boolean sweep(int k1, int k2, byte[][] known, float[][] times) {
    boolean again = false;
    int i1b = (k1<0)? 0:_n1-1;
    int i1e = (k1>0)?-1:_n1;
    int i2b = (k2<0)? 0:_n2-1;
    int i2e = (k2>0)?-1:_n2;
    float[] d = new float[3];
    for (int i2=i2b; i2!=i2e; i2-=k2) {
      int j2 = i2+k2;
      boolean b2 = 0<=j2 && j2<_n2;
      for (int i1=i1b; i1!=i1e; i1-=k1) {
        int j1 = i1+k1;
        boolean b1 = 0<=j1 && j1<_n1;
        if (known[i2][i1]!=0)
          continue;
        _tensors.getTensor(i1,i2,d);
        float d11 = d[0];
        float d12 = d[1];
        float d22 = d[2];
        float e12 = 1.0f/(d11*d22-d12*d12);
        float t0 = times[i2][i1];
        float t1 = b1?times[i2][j1]:INFINITY;
        float t2 = b2?times[j2][i1]:INFINITY;
        float tc = computeTime(d11,d12,d22,k1,k2,t1,t2);
        if (tc<t0) {
          again = again || tc<t0*ONE_MINUS_EPSILON;
          times[i2][i1] = tc;
        }
      }
    }
    return again;
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
    float t0 = INFINITY;
    double ds11 = d11*s1*s1;
    double ds12 = d12*s1*s2;
    double ds22 = d22*s2*s2;
    float e12 = 1.0f/(d11*d22-d12*d12);
    if (t1<INFINITY && t2<INFINITY) {
      double t12 = t1-t2;
      double a = ds11+2.0*ds12+ds22;
      double b = 2.0*(ds12+ds22)*t12;
      double c = ds22*t12*t12-1.0;
      double d = b*b-4.0*a*c;
      if (d>=0.0) {
        double u1 = (-b+sqrt(d))/(2.0*a);
        double u2 = u1+t12;
        double v1 = ds11*u1+ds12*u2;
        double v2 = ds12*u1+ds22*u2;
        if (v1>=0.0 && v2>=0.0)
          t0 = t1+(float)u1;
      }
    }
    if (t1<INFINITY) {
      float tc = t1+sqrt(d22*e12);
      if (tc<t0) t0 = tc;
    }
    if (t2<INFINITY) {
      float tc = t2+sqrt(d11*e12);
      if (tc<t0) t0 = tc;
    }
    return t0;
  }

  private void domarks(byte[][] known, float[][] times, int[][] marks) {
    float[] d = new float[3];
    IntStack stack = new IntStack();
    for (int i2=0,i=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1,++i) {
        if (known[i2][i1]!=0) {
          markFrom(i1,i2,stack,d,known,times,marks);
          plot(marks);
        }
      }
    }
  }
  private void markFrom(
    int i1, int i2, IntStack stack, float[] d, 
    byte[][] known, float[][] times, int[][] marks) 
  {
    // Push the specified sample onto stack.
    int marki = marks[i2][i1];
    stack.clear();
    stack.push(i2);
    stack.push(i1);

    // While stack is not empty, ...
    while (!stack.isEmpty()) {

      // Remove sample from top of stack and get its correct time.
      i1 = stack.pop();
      i2 = stack.pop();
      float tc = times[i2][i1];

      // Has the correct time already been found?
      boolean found = known[i2][i1]!=0;

      // Look for the correct time using only marked neighbors.
      for (int ksweep=0; !found && ksweep<NSWEEP; ++ksweep) {
        int k1 = K1SWEEP[ksweep];
        int k2 = K2SWEEP[ksweep];
        int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
        int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
        boolean m1 = marki==marks[i2][j1];
        boolean m2 = marki==marks[j2][i1];
        _tensors.getTensor(i1,i2,d);
        float d11 = d[0];
        float d12 = d[1];
        float d22 = d[2];
        float t1 = m1?times[i2][j1]:INFINITY;
        float t2 = m2?times[j2][i1]:INFINITY;
        float t0 = computeTime(d11,d12,d22,k1,k2,t1,t2);
        found = t0<=tc;
      }

      // If the correct time was found, ...
      if (found) {

        // Mark this sample and push any unmarked nabors on stack.
        marks[i2][i1] = marki;
        for (int knabor=0; knabor<NNABOR; ++knabor) {
          int k1 = K1NABOR[knabor];
          int k2 = K2NABOR[knabor];
          int j1 = i1+k1;  if (j1<0 || j1>=_n1) continue;
          int j2 = i2+k2;  if (j2<0 || j2>=_n2) continue;
          if (marki!=marks[j2][j1]) {
            stack.push(j2);
            stack.push(j1);
          }
        }
      }
    }
  }

  private void domarksX(byte[][] known, float[][] times, int[][] marks) {
    int n1 = _n1;
    int n2 = _n2;
    float[] d = new float[3];

    // Sort sample indices in order of increasing time.
    int n = n1*n2;
    float[] tsort = new float[n];
    int[] isort = new int[n];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1,++i) {
        tsort[i] = times[i2][i1];
        isort[i] = i;
      }
    }
    Array.quickIndexSort(tsort,isort);
    tsort = null;

    // Flags for samples.
    byte[][] flags = new byte[n2][n1];
    byte marked = 1;

    // Map of border samples.
    BorderMap bmap = new BorderMap(n1,n2);

    // For all samples, in order of increasing time, ...
    for (int i=0; i<n; ++i) {
      int ii = isort[i];
      int i1 = ii%n1;
      int i2 = ii/n1;

      // Skip known samples already marked. We assume that times for all
      // known samples are less than the minimum time for unknown samples.
      if (known[i2][i1]!=0) {
        flags[i2][i1] = marked;
        continue;
      }

      // Determine if marks for all marked nabors are the same. If they are, 
      // then we simply assign that same mark to the current sample.
      boolean marksSame = true;
      boolean markFound = false;
      int markFirst = -1;
      for (int knabor=0; knabor<NNABOR; ++knabor) {
        int k1 = K1NABOR[knabor];
        int k2 = K2NABOR[knabor];
        int j1 = i1+k1; if (j1<0 || j1>=n1) continue;
        int j2 = i2+k2; if (j2<0 || j2>=n2) continue;
        if (flags[j2][j1]==marked) {
          if (markFound) {
            marksSame = marksSame && marks[j2][j1]==markFirst;
          } else {
            markFound = true;
            markFirst = marks[j2][j1];
          }
        }
      }

      // If the marks for marked nabors are not the same, then use the
      // characteristic direction to determine which marked nabor wins.
      // The characteristic direction corresponds to the sweep direction
      // in which the time at the current sample is minimized.
      if (!marksSame) {

        // Compute the group-velocity vector (v1,v2) with least time,
        // and remember the corresponding sweep indices (l1,l2).
        _tensors.getTensor(i1,i2,d);
        float d11 = d[0], d12 = d[1], d22 = d[2];
        int l1 = -1;
        int l2 = -1;
        float v1 = 0.0f;
        float v2 = 0.0f;
        float tmin = INFINITY;
        for (int ksweep=0; ksweep<NSWEEP; ++ksweep) {
          int k1 = K1SWEEP[ksweep];
          int k2 = K2SWEEP[ksweep];
          int j1 = i1+k1; if (j1<0 || j1>=n1) continue;
          int j2 = i2+k2; if (j2<0 || j2>=n2) continue;
          float t1 = times[i2][j1];
          float t2 = times[j2][i1];
          float tk = computeTime(d11,d12,d22,k1,k2,t1,t2);
          if (tk<tmin) {
            float s1 = -k1;
            float s2 = -k2;
            float u1 = s1*(times[i2][i1]-times[i2][j1]);
            float u2 = s2*(times[i2][i1]-times[j2][i1]);
            l1 = j1;
            l2 = j2;
            v1 = d11*u1+d12*u2;
            v2 = d12*u1+d22*u2;
            tmin = tk;
          }
        }

        // Propagate marks for border samples.
        Border b1 = bmap.get(l1,i2);
        Border b2 = bmap.get(i1,l2);
        int m1 = marks[i2][l1];
        int m2 = marks[l2][i1];
        float a1 = abs(v1);
        float a2 = abs(v2);
        float w1 = a1/(a1+a2);
        float w2 = a2/(a1+a2);
        Border b = null;
        if (b1!=null && b2!=null) {
          b = b2.merge(b1);
        } else if (b1!=null) {
          b = b1.merge(m2,w2);
        } else if (b2!=null) {
          b = b2.merge(m1,w1);
        } else {
          b = new Border(m1,w1,m2,w2);
        }
        bmap.set(i1,i2,b);
        markFirst = b.chooseMark();
        markFound = true;
        if (i1==300 && 350<i2 && i2<370) {
          trace("i1="+i1+" i2="+i2);
          trace("  l1="+l1+" l2="+l2);
          trace("  m1="+m1+" m2="+m2);
          trace("  v1="+v1+" v2="+v2);
        }
      }

      if (markFound) {
        marks[i2][i1] = markFirst;
        flags[i2][i1] = marked;
      }
    }
  }
  private static class BorderMap {
    Border[][] bm;
    BorderMap(int n1, int n2) {
      bm = new Border[n2][n1];
    }
    Border get(int i1, int i2) {
      return bm[i2][i1];
    }
    void set(int i1, int i2, Border b) {
      bm[i2][i1] = b;
    }
  }
  private static class Border {
    int ma,mb;
    float wa,wb;
    Border(int ma, float wa, int mb, float wb) {
      this.ma = ma;
      this.wa = wa;
      this.mb = mb;
      this.wb = wb;
    }
    Border merge(Border b) {
      Border c = new Border(ma,wa,mb,wb);
      c = c.merge(b.ma,b.wa);
      c = c.merge(b.mb,b.wb);
      return c;
    }
    Border merge(int m, float w) {
      Border c = new Border(ma,wa,mb,wb);
      if (m==ma) {
        c.wa += w;
        c.wb += 1.0f-w;
      } else if (m==mb) {
        c.wa += 1.0f-w;
        c.wb += w;
      } else if (wa>=wb) {
        c.mb = m;
        c.wb = w;
      } else {
        c.ma = m;
        c.wa = w;
      }
      return c;
    }
    int chooseMark() {
      if (wa>=wb) {
        wa -= 1.0f;
        return ma;
      } else {
        wb -= 1.0f;
        return mb;
      }
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static float[][] toFloat(byte[][] b) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        f[i2][i1] = (float)b[i2][i1];
    return f;
  }
  private static float[][] toFloat(int[][] i) {
    int n1 = i[0].length;
    int n2 = i.length;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        f[i2][i1] = (float)i[i2][i1];
    return f;
  }
  private static void plot(byte[][] b) {
    plot(toFloat(b));
  }
  private static void plot(int[][] i) {
    plot(toFloat(i));
  }
  private static void plot(float[][] f) {
    edu.mines.jtk.mosaic.SimplePlot sp =
      new edu.mines.jtk.mosaic.SimplePlot(
        edu.mines.jtk.mosaic.SimplePlot.Origin.UPPER_LEFT);
    sp.setSize(920,900);
    edu.mines.jtk.mosaic.PixelsView pv = sp.addPixels(f);
    pv.setColorModel(edu.mines.jtk.awt.ColorMap.JET);
    pv.setInterpolation(edu.mines.jtk.mosaic.PixelsView.Interpolation.NEAREST);
  }
}
