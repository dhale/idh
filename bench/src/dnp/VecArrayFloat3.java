/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.concurrent.atomic.AtomicInteger;
import edu.mines.jtk.util.*;

/**
 * A vector represented by a 3D array[n3][n2][n1] of floats.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public class VecArrayFloat3 implements Vec {

  /**
   * Constructs a zero vector with specified dimensions.
   * @param n1 the number of floats in the 1st dimension.
   * @param n2 the number of floats in the 2nd dimension.
   * @param n3 the number of floats in the 3rd dimension.
   */
  public VecArrayFloat3(int n1, int n2, int n3) {
    _a = new float[n3][n2][n1];
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
  }

  /**
   * Constructs a vector that wraps the specified array of floats.
   * @param a the array of floats; by reference, not by copy.
   */
  public VecArrayFloat3(float[][][] a) {
    _a = a;
    _n1 = a[0][0].length;
    _n2 = a[0].length;
    _n3 = a.length;
  }

  /**
   * Gets the array of floats wrapped by this vector.
   * @return the array of floats; by reference, not by copy.
   */
  public float[][][] getArray() {
    return _a;
  }

  /**
   * Gets the number of floats in the 1st array dimension.
   * @return the number of floats in the 1st dimension.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of floats in the 2nd array dimension.
   * @return the number of floats in the 2nd dimension.
   */
  public int getN2() {
    return _n2;
  }

  /**
   * Gets the number of floats in the 3rd array dimension.
   * @return the number of floats in the 3rd dimension.
   */
  public int getN3() {
    return _n3;
  }

  public double epsilon() {
    return Math.ulp(1.0f);
  }

  public VecArrayFloat3 clone() {
    VecArrayFloat3 v = new VecArrayFloat3(_n1,_n2,_n3);
    scopy(_a,v._a);
    return v;
  }

  public double dot(Vec vthat) {
    float[][][] athis = _a;
    float[][][] athat = ((VecArrayFloat3)vthat)._a;
    return sdot(athis,athat);
  }

  public double norm2() {
    return Math.sqrt(sdot(_a,_a));
  }

  public void zero() {
    szero(_a);
  }

  public void scale(double s) {
    sscal((float)s,_a);
  }

  public void add(double sthis, Vec vthat, double sthat) {
    float fthis = (float)sthis;
    float fthat = (float)sthat;
    float[][][] athis = _a;
    float[][][] athat = ((VecArrayFloat3)vthat)._a;
    if (fthis==1.0f) {
      saxpy(fthat,athat,athis);
    } else if (fthat==1.0f) {
      sxpay(fthis,athat,athis);
    } else {
      saxpby(fthat,athat,fthis,athis);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _a;
  private int _n1,_n2,_n3;
  private int _nthread = Threads.getAvailableProcessors();
  //private int _nthread = 0;

  // Zeros array x.
  private static void szero(float[] x) {
    ArrayMath.zero(x);
  }
  private static void szero(float[][] x) {
    ArrayMath.zero(x);
  }
  private void szero(float[][][] x) {
    if (_nthread>1) {
      szeroP(x);
    } else {
      szeroS(x);
    }
  }
  private void szeroS(float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      szero(x[i3]);
  }
  private void szeroP(final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            szero(x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Copys array x to array y.
  private void scopy(float[] x, float[] y) {
    ArrayMath.copy(x,y);
  }
  private void scopy(float[][] x, float[][] y) {
    ArrayMath.copy(x,y);
  }
  private void scopy(float[][][] x, float[][][] y) {
    if (_nthread>1) {
      scopyP(x,y);
    } else {
      scopyS(x,y);
    }
  }
  private void scopyS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      scopy(x[i3],y[i3]);
  }
  private void scopyP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            scopy(x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Returns the dot product x'y.
  private double sdot(float[] x, float[] y) {
    int n1 = x.length;
    double d = 0.0;
    for (int i1=0; i1<n1; ++i1)
      d += x[i1]*y[i1];
    return d;
  }
  private double sdot(float[][] x, float[][] y) {
    int n2 = x.length;
    double d = 0.0;
    for (int i2=0; i2<n2; ++i2)
      d += sdot(x[i2],y[i2]);
    return d;
  }
  private double sdot(float[][][] x, float[][][] y) {
    if (_nthread>1) {
      return sdotP(x,y);
    } else {
      return sdotS(x,y);
    }
  }
  private double sdotS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    double d = 0.0;
    for (int i3=0; i3<n3; ++i3)
      d += sdot(x[i3],y[i3]);
    return d;
  }
  private double sdotP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicDouble ad = new AtomicDouble(0.0);
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          double d = 0.0;
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            d += sdot(x[i3],y[i3]);
          ad.getAndAdd(d);
        }
      });
    }
    Threads.startAndJoin(threads);
    return ad.get();
  }

  // Computes x = a*x.
  private void sscal(float a, float[] x) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      x[i1] *= a;
  }
  private void sscal(float a, float[][] x) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      sscal(a,x[i2]);
  }
  private void sscal(float a, float[][][] x) {
    if (_nthread>1) {
      sscalP(a,x);
    } else {
      sscalS(a,x);
    }
  }
  private void sscalS(float a, float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sscal(a,x[i3]);
  }
  private void sscalP(final float a, final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            sscal(a,x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Computes y = y + a*x.
  private void saxpy(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += a*x[i1];
  }
  private void saxpy(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      saxpy(a,x[i2],y[i2]);
  }
  private void saxpy(float a, float[][][] x, float[][][] y) {
    if (_nthread>1) {
      saxpyP(a,x,y);
    } else {
      saxpyS(a,x,y);
    }
  }
  private void saxpyS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpy(a,x[i3],y[i3]);
  }
  private void saxpyP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            saxpy(a,x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Computes y = x + a*y.
  private void sxpay(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = a*y[i1]+x[i1];
  }
  private void sxpay(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      sxpay(a,x[i2],y[i2]);
  }
  private void sxpay(float a, float[][][] x, float[][][] y) {
    if (_nthread>1) {
      sxpayP(a,x,y);
    } else {
      sxpayS(a,x,y);
    }
  }
  private void sxpayS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sxpay(a,x[i3],y[i3]);
  }
  private void sxpayP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            sxpay(a,x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Computes y = a*x + b*y.
  private void saxpby(float a, float[] x, float b, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = a*x[i1]+b*y[i1];
  }
  private void saxpby(float a, float[][] x, float b, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      saxpby(a,x[i2],b,y[i2]);
  }
  private void saxpby(float a, float[][][] x, float b, float[][][] y) {
    if (_nthread>1) {
      saxpbyP(a,x,b,y);
    } else {
      saxpbyS(a,x,b,y);
    }
  }
  private void saxpbyS(float a, float[][][] x, float b, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpby(a,x[i3],b,y[i3]);
  }
  private void saxpbyP(
    final float a, final float[][][] x, final float b, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray(_nthread);
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            saxpby(a,x[i3],b,y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }
}
