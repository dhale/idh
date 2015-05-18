/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

/**
 * A vector represented by a 2D array[n2][n1] of doubles.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public class VecArrayDouble2 implements Vec {

  /**
   * Constructs a zero vector with specified dimensions.
   * @param n1 the number of doubles in the 1st dimension.
   * @param n2 the number of doubles in the 2nd dimension.
   */
  public VecArrayDouble2(int n1, int n2) {
    _a = new double[n2][n1];
    _n1 = n1;
    _n2 = n2;
  }

  /**
   * Constructs a vector that wraps the specified array of doubles.
   * @param a the array of doubles; by reference, not by copy.
   */
  public VecArrayDouble2(double[][] a) {
    _a = a;
    _n1 = a[0].length;
    _n2 = a.length;
  }

  /**
   * Gets the array of doubles wrapped by this vector.
   * @return the array of doubles; by reference, not by copy.
   */
  public double[][] getArray() {
    return _a;
  }

  /**
   * Gets the number of doubles in the 1st array dimension.
   * @return the number of doubles in the 1st dimension.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of doubles in the 2nd array dimension.
   * @return the number of doubles in the 2nd dimension.
   */
  public int getN2() {
    return _n2;
  }

  public double epsilon() {
    return Math.ulp(1.0f);
  }

  public VecArrayDouble2 clone() {
    VecArrayDouble2 v = new VecArrayDouble2(_n1,_n2);
    for (int i2=0; i2<_n2; ++i2)
      System.arraycopy(_a[i2],0,v._a[i2],0,_n1);
    return v;
  }

  public double dot(Vec vthat) {
    double[][] athis = _a;
    double[][] athat = ((VecArrayDouble2)vthat)._a;
    double sum = 0.0;
    for (int i2=0; i2<_n2; ++i2) {
      double[] athis2 = athis[i2];
      double[] athat2 = athat[i2];
      for (int i1=0; i1<_n1; ++i1)
        sum += athis2[i1]*athat2[i1];
    }
    return sum;
  }

  public double norm2() {
    double sum = 0.0;
    for (int i2=0; i2<_n2; ++i2) {
      double[] a2 = _a[i2];
      for (int i1=0; i1<_n1; ++i1) {
        double ai = a2[i1];
        sum += ai*ai;
      }
    }
    return Math.sqrt(sum);
  }

  public void zero() {
    for (int i2=0; i2<_n2; ++i2) {
      double[] a2 = _a[i2];
      for (int i1=0; i1<_n1; ++i1) {
        a2[i1] = 0.0f;
      }
    }
  }

  public void scale(double s) {
    for (int i2=0; i2<_n2; ++i2) {
      double[] a2 = _a[i2];
      for (int i1=0; i1<_n1; ++i1) {
        a2[i1] *= s;
      }
    }
  }

  public void add(double sthis, Vec vthat, double sthat) {
    double[][] athis = _a;
    double[][] athat = ((VecArrayDouble2)vthat)._a;
    for (int i2=0; i2<_n2; ++i2) {
      double[] athis2 = athis[i2];
      double[] athat2 = athat[i2];
      for (int i1=0; i1<_n1; ++i1) {
        athis2[i1] = athis2[i1]*sthis+athat2[i1]*sthat;
      }
    }
  }

  private double[][] _a;
  private int _n1,_n2;
}
