/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

/**
 * A vector represented by a 1D array of doubles.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public class VecArrayDouble1 implements Vec {

  /**
   * Constructs a zero vector with specified dimensions.
   * @param n1 the number of doubles in the 1st dimension.
   */
  public VecArrayDouble1(int n1) {
    _a = new double[n1];
    _n1 = n1;
  }

  /**
   * Constructs a vector that wraps the specified array of doubles.
   * @param a the array of doubles; by reference, not by copy.
   */
  public VecArrayDouble1(double[] a) {
    _a = a;
    _n1 = a.length;
  }

  /**
   * Gets the array of doubles wrapped by this vector.
   * @return the array of doubles; by reference, not by copy.
   */
  public double[] getArray() {
    return _a;
  }

  /**
   * Gets the number of doubles in the 1st array dimension.
   * @return the number of doubles in the 1st dimension.
   */
  public int getN1() {
    return _n1;
  }

  public double epsilon() {
    return Math.ulp(1.0d);
  }

  public VecArrayDouble1 clone() {
    VecArrayDouble1 v = new VecArrayDouble1(_n1);
    System.arraycopy(_a,0,v._a,0,_n1);
    return v;
  }

  public double dot(Vec vthat) {
    double[] athis = _a;
    double[] athat = ((VecArrayDouble1)vthat)._a;
    double sum = 0.0;
    for (int i1=0; i1<_n1; ++i1)
      sum += athis[i1]*athat[i1];
    return sum;
  }

  public double norm2() {
    double sum = 0.0;
    for (int i1=0; i1<_n1; ++i1) {
      double ai = _a[i1];
      sum += ai*ai;
    }
    return Math.sqrt(sum);
  }

  public void zero() {
    for (int i1=0; i1<_n1; ++i1)
      _a[i1] = 0.0f;
  }

  public void scale(double s) {
    for (int i1=0; i1<_n1; ++i1)
      _a[i1] *= s;
  }

  public void add(double sthis, Vec vthat, double sthat) {
    double[] athis = _a;
    double[] athat = ((VecArrayDouble1)vthat)._a;
    for (int i1=0; i1<_n1; ++i1)
      athis[i1] = athis[i1]*sthis+athat[i1]*sthat;
  }

  private double[] _a;
  private int _n1;
}
