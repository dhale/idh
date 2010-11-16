package lss;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Facilitates loops for 2D filter stencils. In effect, this class pads 
 * 2-D arrays with extra values on all sides, without making a complete 
 * copy of the entire array. This buffering and padding simplifies 
 * common loops over array values when filtering 2D arrays.
 */
public class StencilBuffer2 {

  /**
   * The method used to extrapolate values beyond the ends of input arrays.
   * The default is extrapolation with zero values.
   */
  public enum Extrapolation {
    /**
     * Extrapolate with zero values.
     */
    ZERO_VALUE,
    /**
     * Extrapolate values at the ends with zero slope.
     */
    ZERO_SLOPE
  };

  /**
   * Constructs a buffer for specified parameters. The buffer will have 
   * storage for l2+1+m2 arrays, each with l1+n1+m1 values.
   * @param l1 number of extra samples at beginning in 1st dimension.
   * @param m1 number of extra samples at end in 1st dimension.
   * @param n1 number of samples (not counting extras) in 1st dimension.
   * @param l2 number of extra samples at beginning in 2nd dimension.
   * @param m2 number of extra samples at end in 2nd dimension.
   * @param n2 number of samples (not counting extras) in 2nd dimension.
   */
  StencilBuffer2(int l1, int m1, int n1, int l2, int m2, int n2) {
    _l1 = l1;
    _m1 = m1;
    _n1 = n1;
    _l2 = l2;
    _m2 = m2;
    _n2 = n2;
    _nb2 = _l2+1+_m2;
    _nb1 = _l1+_n1+_m1;
    _i = fillint(-_nb2,_nb2);
    _b = new float[_nb2][_nb1];
  }

  /**
   * Gets buffered values from the specified array.
   * Extrapolates values if the specified index is out of bounds.
   * <p>
   * The returned buffered array has l1+m1 extra values at the ends;
   * the first l1 values and the last m1 values are extrapolated.
   * @param i2 index in 2nd dimension of the buffered array to get.
   * @param a the array from which to copy values to this buffer.
   * @return reference to a buffered array of values.
   */
  public float[] get(int i2, float[][] a) {
    int j2 = j2(i2);
    Check.argument(-_l2<=i2 && i2<=_n2+_m2-1,"index i2="+i2+" is in bounds");
    if (_i[j2]!=i2) {
      if (_extrapolation==Extrapolation.ZERO_SLOPE) {
        i2 = max(0,min(_n2-1,i2));
        copy(_n1,0,a[i2],_l1,_b[j2]);
        fill(a[i2][0],_l1,0,_b[j2]);
        fill(a[i2][_n1-1],_m1,_l1+_n1,_b[j2]);
      } else if (_extrapolation==Extrapolation.ZERO_VALUE) {
        if (0<=i2 && i2<_n2) {
          copy(_n1,0,a[i2],_l1,_b[j2]);
          fill(0.0f,_l1,0,_b[j2]);
          fill(0.0f,_m1,_l1+_n1,_b[j2]);
        } else {
          zero(_b[j2]);
        }
      }
      _i[j2] = i2;
    }
    return _b[j2];
  }

  /**
   * Copies buffered values into the specified array.
   * Values for the specified index must be in this buffer.
   * @param i2 index in 2nd dimension of the array to set.
   * @param a the array into which to copy values from this buffer.
   */
  public void set(int i2, float[][] a) {
    int j2 = j2(i2);
    Check.argument(0<=i2 && i2<=_n2,"index i2="+i2+" is in bounds");
    Check.state(_i[j2]==i2,"array with index i2="+i2+" is in buffer");
    float[] bj = _b[j2];
    float[] ai = a[i2];
    for (int i1=0,j1=_l1; i1<_n1; ++i1,++j1)
      ai[i1] = bj[j1]; 
  }

  /**
   * Adds buffered values into the specified array.
   * Values for the specified index must be in this buffer.
   * @param i2 index in 2nd dimension of the array to which to add.
   * @param a the array into which to add values from this buffer.
   */
  public void add(int i2, float[][] a) {
    int j2 = j2(i2);
    Check.argument(0<=i2 && i2<=_n2,"index i2="+i2+" is in bounds");
    Check.state(_i[j2]==i2,"array with index i2="+i2+" is in buffer");
    float[] bj = _b[j2];
    float[] ai = a[i2];
    for (int i1=0,j1=_l1; i1<_n1; ++i1,++j1)
      ai[i1] += bj[j1]; 
  }

  ////////////////////////////////////////////////////////////////////////////
  // private

  private int _l1,_l2;
  private int _m1,_m2;
  private int _n1,_n2;
  private int _nb1,_nb2;
  private int[] _i;
  private float[][] _b;
  private Extrapolation _extrapolation;

  private int j2(int i2) {
    return (i2+_l2)%_nb1;
  }

  private static void fill(float a, int n, int j, float[] b) {
    while (n>0)
      b[j++] = a;
  }
}
