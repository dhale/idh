package vs30;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Utilities for Vs30.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.01.08
 */
public class Util {

  /**
   * Converts signed short (16-bit) integers to floats.
   * Replaces any non-positive (null) values with zeros.
   * @param x array of shorts.
   * @return array of floats.
   */
  public static float[][] floatsFromShorts(short[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        short xi = x[i2][i1];
        y[i2][i1] = (xi>0)?xi:0.0f;
      }
    }
    return y;
  }
}
