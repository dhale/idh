package tp;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Simple griding for well logs.
 * Assumes that wells are predominately aligned with the x1 axis.
 * @author Dave Hale, Colorado School of Mines
 * @version 2010.04.20
 */
public class WellLogGridder {

  public WellLogGridder(Sampling s1, Sampling s2, Sampling s3, float fnull) {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    _g = fillfloat(fnull,n1,n2,n3);
  }

  /**
   * Grids one well log with specified samples.
   * Replaces any gridded values intersected by this log.
   * @param f array of log sample values f(x1,x2,x3).
   * @param x1 array of log sample x1 coordinates.
   * @param x2 array of log sample x2 coordinates.
   * @param x3 array of log sample x3 coordinates.
    */
  public void gridWellLog(float[] f, float[] x1, float[] x2, float[] x3) {
    int n = f.length;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    double d3 = _s3.getDelta();
    double f1 = _s1.getFirst()-0.5*d1;
    double f2 = _s2.getFirst()-0.5*d2;
    double f3 = _s3.getFirst()-0.5*d3;
    double l1 = _s1.getLast()+0.5*d1;
    double l2 = _s2.getLast()+0.5*d2;
    double l3 = _s3.getLast()+0.5*d3;
    
    // Bin 1st dimension of all log samples that lie inside the sampling grid.
    // For each bin, count and sum the sample values and coordinates (x2,x3).
    double[] nsum = new double[n1];
    double[] fsum = new double[n1];
    double[] x2sum = new double[n1];
    double[] x3sum = new double[n1];
    for (int i=0; i<n; ++i) {
      float x1i = x1[i];
      float x2i = x2[i];
      float x3i = x3[i];
      if (f1<=x1i && x1i<=l1 && 
          f2<=x2i && x2i<=l2 && 
          f3<=x3i && x3i<=l3) {
        int i1 = _s1.indexOfNearest(x1i);
        nsum[i1] += 1;
        fsum[i1] += f[i];  
        x2sum[i1] += x2i;
        x3sum[i1] += x3i;
      }
    }

    // Ignore any bins adjacent to empty bins.
    boolean[] ignore = new boolean[n1];
    if (nsum[1]==0) 
      ignore[0] = true;
    for (int i1=1; i1<n1-1; ++i1) {
      if (nsum[i1-1]==0 || nsum[i1+1]==0)
        ignore[i1] = true;
    }
    if (nsum[n1-2]==0) 
      ignore[n1-1] = true;

    // Assign averages to gridded samples.
    for (int i1=0; i1<n1; ++i1) {
      if (!ignore[i1] && nsum[i1]>0) {
        double scale = 1.0/nsum[i1];
        double x2avg = x2sum[i1]*scale;
        double x3avg = x3sum[i1]*scale;
        double favg = fsum[i1]*scale;
        int i2 = _s2.indexOfNearest(x2avg);
        int i3 = _s3.indexOfNearest(x3avg);
        _g[i3][i2][i1] = (float)favg;
      }
    }
  }

  /**
   * Gets the gridded values.
   * @return array of gridded values; by reference, not by copy.
   */
  public float[][][] getGriddedValues() {
    return _g;
  }

  /**
   * Gets the non-null samples from this gridder.
   * @return array {f,x1,x2,x3} of arrays of non-null samples.
   */
  public float[][] getGriddedSamples() {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    int n = 0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (_g[i3][i2][i1]!=_fnull)
            ++n;
        }
      }
    }
    float[] f = new float[n];
    float[] x1 = new float[n];
    float[] x2 = new float[n];
    float[] x3 = new float[n];
    for (int i3=0,i=0; i3<n3; ++i3) {
      float x3i = (float)_s3.getValue(i3);
      for (int i2=0; i2<n2; ++i2) {
        float x2i = (float)_s2.getValue(i2);
        for (int i1=0; i1<n1; ++i1) {
          if (_g[i3][i2][i1]!=_fnull) {
            float x1i = (float)_s1.getValue(i1);
            x1[i] = x1i;
            x2[i] = x2i;
            x3[i] = x3i;
            f[i] = _g[i3][i2][i1];
            ++i;
          }
        }
      }
    }
    return new float[][]{f,x1,x2,x3};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _s1,_s2,_s3;
  private float _fnull;
  private float[][][] _g;
}

