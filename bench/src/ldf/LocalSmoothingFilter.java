/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic smoothing filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.21
 */
public class LocalSmoothingFilter {

  /**
   * Constructs a local smoothing filter.
   * @param sigma the filter half-width.
   */
  public LocalSmoothingFilter(double sigma) {
    _sigma = (float)sigma;
    makeTable();
  }

  public void applyPass(float[] ds, float[] x, float[] y) {
    int n1 = x.length;

    // Sub-diagonal e of SPD tridiagonal matrix A.
    float[] e = new float[n1];
    if (ds==null) {
      float ss = 0.50f*_sigma*_sigma;
      for (int i1=1; i1<n1; ++i1)
        e[i1] = -ss;
    } else {
      float ss = 0.25f*_sigma*_sigma;
      for (int i1=1; i1<n1; ++i1)
        e[i1] = -ss*(ds[i1-1]+ds[i1]);
    }

    // Diagonal d of SPD tridiagonal matrix.
    float[] d = new float[n1];
    d[0] = 1.0f-e[1];
    for (int i1=1; i1<n1-1; ++i1) {
      d[i1] = 1.0f-e[i1]-e[i1+1];
    }
    d[n1-1] = 1.0f-e[n1-1];

    // Factor A = L * inv(D) * L', where L is lower unit bidiagonal matrix.
    d[0] = 1.0f/d[0];
    for (int i1=1; i1<n1; ++i1) {
      float t = e[i1];
      e[i1] = t*d[i1-1];
      d[i1] = 1.0f/(d[i1]-t*e[i1]);
    }

    // y = inv(L) * x.
    y[0] = x[0];
    for (int i1=1; i1<n1; ++i1) {
      y[i1] = x[i1]-e[i1]*y[i1-1];
    }

    // y = D * inv(L) * x.
    for (int i1=0; i1<n1; ++i1) {
      y[i1] *= d[i1];
    }

    // y = L' * D * inv(L) * x.
    for (int i1=n1-1; i1>0; --i1) {
      y[i1-1] -= e[i1]*y[i1];
    }
  }

  public void applyKill(float[] ds, float[] x, float[] y) {
    applyPass(ds,x,y);
    Array.sub(x,y,y);
  }

  /**
   * Applies a filter that enhances (passes) features that are locally 
   * linear with inline vectors v.
   * @param ds diffusivity scale factors; null, for no scaling.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  public void applyPass(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] u1 = Array.rampfloat(0.0f,1.0f,n1);
    float[] t1 = new float[n1];
    float[] y1 = new float[n1];
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0f,0.0f);
    for (int i2=n2-1; i2>=1; --i2) {
      float[] yi2m1 = (i2>0)?y[i2-1]:y[0];
      for (int i1=0; i1<n1; ++i1) {
        float dsi = (ds!=null)?ds[i2][i1]:1.0f;
        float a = atable(dsi);
        float b = 1.0f-a;
        y[i2][i1] += sqrt(b)*x[i2][i1];
        y1[i1] = a*y[i2][i1];
      }
      maket(v1[i2],u1,t1);
      si.setUniformSamples(yi2m1);
      si.accumulate(n1,t1,y1);
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float dsi = (ds!=null)?ds[i2][i1]:1.0f;
        float a = atable(dsi);
        float b = 1.0f-a;
        y[i2][i1] *= b;
      }
    }
    for (int i2=1; i2<n2; ++i2) {
      float[] yi2m1 = (i2>0)?y[i2-1]:y[0];
      maket(v1[i2],u1,t1);
      si.setUniformSamples(yi2m1);
      si.interpolate(n1,t1,y1);
      for (int i1=0; i1<n1; ++i1) {
        float dsi = (ds!=null)?ds[i2][i1]:1.0f;
        float a = atable(dsi);
        y[i2][i1] += a*y1[i1];
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float dsi = (ds!=null)?ds[i2][i1]:1.0f;
        float a = atable(dsi);
        float b = 1.0f-a;
        y[i2][i1] *= sqrt(b);
      }
    }
  }

  // y = inv(I - A S )*(I - A )x
  // z = (I - A')*inv(I - S'A')y
  //
  // y = (I - A') * inv(I - S' * A')x
  // z = inv(I - A  * S ) * (I - A )y
  //
  // Solve: y[n-1] = x[n] + S'*A'*y[n]
  // Scale: y[n] = (I-A)*(I-A')*y[n]
  // Solve: z[n] = A*S*z[n-1]+y[n]

  /**
   * Applies a filter that attenuates (kills) features that are locally 
   * linear with inline vectors v.
   * @param ds diffusivity scale factors; null, for no scaling.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  public void applyKill(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    applyPass(ds,v1,x,y);
    Array.sub(x,y,y);
  }

  /**
   * Encodes specified fractions as 8-bit byte percentages.
   * Fractions are clipped to lie in the range [0,1].
   * @param s array of fractions.
   * @return array of 8-bit (byte) percentages.
   */
  public static byte[] encodeFractions(float[] s) {
    int n = s.length;
    byte[] b = new byte[n];
    for (int i=0; i<n; ++i) {
      float si = s[i];
      if (si<0.0f) {
        b[i] = 0;
      } else if (si>1.0f) {
        b[i] = 100;
      } else {
        b[i] = (byte)(si*100+0.5f);
      }
    }
    return b;
  }

  /**
   * Encodes specified fractions as 8-bit byte percentages.
   * Fractions are clipped to lie in the range [0,1].
   * @param s array of fractions.
   * @return array of 8-bit (byte) percentages.
   */
  public static byte[][] encodeFractions(float[][] s) {
    int n = s.length;
    byte[][] b = new byte[n][];
    for (int i=0; i<n; ++i) {
      b[i] = encodeFractions(s[i]);
    }
    return b;
  }

  /**
   * Encodes specified fractions as 8-bit byte percentages.
   * Fractions are clipped to lie in the range [0,1].
   * @param s array of fractions.
   * @return array of 8-bit (byte) percentages.
   */
  public static byte[][][] encodeFractions(float[][][] s) {
    int n = s.length;
    byte[][][] b = new byte[n][][];
    for (int i=0; i<n; ++i) {
      b[i] = encodeFractions(s[i]);
    }
    return b;
  }

  /**
   * Encodes specified unit vectors as 16-bit (short) indices.
   * @param u1 array of u1-components of unit vectors.
   * @param u2 array of u2-components of unit vectors.
   * @param u3 array of u3-components of unit vectors.
   * @return array of 16-bit (short) indices.
   */
  public static short[] encodeUnitVectors(float[] u1, float[] u2, float[] u3) {
    return UnitSphereSampling.encode16(u3,u2,u1);
  }

  /**
   * Encodes specified unit vectors as 16-bit (short) indices.
   * @param u1 array of u1-components of unit vectors.
   * @param u2 array of u2-components of unit vectors.
   * @param u3 array of u3-components of unit vectors.
   * @return array of 16-bit (short) indices.
   */
  public static short[][] encodeUnitVectors(
    float[][] u1, float[][] u2, float[][] u3) 
  {
    return UnitSphereSampling.encode16(u3,u2,u1);
  }

  /**
   * Encodes specified unit vectors as 16-bit (short) indices.
   * @param u1 array of u1-components of unit vectors.
   * @param u2 array of u2-components of unit vectors.
   * @param u3 array of u3-components of unit vectors.
   * @return array of 16-bit (short) indices.
   */
  public static short[][][] encodeUnitVectors(
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    return UnitSphereSampling.encode16(u3,u2,u1);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final int _ntable = 101;
  private static final float _stable = (float)(_ntable-1);
  private float[] _atable = new float[_ntable];
  private float _sigma; // filter half-width

  private void makeTable() {
    for (int itable=0; itable<_ntable; ++itable) {
      float sigmai = _sigma*(float)itable/(float)(_ntable-1);
      _atable[itable] = sigmai/(sigmai+sqrt(2.0f));
    }
  }
  private float atable(float ds) {
    int itable = (int)(ds*_stable);
    return _atable[itable];
  }

  private static void maket(float[] v1, float[] u1, float[] t1) {
    int n1 = v1.length;
    float t1min = 0.0f;
    float t1max = (float)(n1-1);
    for (int i1=0; i1<n1; ++i1) {
      float v1i = v1[i1];
      float v2i = sqrt(1.0f-v1i*v1i);
      float vi = v1i/v2i;
      t1[i1] = max(t1min,min(t1max,u1[i1]-vi));
    }
  }
}
