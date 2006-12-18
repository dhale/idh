/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local plane prediction error filtering of images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.12.12
 */
public class LocalPlaneFilter {

  public enum Type {
    QUAD,
    HALE,
    HALE0,
    HALE1,
    FOMEL,
    FOMEL2,
  };

  public LocalPlaneFilter(double sigma) {
    this(sigma,Type.HALE);
  }

  public LocalPlaneFilter(double sigma, Type type) {
    _rgfGradient = new RecursiveGaussianFilter(1.0);
    _rgfSmoother = new RecursiveGaussianFilter(sigma);
    if (type==Type.QUAD) {
      _filter = new QuadFilter();
    } else if (type==Type.HALE) {
      _filter = new HaleFilter();
    } else if (type==Type.HALE0) {
      _filter = new Hale0Filter();
    } else if (type==Type.HALE1) {
      _filter = new Hale1Filter();
    } else if (type==Type.FOMEL) {
      _filter = new FomelFilter();
    } else if (type==Type.FOMEL2) {
      _filter = new Fomel2Filter();
    }
  }

  public float[][][] find(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;

    // Gradient.
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    _rgfGradient.apply1X(x,g1);
    _rgfGradient.applyX1(x,g2);

    // Gradient products.
    float[][] g11 = new float[n2][n1];
    float[][] g12 = g1;
    float[][] g22 = g2;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        g11[i2][i1] = g1i*g1i;
        g12[i2][i1] = g1i*g2i;
        g22[i2][i1] = g2i*g2i;
      }
    }
    
    // Smoothed gradient products comprise the structure tensor.
    float[][] gtt = new float[n2][n1];
    _rgfSmoother.apply0X(g11,gtt);
    _rgfSmoother.applyX0(gtt,g11);
    _rgfSmoother.apply0X(g12,gtt);
    _rgfSmoother.applyX0(gtt,g12);
    _rgfSmoother.apply0X(g22,gtt);
    _rgfSmoother.applyX0(gtt,g22);

    // For each sample, the eigenvector corresponding to the largest
    // eigenvalue is normal to the plane. The size of that eigenvalue
    // is a measure of planarity in the interval [0,1].
    float[][][] p = new float[3][][];
    p[0] = g11;
    p[1] = g12;
    p[2] = g22;
    float[][] a = new float[2][2];
    float[][] v = new float[2][2];
    float[] d = new float[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = g11[i2][i1];
        a[0][1] = g12[i2][i1];
        a[1][0] = g12[i2][i1];
        a[1][1] = g22[i2][i1];
        Eigen.solveSymmetric22(a,v,d);
        p[0][i2][i1] = d[0]/sqrt(d[0]*d[0]+d[1]*d[1]);
        float v1 = v[0][0];
        float v2 = v[0][1];
        if (v1<0.0f) {
          v1 = -v1;
          v2 = -v2;
        }
        p[1][i2][i1] = v1;
        p[2][i2][i1] = v2;
      }
    }
    return p;
  }

  public void applyForward(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float[] c = new float[6];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1];
    for (i1=1; i1<n1-1; ++i1) {
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1]+c[2]*x2[i1-1];
    }
    _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
    y2[i1] = c[1]*x2[i1]+c[2]*x2[i1-1];

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]
             + c[3]*x2m[i1+1]+c[4]*x2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
        y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]+c[2]*x2 [i1-1]
               + c[3]*x2m[i1+1]+c[4]*x2m[i1]+c[5]*x2m[i1-1];
      }
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      y2[i1] = c[1]*x2 [i1]+c[2]*x2 [i1-1]
             + c[4]*x2m[i1]+c[5]*x2m[i1-1];
    }
  }

  public void applyInverse(float[][][] p, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    TridiagonalFMatrix tm = new TridiagonalFMatrix(n1);
    float[] ta = tm.a();
    float[] tb = tm.b();
    float[] tc = tm.c();
    float[] r = new float[n1];
    float[] c = new float[6];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
        _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
        r[i1] = x2[i1];
        ta[i1] = c[2];
        tb[i1] = c[1];
        tc[i1] = c[0];
    }
    tm.solve(r,y2);

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      int i1 = 0;
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      tb[i1] = c[1];
      tc[i1] = c[0];
      r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
        r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
        ta[i1] = c[2];
        tb[i1] = c[1];
        tc[i1] = c[0];
      }
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      ta[i1] = c[2];
      tb[i1] = c[1];
      r[i1] = x2[i1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      tm.solve(r,y2);
    }
  }

  public void smooth(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i1=n1-2; i1>=0; --i1) {
      y[n2-1][i1] = 0.0f;
    }
    for (int i2=n2-2; i2>=0; --i2) {
      y[i2][n1-1] = 0.0f;
      for (int i1=n1-2; i1>=0; --i1) {
        y[i2][i1] = 0.25f*(x[i2][i1]+x[i2][i1+1]+x[i2+1][i1]+x[i2+1][i1+1]);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float P99 = 0.999f;
  private static final float P98 = P99*P99;
  private static final float P49 = P99/2.0f;
  private static final float P48 = P98/2.0f;

  private RecursiveGaussianFilter _rgfGradient;
  private RecursiveGaussianFilter _rgfSmoother;
  private Filter _filter;

  /**
   * Different plane filters have different coefficients.
   * For forward filtering, the coefficients are used as follows:
   * <pre><code>
   * y[i2][i1] = c[0]*x[i2  ][i1+1]+c[1]*x[i2  ][i1]+c[2]*x[i2  ][i1-1]
   *           + c[3]*x[i2-1][i1+1]+c[4]*x[i2-1][i1]+c[5]*x[i2-1][i1-1]
   * </code></pre>
   * For inverse filtering, a corresponding tridiagonal system of 
   * equations is solved. Note that some filters may not be invertible. 
   * The corresponding finite-difference stencil is
   * <pre><code>
   *        i2-1  i2
   *   i1-1 c[5] c[2]
   *   i1   c[4] c[1]
   *   i1+1 c[3] c[0]
   * </code></pre>
   */
  private interface Filter {
    /**
     * Gets coefficients for specified normal unit-vector.
     * @param n1 component of normal for 1st dimension.
     * @param n2 component of normal for 2nd dimension.
     * @param c array of computed coefficients.
     */
    public void getCoefficients(float n1, float n2, float[] c);
  }
  private static class QuadFilter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float np = n1+n2;
      float nm = n1-n2;
      c[0] = 0.0f;
      c[1] = nm;
      c[2] = np;
      c[3] = 0.0f;
      c[4] = -np;
      c[5] = -nm;
    }
  }
  private static class Hale0Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      if (n2<0.0f) {
        c[0] = 0.0f;
        c[1] = n1-n2;
        c[2] = P99*n2;
      } else {
        c[0] = -P99*n2;
        c[1] = n1+n2;
        c[2] = 0.0f;
      }
      c[3] = 0.0f;
      c[4] = -P99*n1;
      c[5] = 0.0f;
    }
  }
  private static class Hale1Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      c[0] = -P49*n2*(n1+n2);
      c[1] =  1.0f;
      c[2] =  P49*n2*(n1-n2);
      c[3] =  0.0f;
      c[4] = -P99*n1*n1;
      c[5] =  0.0f;
    }
  }
  private static class HaleFilter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float np = n1+n2;
      float nm = n1-n2;
      c[0] = P49*nm*np;
      c[1] = 1.0f;
      c[2] = c[0];
      c[3] = -P48*np*np;
      c[4] = -P99*nm*np;
      c[5] = -P48*nm*nm;
    }
  }
  private static class FomelFilter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      if (n1<N1MIN)
        n1 = N1MIN;
      float s = n2/n1;
      c[0] = (1.0f-s)*(2.0f-s)/6.0f;
      c[1] = (2.0f+s)*(2.0f-s)/3.0f;
      c[2] = (1.0f+s)*(2.0f+s)/6.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
    private static final float N1MIN = 0.001f;
  }
  private static class Fomel2Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float t1 = 2.0f*n1;
      c[0] = (n1-n2)*(t1-n2)/6.0f;
      c[1] = (t1+n2)*(t1-n2)/3.0f;
      c[2] = (n1+n2)*(t1+n2)/6.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
  }
}
