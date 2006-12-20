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
    HALE0,
    HALE1,
    HALE2,
    FOMEL1,
    FOMEL2,
  };

  public LocalPlaneFilter(double sigma) {
    this(sigma,Type.HALE2);
  }

  public LocalPlaneFilter(double sigma, Type type) {
    _rgfGradient = new RecursiveGaussianFilter(1.0);
    _rgfSmoother = new RecursiveGaussianFilter(sigma);
    if (type==Type.QUAD) {
      _filter = new QuadFilter();
    } else if (type==Type.HALE0) {
      _filter = new Hale0Filter();
    } else if (type==Type.HALE1) {
      _filter = new Hale1Filter();
    } else if (type==Type.HALE2) {
      _filter = new Hale2Filter();
    } else if (type==Type.FOMEL1) {
      _filter = new Fomel1Filter();
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
        ta[i1] = c[2];
        tb[i1] = c[1];
        tc[i1] = c[0];
        r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      }
      _filter.getCoefficients(p1[i2][i1],p2[i2][i1],c);
      ta[i1] = c[2];
      tb[i1] = c[1];
      r[i1] = x2[i1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      tm.solve(r,y2);
    }
  }

  public void applyForwardX(float[][][] p, float[][] x, float[][] y) {
    Hale2Filter filter = new Hale2Filter();
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
    filter.getCoefficients(i1,i2,p1,p2,c);
    y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1];
    for (i1=1; i1<n1-1; ++i1) {
      filter.getCoefficients(i1,i2,p1,p2,c);
      y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1]+c[2]*x2[i1-1];
    }
    filter.getCoefficients(i1,i2,p1,p2,c);
    y2[i1] = c[1]*x2[i1]+c[2]*x2[i1-1];

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      filter.getCoefficients(i1,i2,p1,p2,c);
      y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]
             + c[3]*x2m[i1+1]+c[4]*x2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        filter.getCoefficients(i1,i2,p1,p2,c);
        y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]+c[2]*x2 [i1-1]
               + c[3]*x2m[i1+1]+c[4]*x2m[i1]+c[5]*x2m[i1-1];
      }
      filter.getCoefficients(i1,i2,p1,p2,c);
      y2[i1] = c[1]*x2 [i1]+c[2]*x2 [i1-1]
             + c[4]*x2m[i1]+c[5]*x2m[i1-1];
    }
  }

  public void applyInverseX(float[][][] p, float[][] x, float[][] y) {
    Hale2Filter filter = new Hale2Filter();
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
      filter.getCoefficients(i1,i2,p1,p2,c);
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
      filter.getCoefficients(i1,i2,p1,p2,c);
      tb[i1] = c[1];
      tc[i1] = c[0];
      r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        filter.getCoefficients(i1,i2,p1,p2,c);
        ta[i1] = c[2];
        tb[i1] = c[1];
        tc[i1] = c[0];
        r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      }
      filter.getCoefficients(i1,i2,p1,p2,c);
      ta[i1] = c[2];
      tb[i1] = c[1];
      r[i1] = x2[i1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      tm.solve(r,y2);
    }
  }

  public void smooth(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int i2 = 0;
    for (i2=0; i2<n2-1; ++i2) {
      float[] yn = y[i2];
      float[] xn = x[i2];
      float[] xp = x[i2+1];
      int i1 = 0;
      float xnm = 0.0f;
      float xpm = 0.0f;
      float xnn = xn[i1];
      float xpn = xp[i1];
      float xnp = xn[i1+1];
      float xpp = xp[i1+1];
      yn[i1] = 0.25f*(xnn+xpn+xnp+xpp);
      for (i1=1; i1<n1-1; ++i1) {
        xnm = xnn;
        xpm = xpn;
        xnn = xnp;
        xpn = xpp;
        xnp = xn[i1+1];
        xpp = xp[i1+1];
        yn[i1] = 0.125f*(xnm+xpm+2.0f*(xnn+xpn)+xnp+xpp);
      }
      xnm = xnn;
      xpm = xpn;
      xnn = xnp;
      xpn = xpp;
      yn[i1] = 0.25f*(xnm+xpm+xnn+xpn);
    }
    float[] yn = y[i2];
    float[] xn = x[i2];
    int i1 = 0;
    float xnm = 0.0f;
    float xnn = xn[i1];
    float xnp = xn[i1+1];
    yn[i1] = 0.5f*(xnn+xnp);
    for (i1=1; i1<n1-1; ++i1) {
      xnm = xnn;
      xnn = xnp;
      xnp = xn[i1+1];
      yn[i1] = 0.25f*(xnm+2.0f*xnn+xnp);
    }
    xnm = xnn;
    xnn = xnp;
    yn[i1] = 0.5f*(xnm+xnn);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float P99 = 1.0000f;
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
  private static class Hale2Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float np = n1+n2;
      float nm = n1-n2;
      c[0] = 0.25f*nm*np;
      c[1] = 0.50f;
      c[2] = c[0];
      c[3] = -0.25f*np*np;
      c[4] = -0.50f*nm*np;
      c[5] = -0.25f*nm*nm;
    }
    // Experiment with averaging normal vectors for symmetric tridi system.
    public void getCoefficients(
      int i1, int i2, float[][] v1, float[][] v2, float[] c) {
      int n1 = v1[0].length;
      int n2 = v1.length;
      int i1n = i1;
      int i2n = i2;
      int i1m = max(i1-1,0);
      int i2m = max(i2-1,0);
      int i1p = min(i1+1,n1-1);
      float v1mm = v1[i2m][i1m];
      float v1mn = v1[i2m][i1n];
      float v1mp = v1[i2m][i1p];
      float v1nm = v1[i2n][i1m];
      float v1nn = v1[i2n][i1n];
      float v1np = v1[i2n][i1p];
      float v2mm = v2[i2m][i1m];
      float v2mn = v2[i2m][i1n];
      float v2mp = v2[i2m][i1p];
      float v2nm = v2[i2n][i1m];
      float v2nn = v2[i2n][i1n];
      float v2np = v2[i2n][i1p];
      float v1m = 0.25f*(v1mm+v1nm+v1mn+v1nn);
      float v1p = 0.25f*(v1mp+v1np+v1mn+v1nn);
      float v2m = 0.25f*(v2mm+v2nm+v2mn+v2nn);
      float v2p = 0.25f*(v2mp+v2np+v2mn+v2nn);
      float vsm = sqrt(v1m*v1m+v2m*v2m);
      v1m /= vsm;
      v2m /= vsm;
      float vsp = sqrt(v1p*v1p+v2p*v2p);
      v1p /= vsp;
      v2p /= vsp;
      float vmm = v1m-v2m;
      float vmp = v1m+v2m;
      float vpm = v1p-v2p;
      float vpp = v1p+v2p;
      c[0] = 0.5f*vpm*vpp;
      c[1] = 0.5f*(vmm*vmm+vpp*vpp);
      c[2] = 0.5f*vmm*vmp;
      c[3] = -0.5f*vpp*vpp;
      c[4] = -0.5f*(vmm*vmp+vpm*vpp);
      c[5] = -0.5f*vmm*vmm;
    }
  }
  private static class Fomel1Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      if (n1<N1MIN)
        n1 = N1MIN;
      float s = n2/n1;
      c[0] = (1.0f-s)*(2.0f-s)/12.0f;
      c[1] = (2.0f+s)*(2.0f-s)/6.0f;
      c[2] = (1.0f+s)*(2.0f+s)/12.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
    private static final float N1MIN = 0.001f;
  }
  private static class Fomel2Filter implements Filter {
    public void getCoefficients(float n1, float n2, float[] c) {
      float t1 = 2.0f*n1;
      c[0] = (n1-n2)*(t1-n2)/12.0f;
      c[1] = (t1+n2)*(t1-n2)/6.0f;
      c[2] = (n1+n2)*(t1+n2)/12.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
  }
}
