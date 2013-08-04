/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Flattens and unflattens locally planar features in a 3D image.
 * <p>
 * We assume that a 3D image f(x1,x2,x3) has coherent features that are
 * locally planar with finite (non-vertical) slope. This class computes and
 * applies coordinate mappings for flattening and unflattening such features.
 * In this version, the mappings change only the vertical coordinate.
 * <p> 
 * Let x1 denote the vertical coordinate for the image f(x1,x2,x3). This class
 * computes a mapping function x1(u1,u2,u3) such that features in a
 * transformed image g(u1,u2,u3) = f(x1(u1,u2),u2,u3) are approximately
 * horizontal. This process is often called "flattening", and the transformed
 * image g(u1,u2,u3) is "the flattened image." Likewise, for any constant u1,
 * the curve defined by the function x1(u1,x2,x3) is called a "horizon".
 * <p>
 * If the coordinates u1, u2 and u3 are sampled finely enough, then the
 * mapping function x1(u1,u2,u3) is invertible. Specifically, there exists an
 * inverse mapping u1(x1,x2,x3) such that x1 = x1(u1(x1,x2,x3),x2,u3) for all
 * x2 and x3. Currently, samplings of u1, u2 and u3 are set to be identical to
 * those for x1, x2 and x3. If used directly, this sampling of the mapping
 * x1(u1,u2,u3) may cause aliasing of the flattened image g(u1,u2,u3), so that
 * f(x1,x2,x3) cannot be recovered by unflattening.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.02
 */
public class Flattener3X {

  /** Coordinate mappings u1(x1,x2,x3) and x1(u1,u2,u3). */
  public static class Mappings {
    
    /** Sampling for the 1st dimension (the vertical coordinate). */
    public Sampling s1;
    
    /** Sampling for the 2nd dimension. */
    public Sampling s2;
    
    /** Sampling for the 3rd dimension. */
    public Sampling s3;

    /** Array of sampled u1(x1,x2,x3). */
    public float[][][] u1;
    
    /** Array of sampled x1(u1,u2,u3). */
    public float[][][] x1;

    /**
     * Uses these mappings to flatten the specified image.
     * @param f the image to flatten.
     * @return the flattened image.
     */
    public float[][][] flatten(float[][][] f) {
      return apply(x1,f);
    }

    /**
     * Uses these mappings to unflatten the specified image.
     * @param f the image to unflatten.
     * @return the unflattened image.
     */
    public float[][][] unflatten(float[][][] f) {
      return apply(u1,f);
    }

    /**
     * Gets the flattening shifts s(u1,u2,u3) = u1 - x1(u1,u2,u3).
     * @return the flattening shifts.
     */
    public float[][][] getShiftsS() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      float[][][] s = new float[n3][n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float u1i = f1+i1*d1;
            s[i3][i2][i1] = u1i-x1[i3][i2][i1];
          }
        }
      }
      return s;
    }

    /**
     * Gets the unflattening shifts r(x1,x2,x3) = u1(x1,x2,x3) - x1.
     * @return the unflattening shifts.
     */
    public float[][][] getShiftsR() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      float[][][] r = new float[n3][n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float x1i = f1+i1*d1;
            r[i3][i2][i1] = u1[i3][i2][i1]-x1i;
          }
        }
      }
      return r;
    }

    private Mappings(
      Sampling s1, Sampling s2, Sampling s3, 
      float[][][] u1, float[][][] x1) 
    {
      this.s1 = s1;
      this.s2 = s2;
      this.s3 = s3;
      this.u1 = u1;
      this.x1 = x1;
    }

    private float[][][] apply(final float[][][] ux, final float[][][] f) {
      final int n1 = s1.getCount();
      final int n2 = s2.getCount();
      final int n3 = s3.getCount();
      final double d1 = s1.getDelta();
      final double f1 = s1.getFirst();
      final SincInterp si = new SincInterp();
      final float[][][] g = new float[n3][n2][n1];
      Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2)
          si.interpolate(n1,d1,f1,f[i3][i2],n1,ux[i3][i2],g[i3][i2]);
      }});
      return g;
    }
  }

  /**
   * Sets half-widths for smoothings in 1st, 2nd and 3rd dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma1 half-width for smoothing in 1st dimension, in samples.
   * @param sigma2 half-width for smoothing in 2nd dimension, in samples.
   * @param sigma3 half-width for smoothing in 3rd dimension, in samples.
   */
  public void setSmoothings(double sigma1, double sigma2, double sigma3) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma3;
  }

  /**
   * Sets parameters that control the number of solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of solver iterations.
   */
  public void setIterations(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  /**
   * Gets mappings computed from specified slopes and planarities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param s3 sampling of 3rd dimension.
   * @param p2 array of slopes of image features in 2nd dimension.
   * @param p3 array of slopes of image features in 3rd dimension.
   * @param ep array of planarities of image features.
   */
  public Mappings getMappingsFromSlopes(
    Sampling s1, Sampling s2, Sampling s3,
    float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    // Sampling parameters.
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    float d3 = (float)s3.getDelta();
    float f1 = (float)s1.getFirst();

    // If necessary, convert units for slopes to samples per sample.
    if (d1!=d2)
      p2 = mul(d2/d1,p2);
    if (d1!=d3)
      p3 = mul(d3/d1,p3);

    // Compute u(x1,x2,x3), in samples.
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] u = rampfloat(0.0f,1.0f,0.0f,0.0f,n1,n2,n3); // u, in samples
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vu = new VecArrayFloat3(u);
    Smoother3 smoother3 = new Smoother3(n1,n2,n3,_sigma1,_sigma2,_sigma3,ep);
    A3 a3 = new A3(smoother3,ep,p2,p3);
    double anorm = 1.0/vu.norm2();
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(ep,p2,p3,b);
    smoother3.applyTranspose(b);
    cs.solve(anorm,a3,vb,vu);
    smoother3.apply(u);
    cleanShifts(u);

    // Compute u1(x1,x2,x3) from u in samples.
    final float[][][] u1 = u;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          u1[i3][i2][i1] = f1+u[i3][i2][i1]*d1;
        }
      }
    }

    // Compute x1(u1,u2).
    final float[][][] x1 = b;
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});

    return new Mappings(s1,s2,s3,u1,x1);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][] wp, float[][][] p2, float[][][] p3) 
    {
      _s3 = s3;
      _wp = wp;
      _p2 = p2;
      _p3 = p3;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      _s3.apply(z);
      applyLhs(_wp,_p2,_p3,z,y);
      _s3.applyTranspose(y);
    }
    private Smoother3 _s3;
    private float[][][] _wp;
    private float[][][] _p2;
    private float[][][] _p3;
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _wp[0][0].length;
      int n2 = _wp[0].length;
      int n3 = _wp.length;
      float[][][] x = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] y = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] ax = zerofloat(n1,n2,n3);
      float[][][] ay = zerofloat(n1,n2,n3);
      VecArrayFloat3 vx = new VecArrayFloat3(x);
      VecArrayFloat3 vy = new VecArrayFloat3(y);
      VecArrayFloat3 vax = new VecArrayFloat3(ax);
      VecArrayFloat3 vay = new VecArrayFloat3(ay);
      apply(vx,vax);
      apply(vy,vay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("A3: yax="+yax+" xay="+xay);
      System.out.println("A3: xax="+xax+" yay="+yay);
    }
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother3 {
    public Smoother3(
      int n1, int n2, int n3, 
      float sigma1, float sigma2, float sigma3, float[][][] ep) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _ep = ep;
      //testSpd();
    }
    public void apply(float[][][] x) {
      smooth3(_sigma3,_ep,x);
      smooth2(_sigma2,_ep,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,_ep,x);
      smooth3(_sigma3,_ep,x);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _ep;
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _ep[0][0].length;
      int n2 = _ep[0].length;
      int n3 = _ep.length;
      float[][][] x = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] y = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] ax = copy(x);
      float[][][] ay = copy(y);
      VecArrayFloat3 vx = new VecArrayFloat3(x);
      VecArrayFloat3 vy = new VecArrayFloat3(y);
      VecArrayFloat3 vax = new VecArrayFloat3(ax);
      VecArrayFloat3 vay = new VecArrayFloat3(ay);
      apply(ax);
      apply(ay);
      applyTranspose(ax);
      applyTranspose(ay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("S3: yax="+yax+" xay="+xay);
      System.out.println("S3: xax="+xax+" yay="+yay);
    }
  }

  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }
  private static void smooth2(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }});
  }

  // Smoothing for dimension 3.
  private static void smooth3(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] s2 = (s!=null)?new float[n3][]:null;
      float[][] x2 = new float[n3][];
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }});
  }

  private static void makeRhs(
    float[][][] wp, float[][][] p2, float[][][] p3, float[][][] y) 
  {
    zero(y);
  }
  private static void applyLhs(
    final float[][][] wp, final float[][][] p2, final float[][][] p3,
    final float[][][] x, final float[][][] y) 
  {
    zero(y);
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,p2,p3,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,p2,p3,x,y);
    }});
  }
  private static void applyLhsSlice3(
    int i3,
    float[][][] wp, float[][][] p2, float[][][] p3,
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float p2i = p2[i3][i2][i1];
        float p3i = p3[i3][i2][i1];
        float wps = wpi*wpi;
        float p2s = p2i*p2i;
        float p3s = p3i*p3i;
        float d11 = wps*(p2s+p3s);
        float d12 = wps*p2i;
        float d13 = wps*p3i;
        float d22 = wps;
        float d33 = wps;
        float xa = 0.0f;
        float xb = 0.0f;
        float xc = 0.0f;
        float xd = 0.0f;
        xa += x00[i1 ];
        xd -= x00[i1m];
        xb += x01[i1 ];
        xc -= x01[i1m];
        xc += x10[i1 ];
        xb -= x10[i1m];
        xd += x11[i1 ];
        xa -= x11[i1m];
        float x1 = 0.25f*(xa+xb+xc+xd);
        float x2 = 0.25f*(xa-xb+xc-xd);
        float x3 = 0.25f*(xa+xb-xc-xd);
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2       ;
        float y3 = d13*x1       +d33*x3;
        float ya = 0.25f*(y1+y2+y3);
        float yb = 0.25f*(y1-y2+y3);
        float yc = 0.25f*(y1+y2-y3);
        float yd = 0.25f*(y1-y2-y3);
        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;
      }
    }
  }

  // Post-processing of computed u (in samples) to ensure monotonic u1.
  private static void cleanShifts(float[][][] u) {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          if (u[i3][i2][i1]<=u[i3][i2][i1-1]+0.01f)
            u[i3][i2][i1] = u[i3][i2][i1-1]+0.01f;
        }
      }
    }
    printStats("u raw",u);
    shiftAndScale(u);
    printStats("u scl",u);
  }
  private static void shiftAndScale(float[][][] u) {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    float umin = min(u);
    float umax = max(u);
    float scale = (n1-1.0f)/(umax-umin);
    float shift = -umin*scale;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          u[i3][i2][i1] = shift+scale*u[i3][i2][i1];
        }
      }
    }
  }
  private static void printStats(String s, float[][][] a) {
    trace(s+": min="+min(a)+" max="+max(a));
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
