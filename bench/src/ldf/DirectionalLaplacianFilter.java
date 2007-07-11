/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Directional Laplacian filter. The Laplacian operator appears in 
 * diffusion equations, and this filter is useful in the context of 
 * anisotropic diffusion filtering.
 * <p>
 * This filter computes y = y+G'DGx where G is the gradient operator, 
 * G' is its adjoint, and D is a local diffusion tensor that determines 
 * for each image sample the direction of the Laplacian filter.
 * <p>
 * For example, if D = dvv' for local diffusivities d and unit vectors v,
 * then G'DGx is zero for image features that are constant in the direction 
 * of v. The diffusivities d depend on local scale factors ds multiplied by 
 * a nominal filter half-width sigma. Specifically, for each sample, 
 * diffusivities d = 0.5*(ds*sigma)*(ds*sigma).
 * <p>
 * Alternatively, if D = d(I-UU'), then the right-hand side G'DGx is zero
 * for image features that are constant in all directions orthogonal to the 
 * unit vectors u. 
 * <p>
 * Directional Laplacian filters are rarely used alone. While zeroing some 
 * features in images, they tend to attenuate many other features as well. 
 * Therefore, these filters are typically used in combinations with others.
 * <p>
 * For example, the filter implied by (I+G'DG)y = G'DGx acts as a notch
 * filter. It attenuates features for which G'DGx is zero while preserving 
 * other features. Diffusivities d (inside D) control the width of the notch.
 * Note that application of this filter requires solving a sparse symmetric 
 * positive-definite system of equations.
 * <p>
 * An even simpler example is the filter implied by (I+G'DG)y = x. This
 * filter smooths features in the directions implied by the tensors D.
 * Again, application of this filter requires solving a sparse symmetric 
 * positive-definite system of equations.
 * <p>
 * The accumulation of the filter output in y = y+G'DGx is useful when
 * constructing such combination filters. Given y = 0, this filter 
 * computes y = G'DGx. Given y = x, this filter computes y = (I+G'DG)x.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.05
 */
public class DirectionalLaplacianFilter {

  /**
   * A 2-D directional Laplacian filters with a 3 x 3 (9-point) stencil.
   * Note that filters need not use these stencils when applied. Rather,
   * stencils are used in certain iterative methods for solving sparse 
   * positive-definite systems of equations such as (I+G'DG)y = x.
   */
  public interface Stencil33 {

    /**
     * Gets an array of stencil coefficients. With care taken near the ends 
     * of arrays, filters may be applied using this stencil as follows:
     * <pre><code>
     * y[i2][i1] += 
     *   a[i1][0]*x[i2-1][i1-1]+a[i1][3]*x[i2  ][i1-1]+a[i1][6]*x[i2+1][i1-1]
     *   a[i1][1]*x[i2-1][i1  ]+a[i1][4]*x[i2  ][i1  ]+a[i1][7]*x[i2+1][i1  ]
     *   a[i1][2]*x[i2-1][i1+1]+a[i1][5]*x[i2  ][i1+1]+a[i1][8]*x[i2+1][i1+1]
     * </code></pre>
     * Coefficients corresponding to samples off the ends of arrays are zero.
     * @param i2 sample index in 2nd dimension.
     * @param a array[n1][9] in which to get coefficients.
     */
    public void get(int i2, float[][] a);
  }

  /**
   * Constructs a directional Laplacian filter with nominal half-width sigma.
   * @param sigma the nominal half-width for this filter.
   */
  public DirectionalLaplacianFilter(double sigma) {
    _sigma = (float)sigma;
  }

  /**
   * Computes y = y+G'DGx, where D = dvv' and G is the gradient operator. 
   * The right-hand-side G'DGx is zero in the direction of the unit vectors v.
   * Diffusivities d are scale factors that multiply the nominal half-width 
   * sigma for this filter.
   * <p>
   * Only components v1 of the inline vectors v are specified; all components 
   * v2 = sqrt(1-v1*v1) are assumed to be non-negative.
   * @param ds scale factors for diffusivity in direction of unit vectors v;
   *  if null, this method uses constant ds = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float dsi = (ds!=null)?_sigma*ds[i2][i1]:_sigma;
        float svi = 0.5f*dsi*dsi;
        float v1i = v1[i2][i1];
        float v2i = sqrt(1.0f-v1i*v1i);
        float d11 = svi*v1i*v1i;
        float d12 = svi*v1i*v2i;
        float d22 = svi*v2i*v2i;
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = d(I-uu') and G is the gradient operator. 
   * The right-hand-side G'DGx is zero in all directions orthogonal to the 
   * unit vectors u. Diffusivities d are scale factors that multiply the 
   * nominal half-width sigma for this filter.
   * <p>
   * Only components u2 of the normal vectors u are specified; all components 
   * u1 = sqrt(1-u2*u2) are assumed to be non-negative.
   * @param ds scale factors for diffusivity in directions orthogonal to
   *  unit vectors u; if null, this method uses constant ds = 1.
   * @param u2 array of 2nd components of normal unit vectors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyNormal(
    float[][] ds, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float dsi = (ds!=null)?_sigma*ds[i2][i1]:_sigma;
        float sui = 0.5f*dsi*dsi;
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float d11 = 1.0f-sui*u1i*u1i;
        float d12 =     -sui*u1i*u2i;
        float d22 = 1.0f-sui*u2i*u2i;
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  /**
   * Makes filter stencil for G'DG.
   * @param ds scale factors for diffusivity in direction of unit vectors v;
   *  if null, this method uses constant ds = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @return the stencil.
   */
  public Stencil33 makeInlineStencil33(float[][] ds, float[][] v1) {
    return new InlineStencil33(_sigma,ds,v1);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // nominal filter half-width

  // Computes y = y+G'DGx for one sample.
  private void apply(
   float d11, float d12, float d22,
   int i1, int i2, float[][] x, float[][] y) 
  {
    float x00 = x[i2  ][i1  ];
    float x01 = x[i2  ][i1-1];
    float x10 = x[i2-1][i1  ];
    float x11 = x[i2-1][i1-1];
    float xa = x00-x11;
    float xb = x01-x10;
    float x1 = 0.5f*(xa-xb);
    float x2 = 0.5f*(xa+xb);
    float y1 = d11*x1+d12*x2;
    float y2 = d12*x1+d22*x2;
    float ya = 0.5f*(y1+y2);
    float yb = 0.5f*(y1-y2);
    y[i2  ][i1  ] += ya;
    y[i2  ][i1-1] -= yb;
    y[i2-1][i1  ] += yb;
    y[i2-1][i1-1] -= ya;
  }

  /**
   * Implements 3 x 3 stencil for inline filters.
   */
  private static class InlineStencil33 implements Stencil33 {
    InlineStencil33(float sigma, float[][] ds, float[][] v1) {
      _i2 = new int[]{-1,-1};
      _n1 = v1[0].length;
      _n2 = v1.length;
      _ds = ds;
      _v1 = v1;
      _cmm = new float[2][_n1+1];
      _cmp = new float[2][_n1+1];
      _cpp = new float[2][_n1+1];
      _scale = sigma/(2.0f*sqrt(2.0f));
    }
    public void get(int i2, float[][] a) {
      updateCache(i2);
      int k0 = (i2  )%2;
      int k1 = (i2+1)%2;
      float[] mm0 = _cmm[k0], mm1 = _cmm[k1];
      float[] mp0 = _cmp[k0], mp1 = _cmp[k1];
      float[] pp0 = _cpp[k0], pp1 = _cpp[k1];
      for (int i1=0,j1=1; i1<_n1; ++i1,++j1) {
        float mm01 = mm0[j1];
        float mm10 = mm1[i1];
        float pp00 = pp0[i1];
        float pp11 = pp1[j1];
        float mp00 = mp0[i1];
        float mp10 = mp1[i1];
        float mp01 = mp0[j1];
        float mp11 = mp1[j1];
        float[] ai = a[i1];
        ai[0] = -pp00;
        ai[1] =  mp00+mp01;
        ai[2] = -mm01;
        ai[3] = -mp00-mp10;
        ai[4] =  pp00+pp11+mm01+mm10;
        ai[5] = -mp01-mp11;
        ai[6] = -mm10;
        ai[7] =  mp10+mp11;
        ai[8] = -pp11;
      }
    }
    private int _n1,_n2; // dimensions of arrays to be filtered
    private float _scale; // constant scale factor to compute coefficients
    private float[][] _ds,_v1; // arrays[n2][n1] of filter coefficients
    private int[] _i2; // array[2] of indices i2 in cache cmm, cmp, and cpp
    private float[][] _cmm,_cmp,_cpp; // array[2][n1+1] for cached products
    private void updateCache(int i2) {
      if (i2!=_i2[0] && i2!=_i2[1])
        computeProducts(i2);
      ++i2;
      if (i2!=_i2[0] && i2!=_i2[1])
        computeProducts(i2);
    }
    private void computeProducts(int i2) {
      if (i2<=0 || i2>=_n2) {
        Array.zero(_cmm[i2%2]);
        Array.zero(_cmp[i2%2]);
        Array.zero(_cpp[i2%2]);
      } else {
        float[] cmm = _cmm[i2%2];
        float[] cmp = _cmp[i2%2];
        float[] cpp = _cpp[i2%2];
        if (_ds!=null) {
          float[] ds = _ds[i2];
          float[] v1 = _v1[i2];
          for (int i1=1; i1<_n1; ++i1) {
            float dsi = ds[i1]*_scale;
            float v1i = v1[i1];
            float v2i = sqrt(1.0f-v1i*v1i);
            float vmi = dsi*(v1i-v2i);
            float vpi = dsi*(v1i+v2i);
            cmm[i1] = vmi*vmi;
            cmp[i1] = vmi*vpi;
            cpp[i1] = vpi*vpi;
          }
        } else {
          float[] v1 = _v1[i2];
          float dsi = _scale;
          for (int i1=1; i1<_n1; ++i1) {
            float v1i = v1[i1];
            float v2i = sqrt(1.0f-v1i*v1i);
            float vmi = dsi*(v1i-v2i);
            float vpi = dsi*(v1i+v2i);
            cmm[i1] = vmi*vmi;
            cmp[i1] = vmi*vpi;
            cpp[i1] = vpi*vpi;
          }
        }
      }
      _i2[i2%2] = i2;
    }
  }

  private void applyInlineStencil33(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] a = new float[n1][9];
    Stencil33 s33 = new InlineStencil33(_sigma,ds,v1);
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2m);
      float[] xm = x[i2m];
      float[] x0 = x[i2 ];
      float[] xp = x[i2p];
      s33.get(i2,a);
      for (int i1=0; i1<n1; ++i1) {
        int i1m = max(i1-1,0);
        int i1p = min(i1+1,n1m);
        float[] ai = a[i1];
        y[i2][i1] += ai[0]*xm[i1m] + ai[3]*x0[i1m] + ai[6]*xp[i1m] +
                     ai[1]*xm[i1 ] + ai[4]*x0[i1 ] + ai[7]*xp[i1 ] +
                     ai[2]*xm[i1p] + ai[5]*x0[i1p] + ai[8]*xp[i1p];
      }
    }
  }

  private void applyInlineStencil33Reverse(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] a = new float[n1][9];
    Stencil33 s33 = new InlineStencil33(_sigma,ds,v1);
    for (int i2=n2-1; i2>=0; --i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2m);
      float[] xm = x[i2m];
      float[] x0 = x[i2 ];
      float[] xp = x[i2p];
      s33.get(i2,a);
      for (int i1=n1-1; i1>=0; --i1) {
        int i1m = max(i1-1,0);
        int i1p = min(i1+1,n1m);
        float[] ai = a[i1];
        y[i2][i1] += ai[0]*xm[i1m] + ai[3]*x0[i1m] + ai[6]*xp[i1m] +
                     ai[1]*xm[i1 ] + ai[4]*x0[i1 ] + ai[7]*xp[i1 ] +
                     ai[2]*xm[i1p] + ai[5]*x0[i1p] + ai[8]*xp[i1p];
      }
    }
  }

  // Tests y'(Ax) = (y'Ax)' = x'(A'y) = x'(Ay)
  public static void main(String[] args) {
    testStencil();
    benchStencil();
  }
  private static void testStencil() {
    int n1 = 5;
    int n2 = 7;
    float[][] ds = Array.randfloat(n1,n2);
    float[][] v1 = Array.randfloat(n1,n2);
    float[][] x = Array.randfloat(n1,n2);
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    DirectionalLaplacianFilter dlf = new DirectionalLaplacianFilter(1.0);
    dlf.applyInline(ds,v1,x,y);
    dlf.applyInlineStencil33Reverse(ds,v1,x,z);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(y);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(z);
    //Array.dump(x);
    //Array.dump(y);
    //Array.dump(z);
    float e = Array.max(Array.abs(Array.sub(z,y)));
    System.out.println("error = "+e);
  }
  private static void benchStencil() {
    int napply;
    double maxtime = 2.0;
    int n1 = 1000;
    int n2 = 1000;
    float[][] ds = Array.randfloat(n1,n2);
    float[][] v1 = Array.randfloat(n1,n2);
    float[][] x = Array.randfloat(n1,n2);
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    DirectionalLaplacianFilter dlf = new DirectionalLaplacianFilter(1.0);
    Stopwatch sw = new Stopwatch();
    sw.restart();
    for (napply=0; sw.time()<maxtime; ++napply) {
      Array.zero(y);
      dlf.applyInline(ds,v1,x,y);
    }
    sw.stop();
    System.out.println(" simple rate = "+napply/sw.time());
    sw.restart();
    for (napply=0; sw.time()<maxtime; ++napply) {
      Array.zero(z);
      dlf.applyInlineStencil33(ds,v1,x,z);
    }
    sw.stop();
    System.out.println("stencil rate = "+napply/sw.time());
    edu.mines.jtk.mosaic.SimplePlot.asPixels(Array.sub(z,y));
    float e = Array.max(Array.abs(Array.sub(z,y)));
    System.out.println("error = "+e);
  }
}
