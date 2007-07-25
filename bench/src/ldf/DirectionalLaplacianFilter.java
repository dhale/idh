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
 * diffusivities d = 0.5*(ds*sigma)*(ds*sigma). The scale factor 0.5 makes
 * the Fourier transform of this filter approximate that of a Gaussian 
 * for small wavenumbers.
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
   * Constructs a directional Laplacian filter with nominal half-width sigma.
   * @param sigma the nominal half-width for this filter.
   */
  public DirectionalLaplacianFilter(double sigma) {
    _sigma = (float)sigma;
  }

  ///////////////////////////////////////////////////////////////////////////
  // 2-D

  /**
   * Computes y = y+G'DGx, where D = dvv' and G is the gradient operator. 
   * @param ds scale factor for diffusivity in direction of unit vector v.
   * @param v1 1st component of unit vector v.
   * @param v2 2nd component of unit vector v.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyInline(
    float ds, float v1, float v2,
    float[][] x, float[][] y) 
  {
    float ss = ds*_sigma;
    float sv = 0.5f*ss*ss;
    float d11 = sv*v1*v1;
    float d12 = sv*v1*v2;
    float d22 = sv*v2*v2;
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = d(I-uu') and G is the gradient operator. 
   * @param ds scale factor for diffusivity orthogonal to unit vector u.
   * @param u1 1st component of unit vector u.
   * @param u2 2nd component of unit vector u.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyNormal(
    float ds, float u1, float u2,
    float[][] x, float[][] y) 
  {
    float ss = ds*_sigma;
    float su = 0.5f*ss*ss;
    float d11 = su*(1.0f-u1*u1);
    float d22 = su*(1.0f-u2*u2);
    float d12 = su*u1*u2;
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = dvv' and G is the gradient operator. 
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
        float ssi = (ds!=null)?_sigma*ds[i2][i1]:_sigma;
        float svi = 0.5f*ssi*ssi;
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
        float ssi = (ds!=null)?_sigma*ds[i2][i1]:_sigma;
        float sui = 0.5f*ssi*ssi;
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float d11 = 1.0f-sui*u1i*u1i;
        float d12 =     -sui*u1i*u2i;
        float d22 = 1.0f-sui*u2i*u2i;
        apply(d11,d12,d22,i1,i2,x,y);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // 3-D

  /**
   * Computes y = y+G'DGx, where D = dww' and G is the gradient operator. 
   * @param ds scale factor for diffusivity in direction of unit vector w.
   * @param w1 1st component of unit vector w.
   * @param w2 2nd component of unit vector w.
   * @param w3 3rd component of unit vector w.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyInline(
    float ds, float w1, float w2, float w3, 
    float[][][] x, float[][][] y) 
  {
    float ss = ds*_sigma;
    float sw = 0.5f*ss*ss;
    float d11 = sw*w1*w1;
    float d12 = sw*w1*w2;
    float d13 = sw*w1*w3;
    float d22 = sw*w2*w2;
    float d23 = sw*w2*w3;
    float d33 = sw*w3*w3;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          apply(d11,d12,d13,d22,d23,d33,i1,i2,i3,x,y);
        }
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = d(I-uu') and G is the gradient operator. 
   * @param ds scale factor for diffusivity orthogonal to unit vector u.
   * @param u1 1st component of unit vector u.
   * @param u2 2nd component of unit vector u.
   * @param u3 3rd component of unit vector u.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyNormal(
    float ds, float u1, float u2, float u3, 
    float[][][] x, float[][][] y) 
  {
    float ss = ds*_sigma;
    float su = 0.5f*ss*ss;
    float d11 = su*(1.0f-u1*u1);
    float d22 = su*(1.0f-u2*u2);
    float d33 = su*(1.0f-u3*u3);
    float d12 = su*u1*u2;
    float d13 = su*u1*u3;
    float d23 = su*u2*u3;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          apply(d11,d12,d13,d22,d23,d33,i1,i2,i3,x,y);
        }
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = dww' and G is the gradient operator. 
   * Unit vectors w are specified by short indices iw that correspond to a 
   * 16-bit sampling of the unit-sphere.
   * @param ds scale factors for diffusivity in direction of unit vectors w;
   *  if null, this method uses constant ds = 1.
   * @param iw unit-sphere sample indices of unit vectors w.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyInline(
    float[][][] ds, short[][][] iw, 
    float[][][] x, float[][][] y) 
  {
    if (_uss16==null)
      _uss16 = new UnitSphereSampling(16);
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ssi = (ds!=null)?_sigma*ds[i3][i2][i1]:_sigma;
          float swi = 0.5f*ssi*ssi;
          float[] wi = _uss16.getPoint(iw[i3][i2][i1]); // {wx,wy,wz}
          float w1i = wi[2]; // w1 = wz
          float w2i = wi[1]; // w2 = wy
          float w3i = wi[0]; // w3 = wx
          float d11 = swi*w1i*w1i;
          float d12 = swi*w1i*w2i;
          float d13 = swi*w1i*w3i;
          float d22 = swi*w2i*w2i;
          float d23 = swi*w2i*w3i;
          float d33 = swi*w3i*w3i;
          apply(d11,d12,d13,d22,d23,d33,i1,i2,i3,x,y);
        }
      }
    }
  }

  /**
   * Computes y = y+G'DGx, where D = d(I-uu') and G is the gradient operator. 
   * Unit vectors u are specified by short indices iu that correspond to a 
   * 16-bit sampling of the unit-sphere.
   * @param ds scale factor for diffusivity orthogonal to unit vectors u;
   *  if null, this method uses constant ds = 1.
   * @param iu unit-sphere sample indices of unit vectors u.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyNormal(
    float[][][] ds, short[][][] iu, 
    float[][][] x, float[][][] y) 
  {
    if (_uss16==null)
      _uss16 = new UnitSphereSampling(16);
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ssi = (ds!=null)?_sigma*ds[i3][i2][i1]:_sigma;
          float sui = 0.5f*ssi*ssi;
          float[] ui = _uss16.getPoint(iu[i3][i2][i1]); // {ux,uy,uz}
          float u1i = ui[2]; // u1 = uz
          float u2i = ui[1]; // u2 = uy
          float u3i = ui[0]; // u3 = ux
          float d11 = sui*(1.0f-u1i*u1i);
          float d22 = sui*(1.0f-u2i*u2i);
          float d33 = sui*(1.0f-u3i*u3i);
          float d12 = sui*u1i*u2i;
          float d13 = sui*u1i*u3i;
          float d23 = sui*u2i*u3i;
          apply(d11,d12,d13,d22,d23,d33,i1,i2,i3,x,y);
        }
      }
    }
  }

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

  private static UnitSphereSampling _uss16; // maps indices to unit-vectors

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

  // Computes y = y+G'DGx for one sample.
  // Operations per sample for this method:
  //    16 loads + 8 stores + 4+12+6+8 adds + 3+9+4 muls
  //  = 16 loads + 8 stores +       30 adds +    16 muls
  // For alternative (more complicated) method with 27-point stencil:
  //    28 loads + 1 store  +       27 adds +    27 muls
  // This does not include the cost of computing the 27 coefficients!
  private void apply(
   float d11, float d12, float d13, float d22, float d23, float d33,
   int i1, int i2, int i3, float[][][] x, float[][][] y) 
  {
    float x000 = x[i3  ][i2  ][i1  ];
    float x001 = x[i3  ][i2  ][i1-1];
    float x010 = x[i3  ][i2-1][i1  ];
    float x100 = x[i3-1][i2  ][i1  ];
    float x011 = x[i3  ][i2-1][i1-1];
    float x101 = x[i3-1][i2  ][i1-1];
    float x110 = x[i3-1][i2-1][i1  ];
    float x111 = x[i3-1][i2-1][i1-1];
    //float x1 = 0.25f*(x000+x010+x100+x110-x001-x011-x101-x111);
    //float x2 = 0.25f*(x000+x001+x100+x101-x010-x011-x110-x111);
    //float x3 = 0.25f*(x000+x001+x010+x011-x100-x101-x110-x111);
    float xa = x000-x111;
    float xb = x001-x110;
    float xc = x010-x101;
    float xd = x100-x011;
    float x1 = 0.25f*(xa-xb+xc+xd);
    float x2 = 0.25f*(xa+xb-xc+xd);
    float x3 = 0.25f*(xa+xb+xc-xd);
    float y1 = d11*x1+d12*x2+d13*x3;
    float y2 = d12*x1+d22*x2+d23*x3;
    float y3 = d13*x1+d23*x2+d33*x3;
    float ya = 0.25f*(y1+y2+y3);
    float yb = 0.25f*(y1-y2+y3);
    float yc = 0.25f*(y1+y2-y3);
    float yd = 0.25f*(y1-y2-y3);
    y[i3  ][i2  ][i1  ] += ya;
    y[i3  ][i2  ][i1-1] -= yd;
    y[i3  ][i2-1][i1  ] += yb;
    y[i3-1][i2  ][i1  ] += yc;
    y[i3  ][i2-1][i1-1] -= yc;
    y[i3-1][i2  ][i1-1] -= yb;
    y[i3-1][i2-1][i1  ] += yd;
    y[i3-1][i2-1][i1-1] -= ya;
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
