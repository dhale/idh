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

  private void applyInline33(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;

    // Diffusivity scale factor = square-root of (sigma*sigma)/2.
    float scale = _sigma/sqrt(2.0f);

    // Initialize cache of products mm, mp, pp to compute filter stencil.
    float[][] cmm = new float[2][n1+1];
    float[][] cmp = new float[2][n1+1];
    float[][] cpp = new float[2][n1+1];
    float[] dsi2p = (ds!=null)?ds[0]:null;
    //updateInlineCache(scale,dsi2p,v1[0],cmm,cmp,cpp);

    // For all i2, ...
    for (int i2=0; i2<n2; ++i2) {

      // Update cache.
      int i2m = max(i2-1,0);
      int i20 = i2;
      int i2p = min(i2+1,n2m);
      dsi2p = (ds!=null)?ds[i2p]:null;
      if (i2<n2m) {
        updateInlineCache(scale,dsi2p,v1[i2p],cmm,cmp,cpp);
      } else {
        Array.zero(cmm[1]);
        Array.zero(cmp[1]);
        Array.zero(cpp[1]);
      }

      float[] mm0 = cmm[0], mm1 = cmm[1];
      float[] mp0 = cmp[0], mp1 = cmp[1];
      float[] pp0 = cpp[0], pp1 = cpp[1];

      // Compute and apply 3 x 3 stencil.
      float[] xm = x[i2m];
      float[] x0 = x[i20];
      float[] xp = x[i2p];
      for (int i1=0,j1=1; i1<n1; ++i1,++j1) {
        int i1m = max(i1-1,0);
        int i10 = i1;
        int i1p = min(i1+1,n1m);
        float mm01 = mm0[j1];
        float mm10 = mm1[i1];
        float pp00 = pp0[i1];
        float pp11 = pp1[j1];
        float mp00 = mp0[i1];
        float mp10 = mp1[i1];
        float mp01 = mp0[j1];
        float mp11 = mp1[j1];
        float amm = -pp00;
        float am0 =  mp00+mp01;
        float amp = -mm01;
        float a0m = -mp00-mp10;
        float a00 =  pp00+pp11+mm01+mm10;
        float a0p = -mp01-mp11;
        float apm = -mm10;
        float ap0 =  mp10+mp11;
        float app = -pp11;
        y[i2][i1] = amm*xm[i1m] + a0m*x0[i1m] + apm*xp[i1m] +
                    am0*xm[i10] + a00*x0[i10] + ap0*xp[i10] +
                    amp*xm[i1p] + a0p*x0[i1p] + app*xp[i1p];
      }
    }
  }
  private static void updateInlineCache(
    float scale, float[] ds, float[] v1,
    float[][] cmm, float[][] cmp, float[][] cpp) 
  {
    int n1 = v1.length;
    float[] cmm1 = cmm[0];  cmm[0] = cmm[1];  cmm[1] = cmm1;
    float[] cmp1 = cmp[0];  cmp[0] = cmp[1];  cmp[1] = cmp1;
    float[] cpp1 = cpp[0];  cpp[0] = cpp[1];  cpp[1] = cpp1;
    scale *= 0.5f; // for finite-difference approximation to gradient
    for (int i1=1; i1<n1; ++i1) {
      float dsi = (ds!=null)?scale*ds[i1]:scale;
      float v1i = v1[i1];
      float v2i = sqrt(1.0f-v1i*v1i);
      float vmi = dsi*(v1i-v2i);
      float vpi = dsi*(v1i+v2i);
      cmm1[i1] = vmi*vmi;
      cmp1[i1] = vmi*vpi;
      cpp1[i1] = vpi*vpi;
    }
  }

  // Tests y'(Ax) = (y'Ax)' = x'(A'y) = x'(Ay)
  public static void main(String[] args) {
    testStencil();
  }
  private static void testStencil() {
    int n1 = 5;
    int n2 = 5;
    //float[][] x = Array.rampfloat(1.0f,1.0f,1.0f,n1,n2);
    float[][] x = Array.randfloat(n1,n2);
    //float[][] x = Array.zerofloat(n1,n2);
    //x[n2-1][n1-1] = 1.0f;
    //x[n2/2][n1/2] = 1.0f;
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    float[][] ds = Array.randfloat(n1,n2);
    float[][] v1 = Array.randfloat(n1,n2);
    ds = null;
    v1 = Array.fillfloat(0.5f,n1,n2);
    DirectionalLaplacianFilter dlf = new DirectionalLaplacianFilter(1.0);
    dlf.applyInline(ds,v1,x,y);
    dlf.applyInline33(ds,v1,x,z);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(y);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(z);
    float e = Array.sum(Array.abs(Array.sub(z,y)));
    Array.dump(x);
    Array.dump(y);
    Array.dump(z);
    System.out.println("error = "+e);
  }
}
