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
 * Local anisotropic diffusion filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.07
 */
public class LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   */
  public LocalDiffusionFilter(double sigma) {
    _sigma = (float)sigma;
  }

  /**
   * Applies a filter that enhances (passes) features that are locally 
   * linear with inline vectors v.
   * @param ds diffusivity scale factors; null, for no scaling.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  public void applyLinearPass(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveLinear(ds,v1,x,y);
  }

  /**
   * Applies a filter that attenuates (kills) features that are locally 
   * linear with inline vectors v.
   * @param ds diffusivity scale factors; null, for no scaling.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  public void applyLinearKill(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveLinear(ds,v1,x,y);
    Array.sub(x,y,y);
  }

  /**
   * Applies a filter that enhances (passes) features that are locally 
   * linear with inline vectors w.
   * Diffusivities d depend on a percentage of the nominal filter half-width 
   * sigma; these percentages are specified by byte values in the array is.
   * Inline vectors w are specified by short indices in the array iw that 
   * correspond to a 16-bit sampling of the unit-sphere.
   * @param is diffusivity scaling percentages; null for no scaling.
   * @param iw unit-sphere 16-bit sample indices for unit vectors w.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyLinearPass(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    solveLinear(is,iw,x,y);
  }

  /**
   * Applies a filter that attenuates (kills) features that are locally 
   * linear with inline vectors w.
   * Diffusivities d depend on a percentage of the nominal filter half-width 
   * sigma; these percentages are specified by byte values in the array is.
   * Inline vectors w are specified by short indices in the array iw that 
   * correspond to a 16-bit sampling of the unit-sphere.
   * @param is diffusivity scaling percentages; null for no scaling.
   * @param iw unit-sphere 16-bit sample indices for unit vectors w.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyLinearKill(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    solveLinear(is,iw,x,y);
    Array.sub(x,y,y);
  }

  /**
   * Applies a filter that enhances (passes) features that are locally 
   * planar with normal vectors u.
   * Diffusivities d depend on a percentage of the nominal filter half-width 
   * sigma; these percentages are specified by byte values in the array is.
   * Normal vectors u are specified by short indices in the array iu that 
   * correspond to a 16-bit sampling of the unit-sphere.
   * @param is diffusivity scaling percentages; null for no scaling.
   * @param iu unit-sphere 16-bit sample indices for unit vectors u.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyPlanarPass(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    solvePlanar(is,iu,x,y);
  }

  /**
   * Applies a filter that attenuates (kills) features that are locally 
   * planar with normal vectors u.
   * Diffusivities d depend on a percentage of the nominal filter half-width 
   * sigma; these percentages are specified by byte values in the array is.
   * Normal vectors u are specified by short indices in the array iu that 
   * correspond to a 16-bit sampling of the unit-sphere.
   * @param is diffusivity scaling percentages; null for no scaling.
   * @param iu unit-sphere 16-bit sample indices for unit vectors u.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyPlanarKill(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    solvePlanar(is,iu,x,y);
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
  // protected

  /**
   * Solves (I+G'DG)y = x for inline diffusion tensors D = dvv'.
   * @param ds diffusivity scale factors; null, for constant = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  protected void solveLinear(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float dsmax = (ds!=null)?Array.max(ds):1.0f;
    float ssmax = dsmax*_sigma;
    int ns = 1+(int)(ssmax*ssmax);
    float sigma = _sigma/sqrt(ns);
    Array.copy(x,y);
    //float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    //float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    float r = 0.5f;
    float s = 0.5f;
    for (int is=0; is<ns; ++is,x=y) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float dsi = (ds!=null)?sigma*ds[i2][i1]:sigma;
          float svi = 0.5f*dsi*dsi;
          float v1i = v1[i2][i1];
          float v2i = sqrt(1.0f-v1i*v1i);
          float d11 = svi*v1i*v1i;
          float d12 = svi*v1i*v2i;
          float d22 = svi*v2i*v2i;
          float x00 = x[i2  ][i1  ];
          float x01 = x[i2  ][i1-1];
          float x10 = x[i2-1][i1  ];
          float x11 = x[i2-1][i1-1];
          float x1 = r*(x00-x01)+s*(x10-x11);
          float x2 = r*(x00-x10)+s*(x01-x11);
          float y1 = d11*x1+d12*x2;
          float y2 = d12*x1+d22*x2;
          y[i2  ][i1  ] -= r*y1+r*y2;
          y[i2  ][i1-1] += r*y1-s*y2;
          y[i2-1][i1  ] -= s*y1-r*y2;
          y[i2-1][i1-1] += s*y1+s*y2;
        }
      }
    }
  }

  protected void solveLinear(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    Check.state(false,"method implemented in subclass");
  }

  protected void solvePlanar(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    Check.state(false,"method implemented in subclass");
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // filter half-width
} 
