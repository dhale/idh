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
   * @param ds scale factors for diffusivity inline with unit vectors v;
   *  if null, this method uses constant ds = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  public void applyInlinePass(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveInline(ds,v1,x,y);
  }

  /**
   * Applies a filter that attenuates (kills) features that are locally 
   * linear with inline vectors v.
   * @param ds scale factors for diffusivity inline with unit vectors v;
   *  if null, this method uses constant ds = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  public void applyInlineKill(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    solveInline(ds,v1,x,y);
    Array.sub(x,y,y);
  }

  /**
   * Applies a filter that enhances (passes) features that are locally 
   * linear with inline vectors w.
   * The unit vectors w are specified by short indices iw that correspond 
   * to a 16-bit sampling of the unit-sphere.
   * @param ds scale factors for diffusivity inline with unit vectors w;
   *  if null, this method uses constant ds = 1.
   * @param iw unit-sphere sample indices of unit vectors w.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyInlinePass(
    float[][][] ds, short[][][] iw, float[][][] x, float[][][] y) 
  {
    solveInline(ds,iw,x,y);
  }

  /**
   * Applies a filter that enhances (passes) features that are locally 
   * planar with normal vectors u.
   * The unit vectors u are specified by short indices iu that correspond 
   * to a 16-bit sampling of the unit-sphere.
   * @param ds scale factors for diffusivity normal to unit vectors u;
   *  if null, this method uses constant ds = 1.
   * @param iu unit-sphere sample indices of unit vectors u.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void applyNormalPass(
    float[][][] ds, short[][][] iu, float[][][] x, float[][][] y) 
  {
    solveNormal(ds,iu,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  /**
   * Solves (I+G'DG)y = x for inline diffusion tensors D = dvv'.
   * @param ds scale factors for diffusivity in direction of unit vectors v;
   *  if null, this method uses constant ds = 1.
   * @param v1 array of 1st components of inline unit vectors.
   * @param x array with input image; must be distinct from y.
   * @param y array with output image; must be distinct from x.
   */
  protected void solveInline(
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

  protected void solveInline(
    float[][][] ds, short[][][] iw, float[][][] x, float[][][] y) 
  {
    Check.state(false,"method implemented in subclass");
  }

  protected void solveNormal(
    float[][][] ds, short[][][] iu, float[][][] x, float[][][] y) 
  {
    Check.state(false,"method implemented in subclass");
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // filter half-width
} 
