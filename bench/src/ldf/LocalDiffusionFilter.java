/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic diffusion filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.04.08
 */
public class LocalDiffusionFilter {

  public LocalDiffusionFilter(double sigma) {
    _sigma = (float)sigma;
  }

  /**
   * Applies this local anisotropic diffusion filter.
   * Input and output arrays must be distinct.
   * @param su diffusivities in direction of normal vector u.
   * @param sv diffusivities in direction of normal vector v.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */

  public void apply(
    float[][] su, float[][] sv, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    int ns = 1+(int)(_sigma*_sigma);
    float ss = 0.5f*(_sigma*_sigma)/(float)ns;
    float[] xi20 = new float[n1];
    float[] xi21 = new float[n1];
    for (int is=0; is<ns; ++is,x=y) {
      for (int i2=0; i2<n2; ++i2) {
        float[] xtmp = xi20;  xi20 = xi21;  xi21 = xtmp;
        for (int i1=0; i1<n1; ++i1)
          xi20[i1] = x[i2][i1];
        for (int i1=0; i1<n1; ++i1) {
          float sui = su[i2][i1];
          float svi = sv[i2][i1];
          float u2i = u2[i2][i1];
          float u1i = sqrt(1.0f-u2i*u2i);
          float v2i =  u1i;
          float v1i = -u2i;
          float a11 = ss*(sui*u1i*u1i+svi*v1i*v1i);
          float a12 = ss*(sui*u1i*u2i+svi*v1i*v2i);
          float a22 = ss*(sui*u2i*u2i+svi*v2i*v2i);
          float x00 = xi20[i1];
          float x10 = xi21[i1];
          float x01 = (i1>0)?xi20[i1-1]:0.0f;
          float x11 = (i1>0)?xi21[i1-1]:0.0f;
          float xa = x00-x11;
          float xb = x01-x10;
          float x1 = 0.5f*(xa-xb);
          float x2 = 0.5f*(xa+xb);
          float y1 = a11*x1+a12*x2;
          float y2 = a12*x1+a22*x2;
          float ya = 0.5f*(y1+y2);
          float yb = 0.5f*(y1-y2);
          y[i2][i1] = x00-ya;
          if (i1>0) y[i2][i1-1] += yb;
          if (i2>0) y[i2-1][i1] -= yb;
          if (i2>0 && i1>0) y[i2-1][i1-1] += ya;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  float _sigma;
}
