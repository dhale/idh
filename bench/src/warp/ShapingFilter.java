/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A shaping filter h, such that convolution h * x ~ y for specified x and y.
 * The filter h shapes (approximately) a specified input sequence x into a
 * specified output sequence y. The approximation to minimize the sum of
 * squared differences between the sequences h*x and y.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.25
 */
public class ShapingFilter {

  public static float[] design(
    int nh, int kh, 
    int nx, int kx, float[] x, 
    int ny, int ky, float[] y) 
  {
    float[] cxx = new float[nh];
    float[] cxy = new float[nh];
    xcor(nx,kx,x,nx,kx,x,nh, 0,cxx);
    xcor(nx,kx,x,ny,ky,y,nh,kh,cxy);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cxx);
    return stm.solve(cxy);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing
  public static void main(String[] args) {
    int nx = 3;
    int kx = 0;
    float[] x = {1.0f,-1.8f,0.81f};
    int ny = 1;
    int ky = 0;
    float[] y = {1.0f};
    int nh = 101;
    int kh = -5;
    float[] h = design(nh,kh,nx,kx,x,ny,ky,y);
    int nz = 101;
    int kz = -5;
    float[] z = new float[nz];
    conv(nh,kh,h,nx,kx,x,nz,kz,z);
    dump(h);
    dump(z);
  }
}
