/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * A local diffusion kernel for use in anisotropic diffusion filtering.
 * <p>
 * This kernel is a filter that computes y = y+G'DGx where G is the 
 * gradient operator, G' is its adjoint, and D is a local diffusion 
 * tensor that determines for each image sample the filter coefficients.
 * <p>
 * A local diffusion kernel is typically used in combinations with others.
 * For example, the filter implied by (I+G'DG)y = G'DGx acts as a notch
 * filter. It attenuates features for which G'DG is zero while preserving 
 * other features. The diffusion tensors D control the width, orientation,
 * and anisotropy of the spectral notch. Note that application of this filter 
 * requires solving a sparse symmetric positive-definite system of equations.
 * <p>
 * An even simpler example is the filter implied by (I+G'DG)y = x. This
 * filter smooths features in the directions implied by the tensors D.
 * Again, application of this filter requires solving a sparse symmetric 
 * positive-definite system of equations.
 * <p>
 * The accumulation of the kernel output in y = y+G'DGx is useful when
 * constructing such combination filters. Given y = 0, this kernel
 * computes y = G'DGx. Given y = x, it computes y = (I+G'DG)x.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.10
 */
public class LocalDiffusionKernel {

  /**
   * Constructs a local diffusion kernel.
   */
  public LocalDiffusionKernel() {
    this(1.0/12.0);
  }

  /**
   * Constructs a local diffusion kernel with an experimental factor.
   * @param rs the experimental factor for tuning gradient approximations.
   */
  public LocalDiffusionKernel(double rs) {
    _rs = (float)rs;
  }

  ///////////////////////////////////////////////////////////////////////////
  // 2-D

  /**
   * Computes y = y+G'DGx, for specified local diffusion tensors D.
   * @param ldt local diffusion tensors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void apply(LocalDiffusionTensors2 ldt, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] d = new float[3];
    for (int i2=0; i2<=n2; ++i2) {
      int j2 = (i2<n2)?i2:n2-1;
      for (int i1=0; i1<=n1; ++i1) {
        int j1 = (i1<n1)?i1:n1-1;
        ldt.getTensor(j1,j2,d);
        apply(_rs,d[0],d[1],d[2],i1,i2,n1,n2,x,y);
      }
    }
  }

  /**
   * Gets filter coefficients for this kernel.
   * @param ldt local diffusion 2x2 tensors.
   * @return arrays[5][n2][n1] of coefficients.
   *  The five arrays are {s00,s0p,spm,sp0,spp}, where s00 corresponds to 
   *  lag[0][0], s0p to lag[0][1], spm to lag[1][-1], and so on. 
   *  (Here, m, 0, and p denote minus, zero, and plus, respectively.)
   * @return the filter.
   */
  public float[][][] getCoefficients(LocalDiffusionTensors2 ldt) {
    int n1 = ldt.getN1();
    int n2 = ldt.getN2();
    float[][][] s = new float[5][n2][n1];
    float[][] s00 = s[0];
    float[][] s0p = s[1];
    float[][] spm = s[2];
    float[][] sp0 = s[3];
    float[][] spp = s[4];
    float[] d = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
        ldt.getTensor(i1p,i2p,d);
        float a = 0.5f*d[0];
        float b = 0.5f*d[1];
        float c = 0.5f*d[2];
        float t = 2.0f*_rs*(a+c);
        //float t = abs(b);
        s00[i2p][i1p] += a+b+c-t;
        s00[i2p][i1m] += a-b+c-t;
        s00[i2m][i1p] += a-b+c-t;
        s00[i2m][i1m] += a+b+c-t;
        s0p[i2p][i1m] -= a-t;
        s0p[i2m][i1m] -= a-t;
        spp[i2m][i1m] -= b+t;
        spm[i2m][i1p] += b-t;
        sp0[i2m][i1p] -= c-t;
        sp0[i2m][i1m] -= c-t;
        if (i1m==0) {
          c *= 2.0f;
          s00[i2p][i1m] += c;
          s00[i2m][i1m] += c;
          sp0[i2m][i1m] -= c;
        } else if (i1p==n1-1) {
          c *= 2.0f;
          s00[i2p][i1p] += c;
          s00[i2m][i1p] += c;
          sp0[i2m][i1p] -= c;
        }
        if (i2m==0) {
          a *= 2.0f;
          s00[i2m][i1p] += a;
          s00[i2m][i1m] += a;
          s0p[i2m][i1m] -= a;
        } else if (i2p==n2-1) {
          a *= 2.0f;
          s00[i2p][i1p] += a;
          s00[i2p][i1m] += a;
          s0p[i2p][i1m] -= a;
        }
      }
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////
  // 3-D


  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _rs; // experimental parameter for gradient approximation

  // Computes y = y+G'DGx for one sample.
  private static void apply(
   float rs, float d11, float d12, float d22,
   int i1, int i2, int n1, int n2, float[][] x, float[][] y) 
  {
    //float dd = 0.5f*abs(d12);
    float dd = rs*(d11+d22);
    float aa = (0.5f*d11-dd);
    float bm = (0.5f*d12-dd);
    float bp = (0.5f*d12+dd);
    float cc = (0.5f*d22-dd);
    int i1m = ( 0<i1)?i1-1:0;
    int i2m = ( 0<i2)?i2-1:0;
    int i1p = (i1<n1)?i1:n1-1;
    int i2p = (i2<n2)?i2:n2-1;
    float xpp = x[i2p][i1p];
    float xpm = x[i2p][i1m];
    float xmp = x[i2m][i1p];
    float xmm = x[i2m][i1m];
    float apppm = aa*(xpp-xpm);
    float ampmm = aa*(xmp-xmm);
    float bppmm = bp*(xpp-xmm);
    float bpmmp = bm*(xpm-xmp);
    float cppmp = cc*(xpp-xmp);
    float cpmmm = cc*(xpm-xmm);
    y[i2p][i1p] += apppm+bppmm+cppmp;
    y[i2p][i1m] -= apppm+bpmmp-cpmmm;
    y[i2m][i1p] += ampmm+bpmmp-cppmp;
    y[i2m][i1m] -= ampmm+bppmm+cpmmm;
  }

  // Computes y = y+G'DGx for one sample.
  // This old method assumes that rs = 0.25.
  private static void applyOld(
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

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void testApplyMethods() {
    int n1 = 3;
    int n2 = 3;
    float[][] x = Array.zerofloat(n1,n2); x[0][0] = 1.0f;
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    float theta = FLT_PI*2.0f/8.0f;
    float v1 = sin(theta);
    float v2 = cos(theta);
    float d0 = 0.0f;
    float d1 = 1.0f;
    float d11 = d0+d1*v1*v1;
    float d12 =    d1*v2*v1;
    float d22 = d0+d1*v2*v2;
    float rs = 1.0f/12.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        apply(rs,d11,d12,d22,i1,i2,n1,n2,x,y);
        applyOld(d11,d12,d22,i1,i2,x,z);
      }
    }
    Array.dump(y);
    Array.dump(z);
  }

  public static void testGetCoefficients() {
    int n1 = 5;
    int n2 = 5;
    float[][] x = Array.zerofloat(n1,n2); 
    x[0][0] = x[n2-1][0] = x[0][n1-1] = x[n2-1][n1-1] = 1.0f;
    //x[n2/2][n1/2] = 1.0f;
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    float theta = -FLT_PI*2.0f/8.0f;
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sin(theta),n1,n2);
    LocalDiffusionTensors2 ldt = new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    float rs = 1.0f/12.0f;
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(rs);
    LocalSpd9Filter lsf = new LocalSpd9Filter(ldk.getCoefficients(ldt));
    ldk.apply(ldt,x,y);
    lsf.apply(x,z);
    Array.dump(x);
    Array.dump(y);
    Array.dump(z);
  }

  public static void main(String[] args) {
    //testApplyMethods();
    testGetCoefficients();
  }
}
