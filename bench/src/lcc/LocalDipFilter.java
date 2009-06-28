/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local dip filtering of images. Local dip filters are linear filters
 * designed to attenuate events with a specified dip that may vary from
 * sample to sample. Because their coefficients are not constant, local
 * dip filters are not shift-invariant and therefore cannot be implemented
 * efficiently with fast Fourier transforms.
 * <p>
 * Local dip filters are implemented with finite-difference approximations 
 * to derivatives. These approximations are best for low wavenumbers (or
 * low frequencies), and are used to construct directional derivatives.
 * A basic filter is simply a directional Laplacian that zeros specified
 * dips, but also attenuates a rather broad range of nearby dips.
 * <p>
 * More useful dip and notch filters use this basic filter in combinations
 * with its inverse to limit attenuation primarily to the specified dips.
 * Small parameters may be specified to control the range of dips attenuated
 * by dip or notch combination filters.
 * <p>
 * Contours of constant amplitude for notch filters are parallel lines
 * near the origin in the wavenumber domain. For dip filters, these
 * contours are radial lines.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.03.20
 */
public class LocalDipFilter {

  /**
   * Use of minimum-phase filter factorization.
   * <em>Factorization is not yet implemented for 3-D filters.</em>
   * <dl>
   * <dt>NOT<dd>
   * not used
   * <dt>PCG<dd>
   * used to precondition inverse filter by conjugate gradients
   * <dt>INV<dd>
   * used directly to implement inverse filter only
   * <dt>ALL<dd>
   * used for all forward and inverse filters
   * </dl>
   */
  public enum Factor {
    NOT,
    PCG,
    INV,
    ALL
  };

  /**
   * Constructs a local dip filter that does not use factorization. 
   */
  public LocalDipFilter() {
    this(Factor.NOT);
  }

  /**
   * Constructs a local dip filter with specified factorization.
   * @param factor filter factorization.
   */
  public LocalDipFilter(Factor factor) {
    _factor = factor;
  }

  /**
   * Applies a line filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param vw array of notch filter widths.
   * @param v1 array of 1st components of line vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyLineFilter(
    float[][] vw, float[][] v1, float[][] x, float[][] y) 
  {
    float[][] t = new float[x.length][x[0].length];
    applyForwardNot(null,v1,x,t);
    applyInverseFac(vw,v1,t,y);
  }

  /**
   * Applies this local dip filter.
   * Input and output arrays must be distinct.
   * @param sd small parameter that controls width of dip filter.
   * @param sn small parameter that controls width of notch filter.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyForward(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    if (_factor!=Factor.ALL) {
      applyForwardNot(sd,sn,u2,x,y);
    } else {
      applyForwardFac(sd,sn,u2,x,y);
    }
  }

  /**
   * Applies the inverse of this local dip filter.
   * Input and output arrays must be distinct.
   * @param sd small parameter that controls width of dip filter.
   * @param sn small parameter that controls width of notch filter.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyInverse(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    if (_factor==Factor.NOT) {
      applyInverseNot(sd,sn,u2,x,y);
    } else if (_factor==Factor.PCG) {
      applyInversePcg(sd,sn,u2,x,y);
    } else {
      applyInverseFac(sd,sn,u2,x,y);
    }
  }

  /**
   * Applies a dip filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param sd small parameter that controls width of dip filter.
   *  This parameter is used for the inverse filter only.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyDip(
    float sd, float[][] u2, float[][] x, float[][] y) 
  {
    float[][] t = new float[x.length][x[0].length];
    applyForward(0.0f,0.0f,u2,x,t);
    applyInverse(sd,0.01f,u2,t,y);
  }

  /**
   * Applies a notch filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param sn small parameter that controls width of notch filter.
   *  This parameter is used for the inverse filter only.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyNotch(
    float sn, float[][] u2, float[][] x, float[][] y) 
  {
    float[][] t = new float[x.length][x[0].length];
    applyForward(0.0f,0.0f,u2,x,t);
    applyInverse(0.0f,sn,u2,t,y);
  }

  /**
   * Applies this local dip filter.
   * Input and output arrays must be distinct.
   * @param sd small parameter that controls width of dip filter.
   * @param sn small parameter that controls width of notch filter.
   * @param u2 array of 2nd components of normal vectors.
   * @param u3 array of 3rd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyForward(
    float sd, float sn, 
    float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    if (_factor!=Factor.ALL) {
      applyForwardNot(sd,sn,u2,u3,x,y);
    } else {
      applyForwardNot(sd,sn,u2,u3,x,y);
    }
  }

  /**
   * Applies the inverse of this local dip filter.
   * Input and output arrays must be distinct.
   * @param sd small parameter that controls width of dip filter.
   * @param sn small parameter that controls width of notch filter.
   * @param u2 array of 2nd components of normal vectors.
   * @param u3 array of 3rd components of normal vectors.
   * @param x array with input image
   * @param y array with output image
   */
  public void applyInverse(
    float sd, float sn, 
    float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    if (_factor==Factor.NOT) {
      applyInverseNot(sd,sn,u2,u3,x,y);
    } else if (_factor==Factor.PCG) {
      applyInverseNot(sd,sn,u2,u3,x,y);
    } else {
      applyInverseNot(sd,sn,u2,u3,x,y);
    }
  }

  /**
   * Applies a dip filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param sd small parameter that controls width of dip filter.
   *  This parameter is used for the inverse filter only.
   * @param u2 array of 2nd components of normal vectors.
   * @param u3 array of 3rd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyDip(
    float sd, float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    float[][][] t = new float[x.length][x[0].length][x[0][0].length];
    applyForward(0.0f,0.0f,u2,u3,x,t);
    applyInverse(sd,0.01f,u2,u3,t,y);
  }

  /**
   * Applies a notch filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param sn small parameter that controls width of notch filter.
   *  This parameter is used for the inverse filter only.
   * @param u2 array of 2nd components of normal vectors.
   * @param u3 array of 3rd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyNotch(
    float sn, float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    float[][][] t = new float[x.length][x[0].length][x[0][0].length];
    applyForward(0.0f,0.0f,u2,u3,x,t);
    applyInverse(0.0f,sn,u2,u3,t,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float CG_SMALL = 0.000001f;

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  private Factor _factor; // filter factorization
  private Map<Small,FactoredFilter> _ffmap; // minimum-phase factors
  private static FactoredFilter2 _ff2;

  private static class Small {
    float sd,sn;
    Small(float sd, float sn) {
      this.sd = sd;
      this.sn = sn;
    }
    public boolean equals(Object o) {
      Small that = (Small)o;
      return this.sd==that.sd && this.sn==that.sn;
    }
    public int hashCode() {
      return Float.floatToIntBits(sd)^Float.floatToIntBits(sn);
    }
  }

  private FactoredFilter getFactoredFilter(float sd, float sn) {
    Small small = new Small(sd,sn);
    if (_ffmap==null)
      _ffmap = new HashMap<Small,FactoredFilter>();
    FactoredFilter ff = _ffmap.get(small);
    if (ff==null) {
      ff = new FactoredFilter(small);
      _ffmap.put(small,ff);
    }
    return ff;
  }

  // Forward (not factored) filter.
  private void applyForwardNot(
    float[][] vw, float[][] v1, float[][] x, float[][] y) 
  {
    ArrayMath.zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float vwi = vw!=null?vw[i2][i1]:0.0f;
        float v1i = v1[i2][i1];
        float v2i = sqrt(1.0f-v1i*v1i);
        float d11 = v1i*v1i;
        float d12 = v1i*v2i;
        float d22 = v2i*v2i;
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
        y[i2  ][i1  ]  = ya+vwi*x[i2][i1];
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Inverse factored filter.
  private void applyInverseFac(
    float[][] vw, float[][] v1, float[][] x, float[][] y) 
  {
    if (_ff2==null)
      _ff2 = new FactoredFilter2();
    _ff2.applyInverse(vw,v1,x,y);
  }

  // Forward (not factored) filter.
  private void applyForwardNot(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    ArrayMath.zero(y);
    float aone = 1.0f+sd;
    float aeps = sn;
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float a11 = aone-u1i*u1i;
        float a12 =     -u1i*u2i;
        float a22 = aone-u2i*u2i;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.5f*(xa-xb);
        float x2 = 0.5f*(xa+xb);
        float y1 = a11*x1+a12*x2;
        float y2 = a12*x1+a22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ]  = ya+aeps*x[i2][i1];
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Directional isotropic Laplacian filter.
  // The first-derivative approximation used here is consistent
  // with an O(h^2) isotropic approximation to the Laplacian.
  // NOTE: if we use this, must also change preconditioner.
  private void applyForwardIso(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    float aone = 1.0f+sd;
    float aeps = sn;
    float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float a11 = aone-u1i*u1i;
        float a12 =     -u1i*u2i;
        float a22 = aone-u2i*u2i;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float x1 = r*(x00-x01)+s*(x10-x11);
        float x2 = r*(x00-x10)+s*(x01-x11);
        float y1 = a11*x1+a12*x2;
        float y2 = a12*x1+a22*x2;
        y[i2  ][i1  ]  = r*(y1+y2)+aeps*x[i2][i1];
        y[i2  ][i1-1] -= r*y1-s*y2;
        y[i2-1][i1  ] += s*y1-r*y2;
        y[i2-1][i1-1] -= s*(y1+y2);
      }
    }
  }

  // Forward factored filter.
  private void applyForwardFac(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    FactoredFilter ff = getFactoredFilter(sd,sn);
    ff.applyForward(u2,x,y);
  }

  // Inverse filter via conjugate gradients without preconditioning.
  private void applyInverseNot(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    trace("applyInverseNot: sd="+sd+" sn="+sn);
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    ArrayMath.zero(y);
    ArrayMath.copy(x,r);
    ArrayMath.copy(r,s);
    float rr = dot(r,r);
    float stop = rr*CG_SMALL;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForwardNot(sd,sn,u2,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      float rrold = rr;
      rr = dot(r,r);
      float beta = rr/rrold;
      for (int i2=0; i2<n2; ++i2) {
        float[] r2 = r[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = r2[i1]+beta*s2[i1];
      }
      //trace("niter="+niter+" rr="+rr);
    }
    trace("niter="+niter+" rr="+rr);
  }

  // Inverse filter with pre-conditioned conjugate gradients.
  private void applyInversePcg(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    trace("applyInversePcg: sd="+sd+" sn="+sn);
    FactoredFilter ff = getFactoredFilter(sd,sn);
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    float[][] w = new float[n2][n1]; // s
    ArrayMath.zero(y);
    ArrayMath.copy(x,r);
    ff.applyInverse(u2,r,s);
    float rr = dot(r,s);
    float stop = rr*CG_SMALL;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<100 && rr>stop; ++niter) {
      applyForwardNot(sd,sn,u2,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      ff.applyInverse(u2,r,w);
      float rrold = rr;
      rr = dot(r,w);
      float beta = rr/rrold;
      for (int i2=0; i2<n2; ++i2) {
        float[] w2 = w[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = w2[i1]+beta*s2[i1];
      }
      //trace("niter="+niter+" rr="+rr);
    }
    trace("niter="+niter+" rr="+rr);
  }

  // Inverse factored filter.
  private void applyInverseFac(
    float sd, float sn, float[][] u2, float[][] x, float[][] y) 
  {
    trace("applyInverseFac: sd="+sd+" sn="+sn);
    FactoredFilter ff = getFactoredFilter(sd,sn);
    ff.applyInverse(u2,x,y);
  }

  // Forward (not factored) filter.
  private void applyForwardNot(
    float sd, float sn, 
    float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    float aone = 1.0f+sd;
    float aeps = sn;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float u3i = u3[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u1i = sqrt(1.0f-u2i*u2i-u3i*u3i);
          float a11 = aone-u1i*u1i;
          float a22 = aone-u2i*u2i;
          float a33 = aone-u3i*u3i;
          float a12 =     -u1i*u2i;
          float a13 =     -u1i*u3i;
          float a23 =     -u2i*u3i;
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
          float y1 = a11*x1+a12*x2+a13*x3;
          float y2 = a12*x1+a22*x2+a23*x3;
          float y3 = a13*x1+a23*x2+a33*x3;
          float ya = 0.25f*(y1+y2+y3);
          float yb = 0.25f*(y1-y2+y3);
          float yc = 0.25f*(y1+y2-y3);
          float yd = 0.25f*(y1-y2-y3);
          y[i3  ][i2  ][i1  ]  = ya + aeps*x[i3][i2][i1];
          y[i3  ][i2  ][i1-1] -= yd;
          y[i3  ][i2-1][i1  ] += yb;
          y[i3-1][i2  ][i1  ] += yc;
          y[i3  ][i2-1][i1-1] -= yc;
          y[i3-1][i2  ][i1-1] -= yb;
          y[i3-1][i2-1][i1  ] += yd;
          y[i3-1][i2-1][i1-1] -= ya;
        }
      }
    }
  }

  // Inverse filter via conjugate gradients without preconditioning.
  private void applyInverseNot(
    float sd, float sn, 
    float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] r = new float[n3][n2][n1]; // r
    float[][][] s = new float[n3][n2][n1]; // d
    float[][][] t = new float[n3][n2][n1]; // q
    ArrayMath.zero(y);
    ArrayMath.copy(x,r);
    ArrayMath.copy(r,s);
    float rr = dot(r,r);
    float stop = rr*CG_SMALL;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForwardNot(sd,sn,u2,u3,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      float rrold = rr;
      rr = dot(r,r);
      float beta = rr/rrold;
      for (int i3=0; i3<n3; ++i3) {
        float[][] r3 = r[i3];
        float[][] s3 = s[i3];
        for (int i2=0; i2<n2; ++i2) {
          float[] r32 = r3[i2];
          float[] s32 = s3[i2];
          for (int i1=0; i1<n1; ++i1)
            s32[i1] = r32[i1]+beta*s32[i1];
        }
      }
      trace("niter="+niter+" rr="+rr);
    }
  }

  // A local line filter approximated with minimum-phase factors.
  // Factors are tabulated as a function of sigma and theta.
  private static class FactoredFilter2 {
    FactoredFilter2() {
      int maxlag = 4;
      /*
      int nlag = maxlag+2+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
        lag2[ilag] = (ilag<=maxlag)?0:1;
      }
      */
      int nlag = 4+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=1)?ilag:ilag-2-maxlag;
        lag2[ilag] = (ilag<=1)?0:1;
      }
      float[][] w = new float[3][3];
      float[][] v = new float[3][3];
      float[][] t = new float[3][3];
      float[][] r = new float[3][3];
      t[1][1] = 1.0f;
      LocalDipFilter ldf = new LocalDipFilter(Factor.NOT);
      CausalFilter cf = new CausalFilter(lag1,lag2);
      for (int itheta=0; itheta<NTHETA; ++itheta) {
        float theta = FTHETA+itheta*DTHETA;
        for (int iwidth=0; iwidth<NWIDTH; ++iwidth) {
          float width = FWIDTH+iwidth*DWIDTH;
          ArrayMath.fill(width,w);
          ArrayMath.fill(-sin(theta),v);
          ldf.applyForwardNot(w,v,t,r);
          cf.factorWilsonBurg(100,0.000001f,r);
          _atable[itheta][iwidth] = cf.getA();
        }
      }
      _lcf = new LocalCausalFilter(lag1,lag2);
    }
    void applyForward(float[][] vw, float[][] v1, float[][] x, float[][] y) {
      A2 a2 = new A2(_atable,vw,v1);
      _lcf.apply(a2,x,y);
      _lcf.applyTranspose(a2,y,y);
    }
    void applyInverse(float[][] vw, float[][] v1, float[][] x, float[][] y) {
      A2 a2 = new A2(_atable,vw,v1);
      _lcf.applyInverseTranspose(a2,x,y);
      _lcf.applyInverse(a2,y,y);
    }
    private static final float WIDTH_MIN = 0.01f;
    private static final float WIDTH_MAX = 1.00f;
    private static final int NWIDTH = 11;
    private static final float FWIDTH = WIDTH_MIN;
    private static final float DWIDTH = (WIDTH_MAX-WIDTH_MIN)/(float)(NWIDTH-1);
    private static final float SWIDTH = 0.9999f/DWIDTH;
    private static final int NTHETA = 21;
    private static final float FTHETA = -0.5f*FLT_PI;
    private static final float DTHETA = FLT_PI/(float)(NTHETA-1);;
    private static final float STHETA = 0.9999f/DTHETA;
    private static class A2 implements LocalCausalFilter.A2 {
      A2(float[][][] atable, float[][] vw, float[][] v1) {
        _at = atable;
        _vw = vw;
        _v1 = v1;
      }
      public void get(int i1, int i2, float[] a) {
        float theta = -asin(_v1[i2][i1]);
        float sigma = min(max(_vw[i2][i1],WIDTH_MIN),WIDTH_MAX);
        float s = (sigma-FWIDTH)*SWIDTH;
        float t = (theta-FTHETA)*STHETA;
        int is = (int)s;
        int it = (int)t;
        float s1 = s-(float)is;
        float t1 = t-(float)it;
        float s0 = 1.0f-s1;
        float t0 = 1.0f-t1;
        float[] a00 = _at[it  ][is  ];
        float[] a01 = _at[it  ][is+1];
        float[] a10 = _at[it+1][is  ];
        float[] a11 = _at[it+1][is+1];
        int n = a00.length;
        for (int j=0; j<n; ++j)
          a[j] = t0*(s0*a00[j]+s1*a01[j])+t1*(s0*a10[j]+s1*a11[j]);
      }
      private float[][][] _at;
      private float[][] _vw;
      private float[][] _v1;
    }
    private LocalCausalFilter _lcf;
    private float[][][] _atable = new float[NTHETA][NWIDTH][];
  }

  // A local dip filter approximated with minimum-phase factors.
  private static class FactoredFilter {
    FactoredFilter(Small small) {
      LocalDipFilter ldf = new LocalDipFilter(Factor.NOT);
      int maxlag = 4;
      /*
      int nlag = maxlag+2+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
        lag2[ilag] = (ilag<=maxlag)?0:1;
      }
      */
      int nlag = 4+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=1)?ilag:ilag-2-maxlag;
        lag2[ilag] = (ilag<=1)?0:1;
      }
      float[][] u = new float[3][3];
      float[][] t = new float[3][3];
      float[][] r = new float[3][3];
      t[1][1] = 1.0f;
      CausalFilter cf = new CausalFilter(lag1,lag2);
      for (int itheta=0; itheta<NTHETA; ++itheta) {
        float theta = FTHETA+itheta*DTHETA;
        trace("theta = "+theta*180.0f/FLT_PI);
        ArrayMath.fill(-sin(theta),u);
        ldf.applyForward(small.sd,small.sn,u,t,r);
        //ArrayMath.dump(r);
        cf.factorWilsonBurg(100,0.000001f,r);
        _atable[itheta] = cf.getA();
        ArrayMath.dump(cf.getA());
      }
      _lcf = new LocalCausalFilter(lag1,lag2);
    }
    void applyForward(float[][] u2, float[][] x, float[][] y) {
      A2 a2 = new A2(_atable,u2);
      _lcf.apply(a2,x,y);
      _lcf.applyTranspose(a2,y,y);
    }
    void applyInverse(float[][] u2, float[][] x, float[][] y) {
      A2 a2 = new A2(_atable,u2);
      _lcf.applyInverseTranspose(a2,x,y);
      _lcf.applyInverse(a2,y,y);
    }
    private static int NTHETA = 13;
    private static float FTHETA = -0.5f*FLT_PI;
    private static float DTHETA = FLT_PI/(float)(NTHETA-1);;
    private static float STHETA = 0.9999f/DTHETA;
    private static class A2 implements LocalCausalFilter.A2 {
      A2(float[][] atable, float[][] u2) {
        _at = atable;
        _u2 = u2;
      }
      public void xget(int i1, int i2, float[] a) {
        float theta = -asin(_u2[i2][i1]);
        int i = (int)(0.5f+(theta-FTHETA)*STHETA);
        float[] ai = _at[i];
        int n = ai.length;
        for (int j=0; j<n; ++j)
          a[j] = ai[j];
      }
      public void get(int i1, int i2, float[] a) {
        float theta = -asin(_u2[i2][i1]);
        float t = (theta-FTHETA)*STHETA;
        int i = (int)t;
        float w = t-(float)i;
        float[] ai = _at[i];
        float[] aj = _at[i+1];
        int n = ai.length;
        for (int j=0; j<n; ++j)
          a[j] = (1.0f-w)*ai[j]+w*aj[j];
      }
      private float[][] _at;
      private float[][] _u2;
    }
    private LocalCausalFilter _lcf;
    private float[][] _atable = new float[NTHETA][];
  }

  private static float dot(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    double s = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        s += x2[i1]*y2[i1];
    }
    return (float)s;
  }
  private static float dot(float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    double s = 0.0;
    for (int i3=0; i3<n3; ++i3) {
      float[][] x3 = x[i3];
      float[][] y3 = y[i3];
      for (int i2=0; i2<n2; ++i2) {
        float[] x32 = x3[i2];
        float[] y32 = y3[i2];
        for (int i1=0; i1<n1; ++i1)
          s += x32[i1]*y32[i1];
      }
    }
    return (float)s;
  }

  private static void saxpy(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        y2[i1] += a*x2[i1];
    }
  }
  private static void saxpy(float a, float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      float[][] x3 = x[i3];
      float[][] y3 = y[i3];
      for (int i2=0; i2<n2; ++i2) {
        float[] x32 = x3[i2];
        float[] y32 = y3[i2];
        for (int i1=0; i1<n1; ++i1)
          y32[i1] += a*x32[i1];
      }
    }
  }
}
