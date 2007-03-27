/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local line filtering of images.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.03.25
 */
public class LocalLineFilter {

  /**
   * Use of minimum-phase filter factorization.
   * <em>Factorization is not yet implemented for 3-D filters.</em>
   * <dl>
   * <dt>NOT<dd>
   * not used
   * <dt>PCG<dd>
   * used to precondition inverse filter by conjugate gradients
   * <dt>INV<dd>
   * used for inverse filter only
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
   * Constructs a local line filter that does not use factorization. 
   */
  public LocalLineFilter() {
    this(Factor.NOT);
  }

  /**
   * Constructs a local line filter with specified factorization.
   * @param factor filter factorization.
   */
  public LocalLineFilter(Factor factor) {
    //_factor = factor;
  }

  /**
   * Applies this local line filter.
   * Input and output arrays must be distinct.
   * @param sl small parameter that controls width of line filter.
   * @param sn small parameter that controls width of notch filter.
   * @param w1 array of 1st components of line vectors.
   * @param w2 array of 2nd components of line vectors.
   * @param w3 array of 3rd components of line vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyForward(
    float sl, float sn, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    applyForwardNot(sl,sn,w1,w2,w3,x,y);
  }
  public void applyForward(
    float sl, float sn, float[][][] e1,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    applyForwardNot(sl,sn,e1,w1,w2,w3,x,y);
  }

  /**
   * Applies the inverse of this local line filter.
   * Input and output arrays must be distinct.
   * @param sl small parameter that controls width of line filter.
   * @param sn small parameter that controls width of notch filter.
   * @param w1 array of 1st components of line vectors.
   * @param w2 array of 2nd components of line vectors.
   * @param w3 array of 3rd components of line vectors.
   * @param x array with input image
   * @param y array with output image
   */
  public void applyInverse(
    float sl, float sn, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    applyInverseNot(sl,sn,w1,w2,w3,x,y);
  }
  public void applyInverse(
    float sl, float sn, float[][][] e1,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    applyInverseNot(sl,sn,e1,w1,w2,w3,x,y);
  }

  /**
   * Applies a line filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param sl small parameter that controls width of line filter.
   *  This parameter is used for the inverse filter only.
   * @param w1 array of 1st components of line vectors.
   * @param w2 array of 2nd components of line vectors.
   * @param w3 array of 3rd components of line vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyLine(
    float sl,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    float[][][] t = new float[x.length][x[0].length][x[0][0].length];
    applyForward(0.0f,0.0f,w1,w2,w3,x,t);
    applyInverse(  sl,0.0f,w1,w2,w3,t,y);
  }
  public void applyLine(
    float sl, float[][][] e1,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    float[][][] t = new float[x.length][x[0].length][x[0][0].length];
    applyForward(0.0f,0.0f,e1,w1,w2,w3,x,t);
    applyInverse(  sl,0.0f,e1,w1,w2,w3,t,y);
  }

  /**
   * Applies a notch filter comprised of forward and inverse filters.
   * Input and output arrays must be distinct.
   * @param sn small parameter that controls width of notch filter.
   *  This parameter is used for the inverse filter only.
   * @param w1 array of 1st components of line vectors.
   * @param w2 array of 2nd components of line vectors.
   * @param w3 array of 3rd components of line vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyNotch(
    float sn,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    float[][][] t = new float[x.length][x[0].length][x[0][0].length];
    applyForward(0.0f,0.0f,w1,w2,w3,x,t);
    applyInverse(0.0f,  sn,w1,w2,w3,t,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float CG_SMALL = 0.000001f;

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  // Forward (not factored) filter.
  private void applyForwardNot(
    float sl, float sn, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float w1i = w1[i3][i2][i1];
          float w2i = w2[i3][i2][i1];
          float w3i = w3[i3][i2][i1];
          float a11 = sl+w1i*w1i;
          float a22 = sl+w2i*w2i;
          float a33 = sl+w3i*w3i;
          float a12 =    w1i*w2i;
          float a13 =    w1i*w3i;
          float a23 =    w2i*w3i;
          float x000 = x[i3][i2][i1];
          float x001 = (i1>0)?x[i3][i2][i1-1]:0.0f;
          float x010 = (i2>0)?x[i3][i2-1][i1]:0.0f;
          float x100 = (i3>0)?x[i3-1][i2][i1]:0.0f;
          float x011 = (i2>0 && i1>0)?x[i3][i2-1][i1-1]:0.0f;
          float x101 = (i3>0 && i1>0)?x[i3-1][i2][i1-1]:0.0f;
          float x110 = (i3>0 && i2>0)?x[i3-1][i2-1][i1]:0.0f;
          float x111 = (i3>0 && i2>0 && i1>0)?x[i3-1][i2-1][i1-1]:0.0f;
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
          y[i3][i2][i1] = ya + sn*x[i3][i2][i1];
          if (i1>0) y[i3][i2][i1-1] -= yd;
          if (i2>0) y[i3][i2-1][i1] += yb;
          if (i3>0) y[i3-1][i2][i1] += yc;
          if (i2>0 && i1>0) y[i3][i2-1][i1-1] -= yc;
          if (i3>0 && i1>0) y[i3-1][i2][i1-1] -= yb;
          if (i3>0 && i2>0) y[i3-1][i2-1][i1] += yd;
          if (i3>0 && i2>0 && i1>0) y[i3-1][i2-1][i1-1] -= ya;
        }
      }
    }
  }
  private void applyForwardNot(
    float sl, float sn, float[][][] e1,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float w1i = w1[i3][i2][i1];
          float w2i = w2[i3][i2][i1];
          float w3i = w3[i3][i2][i1];
          float e1i = e1[i3][i2][i1];
          float sw  = e1i*e1i*e1i;
          float a11 = sl+sw*w1i*w1i;
          float a22 = sl+sw*w2i*w2i;
          float a33 = sl+sw*w3i*w3i;
          float a12 =    sw*w1i*w2i;
          float a13 =    sw*w1i*w3i;
          float a23 =    sw*w2i*w3i;
          float x000 = x[i3][i2][i1];
          float x001 = (i1>0)?x[i3][i2][i1-1]:0.0f;
          float x010 = (i2>0)?x[i3][i2-1][i1]:0.0f;
          float x100 = (i3>0)?x[i3-1][i2][i1]:0.0f;
          float x011 = (i2>0 && i1>0)?x[i3][i2-1][i1-1]:0.0f;
          float x101 = (i3>0 && i1>0)?x[i3-1][i2][i1-1]:0.0f;
          float x110 = (i3>0 && i2>0)?x[i3-1][i2-1][i1]:0.0f;
          float x111 = (i3>0 && i2>0 && i1>0)?x[i3-1][i2-1][i1-1]:0.0f;
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
          y[i3][i2][i1] = ya + sn*x[i3][i2][i1];
          if (i1>0) y[i3][i2][i1-1] -= yd;
          if (i2>0) y[i3][i2-1][i1] += yb;
          if (i3>0) y[i3-1][i2][i1] += yc;
          if (i2>0 && i1>0) y[i3][i2-1][i1-1] -= yc;
          if (i3>0 && i1>0) y[i3-1][i2][i1-1] -= yb;
          if (i3>0 && i2>0) y[i3-1][i2-1][i1] += yd;
          if (i3>0 && i2>0 && i1>0) y[i3-1][i2-1][i1-1] -= ya;
        }
      }
    }
  }

  // Inverse filter via conjugate gradients without preconditioning.
  private void applyInverseNot(
    float sl, float sn, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] r = new float[n3][n2][n1]; // r
    float[][][] s = new float[n3][n2][n1]; // d
    float[][][] t = new float[n3][n2][n1]; // q
    Array.zero(y);
    Array.copy(x,r);
    Array.copy(r,s);
    float rr = dot(r,r);
    float stop = rr*CG_SMALL;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForwardNot(sl,sn,w1,w2,w3,s,t);
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
  private void applyInverseNot(
    float sl, float sn, float[][][] e1,
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] r = new float[n3][n2][n1]; // r
    float[][][] s = new float[n3][n2][n1]; // d
    float[][][] t = new float[n3][n2][n1]; // q
    Array.zero(y);
    Array.copy(x,r);
    Array.copy(r,s);
    float rr = dot(r,r);
    float stop = rr*CG_SMALL;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForwardNot(sl,sn,e1,w1,w2,w3,s,t);
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
