/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local dip filtering of images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.03.20
 */
public class LocalDipFilter {

  /**
   * Filter type.
   */
  public enum Type {
    SIMPLE,
    NOTCH,
    DIP,
  };

  public LocalDipFilter(Type type) {
    this(type,0.01);
  }

  public LocalDipFilter(Type type, double small) {
    _small = (float)small;
    _tiny = 0.01f*_small;
    _type = type;
  }

  /**
   * Applies this local dip filter.
   * Input and output arrays must be distinct.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyForward(float[][] u2, float[][] x, float[][] y) {
    if (_type==Type.SIMPLE) {
      applyForward(_small,u2,x,y);
    } else {
      applyForward(_tiny,u2,x,y);
      applyInverse(_small,u2,x,y);
    }
  }

  /**
   * Applies the inverse of this local dip filter.
   * Input and output arrays must be distinct.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image
   * @param y array with output image
   */
  public void applyInverse(
    float[][] u1, float[][] u2, 
    float[][] x, float[][] y) 
  {
    if (_type==Type.SIMPLE) {
      applyInverse(_small,u2,x,y);
    } else {
      applyForward(_small,u2,x,y);
      applyInverse(_tiny,u2,x,y);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Type _type; // filter type
  private float _small; // small non-negative number
  private float _tiny; // really small non-negative number

  private static final boolean TRACE = false;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  // Directional Laplacian filter.
  private void applyForward(
    float small, float[][] u2, float[][] x, float[][] y) 
  {
    float aone = (_type==Type.DIP)?1.0f+small:1.0f;
    float anotch = (_type==Type.NOTCH)?small:0.0f;
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float a11 = aone-u1i*u1i;
        float a12 =     -u1i*u2i;
        float a22 = aone-u2i*u2i;
        float x00 = x[i2][i1];
        float x01 = (i1>0)?x[i2][i1-1]:0.0f;
        float x10 = (i2>0)?x[i2-1][i1]:0.0f;
        float x11 = (i2>0 && i1>0)?x[i2-1][i1-1]:0.0f;
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.5f*(xa-xb);
        float x2 = 0.5f*(xa+xb);
        float y1 = a11*x1+a12*x2;
        float y2 = a12*x1+a22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2][i1] = ya+anotch*x[i2][i1];
        if (i1>0) y[i2][i1-1] -= yb;
        if (i2>0) y[i2-1][i1] += yb;
        if (i2>0 && i1>0) y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Directional isotropic Laplacian filter.
  // The first-derivative approximation used here is consistent
  // with an O(h^2) isotropic approximation to the Laplacian.
  // NOTE: if we use this, must also change preconditioner.
  private void applyForwardIsotropic(
    float small, float[][] u2, float[][] x, float[][] y) 
  {
    float aone = (_type==Type.DIP)?1.0f+small:1.0f;
    float anotch = (_type==Type.NOTCH)?small:0.0f;
    float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u2i = u2[i2][i1];
        float u1i = sqrt(1.0f-u2i*u2i);
        float a11 = aone-u1i*u1i;
        float a12 =     -u1i*u2i;
        float a22 = aone-u2i*u2i;
        float x00 = x[i2][i1];
        float x01 = (i1>0)?x[i2][i1-1]:0.0f;
        float x10 = (i2>0)?x[i2-1][i1]:0.0f;
        float x11 = (i2>0 && i1>0)?x[i2-1][i1-1]:0.0f;
        float x1 = r*(x00-x01)+s*(x10-x11);
        float x2 = r*(x00-x10)+s*(x01-x11);
        float y1 = a11*x1+a12*x2;
        float y2 = a12*x1+a22*x2;
        y[i2][i1] = r*(y1+y2)+anotch*x[i2][i1];
        if (i1>0) y[i2][i1-1] -= r*y1-s*y2;
        if (i2>0) y[i2-1][i1  ] += s*y1-r*y2;
        if (i2>0 && i1>0) y[i2-1][i1-1] -= s*(y1+y2);
      }
    }
  }

  // Inverse of symmetric positive-definite filter without pre-conditioning.
  private void applyInverse(
    float small, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    Array.zero(y);
    Array.copy(x,r);
    Array.copy(r,s);
    float rr = dot(r,r);
    float stop = rr*0.00001f;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForward(small,u2,s,t);
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
      trace("niter="+niter+" rr="+rr);
    }
  }

  // Inverse of symmetric positive-definite filter with pre-conditioning.
  /*
  private void applyInverseSpdPc(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    SpdMpFilter smf = new SpdMpFilter(_stability,u1,u2);
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    float[][] w = new float[n2][n1]; // s
    Array.zero(y);
    Array.copy(x,r);
    smf.applyInverse(r,s);
    float rr = dot(r,s);
    float stop = rr*0.00001f;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<100 && rr>stop; ++niter) {
      applyForwardSpd(u1,u2,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      smf.applyInverse(r,w);
      float rrold = rr;
      rr = dot(r,w);
      float beta = rr/rrold;
      for (int i2=0; i2<n2; ++i2) {
        float[] w2 = w[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = w2[i1]+beta*s2[i1];
      }
      trace("niter="+niter+" rr="+rr);
    }
  }
  */
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

  // Directional Laplacian filter.
  private void applyForward(
    float small, float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
  {
    float aone = (_type==Type.DIP)?1.0f+small:1.0f;
    float anotch = (_type==Type.NOTCH)?small:0.0f;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u3i = u3[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u1i = sqrt(1.0f-u2i*u2i+u3i*u3i);
          float a11 = aone-u1i*u1i;
          float a22 = aone-u2i*u2i;
          float a33 = aone-u3i*u3i;
          float a12 =     -u1i*u2i;
          float a13 =     -u1i*u3i;
          float a23 =     -u2i*u3i;
          float x000 = x[i3][i2][i1];
          float x001 = (i1>0)?x[i3][i2][i1-1]:0.0f;
          float x010 = (i2>0)?x[i3][i2-1][i1]:0.0f;
          float x100 = (i3>0)?x[i3-1][i2][i1]:0.0f;
          float x011 = (i2>0 && i1>0)?x[i3][i2-1][i1-1]:0.0f;
          float x101 = (i3>0 && i1>0)?x[i3-1][i2][i1-1]:0.0f;
          float x110 = (i3>0 && i2>0)?x[i3-1][i2-1][i1]:0.0f;
          float x111 = (i3>0 && i2>0 && i1>0)?x[i3-1][i2-1][i1-1]:0.0f;
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
          y[i3][i2][i1] = ya + anotch*x[i3][i2][i1];
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

  // Inverse of symmetric positive-definite filter without pre-conditioning.
  private void applyInverse(
    float small, float[][][] u2, float[][][] u3, float[][][] x, float[][][] y) 
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
    float stop = rr*0.00001f;
    trace("stop="+stop);
    int niter;
    for (niter=0; niter<200 && rr>stop; ++niter) {
      applyForward(small,u2,u3,s,t);
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

  // A local dip filter approximated with minimum-phase factors.
  private static class FactoredFilter {
    FactoredFilter(LocalDipFilter ldf) {
      makeLcfTheta(ldf);
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
    private static int NTHETA = 33;
    private static float FTHETA = -0.5f*FLT_PI;
    private static float DTHETA = FLT_PI/(float)(NTHETA-1);;
    private static float STHETA = 0.9999f/DTHETA;
    private static class A2 implements LocalCausalFilter.A2 {
      A2(float[][] atable, float[][] u2) {
        _at = atable;
        _u2 = u2;
      }
      public void get(int i1, int i2, float[] a) {
        float theta = -asin(_u2[i2][i1]);
        int i = (int)(0.5f+(theta-FTHETA)*STHETA);
        float[] ai = _at[i];
        int n = ai.length;
        for (int j=0; j<n; ++j)
          a[j] = ai[j];
      }
      private float[][] _at;
      private float[][] _u2;
    }
    private float[][] _atable = new float[NTHETA][];
    private LocalCausalFilter _lcf;
    private void makeLcfTheta(LocalDipFilter ldf) {
      int maxlag = 6;
      int nlag = maxlag+2+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
        lag2[ilag] = (ilag<=maxlag)?0:1;
      }
      float[][] u = new float[3][3];
      float[][] t = new float[3][3];
      float[][] r = new float[3][3];
      t[1][1] = 1.0f;
      CausalFilter cf = new CausalFilter(lag1,lag2);
      for (int itheta=0; itheta<NTHETA; ++itheta) {
        float theta = FTHETA+itheta*DTHETA;
        Array.fill(-sin(theta),u);
        ldf.applyForward(u,t,r);
        cf.factorWilsonBurg(100,0.000001f,r);
        _atable[itheta] = cf.getA();
      }
      _lcf = new LocalCausalFilter(lag1,lag2);
    }
  }
}
