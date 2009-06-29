/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.CausalFilter;
import edu.mines.jtk.dsp.LocalCausalFilter;
import edu.mines.jtk.la.TridiagonalFMatrix;
import edu.mines.jtk.util.ArrayMath;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local plane filtering of images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.02.10
 */
public class LocalPlaneFilter {

  /**
   * Filter type.
   */
  public enum Type {
    CLAERBOUT1, // Jon Claerbout's "wavekill" filter A
    FOMEL1, // Fomel's 2002 filter with coefficients a function of slope
    FOMEL2, // Fomel's filter modified with coefficients a function of angle
    HALE1, // simplest finite-difference operators
    HALE2, // folded filter B
    HALE3, // SPD A'A with no preconditioner
    HALE4, // SPD A'A with minimum-phase factors
    HALE5, // SPD A'A with SSOR preconditioner
    HALE6, // SPD SSOR A'A
    HALE7, // folded filter B with conjugate gradients
  };

  public LocalPlaneFilter(Type type) {
    this(type,0.0);
  }

  public LocalPlaneFilter(Type type, double small) {
    _stability = (float)(1.0+small);
    _type = type;
    if (type==Type.CLAERBOUT1) {
      _filter6 = new Claerbout1Filter();
    } else if (type==Type.FOMEL1) {
      _filter6 = new Fomel1Filter();
    } else if (type==Type.FOMEL2) {
      _filter6 = new Fomel2Filter();
    } else if (type==Type.HALE1) {
      _filter6 = new Hale1Filter();
    } else if (type==Type.HALE2) {
      _filter6 = new Hale2Filter();
    }
  }

  /**
   * Applies this local plane filter.
   * Input and output arrays must be distinct.
   * @param u1 array of 1st components of normal vectors.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image.
   * @param y array with output image.
   */
  public void applyForward(
    float[][] u1, float[][] u2, 
    float[][] x, float[][] y) 
  {
    if (_filter6!=null) {
      applyForward(_filter6,u1,u2,x,y);
    } else if (_type==Type.HALE3) {
      applyForwardSpd(u1,u2,x,y);
    } else if (_type==Type.HALE4) {
      SpdMpFilter.applyForward(_stability,u1,u2,x,y);
    } else if (_type==Type.HALE5) {
      applyForwardSpd(u1,u2,x,y);
    } else if (_type==Type.HALE6) {
      applyForwardSsor(u1,u2,x,y);
    } else if (_type==Type.HALE7) {
      applyForwardFolded(u1,u2,x,y);
    }
  }

  /**
   * Applies the inverse of this local plane filter.
   * Input and output arrays must be distinct.
   * @param u1 array of 1st components of normal vectors.
   * @param u2 array of 2nd components of normal vectors.
   * @param x array with input image
   * @param y array with output image
   */
  public void applyInverse(
    float[][] u1, float[][] u2, 
    float[][] x, float[][] y) 
  {
    if (_filter6!=null) {
      applyInverse(_filter6,u1,u2,x,y);
    } else if (_type==Type.HALE3) {
      applyInverseSpd(u1,u2,x,y);
    } else if (_type==Type.HALE4) {
      SpdMpFilter.applyInverse(_stability,u1,u2,x,y);
    } else if (_type==Type.HALE5) {
      applyInverseSpdSsor(u1,u2,x,y);
    } else if (_type==Type.HALE6) {
      applyInverseSsor(u1,u2,x,y);
    } else if (_type==Type.HALE7) {
      applyInverseFolded(u1,u2,x,y);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Type _type; // filter type
  private float _stability; // 1 plus a small number
  private Filter6 _filter6; // for 6-coefficient filters

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  /**
   * Filters defined by six coefficients for each image sample.
   * Different plane filters have different coefficients.
   * For forward filtering, the coefficients are used as follows:
   * <pre><code>
   * y[i2][i1] = c[0]*x[i2  ][i1+1]+c[1]*x[i2  ][i1]+c[2]*x[i2  ][i1-1]
   *           + c[3]*x[i2-1][i1+1]+c[4]*x[i2-1][i1]+c[5]*x[i2-1][i1-1]
   * </code></pre>
   * For inverse filtering, a corresponding tridiagonal system of 
   * equations is solved. Note that some filters may not be invertible. 
   * The corresponding finite-difference stencil is
   * <pre><code>
   *        i2-1  i2
   *   i1-1 c[5] c[2]
   *   i1   c[4] c[1]
   *   i1+1 c[3] c[0]
   * </code></pre>
   */
  private interface Filter6 {
    /**
     * Gets coefficients for specified normal unit-vector.
     * @param i1 sample index in 1st dimension.
     * @param i2 sample index in 2nd dimension.
     * @param u1 array of 1st components of normal vectors.
     * @param u2 array of 2nd components of normal vectors.
     * @param c array of computed coefficients.
     */
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c);
  }

  // Simple Claerbout-Nichols wavekill filter (only 4 non-zero coefficients)
  // This quarter-plane filter is not zero-phase and not invertible.
  private static class Claerbout1Filter implements Filter6 {
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c) 
    {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float um = (u1i-u2i);
      float up = (u1i+u2i);
      c[0] = 0.0f;
      c[1] = um;
      c[2] = up;
      c[3] = 0.0f;
      c[4] = -up;
      c[5] = -um;
    }
  }

  // Fomel's 2002 filter. Coefficients are infinite for vertical planes.
  // Coefficients are multipled by 2 for consistency with other filters. 
  // Not invertible for steeper planes.
  private static class Fomel1Filter implements Filter6 {
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c) 
    {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      if (u1i<U1MIN)
        u1i = U1MIN;
      float s = u2i/u1i;
      c[0] = 2.0f*(1.0f-s)*(2.0f-s)/12.0f;
      c[1] = 2.0f*(2.0f+s)*(2.0f-s)/6.0f;
      c[2] = 2.0f*(1.0f+s)*(2.0f+s)/12.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
    private static final float U1MIN = 0.001f;
  }

  // Fomel's 2002 filter with coefficients multiplied by u1.
  // This modification makes filter coefficients finite for all planes.
  // Coefficients are multipled by 2 for consistency with other filters. 
  // Not invertible for steeper planes.
  private static class Fomel2Filter implements Filter6 {
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c) 
    {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float t1i = 2.0f*u1i;
      c[0] = 2.0f*(u1i-u2i)*(t1i-u2i)/12.0f;
      c[1] = 2.0f*(t1i+u2i)*(t1i-u2i)/6.0f;
      c[2] = 2.0f*(u1i+u2i)*(t1i+u2i)/12.0f;
      c[3] = -c[2];
      c[4] = -c[1];
      c[5] = -c[0];
    }
  }

  // Simple half-plane four-coefficient filter with stable inverse.
  private static class Hale1Filter implements Filter6 {
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c) 
    {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      c[0] = -0.5f*u2i*(u1i+u2i);
      c[1] =  1.0f;
      c[2] =  0.5f*u2i*(u1i-u2i);
      c[3] =  0.0f;
      c[4] = -u1i*u1i;
      c[5] =  0.0f;
    }
  }

  // Half-plane filter obtained by symmetric folding of zero-phase filter.
  // Coefficients are simply computed for each sample.
  // Inverse is barely stable for constant coefficients.
  private static class Hale2Filter implements Filter6 {
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c) 
    {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float um = 0.5f*(u1i-u2i);
      float up = 0.5f*(u1i+u2i);
      c[0] = 2.0f*um*up;
      c[1] = 1.0f;
      c[2] = c[0];
      c[3] = -2.0f*up*up;
      c[4] = -4.0f*um*up;
      c[5] = -2.0f*um*um;
    }
  }

  // Half-plane filter obtained by symmetric folding of the impulse
  // response of the cascade of a quarter-plane wavekill filter and 
  // its transpose. Coefficients are tricky and costly to compute.
  // The inverse is barely stable for constant-coefficients.
  private static class Hale2bFilter implements Filter6 {
    public void getCoefficients(int i1, int i2, 
      float[][] u1, float[][] u2, float[] c) 
    {
      int n1 = u1[0].length;
      int n2 = u1.length;
      int j1 = min(i1+1,n1-1);
      int j2 = min(i2+1,n2-1);
      float u1ii = u1[i2][i1];
      float u1ij = u1[i2][j1];
      float u1ji = u1[j2][i1];
      float u1jj = u1[j2][j1];
      float u2ii = u2[i2][i1];
      float u2ij = u2[i2][j1];
      float u2ji = u2[j2][i1];
      float u2jj = u2[j2][j1];
      float umii = 0.5f*(u1ii-u2ii);
      float umij = 0.5f*(u1ij-u2ij);
      float umji = 0.5f*(u1ji-u2ji);
      float umjj = 0.5f*(u1jj-u2jj);
      float upii = 0.5f*(u1ii+u2ii);
      float upij = 0.5f*(u1ij+u2ij);
      float upji = 0.5f*(u1ji+u2ji);
      float upjj = 0.5f*(u1jj+u2jj);
      c[0] = umii*upii+umji*upji;
      c[1] = umii*umii+upij*upij+upji*upji+umjj*umjj;
      c[2] = umij*upij+umjj*upjj;
      c[3] = -upij*upij-upji*upji;
      c[4] = -umii*upii-umji*upji-umij*upij-umjj*upjj;
      c[5] = -umii*umii-umjj*umjj;
    }
  }

  // Applys a 6-coefficient filter.
  private void applyForward(Filter6 filter, 
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] c = new float[6];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    int i1 = 0;
    filter.getCoefficients(i1,i2,u1,u2,c);
    c[1] *= _stability;
    y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1];
    for (i1=1; i1<n1-1; ++i1) {
      filter.getCoefficients(i1,i2,u1,u2,c);
      c[1] *= _stability;
      y2[i1] = c[0]*x2[i1+1]+c[1]*x2[i1]+c[2]*x2[i1-1];
    }
    filter.getCoefficients(i1,i2,u1,u2,c);
    c[1] *= _stability;
    y2[i1] = c[1]*x2[i1]+c[2]*x2[i1-1];

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] x2m = x[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      i1 = 0;
      filter.getCoefficients(i1,i2,u1,u2,c);
      c[1] *= _stability;
      y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]
             + c[3]*x2m[i1+1]+c[4]*x2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        filter.getCoefficients(i1,i2,u1,u2,c);
        c[1] *= _stability;
        y2[i1] = c[0]*x2 [i1+1]+c[1]*x2 [i1]+c[2]*x2 [i1-1]
               + c[3]*x2m[i1+1]+c[4]*x2m[i1]+c[5]*x2m[i1-1];
      }
      filter.getCoefficients(i1,i2,u1,u2,c);
      c[1] *= _stability;
      y2[i1] = c[1]*x2 [i1]+c[2]*x2 [i1-1]
             + c[4]*x2m[i1]+c[5]*x2m[i1-1];
    }
  }

  // Applys inverse of a 6-coefficient filter.
  private void applyInverse(Filter6 filter,
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    TridiagonalFMatrix tm = new TridiagonalFMatrix(n1);
    float[] ta = tm.a();
    float[] tb = tm.b();
    float[] tc = tm.c();
    float[] r = new float[n1];
    float[] c = new float[6];

    // For i2=0, x[i2-1][i1] = y[i2-1][i1] = x[i2][-1] = x[i2][n1] = 0.
    int i2 = 0;
    float[] x2 = x[i2];
    float[] y2 = y[i2];
    for (int i1=0; i1<n1; ++i1) {
      filter.getCoefficients(i1,i2,u1,u2,c);
      c[1] *= _stability;
      r[i1] = x2[i1];
      ta[i1] = c[2];
      tb[i1] = c[1];
      tc[i1] = c[0];
    }
    tm.solve(r,y2);

    // For all i2>0, assume that x[i2][-1] = x[i2][n1] = 0.
    for (i2=1; i2<n2; ++i2) {
      float[] y2m = y[i2-1];
      x2 = x[i2];
      y2 = y[i2];
      int i1 = 0;
      filter.getCoefficients(i1,i2,u1,u2,c);
      c[1] *= _stability;
      tb[i1] = c[1];
      tc[i1] = c[0];
      r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1];
      for (i1=1; i1<n1-1; ++i1) {
        filter.getCoefficients(i1,i2,u1,u2,c);
        c[1] *= _stability;
        ta[i1] = c[2];
        tb[i1] = c[1];
        tc[i1] = c[0];
        r[i1] = x2[i1]-c[3]*y2m[i1+1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      }
      filter.getCoefficients(i1,i2,u1,u2,c);
      c[1] *= _stability;
      ta[i1] = c[2];
      tb[i1] = c[1];
      r[i1] = x2[i1]-c[4]*y2m[i1]-c[5]*y2m[i1-1];
      tm.solve(r,y2);
    }
  }

  // Directional Laplacian filter.
  private void applyForwardSpd(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float epsilon = _stability-1.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float a11 = 1.0f-u1i*u1i;
        float a12 =     -u1i*u2i;
        float a22 = 1.0f-u2i*u2i;
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
        y[i2][i1] = ya+epsilon*x00;
        if (i1>0) y[i2][i1-1] -= yb;
        if (i2>0) y[i2-1][i1] += yb;
        if (i2>0 && i1>0) y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Folded directional Laplacian filter.
  private void applyForwardFolded(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float epsilon = _stability-1.0f;
    float[] xi20 = new float[n1];
    float[] xi21 = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] xtmp = xi20;  xi20 = xi21;  xi21 = xtmp;
      float[] yi2 = y[i2];
      ArrayMath.copy(x[i2],xi20);
      float u1i = u1[i2][0];
      float u2i = u2[i2][0];
      float v11 = 1.0f-u1i*u1i;
      float v12 =     -u1i*u2i;
      float v22 = 1.0f-u2i*u2i;
      float x00 = xi20[0];
      float x10 = xi21[0];
      float x01 = 0.0f;
      float x11 = 0.0f;
      float xa = x00-x11;
      float xb = x01-x10;
      float x1 = 0.5f*(xa-xb);
      float x2 = 0.5f*(xa+xb);
      float y1 = v11*x1+v12*x2;
      float y2 = v12*x1+v22*x2;
      float ya = y1+y2;
      float yb = 0.0f;
      yi2[0] = ya+epsilon*x00;
      for (int i1=1; i1<n1; ++i1) {
        u1i = u1[i2][i1];
        u2i = u2[i2][i1];
        v11 = 1.0f-u1i*u1i;
        v12 =     -u1i*u2i;
        v22 = 1.0f-u2i*u2i;
        x00 = xi20[i1];
        x10 = xi21[i1];
        x01 = xi20[i1-1];
        x11 = xi21[i1-1];
        xa = x00-x11;
        xb = x01-x10;
        x1 = 0.5f*(xa-xb);
        x2 = 0.5f*(xa+xb);
        y1 = v11*x1+v12*x2;
        y2 = v12*x1+v22*x2;
        yi2[i1 ]   = y1+y2+epsilon*x00;
        yi2[i1-1] -= y1-y2;
      }
    }
  }
  private void applyForwardFoldedB0(
    float[] u1, float[] u2, float[] x, float[] y) 
  {
    int n1 = x.length;
    float epsilon = _stability-1.0f;
    float u1i = u1[0];
    float u2i = u2[0];
    float v11 = 1.0f-u1i*u1i;
    float v12 =     -u1i*u2i;
    float v22 = 1.0f-u2i*u2i;
    float x00 = x[0];
    float x01 = 0.0f;
    float x1 = 0.5f*(x00-x01);
    float x2 = 0.5f*(x00+x01);
    float y1 = v11*x1+v12*x2;
    float y2 = v12*x1+v22*x2;
    y[0]  = y1+y2+epsilon*x00;
    for (int i1=1; i1<n1; ++i1) {
      u1i = u1[i1];
      u2i = u2[i1];
      v11 = 1.0f-u1i*u1i;
      v12 =     -u1i*u2i;
      v22 = 1.0f-u2i*u2i;
      x00 = x[i1];
      x01 = x[i1-1];
      x1 = 0.5f*(x00-x01);
      x2 = 0.5f*(x00+x01);
      y1 = v11*x1+v12*x2;
      y2 = v12*x1+v22*x2;
      y[i1  ]  = y1+y2+epsilon*x00;
      y[i1-1] -= y1-y2;
    }
  }
  private void applyForwardFoldedB1(
    float[] u1, float[] u2, float[] x, float[] y) 
  {
    int n1 = x.length;
    float u1i = u1[0];
    float u2i = u2[0];
    float v11 = 1.0f-u1i*u1i;
    float v12 =     -u1i*u2i;
    float v22 = 1.0f-u2i*u2i;
    float x10 = x[0];
    float x11 = 0.0f;
    float x1 = 0.5f*(x10-x11);
    float x2 = 0.5f*(x10+x11);
    float y1 = v11*x1-v12*x2;
    float y2 = v12*x1-v22*x2;
    y[0]  = y1+y2;
    for (int i1=1; i1<n1; ++i1) {
      u1i = u1[i1];
      u2i = u2[i1];
      v11 = 1.0f-u1i*u1i;
      v12 =     -u1i*u2i;
      v22 = 1.0f-u2i*u2i;
      x10 = x[i1];
      x11 = x[i1-1];
      x1 = 0.5f*(x10-x11);
      x2 = 0.5f*(x10+x11);
      y1 = v11*x1-v12*x2;
      y2 = v12*x1-v22*x2;
      y[i1  ]  = y1+y2;
      y[i1-1] -= y1-y2;
    }
  }
  private void applyInverseFolded(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] r = new float[n1];
    applyInverseFoldedB0(u1[0],u2[0],x[0],y[0]);
    for (int i2=1; i2<n2; ++i2) {
      applyForwardFoldedB1(u1[i2],u2[i2],y[i2-1],r);
      ArrayMath.sub(x[i2],r,r);
      ArrayMath.copy(y[i2-1],y[i2]);
      applyInverseFoldedB0(u1[i2],u2[i2],r,y[i2]);
    }
  }
  private void applyInverseFoldedB0(
    float[] u1, float[] u2, float[] x, float[] y) 
  {
    int n1 = x.length;
    float[] r = new float[n1]; // r
    float[] s = new float[n1]; // d
    float[] t = new float[n1]; // q
    //ArrayMath.zero(y);
    //ArrayMath.copy(x,r);
    applyForwardFoldedB0(u1,u2,y,t);
    ArrayMath.sub(x,t,r);
    ArrayMath.copy(r,s);
    float rr = dot(r,r);
    float small = rr*1.0e-12f;
    int niter;
    for (niter=0; niter<200 && rr>small; ++niter) {
      applyForwardFoldedB0(u1,u2,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      float rrold = rr;
      rr = dot(r,r);
      float beta = rr/rrold;
      for (int i1=0; i1<n1; ++i1)
        s[i1] = r[i1]+beta*s[i1];
    }
    trace("small="+small+" rr="+rr+" niter="+niter);
  }
  private static float dot(float[] x, float[] y) {
    int n1 = x.length;
    double s = 0.0;
    for (int i1=0; i1<n1; ++i1)
      s += x[i1]*y[i1];
    return (float)s;
  }
  private static void saxpy(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += a*x[i1];
  }

  // Directional isotropic Laplacian filter.
  // The first-derivative approximation used here is consistent
  // with an O(h^2) isotropic approximation to the Laplacian.
  // NOTE: if we use this, must also change preconditioner.
  private void applyForwardSpdIsotropic(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    float r = 0.5f*(1.0f+sqrt(2.0f/3.0f));
    float s = 0.5f*(1.0f-sqrt(2.0f/3.0f));
    int n1 = x[0].length;
    int n2 = x.length;
    float epsilon = _stability-1.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float a11 = 1.0f-u1i*u1i;
        float a12 =     -u1i*u2i;
        float a22 = 1.0f-u2i*u2i;
        float x00 = x[i2][i1];
        float x01 = (i1>0)?x[i2][i1-1]:0.0f;
        float x10 = (i2>0)?x[i2-1][i1]:0.0f;
        float x11 = (i2>0 && i1>0)?x[i2-1][i1-1]:0.0f;
        float x1 = r*(x00-x01)+s*(x10-x11);
        float x2 = r*(x00-x10)+s*(x01-x11);
        float y1 = a11*x1+a12*x2;
        float y2 = a12*x1+a22*x2;
        y[i2][i1] = r*y1+r*y2+epsilon*x00;
        if (i1>0) y[i2][i1-1] -= r*y1-s*y2;
        if (i2>0) y[i2-1][i1  ] += s*y1-r*y2;
        if (i2>0 && i1>0) y[i2-1][i1-1] -= s*y1+s*y2;
      }
    }
  }

  // Inverse of symmetric positive-definite filter without pre-conditioning.
  private void applyInverseSpd(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    ArrayMath.zero(y);
    ArrayMath.copy(x,r);
    ArrayMath.copy(r,s);
    float rr = dot(r,r);
    float small = rr*0.00001f;
    trace("small="+small);
    int niter;
    for (niter=0; niter<200 && rr>small; ++niter) {
      applyForwardSpd(u1,u2,s,t);
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
    ArrayMath.zero(y);
    ArrayMath.copy(x,r);
    smf.applyInverse(r,s);
    float rr = dot(r,s);
    float small = rr*0.00001f;
    trace("small="+small);
    int niter;
    for (niter=0; niter<100 && rr>small; ++niter) {
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
  // Inverse of symmetric positive-definite filter with SSOR pre-conditioning.
  private void applyInverseSpdSsor(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1]; // r
    float[][] s = new float[n2][n1]; // d
    float[][] t = new float[n2][n1]; // q
    float[][] w = new float[n2][n1]; // s
    ArrayMath.zero(y);
    ArrayMath.copy(x,r);
    applyInverseSsor(u1,u2,r,s);
    float rr = dot(r,s);
    float small = rr*0.00001f;
    trace("small="+small);
    int niter;
    for (niter=0; niter<100 && rr>small; ++niter) {
      applyForwardSpd(u1,u2,s,t);
      float alpha = rr/dot(s,t);
      saxpy( alpha,s,y);
      saxpy(-alpha,t,r);
      applyInverseSsor(u1,u2,r,w);
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
  private void applyForwardSsor(
    final float[][] u1, final float[][] u2, float[][] x, float[][] y) 
  {
    int[] lag1 = {0, 1,-1, 0, 1};
    int[] lag2 = {0, 0, 1, 1, 1};
    LocalCausalFilter lcf = new LocalCausalFilter(lag1,lag2);
    LocalCausalFilter.A2 a2 = new LocalCausalFilter.A2() {
      public void get(int i1, int i2, float[] a) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float um = 0.5f*(u1i-u2i);
        float up = 0.5f*(u1i+u2i);
        a[0] = _stability;
        a[1] = 2.0f*um*up;
        a[2] = -up*up;
        a[3] = -a[1];
        a[4] = -um*um;
      }
    };
    lcf.apply(a2,x,y);
    lcf.applyTranspose(a2,y,y);
  }
  private void applyInverseSsor(
    final float[][] u1, final float[][] u2, float[][] x, float[][] y) 
  {
    int[] lag1 = {0, 1,-1, 0, 1};
    int[] lag2 = {0, 0, 1, 1, 1};
    LocalCausalFilter lcf = new LocalCausalFilter(lag1,lag2);
    LocalCausalFilter.A2 a2 = new LocalCausalFilter.A2() {
      public void get(int i1, int i2, float[] a) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float um = 0.5f*(u1i-u2i);
        float up = 0.5f*(u1i+u2i);
        a[0] = _stability;
        a[1] = 2.0f*um*up;
        a[2] = -up*up;
        a[3] = -a[1];
        a[4] = -um*um;
      }
    };
    lcf.applyInverseTranspose(a2,x,y);
    lcf.applyInverse(a2,y,y);
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
  private void applyForwardSpd(
    float[][][] u2, float[][][] u3,
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float epsilon = _stability-1.0f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u3i = u3[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u1i = sqrt(1.0f-u2i*u2i+u3i*u3i);
          float a11 = 1.0f-u1i*u1i;
          float a22 = 1.0f-u2i*u2i;
          float a33 = 1.0f-u3i*u3i;
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
          y[i3][i2][i1] = ya + epsilon*x[i3][i2][i1];
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
  private void applyInverseSpd(
    float[][][] u2, float[][][] u3,
    float[][][] x, float[][][] y) 
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
    float small = rr*0.00001f;
    trace("small="+small);
    int niter;
    for (niter=0; niter<200 && rr>small; ++niter) {
      applyForwardSpd(u2,u3,s,t);
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

  // Symmetric-positive definite filter. Applies the quarter-plane
  // wavekill filter followed by the transpose of that filter. Care
  // is taken to get the transpose correct for variable coefficients.
  private void applyForwardSpdOld(
    float[][] u1, float[][] u2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float u1i,u2i,um,up,umy,upy;
    ArrayMath.copy(x,y);

    // Apply plane filter.
    for (int i2=n2-1; i2>0; --i2) {
      float[] y2 = y[i2];
      float[] y2m = y[i2-1];
      for (int i1=n1-1; i1>0; --i1) {
        u1i = u1[i2][i1]; u2i = u2[i2][i1]; 
        um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
        y2[i1] = um*(y2[i1]-y2m[i1-1])+up*(y2[i1-1]-y2m[i1]);
      }
      u1i = u1[i2][0]; u2i = u2[i2][0]; 
      um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
      y2[0] = um*y2[0]-up*y2m[0];
    }
    for (int i1=n1-1; i1>0; --i1) {
      u1i = u1[0][i1]; u2i = u2[0][i1]; 
      um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
      y[0][i1] = um*y[0][i1]+up*y[0][i1-1];
    }
    u1i = u1[0][0]; u2i = u2[0][0]; 
    um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
    y[0][0] = um*y[0][0];

    // Apply transpose of plane filter (with care for variable coefficients).
    u1i = u1[0][0]; u2i = u2[0][0]; 
    um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
    umy = um*y[0][0]; upy = up*y[0][0];
    y[0][0] = umy;
    for (int i1=1; i1<n1; ++i1) {
      u1i = u1[0][i1]; u2i = u2[0][i1]; 
      um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
      umy = um*y[0][i1]; upy = up*y[0][i1];
      y[0][i1] = umy;
      y[0][i1-1] += upy;
    }
    for (int i2=1; i2<n2; ++i2) {
      float[] y2 = y[i2];
      float[] y2m = y[i2-1];
      u1i = u1[i2][0]; u2i = u2[i2][0]; 
      um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
      umy = um*y2[0]; upy = up*y2[0];
      y2[0] = umy;
      y2m[0] -= upy;
      for (int i1=1; i1<n1; ++i1) {
        u1i = u1[i2][i1]; u2i = u2[i2][i1]; 
        um = 0.5f*(u1i-u2i); up = 0.5f*(u1i+u2i);
        umy = um*y2[i1]; upy = up*y2[i1];
        y2[i1] = umy;
        y2[i1-1] += upy;
        y2m[i1] -= upy;
        y2m[i1-1] -= umy;
      }
    }

    // Stability.
    float small = _stability-1.0f;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      for (int i1=0; i1<n1; ++i1)
        y2[i1] += small*x2[i1];
    }
  }

  // Symmetric positive-definite minimum-phase filter.
  private static class SpdMpFilter {
    SpdMpFilter(float stability, float[][] u1, float[][] u2) {
      _a2 = new A2(stability,u1,u2);
    }
    void applyForward(float[][] x, float[][] y) {
      _lcf.apply(_a2,x,y);
      _lcf.applyTranspose(_a2,y,y);
    }
    void applyInverse(float[][] x, float[][] y) {
      _lcf.applyInverseTranspose(_a2,x,y);
      _lcf.applyInverse(_a2,y,y);
    }
    static void applyForward(
      float stability, float[][] u1, float[][] u2, float[][] x, float[][] y) 
    {
      A2 a2 = new A2(stability,u1,u2);
      _lcf.apply(a2,x,y);
      _lcf.applyTranspose(a2,y,y);
    }
    static void applyInverse(
      float stability, float[][] u1, float[][] u2, float[][] x, float[][] y) 
    {
      A2 a2 = new A2(stability,u1,u2);
      _lcf.applyInverseTranspose(a2,x,y);
      _lcf.applyInverse(a2,y,y);
    }
    private static class A2 implements LocalCausalFilter.A2 {
      A2(float stability, float[][] u1, float[][] u2) {
        if (SAMPLE_THETA) {
          _index = indexTheta(u1,u2);
        } else {
          _index = indexU2(u1,u2);
        }
        _istab = 0;
        for (int istab=1; istab<NSTAB; ++istab) {
          if (abs(STABS[istab]-stability)<abs(STABS[_istab]-stability))
            _istab = istab;
        }
      }
      public void get(int i1, int i2, float[] a) {
        ArrayMath.copy(_aTable[_istab][_index[i2][i1]],a);
      }
      private int[][] _index;
      private int _istab;
    }
    private A2 _a2;
    private static boolean SAMPLE_THETA = true;
    private static int NSTAB = 3;
    private static float[] STABS = {1.1f,1.01f,1.001f};
    private static int NTHETA = 33;
    private static float FTHETA = -0.5f*FLT_PI;
    private static float DTHETA = FLT_PI/(float)(NTHETA-1);;
    private static float STHETA = 0.9999f/DTHETA;
    private static int NU2 = NTHETA;
    private static float FU2 = -1.0f;
    private static float DU2 = 2.0f/(float)(NU2-1);;
    private static float SU2 = 0.9999f/DU2;
    private static LocalCausalFilter _lcf;
    private static float[][][] _aTable = new float[NSTAB][NTHETA][];
    static {
      if (SAMPLE_THETA) {
        makeLcfTheta();
      } else {
        makeLcfU2();
      }
    }
    private static int[][] indexTheta(float[][] u1, float[][] u2) {
      int n1 = u1[0].length;
      int n2 = u1.length;
      int[][] i = new int[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float theta = asin(u2[i2][i1]);
          i[i2][i1] = (int)(0.5f+(theta-FTHETA)*STHETA);
        }
      }
      return i;
    }
    private static void makeLcfTheta() {
      // Stencil lag1:
      // lag2 = 0: 6  5  4  3  2  1  0
      // lag2 = 1:                1  0 -1 -2 -3 -4 -5
      int maxlag = 6;
      int nlag = maxlag+2+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
        lag2[ilag] = (ilag<=maxlag)?0:1;
      }
      for (int istab=0; istab<NSTAB; ++istab) {
        float stability = STABS[istab];
        for (int itheta=0; itheta<NTHETA; ++itheta) {
          float theta = FTHETA+itheta*DTHETA;
          float c = cos(theta);
          float s = sin(theta);
          float m12 = 0.5f*(c-s);
          float p12 = 0.5f*(c+s);
          float[][] r = {
            {    -m12*m12,   -2.0f*m12*p12,      -p12*p12},
            {2.0f*m12*p12,       stability,  2.0f*m12*p12},
            {    -p12*p12,   -2.0f*m12*p12,      -m12*m12}
          };
          CausalFilter cf = new CausalFilter(lag1,lag2);
          cf.factorWilsonBurg(100,0.000001f,r);
          _aTable[istab][itheta] = cf.getA();
        }
      }
      _lcf = new LocalCausalFilter(lag1,lag2);
    }
    private static int[][] indexU2(float[][] u1, float[][] u2) {
      int n1 = u1[0].length;
      int n2 = u1.length;
      int[][] i = new int[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u2i = u2[i2][i1];
          i[i2][i1] = (int)(0.5f+(u2i-FU2)*SU2);
        }
      }
      return i;
    }
    private static void makeLcfU2() {
      int maxlag = 6;
      int nlag = maxlag+2+maxlag;
      int[] lag1 = new int[nlag];
      int[] lag2 = new int[nlag];
      for (int ilag=0; ilag<nlag; ++ilag) {
        lag1[ilag] = (ilag<=maxlag)?ilag:ilag-2*maxlag;
        lag2[ilag] = (ilag<=maxlag)?0:1;
      }
      for (int istab=0; istab<NSTAB; ++istab) {
        float stability = STABS[istab];
        for (int iu2=0; iu2<NU2; ++iu2) {
          float u2 = FU2+iu2*DU2;
          float u1 = sqrt(1.0f-u2*u2);
          float m12 = 0.5f*(u1-u2);
          float p12 = 0.5f*(u1+u2);
          float[][] r = {
            {    -m12*m12,   -2.0f*m12*p12,      -p12*p12},
            {2.0f*m12*p12,       stability,  2.0f*m12*p12},
            {    -p12*p12,   -2.0f*m12*p12,      -m12*m12}
          };
          CausalFilter cf = new CausalFilter(lag1,lag2);
          cf.factorWilsonBurg(100,0.000001f,r);
          _aTable[istab][iu2] = cf.getA();
        }
      }
      _lcf = new LocalCausalFilter(lag1,lag2);
    }
  }
}
