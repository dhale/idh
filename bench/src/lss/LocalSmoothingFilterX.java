/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lss;

import java.util.logging.Logger;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.AtomicFloat;
import edu.mines.jtk.util.Threads;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Local smoothing of images with tensor filter coefficients.
 * <p>
 * EXPERIMENTAL VERSION
 * <p>
 * Smoothing is performed by solving a sparse symmetric positive-definite
 * system of equations: (I+G'DG)y = x, where G is a matrix of gradient 
 * operators, D is a matrix of tensor filter coefficients, x is an input 
 * image, and y is an output image.
 * <p>
 * This sparse system of filter equations is solved iteratively, beginning
 * with y = x. Iterations continue until either the error in the solution 
 * y is below a specified threshold or the number of iterations exceeds a 
 * specified limit.
 * <p>
 * For low wavenumbers the output of this filter approximates the solution 
 * to an anisotropic inhomogeneous diffusion equation, where the filter 
 * input x corresponds to the initial condition at time t = 0 and filter 
 * output y corresponds to the solution at some later time t.
 * <p>
 * Derivatives in the gradient operator G are approximated by small
 * finite difference stencils, and for two (and higher) dimensions
 * the frequency responses of these approximations have zeros in the 
 * Nyquist corners. Therefore, for these high corner frequencies, 
 * (I+G'DG)y = y = x. In other words, if these frequencies are present 
 * in the input image x, then they will appear unattenuated in the 
 * output image y. In 2D output images y, these frequencies appear as 
 * a checkerboard pattern.
 * <p>
 * To attenuate these high frequencies, an optional smoothing filter 
 * may be applied to the input x, before solving for the output y. 
 * However, if specified, this optional smoothing is applied to all
 * input samples, even those for which the tensor field D may be zero.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.12.08
 */
public class LocalSmoothingFilterX {

  /**
   * Constructs a local smoothing filter with default parameters.
   * The default parameter small is 0.01 and the default maximum 
   * number of iterations is 100.
   */
  public LocalSmoothingFilterX() {
    this(0.01,100);
  }

  /**
   * Constructs a local smoothing filter with specified iteration parameters.
   * @param small stop when norm of residuals is less than this factor times
   *  the norm of the input array.
   * @param niter stop when number of iterations exceeds this limit.
   */
  public LocalSmoothingFilterX(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  /**
   * Applies this filter for specified constant scale factor.
   * Local smoothing for 1D arrays is a special case that requires no tensors. 
   * All tensors are implicitly scalar values equal to one, so that filtering 
   * is determined entirely by the specified constant scale factor.
   * @param c constant scale factor.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float c, float[] x, float[] y) {
    apply(c,null,x,y);
  }

  /**
   * Applies this filter for specified scale factors.
   * Local smoothing for 1D arrays is a special case that requires no tensors. 
   * All tensors are implicitly scalar values equal to one, so that filtering 
   * is determined entirely by the specified scale factors.
   * @param c constant scale factor.
   * @param s array of scale factors.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float c, float[] s, float[] x, float[] y) {
    int n1 = x.length;

    // Sub-diagonal e of SPD tridiagonal matrix I+G'DG; e[0] = e[n1] = 0.0.
    float[] e = new float[n1+1];
    if (s!=null) {
      c = -0.5f*c;
      for (int i1=1; i1<n1; ++i1)
        e[i1] = c*(s[i1]+s[i1-1]);
    } else {
      c = -c;
      for (int i1=1; i1<n1; ++i1)
        e[i1] = c;
    }

    // Work array w overwrites sub-diagonal array e.
    float[] w = e;

    // Solve tridiagonal system of equations (I+G'DG)y = x.
    float t = 1.0f-e[0]-e[1];
    y[0] = x[0]/t;
    for (int i1=1; i1<n1; ++i1) {
      float di = 1.0f-e[i1]-e[i1+1]; // diagonal element
      float ei = e[i1]; // sub-diagonal element
      w[i1] = ei/t;
      t = di-ei*w[i1];
      y[i1] = (x[i1]-ei*y[i1-1])/t;
    }
    for (int i1=n1-1; i1>0; --i1)
      y[i1-1] -= w[i1]*y[i1];
  }

  /**
   * Applies this filter for specified tensor coefficients.
   * @param d tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(Tensors2 d, float[][] x, float[][] y) 
  {
    apply(d,1.0f,x,y);
  }

  /**
   * Applies this filter for specified tensor coefficients and scale factor.
   * @param d tensor coefficients.
   * @param c constant scale factor for tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(Tensors2 d, float c, float[][] x, float[][] y) {
    apply(d,c,null,x,y);
  }

  /**
   * Applies this filter for specified tensor coefficients and scale factors.
   * @param d tensor coefficients.
   * @param c constant scale factor for tensor coefficients.
   * @param s array of scale factors for tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    Operator2 a = new LhsOperator2(d,c,s);
    //float[][] r = like(x);
    //applyRhs(x,r);
    //smooth(x,r);
    if (METHOD==2222) {
      float[][] r = like(x);
      smooth(x,r);
      x = r;
    }
    solve(a,x,y);
  }
  public void applyTranspose(Tensors2 d, float[][] x, float[][] y) {
    applyTranspose(d,1.0f,x,y);
  }
  public void applyTranspose(
    Tensors2 d, float c, float[][] x, float[][] y) 
  {
    applyTranspose(d,c,null,x,y);
  }
  public void applyTranspose(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    Operator2 a = new LhsOperator2(d,c,s);
    //float[][] r = like(x);
    solve(a,x,y);
    //smooth(r,y);
    //applyRhs(r,y);
  }
  public static float[][] like(float[][] x) {
    return new float[x.length][x[0].length];
  }

  /**
   * Applies this filter for specified tensor coefficients.
   * @param d tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(Tensors3 d, float[][][] x, float[][][] y) 
  {
    apply(d,1.0f,x,y);
  }

  /**
   * Applies this filter for specified tensor coefficients and scale factor.
   * @param d tensor coefficients.
   * @param c constant scale factor for tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(Tensors3 d, float c, float[][][] x, float[][][] y) {
    apply(d,c,null,x,y);
  }

  /**
   * Applies this filter for specified tensor coefficients and scale factors.
   * @param d tensor coefficients.
   * @param c constant scale factor for tensor coefficients.
   * @param s array of scale factors for tensor coefficients.
   * @param x input array.
   * @param y output array.
   */
  public void apply(
    Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    Operator3 a = new LhsOperator3(d,c,s);
    float[][][] r = applyRhs(x);
    solve(a,r,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static Logger log = 
    Logger.getLogger(LocalSmoothingFilter.class.getName());

  private static final boolean PARALLEL = true; // false for single-threaded
  private static final boolean SMOOTH = false; // false for I instead of S'S

  private float _small; // stop iterations when residuals are small
  private int _niter; // number of iterations

  public static int METHOD = 22; // smoothing method

  /**
   * A symmetric positive-definite operator.
   */
  private static interface Operator2 {
    public void apply(float[][] x, float[][] y);
  }
  private static interface Operator3 {
    public void apply(float[][][] x, float[][][] y);
  }

  private static class LhsOperator2 implements Operator2 {
    LhsOperator2(Tensors2 d, float c, float[][] s) {
      _d = d;
      _c = c;
      _s = s;
    }
    public void apply(float[][] x, float[][] y) {
      applyLhs(_d,_c,_s,x,y);
    }
    private Tensors2 _d;
    private float _c;
    private float[][] _s;
  }

  private static class LhsOperator3 implements Operator3 {
    LhsOperator3(Tensors3 d, float c, float[][][] s) {
      _d = d;
      _c = c;
      _s = s;
    }
    public void apply(float[][][] x, float[][][] y) {
      applyLhs(_d,_c,_s,x,y);
    }
    private Tensors3 _d;
    private float _c;
    private float[][][] _s;
  }

  /**
   * Computes y = S'Sx.
   */
  private static void applyRhsX(float[][] x, float[][] y) {
    if (SMOOTH) {
      smooth(x,y);
    } else {
      scopy(x,y);
    }
  }
  private static void smooth(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float a = 0.2500f;
    float b = 0.1250f;
    float c = 0.0625f;
    /*
    float p = 0.3f;
    float q = 1.0f+p*p;
    float s = 1.0f/pow(1.0f+p,4.0f);
    float a = s*q*q;
    float b = s*p*q;
    float c = s*p*p;
    */
    float[] xm = new float[1+n1+1];
    float[] x0 = new float[1+n1+1];
    float[] xp = new float[1+n1+1];
    float[] xt;
    copy(n1,0,x[0],1,xm); xm[0] = xm[1]; xm[n1+1] = xm[n1];
    copy(n1,0,x[0],1,x0); x0[0] = x0[1]; x0[n1+1] = x0[n1];
    for (int j2=0; j2<n2; ++j2,xt=xm,xm=x0,x0=xp,xp=xt) {
      int i2p = min(j2+1,n2-1);
      copy(n1,0,x[i2p],1,xp); xp[0] = xp[1]; xp[n1+1] = xp[n1];
      float xmm, xm0 = xm[0], xmp = xm[1];
      float x0m, x00 = x0[0], x0p = x0[1]; 
      float xpm, xp0 = xp[0], xpp = xp[1]; 
      for (int j1=0,i1p=2; j1<n1; ++j1,++i1p) {
        xmm = xm0; xm0 = xmp; xmp = xm[i1p]; 
        x0m = x00; x00 = x0p; x0p = x0[i1p]; 
        xpm = xp0; xp0 = xpp; xpp = xp[i1p]; 
        y[j2][j1] = a*x00+b*(x0m+xm0+xp0+x0p)+c*(xmm+xpm+xmp+xpp);
      }
    }
  }
  private static void smoothX(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float kmax = 1.0f, kmaxs = kmax*kmax;
    int npad = 29;
    int n1fft = FftReal.nfftSmall(n1+npad);
    int n2fft = FftComplex.nfftSmall(n2+npad);
    float[][] xfft = new float[n2fft][n1fft+2];
    float[][] hfft = new float[n2fft][n1fft+2];
    float j1 = n1fft/2;
    float j2 = n2fft/2;
    copy(n1,n2,x,xfft);
    KaiserWindow kw = KaiserWindow.fromErrorAndLength(0.005,npad);
    for (int i2=0; i2<n2fft; ++i2) {
      float x2 = i2-j2;
      for (int i1=0; i1<n1fft; ++i1) {
        float x1 = i1-j1;
        float ri = sqrt(x1*x1+x2*x2);
        float wi = (float)(kw.evaluate(x1)*kw.evaluate(x2));
        float hi = wi*kmaxs*h2(kmax*ri);
        hfft[i2][i1] = hi;
      }
    }
    FftReal fft1 = new FftReal(n1fft);
    FftComplex fft2 = new FftComplex(n2fft);
    int m1 = n1fft/2+1;
    fft1.realToComplex1(-1,n2,xfft,xfft);
    fft2.complexToComplex2(-1,m1,xfft,xfft);
    fft1.realToComplex1(-1,n2,hfft,hfft);
    fft2.complexToComplex2(-1,m1,hfft,hfft);
    for (int i2=0; i2<n2fft; ++i2) {
      for (int i1=0,ir=0,ii=1; i1<m1; ++i1,ir+=2,ii+=2) {
        float hr = hfft[i2][ir];
        float hi = hfft[i2][ii];
        float hh = sqrt(hr*hr+hi*hi);
        xfft[i2][ir] *= hh;
        xfft[i2][ii] *= hh;
      }
    }
    fft2.complexToComplex2(1,m1,xfft,xfft);
    fft1.complexToReal1(1,n2,xfft,xfft);
    copy(n1,n2,xfft,y);
    mul(1.0f/(float)n1fft/(float)n2fft,y,y);
  }

  private static float h2(float r) {
    if (r==0.0) {
      return FLT_PI/4.0f;
    } else {
      float pir = FLT_PI*r;
      return j1(pir)/(2.0f*r);
    }
  }
  private static float j1(float x) {
    return (float)j1((double)x);
  }
  private static double j1(double x) {
    double ax = abs(x);
    if (ax<8.0) {
      double xx = x*x;
      double num = x*(72362614232.0 + 
        xx*(-7895059235.0 +
        xx*(242396853.1 +
        xx*(-2972611.439 +
        xx*(15704.48260+
        xx*(-30.16036606))))));
      double den = 144725228442.0 + 
        xx*(2300535178.0 +
        xx*(18583304.74 +
        xx*(99447.43394 +
        xx*(376.9991397 +
        xx))));
      return num/den;
    } else {
      double z = 8.0/ax;
      double zz = z*z;
      double t1 = 1.0 + 
        zz*(0.183105e-2 +
        zz*(-0.3516396496e-4 +
        zz*(0.2457520174e-5 +
        zz*(-0.240337019e-6))));
      double t2 = 0.04687499995 + 
        zz*(-0.2002690873e-3 +
        zz*(0.8449199096e-5 +
        zz*(-0.88228987e-6 +
        zz*0.105787412e-6)));
      double am = ax-2.356194491;
      double y = sqrt(0.636619772/ax)*(cos(am)*t1-z*sin(am)*t2);
      return (x<0.0)?-y:y;
    }
  }

  /**
   * Returns y = S'Sx.
   */
  private static void applyRhs(float[][] x, float[][] y) {
    if (SMOOTH) {
      szero(y);
      int n1 = x[0].length;
      int n2 = x.length;
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float x00 = x[i2  ][i1  ];
          float x01 = x[i2  ][i1-1];
          float x10 = x[i2-1][i1  ];
          float x11 = x[i2-1][i1-1];
          //         0.0625 = 1/16
          float xs = 0.0625f*(x00+x01+x10+x11);
          y[i2  ][i1  ] += xs;
          y[i2  ][i1-1] += xs;
          y[i2-1][i1  ] += xs;
          y[i2-1][i1-1] += xs;
        }
      }
    } else {
      scopy(x,y);
    }
  }

  /**
   * Returns y = S'Sx.
   */
  private static float[][][] applyRhs(float[][][] x) {
    if (!SMOOTH)
      return x;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] y = new float[n3][n2][n1];
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        float[] x00 = x[i3  ][i2  ];
        float[] x01 = x[i3  ][i2-1];
        float[] x10 = x[i3-1][i2  ];
        float[] x11 = x[i3-1][i2-1];
        float[] y00 = y[i3  ][i2  ];
        float[] y01 = y[i3  ][i2-1];
        float[] y10 = y[i3-1][i2  ];
        float[] y11 = y[i3-1][i2-1];
        for (int i1=1; i1<n1; ++i1) {
          int i1m = i1-1;
          float x000 = x00[i1 ];
          float x001 = x00[i1m];
          float x010 = x01[i1 ];
          float x011 = x01[i1m];
          float x100 = x10[i1 ];
          float x101 = x10[i1m];
          float x110 = x11[i1 ];
          float x111 = x11[i1m];
          //         0.015625 = 1/64
          float xs = 0.015625f*(x000+x001+x010+x011+x100+x101+x110+x111);
          y00[i1 ] += xs;
          y00[i1m] += xs;
          y01[i1 ] += xs;
          y01[i1m] += xs;
          y10[i1 ] += xs;
          y10[i1m] += xs;
          y11[i1 ] += xs;
          y11[i1m] += xs;
        }
      }
    }
    return y;
  }


  public static void applyLhs(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    if (METHOD==2222) { // 2x2 with 2x2 smoothing
      applyLhs2222(d,c,s,x,y);
    } else if (METHOD==22) { // 2x2
      applyLhs22(d,c,s,x,y);
    } else if (METHOD==24) { // 2x4
      applyLhs24(d,c,s,x,y);
    } else if (METHOD==242) { // 2x4, with filter buffer
      applyLhs242(d,c,s,x,y);
    } else if (METHOD==33) { // Scharr's 3x3
      applyLhs33(d,c,s,x,y);
    } else if (METHOD==332) { // Scharr's 3x3, with filter buffer
      applyLhs332(d,c,s,x,y);
    } else if (METHOD==44) {
      applyLhs44(d,c,s,x,y);
    } else if (METHOD==51) {
      applyLhs51(d,c,s,x,y);
    } else if (METHOD==71) {
      applyLhs71(d,c,s,x,y);
    } else if (METHOD==712) {
      applyLhs712(d,c,s,x,y);
    } else if (METHOD==91) {
      applyLhs91(d,c,s,x,y);
    } else if (METHOD==131) {
      applyLhs131(d,c,s,x,y);
    }
  }

  private static void applyLhs2222(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    smooth(x,y);
    c *= 0.25f;
    float[] di = new float[3];
    for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
      float[] x0 = x[i2 ];
      float[] xm = x[i2m];
      float[] y0 = y[i2 ];
      float[] ym = y[i2m];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x00 = x0[i1 ];
        float x01 = x0[i1m];
        float x10 = xm[i1 ];
        float x11 = xm[i1m];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = xa-xb;
        float x2 = xa+xb;
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y0[i1 ] += ya;
        y0[i1m] -= yb;
        ym[i1 ] += yb;
        ym[i1m] -= ya;
      }
    }
  }

  private static void applyLhs22(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    scopy(x,y);
    c *= 0.25f;
    float[] di = new float[3];
    for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
      float[] x0 = x[i2 ];
      float[] xm = x[i2m];
      float[] y0 = y[i2 ];
      float[] ym = y[i2m];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x00 = x0[i1 ];
        float x01 = x0[i1m];
        float x10 = xm[i1 ];
        float x11 = xm[i1m];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = xa-xb;
        float x2 = xa+xb;
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y0[i1 ] += ya;
        y0[i1m] -= yb;
        ym[i1 ] += yb;
        ym[i1m] -= ya;
      }
    }
  }

  public static void applyLhs24(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    float p = 0.18f;
    float a = 0.5f*(1.0f+p);
    float b = 0.5f*(    -p);
    c *= a*a;
    b /= a;
    scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    int i2m2, i2m1 = 0, i2p0 = 0, i2p1 = 1;
    for (int i2=1; i2<n2; ++i2) {
      i2m2 = i2m1; i2m1 = i2p0; i2p0 = i2p1; ++i2p1; 
      if (i2p1>=n1) i2p1 = n1-1;
      float[] xm2=x[i2m2], xm1=x[i2m1], xp0=x[i2p0], xp1=x[i2p1];
      float[] ym2=y[i2m2], ym1=y[i2m1], yp0=y[i2p0], yp1=y[i2p1];
      int m2, m1 = 0, p0 = 0, p1 = 1;
      for (int i1=1; i1<n1; ++i1) {
        m2 = m1; m1 = p0; p0 = p1; ++p1; 
        if (p1>=n1) p1 = n1-1;
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float xa = xp0[p0]-xm1[m1];
        float xb = xm1[p0]-xp0[m1];
        float x1 = xa+xb+b*(xp1[p0]+xm2[p0]-xp1[m1]-xm2[m1]);
        float x2 = xa-xb+b*(xp0[p1]+xp0[m2]-xm1[p1]-xm1[m2]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        float yc = b*y1;
        float yd = b*y2;
        yp0[p0] += ya; ym1[m1] -= ya;
        ym1[p0] += yb; yp0[m1] -= yb;
        yp1[p0] += yc; ym2[m1] -= yc;
        ym2[p0] += yc; yp1[m1] -= yc;
        yp0[p1] += yd; ym1[m2] -= yd;
        yp0[m2] += yd; ym1[p1] -= yd;
      }
    }
  }
  public static void applyLhs242(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    float p = 0.18f;
    float a = 0.5f*(1.0f+p);
    float b = 0.5f*(    -p);
    c *= a*a;
    b /= a;
    if (x!=y) 
      scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    FilterBuffer2 fbx = new FilterBuffer2(2,1,2,1,x);
    FilterBuffer2 fby = new FilterBuffer2(2,1,2,1,y);
    fbx.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fby.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fbx.setMode(FilterBuffer2.Mode.INPUT);
    fby.setMode(FilterBuffer2.Mode.INPUT_OUTPUT);
    for (int i2=1; i2<n2; ++i2) {
      float[] xm2 = fbx.get(i2-2);
      float[] xm1 = fbx.get(i2-1);
      float[] xp0 = fbx.get(i2  );
      float[] xp1 = fbx.get(i2+1);
      float[] ym2 = fby.get(i2-2);
      float[] ym1 = fby.get(i2-1);
      float[] yp0 = fby.get(i2  );
      float[] yp1 = fby.get(i2+1);
      /*
      float xm2m2, xm2m1=xm2[1], xm2p0=xm2[2], xm2p1=xm2[3];
      float xm1m2, xm1m1=xm1[1], xm1p0=xm1[2], xm1p1=xm1[3];
      float xp0m2, xp0m1=xp0[1], xp0p0=xp0[2], xp0p1=xp0[3];
      float xp1m2, xp1m1=xp1[1], xp1p0=xp1[2], xp1p1=xp1[3];
      float ym2m2, ym2m1=ym2[1], ym2p0=ym2[2], ym2p1=ym2[3];
      float ym1m2, ym1m1=ym1[1], ym1p0=ym1[2], ym1p1=ym1[3];
      float yp0m2, yp0m1=yp0[1], yp0p0=yp0[2], yp0p1=yp0[3];
      float yp1m2, yp1m1=yp1[1], yp1p0=yp1[2], yp1p1=yp1[3];
      for (int i1=1,m2=1,p1=4; i1<n1; ++i1,++m2,++p1) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        xm2m2=xm2m1; xm2m1=xm2p0; xm2p0=xm2p1; xm2p1 = xm2[p1];
        xm1m2=xm1m1; xm1m1=xm1p0; xm1p0=xm1p1; xm1p1 = xm1[p1];
        xp0m2=xp0m1; xp0m1=xp0p0; xp0p0=xp0p1; xp0p1 = xp0[p1];
        xp1m2=xp1m1; xp1m1=xp1p0; xp1p0=xp1p1; xp1p1 = xp1[p1];
        ym2m2=ym2m1; ym2m1=ym2p0; ym2p0=ym2p1; ym2p1 = ym2[p1];
        ym1m2=ym1m1; ym1m1=ym1p0; ym1p0=ym1p1; ym1p1 = ym1[p1];
        yp0m2=yp0m1; yp0m1=yp0p0; yp0p0=yp0p1; yp0p1 = yp0[p1];
        yp1m2=yp1m1; yp1m1=yp1p0; yp1p0=yp1p1; yp1p1 = yp1[p1];
        float xa = xp0p0-xm1m1;
        float xb = xm1p0-xp0m1;
        float x1 = xa+xb+b*(xp1p0+xm2p0-xp1m1-xm2m1);
        float x2 = xa-xb+b*(xp0p1+xp0m2-xm1p1-xm1m2);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        float yc = b*y1;
        float yd = b*y2;
        yp0p0 += ya; ym1m1 -= ya;
        ym1p0 += yb; yp0m1 -= yb;
        yp1p0 += yc; ym2m1 -= yc;
        ym2p0 += yc; yp1m1 -= yc;
        yp0p1 += yd; ym1m2 -= yd;
        yp0m2 += yd; ym1p1 -= yd;
        ym2[m2] = ym2m2;
        ym1[m2] = ym1m2;
        yp0[m2] = yp0m2;
        yp1[m2] = yp1m2;
      }
      ym2[n1] = ym2m1; ym2[n1+1] = ym2p0; ym2[n1+2] = ym2p1;
      ym1[n1] = ym1m1; ym1[n1+1] = ym1p0; ym1[n1+2] = ym1p1;
      yp0[n1] = yp0m1; yp0[n1+1] = yp0p0; yp0[n1+2] = yp0p1;
      yp1[n1] = yp1m1; yp1[n1+1] = yp1p0; yp1[n1+2] = yp1p1;
      */
      for (int i1=1,m2=1,m1=2,p0=3,p1=4; i1<n1; ++i1,++m2,++m1,++p0,++p1) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float xa = xp0[p0]-xm1[m1];
        float xb = xm1[p0]-xp0[m1];
        float x1 = xa+xb+b*(xp1[p0]+xm2[p0]-xp1[m1]-xm2[m1]);
        float x2 = xa-xb+b*(xp0[p1]+xp0[m2]-xm1[p1]-xm1[m2]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        float yc = b*y1;
        float yd = b*y2;
        yp0[p0] += ya; ym1[m1] -= ya;
        ym1[p0] += yb; yp0[m1] -= yb;
        yp1[p0] += yc; ym2[m1] -= yc;
        ym2[p0] += yc; yp1[m1] -= yc;
        yp0[p1] += yd; ym1[m2] -= yd;
        yp0[m2] += yd; ym1[p1] -= yd;
      }
    }
    fby.flush();
  }

  public static void applyLhs33(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    float a =  10.0f/32.0f;
    float b =   3.0f/32.0f;
    c *= a*a;
    b /= a;
    scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    for (int i2=1; i2<n2-1; ++i2) {
      float[] xm = x[i2-1], x0 = x[i2], xp = x[i2+1];
      float[] ym = y[i2-1], y0 = y[i2], yp = y[i2+1];
      float xmm, xm0 = xm[0], xmp = xm[1];
      float x0m, x00 = x0[0], x0p = x0[1];
      float xpm, xp0 = xp[0], xpp = xp[1];
      float ymm, ym0 = ym[0], ymp = ym[1];
      float y0m, y00 = y0[0], y0p = y0[1];
      float ypm, yp0 = yp[0], ypp = yp[1];
      for (int i1m=0,i1=1,i1p=2; i1p<n1; ++i1m,++i1,++i1p) {
        xmm = xm0; xm0 = xmp; xmp = xm[i1p];
        x0m = x00; x00 = x0p; x0p = x0[i1p];
        xpm = xp0; xp0 = xpp; xpp = xp[i1p];
        ymm = ym0; ym0 = ymp; ymp = ym[i1p];
        y0m = y00; y00 = y0p; y0p = y0[i1p];
        ypm = yp0; yp0 = ypp; ypp = yp[i1p];
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float xa = b*(xpp-xmm);
        float xb = b*(xmp-xpm);
        float x1 = x0p-x0m+xa+xb;
        float x2 = xp0-xm0+xa-xb;
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = b*(y1+y2);
        float yb = b*(y1-y2);
        y0p += y1; y0m -= y1;
        ypp += ya; ymm -= ya;
        ymp += yb; ypm -= yb;
        yp0 += y2; ym0 -= y2;
        ym[i1m] = ymm;
        y0[i1m] = y0m;
        yp[i1m] = ypm;
      }
      ym[n1-2] = ym0; ym[n1-1] = ymp;
      y0[n1-2] = y00; y0[n1-1] = y0p;
      yp[n1-2] = yp0; yp[n1-1] = ypp;
    }
  }
  public static void applyLhs33X(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    float a = 10.0f/32.0f;
    float b =  3.0f/32.0f;
    c *= a*a;
    b /= a;
    scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    for (int i2=1; i2<n2-1; ++i2) {
      float[] xm1 = x[i2-1], xp0 = x[i2], xp1 = x[i2+1];
      float[] ym1 = y[i2-1], yp0 = y[i2], yp1 = y[i2+1];
      for (int i1=1,m1=0,p0=1,p1=2; p1<n1; ++i1,++m1,++p0,++p1) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float xa = b*(xp1[p1]-xm1[m1]);
        float xb = b*(xm1[p1]-xp1[m1]);
        float x1 = xp0[p1]-xp0[m1]+xa+xb;
        float x2 = xp1[p0]-xm1[p0]+xa-xb;
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = b*(y1+y2);
        float yb = b*(y1-y2);
        ym1[m1] -= ya; yp0[m1] -= y1; yp1[m1] -= yb;
        ym1[p0] -= y2;                yp1[p0] += y2;
        ym1[p1] += yb; yp0[p1] += y1; yp1[p1] += ya;
      }
    }
  }
  public static void applyLhs332(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    float a = 10.0f/32.0f;
    float b =  3.0f/32.0f;
    c *= a*a;
    b /= a;
    if (x!=y) 
      scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    FilterBuffer2 fbx = new FilterBuffer2(1,1,1,1,x);
    FilterBuffer2 fby = new FilterBuffer2(1,1,1,1,y);
    fbx.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fby.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fbx.setMode(FilterBuffer2.Mode.INPUT);
    fby.setMode(FilterBuffer2.Mode.INPUT_OUTPUT);
    for (int i2=0; i2<n2; ++i2) {
      float[] xm = fbx.get(i2-1), x0 = fbx.get(i2), xp = fbx.get(i2+1);
      float[] ym = fby.get(i2-1), y0 = fby.get(i2), yp = fby.get(i2+1);
      for (int i1=0,j1m=0,j1=1,j1p=2; i1<n1; ++i1,++j1m,++j1,++j1p) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float xmm = xm[j1m];
        float xm0 = xm[j1 ];
        float xmp = xm[j1p];
        float x0m = x0[j1m];
        float x0p = x0[j1p];
        float xpm = xp[j1m];
        float xp0 = xp[j1 ];
        float xpp = xp[j1p];
        float xa = xpp-xmm;
        float xb = xmp-xpm;
        float x1 = x0p-x0m+b*(xa+xb);
        float x2 = xp0-xm0+b*(xa-xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = b*(y1+y2);
        float yb = b*(y1-y2);
        y0[j1p] += y1;
        y0[j1m] -= y1;
        yp[j1p] += ya;
        ym[j1m] -= ya;
        ym[j1p] += yb;
        yp[j1m] -= yb;
        yp[j1 ] += y2;
        ym[j1 ] -= y2;
      }
    }
    fby.flush();
  }

  public static void applyLhs44(
    Tensors2 t, float cs, float[][] s, float[][] x, float[][] y) 
  {
    float p = -0.087584f;
    float a =  (0.5f-p)*27.0f/24.0f;
    float b =  (     p)*27.0f/24.0f;
    float c = -(0.5f-p)/24.0f;
    float d = -(     p)/24.0f;
    if (x!=y) 
      scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    FilterBuffer2 fbx = new FilterBuffer2(2,1,2,1,x);
    FilterBuffer2 fby = new FilterBuffer2(2,1,2,1,y);
    fbx.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fby.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fbx.setMode(FilterBuffer2.Mode.INPUT);
    fby.setMode(FilterBuffer2.Mode.INPUT_OUTPUT);
    for (int i2=1; i2<n2; ++i2) {
      float[] xm2 = fbx.get(i2-2);
      float[] xm1 = fbx.get(i2-1);
      float[] xp0 = fbx.get(i2  );
      float[] xp1 = fbx.get(i2+1);
      float[] ym2 = fby.get(i2-2);
      float[] ym1 = fby.get(i2-1);
      float[] yp0 = fby.get(i2  );
      float[] yp1 = fby.get(i2+1);
      for (int i1=1,j1=2+i1; i1<n1; ++i1,++j1) {
        int m2 = j1-2, m1 = j1-1, p0 = j1, p1 = j1+1;
        t.getTensor(i1,i2,di);
        float csi = (s!=null)?cs*s[i2][i1]:cs;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float xm2m2=xm2[m2], xm1m2=xm1[m2], xp0m2=xp0[m2], xp1m2=xp1[m2];
        float xm2m1=xm2[m1], xm1m1=xm1[m1], xp0m1=xp0[m1], xp1m1=xp1[m1];
        float xm2p0=xm2[p0], xm1p0=xm1[p0], xp0p0=xp0[p0], xp1p0=xp1[p0];
        float xm2p1=xm2[p1], xm1p1=xm1[p1], xp0p1=xp0[p1], xp1p1=xp1[p1];
        float x1 = a*(xp0p0+xm1p0-xp0m1-xm1m1) +
                   b*(xp1p0+xm2p0-xp1m1-xm2m1) +
                   c*(xp0m2+xm1m2-xp0p1-xm1p1) +
                   d*(xp1m2+xm2m2-xp1p1-xm2p1);
        float x2 = a*(xp0p0+xp0m1-xm1p0-xm1m1) +
                   b*(xp0p1+xp0m2-xm1p1-xm1m2) +
                   c*(xm2p0+xm2m1-xp1p0-xp1m1) +
                   d*(xm2p1+xm2m2-xp1p1-xp1m2);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ay1 = a*y1, ay2 = a*y2;
        float by1 = b*y1, by2 = b*y2;
        float cy1 = c*y1, cy2 = c*y2;
        float dy1 = d*y1, dy2 = d*y2;
        yp0[p0] += ay1+ay2;
        ym1[p0] += ay1-ay2;
        yp0[m1] -= ay1-ay2;
        ym1[m1] -= ay1+ay2;
        yp1[p0] += by1-cy2;
        ym2[p0] += by1+cy2;
        yp1[m1] -= by1+cy2;
        ym2[m1] -= by1-cy2;
        yp0[m2] += cy1+by2;
        ym1[m2] += cy1-by2;
        yp0[p1] -= cy1-by2;
        ym1[p1] -= cy1+by2;
        yp1[m2] += dy1-dy2;
        ym2[m2] += dy1+dy2;
        yp1[p1] -= dy1+dy2;
        ym2[p1] -= dy1-dy2;
      }
    }
    fby.flush();
  }

  private static float D21 =  0.6908013f;
  private static float D22 = -0.0899626f;
  public static void applyLhs51(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    if (x!=y) 
      scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    FilterBuffer2 fbx = new FilterBuffer2(2,2,2,2,x);
    FilterBuffer2 fby = new FilterBuffer2(2,2,2,2,y);
    fbx.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fby.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fbx.setMode(FilterBuffer2.Mode.INPUT);
    fby.setMode(FilterBuffer2.Mode.INPUT_OUTPUT);
    for (int i2=0; i2<n2; ++i2) {
      float[] x00 = fbx.get(i2);
      float[] xp1 = fbx.get(i2+1), xm1 = fbx.get(i2-1);
      float[] xp2 = fbx.get(i2+2), xm2 = fbx.get(i2-2);
      float[] y00 = fby.get(i2);
      float[] yp1 = fby.get(i2+1), ym1 = fby.get(i2-1);
      float[] yp2 = fby.get(i2+2), ym2 = fby.get(i2-2);
      for (int i1=0,j1=2; i1<n1; ++i1,++j1) {
        int j1m2 = j1-2, j1m1 = j1-1, j1p1 = j1+1, j1p2 = j1+2;
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x1 = D21*(x00[j1p1]-x00[j1m1]) +
                   D22*(x00[j1p2]-x00[j1m2]);
        float x2 = D21*(xp1[j1  ]-xm1[j1  ]) +
                   D22*(xp2[j1  ]-xm2[j1  ]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float d1y1 = D21*y1, d1y2 = D21*y2;
        float d2y1 = D22*y1, d2y2 = D22*y2;
        y00[j1p1] += d1y1; y00[j1m1] -= d1y1;
        y00[j1p2] += d2y1; y00[j1m2] -= d2y1;
        yp1[j1  ] += d1y2; ym1[j1  ] -= d1y2;
        yp2[j1  ] += d2y2; ym2[j1  ] -= d2y2;
      }
    }
    fby.flush();
  }

  //private static float D31 =  0.8630447f;
  //private static float D32 = -0.2677244f;
  //private static float D33 =  0.0670329f;
  private static float D31 =  0.830893f;
  private static float D32 = -0.227266f;
  private static float D33 =  0.042877f;
  //private static float D31 =  45.0f/60.0f;
  //private static float D32 =  -9.0f/60.0f;
  //private static float D33 =   1.0f/60.0f;

  public static void d7f(float[] x, float[] y) {
    int n1 = x.length;
    int n1m1 = n1-1, n1m2 = n1-2, n1m3 = n1-3, 
        n1m4 = n1-4, n1m5 = n1-5, n1m6 = n1-6;
    y[0] = D31*(x[1]-x[0]) +
           D32*(x[2]-x[0]) +
           D33*(x[3]-x[0]);
    y[1] = D31*(x[2]-x[0]) +
           D32*(x[3]-x[0]) +
           D33*(x[4]-x[0]);
    y[2] = D31*(x[3]-x[1]) +
           D32*(x[4]-x[0]) +
           D33*(x[5]-x[0]);
    float xm3, xm2 = x[0], xm1 = x[1], xp0 = x[2], 
               xp1 = x[3], xp2 = x[4], xp3 = x[5];
    for (int i1=3; i1<n1m3; ++i1) {
      xm3 = xm2; xm2 = xm1; xm1 = xp0;
      xp0 = xp1; xp1 = xp2; xp2 = xp3;
      xp3 = x[i1+3];
      y[i1] = D31*(xp1-xm1) +
              D32*(xp2-xm2) +
              D33*(xp3-xm3);
    }
    y[n1m3] = D31*(x[n1m2]-x[n1m4]) +
              D32*(x[n1m1]-x[n1m5]) +
              D33*(x[n1m1]-x[n1m6]);
    y[n1m2] = D31*(x[n1m1]-x[n1m3]) +
              D32*(x[n1m1]-x[n1m4]) +
              D33*(x[n1m1]-x[n1m5]);
    y[n1m1] = D31*(x[n1m1]-x[n1m2]) +
              D32*(x[n1m1]-x[n1m3]) +
              D33*(x[n1m1]-x[n1m4]);
  }
  public static void d7t(float[] x, float[] y) {
    int n1 = x.length;
    int n1m1 = n1-1, n1m2 = n1-2, n1m3 = n1-3, 
        n1m4 = n1-4, n1m5 = n1-5, n1m6 = n1-6;
    y[0] += D31*(-x[0]-x[1]) +
            D32*(-x[0]-x[1]-x[2]) +
            D33*(-x[0]-x[1]-x[2]-x[3]);
    y[1] += D31*( x[0]-x[2]) +
            D32*(     -x[3]) +
            D33*(     -x[4]);
    y[2] += D31*( x[1]-x[3]) +
            D32*( x[0]-x[4]) +
            D33*(     -x[5]);
    float xm3, xm2 = x[0], xm1 = x[1], xp0 = x[2], 
               xp1 = x[3], xp2 = x[4], xp3 = x[5];
    for (int i1=3; i1<n1m3; ++i1) {
      xm3 = xm2; xm2 = xm1; xm1 = xp0;
      xp0 = xp1; xp1 = xp2; xp2 = xp3;
      xp3 = x[i1+3];
      y[i1] += D31*(xm1-xp1) +
               D32*(xm2-xp2) +
               D33*(xm3-xp3);
    }
    y[n1m3] += D31*(x[n1m4]-x[n1m2]) +
               D32*(x[n1m5]-x[n1m1]) +
               D33*(x[n1m6]        );
    y[n1m2] += D31*(x[n1m3]-x[n1m1]) +
               D32*(x[n1m4]        ) +
               D33*(x[n1m5]        );
    y[n1m1] += D31*(x[n1m1]+x[n1m2]) +
               D32*(x[n1m1]+x[n1m2]+x[n1m3]) +
               D33*(x[n1m1]+x[n1m2]+x[n1m3]+x[n1m4]);
  }
  public static void g7f1(float[][] x, float[][] g1) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      d7f(x[i2],g1[i2]);
  }
  public static void g7f2(float[][] x, float[][] g2) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n2m1 = n2-1, n2m2 = n2-2, n2m3 = n2-3, 
        n2m4 = n2-4, n2m5 = n2-5, n2m6 = n2-6;
    for (int i2=0; i2<n2; ++i2) {
      float[] xm3 = (i2>2)?x[i2-3]:x[0];
      float[] xm2 = (i2>1)?x[i2-2]:x[0];
      float[] xm1 = (i2>0)?x[i2-1]:x[0];
      float[] xp1 = (i2<n2m1)?x[i2+1]:x[n2m1];
      float[] xp2 = (i2<n2m2)?x[i2+2]:x[n2m1];
      float[] xp3 = (i2<n2m3)?x[i2+3]:x[n2m1];
      float[] g2i = g2[i2];
      for (int i1=0; i1<n1; ++i1) {
        g2i[i1] = D31*(xp1[i1]-xm1[i1]) +
                  D32*(xp2[i1]-xm2[i1]) +
                  D33*(xp3[i1]-xm3[i1]);
      }
    }
  }
  public static void g7t1(float[][] g1, float[][] x) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      d7t(g1[i2],x[i2]);
  }
  public static void g7t2(float[][] g2, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n2m1 = n2-1, n2m2 = n2-2, n2m3 = n2-3, 
        n2m4 = n2-4, n2m5 = n2-5, n2m6 = n2-6;
    for (int i1=0; i1<n1; ++i1) {
      x[0][i1] += D31*(-g2[0][i1]-g2[1][i1]) +
                  D32*(-g2[0][i1]-g2[1][i1]-g2[2][i1]) +
                  D33*(-g2[0][i1]-g2[1][i1]-g2[2][i1]-g2[3][i1]);
      x[1][i1] += D31*( g2[0][i1]-g2[2][i1]) +
                  D32*(          -g2[3][i1]) +
                  D33*(          -g2[4][i1]);
      x[2][i1] += D31*( g2[1][i1]-g2[3][i1]) +
                  D32*( g2[0][i1]-g2[4][i1]) +
                  D33*(          -g2[5][i1]);
    }
    for (int i2=3; i2<n2m3; ++i2) {
      float[] gm3 = g2[i2-3];
      float[] gm2 = g2[i2-2];
      float[] gm1 = g2[i2-1];
      float[] gp1 = g2[i2+1];
      float[] gp2 = g2[i2+2];
      float[] gp3 = g2[i2+3];
      float[] x2i = x[i2];
      for (int i1=0; i1<n1; ++i1) {
        x2i[i1] += D31*(gm1[i1]-gp1[i1]) +
                   D32*(gm2[i1]-gp2[i1]) +
                   D33*(gm3[i1]-gp3[i1]);
      }
    }
    for (int i1=0; i1<n1; ++i1) {
      x[n2m3][i1] += D31*(g2[n2m4][i1]-g2[n2m2][i1]) +
                     D32*(g2[n2m5][i1]-g2[n2m1][i1]) +
                     D33*(g2[n2m6][i1]             );
      x[n2m2][i1] += D31*(g2[n2m3][i1]-g2[n2m1][i1]) +
                     D32*(g2[n2m4][i1]             ) +
                     D33*(g2[n2m5][i1]             );
      x[n2m1][i1] += D31*(g2[n2m1][i1]+g2[n2m2][i1]) +
                     D32*(g2[n2m1][i1]+g2[n2m2][i1]+g2[n2m3][i1]) +
                     D33*(g2[n2m1][i1]+g2[n2m2][i1]+g2[n2m3][i1]+g2[n2m4][i1]);
    }
  }
  public static void g7f(float[][] x, float[][] g1, float[][] g2) {
    g7f1(x,g1);
    g7f2(x,g2);
  }
  public static void g7t(float[][] g1, float[][] g2, float[][] x) {
    g7t1(g1,x);
    g7t2(g2,x);
  }
  public static void testGrad1() {
    int n = 21;
    float[] x = randfloat(n);
    float[] y = randfloat(n);
    //float[] x = zerofloat(n); x[0] = x[n/2] = x[n-1] = 1.0f;
    //float[] y = zerofloat(n); y[0] = y[n/2] = y[n-1] = 1.0f;
    float[] gx = zerofloat(n);
    float[] gy = zerofloat(n);
    d7f(x,gx); // Gx
    d7t(y,gy); // G'y
    dump(gx);
    dump(gy);
    float ygx = sum(mul(y,gx)); // y'Gx
    float xgy = sum(mul(x,gy)); // x'G'y
    trace("ygx="+ygx);
    trace("xgy="+xgy);
  }
  public static void testGrad2() {
    int n1 = 21;
    int n2 = 22;
    float[][] x = randfloat(n1,n2);
    float[][] y1 = randfloat(n1,n2);
    float[][] y2 = randfloat(n1,n2);
    float[][] y = zerofloat(n1,n2);
    float[][] x1 = zerofloat(n1,n2);
    float[][] x2 = zerofloat(n1,n2);
    g7f(x,x1,x2);
    g7t(y1,y2,y);
    float ygx = sum(add(mul(y1,x1),mul(y2,x2)));
    float xgy = sum(mul(x,y));
    trace("ygx="+ygx);
    trace("xgy="+xgy);
  }
  public static void main(String[] args) {
    testGrad1();
    testGrad2();
  }
  public static void applyLhs712(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    scopy(x,y);
    float[] di = new float[3];
    float[][] x1 = new float[n2][n1];
    float[][] x2 = new float[n2][n1];
    float[][] y1 = x1;
    float[][] y2 = x2;
    g7f(x,x1,x2);
    for (int i2=0; i2<n2; ++i2) {
      float[] x1i2 = x1[i2];
      float[] x2i2 = x2[i2];
      float[] y1i2 = y1[i2];
      float[] y2i2 = y2[i2];
      for (int i1=0; i1<n1; ++i1) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x1i = x1i2[i1];
        float x2i = x2i2[i1];
        y1i2[i1] = d11*x1i+d12*x2i;
        y2i2[i1] = d12*x1i+d22*x2i;
      }
    }
    g7t(y1,y2,y);
  }
  public static void applyLhs71(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    scopy(x,y);
    float[] di = new float[3];
    int i2m3,i2m2=0,i2m1=0,i2p0=0,i2p1=0,i2p2=1,i2p3=2;
    for (int i2=0; i2<n2; ++i2) {
      i2m3 = i2m2; i2m2 = i2m1; i2m1 = i2p0;
      i2p0 = i2p1; i2p1 = i2p2; i2p2 = i2p3; ++i2p3;
      if (i2p1>=n2) i2p1 = n2-1;
      if (i2p2>=n2) i2p2 = n2-1;
      if (i2p3>=n2) i2p3 = n2-1;
      float[] xm3 = x[i2m3], xm2 = x[i2m2], xm1 = x[i2m1];
      float[] xp3 = x[i2p3], xp2 = x[i2p2], xp1 = x[i2p1];
      float[] xp0 = x[i2p0];
      float[] ym3 = y[i2m3], ym2 = y[i2m2], ym1 = y[i2m1];
      float[] yp3 = y[i2p3], yp2 = y[i2p2], yp1 = y[i2p1];
      float[] yp0 = y[i2p0];
      int m3,m2=0,m1=0,p0=0,p1=0,p2=1,p3=2;
      for (int i1=0; i1<n1; ++i1) {
        m3 = m2; m2 = m1; m1 = p0;
        p0 = p1; p1 = p2; p2 = p3; ++p3;
        if (p1>=n1) p1 = n1-1;
        if (p2>=n1) p2 = n1-1;
        if (p3>=n1) p3 = n1-1;
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x1 = D31*(xp0[p1]-xp0[m1]) +
                   D32*(xp0[p2]-xp0[m2]) +
                   D33*(xp0[p3]-xp0[m3]);
        float x2 = D31*(xp1[p0]-xm1[p0]) +
                   D32*(xp2[p0]-xm2[p0]) +
                   D33*(xp3[p0]-xm3[p0]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float d1y1 = D31*y1, d1y2 = D31*y2;
        float d2y1 = D32*y1, d2y2 = D32*y2;
        float d3y1 = D33*y1, d3y2 = D33*y2;
        yp0[p1] += d1y1; yp0[m1] -= d1y1;
        yp1[p0] += d1y2; ym1[p0] -= d1y2;
        yp0[p2] += d2y1; yp0[m2] -= d2y1;
        yp2[p0] += d2y2; ym2[p0] -= d2y2;
        yp0[p3] += d3y1; yp0[m3] -= d3y1;
        yp3[p0] += d3y2; ym3[p0] -= d3y2;
      }
    }
  }

  private static float D41 =  0.8947167f;
  private static float D42 = -0.3153471f;
  private static float D43 =  0.1096895f;
  private static float D44 = -0.0259358f;
  public static void applyLhs91(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    scopy(x,y);
    float[] di = new float[3];
    int i2m4,i2m3=0,i2m2=0,i2m1=0,i2p0=0,i2p1=0,i2p2=1,i2p3=2,i2p4=3;
    for (int i2=0; i2<n2; ++i2) {
      i2m4 = i2m3; i2m3 = i2m2; i2m2 = i2m1; i2m1 = i2p0;
      i2p0 = i2p1; i2p1 = i2p2; i2p2 = i2p3; i2p3 = i2p4; ++i2p4;
      if (i2p1>=n2) i2p1 = n2-1;
      if (i2p2>=n2) i2p2 = n2-1;
      if (i2p3>=n2) i2p3 = n2-1;
      if (i2p4>=n2) i2p4 = n2-1;
      float[] xm4 = x[i2m4], xm3 = x[i2m3], xm2 = x[i2m2], xm1 = x[i2m1];
      float[] xp4 = x[i2p4], xp3 = x[i2p3], xp2 = x[i2p2], xp1 = x[i2p1];
      float[] xp0 = x[i2p0];
      float[] ym4 = y[i2m4], ym3 = y[i2m3], ym2 = y[i2m2], ym1 = y[i2m1];
      float[] yp4 = y[i2p4], yp3 = y[i2p3], yp2 = y[i2p2], yp1 = y[i2p1];
      float[] yp0 = y[i2p0];
      int m4,m3=0,m2=0,m1=0,p0=0,p1=0,p2=1,p3=2,p4=3;
      for (int i1=0; i1<n1; ++i1) {
        m4 = m3; m3 = m2; m2 = m1; m1 = p0;
        p0 = p1; p1 = p2; p2 = p3; p3 = p4; ++p4;
        if (p1>=n1) p1 = n1-1;
        if (p2>=n1) p2 = n1-1;
        if (p3>=n1) p3 = n1-1;
        if (p4>=n1) p4 = n1-1;
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x1 = D41*(xp0[p1]-xp0[m1]) +
                   D42*(xp0[p2]-xp0[m2]) +
                   D43*(xp0[p3]-xp0[m3]) +
                   D44*(xp0[p4]-xp0[m4]);
        float x2 = D41*(xp1[p0]-xm1[p0]) +
                   D42*(xp2[p0]-xm2[p0]) +
                   D43*(xp3[p0]-xm3[p0]) +
                   D44*(xp4[p0]-xm4[p0]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float d1y1 = D41*y1, d1y2 = D41*y2;
        float d2y1 = D42*y1, d2y2 = D42*y2;
        float d3y1 = D43*y1, d3y2 = D43*y2;
        float d4y1 = D44*y1, d4y2 = D44*y2;
        yp0[p1] += d1y1; yp0[m1] -= d1y1;
        yp0[p2] += d2y1; yp0[m2] -= d2y1;
        yp0[p3] += d3y1; yp0[m3] -= d3y1;
        yp0[p4] += d4y1; yp0[m4] -= d4y1;
        yp1[p0] += d1y2; ym1[p0] -= d1y2;
        yp2[p0] += d2y2; ym2[p0] -= d2y2;
        yp3[p0] += d3y2; ym3[p0] -= d3y2;
        yp4[p0] += d4y2; ym4[p0] -= d4y2;
      }
    }
  }

  public static void applyLhs912(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    if (x!=y) 
      scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    FilterBuffer2 fbx = new FilterBuffer2(4,4,4,4,x);
    FilterBuffer2 fby = new FilterBuffer2(4,4,4,4,y);
    fbx.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fby.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fbx.setMode(FilterBuffer2.Mode.INPUT);
    fby.setMode(FilterBuffer2.Mode.INPUT_OUTPUT);
    for (int i2=0; i2<n2; ++i2) {
      float[] x00 = fbx.get(i2);
      float[] xp1 = fbx.get(i2+1), xm1 = fbx.get(i2-1);
      float[] xp2 = fbx.get(i2+2), xm2 = fbx.get(i2-2);
      float[] xp3 = fbx.get(i2+3), xm3 = fbx.get(i2-3);
      float[] xp4 = fbx.get(i2+4), xm4 = fbx.get(i2-4);
      float[] y00 = fby.get(i2);
      float[] yp1 = fby.get(i2+1), ym1 = fby.get(i2-1);
      float[] yp2 = fby.get(i2+2), ym2 = fby.get(i2-2);
      float[] yp3 = fby.get(i2+3), ym3 = fby.get(i2-3);
      float[] yp4 = fby.get(i2+4), ym4 = fby.get(i2-4);
      for (int i1=0,j1=4; i1<n1; ++i1,++j1) {
        int j1m2 = j1-2, j1m1 = j1-1, j1p1 = j1+1, j1p2 = j1+2;
        int j1m4 = j1-4, j1m3 = j1-3, j1p3 = j1+3, j1p4 = j1+4;
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x1 = D41*(x00[j1p1]-x00[j1m1]) +
                   D42*(x00[j1p2]-x00[j1m2]) +
                   D43*(x00[j1p3]-x00[j1m3]) +
                   D44*(x00[j1p4]-x00[j1m4]);
        float x2 = D41*(xp1[j1  ]-xm1[j1  ]) +
                   D42*(xp2[j1  ]-xm2[j1  ]) +
                   D43*(xp3[j1  ]-xm3[j1  ]) +
                   D44*(xp4[j1  ]-xm4[j1  ]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float d1y1 = D41*y1, d1y2 = D41*y2;
        float d2y1 = D42*y1, d2y2 = D42*y2;
        float d3y1 = D43*y1, d3y2 = D43*y2;
        float d4y1 = D44*y1, d4y2 = D44*y2;
        y00[j1p1] += d1y1; y00[j1m1] -= d1y1;
        y00[j1p2] += d2y1; y00[j1m2] -= d2y1;
        y00[j1p3] += d3y1; y00[j1m3] -= d3y1;
        y00[j1p4] += d4y1; y00[j1m4] -= d4y1;
        yp1[j1  ] += d1y2; ym1[j1  ] -= d1y2;
        yp2[j1  ] += d2y2; ym2[j1  ] -= d2y2;
        yp3[j1  ] += d3y2; ym3[j1  ] -= d3y2;
        yp4[j1  ] += d4y2; ym4[j1  ] -= d4y2;
      }
    }
    fby.flush();
  }

  private static float D61 =  0.9483655f; 
  private static float D62 = -0.4030691f;
  private static float D63 =  0.2023745f;
  private static float D64 = -0.0988080f;
  private static float D65 =  0.0421765f;
  private static float D66 = -0.0133106f;
  public static void applyLhs131(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    if (x!=y) 
      scopy(x,y);
    int n1 = x[0].length;
    int n2 = x.length;
    float[] di = new float[3];
    FilterBuffer2 fbx = new FilterBuffer2(6,6,6,6,x);
    FilterBuffer2 fby = new FilterBuffer2(6,6,6,6,y);
    fbx.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fby.setExtrapolation(FilterBuffer2.Extrapolation.ZERO_SLOPE);
    fbx.setMode(FilterBuffer2.Mode.INPUT);
    fby.setMode(FilterBuffer2.Mode.INPUT_OUTPUT);
    for (int i2=0; i2<n2; ++i2) {
      float[] x00 = fbx.get(i2);
      float[] xp1 = fbx.get(i2+1), xm1 = fbx.get(i2-1);
      float[] xp2 = fbx.get(i2+2), xm2 = fbx.get(i2-2);
      float[] xp3 = fbx.get(i2+3), xm3 = fbx.get(i2-3);
      float[] xp4 = fbx.get(i2+4), xm4 = fbx.get(i2-4);
      float[] xp5 = fbx.get(i2+5), xm5 = fbx.get(i2-5);
      float[] xp6 = fbx.get(i2+6), xm6 = fbx.get(i2-6);
      float[] y00 = fby.get(i2);
      float[] yp1 = fby.get(i2+1), ym1 = fby.get(i2-1);
      float[] yp2 = fby.get(i2+2), ym2 = fby.get(i2-2);
      float[] yp3 = fby.get(i2+3), ym3 = fby.get(i2-3);
      float[] yp4 = fby.get(i2+4), ym4 = fby.get(i2-4);
      float[] yp5 = fby.get(i2+5), ym5 = fby.get(i2-5);
      float[] yp6 = fby.get(i2+6), ym6 = fby.get(i2-6);
      for (int i1=0,j1=6; i1<n1; ++i1,++j1) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x1 = D61*(x00[j1+1]-x00[j1-1]) +
                   D62*(x00[j1+2]-x00[j1-2]) +
                   D63*(x00[j1+3]-x00[j1-3]) +
                   D64*(x00[j1+4]-x00[j1-4]) +
                   D65*(x00[j1+5]-x00[j1-5]) +
                   D66*(x00[j1+6]-x00[j1-6]);
        float x2 = D61*(xp1[j1  ]-xm1[j1  ]) +
                   D62*(xp2[j1  ]-xm2[j1  ]) +
                   D63*(xp3[j1  ]-xm3[j1  ]) +
                   D64*(xp4[j1  ]-xm4[j1  ]) +
                   D65*(xp5[j1  ]-xm5[j1  ]) +
                   D66*(xp6[j1  ]-xm6[j1  ]);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float y11 = D61*y1, y21 = D61*y2;
        float y12 = D62*y1, y22 = D62*y2;
        float y13 = D63*y1, y23 = D63*y2;
        float y14 = D64*y1, y24 = D64*y2;
        float y15 = D65*y1, y25 = D65*y2;
        float y16 = D66*y1, y26 = D66*y2;
        y00[j1+1] += y11; y00[j1-1] -= y11;
        y00[j1+2] += y12; y00[j1-2] -= y12;
        y00[j1+3] += y13; y00[j1-3] -= y13;
        y00[j1+4] += y14; y00[j1-4] -= y14;
        y00[j1+5] += y15; y00[j1-5] -= y15;
        y00[j1+6] += y16; y00[j1-6] -= y16;
        yp1[j1  ] += y21; ym1[j1  ] -= y21;
        yp2[j1  ] += y22; ym2[j1  ] -= y22;
        yp3[j1  ] += y23; ym3[j1  ] -= y23;
        yp4[j1  ] += y24; ym4[j1  ] -= y24;
        yp5[j1  ] += y25; ym5[j1  ] -= y25;
        yp6[j1  ] += y26; ym6[j1  ] -= y26;
      }
    }
    fby.flush();
  }

  /**
   * Computes y = (S'S+G'DG)x. Arrays x and y must be distinct.
   */
  private static void applyLhsX(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    applyRhs(x,y);
    float[] di = new float[3];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        d.getTensor(i1,i2,di);
        float csi = (s!=null)?c*s[i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.25f*(xa-xb);
        float x2 = 0.25f*(xa+xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void applyLhsXX(
    Tensors2 d, float c, float[][] s, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    scopy(x,y);
    float[] di = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
        d.getTensor(i1p,i2p,di);
        float csi = 0.5f*((s!=null)?c*s[i2p][i1p]:c);
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d22 = di[2]*csi;
        float t = min(d11,d22,abs(d12));
        //float t = 0.1666667f*(d11+d22);
        //float t = 0.3333333f*(d11+d22);
        //float t = 0.4f*(d11+d22);
        //float t = 0.5f*(d11+d22);
        float xpp = x[i2p][i1p];
        float xpm = x[i2p][i1m];
        float xmp = x[i2m][i1p];
        float xmm = x[i2m][i1m];
        float apppm = (d11-t)*(xpp-xpm);
        float ampmm = (d11-t)*(xmp-xmm);
        float bppmm = (d12+t)*(xpp-xmm);
        float bpmmp = (d12-t)*(xpm-xmp);
        float cppmp = (d22-t)*(xpp-xmp);
        float cpmmm = (d22-t)*(xpm-xmm);
        y[i2p][i1p] += apppm+bppmm+cppmp;
        y[i2p][i1m] -= apppm+bpmmp-cpmmm;
        y[i2m][i1p] += ampmm+bpmmp-cppmp;
        y[i2m][i1m] -= ampmm+bppmm+cpmmm;
      }
    }
  }

  private static void applyLhs(
    Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    if (SMOOTH) {
      szero(y);
    } else {
      scopy(x,y);
    }
    if (PARALLEL) {
      applyLhsParallel(d,c,s,x,y);
    } else {
      applyLhsSerial(d,c,s,x,y);
    }
  }

  private static void applyLhsSerial(
    Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3)
      applyLhsSlice3(i3,d,c,s,x,y);
  }

  private static void applyLhsParallel(
    final Tensors3 d, final float c, final float[][][] s, 
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;

    // i3 = 1, 3, 5, ...
    final AtomicInteger a1 = new AtomicInteger(1);
    Thread[] thread1 = Threads.makeArray();
    for (int ithread=0; ithread<thread1.length; ++ithread) {
      thread1[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a1.getAndAdd(2); i3<n3; i3=a1.getAndAdd(2))
            applyLhsSlice3(i3,d,c,s,x,y);
        }
      });
    }
    Threads.startAndJoin(thread1);

    // i3 = 2, 4, 6, ...
    final AtomicInteger a2 = new AtomicInteger(2);
    Thread[] thread2 = Threads.makeArray();
    for (int ithread=0; ithread<thread2.length; ++ithread) {
      thread2[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a2.getAndAdd(2); i3<n3; i3=a2.getAndAdd(2))
            applyLhsSlice3(i3,d,c,s,x,y);
        }
      });
    }
    Threads.startAndJoin(thread2);
  }


  /**
   * Computes y = (S'S+G'DG)x for one constant-i3 slice.
   */
  private static void applyLhsSlice3(
    int i3, Tensors3 d, float c, float[][][] s, float[][][] x, float[][][] y) 
  {
    float[] di = new float[6];
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        d.getTensor(i1,i2,i3,di);
        float csi = (s!=null)?c*s[i3][i2][i1]:c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d13 = di[2]*csi;
        float d22 = di[3]*csi;
        float d23 = di[4]*csi;
        float d33 = di[5]*csi;
        applyLhs(i1,d11,d12,d13,d22,d23,d33,x00,x01,x10,x11,y00,y01,y10,y11);
      }
    }
  }

  /**
   * Computes y = (S'S+G'DG)x for one sample.
   */
  private static void applyLhs(int i1,
   float d11, float d12, float d13, float d22, float d23, float d33,
   float[] x00, float[] x01, float[] x10, float[] x11,
   float[] y00, float[] y01, float[] y10, float[] y11)
  {
    int i1m = i1-1;
    float x000 = x00[i1 ];
    float x001 = x00[i1m];
    float x010 = x01[i1 ];
    float x011 = x01[i1m];
    float x100 = x10[i1 ];
    float x101 = x10[i1m];
    float x110 = x11[i1 ];
    float x111 = x11[i1m];
    //float x1 = 0.0625f*(x000+x010+x100+x110-x001-x011-x101-x111);
    //float x2 = 0.0625f*(x000+x001+x100+x101-x010-x011-x110-x111);
    //float x3 = 0.0625f*(x000+x001+x010+x011-x100-x101-x110-x111);
    float xa = x000-x111;
    float xb = x001-x110;
    float xc = x010-x101;
    float xd = x100-x011;
    float x1 = 0.0625f*(xa-xb+xc+xd);
    float x2 = 0.0625f*(xa+xb-xc+xd);
    float x3 = 0.0625f*(xa+xb+xc-xd);
    float y1 = d11*x1+d12*x2+d13*x3;
    float y2 = d12*x1+d22*x2+d23*x3;
    float y3 = d13*x1+d23*x2+d33*x3;
    float ya = y1+y2+y3;
    float yb = y1-y2+y3;
    float yc = y1+y2-y3;
    float yd = y1-y2-y3;
    if (SMOOTH) {
      float xs = 0.015625f*(x000+x001+x010+x011+x100+x101+x110+x111);
      y00[i1 ] += ya+xs;
      y00[i1m] -= yd-xs;
      y01[i1 ] += yb+xs;
      y01[i1m] -= yc-xs;
      y10[i1 ] += yc+xs;
      y10[i1m] -= yb-xs;
      y11[i1 ] += yd+xs;
      y11[i1m] -= ya-xs;
    } else {
      y00[i1 ] += ya;
      y00[i1m] -= yd;
      y01[i1 ] += yb;
      y01[i1m] -= yc;
      y10[i1 ] += yc;
      y10[i1m] -= yb;
      y11[i1 ] += yd;
      y11[i1m] -= ya;
    }
  }

  /**
   * Solves Ax = b via conjugate gradient iterations. (No preconditioner.)
   * Uses the initial values of x; does not assume they are zero.
   */
  private void solve(Operator2 a, float[][] b, float[][] x) {
    int n1 = b[0].length;
    int n2 = b.length;
    float[][] d = new float[n2][n1];
    float[][] q = new float[n2][n1];
    float[][] r = new float[n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    scopy(r,d);
    float delta = sdot(r,r);
    float deltaBegin = delta;
    float deltaSmall = sdot(b,b)*_small*_small;
    log.fine("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      log.finer("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
      a.apply(d,q);
      float dq = sdot(d,q);
      float alpha = delta/dq;
      saxpy( alpha,d,x);
      saxpy(-alpha,q,r);
      float deltaOld = delta;
      delta = sdot(r,r);
      float beta = delta/deltaOld;
      sxpay(beta,r,d);
    }
    log.fine("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }
  private void solve(Operator3 a, float[][][] b, float[][][] x) {
    int n1 = b[0][0].length;
    int n2 = b[0].length;
    int n3 = b.length;
    float[][][] d = new float[n3][n2][n1];
    float[][][] q = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    scopy(b,r);
    a.apply(x,q);
    saxpy(-1.0f,q,r); // r = b-Ax
    scopy(r,d);
    float delta = sdot(r,r);
    float deltaBegin = delta;
    float deltaSmall = sdot(b,b)*_small*_small;
    log.fine("solve: delta="+delta);
    int iter;
    for (iter=0; iter<_niter && delta>deltaSmall; ++iter) {
      log.finer("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
      a.apply(d,q);
      float dq = sdot(d,q);
      float alpha = delta/dq;
      saxpy( alpha,d,x);
      saxpy(-alpha,q,r);
      float deltaOld = delta;
      delta = sdot(r,r);
      float beta = delta/deltaOld;
      sxpay(beta,r,d);
    }
    log.fine("  iter="+iter+" delta="+delta+" ratio="+delta/deltaBegin);
  }

  // Zeros array x.
  private static void szero(float[][] x) {
    zero(x);
  }
  private static void szero(float[][][] x) {
    if (PARALLEL) {
      szeroP(x);
    } else {
      szeroS(x);
    }
  }
  private static void szeroS(float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      szero(x[i3]);
  }
  private static void szeroP(final float[][][] x) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            szero(x[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Copys array x to array y.
  private static void scopy(float[][] x, float[][] y) {
    copy(x,y);
  }
  private static void scopy(float[][][] x, float[][][] y) {
    if (PARALLEL) {
      scopyP(x,y);
    } else {
      scopyS(x,y);
    }
  }
  private static void scopyS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      scopy(x[i3],y[i3]);
  }
  private static void scopyP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            scopy(x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Returns the dot product x'y.
  private static float sdot(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float d = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        d += x2[i1]*y2[i1];
      }
    }
    return d;
  }
  private static float sdot(float[][][] x, float[][][] y) {
    if (PARALLEL) {
      return sdotP(x,y);
    } else {
      return sdotS(x,y);
    }
  }
  private static float sdotS(float[][][] x, float[][][] y) {
    int n3 = x.length;
    float d = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      d += sdot(x[i3],y[i3]);
    return d;
  }
  private static float sdotP(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final AtomicFloat ad = new AtomicFloat(0.0f);
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          float d = 0.0f;
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            d += sdot(x[i3],y[i3]);
          ad.getAndAdd(d);
        }
      });
    }
    Threads.startAndJoin(threads);
    return ad.get();
  }

  // Computes y = y + a*x.
  private static void saxpy(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        y2[i1] += a*x2[i1];
      }
    }
  }
  private static void saxpy(float a, float[][][] x, float[][][] y) {
    if (PARALLEL) {
      saxpyP(a,x,y);
    } else {
      saxpyS(a,x,y);
    }
  }
  private static void saxpyS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      saxpy(a,x[i3],y[i3]);
  }
  private static void saxpyP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            saxpy(a,x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Computes y = x + a*y.
  private static void sxpay(float a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2], y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        y2[i1] = a*y2[i1]+x2[i1];
      }
    }
  }
  private static void sxpay(float a, float[][][] x, float[][][] y) {
    if (PARALLEL) {
      sxpayP(a,x,y);
    } else {
      sxpayS(a,x,y);
    }
  }
  private static void sxpayS(float a, float[][][] x, float[][][] y) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      sxpay(a,x[i3],y[i3]);
  }
  private static void sxpayP(
    final float a, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    final AtomicInteger a3 = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=a3.getAndIncrement(); i3<n3; i3=a3.getAndIncrement())
            sxpay(a,x[i3],y[i3]);
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  public static void xmain(String[] args) {
    LocalSmoothingFilterX lsf = new LocalSmoothingFilterX(0.00001,100);
    Tensors2 t = new Tensors2() {
      public void getTensor(int i1, int i2, float[] a) {
        a[0] = 1.0f; a[1] = 0.0f; a[2] = 0.0f;
      }
    };
    int n1=5,n2=5;
    //float[][] x = randfloat(n1,n2);
    //float[][] y = randfloat(n1,n2);
    //float[][] x = fillfloat(1.0f,n1,n2);
    //float[][] y = fillfloat(1.0f,n1,n2);
    float[][] x = zerofloat(n1,n2);
    float[][] y = zerofloat(n1,n2);
    x[2][2] = 1.0f;
    y[2][2] = 1.0f;
    float[][] sx = like(x);
    float[][] sy = like(y);
    lsf.apply(t,x,sx); // Sx
    lsf.applyTranspose(t,y,sy); // S'y
    //applyRhs(x,sx);
    //applyLhs(t,1.0f,null,y,sy);
    dump(sx);
    dump(sy);
    float ysx = sdot(y,sx); // y'(Sx)
    float xsy = sdot(x,sy); // x'(S'y)
    trace("ysx="+ysx+" xsy="+xsy);
  }
}
