/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import edu.mines.jtk.dsp.FftComplex;
import edu.mines.jtk.dsp.FftReal;
import edu.mines.jtk.sgl.SimpleFrame; // FOR DEVELOPMENT ONLY
import edu.mines.jtk.util.Check;

import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

/**
 * 3D smoothing along fault planes with a rotated Gaussian filter.
 * This smoothing filter is implemented with FFTs, and is optimized 
 * for the case of repeated applications for different rotations to 
 * multiple arrays.
 * <p>
 * A fault plane is defined by two angles phi and theta in the range 
 * [-90,90] degrees. The angle phi is the strike of the fault, and 
 * the angle theta is its deviation from vertical. 
 * <p>
 * To define these angles more precisely, let vectors u, v and w be
 * initially aligned with axes 3, 2 and 1 of a 3D array, respectively.
 * (Axes 3, 2, and 1 typically correspond to axes X, Y, and Z, where
 * axis 1 and axis Z are vertical coordinates increasing downward.)
 * First rotate the vectors (u,v,w) by angle phi about vector w (which
 * is aligned with axis 1), where angle phi is positive clockwise when
 * viewed from the head of vector u towards its tail. The vectors u and
 * v may now point in new directions, but both still lie in the 23
 * plane. Then rotate the vectors u and w by angle theta about the
 * vector v, where angle theta is positive counterclockwise when viewed
 * from the head of vector v towards its tail. After these rotations,
 * vector v remains in the 23 plane. The vectors v and w lie within the
 * fault plane (w points down the fault, and v is horizontal), and the
 * vector u is normal to that plane.
 * <p>
 * Let ct = cos(theta), st = sin(theta), cp = cos(phi), sp = sin(phi),
 * where both angles theta and phi are in the range [-90,90] degrees. 
 * The three orthonormal vectors (u,v,w) after the two rotations are 
 * u = {-st,-ct*sp, ct*cp}.
 * v = {  0,    cp,    sp},
 * w = { ct,-st*sp, st*cp}, 
 * <p>
 * Gaussian smoothing is specified by half-widths sigma for each of 
 * the principle directions of the rotated vectors (u,v,w).
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.06.24
 */
public class FaultPlaneSmoother {

  /**
   * Constructs a smoother for specified parameters and input.
   * @param sigmau half-width in u direction down the fault.
   * @param sigmav half-width in v direction along fault strike.
   * @param sigmaw half-width in w direction normal to fault.
   * @param f input 3D array to be smoothed.
   */
  public FaultPlaneSmoother(
    double sigmau, double sigmav, double sigmaw, float[][][] f)
  {
    this(sigmau,sigmav,sigmaw,new float[][][][]{f});
  }

  /**
   * Constructs a smoother for specified parameters and input.
   * @param sigmau half-width in u direction down the fault.
   * @param sigmav half-width in v direction along fault strike.
   * @param sigmaw half-width in w direction normal to fault.
   * @param f input array of 3D arrays to be smoothed.
   */
  public FaultPlaneSmoother(
    double sigmau, double sigmav, double sigmaw, float[][][][] f)
  {
    _sigmau = (float)sigmau;
    _sigmav = (float)sigmav;
    _sigmaw = (float)sigmaw;
    _n1 = f[0][0][0].length;
    _n2 = f[0][0].length;
    _n3 = f[0].length;
    _nf = f.length;

     // Construct FFTs for all three array dimensions.
    int npad = (int)(3.0f*max(_sigmau,_sigmav,_sigmaw));
    int n1pad = _n1+npad;
    int n2pad = _n2+npad;
    int n3pad = _n3+npad;
    _n1fft = FftComplex.nfftSmall(n1pad);
    _n2fft = FftComplex.nfftSmall(n2pad);
    _n3fft = FftReal.nfftSmall(n3pad);
    _nk1 = _n1fft;
    _nk2 = _n2fft;
    _nk3 = _n3fft/2+1;
    _fft1 = new FftComplex(_n1fft);
    _fft2 = new FftComplex(_n2fft);
    _fft3 = new FftReal(_n3fft);

    // Space for inputs after FFT for only the 3rd dimension.
    // The 3rd dimension is slowest, and we save space by not 
    // padding and transforming over all three dimensions.
    _fk = new float[_nf][_nk3][_n2][_n1*2];

    // Perform forward FFT over 3rd dimension only.
    for (int i=0; i<_nf; ++i)
      fft3Forward(f[i],_fk[i]);
  }
  /**
   * Applies this filter for specified fault angles.
   * @param phi fault strike angle in degrees [-90,90].
   * @param theta fault angle from vertical in degrees [-90,90].
   * @param g output smoothed 3D array.
   */
  public void apply(double phi, double theta, float[][][] g) {
    apply(phi,theta,new float[][][][]{g});
  }

  /**
   * Applies this filter for specified fault angles.
   * @param phi fault strike angle in degrees [-90,90].
   * @param theta fault angle from vertical in degrees [-90,90].
   * @param g output array of smoothed 3D arrays.
   */
  public void apply(double phi, double theta, float[][][][] g) {
    Check.argument(_n1==g[0][0][0].length,"dimension n1 of g is valid");
    Check.argument(_n2==g[0][0].length,"dimension n2 of g is valid");
    Check.argument(_n3==g[0].length,"dimension n3 of g is valid");
    Check.argument(_nf==g.length,"number of arrays in g is valid");
    //float[][][] h = makeFilter(phi,theta);
    float[][][] gk = new float[_nk3][_n2][_n1*2];
    for (int i=0; i<_nf; ++i) { // for all arrays, ...
      //applyFilter(h,_fk[i],gk); // multiply by transform of filter
      applyFilter(phi,theta,_fk[i],gk); // multiply by transform of filter
      fft3Inverse(gk,g[i]); // inverse FFT over 3rd dimension only
    }
  }

  /**
   * Applies this filter for specified fault angles.
   * @param phi fault strike angle in degrees [-90,90].
   * @param theta fault angle from vertical in degrees [-90,90].
   * @return output array of smoothed 3D arrays.
   */
  public float[][][][] apply(double phi, double theta) {
    float[][][][] g = new float[_nf][_n3][_n2][_n1];
    apply(phi,theta,g);
    return g;
  }

  private float _sigmau,_sigmav,_sigmaw;
  private int _n1,_n2,_n3,_nf;
  private int _n1fft,_n2fft,_n3fft;
  private int _nk1,_nk2,_nk3;
  private FftComplex _fft1;
  private FftComplex _fft2;
  private FftReal _fft3;
  private float[][][][] _fk; // inputs after FFT over axis 3

  private float[][][] makeFilter(double phi, double theta) {
    float p = (float)toRadians(phi);
    float t = (float)toRadians(theta);
    float cp = cos(p), sp = sin(p);
    float ct = cos(t), st = sin(t);
    final float u1 =  -st, u2 = -ct*sp, u3 = ct*cp; // u normal to fault
    final float v1 = 0.0f, v2 =     cp, v3 =    sp; // v along strike
    final float w1 =   ct, w2 = -st*sp, w3 = st*cp; // w down the fault
    final float twopi = 2.0f*FLT_PI;
    final float dk1 = twopi/_n1fft;
    final float dk2 = twopi/_n2fft;
    final float dk3 = twopi/_n3fft;
    final float sigmaus = _sigmau*_sigmau;
    final float sigmavs = _sigmav*_sigmav;
    final float sigmaws = _sigmaw*_sigmaw;
    final float hscale = 1.0f/_n1fft/_n2fft/_n3fft;
    final float[][][] h = new float[_nk3][_nk2][_nk1];
    loop(_nk3,new LoopInt() {
    public void compute(int i3) {
      float k3 = i3*dk3;
      float uk3 = u3*k3;
      float vk3 = v3*k3;
      float wk3 = w3*k3;
      for (int i2=0; i2<_nk2; ++i2) {
        float k2 = i2*dk2;
        if (i2*2>_nk2) k2 -= twopi;
        float uk23 = u2*k2+uk3;
        float vk23 = v2*k2+vk3;
        float wk23 = w2*k2+wk3;
        float[] h32 = h[i3][i2];
        for (int i1=0; i1<_nk1; ++i1) {
          float k1 = i1*dk1;
          if (i1*2>_nk1) k1 -= twopi;
          float uk = u1*k1+uk23;
          float vk = v1*k1+vk23;
          float wk = w1*k1+wk23;
          float s = sigmaus*uk*uk+sigmavs*vk*vk+sigmaws*wk*wk;
          if (s<10.0f)
            h32[i1] = exp(-0.5f*s)*hscale;
        }
      }
    }});
    return h;
  }

  private void applyFilter(float[][][] h, float[][][] f, float[][][] g) {
    // arrays h[nk3][nk2][nk1], f[nk3][n2][2*n1], g[nk3][n2][2*n1]
    final float[][][] hh = h;
    final float[][][] ff = f;
    final float[][][] gg = g;
    loop(_nk3,new LoopInt() {
    public void compute(int i3) {
      float[][] gk = new float[_nk2][_nk1*2];
      copy(2*_n1,_n2,ff[i3],gk);
      _fft2.complexToComplex2(-1,_n1,gk,gk);
      _fft1.complexToComplex1(-1,_nk2,gk,gk);
      for (int i2=0; i2<_nk2; ++i2) {
        float[] g32 = gk[i2];
        float[] h32 = hh[i3][i2];
        for (int i1=0,i1r=0,i1i=1; i1<_nk1; ++i1,i1r+=2,i1i+=2) {
          float hi = h32[i1];
          g32[i1r] *= hi;
          g32[i1i] *= hi;
        }
      }
      _fft1.complexToComplex1(1,_nk2,gk,gk);
      _fft2.complexToComplex2(1,_n1,gk,gk);
      copy(2*_n1,_n2,gk,gg[i3]);
    }});
  }

  private void applyFilter(
      double phi, double theta, float[][][] f,  float[][][] g) {
    final float[][][] ff = f;
    final float[][][] gg = g;
    float p = (float)toRadians(phi);
    float t = (float)toRadians(theta);
    float cp = cos(p), sp = sin(p);
    float ct = cos(t), st = sin(t);
    final float u1 =  -st, u2 = -ct*sp, u3 = ct*cp; // u normal to fault
    final float v1 = 0.0f, v2 =     cp, v3 =    sp; // v along strike
    final float w1 =   ct, w2 = -st*sp, w3 = st*cp; // w down the fault
    final float twopi = 2.0f*FLT_PI;
    final float dk1 = twopi/_n1fft;
    final float dk2 = twopi/_n2fft;
    final float dk3 = twopi/_n3fft;
    final float sigmaus = _sigmau*_sigmau;
    final float sigmavs = _sigmav*_sigmav;
    final float sigmaws = _sigmaw*_sigmaw;
    final float hscale = 1.0f/_n1fft/_n2fft/_n3fft;
    loop(_nk3,new LoopInt() {
    public void compute(int i3) {
      float[][] gk = new float[_nk2][_nk1*2];
      copy(2*_n1,_n2,ff[i3],gk);
      _fft2.complexToComplex2(-1,_n1,gk,gk);
      _fft1.complexToComplex1(-1,_nk2,gk,gk);
      float k3 = i3*dk3;
      float uk3 = u3*k3;
      float vk3 = v3*k3;
      float wk3 = w3*k3;
      for (int i2=0; i2<_nk2; ++i2) {
        float[] g32 = gk[i2];
        float k2 = i2*dk2;
        if (i2*2>_nk2) k2 -= twopi;
        float uk23 = u2*k2+uk3;
        float vk23 = v2*k2+vk3;
        float wk23 = w2*k2+wk3;
        for (int i1=0,i1r=0,i1i=1; i1<_nk1; ++i1,i1r+=2,i1i+=2) {
          float k1 = i1*dk1;
          if (i1*2>_nk1) k1 -= twopi;
          float uk = u1*k1+uk23;
          float vk = v1*k1+vk23;
          float wk = w1*k1+wk23;
          float s = sigmaus*uk*uk+sigmavs*vk*vk+sigmaws*wk*wk;
          float hi = (s<10.0f) ? exp(-0.5f*s)*hscale : 0.0f;
          g32[i1r] *= hi;
          g32[i1i] *= hi;
        }
      }
      _fft1.complexToComplex1(1,_nk2,gk,gk);
      _fft2.complexToComplex2(1,_n1,gk,gk);
      ccopy(_n1,_n2,gk,gg[i3]);
    }});
  }

  // Forward FFT over 3rd dimension only.
  private void fft3Forward(final float[][][] fx, final float[][][] fk) {
    loop(_n2,new LoopInt() {
    public void compute(int i2) {
      float[][] fxpad = new float[_n3fft][_n1];
      float[][] fkpad = new float[_nk3][_n1*2];
      for (int i3=0; i3<_n3; ++i3)
        copy(fx[i3][i2],fxpad[i3]);
      _fft3.realToComplex2(-1,_n1,fxpad,fkpad);
      for (int i3=0; i3<_nk3; ++i3)
        ccopy(fkpad[i3],fk[i3][i2]);
    }});
  }

  // Inverse FFT over 3rd dimension only.
  private void fft3Inverse(final float[][][] fk, final float[][][] fx) {
    loop(_n2,new LoopInt() {
    public void compute(int i2) {
      float[][] fkpad = new float[_nk3][_n1*2];
      float[][] fxpad = new float[_n3fft][_n1];
      for (int i3=0; i3<_nk3; ++i3)
        ccopy(fk[i3][i2],fkpad[i3]);
      _fft3.complexToReal2(1,_n1,fkpad,fxpad);
      for (int i3=0; i3<_n3; ++i3)
        copy(fxpad[i3],fx[i3][i2]);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    testSmoother();
  }
  public static void testSmoother() {
    int n1 = 201;
    int n2 = 202;
    int n3 = 203;
    float[][][] f1 = zerofloat(n1,n2,n3);
    float[][][] f2 = zerofloat(n1,n2,n3);
    float[][][] g1 = zerofloat(n1,n2,n3);
    float[][][] g2 = zerofloat(n1,n2,n3);
    f1[n3/2][n2/2][n1/2] = 1.0f;
    f2[n3/4][n2/4][n1/4] = 1.0f;
    float[][][][] f = {f1,f2};
    float[][][][] g = {g1,g2};
    double sigmau = 1.0;
    double sigmav = 15.0;
    double sigmaw = 30.0;
    double phi = 30.0;
    double theta = 45.0;
    FaultPlaneSmoother fps = new FaultPlaneSmoother(sigmau,sigmav,sigmaw,f);
    fps.apply(phi,theta,g);
    SimpleFrame.asImagePanels(g1);
    SimpleFrame.asImagePanels(g2);
    System.out.println("g1 max="+max(g1));
    System.out.println("g2 max="+max(g2));
  }
}
