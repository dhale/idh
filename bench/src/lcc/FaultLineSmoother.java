/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

// FOR DEVELOPMENT ONLY
import edu.mines.jtk.mosaic.*;

/**
 * 2D smoothing along fault lines with a rotated Gaussian filter.
 * This smoothing filter is implemented with FFTs, and is optimized 
 * for the case of repeated applications for different rotations to 
 * multiple arrays.
 * <p>
 * A fault line is defined by an angle theta in the range [-90,90] 
 * degrees. The angle theta is the line's deviation from vertical. 
 * <p>
 * To define this angle more precisely, let vectors u and v be 
 * initially aligned with axes 1 and 2 of a 2D array, respectively. 
 * (Axis 1 is typically the vertical axis.) Rotate the vectors (u,v)
 * counter-clockwise by angle theta. After this rotation, the vector
 * u points down the fault line, and the vector v is normal to that 
 * line.
 * <p>
 * Let ct = cos(theta) and st = sin(theta). After rotation, the two 
 * orthonormal vectors (u,v) are u = { ct, st}, and v = {-st, ct}.
 * <p>
 * Gaussian smoothing is specified by half-widths sigma for each of 
 * the principle directions of the rotated vectors (u,v).
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.06.24
 */
public class FaultLineSmoother {

  /**
   * Constructs a smoother for specified parameters and input.
   * @param sigmau half-width in u direction, down the fault.
   * @param sigmav half-width in v direction, normal to the fault.
   * @param f input 2D array to be smoothed.
   */
  public FaultLineSmoother(
    double sigmau, double sigmav, float[][] f) 
  {
    this(sigmau,sigmav,new float[][][]{f});
  }

  /**
   * Constructs a smoother for specified parameters and input.
   * @param sigmau half-width in u direction, down the fault.
   * @param sigmav half-width in v direction, normal to the fault.
   * @param f array of 2D input arrays to be smoothed.
   */
  public FaultLineSmoother(
    double sigmau, double sigmav, float[][][] f) 
  {
    _sigmau = (float)sigmau;
    _sigmav = (float)sigmav;
    _n1 = f[0][0].length;
    _n2 = f[0].length;
    _nf = f.length;
    int npad = (int)(3.0f*max(_sigmau,_sigmav));
    int n1pad = _n1+npad;
    int n2pad = _n2+npad;
    _n1fft = FftReal.nfftSmall(n1pad);
    _n2fft = FftComplex.nfftSmall(n2pad);
    _nk1 = _n1fft/2+1;
    _nk2 = _n2fft;
    _fft1 = new FftReal(_n1fft);
    _fft2 = new FftComplex(_n2fft);
    _f2 = new float[_nf][_n2fft][_n1fft+2];
    for (int jf=0; jf<_nf; ++jf) {
      zero(_f2[jf]);
      copy(_n1,_n2,f[jf],_f2[jf]);
      _fft1.realToComplex1(-1,_n2,_f2[jf],_f2[jf]);
      _fft2.complexToComplex2(-1,_nk1,_f2[jf],_f2[jf]);
    }
  }

  /**
   * Applies this filter for specified fault angle.
   * @param theta fault angle from vertical in degrees [-90,90].
   * @param g output smoothed 2D array.
   */
  public void apply(double theta, float[][] g) {
    apply(theta,new float[][][]{g});
  }

  /**
   * Applies this filter for specified fault angle.
   * @param theta fault angle from vertical in degrees [-90,90].
   * @param g output array of smoothed 2D arrays.
   */
  public void apply(double theta, float[][][] g) {
    Check.argument(_n1==g[0][0].length,"dimension n1 of g is valid");
    Check.argument(_n2==g[0].length,"dimension n2 of g is valid");
    Check.argument(_nf==g.length,"number of arrays in g is valid");
    float t = (float)toRadians(theta);
    float ct = cos(t), st = sin(t);
    float u1 =  ct, u2 = st;
    float v1 = -st, v2 = ct;
    float twopi = 2.0f*FLT_PI;
    float dk1 = twopi/_n1fft;
    float dk2 = twopi/_n2fft;
    float sigmaus = _sigmau*_sigmau;
    float sigmavs = _sigmav*_sigmav;
    float hscale = 1.0f/_n1fft/_n2fft;
    float[][] h2 = new float[_nk2][_nk1];
    for (int jk2=0; jk2<_nk2; ++jk2) {
      float k2 = jk2*dk2;
      if (jk2*2>_nk2) k2 -= twopi;
      float[] h2j = h2[jk2];
      for (int jk1=0; jk1<_nk1; ++jk1) {
        float k1 = jk1*dk1;
        float uk = u1*k1+u2*k2;
        float vk = v1*k1+v2*k2;
        float s = sigmaus*uk*uk+sigmavs*vk*vk;
        if (s<10.0f)
          h2j[jk1] = exp(-0.5f*s)*hscale;
      }
    }
    float[][] g2 = new float[_n2fft][_n1fft+2];
    for (int jg=0; jg<_nf; ++jg) {
      for (int jk2=0; jk2<_nk2; ++jk2) {
        float[] f2j = _f2[jg][jk2];
        float[] g2j = g2[jk2];
        float[] h2j = h2[jk2];
        for (int jk1=0,jk1r=0,jk1i=1; jk1<_nk1; ++jk1,jk1r+=2,jk1i+=2) {
          float h2jj = h2j[jk1];
          g2j[jk1r] = h2jj*f2j[jk1r];
          g2j[jk1i] = h2jj*f2j[jk1i];
        }
      }
      _fft2.complexToComplex2(1,_nk1,g2,g2);
      _fft1.complexToReal1(1,_n2,g2,g2);
      copy(_n1,_n2,g2,g[jg]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigmau,_sigmav;
  private int _n1,_n2,_nf;
  private int _n1fft,_n2fft;
  private int _nk1,_nk2;
  private FftReal _fft1;
  private FftComplex _fft2;
  private float[][][] _f2;

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    testSmoother();
  }
  public static void testSmoother() {
    int n1 = 357;
    int n2 = 357;
    float[][] f = zerofloat(n1,n2);
    float[][] g = zerofloat(n1,n2);
    f[   0][   0] = 1.0f;
    f[   0][n1-1] = 1.0f;
    f[n2-1][   0] = 1.0f;
    f[n2-1][n1-1] = 1.0f;
    f[n2/2][n1/2] = 1.0f;
    double sigmau = 30.0;
    double sigmav = 1.0;
    double theta = 15.0;
    FaultLineSmoother fls = new FaultLineSmoother(sigmau,sigmav,f);
    fls.apply(theta,g);
    SimplePlot.asPixels(f).addColorBar();
    SimplePlot.asPixels(g).addColorBar();
  }
}
