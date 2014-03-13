/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package util;

import java.util.Random;

import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Generates fake data for use in testing and demonstrations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2010.12.13
 */
public class FakeData {

  /**
   * Plots the specified fake data. Data type is specified by
   * method name. For example, type "seismic3d2010A" corresponds to
   * the method seismic3d2010A().
   */
  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable(){
      public void run() {
        go(args);
      }
    });
  }
  private static void go(String[] args) {
    if (args.length==0 || args[0].equals("seismic3d2010A")) {
      float[][][] f = seismic3d2010A();
      SimpleFrame frame = new SimpleFrame();
      frame.addImagePanels(f);
      frame.getOrbitView().setScale(2.0);
      frame.setSize(900,900);
    } else if (args[0].equals("seismic2d2011A")) {
      float[][] f = seismic2d2011A(251,501,45.0);
        SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
        PixelsView pv = sp.addPixels(f);
        pv.setInterpolation(PixelsView.Interpolation.LINEAR);
        sp.setSize(1200,600);
    } else if (args[0].equals("seismicAndShifts2d2011A")) {
      float[][][] fgsr = seismicAndShifts2d2011A(251,501,45.0);
      for (int i=0; i<6; ++i) {
        SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
        PixelsView pv = sp.addPixels(fgsr[i]);
        pv.setInterpolation(PixelsView.Interpolation.LINEAR);
        if (i>=2) {
          pv.setClips(-35.0f,35.0f);
          pv.setColorModel(ColorMap.JET);
        } else {
          pv.setClips(-1.0f,1.0f);
        }
        sp.addColorBar();
        sp.setSize(1220,600);
      }
    } else if (args[0].equals("seismicAndSlopes2d2014A")) {
      float[][][] fs = seismicAndSlopes2d2014A();
      for (float[][] f:fs) {
        SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
        PixelsView pv = sp.addPixels(f);
        pv.setInterpolation(PixelsView.Interpolation.NEAREST);
        pv.setClips(-3.0f,3.0f);
        sp.addColorBar();
        sp.setSize(710,640);
        sp.getPlotPanel().setColorBarWidthMinimum(50);
      }
    } else {
      System.out.println("unrecognized type of fake data");
    }
  }

  /**
   * Returns a fake 2D seismic image, version 2011A.
   * This version simulates different types of structural deformation.
   * For the specified dip, the first n2/2 traces are rotated,
   * and the last n2/2 traces are vertically sheared. Shifts are
   * tapered and smoothed to reduce shifts near the middle trace
   * with index n2/2 and near edges of the images.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param dip maximum dip angle, in degrees.
   * @return the fake image.
   */
  public static float[][] seismic2d2011A(
    int n1, int n2, double dip)
  {
    float[][][] fgsr = seismicAndShifts2d2011A(n1,n2,dip);
    return fgsr[1];
  }

  /**
   * Returns fake 2D seismic images and shifts, version 2011A.
   * This version simulates different types of structural deformation.
   * For the specified dip, the first n2/2 traces are rotated,
   * and the last n2/2 traces are vertically sheared. Shifts are
   * tapered and smoothed to reduce shifts near the middle trace
   * with index n2/2 and near edges of the images.
   * <p>
   * Elements of the returned array {f,g,s1,s2,r1,r2} are 2D arrays
   * for a transform g(x1,x2) = f(x1+s1(x1,x2),x2+s2(x1,x2)) 
   * and its inverse f(u1,u2) = g(u1-r1(u1,u2),u2-r2(u1,u2)).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param dip maximum dip angle, in degrees.
   * @return array {f,g,s1,s2,r1,r2} of images and shifts.
   */
  public static float[][][] seismicAndShifts2d2011A(
    int n1, int n2, double dip)
  {
    float fpeak = 0.125f; // 1/4 of Nyquist
    float fmax = 2.0f*fpeak;
    float[][] f = makeEvents(n1,n2);
    f = addRickerWavelet(fpeak,f);
    f = mul(1.0f/max(abs(f)),f);
    float[][][] gsr = addRotateAndShear(f,dip,n2/8.0);
    //f = applyShiftsR(new float[][][]{gsr[3],gsr[4]},gsr[0]); // for testing 
    return new float[][][]{f,gsr[0],gsr[1],gsr[2],gsr[3],gsr[4]};
  }

  /**
   * Returns a fake 3D seismic image, version 2010A. The image includes 
   * default structure, a fault, an unconformity, amplitude variations,
   * and random noise.
   * @return the fake 3D seismic image.
   */
  public static float[][][] seismic3d2010A() {
    return seismic3d2010A(201,201,201,20.0,10.0,30.0,0.5,0.5);
  }

  /**
   * Returns a fake 2D seismic image, version 2012A. As options, the image 
   * may include structure, a fault, an unconformity, and random noise.
   * <p>
   * The image is initially a random sequence of horizontal reflections.
   * Random structure may be specified with a maximum vertical shift, 
   * but all shifts are limited to avoid aliasing. If a fault is 
   * specified, displacement along the fault increases with vertical 
   * depth (or time). If an erosional unconformity is specified, events 
   * above the specified sample are replaced with horizontal events.
   * Smooth random variations in signal amplitudes may be specified
   * by a minimum scale factor less than one. After these steps, the 
   * image is convolved vertically with a Ricker wavelet that has
   * a peak frequency equal to 1/5 of the Nyquist frequency. Finally, 
   * bandlimited random noise may be added with a specified rms
   * amplitude, relative to the rms signal amplitude.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param structure maximum vertical shift when adding structure.
   * @param fault maximum displacement when adding a fault.
   * @param erosion sample at which to create an erosional unconformity.
   * @param amplitude minimum scale factor when scaling amplitudes.
   * @param noise rms of added noise, relative to rms signal.
   * @return the fake 2D seismic image.
   */
  public static float[][] seismic2d2012A(
    int n1, int n2,
    double structure, double fault, double erosion, 
    double amplitude, double noise) 
  {
    float fpeak = 0.2f;
    float fmax = 2.0f*fpeak;
    float smaxStructure = (float)structure;
    float smaxFault = (float)fault;
    int kerosion = (int)erosion;
    float aminScale = (float)amplitude;
    float rmsNoise = (float)noise;
    float[][] f = makeEvents(n1,n2);
    if (smaxStructure>0.0f) 
      f = addStructure(smaxStructure,fmax,f);
    if (smaxFault>0.0f)
      f = addFault(smaxFault,f);
    if (kerosion>0)
      f = addErosion(kerosion,f);
    if (aminScale<1.0f)
      f = addAmplitude(aminScale,f);
    f = addRickerWavelet(fpeak,f);
    if (rmsNoise>0.0f)
      f = addNoise(rmsNoise,f);
    f = mul(1.0f/max(abs(f)),f);
    return f;
  }

  /**
   * Returns a fake 3D seismic image, version 2010A. As options, the image 
   * may include structure, a fault, an unconformity, and random noise.
   * <p>
   * The image is initially a random sequence of horizontal reflections.
   * Random structure may be specified with a maximum vertical shift, 
   * but all shifts are limited to avoid aliasing. If a fault is 
   * specified, displacement along the fault increases with vertical 
   * depth (or time). If an erosional unconformity is specified, events 
   * above the specified sample are replaced with horizontal events.
   * Smooth random variations in signal amplitudes may be specified
   * by a minimum scale factor less than one. After these steps, the 
   * image is convolved vertically with a Ricker wavelet that has
   * a peak frequency equal to 1/5 of the Nyquist frequency. Finally, 
   * bandlimited random noise may be added with a specified rms
   * amplitude, relative to the rms signal amplitude.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param structure maximum vertical shift when adding structure.
   * @param fault maximum displacement when adding a fault.
   * @param erosion sample at which to create an erosional unconformity.
   * @param amplitude minimum scale factor when scaling amplitudes.
   * @param noise rms of added noise, relative to rms signal.
   * @return the fake 3D seismic image.
   */
  public static float[][][] seismic3d2010A(
    int n1, int n2, int n3,
    double structure, double fault, double erosion, 
    double amplitude, double noise) 
  {
    float fpeak = 0.2f;
    float fmax = 2.0f*fpeak;
    float smaxStructure = (float)structure;
    float smaxFault = (float)fault;
    int kerosion = (int)erosion;
    float aminScale = (float)amplitude;
    float rmsNoise = (float)noise;
    float[][][] f = makeEvents(n1,n2,n3);
    if (smaxStructure>0.0f) 
      f = addStructure(smaxStructure,fmax,f);
    if (smaxFault>0.0f)
      f = addFault(smaxFault,f);
    if (kerosion>0)
      f = addErosion(kerosion,f);
    if (aminScale<1.0f)
      f = addAmplitude(aminScale,f);
    f = addRickerWavelet(fpeak,f);
    if (rmsNoise>0.0f)
      f = addNoise(rmsNoise,f);
    f = mul(1.0f/max(abs(f)),f);
    return f;
  }

  /**
   * Returns a fake noisy 2D seismic image with noise-free slopes.
   * The image and slopes are designed to be used to test methods for
   * estimating slopes of locally linear features apparent in 2D seismic
   * images. The image contains sinusoidal folding, horizontal and dipping
   * layers, two unconformities, and two intersecting faults with throws that
   * increase linearly with depth. The rms noise-to-signal ratio is 0.5.
   */
  public static float[][][] seismicAndSlopes2d2014A() {
    return seismicAndSlopes2d2014A(0.5f);
  }

  /**
   * Returns a fake 2D seismic image with slopes.
   * The fake seismic image contains sinusoidal folding, horizontal and
   * dipping layers, two unconformities, and two intersecting faults with
   * throws that increase linearly with depth. While the image may have
   * a specified amount of additive noise, the slopes are noise-free.
   * @param noise rms of noise (relative to signal) added to the image.
   */
  public static float[][][] seismicAndSlopes2d2014A(double noise) {
    int n1 = 501;
    int n2 = 501;
    float[][][] p = makeReflectivityWithNormals(n1,n2);
    float[][][] q = makeReflectivityWithNormals(n1,n2);
    float[][][] r = makeReflectivityWithNormals(n1,n2);
    Linear1 throw1 = new Linear1(0.0f,0.10f);
    Linear1 throw2 = new Linear1(0.0f,0.10f);
    LinearFault2 fault1 = new LinearFault2(n2*0.2f,15.0f,throw1);
    LinearFault2 fault2 = new LinearFault2(n2*0.4f,-15.0f,throw2);
    Sinusoidal2 fold = new Sinusoidal2(0.0f,0.05f,1.0e-4f,2.0e-4f);
    VerticalShear2 shear = new VerticalShear2(new Linear1(0.0f,0.05f));
    p = apply(fold,p);
    p = combine(n1/3,q,p);
    p = apply(shear,p);
    p = combine(n1/6,r,p);
    p = apply(fault1,p);
    p = apply(fault2,p);
    p = addWavelet(0.1,p);
    p[0] = addNoise((float)noise,p[0]);
    p[1] = neg(div(p[2],p[1]));
    return new float[][][]{p[0],p[1]};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static SincInterp _si = new SincInterp();
  static {
    _si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
  }

  /**
   * Adds both rotation and shear to the specified image.
   * Returns an array {g,s1,s2,r1,r2} such that the transform is
   * g(x1,x2) = f(x1+s1(x1,x2),x2+s2(x1,x2)) and its inverse is
   * f(x1,x2) = g(u1-r1(u1,u2),u2-r2(u1,u2)).
   */
  private static float[][][] addRotateAndShear(
    float[][] f, double dip, double sigma) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][][] r = makeRotateAndShear(n1,n2,dip,sigma);
    float[][][] s = makeShiftsSFromR(r);
    float[][] g = applyShiftsS(s,f);
    return new float[][][]{g,s[0],s[1],r[0],r[1]};
  }

  /**
   * Returns shifts r for both rotation and shear.
   */
  private static float[][][] makeRotateAndShear(
    int n1, int n2, double dip, double sigma)
  {
    double c1r = 2*n1/4;
    double c2r = 1*n2/4;
    double c1s = 2*n1/4;
    double c2s = 3*n2/4;
    int i2lo = n2/2-n2/16;
    int i2hi = n2/2+n2/16;
    float[][][] rr = makeRotateOrShear(n1,n2,c1r,c2r,dip,sigma,true);
    float[][][] rs = makeRotateOrShear(n1,n2,c1s,c2s,dip,sigma,false);
    float[][] r1r = rr[0], r2r = rr[1];
    float[][] r1s = rs[0], r2s = rs[1];
    float[][] r1 = new float[n2][n1];
    float[][] r2 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float wr,ws;
      if (i2<=i2lo) {
        wr = 1.0f;
        ws = 0.0f;
      } else if (i2>=i2hi) {
        wr = 0.0f;
        ws = 1.0f;
      } else {
        wr = 0.5f+0.5f*cos(FLT_PI*(float)(i2-i2lo)/(float)(i2hi-i2lo));
        ws = 1.0f-wr;
      }
      for (int i1=0; i1<n1; ++i1) {
        r1[i2][i1] = wr*r1r[i2][i1]+ws*r1s[i2][i1];
        r2[i2][i1] = wr*r2r[i2][i1]+ws*r2s[i2][i1];
      }
    }
    return new float[][][]{r1,r2};
  }

  /**
   * Returns shifts r for either rotation or shear, centered at (c1,c2).
   */
  private static float[][][] makeRotateOrShear(
    int n1, int n2, double c1, double c2, 
    double dip, double sigma, boolean rotate)
  {
    float[][] r1 = new float[n2][n1];
    float[][] r2 = new float[n2][n1];
    double cdip = cos(toRadians(dip));
    double sdip = sin(toRadians(dip));
    for (int i2=0; i2<n2; ++i2) {
      double u2 = i2-c2;
      for (int i1=0; i1<n1; ++i1) {
        double u1 = i1-c1;
        double x1,x2;
        if (rotate) {
          x1 = u1*cdip+u2*sdip;
          x2 = u2*cdip-u1*sdip;
        } else {
          x1 = u1+u2*sdip/cdip;
          x2 = u2;
        }
        double e = exp(-0.5*(u1*u1+u2*u2)/(sigma*sigma));
        r1[i2][i1] = (float)(e*(u1-x1));
        r2[i2][i1] = (float)(e*(u2-x2));
      }
    }
    return new float[][][]{r1,r2};
  }

  /**
   * Applies shifts s(x) via g(x) = f(x+s(x)).
   * s[0] and s[1] contain the shifts s1(x1,x2) and s2(x1,x2).
   */
  private static float[][] applyShiftsS(float[][][] s, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] s1 = s[0];
    float[][] s2 = s[1];
    SincInterp si = new SincInterp();
    si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
    float[][] g = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      double x2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        double x1 = i1;
        g[i2][i1] = si.interpolate(
          n1,1.0,0.0,n2,1.0,0.0,f,x1+s1[i2][i1],x2+s2[i2][i1]);
      }
    }
    return g;
  }

  /**
   * Applies shifts r(u) via f(u) = g(u-r(u)).
   * r[0] and r[1] contain the shifts r1(u1,u2) and r2(u1,u2).
   */
  private static float[][] applyShiftsR(float[][][] r, float[][] g) {
    int n1 = g[0].length;
    int n2 = g.length;
    float[][] r1 = r[0];
    float[][] r2 = r[1];
    SincInterp si = new SincInterp();
    si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      double u2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        double u1 = i1;
        f[i2][i1] = si.interpolate(
          n1,1.0,0.0,n2,1.0,0.0,g,u1-r1[i2][i1],u2-r2[i2][i1]);
      }
    }
    return f;
  }

  /**
   * Returns shifts s such that s(x) = r(x-s(x)).
   */
  private static float[][][] makeShiftsSFromR(float[][][] r) {
    float[][] r1 = r[0];
    float[][] r2 = r[1];
    int n1 = r1[0].length;
    int n2 = r1.length;
    float[][] s1 = copy(r1);
    float[][] s2 = copy(r2);
    LinearInterpolator li1 = new LinearInterpolator();
    li1.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li1.setUniform(n1,1.0,0.0,n2,1.0,0.0,r1);
    LinearInterpolator li2 = new LinearInterpolator();
    li2.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li2.setUniform(n1,1.0,0.0,n2,1.0,0.0,r2);
    for (int i2=0; i2<n2; ++i2) {
      double x2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        double x1 = i1;
        double s1i = s1[i2][i1];
        double s2i = s2[i2][i1];
        double s1p,s2p,ds1,ds2;
        do {
          double u1 = x1+s1i;
          double u2 = x2+s2i;
          s1p = s1i;
          s2p = s2i;
          s1i = li1.interpolate(u1,u2);
          s2i = li2.interpolate(u1,u2);
          ds1 = s1i-s1p;
          ds2 = s2i-s2p;
        } while (ds1*ds1+ds2*ds2>0.0001);
        s1[i2][i1] = (float)s1i;
        s2[i2][i1] = (float)s2i;
      }
    }
    return new float[][][]{s1,s2};
  }

  private static float[] smoothRandomCurve(
    int seed, double sigma, double smin, double smax, int n2) {
    Random r = new Random(seed);
    float[] s = new float[n2];
    for (int i2=0; i2<n2; ++i2) {
      if (r.nextFloat()>0.9f)
        s[i2] = r.nextFloat()-0.5f;
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply0(s,s);
    s = add((float)smin,mul(s,(float)(smax-smin)/(max(s)-min(s))));
    return s;
  }

  private static float[][] smoothRandomSurface(
    int seed, double sigma, double smin, double smax, int n2, int n3) {
    Random r = new Random(seed);
    float[][] s = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        if (r.nextFloat()>0.99f)
          s[i3][i2] = r.nextFloat()-0.5f;
      }
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply00(s,s);
    s = add((float)smin,mul(s,(float)(smax-smin)/(max(s)-min(s))));
    return s;
  }

  private static float[][] addAmplitude(double amin, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = copy(f);
    float[] a1 = smoothRandomCurve(140,20.0,amin,1.0,n2);
    float[] a2 = smoothRandomCurve(240,20.0,amin,1.0,n2);
    for (int i2=0; i2<n2; ++i2) {
      float fa = a1[i2];
      float da = (a2[i2]-fa)/n1;
      float[] a = rampfloat(fa,da,n1);
      mul(a,g[i2],g[i2]);
    }
    return g;
  }

  private static float[][][] addAmplitude(double amin, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] g = copy(f);
    float[][] a1 = smoothRandomSurface(140,20.0,amin,1.0,n2,n3);
    float[][] a2 = smoothRandomSurface(240,20.0,amin,1.0,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float fa = a1[i3][i2];
        float da = (a2[i3][i2]-fa)/n1;
        float[] a = rampfloat(fa,da,n1);
        mul(a,g[i3][i2],g[i3][i2]);
      }
    }
    return g;
  }

  private static float[][] addStructure(
    float smax, float fmax, float[][] f) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[] s1 = smoothRandomCurve(393,20.0,-smax,smax,n2);
    float[] s2 = smoothRandomCurve(394,20.0,-smax,smax,n2);
    s1 = limitShifts(smax,fmax,s1);
    s2 = limitShifts(smax,fmax,s2);
    float[][] g = new float[n2][n1];
    SincInterp si = new SincInterp();
    for (int i2=0; i2<n2; ++i2) {
      float fs = s1[i2];
      float ds = 1.0f+(s2[i2]-fs)/n1;;
      float[] s = rampfloat(fs,ds,n1);
      si.interpolate(n1,1.0,0.0,f[i2],n1,s,g[i2]);
    }
    return g;
  }
  private static float[] limitShifts(double smax, double fmax, float[] s) {
    int n2 = s.length;
    smax = min(smax,0.5f/(float)fmax);
    float sabs = 0.0f;
    for (int i2=1; i2<n2; ++i2) {
      float s00 = s[i2  ];
      float s01 = s[i2-1];
      float s2 = 0.5f*(s00-s01);
      if (abs(s2)>sabs) sabs = abs(s2);
    }
    if (sabs>smax)
      s = mul((float)smax/sabs,s);
    return s;
  }

  private static float[][][] addStructure(
    float smax, float fmax, float[][][] f) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][] s1 = smoothRandomSurface(393,20.0,-smax,smax,n2,n3);
    float[][] s2 = smoothRandomSurface(394,20.0,-smax,smax,n2,n3);
    s1 = limitShifts(smax,fmax,s1);
    s2 = limitShifts(smax,fmax,s2);
    float[][][] g = new float[n3][n2][n1];
    SincInterp si = new SincInterp();
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float fs = s1[i3][i2];
        float ds = 1.0f+(s2[i3][i2]-fs)/n1;;
        float[] s = rampfloat(fs,ds,n1);
        si.interpolate(n1,1.0,0.0,f[i3][i2],n1,s,g[i3][i2]);
      }
    }
    return g;
  }
  private static float[][] limitShifts(double smax, double fmax, float[][] s) {
    int n2 = s[0].length;
    int n3 = s.length;
    smax = min(smax,0.5f/(float)fmax);
    float sabs = 0.0f;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        float s00 = s[i3  ][i2  ];
        float s01 = s[i3  ][i2-1];
        float s10 = s[i3-1][i2  ];
        float s11 = s[i3-1][i2-1];
        float s2 = 0.5f*(s00-s01+s10-s11);
        float s3 = 0.5f*(s00-s10+s01-s11);
        if (abs(s2)>sabs) sabs = abs(s2);
        if (abs(s3)>sabs) sabs = abs(s3);
      }
    }
    if (sabs>smax)
      s = mul((float)smax/sabs,s);
    return s;
  }

  private static float[][] addFault(float smax, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = new float[n2][n1];
    float c1 = 0.0f;
    float c2 = -n2;
    float d2 = n2/2-c2;
    float rs = d2*d2;
    SincInterp si = new SincInterp();
    for (int i2=0; i2<n2; ++i2) {
      d2 = i2-c2;
      float ds = d2*d2-rs;
      if (ds<0.0) {
        float[] x1 = rampfloat(0.0f,1.0f-smax/n1,n1);
        si.interpolate(n1,1.0,0.0,f[i2],n1,x1,g[i2]);
      } else {
        copy(f[i2],g[i2]);
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      d2 = i2-c2;
      for (int i1=0; i1<n1; ++i1) {
        float d1 = i1-c1;
        float ds = d1*d1+d2*d2-rs;
        if (ds>0.0)
          g[i2][i1] = f[i2][i1];
      }
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.applyX0(g,g);
    return g;
  }

  private static float[][][] addFault(float smax, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] g = new float[n3][n2][n1];
    float c1 = 0.0f;
    float c2 = -n2;
    float c3 = -n3*3;
    float d2 = n2/2-c2;
    float d3 = n3/2-c3;
    float rs = d2*d2+d3*d3;
    SincInterp si = new SincInterp();
    for (int i3=0; i3<n3; ++i3) {
      d3 = i3-c3;
      for (int i2=0; i2<n2; ++i2) {
        d2 = i2-c2;
        float ds = d2*d2+d3*d3-rs;
        if (ds<0.0) {
          float[] x1 = rampfloat(0.0f,1.0f-smax/n1,n1);
          si.interpolate(n1,1.0,0.0,f[i3][i2],n1,x1,g[i3][i2]);
        } else {
          copy(f[i3][i2],g[i3][i2]);
        }
      }
    }
    for (int i3=0; i3<n3; ++i3) {
      d3 = i3-c3;
      for (int i2=0; i2<n2; ++i2) {
        d2 = i2-c2;
        for (int i1=0; i1<n1; ++i1) {
          float d1 = i1-c1;
          float ds = d1*d1+d2*d2+d3*d3-rs;
          if (ds>0.0)
            g[i3][i2][i1] = f[i3][i2][i1];
        }
      }
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.applyX0X(g,g);
    rgf.applyXX0(g,g);
    return g;
  }

  private static float[][] addErosion(int kerosion, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = makeEvents(kerosion,n2);
    float[][] h = copy(f);
    copy(kerosion,n2,g,h);
    return h;
  }

  private static float[][][] addErosion(int kerosion, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] g = makeEvents(kerosion,n2,n3);
    float[][][] h = copy(f);
    copy(kerosion,n2,n3,g,h);
    return h;
  }

  private static float[] makeEvents(int n1) {
    Random r = new Random(31415);
    return pow(mul(2.0f,sub(randfloat(r,n1),0.5f)),7.0f);
  }
  private static float[][] makeEvents(int n1, int n2) {
    float[][] f = new float[n2][n1];
    f[0] = makeEvents(n1);
    for (int i2=0; i2<n2; ++i2)
      copy(f[0],f[i2]);
    return f;
  }
  private static float[][][] makeEvents(int n1, int n2, int n3) {
    float[][][] f = new float[n3][n2][n1];
    f[0] = makeEvents(n1,n2);
    for (int i3=0; i3<n3; ++i3)
      copy(f[0],f[i3]);
    return f;
  }

  private static float[] addRickerWavelet(double fpeak, float[] f) {
    int n1 = f.length;
    int ih = (int)(3.0/fpeak);
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0; jh<nh; ++jh)
      h[jh] = ricker(fpeak,jh-ih);
    float[] g = new float[n1];
    Conv.conv(nh,-ih,h,n1,0,f,n1,0,g);
    return g;
  }
  private static float[][] addRickerWavelet(double fpeak, float[][] f) {
    int n2 = f.length;
    float[][] g = new float[n2][];
    for (int i2=0; i2<n2; ++i2)
      g[i2] = addRickerWavelet(fpeak,f[i2]);
    return g;
  }
  private static float[][][] addRickerWavelet(double fpeak, float[][][] f) {
    int n3 = f.length;
    float[][][] g = new float[n3][][];
    for (int i3=0; i3<n3; ++i3)
      g[i3] = addRickerWavelet(fpeak,f[i3]);
    return g;
  }
  private static float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }

  private static float[][] addNoise(float nrms, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    Random r = new Random(31415);
    //nrms *= max(abs(f));
    float[][] g = mul(2.0f,sub(randfloat(r,n1,n2),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply10(g,g); // 1st derivative enhances high-frequencies
    float frms = sqrt(sum(mul(f,f))/n1/n2);
    float grms = sqrt(sum(mul(g,g))/n1/n2);
    g = mul(g,nrms*frms/grms);
    return add(f,g);
  }

  private static float[][][] addNoise(float nrms, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    Random r = new Random(31415);
    nrms *= max(abs(f));
    float[][][] g = mul(2.0f,sub(randfloat(r,n1,n2,n3),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply100(g,g); // 1st derivative enhances high-frequencies
    float frms = sqrt(sum(mul(f,f))/n1/n2/n3);
    float grms = sqrt(sum(mul(g,g))/n1/n2/n3);
    g = mul(g,nrms*frms/grms);
    return add(f,g);
  }

  /**
   * A 1D coordinate mapping, with derivative.
   */
  private interface F1 {
    /** Returns the value f(x). */
    public float f(float x);
    /** Returns the derivative df/dx(x). */
    public float df(float x);
  }

  /**
   * A 2D coordinate mapping, with derivatives (Jacobian).
   */
  private interface F2 {
    /** Returns the value f1(x1,x2). */
    public float f1(float x1, float x2);
    /** Returns the value f2(x1,x2). */
    public float f2(float x1, float x2);
    /** Returns the partial derivative df1/dx1(x1,x2). */
    public float df11(float x1, float x2);
    /** Returns the partial derivative df1/dx2(x1,x2). */
    public float df12(float x1, float x2);
    /** Returns the partial derivative df2/dx1(x1,x2). */
    public float df21(float x1, float x2);
    /** Returns the partial derivative df2/dx2(x1,x2). */
    public float df22(float x1, float x2);
  }

  /**
   * Applies a 2D coordinate transform to an image.
   * @param f coordinate transform f(x).
   * @param p input image p(x).
   * @return transformed image q(x) = p(f(x)).
   */
  private static float[][] apply(F2 f, float[][] p) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[][] q = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float f1 = f.f1(i1,i2);
        float f2 = f.f2(i1,i2);
        q[i2][i1] = _si.interpolate(n2,1.0,0.0,n1,1.0,0.0,p,f1,f2);
      }
    }
    return q;
  }
  /**
   * Applies a 2D coordinate transform to an image and normal vectors.
   * The input array {p0,p1,p2} contains the input image p0 and 1st and 2nd
   * components of normal vectors, p1 and p2. The returned array {q0,q1,q2}
   * contains the corresponding transformed image and normal vectors. All
   * normal vectors are unit vectors.
   * @param f coordinate transform f(x).
   * @param p input image p(x) and normal vectors.
   * @return transformed image q(x) = p(f(x) and normal vectors.
   */
  private static float[][][] apply(F2 f, float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] q0 = apply(f,p[0]);
    float[][] q1 = apply(f,p[1]);
    float[][] q2 = apply(f,p[2]);
    for (int i2=0; i2<n2; ++i2) {
      float x2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        float x1 = i1;
        float q1i = f.df11(x1,x2)*q1[i2][i1]+f.df21(x1,x2)*q2[i2][i1];
        float q2i = f.df12(x1,x2)*q1[i2][i1]+f.df22(x1,x2)*q2[i2][i1];
        float qsi = 1.0f/sqrt(q1i*q1i+q2i*q2i);
        q1[i2][i1] = q1i*qsi;
        q2[i2][i1] = q2i*qsi;
      }
    }
    return new float[][][]{q0,q1,q2};
  }

  /**
   * A linear 1D coordinate mapping f(x) = a+b*x.
   */
  private static class Linear1 implements F1 {

    /**
     * Constructs a linear 1D coordinate mapping.
     * @param a the intercept.
     * @param b the slope.
     */
    public Linear1(float a, float b) {
      _a = a;
      _b = b;
    }
    public float f(float x) {
      return _a+_b*x;
    }
    public float df(float x) {
      return _b;
    }
    private float _a; // intercept
    private float _b; // slope
  }

  /**
   * A linear fault in a 2D image.
   */
  private static class LinearFault2 implements F2 {

    /**
     * Constructs a linear fault.
     * @param ftrace coordinate x2 of the fault where x1 = 0.
     * @param ftheta fault dip, measured in degrees from vertical.
     * @param fthrow fault throw, a function of coordinate x1.
     */
    public LinearFault2(float ftrace, float ftheta, F1 fthrow) {
      float rtheta = toRadians(ftheta);
      float ctheta = cos(rtheta);
      float stheta = sin(rtheta);
      _a0 = ftrace*ctheta;
      _a1 = stheta;
      _a2 = -ctheta;
      _d1 = fthrow;
      if (_a1/_a2>0.0f) {
        _a0 = -_a0;
        _a1 = -_a1;
        _a2 = -_a2;
      }
    }
    public float f1(float x1, float x2) {
      float f = x1;
      if (faulted(x1,x2))
        f -= _d1.f(x1);
      return f;
    }
    public float f2(float x1, float x2) {
      float f = x2;
      if (_a0+_a1*x1+_a2*x2<=0.0f)
        f += _d1.f(x1)*_a1/_a2;
      return f;
    }
    public float df11(float x1, float x2) {
      float d = 1.0f;
      if (faulted(x1,x2))
        d -= _d1.df(x1);
      return d;
    }
    public float df12(float x1, float x2) {
      return 0.0f;
    }
    public float df21(float x1, float x2) {
      float d = 0.0f;
      if (faulted(x1,x2))
        d += _d1.df(x1)*_a1/_a2;
      return d;
    }
    public float df22(float x1, float x2) {
      return 1.0f;
    }
    private float _a0,_a1,_a2;
    private F1 _d1;
    private boolean faulted(float x1, float x2) {
      return _a0+_a1*x1+_a2*x2<=0.0f;
    }
  }

  /**
   * Vertical shear of a 2D image.
   */
  private static class VerticalShear2 implements F2 {
    public VerticalShear2(Linear1 s1) {
      _s1 = s1;
    }
    public float f1(float x1, float x2) {
      return x1-_s1.f(x2);
    }
    public float f2(float x1, float x2) {
      return x2;
    }
    public float df11(float x1, float x2) {
      return 1.0f;
    }
    public float df12(float x1, float x2) {
      return -_s1.df(x2);
    }
    public float df21(float x1, float x2) {
      return 0.0f;
    }
    public float df22(float x1, float x2) {
      return 1.0f;
    }
    private Linear1 _s1;
  }

  /**
   * Sinusoidal folding in a 2D image.
   */
  private static class Sinusoidal2 implements F2 {
    public Sinusoidal2(float a1, float b1, float a2, float b2) {
      _a1 = a1;
      _b1 = b1;
      _a2 = a2;
      _b2 = b2;
    }
    public float f1(float x1, float x2) {
      return x1-(_a1+_b1*x1)*sin((_a2+_b2*x2)*x2);
    }
    public float f2(float x1, float x2) {
      return x2;
    }
    public float df11(float x1, float x2) {
      return 1.0f-_b1*sin((_a2+_b2*x2)*x2);
    }
    public float df12(float x1, float x2) {
      return -(_a1+_b1*x1)*(_a2+2.0f*_b2*x2)*cos((_a2+_b2*x2)*x2);
    }
    public float df21(float x1, float x2) {
      return 0.0f;
    }
    public float df22(float x1, float x2) {
      return 1.0f;
    }
    private float _a1,_b1,_a2,_b2;
  }

  private static float[][][] makeReflectivityWithNormals(int n1, int n2) {
    Random random = new Random(31);
    float[] r = pow(mul(2.0f,sub(randfloat(random,n1),0.5f)),5.0f);
    float[][][] p = new float[3][n2][n1];
    for (int i2=0; i2<n2; ++i2)
      copy(r,p[0][i2]);
    p[1] = fillfloat(1.0f,n1,n2);
    p[2] = fillfloat(0.0f,n1,n2);
    return p;
  }

  private static float[][][] combine(
    float depth, float[][][] pa, float[][][] pb) 
  {
    int n1 = pa[0][0].length;
    int n2 = pa[0].length;
    int nc = pa.length;
    float[][][] pc = new float[nc][n2][n1];
    for (int ic=0; ic<nc; ++ic) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          pc[ic][i2][i1] = (i1<depth)?pa[ic][i2][i1]:pb[ic][i2][i1];
        }
      }
    }
    return pc;
  }

  private static float[][][] addWavelet(double fpeak, float[][][] p) {
    double sigma = max(1.0,1.0/(2.0*PI*fpeak));
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] p0 = p[0];
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float[][] q = copy(p0);
    float[][] q1 = new float[n2][n1];
    float[][] q2 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    for (int id=0; id<2; ++id) { // 2nd directional derivative of Gaussian
      rgf.apply10(q,q1);
      rgf.apply01(q,q2);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          q[i2][i1] = p1[i2][i1]*q1[i2][i1]+p2[i2][i1]*q2[i2][i1];
        }
      }
    }
    q = mul(q,-1.0f/sqrt(sum(mul(q,q))/n1/n2));
    return new float[][][]{q,p1,p2};
  }
}
