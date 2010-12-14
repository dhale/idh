/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package util;

import java.util.Random;

import edu.mines.jtk.dsp.*;
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
  public static void main(String[] args) {
    if (args.length==0 || args[0]=="seismic3d2010A") {
      float[][][] f = seismic3d2010A();
      SimpleFrame frame = new SimpleFrame();
      frame.addImagePanels(f);
      frame.getOrbitView().setScale(2.0);
      frame.setSize(900,900);
    } else {
      System.out.println("unrecognized type of fake data");
    }
  }

  /**
   * Returns a fake 3D seismic image, version 2010A. The image includes 
   * default structure, a fault, an unconformity, amplitude variations,
   * and random noise.
   * @return the fake 3d seismic image.
   */
  public static float[][][] seismic3d2010A() {
    return seismic3d2010A(201,201,201,20.0,10.0,30.0,0.5,0.5);
  }

  /**
   * Returns a fake 3D seismic image, version A. As options, this image 
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
   * @return the fake 3d seismic image.
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

  ///////////////////////////////////////////////////////////////////////////
  // private

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
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float fs = s1[i3][i2];
        float ds = 1.0f+(s2[i3][i2]-fs)/n1;;
        float[] s = rampfloat(fs,ds,n1);
        si.setUniformSamples(f[i3][i2]);
        si.interpolate(n1,s,g[i3][i2]);
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
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    for (int i3=0; i3<n3; ++i3) {
      d3 = i3-c3;
      for (int i2=0; i2<n2; ++i2) {
        d2 = i2-c2;
        float ds = d2*d2+d3*d3-rs;
        if (ds<0.0) {
          float[] x1 = rampfloat(0.0f,1.0f-smax/n1,n1);
          si.setUniformSamples(f[i3][i2]);
          si.interpolate(n1,x1,g[i3][i2]);
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

  private static float[][][] addErosion(
    int kerosion, float[][][] f) 
  {
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
    return pow(mul(2.0f,sub(randfloat(r,n1),0.5f)),8.0f);
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
}
