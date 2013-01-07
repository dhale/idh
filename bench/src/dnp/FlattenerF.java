/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.util.logging.Logger;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates and applies shifts to flatten features in 2D and 3D images.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.03
 */
public class FlattenerF {

  public FlattenerF() {
    this(8.0,0.1);
  }

  public FlattenerF(double sigma, double epsilon) {
    _sigma = (float)sigma;
    _epsilon = (float)epsilon;
  }

  public float[][] findShifts(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;

    // Estimate slopes.
    float[][] p2 = new float[n2][n1];
    findSlopes(_sigma,f,p2);
    //fill(-0.04f,p2);
    float a2 = sum(p2)/n1/n2;

    // Frequency sampling.
    float pi = FLT_PI;
    float twopi = 2.0f*pi;
    int npad = 100;
    int n1fft = FftReal.nfftFast(n1+npad);
    int n2fft = FftComplex.nfftFast(n2+npad);
    int nk1 = n1fft/2+1;
    int nk2 = n2fft;
    float dk1 = twopi/n1fft;
    float dk2 = twopi/n2fft;

    // Pad image for FFT.
    float[][] p2fft = new float[nk2][2*nk1];
    //pad(n1fft,n2fft,p2,p2fft);
    copy(n1,n2,p2,p2fft);
    p2 = null;

    // Forward FFT p2(x1,x2) => p2(k1,k2).
    FftReal fft1 = new FftReal(n1fft);
    FftComplex fft2 = new FftComplex(n2fft);
    fft1.realToComplex1(-1,n2fft,p2fft,p2fft);
    fft2.complexToComplex2(-1,nk1,p2fft,p2fft);

    // Compute shifts s(k1,k2) from p2(k1,k2).
    float es = _epsilon*_epsilon;
    float ek1 = 0.01f*dk1;
    float ek2 = 0.01f*dk2;
    float tiny = es*ek1*ek1+ek2*ek2;
    float[][] sfft = p2fft;
    float scale = 1.0f/n1fft/n2fft;
    for (int ik2=0; ik2<nk2; ++ik2) {
      float k2 = ik2*dk2;
      if (k2>pi) k2 -= twopi;
      float k2s = k2*k2;
      for (int ik1=0,ikr=0,iki=1; ik1<nk1; ++ik1,ikr+=2,iki+=2) {
        float k1 = ik1*dk1;
        float k1s = k1*k1;
        float p2r = p2fft[ik2][ikr];
        float p2i = p2fft[ik2][iki];
        float s2 = scale*k2/(es*k1s+k2s+tiny);
        sfft[ik2][ikr] = -s2*p2i;
        sfft[ik2][iki] =  s2*p2r;
      }
    }

    // Inverse FFT s(k1,k2) => s(x1,x2).
    fft2.complexToComplex2(1,nk1,sfft,sfft);
    fft1.complexToReal1(1,nk2,sfft,sfft);
    float[][] s = copy(n1,n2,sfft);
    
    float avg = 0.0f;
    float ds1 = 0.0f;
    float ds2 = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        avg += s[i2][i1];
        if (i1>0) ds1 += s[i2][i1]-s[i2][i1-1];
        if (i2>0) ds2 += s[i2][i1]-s[i2-1][i1];
      }
    }
    avg = avg/n1/n2;
    ds1 = ds1/(n1-1)/n2;
    ds2 = ds2/n1/(n2-1);
    float b0 = avg;
    float b1 = ds1;
    float b2 = ds2+a2;
    System.out.println("b0="+b0+" b1="+b1+" b2="+b2);
    for (int i2=0,j2=n2/2; i2<n2; ++i2) {
      for (int i1=0,j1=n1/2; i1<n1; ++i1) {
        s[i2][i1] -= b0+b1*(i1-j1)+b2*(i2-j2);
      }
    }

    return s;
  }

  public float[][] applyShifts(float[][] f, float[][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    SincInterp si = new SincInterp();
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      sub(r,s[i2],t);
      si.interpolate(n1,1.0,0.0,f[i2],n1,t,g[i2]);
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // smoothing half-width in 1st dimension for slopes
  private float _epsilon; // smoothing in 1st dimension for shifts
 
  private static void findSlopes(float sigma, float[][] f, float[][] p2) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = p2;
    float sigma1 = sigma;
    float sigma2 = 1.0f;
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    lof.applyForNormal(f,u1,u2);
    float p2min = -10.0f;
    float p2max =  10.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        if (u2i<0.0f) {
          p2[i2][i1] = (u2i>p2min*u1i)?-u2i/u1i:p2min;
        } else {
          p2[i2][i1] = (u2i<p2max*u1i)?-u2i/u1i:p2max;
        }
      }
    }
    ZeroMask zm = new ZeroMask(0.1,sigma,1,f);
    zm.apply(0.0f,p2);
  }

  private static void pad(int n1pad, int n2pad, float[][] p, float[][] q) {
    int n1 = p[0].length;
    int n2 = p.length;
    copy(n1,n2,p,q);
    float s1 = 1.0f/(n1pad-n1+1);
    float s2 = 1.0f/(n2pad-n2+1);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=n1; i1<n1pad; ++i1) {
        float w = (i1-n1+1)*s1;
          q[i2][i1] = w*q[i2][0]+(1.0f-w)*q[i2][n1-1];
      }
    }
    for (int i2=n2; i2<n2pad; ++i2) {
      float w = (i2-n2+1)*s2;
      for (int i1=0; i1<n1pad; ++i1) {
        q[i2][i1] = w*q[0][i1]+(1.0f-w)*q[n2-1][i1];
      }
    }
  }

  private static void simplify(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    double a00 = 0.0;
    double a01 = 0.0;
    double a02 = 0.0;
    double a11 = 0.0;
    double a12 = 0.0;
    double a22 = 0.0;
    double b0 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a00 += 1.0;
        a01 += i1;
        a02 += i2;
        a11 += i1*i1;
        a12 += i1*i2;
        a22 += i2*i2;
        b0 += s[i2][i1];
        b1 += i1*s[i2][i1];
        b2 += i2*s[i2][i1];
      }
    }
    double[][] a = {
      {a00,a01,a02},
      {a01,a11,a12},
      {a02,a12,a22}
    };
    DMatrix am = new DMatrix(new double[][]{
      {a00,a01,a02},
      {a01,a11,a12},
      {a02,a12,a22}
    });
    DMatrix bm = new DMatrix(new double[][]{{b0},{b1},{b2}});
    DMatrixLud lud = new DMatrixLud(am);
    am = lud.solve(bm);
    float a0 = (float)am.get(0,0);
    float a1 = (float)am.get(1,0);
    float a2 = (float)am.get(2,0);
    System.out.println("a0="+a0+" a1="+a1+" a2="+a2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        s[i2][i1] -= a0+a1*i1+a2*i2;
      }
    }
  }
} 
