/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

import dnp.LocalSlopeFinder;
import het.RecursiveExponentialFilter;

// FOR DEVELOPMENT ONLY
import edu.mines.jtk.sgl.*;

/*
Algorithm:
given image f
find slopes p2,p3
use slopes p2,p3 to align images fm and fp
compute correlation products cp = {cmp,cmm,cpp}
Fourier transform products cp
for all phis
  for all thetas
    smooth and inverse transform products cp
    c = (<fmp>*<fmp>)/(<fmm><fpp>), if <fmp> > 0; 0, otherwise
    map c to fault likelihood c = 1-c
    remember c max and corresponding phi and theta max
thin cmax by smoothing laterally and picking peaks
use p2,p3,cmax for structure-oriented smoothing
for all fault surfaces
  gather image samples alongside fault
  cross-correlate to find displacements
*/

/**
 * Finds faults in 3D seismic images.
 * <p>
 * This version computes fault likelihoods from correlation coefficients,
 * which are ratios of smoothed correlation products, and it uses rotated
 * Gaussian filters implemented with FFTs to perform that smoothing.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.06.24
 */
public class FaultFinder3B {

  /**
   * Constructs a fault finder with specified parameters.
   * @param slopeMax maximum slope of seismic reflections.
   * @param shiftMax maximum fault shift, in samples.
   * @param thetaMax maximum fault angle theta, in degrees.
   */
  public FaultFinder3B(double slopeMax, double shiftMax, double thetaMax) {
    _slopeMax = (float)slopeMax;
    _shiftMax = (float)shiftMax;
    _thetaMax = (float)thetaMax;
    _faultLengthMin = 2.0f*_shiftMax;
    _sigma = 2.0f*_shiftMax;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    _st  = makeThetaSampling(-_thetaMax,_thetaMax);
    _sp = makePhiSampling(-90.0,90.0);
  }

  public void setPhiSampling(Sampling sp) {
    _sp = sp;
  }

  public void setThetaSampling(Sampling st) {
    _st = st;
  }

  public Sampling makePhiSampling(double phiMin, double phiMax) {
    return angleSampling(SIGMA2C,phiMin,phiMax);
  }

  public Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma,thetaMin,thetaMax);
  }

  public float[][][][] findSlopes(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(SIGMA1P,_slopeMax);
    lsf.findSlopes(f,p2,p3,null);
    return new float[][][][]{p2,p3};
  }

  public static float[][][] taper(int m, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] g = copy(f);
    float[] t = new float[m];
    for (int i=0; i<m; ++i) {
      t[i] = (float)(0.54+0.46*cos(PI*(m-i)/m));
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0,j1=n1-1; i1<m; ++i1,--j1) {
          float ti = t[i1];
          g[i3][i2][i1] *= ti;
          g[i3][i2][j1] *= ti;
        }
      }
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0,j2=n2-1; i2<m; ++i2,--j2) {
        float ti = t[i2];
        for (int i1=0; i1<n1; ++i1) {
          g[i3][i2][i1] *= ti;
          g[i3][j2][i1] *= ti;
        }
      }
    }
    for (int i3=0,j3=n3-1; i3<m; ++i3,--j3) {
      float ti = t[i3];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          g[i3][i2][i1] *= ti;
          g[j3][i2][i1] *= ti;
        }
      }
    }
    return g;
  }

  public float[][][][] align(float[][][] p2, float[][][] p3, float[][][] f) {
    /*
    int n1 = n1(f);
    int n2 = n2(f);
    int n3 = n3(f);
    float[][][] fm = new float[n3][n2][];
    float[][][] fp = new float[n3][n2][];
    float[] xm = new float[n1];
    float[] xp = new float[n1];
    float phir = (float)toRadians(phi);
    float cosp = cos(phir);
    float sinp = sin(phir);
    _si.setUniformSampling(n1,1.0,0.0);
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        float[] p2m = p2[i3m][i2];
        float[] p2p = p2[i3p][i2];
        float[] p3m = p3[i3m][i2];
        float[] p3p = p3[i3p][i2];
        float[] f2m = f[i3m][i2];
        float[] f2p = f[i3p][i2];
        if (p2m!=null && p2p!=null && 
            p3m!=null && p3p!=null &&
            f2m!=null && f2p!=null) {
          for (int i1=0; i1<n1; ++i1) {
            // for phi =   0, p3 =  p3 (no change)
            // for phi =  90, p3 = -p2
            // for phi = -90, p3 =  p2
            float p3mi = cosp*p3m[i1]-sinp*p2m[i1];
            float p3pi = cosp*p3p[i1]-sinp*p2p[i1];
            xm[i1] = i1-p3mi;
            xp[i1] = i1+p3pi;
          }
          fm[i3][i2] = new float[n1];
          fp[i3][i2] = new float[n1];
          _si.setUniformSamples(f2m);
          _si.interpolate(n1,xm,fm[i3][i2]);
          _si.setUniformSamples(f2p);
          _si.interpolate(n1,xp,fp[i3][i2]);
        }
      }
    }
    return new float[][][][]{fm,fp};
    */
    return null;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final double SIGMA1P = 8.0; // for estimating slopes
  private static final double SIGMA2C = 4.0; // for correlation coeff

  private float _sigma; // for smoothing correlation products vertically 
  private float _slopeMax; // maximum slope used in alignment
  private float _thetaMax; // maximum fault angle theta
  private float _shiftMax; // maximum vertical throw for any fault
  private float _faultLengthMin; // min number of samples in a fault
  private SincInterpolator _si; // for shifting, squeezing and stretching
  private Sampling _st; // sampling of fault angle theta
  private Sampling _sp; // sampling of fault strike phi

  // Gaussian smoothing for correlation products.
  private static class Smoother {
    public Smoother(
      double sigmau, double sigmav, double sigmaw,
      float[][][][] f)
    {
      _sigmau = (float)sigmau;
      _sigmav = (float)sigmav;
      _sigmaw = (float)sigmaw;
      _n1 = f[0][0][0].length;
      _n2 = f[0][0].length;
      _n3 = f[0].length;
      _nf = f.length;
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
      _fk = new float[_nf][_nk3][_n2][_n1*2];
      final float[][][][] fx = f;
      final float[][][][] fk = _fk;
      loop(_n2,new LoopInt() {
      public void compute(int i2) {
        float[][] fxpad = new float[_n3fft][_n1];
        float[][] fkpad = new float[_nk3][_n1*2];
        for (int i=0; i<_nf; ++i) {
          for (int i3=0; i3<_n3; ++i3)
            copy(fx[i][i3][i2],fxpad[i3]);
          _fft3.realToComplex2(-1,_n1,fxpad,fkpad);
          for (int i3=0; i3<_nk3; ++i3)
            ccopy(fkpad[i3],fk[i][i3][i2]);
        }
      }});
    }
    public void apply(double phi, double theta, float[][][][] g) {
      float[][][] h = makeFilter(phi,theta);
      final float[][][] gk = new float[_nk3][_n2][_n1*2];
      final float[][][][] gx = g;
      for (int i=0; i<_nf; ++i) {
        applyFilter(h,_fk[i],gk);
        loop(_n2,new LoopInt() {
        public void compute(int i2) {
          float[][] gkpad = new float[_nk3][_n1*2];
          float[][] gxpad = new float[_n3fft][_n1];
          for (int i=0; i<_nf; ++i) {
            for (int i3=0; i3<_nk3; ++i3)
              ccopy(gk[i3][i2],gkpad[i3]);
            _fft3.complexToReal2(1,_n1,gkpad,gxpad);
            for (int i3=0; i3<_n3; ++i3)
              copy(gxpad[i3],gx[i][i3][i2]);
          }
        }});
      }
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
      final float u1 =   ct, u2 = -st*sp, u3 = st*cp; // u down the fault
      final float v1 = 0.0f, v2 =     cp, v3 =    sp; // v along strike
      final float w1 =  -st, w2 = -ct*sp, w3 = ct*cp; // w normal to fault
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
        for (int i2=0; i2<_nk2; ++i2) {
          float k2 = i2*dk2;
          if (i2*2>_nk2) k2 -= twopi;
          float[] h32 = h[i3][i2];
          for (int i1=0; i1<_nk1; ++i1) {
            float k1 = i1*dk1;
            if (i1*2>_nk1) k1 -= twopi;
            float uk = u1*k1+u2*k2+u3*k3;
            float vk = v1*k1+v2*k2+v3*k3;
            float wk = w1*k1+w2*k2+w3*k3;
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
  }

  // Samples on opposite sides of a fault.
  private static class FaultSamples {
    public static final int KA = 2; // offset used to find fault shifts
    public static final int KB = KA+1; // offset used with KA for slopes
    public int[] k1,k2; // sample indices of fault
    public float[] fmb,fma,fpa,fpb; // image samples at -ka-1, -ka, ka, ka+1
    public FaultSamples(
      int[] k1, int[] k2,
      float[] fmb, float[] fma, float[] fpa, float[] fpb)
    {
      this.k1 = k1;
      this.k2 = k2;
      this.fmb = fmb;
      this.fma = fma;
      this.fpa = fpa;
      this.fpb = fpb;
    }
  }

  private FaultSamples findFaultSamples(
    int i1, int i2, float[][] c, float[][] f, boolean[][] b) 
  {
    if (c[i2][i1]==0.0f)
      return null;
    int n1 = c[0].length;
    int n2 = c.length;
    int k2min = 0;
    int k2max = n2-1;
    int[] k2s = new int[n1];
    int n = 0;
    boolean done = false;
    for (int k1=i1,k2=i2; !done && k1<n1; ++k1) {
      if (c[k2][k1]>0.0f) {
        b[k2][k1] = true;
        k2s[n++] = k2;
      } else if (k2<k2max && c[k2+1][k1]>0.0f) {
        ++k2;
        b[k2][k1] = true;
        k2s[n++] = k2;
      } else if (k2>k2min && c[k2-1][k1]>0.0f) {
        --k2;
        b[k2][k1] = true;
        k2s[n++] = k2;
      } else {
        done = true;
      }
    }
    if (n<_faultLengthMin)
      return null;
    int[] k1 = new int[n];
    int[] k2 = new int[n];
    float[] fmb = new float[n];
    float[] fma = new float[n];
    float[] fpa = new float[n];
    float[] fpb = new float[n];
    int ka = FaultSamples.KA;
    int kb = FaultSamples.KB;
    for (int k=0; k<n; ++k) {
      int k1k = i1+k;
      int k2k = k2s[k];
      int kmb = max(k2min,k2k-kb);
      int kma = max(k2min,k2k-ka);
      int kpa = min(k2max,k2k+ka);
      int kpb = min(k2max,k2k+kb);
      fmb[k] = f[kmb][k1k];
      fma[k] = f[kma][k1k];
      fpa[k] = f[kpa][k1k];
      fpb[k] = f[kpb][k1k];
      k1[k] = k1k;
      k2[k] = k2k;
    }
    return new FaultSamples(k1,k2,fmb,fma,fpa,fpb);
  }

  private Sampling angleSampling(double sigma, double amin, double amax) {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
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
    float[][][] f = zerofloat(n1,n2,n3);
    float[][][] g = zerofloat(n1,n2,n3);
    f[n3/2][n2/2][n1/2] = 1.0f;
    float[][][][] fs = {f};
    float[][][][] gs = {g};
    double sigmau = 30.0;
    double sigmav = 15.0;
    double sigmaw = 1.0;
    double phi = 45.0;
    double theta = 45.0;
    Smoother s = new Smoother(sigmau,sigmav,sigmaw,fs);
    s.apply(phi,theta,gs);
    SimpleFrame.asImagePanels(g);
  }
}
