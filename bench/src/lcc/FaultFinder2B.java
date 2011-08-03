/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

import dnp.LocalSlopeFinder;
import het.RecursiveExponentialFilter;

// FOR DEVELOPMENT ONLY
import java.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * Finds faults in 2D seismic images.
 * <p>
 * This version computes fault likelihoods from correlation coefficients,
 * which are ratios of smoothed correlation products, and it uses rotated
 * Gaussian filters implemented with FFTs to perform that smoothing.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.06.24
 */
public class FaultFinder2B {

  /**
   * Constructs a fault finder with specified parameters.
   * @param slopeMax maximum slope of seismic reflections.
   * @param shiftMax maximum fault shift, in samples.
   * @param thetaMax maximum fault angle theta, in degrees.
   */
  public FaultFinder2B(double slopeMax, double shiftMax, double thetaMax) {
    _slopeMax = (float)slopeMax;
    _shiftMax = (float)shiftMax;
    _thetaMax = (float)thetaMax;
    _faultLengthMin = 2.0f*_shiftMax;
    _sigma = 2.0f*_shiftMax;
    _rgf1 = new RecursiveGaussianFilter(1.0f);
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    _st  = angleSampling(-_thetaMax,_thetaMax);
  }

  public Sampling angleSampling(double angleMin, double angleMax) {
    double fa = angleMin;
    double da = toDegrees(0.5/_sigma);
    int na = 3+(int)((angleMax-angleMin)/da);
    da = (angleMax-angleMin)/(na-1);
    return new Sampling(na,da,fa);
  }

  public float[][] findSlopes(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] p = new float[n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(SIGMA1_SLOPE,_slopeMax);
    lsf.findSlopes(f,p);
    return p;
  }
  private static final double SIGMA1_SLOPE = 8.0;

  public float[][][] findThetas(float[][] p, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][][] fmp = align(1,p,f);
    float[][][] ct = ctScan(fmp);
    return ctThin(ct);
  }

  public float[][][] ctThin(float[][][] ct) {
    float[][] c = ct[0];
    float[][] t = ct[1];
    int n1 = c[0].length;
    int n2 = c.length;
    //pow(c,4.00f,c);
    //_rgf1.apply00(c,c);
    //pow(c,0.25f,c);
    float dt = (float)_st.getDelta();
    float ft = (float)_st.getFirst();
    float lt = (float)_st.getLast();
    float tmin = ft-0.5f*dt;
    float tmax = lt+0.5f*dt;
    float[][] cc = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
      float[] ti = t[i2  ];
      float[] ci = c[i2  ];
      float[] cm = c[i2-1];
      float[] cp = c[i2+1];
      for (int i1=0; i1<n1; ++i1) {
        float cii = ci[i1];
        float tii = ti[i1];
        if (cm[i1]<cii && cp[i1]<cii && tmin<tii && tii<tmax) {
          cc[i2][i1] = cii;
          tt[i2][i1] = tii;
        }
      }
    }
    return new float[][][]{cc,tt};
  }

  public float[][][] ctScan(float[][][] fmp) {
    int n1 = fmp[0][0].length;
    int n2 = fmp[0].length;
    int nt = _st.getCount();
    float dt = (float)_st.getDelta();
    float ft = (float)_st.getFirst();
    float[][] c = new float[n2][n1];
    float[][] t = new float[n2][n1];
    for (int it=0; it<nt; ++it) {
      float ti = ft+it*dt;
      float[][] cx = xcor(ti,fmp);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float ci = 1.0f-cx[i2][i1];
          if (ci>c[i2][i1]) {
            c[i2][i1] = ci;
            t[i2][i1] = ti;
          }
        }
      }
    }
    return new float[][][]{c,t};
  }

  public float[][] smooth(
    double sigma, float[][] p, float[][] c, float[][] f) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    EigenTensors2 d = new EigenTensors2(n1,n2);
    d.setEigenvalues(0.001f,1.00f);
    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        s[i2][i1] = 1.0f-pow(c[i2][i1],0.1f);
        //s[i2][i1] = (c[i2][i1]>0.0)?0.0f:1.0f;
        float dip = atan(p[i2][i1]);
        float u1 =  cos(dip);
        float u2 = -sin(dip);
        d.setEigenvectorU(i1,i2,u1,u2);
      }
    }
    float a = (float)(0.5*sigma*sigma);
    float[][] g = new float[n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,a,s,f,g);
    return g;
  }

  public static float[][] taper(int m, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = copy(f);
    float[] t = new float[m];
    for (int i=0; i<m; ++i) {
      t[i] = (float)(0.54+0.46*cos(PI*(m-i)/m));
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0,j1=n1-1; i1<m; ++i1,--j1) {
        float ti = t[i1];
        g[i2][i1] *= ti;
        g[i2][j1] *= ti;
      }
    }
    for (int i2=0,j2=n2-1; i2<m; ++i2,--j2) {
      float ti = t[i2];
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] *= ti;
        g[j2][i1] *= ti;
      }
    }
    return g;
  }

  /**
   * Returns left (fm) and right (fp) images aligned for specified slopes.
   * Uses the slopes p[i] to compute fm[i] = f[i-k] and fp[i] = f[i+k],
   * such that events in fm[i] and fp[i] are aligned.
   * and fpp[i] = f[i+k]*f[i+k]. These products are used elsewhere to 
   * compute numerators and denominators of local correlation coefficients.
   * @param k half the cross-correlation lag.
   * @param p array of slopes.
   * @param f array of image samples.
   * @return array {fm,fp} of aligned images.
   */
  public float[][][] align(int k, float[][] p, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] fm = new float[n2][n1];
    float[][] fp = new float[n2][n1];
    float[] xm = new float[n1];
    float[] xp = new float[n1];
    _si.setUniformSampling(n1,1.0,0.0);
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-k,0);
      int i2p = min(i2+k,n2-1);
      float[] p2m = p[i2m];
      float[] p2p = p[i2p];
      float[] f2m = f[i2m];
      float[] f2p = f[i2p];
      float[] fm2 = fm[i2];
      float[] fp2 = fp[i2];
      for (int i1=0; i1<n1; ++i1) {
        xm[i1] = i1-k*p2m[i1];
        xp[i1] = i1+k*p2p[i1];
      }
      _si.setUniformSamples(f2m);
      _si.interpolate(n1,xm,fm2);
      _si.setUniformSamples(f2p);
      _si.interpolate(n1,xp,fp2);
    }
    return new float[][][]{fm,fp};
  }

  /**
   * Returns local correlation coefficients for specified angle.
   * @param theta the angle theta, in degrees.
   * @param fmp array {fm,fp} of two slope-aligned images.
   * @return array of local correlation coefficients.
   */
  public float[][] xcor(float theta, float[][][] fmp) {
    int n1 = fmp[0][0].length;
    int n2 = fmp[0].length;
    float[][] fm = fmp[0];
    float[][] fp = fmp[1];
    float[][] cmp = new float[n2][n1];
    float[][] cmm = new float[n2][n1];
    float[][] cpp = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fmi = fm[i2][i1];
        float fpi = fp[i2][i1];
        cmp[i2][i1] = fmi*fpi;
        cmm[i2][i1] = fmi*fmi;
        cpp[i2][i1] = fpi*fpi;
      }
    }
    float[][][] cs = {cmp,cmm,cpp};
    FaultLineSmoother fls = new FaultLineSmoother(_sigma,2.0,cs);
    fls.apply(theta,cs);
    float[][] c = cmp;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float cmpi = cmp[i2][i1];
        float cmmi = cmm[i2][i1];
        float cppi = cpp[i2][i1];
        if (cmpi<0.0f) {
          c[i2][i1] = 0.0f;
        } else {
          float cnum = cmpi*cmpi;
          float cden = cmmi*cppi;
          c[i2][i1] = (cnum<cden)?cnum/cden:1.0f;
        }
      }
    }
    clip(0.0f,1.0f,c,c);
    return c;
  }

  public float[][] findShifts(float[][][] ct, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] c = ct[0];
    //float[][] t = ct[1]; // currently not used
    float[][] u = new float[n2][n1];
    boolean[][] b = new boolean[n2][n1];
    int ka = FaultSamples.KA;
    int kb = FaultSamples.KB;
    int uMinShift = -(int)_shiftMax;
    int uMaxShift =  (int)_shiftMax;
    int pMinShift = -(int)_slopeMax*kb;
    int pMaxShift =  (int)_slopeMax*kb;
    LocalShiftFinderX lsfu = new LocalShiftFinderX(_sigma);
    LocalShiftFinderX lsfp = new LocalShiftFinderX(8.0);
    lsfu.setSmoothShifts(true);
    lsfp.setSmoothShifts(false);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (!b[i2][i1] && c[i2][i1]>0.0f) {
          FaultSamples fs = findFaultSamples(i1,i2,c,f,b);
          if (fs!=null) {
            int[] k1 = fs.k1;
            int[] k2 = fs.k2;
            float[] fmb = fs.fmb;
            float[] fma = fs.fma;
            float[] fpa = fs.fpa;
            float[] fpb = fs.fpb;
            int n = k1.length;
            float[] pm = new float[n];
            float[] pp = new float[n];
            float[] uk = new float[n];
            lsfp.find1(pMinShift,pMaxShift,fmb,fma,pm);
            lsfp.find1(pMinShift,pMaxShift,fpa,fpb,pp);
            lsfu.find1(uMinShift,uMaxShift,fma,fpa,uk);
            if (n>1500) {
              System.out.println("i1="+i1+" i2="+i2+" n="+n);
              SimplePlot sp = new SimplePlot();
              sp.setTitle("i1="+i1+" i2="+i2);
              //sp.addPoints(fm2).setLineColor(Color.GREEN);
              sp.addPoints(fma).setLineColor(Color.RED);
              sp.addPoints(fpa).setLineColor(Color.BLUE);
              sp = new SimplePlot();
              sp.setTitle("i1="+i1+" i2="+i2);
              sp.addPoints(uk).setLineColor(Color.RED);
              //sp.addPoints(fp2).setLineColor(Color.MAGENTA);
            }
            float scale = (float)ka/(float)(kb-ka);
            for (int k=0; k<n; ++k)
              u[k2[k]][k1[k]] = uk[k]-scale*(pm[k]+pp[k]);
          }
        }
      }
    }
    return u;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

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

  private float _sigma;
  private float _slopeMax;
  private float _thetaMax;
  private float _shiftMax;
  private float _faultLengthMin;
  private RecursiveGaussianFilter _rgf1;
  private SincInterpolator _si;
  private Sampling _st;
}
