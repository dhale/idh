/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import java.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * Computes 2D fault likelihoods by scanning over fault dip angles.
 * Fault likelihoods are in the range [0,1], where 0 and 1 denote 
 * lowest and highest likelihoods, respectively.
 * <p>
 * A fault scanner is constructed from precomputed semblance numerators 
 * and denominators. During a scan over fault orientations, these images 
 * are smoothed within fault lines having different dip angles, before 
 * computing the semblance ratios. Fault likelihoods are computed from 
 * semblances. For each sample location, the fault likelihood equals the
 * maximum likelihood computed at that location during a scan over all
 * dip angles.
 * <p>
 * Scan results include the dip angles for which the maximum likelihoods 
 * occured. Fault dip angles in [-90,90] degrees are measured from the 
 * vertical 1st image axis, and are positive when the horizontal projection 
 * of the dip vector that points downward has a positive component in the 
 * direction of the 2nd image axis. Note that fault dip is here defined 
 * such that vertical faults have zero dip.
 * <p>
 * EXPERIMENTAL: this class provides two methods for smoothing
 * semblance numerators and denominators along fault lines.
 * One method first shears those factors before applying fast 
 * recursive axis-aligned smoothing filters. The other method 
 * uses rotated Gaussian filters implemented with FFTs.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.10.31
 */
public class FaultScanner2 {

  /**
   * Method used to smooth semblance numerators and denominators.
   */
  public enum Smoother {
    SHEAR,
    FFT
  }

  /**
   * Constructs a fault scanner with specified parameters.
   * @param sigmaTheta half-width for smoothing along dip of fault lines.
   * @param snd array {snum,sden} of semblance numerators/denominators.
   */
  public FaultScanner2(double sigmaTheta, float[][][] snd) {
    this(sigmaTheta,snd,Smoother.SHEAR);
  }

  /**
   * Constructs a fault scanner with specified parameters.
   * @param sigmaTheta half-width for smoothing up-down dip of fault planes.
   * @param snd array {snum,sden} of semblance numerators/denominators.
   * @param s smoothing method.
   */
  public FaultScanner2(double sigmaTheta, float[][][] snd, Smoother s) {
    _sigmaTheta = sigmaTheta;
    if (s==Smoother.SHEAR) {
      _snd = snd;
    } else {
      _fls = new FaultLineSmoother(_sigmaTheta,0.0,snd);
    }
    _n1 = snd[0][0].length;
    _n2 = snd[0].length;
    _rgf1 = new RecursiveGaussianFilter(1.0);
    _rgf4 = new RecursiveGaussianFilter(4.0);
  }

  /**
   * Returns a sampling for fault dips. This sampling ensures that 
   * fault dips are adequately sampled in scans for fault orientations.
   * @param thetaMin the minimum dip, in degrees.
   * @param thetaMax the maximum dip, in degrees.
   * @return the dip sampling.
   */
  public Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigmaTheta,thetaMin,thetaMax);
  }

  /**
   * Computes and returns fault likelihoods for only one dip.
   * @param theta fault dip, in degrees.
   * @return array of fault likelihoods.
   */
  public float[][] likelihood(double theta) {
    float[][][] ft = scan(theta,theta);
    return ft[0];
  }

  /**
   * Scans within specified bounds for fault dips.
   * @param thetaMin minimum fault dip, in degrees.
   * @param thetaMax maximum fault dip, in degrees.
   * @return array {f,t} of fault likelihoods, strikes, and dips.
   */
  public float[][][] scan(double thetaMin, double thetaMax) {
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scan(st);
  }

  /**
   * Scans with the specified sampling of fault dips.
   * @param thetaSampling sampling of fault dip angles, in degrees.
   * @return array {f,t} of fault likelihoods and dips.
   */
  public float[][][] scan(Sampling thetaSampling) {
    if (_snd!=null) {
      return scanS(thetaSampling);
    } else {
      return scanF(thetaSampling);
    }
  }

  public float[][][] thin(float[][][] ft) {
    int n1 = ft[0][0].length;
    int n2 = ft[0].length;
    float[][] f = ft[0];
    float[][] t = ft[1];
    f = copy(f);
    pow(f,4.00f,f);
    _rgf1.apply00(f,f);
    pow(f,0.25f,f);
    float[][] ff = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fm = f[i2-1][i1];
        float fp = f[i2+1][i1];
        float fi = f[i2  ][i1];
        float ti = t[i2  ][i1];
        if (fm<fi && fp<fi) {
          ff[i2][i1] = fi;
          tt[i2][i1] = ti;
        }
      }
    }
    return new float[][][]{ff,tt};
  }

  public static float[][] smooth(
    double sigma, float[][] p, float[][] f, float[][] g)
  {
    int n1 = g[0].length;
    int n2 = g.length;
    EigenTensors2 d = new EigenTensors2(n1,n2);
    d.setEigenvalues(0.001f,1.00f); // smooth mostly along structure
    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        //s[i2][i1] = 1.0f-pow(f[i2][i1],0.1f); // almost binary
        s[i2][i1] = (f[i2][i1]<0.1f)?1.0f:0.0f; // binary
        float pi = p[i2][i1];
        float u1i = 1.0f/sqrt(1.0f+pi*pi);
        float u2i = -pi*u1i;
        d.setEigenvectorU(i1,i2,u1i,u2i);
      }
    }
    float c = (float)(0.5*sigma*sigma);
    float[][] h = new float[n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,s,g,h);
    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _snd;
  private FaultLineSmoother _fls;
  private RecursiveGaussianFilter _rgf1;
  private RecursiveGaussianFilter _rgf4;
  private double _sigmaTheta;
  private int _n1,_n2;

  // This scan smooths semblance numerators and denominators along fault
  // lines using Gaussian filters implemented with FFTs.
  private float[][][] scanF(Sampling thetaSampling) {
    Sampling st = thetaSampling;
    FaultSemblance fs = new FaultSemblance();
    final int n1 = _n1;
    final int n2 = _n2;
    final float[][] f = new float[n2][n1];
    final float[][] t = new float[n2][n1];
    int nt = st.getCount();
    for (int it=0; it<nt; ++it) {
      final float theta = (float)st.getValue(it);
      float[][][] snd = _fls.apply(theta);
      final float[][] s = fs.semblanceFromNumDen(snd);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float si = s[i2][i1]; // semblance
          si = si*si; // semblance^2
          si = si*si; // semblance^4
          si = si*si; // semblance^8
          float fi = 1.0f-si;
          if (fi>f[i2][i1]) {
            f[i2][i1] = fi;
            t[i2][i1] = theta;
          }
        }
      }
    }
    return new float[][][]{f,t};
  }

  // This scan smooths semblance numerators and denominators along fault
  // lines by shearing those factors before applying fast recursive 
  // axis-aligned smoothing filters.
  private float[][][] scanS(Sampling thetaSampling) {
    Sampling st = thetaSampling;
    int nt = st.getCount();
    int n1 = _n1, n2 = _n2;
    float[][] sn = _snd[0];
    float[][] sd = _snd[1];
    float[][] f = new float[n2][n1];
    float[][] t = new float[n2][n1];
    SincInterpolator si = new SincInterpolator();
    RecursiveExponentialFilter ref = makeRef(_sigmaTheta);
    for (int it=0; it<nt; ++it) {
      float ti = (float)st.getValue(it);
      float theta = toRadians(ti);
      float shear = tan(theta);
      float sigma = (float)_sigmaTheta*cos(theta);
      float[][] sns = shear(si,shear,sn);
      float[][] sds = shear(si,shear,sd);
      ref.apply1(sns,sns);
      ref.apply1(sds,sds);
      float[][] s = unshear(si,shear,semblance(sns,sds));
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float ss = s[i2][i1]; // semblance
          ss = ss*ss; // semblance^2
          ss = ss*ss; // semblance^4
          ss = ss*ss; // semblance^8
          float fi = 1.0f-ss;
          if (fi<0.0f) fi = 0.0f;
          if (fi>1.0f) fi = 1.0f;
          if (fi>f[i2][i1]) {
            f[i2][i1] = fi;
            t[i2][i1] = ti;
          }
        }
      }
    }
    return new float[][][]{f,t};
  }

  // Sampling of angles depends on extent of smoothing.
  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  // Shear horizontally such that q(i1,i2) = p(i1,i2+s*i1).
  private static float[][] shear(
    SincInterpolator si, double s, float[][] p) 
  {
    int n1 = p[0].length;
    int n2p = p.length;
    int n2q = n2p+(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] q = new float[n2q][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    si.setUniform(n2p,1.0,0.0,pp);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2p; ++i2)
        pp[i2] = p[i2][i1];
      double f2q = (s<0.0f)?s*i1:s*i1-dqp;
      si.interpolate(n2q,1.0f,f2q,qq);
      for (int i2=0; i2<n2q; ++i2)
        q[i2][i1] = qq[i2];
    }
    return q;
  }

  // Unshear horizontally such that p(i1,i2) = q(i1,i2-s*i1).
  private static float[][] unshear(
    SincInterpolator si, double s, float[][] q) 
  {
    int n1 = q[0].length;
    int n2q = q.length;
    int n2p = n2q-(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] p = new float[n2p][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    si.setUniform(n2q,1.0,0.0,qq);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2q; ++i2)
        qq[i2] = q[i2][i1];
      double f2p = (s<0.0f)?-s*i1:-s*i1+dqp;
      si.interpolate(n2p,1.0f,f2p,pp);
      for (int i2=0; i2<n2p; ++i2)
        p[i2][i1] = pp[i2];
    }
    return p;
  }

  // Smoothing filter
  private static RecursiveExponentialFilter makeRef(double sigma) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    return ref;
  }

  // Computes semblance from specified numerators and denominators.
  private float[][] semblance(float[][] sn, float[][] sd) {
    int n1 = sn[0].length;
    int n2 = sn.length;
    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] s2 = s[i2];
      float[] sn2 = sn[i2];
      float[] sd2 = sd[i2];
      for (int i1=0; i1<n1; ++i1) {
        float sni = sn2[i1];
        float sdi = sd2[i1];
        if (sdi<=0.0f || sni<0.0f) {
          s2[i1] = 0.0f;
        } else if (sdi<sni) {
          s2[i1] = 1.0f;
        } else {
          s2[i1] = sni/sdi;
        }
      }
    }
    return s;
  }

  // Returns a mapping j2[i2][i1] from indices (i1,i2) in sheared
  // coordinates to indices j2 of (i1,j2) in unsheared uncoordinates.
  private static int[][] j2Map(float s, int n1, int n2) {
    int ksn1 = (int)(s*n1);
    int[][] j2 = new int[n2][n1];
    for (int i1=0; i1<n1; ++i1) {
      int k2 = round((s<0.0f)?s*i1:s*i1-ksn1);
      for (int i2=0; i2<n2; ++i2) {
        int j2i = i2+k2;
        j2[i2][i1] = (0<=j2i && j2i<n2)?j2i:-1;
      }
    }
    return j2;
  }

  // This scan computes shifts along fault lines.
  public float[][] scanForShifts(
    int shiftMin, int shiftMax, Sampling thetaSampling,
    float fmin, float[][][] ft, float[][] p, float[][] g) 
  {
    DynamicWarping dw = new DynamicWarping(shiftMin,shiftMax);
    dw.setStretchMax(0.125);
    int nt = thetaSampling.getCount();
    int n1 = g[0].length;
    int n2 = g.length;
    float[][] f = ft[0];
    float[][] t = ft[1];
    float[][] s = new float[n2][n1];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    int nwarp = 0;
    for (int it=0; it<nt; ++it) {
      float ti = (float)thetaSampling.getValue(it);
      //if (ti!=4.5f) continue;
      float theta = toRadians(ti);
      float shear = tan(theta);
      float[][] gs = shear(si,shear,g);
      int[][] j2 = j2Map(shear,n1,n2);
      for (int i2=2; i2<n2-2; ++i2) {
        boolean mustWarp = false;
        for (int i1=0; i1<n1 && !mustWarp; ++i1) {
          int k2 = j2[i2][i1];
          mustWarp = k2>=0 && ti==t[k2][i1] && fmin<f[k2][i1];
        }
        if (mustWarp) {
          ++nwarp;
          //int i2m = (ti>=0.0)?i2-2:i2+2;
          //int i2p = (ti>=0.0)?i2+2:i2-2;
          int i2m = i2-2; 
          int i2p = i2+2; 
          float sp = 0.5f*(i2p-i2m);
          float[] pm = p[i2m];
          float[] pp = p[i2p];
          float[][] e = dw.computeErrors(gs[i2m],gs[i2p]);
          float[][] d = dw.accumulateForward(e);
          float[] u = dw.findShiftsReverse(d,e);
          //_rgf4.apply0(u,u);
          for (int i1=0; i1<n1; ++i1) {
            int k2 = j2[i2][i1];
            if (k2>=0 && ti==t[k2][i1] && fmin<f[k2][i1]) {
              s[k2][i1] = u[i1]-sp*(pm[i1]+pp[i1]);
              if (k2==149 && i1==30) {
                System.out.println(
                  "f="+f[k2][i1]+" t="+t[k2][i1]+" s="+s[k2][i1]);
                SimplePlot plot = new SimplePlot();
                plot.addPoints(gs[i2m]).setLineColor(Color.BLACK);
                plot.addPoints(gs[i2p]).setLineColor(Color.RED);
                SimplePlot plot2 = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
                float[] x1 = rampfloat(0.0f,1.0f,n1);
                Sampling s1 = new Sampling(n1,1.0f,0.0f);
                Sampling sl =
                  new Sampling(1+(shiftMax-shiftMin),1.0f,shiftMin);
                plot2.addPixels(sl,s1,pow(e,0.25f));
                plot2.addPoints(u,x1).setLineColor(Color.YELLOW);
              }
            }
          }
        }
      }
      /*
      SimplePlot.asPixels(g);
      SimplePlot.asPixels(f);
      SimplePlot.asPixels(t);
      SimplePlot.asPixels(s);
      */
    }
    System.out.println("warp fraction = "+(float)nwarp/n2/nt);
    return s;
  }
}
