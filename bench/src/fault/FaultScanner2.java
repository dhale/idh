/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import java.util.ArrayList;
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
      _fls = new FaultLineSmoother(_sigmaTheta,1.0,snd);
    }
    _n1 = snd[0][0].length;
    _n2 = snd[0].length;
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
    float[][] fs = new float[n2][n1];
    float[][] ff = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    new RecursiveGaussianFilter(1.0).apply00(f,fs);
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fm = fs[i2-1][i1];
        float fp = fs[i2+1][i1];
        float fi = fs[i2  ][i1];
        if (fm<fi && fp<fi) {
          ff[i2][i1] = f[i2][i1];
          tt[i2][i1] = t[i2][i1];
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
    SincInterp si = new SincInterp();
    si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
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
    SincInterp si, double s, float[][] p) 
  {
    int n1 = p[0].length;
    int n2p = p.length;
    int n2q = n2p+(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] q = new float[n2q][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2p; ++i2)
        pp[i2] = p[i2][i1];
      double f2q = (s<0.0f)?s*i1:s*i1-dqp;
      si.interpolate(n2p,1.0,0.0,pp,n2q,1.0f,f2q,qq);
      for (int i2=0; i2<n2q; ++i2)
        q[i2][i1] = qq[i2];
    }
    return q;
  }

  // Unshear horizontally such that p(i1,i2) = q(i1,i2-s*i1).
  private static float[][] unshear(
    SincInterp si, double s, float[][] q) 
  {
    int n1 = q[0].length;
    int n2q = q.length;
    int n2p = n2q-(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] p = new float[n2p][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2q; ++i2)
        qq[i2] = q[i2][i1];
      double f2p = (s<0.0f)?-s*i1:-s*i1+dqp;
      si.interpolate(n2q,1.0,0.0,qq,n2p,1.0f,f2p,pp);
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

  ///////////////////////////////////////////////////////////////////////////
  // fault analysis

  public Faults findFaults(float[][][] ft, int minLength) {
    return new Faults(ft,minLength);
  }

  public static class Faults {
    Faults(float[][][] ft, int minLength) {
      _n1 = ft[0][0].length;
      _n2 = ft[0].length;
      FaultNode[][] fns = findFaultNodes(ft);
      _fls = findFaultLines(fns,minLength);
    }

    public void findShifts(float[][] g, float[][] p, int smin, int smax) {
      for (FaultLine fl:_fls)
        findShiftsFl(fl,g,p,smin,smax);
    }

    public void clean() {
      ArrayList<FaultLine> fls = new ArrayList<FaultLine>();
      for (FaultLine fl:_fls)
        if (shiftsValid(fl))
          fls.add(fl);
      _fls = fls;
    }

    public int getCount() {
      return _fls.size();
    }

    public float[][] getShifts() {
      float[][] s = new float[_n2][_n1];
      for (FaultLine fl:_fls) {
        FaultNode[] fns = fl.getNodes();
        int nn = fns.length;
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns[in];
          int i1 = fn.i1;
          int i2m = fn.i2m;
          int i2p = fn.i2p;
          s[i2m][i1] = fn.smp;
          s[i2p][i1] = fn.spm;
        }
      }
      return s;
    }

    public float[][] getLikelihoods() {
      float[][] f = new float[_n2][_n1];
      for (FaultLine fl:_fls) {
        FaultNode[] fns = fl.getNodes();
        int nn = fns.length;
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns[in];
          int i1 = fn.i1;
          int i2 = fn.i2;
          f[i2][i1] = fn.fl;
        }
      }
      return f;
    }

    public float[][] getDips() {
      float[][] t = new float[_n2][_n1];
      for (FaultLine fl:_fls) {
        FaultNode[] fns = fl.getNodes();
        int nn = fns.length;
        for (int in=0; in<nn; ++in) {
          FaultNode fn = fns[in];
          int i1 = fn.i1;
          int i2 = fn.i2;
          t[i2][i1] = fn.ft;
        }
      }
      return t;
    }

    private int _n1,_n2;
    private ArrayList<FaultLine> _fls;
  }

  // One node in a fault corresponds to one image sample.
  private static class FaultNode {
    FaultNode(int i1, int i2, int i2m, int i2p, float x2, float fl, float ft) {
      this.i1 = i1;
      this.i2 = i2;
      this.i2m = i2m;
      this.i2p = i2p;
      this.x2 = x2;
      this.fl = fl;
      this.ft = ft;
    }
    int i1,i2; // sample indices (i1,i2) for this node
    int i2m,i2p; // indices i2 for minus and plus sides of fault
    float x2; // horizontal fault location in interval [i2m,i2p]
    float fl,ft; // fault likelihood and theta
    float smp,spm; // fault shifts both mp and pm sides of fault
    FaultNode above; // node above
    FaultNode below; // node below
  }

  // A fault line is a linked list of fault nodes.
  private static class FaultLine {
    FaultLine(FaultNode top) {
      this.top = top;
      this.n = 1;
      for (FaultNode fn=top.below; fn!=null; fn=fn.below)
        ++n;
    }
    int n; // number of nodes
    FaultNode top; // top node
    FaultNode[] getNodes() {
      FaultNode[] fns = new FaultNode[n];
      int in = 0;
      for (FaultNode fn=top; fn!=null; fn=fn.below,++in)
        fns[in] = fn;
      return fns;
    }
  }

  private static FaultNode[][] findFaultNodes(float[][][] ft) {
    int n1 = ft[0][0].length;
    int n2 = ft[0].length;
    float[][] f = ft[0];
    float[][] t = ft[1];
    float[][] fs = new float[n2][n1];
    new RecursiveGaussianFilter(1.0).apply00(f,fs);
    FaultNode[][] fns = new FaultNode[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fm = fs[i2-1][i1];
        float fp = fs[i2+1][i1];
        float fi = fs[i2][i1];
        if (fm<fi && fp<fi) {
          float c0 = fi;
          float c1 = 0.5f*(fp-fm);
          float c2 = 0.5f*(fp+fm)-fi;
          float ui = (c2<0.0f)?-0.5f*c1/c2:0.0f;
          float x2 = i2+ui;
          int i2m = i2;
          int i2p = i2;
          if (ui>=0.0f) {
            ++i2p;
          } else {
            --i2m;
          }
          fns[i2][i1] = new FaultNode(i1,i2,i2m,i2p,x2,f[i2][i1],t[i2][i1]);
        }
      }
    }
    return fns;
  }

  private static ArrayList<FaultLine> findFaultLines(
    FaultNode[][] fns, int minNodesInLine) 
  {
    int n1 = fns[0].length;
    int n2 = fns.length;

    // List of fault lines.
    ArrayList<FaultLine> fls = new ArrayList<FaultLine>();

    // For all fault nodes in the array of nodes, ...
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        FaultNode fni = fns[i2][i1];

        // If a node exists and has not yet been linked to others, ...
        if (fni!=null && fni.above==null && fni.below==null) {

          // Count the number of nodes.
          int nn = 1;

          // Search below for nabors to link to this node.
          for (FaultNode fnt=fni; fnt!=null; fnt=fnt.below) {
            FaultNode fnb = naborBelow(fnt,fns);
            if (fnb!=null && fnt==naborAbove(fnb,fns)) {
              fnb.above = fnt;
              fnt.below = fnb;
              ++nn;
            }
          }

          // Search above for nabors to link to this node.
          // Also remember the top node, which may become the
          // head node in the list of nodes for a fault line.
          FaultNode fnr = fni;
          for (FaultNode fnt=fni; fnt!=null; fnt=fnt.above) {
            FaultNode fna = naborAbove(fnt,fns);
            if (fna!=null && fnt==naborBelow(fna,fns)) {
              fna.below = fnt;
              fnt.above = fna;
              ++nn;
            }
            fnr = fnt;
          }

          // If sufficient number of nodes, construct new fault line.
          if (nn>=minNodesInLine)
            fls.add(new FaultLine(fnr));
        }
      }
    }
    return fls;
  }

  private static float[][] getImageSamples(FaultLine fl, int k, float[][] g) {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    int n1 = g[0].length;
    int n2 = g.length;
    int n2m = n2-1;
    float[] gm = new float[nn];
    float[] gp = new float[nn];
    for (int in=0; in<nn; ++in) {
      FaultNode fn = fns[in];
      int i1 = fn.i1;
      int i2 = fn.i2;
      int i2m = max(  0,i2-k);
      int i2p = min(n2m,i2+k);
      gm[in] = g[i2m][i1];
      gp[in] = g[i2p][i1];
    }
    return new float[][]{gm,gp};
  }

  private static boolean shiftsValid(FaultLine fl) {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    int nv = 0;
    for (int in=0; in<nn; ++in) {
      FaultNode fn = fns[in];
      if (fn.smp*fn.spm<=0.0f)
        ++nv;
    }
    return 10*nv>8*nn; // 80% of shifts must have opposite signs
  }

  private static void cleanShifts(FaultLine fl) {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    for (int in=0; in<nn; ++in) {
      FaultNode fn = fns[in];
      if (fn.smp*fn.spm>=0.0f) {
        fn.smp = 0.0f;
        fn.spm = 0.0f;
      }
    }
  }

  private static void findShiftsFl(
    FaultLine fl, float[][] g, float[][] p, int smin, int smax) 
  {
    findShiftsFl( 1,fl,g,p,smin,smax);
    findShiftsFl(-1,fl,g,p,smin,smax);
  }

  private static void findShiftsFl(
    int sign, FaultLine fl, float[][] g, float[][] p, int smin, int smax) 
  {
    FaultNode[] fns = fl.getNodes();
    int nn = fns.length;
    int ssmin = (nn+smin>10)?smin:10-nn;
    int ssmax = (nn-smax>10)?smax:nn-10;
    DynamicWarping dw = new DynamicWarping(ssmin,ssmax);
    dw.setStretchMax(0.5);
    float[][] gmp = getImageSamples(fl,2,g);
    float[][] pmp = getImageSamples(fl,2,p);
    float[] g1,g2,p1,p2;
    float ps;
    if (sign>0) {
      g1 = gmp[0]; g2 = gmp[1];
      p1 = pmp[0]; p2 = pmp[1];
      ps = 2.0f;
    } else {
      g1 = gmp[1]; g2 = gmp[0];
      p1 = pmp[1]; p2 = pmp[0];
      ps = -2.0f;
    }
    float[][] se = dw.computeErrors(g1,g2);
    float[][] sd = dw.accumulateForward(se);
    float[] s = dw.findShiftsReverse(sd,se);
    for (int in=0; in<nn; ++in)
      s[in] -= ps*(p1[in]+p2[in]);
    if (sign>0) {
      for (int in=0; in<nn; ++in)
        fns[in].smp = s[in];
    } else {
      for (int in=0; in<nn; ++in)
        fns[in].spm = s[in];
    }
  }

  // Minimum acceptable dot product of two aligned fault dip vectors.
  // Nodes can be neighbors only if their dip vectors are so aligned.
  private static final float AMIN = 0.9f;

  // Returns dot product of dip vectors for two specified nodes.
  private static float dotDipVectors(FaultNode fn1, FaultNode fn2) {
    float t1 = toRadians(fn1.ft);
    float t2 = toRadians(fn2.ft);
    float c1 = cos(t1);
    float c2 = cos(t2);
    float s1 = sin(t1);
    float s2 = sin(t2);
    return c1*c2+s1*s2;
  }

  // Determines which (if any) of nodes f1,f2,f3 is aligned with fn.
  private static FaultNode alignDip(
    FaultNode fn, FaultNode f1, FaultNode f2, FaultNode f3) 
  {
    float a1 = (f1!=null)?dotDipVectors(fn,f1):0.0f;
    float a2 = (f2!=null)?dotDipVectors(fn,f2):0.0f;
    float a3 = (f3!=null)?dotDipVectors(fn,f3):0.0f;
    float am = max(a1,a2,a3);
    FaultNode fm = null;
    if (am>=AMIN) 
      fm = (a1==am)?f1:(a2==am)?f2:f3;
    return fm;
  }

  // Returns the best nabor (if any) above the specified node.
  private static FaultNode naborAbove(FaultNode fn, FaultNode[][] fns) {
    int n1 = fns[0].length;
    int n2 = fns.length;
    int i1 = fn.i1;
    int i2 = fn.i2;
    FaultNode fm = null;
    if (i1>0) {
      FaultNode f1 = fns[i2][i1-1];
      FaultNode f2 = (i2>0)?fns[i2-1][i1-1]:null;
      FaultNode f3 = (i2<n2-1)?fns[i2+1][i1-1]:null;
      fm = alignDip(fn,f1,f2,f3);
    }
    return fm;
  }

  // Returns the best nabor (if any) below the specified node.
  private static FaultNode naborBelow(FaultNode fn, FaultNode[][] fns) {
    int n1 = fns[0].length;
    int n2 = fns.length;
    int i1 = fn.i1;
    int i2 = fn.i2;
    FaultNode fm = null;
    if (i1<n1-1) {
      FaultNode f1 = fns[i2][i1+1];
      FaultNode f2 = (i2>0)?fns[i2-1][i1+1]:null;
      FaultNode f3 = (i2<n2-1)?fns[i2+1][i1+1]:null;
      fm = alignDip(fn,f1,f2,f3);
    }
    return fm;
  }
}
