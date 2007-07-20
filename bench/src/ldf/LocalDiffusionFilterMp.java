/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.io.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic diffusion filter via minimum-phase factors.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.12
 */
public class LocalDiffusionFilterMp extends LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   */
  public LocalDiffusionFilterMp(double sigma) {
    super(sigma);
    _sigma = (float)sigma;
    _dlf = new DirectionalLaplacianFilter(1.0f);
  }

  /**
   * Constructs a local diffusion filter with pre-computed coefficients.
   * Coefficients are stored in a file with specified name.
   * @param sigma the nominal half-width for this filter.
   * @param fileName name of file containing pre-computed filters.
   */
  public LocalDiffusionFilterMp(double sigma, String fileName) {
    super(sigma);
    _sigma = (float)sigma;
    _dlf = new DirectionalLaplacianFilter(1.0f);
    _file = new File(fileName);
    Check.argument(_file.exists(),"file "+fileName+" exists");
    Check.argument(_file.canRead(),"file "+fileName+" is readable");
    Check.argument(_file.isFile(),"file "+fileName+" is a normal file");
  }

  public static void precomputeFilters(String fileName) {
  }

  /*
   * 2-D NSIGMA,NVECTOR,NLAG
   * 2-D sigma[NSIGMA]
   * 2-D vector[NSIGMA][NVECTOR][2]
   * 2-D inline[NSIGMA][NVECTOR][NLAG]
   * 2-D normal[NSIGMA][NVECTOR][NLAG]
   *
   * 3-D NSIGMA,NVECTOR,NLAG
   * 3-D sigma[NSIGMA]
   * 3-D vector[NSIGMA][NVECTOR][3]
   * 3-D inline[NSIGMA][NVECTOR][NLAG]
   * 3-D normal[NSIGMA][NVECTOR][NLAG]
   */

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void solveInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] t = new float[n2][n1];
    _dlf.applyInline(null,v1,x,t);
    if (_fif2==null)
      _fif2 = new FactoredFilter2(FactoredFilter2.Type.INLINE);
    _fif2.applyInverse(_sigma,ds,v1,t,y);
    Array.sub(x,y,y);
  }

  protected void solveInline(
    float[][][] ds, short[][][] iw, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] t = new float[n3][n2][n1];
    _dlf.applyInline(null,iw,x,t);
    if (_fif3==null)
      _fif3 = new FactoredFilter3(FactoredFilter3.Type.INLINE);
    _fif3.applyInverse(_sigma,ds,iw,t,y);
    Array.sub(x,y,y);
  }

  protected void solveNormal(
    float[][][] ds, short[][][] iu, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] t = new float[n3][n2][n1];
    _dlf.applyInline(null,iu,x,t);
    if (_fif3==null)
      _fif3 = new FactoredFilter3(FactoredFilter3.Type.NORMAL);
    _fif3.applyInverse(_sigma,ds,iu,t,y);
    Array.sub(x,y,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma;
  private File _file;
  private DirectionalLaplacianFilter _dlf;
  private static FactoredFilter2 _fif2; // 2-D inline filter
  private static FactoredFilter2 _fnf2; // 2-D normal filter
  private static FactoredFilter3 _fif3; // 3-D inline filter
  private static FactoredFilter3 _fnf3; // 3-D normal filter
  private static UnitSphereSampling _uss16 = new UnitSphereSampling(16);

  // A local inline filter approximated with minimum-phase factors.
  // Factors are tabulated as a function of sigma and theta.
  private static class FactoredFilter2 {

    static enum Type {
      INLINE,
      NORMAL
    };

    FactoredFilter2(Type type) {
      trace("FactoredFilter2: constructing filters ...");

      // A causal filter to compute tabulated factors.
      CausalFilter cf = new CausalFilter(LAG1,LAG2);

      // Arrays used to get 3x3 auto-correlations to be factored.
      float[][] t = new float[3][3];
      float[][] r = new float[3][3];
      t[1][1] = 1.0f;

      // A filter to compute the auto-correlations to be factored. 
      // With sigma = sqrt(2.0), we compensate for the scaling by 
      // sigma*sigma/2 performed by the directional Laplacian filter.
      DirectionalLaplacianFilter dlf = 
        new DirectionalLaplacianFilter(sqrt(2.0f));

      // Which table of coefficients is being computed?
      float[][][] atable = (type==Type.INLINE)?ATABLE_INLINE:ATABLE_NORMAL;

      // Tabulated factors for all angles theta and half-widths sigma.
      // Note that for normal filters we specify u1 = v2 and u2 = v1.
      // Because vectors u and v are orthogonal, we might more properly
      // use u2 = -v1 = -sin(theta). However, during table lookup we will 
      // use u2 = v1 = sin(theta) for both inline and normal filters, so
      // we need to do the same here.
      for (int itheta=0; itheta<NTHETA; ++itheta) {
        float theta = FTHETA+itheta*DTHETA;
        for (int isigma=0; isigma<NSIGMA; ++isigma) {
          float sigma = FSIGMA+isigma*DSIGMA;
          float scale = 2.0f/(sigma*sigma);
          Array.mul(scale,t,r);
          float v1 = sin(theta);
          float v2 = cos(theta);
          if (type==Type.INLINE) {
            dlf.applyInline(1.0f,v1,v2,t,r);
          } else {
            dlf.applyNormal(1.0f,v2,v1,t,r); // use u2 = v1 (not -v1)!
          }
          cf.factorWilsonBurg(100,0.000001f,r);
          atable[itheta][isigma] = cf.getA();
        }
      }

      _type = type;
      trace("...  done.");
    }

    void applyForward(
      float sigma, float[][] ds, float[][] v1, float[][] x, float[][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?ATABLE_INLINE:ATABLE_NORMAL;
      A2 a2 = new A2(atable,sigma,ds,v1);
      LCF.apply(a2,x,y);
      LCF.applyTranspose(a2,y,y);
    }

    void applyInverse(
      float sigma, float[][] ds, float[][] v1, float[][] x, float[][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?ATABLE_INLINE:ATABLE_NORMAL;
      A2 a2 = new A2(atable,sigma,ds,v1);
      LCF.applyInverseTranspose(a2,x,y);
      LCF.applyInverse(a2,y,y);
    }

    private Type _type;

    // Provides filter coefficients for each sample index via
    // bilinear interpolation of pre-computed filter coefficients.
    // This class looks as if written for only inline filters, but it
    // can also handle normal filters by passing an appropriate table 
    // of coefficients to the constructor.
    private static class A2 implements LocalCausalFilter.A2 {
      A2(float[][][] atable, float sigma, float[][] ds, float[][] v1) 
      {
        _at = atable;
        _sigma = sigma;
        _ds = ds;
        _v1 = v1;
      }
      public void get(int i1, int i2, float[] a) {
        float theta = asin(_v1[i2][i1]);
        float sigma = _sigma;
        if (_ds!=null) sigma *= _ds[i2][i1];
        if (sigma<SIGMA_MIN) sigma = SIGMA_MIN;
        if (sigma>SIGMA_MAX) sigma = SIGMA_MAX;
        float s = (sigma-FSIGMA)*SSIGMA;
        float t = (theta-FTHETA)*STHETA;
        int is = (int)s;
        int it = (int)t;
        float s1 = s-(float)is;
        float t1 = t-(float)it;
        float s0 = 1.0f-s1;
        float t0 = 1.0f-t1;
        float[] a00 = _at[it  ][is  ];
        float[] a01 = _at[it  ][is+1];
        float[] a10 = _at[it+1][is  ];
        float[] a11 = _at[it+1][is+1];
        int n = a00.length;
        for (int j=0; j<n; ++j)
          a[j] = t0*(s0*a00[j]+s1*a01[j])+t1*(s0*a10[j]+s1*a11[j]);
      }
      private float[][][] _at;
      private float _sigma;
      private float[][] _ds;
      private float[][] _v1;
    }

    // Sampling of sigma for tabulated filter coefficients.
    private static final float SIGMA_MIN =  0.1f;
    private static final float SIGMA_MAX = 20.0f;
    private static final int NSIGMA = 20;
    private static final float FSIGMA = SIGMA_MIN;
    private static final float DSIGMA = (SIGMA_MAX-SIGMA_MIN)/(float)(NSIGMA-1);
    private static final float SSIGMA = (1.0f-4.0f*FLT_EPSILON)/DSIGMA;

    // Sampling of theta for tabulated filter coefficients.
    private static final float THETA_MIN = -0.5f*FLT_PI;
    private static final float THETA_MAX =  0.5f*FLT_PI;
    private static final int NTHETA = 21;
    private static final float FTHETA = THETA_MIN;
    private static final float DTHETA = (THETA_MAX-THETA_MIN)/(float)(NTHETA-1);
    private static final float STHETA = (1.0f-4.0f*FLT_EPSILON)/DTHETA;

    // Tables of coefficients for both inline and normal filters.
    private static float[][][] ATABLE_INLINE = new float[NTHETA][NSIGMA][];
    private static float[][][] ATABLE_NORMAL = new float[NTHETA][NSIGMA][];

    // Lags for minimum-phase factors, with the following stencil:
    //   lag1 =    1  0 -1 -2 -3 -4
    //   --------------------------
    //   lag2 = 0: x  x
    //          1: x  x  x  x  x  x
    private static int MLAG = 4;
    private static int NLAG = 4+MLAG;
    private static int[] LAG1 = new int[NLAG];
    private static int[] LAG2 = new int[NLAG];
    private static LocalCausalFilter LCF;
    static {
      for (int ilag=0; ilag<NLAG; ++ilag) {
        LAG1[ilag] = (ilag<=1)?ilag:ilag-2-MLAG;
        LAG2[ilag] = (ilag<=1)?0:1;
      }
      LCF = new LocalCausalFilter(LAG1,LAG2);
    }
  }

  // A filter approximated with minimum-phase factors. Factors are 
  // tabulated as a function of sigma and a unit-sphere sampling index.
  private static class FactoredFilter3 {

    static enum Type {
      INLINE,
      NORMAL
    };

    FactoredFilter3(Type type) {
      trace("FactoredFilter3: constructing filters ...");

      // A causal filter to compute tabulated factors.
      CausalFilter cf = new CausalFilter(LAG1,LAG2,LAG3);

      // Arrays used to compute 3x3x3 auto-correlations to be factored.
      float[][][] t = new float[3][3][3];
      float[][][] r = new float[3][3][3];
      t[1][1][1] = 1.0f;

      // A filter to compute the auto-correlations to be factored. 
      // With sigma = sqrt(2.0), we compensate for the scaling by 
      // sigma*sigma/2 performed by the directional Laplacian filter.
      DirectionalLaplacianFilter dlf = 
        new DirectionalLaplacianFilter(sqrt(2.0f));

      // Which table of coefficients is being computed?
      float[][][] atable = (type==Type.INLINE)?ATABLE_INLINE:ATABLE_NORMAL;

      // Tabulated factors for unit-vector indices and half-widths sigma.
      trace("NVEC="+NVEC+" NSIGMA="+NSIGMA);
      //for (int iw=1; iw<=NVEC; ++iw) {
      for (int iw=1; iw<=NVEC; iw+=10) {
        for (int isigma=NSIGMA-1; isigma<NSIGMA; ++isigma) {
          trace("iw="+iw+" isigma="+isigma);
          float sigma = FSIGMA+isigma*DSIGMA;
          float scale = 2.0f/(sigma*sigma);
          Array.mul(scale,t,r);
          float[] w = USS.getPoint(iw);
          float w1 = w[2];
          float w2 = w[1];
          float w3 = w[0];
          if (type==Type.INLINE) {
            dlf.applyInline(1.0f,w1,w2,w3,t,r);
          } else {
            dlf.applyNormal(1.0f,w1,w2,w3,t,r);
          }
          Array.dump(w);
          Array.dump(r);
          cf.factorWilsonBurg(100,0.000001f,r);
          atable[iw][isigma] = cf.getA();
          dumpA(cf.getA());
        }
      }

      _type = type;
      trace("...  done.");
    }
    static void dumpA(float[] a) {
      int nlag = a.length;
      int n = MLAG+1+MLAG;
      float[][] b = new float[2*n][n];
      for (int ilag=0; ilag<nlag; ++ilag) {
        int i1 = n-1-(LAG1[ilag]+MLAG);
        int i2 = n-1-(LAG2[ilag]+MLAG);
        int i3 = LAG3[ilag];
        b[i3*n+i2][i1] = a[ilag];
      }
      edu.mines.jtk.mosaic.SimplePlot.asPixels(b);
    }

    void applyForward(
      float sigma, float[][][] ds, short[][][] iw, 
      float[][][] x, float[][][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?ATABLE_INLINE:ATABLE_NORMAL;
      A3 a3 = new A3(atable,sigma,ds,iw);
      LCF.apply(a3,x,y);
      LCF.applyTranspose(a3,y,y);
    }

    void applyInverse(
      float sigma, float[][][] ds, short[][][] iw, 
      float[][][] x, float[][][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?ATABLE_INLINE:ATABLE_NORMAL;
      A3 a3 = new A3(atable,sigma,ds,iw);
      LCF.applyInverseTranspose(a3,x,y);
      LCF.applyInverse(a3,y,y);
    }

    private Type _type;

    // Provides filter coefficients for each sample index via
    // linear interpolation of pre-computed filter coefficients.
    private static class A3 implements LocalCausalFilter.A3 {
      A3(float[][][] atable, float sigma, float[][][] ds, short[][][] iw) {
        _sigma = sigma;
        _at = atable;
        _ds = ds;
        _iw = iw;
      }
      public void get(int i1, int i2, int i3, float[] a) {
        int iw = _iw[i3][i2][i1];
        int ia = IA[iw];
        int ib = IB[iw];
        int ic = IC[iw];
        float wa = WA[iw];
        float wb = WB[iw];
        float wc = WC[iw];
        float sigma = _sigma;
        if (_ds!=null) sigma *= _ds[i3][i2][i1];
        if (sigma<SIGMA_MIN) sigma = SIGMA_MIN;
        if (sigma>SIGMA_MAX) sigma = SIGMA_MAX;
        float s = (sigma-FSIGMA)*SSIGMA;
        int is = (int)s;
        float s1 = s-(float)is;
        float s0 = 1.0f-s1;
        float[][] aa = _at[ia];
        float[][] ab = _at[ib];
        float[][] ac = _at[ic];
        float[] aa0 = _at[ia][is  ];
        float[] aa1 = _at[ia][is+1];
        float[] ab0 = _at[ib][is  ];
        float[] ab1 = _at[ib][is+1];
        float[] ac0 = _at[ic][is  ];
        float[] ac1 = _at[ic][is+1];
        int n = aa0.length;
        for (int j=0; j<n; ++j)
          a[j] = s0*(wa*aa0[j]+wb*ab0[j]+wc*ac0[j]) +
                 s1*(wa*aa1[j]+wb*ab1[j]+wc*ac1[j]);
      }
      private float _sigma;
      private float[][][] _at;
      private float[][][] _ds;
      private short[][][] _iw;
    }

    // Sampling of sigma for tabulated filter coefficients.
    private static final float SIGMA_MIN =  0.1f;
    private static final float SIGMA_MAX = 20.0f;
    private static final int NSIGMA = 20;
    private static final float FSIGMA = SIGMA_MIN;
    private static final float DSIGMA = (SIGMA_MAX-SIGMA_MIN)/(float)(NSIGMA-1);
    private static final float SSIGMA = (1.0f-4.0f*FLT_EPSILON)/DSIGMA;

    // Sampling of the unit sphere for tabulated filter coefficients.
    // This sampling must be coarser (using fewer bits) than the 16-bit
    // sampling used to encode unit vectors.
    private static final int NBIT = 8;
    private static final UnitSphereSampling USS =
      new UnitSphereSampling(NBIT);
    private static final int NVEC = USS.getMaxIndex();

    // Precomputed weights for linear interpolation within triangles
    // of the coarse unit-sphere sampling. Because unit-vectors are
    // quantized to 16 bits, only these weights will be needed to
    // interpolate filter coefficients.
    private static int NW = _uss16.getMaxIndex(); // number of weights
    private static int[] IA = new int[1+NW]; // IA[0] unused
    private static int[] IB = new int[1+NW]; // IB[0] unused
    private static int[] IC = new int[1+NW]; // IC[0] unused
    private static float[] WA = new float[1+NW]; // WA[0] unused
    private static float[] WB = new float[1+NW]; // WB[0] unused
    private static float[] WC = new float[1+NW]; // WC[0] unused

    // Tables of coefficients for both inline and normal filters.
    private static float[][][] ATABLE_INLINE = new float[1+NVEC][NSIGMA][];
    private static float[][][] ATABLE_NORMAL = new float[1+NVEC][NSIGMA][];

    // Lags for minimum-phase factors, with the following stencil:
    //                lag1 =  4  3  2  1  0 -1 -2 -3 -4
    //   ----------------------------------------------
    //   lag3 = 0, lag2 =  0:          x  x
    //                     1:          x  x  x  x  x  x
    //   ----------------------------------------------
    //   lag3 = 1, lag2 = -4: x  x  x  x  x  x  x  x  x
    //                    -3: x  x  x  x  x  x  x  x  x
    //                    -2: x  x  x  x  x  x  x  x  x
    //                    -1: x  x  x  x  x  x  x  x  x
    //                     0:          x  x  x  x  x  x
    //                     1:          x  x  x  x  x  x
    private static int MLAG = 4;
    private static int NLAG = 2+3*(2+MLAG)+MLAG*(MLAG+1+MLAG); // = 56
    private static int[] LAG1 = new int[NLAG];
    private static int[] LAG2 = new int[NLAG];
    private static int[] LAG3 = new int[NLAG];
    private static LocalCausalFilter LCF;

    // Initialization of static fields.
    static {
      trace("nlag="+NLAG);
      for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
        int jlag2 = (ilag3==0)?0:-MLAG;
        int klag2 = 1;
        for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
          int jlag1 = (ilag3==0 && ilag2==0)?0:-MLAG;
          int klag1 = (ilag3==0 || ilag2>=0)?1: MLAG;
          for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
            LAG1[ilag] = ilag1;
            LAG2[ilag] = ilag2;
            LAG3[ilag] = ilag3;
            trace("ilag="+ilag+" lag1="+ilag1+" ilag2="+ilag2+" ilag3="+ilag3);
          }
        }
      }
      LCF = new LocalCausalFilter(LAG1,LAG2,LAG3);
      for (int iw=1; iw<NW; ++iw) {
        float[] wi = _uss16.getPoint(iw);
        int[] iabc = USS.getTriangle(wi);
        float[] wabc = USS.getWeights(wi,iabc);
        IA[iw] = iabc[0];
        IB[iw] = iabc[1];
        IC[iw] = iabc[2];
        WA[iw] = wabc[0];
        WB[iw] = wabc[1];
        WC[iw] = wabc[2];
      }
    }
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  public static void main(String[] args) {
    testFactoredFilter3();
  }

  private static void testFactoredFilter3() {
    //FactoredFilter2 ff2 = new FactoredFilter2(FactoredFilter2.Type.INLINE);
    FactoredFilter3 ff3 = new FactoredFilter3(FactoredFilter3.Type.NORMAL);
  }
} 
