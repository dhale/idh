/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.io.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local anisotropic diffusion filter via minimum-phase factors.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.23
 */
public class LocalDiffusionFilterMp extends LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   */
  public LocalDiffusionFilterMp(double sigma) {
    this(sigma,null);
  }

  /**
   * Constructs a local diffusion filter with pre-computed coefficients.
   * Coefficients are stored in a file with specified name.
   * @param sigma the nominal half-width for this filter.
   * @param ffile name of file containing pre-computed filters.
   *  If null, filters will be computed on-the-fly as necessary.
   */
  public LocalDiffusionFilterMp(double sigma, String ffile) {
    super(sigma);
    if (ffile!=null) {
      File file = new File(ffile);
      Check.argument(file.exists(),"file "+file+" exists");
      Check.argument(file.canRead(),"file "+file+" is readable");
      Check.argument(file.isFile(),"file "+file+" is a file");
    }
    _ffile = ffile;
    _sigma = (float)sigma;
    _dlf = new DirectionalLaplacianFilter(1.0f);
  }

  /**
   * Saves filter coefficients to a file with specified name.
   * @param ffile name of file to contain pre-computed filters.
   */
  public void save(String ffile) {
    try {
      ensureInlineFilter3();
      ensureNormalFilter3();
      ArrayFile af = new ArrayFile(ffile,"rw");
      af.writeInt(_fileFormat);
      _fif3.save(af);
      _fnf3.save(af);
      af.close();
    } catch (IOException ioe) {
      Check.state(false,"no exception "+ioe);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void solveInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    ensureInlineFilter2();
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] t = new float[n2][n1];
    _dlf.applyInline(null,v1,x,t);
    _fif2.applyInverse(_sigma,ds,v1,t,y);
    Array.sub(x,y,y);
  }

  protected void solveInline(
    float[][][] ds, short[][][] iw, float[][][] x, float[][][] y) 
  {
    ensureInlineFilter3();
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] t = new float[n3][n2][n1];
    //_dlf.applyInline(null,iw,x,t);
    _fif3.applyInverse(_sigma,ds,iw,x,y);
    //Array.sub(x,y,y);
  }

  protected void solveNormal(
    float[][][] ds, short[][][] iu, float[][][] x, float[][][] y) 
  {
    ensureNormalFilter3();
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] t = new float[n3][n2][n1];
    _dlf.applyNormal(null,iu,x,t);
    _fnf3.applyInverse(_sigma,ds,iu,t,y);
    Array.sub(x,y,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma;
  private String _ffile;
  private DirectionalLaplacianFilter _dlf;
  private static final int _fileFormat = 1;
  private static FactoredFilter2 _fif2; // 2-D inline filter
  private static FactoredFilter2 _fnf2; // 2-D normal filter
  private static FactoredFilter3 _fif3; // 3-D inline filter
  private static FactoredFilter3 _fnf3; // 3-D normal filter
  private static UnitSphereSampling _uss16 = new UnitSphereSampling(16);

  private void ensureInlineFilter2() {
    if (_fif2==null)
      _fif2 = new FactoredFilter2(FactoredFilter2.Type.INLINE);
  }

  private void ensureNormalFilter2() {
    if (_fnf2==null)
      _fnf2 = new FactoredFilter2(FactoredFilter2.Type.NORMAL);
  }

  private void ensureInlineFilter3() {
    if (_fif3==null) {
      if (_ffile==null) {
        _fif3 = new FactoredFilter3(FactoredFilter3.Type.INLINE);
      } else {
        _fif3 = new FactoredFilter3(FactoredFilter3.Type.INLINE,_ffile);
      }
    }
  }

  private void ensureNormalFilter3() {
    if (_fnf3==null)
      if (_ffile==null) {
        _fnf3 = new FactoredFilter3(FactoredFilter3.Type.NORMAL);
      } else {
        _fnf3 = new FactoredFilter3(FactoredFilter3.Type.NORMAL,_ffile);
      }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

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
      CausalFilter cf = new CausalFilter(_lag1,_lag2);

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
      float[][][] atable = (type==Type.INLINE)?_atableInline:_atableNormal;

      // Tabulated factors for all angles theta and half-widths sigma.
      // Note that for normal filters we specify u1 = v2 and u2 = v1.
      // Because vectors u and v are orthogonal, we might more properly
      // use u2 = -v1 = -sin(theta). However, during table lookup we will 
      // use u2 = v1 = sin(theta) for both inline and normal filters, so
      // we need to do the same here.
      for (int isigma=0; isigma<_nsigma; ++isigma) {
        float sigma = _fsigma+isigma*_dsigma;
        float scale = 2.0f/(sigma*sigma);
        for (int itheta=0; itheta<_ntheta; ++itheta) {
          float theta = _ftheta+itheta*_dtheta;
          Array.mul(scale,t,r);
          float v1 = sin(theta);
          float v2 = cos(theta);
          if (type==Type.INLINE) {
            dlf.applyInline(1.0f,v1,v2,t,r);
          } else {
            dlf.applyNormal(1.0f,v2,v1,t,r); // use u2 = v1 (not -v1)!
          }
          cf.factorWilsonBurg(100,0.000001f,r);
          atable[isigma][itheta] = cf.getA();
        }
      }

      _type = type;
      trace("...  done.");
    }

    void applyForward(
      float sigma, float[][] ds, float[][] v1, float[][] x, float[][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?_atableInline:_atableNormal;
      A2 a2 = new A2(atable,sigma,ds,v1);
      _lcf.apply(a2,x,y);
      _lcf.applyTranspose(a2,y,y);
    }

    void applyInverse(
      float sigma, float[][] ds, float[][] v1, float[][] x, float[][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?_atableInline:_atableNormal;
      A2 a2 = new A2(atable,sigma,ds,v1);
      _lcf.applyInverseTranspose(a2,x,y);
      _lcf.applyInverse(a2,y,y);
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
        if (sigma<_sigmaMin) sigma = _sigmaMin;
        if (sigma>_sigmaMax) sigma = _sigmaMax;
        float s = (sigma-_fsigma)*_ssigma;
        float t = (theta-_ftheta)*_stheta;
        int is = (int)s;
        int it = (int)t;
        float s1 = s-(float)is;
        float t1 = t-(float)it;
        float s0 = 1.0f-s1;
        float t0 = 1.0f-t1;
        float[] a00 = _at[is  ][it  ];
        float[] a01 = _at[is  ][it+1];
        float[] a10 = _at[is+1][it  ];
        float[] a11 = _at[is+1][it+1];
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
    private static final float _sigmaMin =  0.1f;
    private static final float _sigmaMax = 20.0f;
    private static final int _nsigma = 20;
    private static final float _fsigma = _sigmaMin;
    private static final float _dsigma = (_sigmaMax-_sigmaMin) /
                                         (float)(_nsigma-1);
    private static final float _ssigma = (1.0f-4.0f*FLT_EPSILON)/_dsigma;

    // Sampling of theta for tabulated filter coefficients.
    private static final float _thetaMin = -0.5f*FLT_PI;
    private static final float _thetaMax =  0.5f*FLT_PI;
    private static final int _ntheta = 21;
    private static final float _ftheta = _thetaMin;
    private static final float _dtheta = (_thetaMax-_thetaMin) /
                                         (float)(_ntheta-1);
    private static final float _stheta = (1.0f-4.0f*FLT_EPSILON)/_dtheta;

    // Tables of coefficients for both inline and normal filters.
    private static float[][][] _atableInline = new float[_nsigma][_ntheta][];
    private static float[][][] _atableNormal = new float[_nsigma][_ntheta][];

    // Lags for minimum-phase factors, with the following stencil:
    //   lag1 =    1  0 -1 -2 -3 -4
    //   --------------------------
    //   lag2 = 0: x  x
    //          1: x  x  x  x  x  x
    private static int _mlag = 4;
    private static int _nlag = 4+_mlag;
    private static int[] _lag1 = new int[_nlag];
    private static int[] _lag2 = new int[_nlag];
    private static LocalCausalFilter _lcf;
    static {
      for (int ilag=0; ilag<_nlag; ++ilag) {
        _lag1[ilag] = (ilag<=1)?ilag:ilag-2-_mlag;
        _lag2[ilag] = (ilag<=1)?0:1;
      }
      _lcf = new LocalCausalFilter(_lag1,_lag2);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  // A filter approximated with minimum-phase factors. Factors are 
  // tabulated as a function of sigma and a unit-sphere sampling index.
  private static class FactoredFilter3 {

    static enum Type {
      INLINE,
      NORMAL
    };

    FactoredFilter3(Type type) {
      trace("FactoredFilter3: begin make ...");
      _type = type;

      // A causal filter to with same lags as the local causal filter.
      CausalFilter cf = new CausalFilter(_lag1,_lag2,_lag3);

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
      float[][][] atable = (type==Type.INLINE)?_atableInline:_atableNormal;

      // Filters for all half-widths sigma and unit-vector indices.
      trace("_nsigma="+_nsigma+" _nvec="+_nvec);
      for (int isigma=0; isigma<_nsigma; ++isigma) {
        float sigma = _fsigma+isigma*_dsigma;
        float scale = 2.0f/(sigma*sigma);
        for (int ivec=0; ivec<_nvec; ++ivec) {
          trace("isigma="+isigma+" ivec="+ivec);
          Array.mul(scale,t,r);
          float[] v = _uss.getPoint(1+ivec);
          float v1 = v[2];
          float v2 = v[1];
          float v3 = v[0];
          if (type==Type.INLINE) {
            dlf.applyInline(1.0f,v1,v2,v3,t,r);
          } else {
            dlf.applyNormal(1.0f,v1,v2,v3,t,r);
          }
          trace("v1="+v1+" v2="+v2+" v3="+v3);
          Array.dump(r);
          cf.factorWilsonBurg(100,0.000001f,r);
          atable[isigma][ivec] = cf.getA();
          //dumpA(atable[isigma][ivec]);
        }
      }
      trace("...  done.");
    }

    FactoredFilter3(Type type, String ffile) {
      trace("FactoredFilter3: begin load ...");
      _type = type;

      try {
        ArrayFile af = new ArrayFile(ffile,"r");

        // Ensure known file format.
        int fileFormat = af.readInt();
        Check.state(_fileFormat==fileFormat,"known file format");

        // Find filter with specified type.
        int itype = itype(type);
        int nbyte = nbyte(type);
        boolean found = false;
        while (!found) {
          int itypeRead = af.readInt();
          int nbyteRead = af.readInt();
          found = (itype==itypeRead && nbyte==nbyteRead);
          if (!found)
            af.skipBytes(nbyteRead);
        }

        // Read arrays of filter coefficients.
        float[][][] atable = (type==Type.INLINE)?_atableInline:_atableNormal;
        for (int isigma=0; isigma<_nsigma; ++isigma) {
          for (int ivec=0; ivec<_nvec; ++ivec) {
            atable[isigma][ivec] = new float[_nlag];
            af.readFloats(atable[isigma][ivec]);
          }
        }

        af.close();
      } catch (IOException ioe) {
        Check.state(false,"loaded filters without exception "+ioe);
      }
      trace("...  done.");
    }

    void save(ArrayFile af) {
      trace("FactoredFilter3: begin save ...");

      try {
        af.writeInt(itype(_type));
        af.writeInt(nbyte(_type));
        float[][][] atable = (_type==Type.INLINE)?_atableInline:_atableNormal;
        for (int isigma=0; isigma<_nsigma; ++isigma) {
          for (int ivec=0; ivec<_nvec; ++ivec) {
            af.writeFloats(atable[isigma][ivec]);
          }
        }
      } catch (IOException ioe) {
        Check.state(false,"saved filters without exception "+ioe);
      }
      trace("...  done.");
    }

    void applyForward(
      float sigma, float[][][] ds, short[][][] iw, 
      float[][][] x, float[][][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?_atableInline:_atableNormal;
      A3 a3 = new A3(atable,sigma,ds,iw);
      _lcf.apply(a3,x,y);
      _lcf.applyTranspose(a3,y,y);
    }

    void applyInverse(
      float sigma, float[][][] ds, short[][][] iw, 
      float[][][] x, float[][][] y) 
    {
      float[][][] atable = (_type==Type.INLINE)?_atableInline:_atableNormal;
      A3 a3 = new A3(atable,sigma,ds,iw);
      _lcf.applyInverseTranspose(a3,x,y);
      _lcf.applyInverse(a3,y,y);
    }

    private Type _type; // filter type, inline or normal

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
        int iw = abs(_iw[i3][i2][i1]);
        int ia = _ia[iw];
        int ib = _ib[iw];
        int ic = _ic[iw];
        float wa = _wa[iw];
        float wb = _wb[iw];
        float wc = _wc[iw];
        float sigma = _sigma;
        if (_ds!=null) sigma *= _ds[i3][i2][i1];
        if (sigma<_sigmaMin) sigma = _sigmaMin;
        if (sigma>_sigmaMax) sigma = _sigmaMax;
        float s = (sigma-_fsigma)*_ssigma;
        int is = (int)s;
        float s1 = s-(float)is;
        float s0 = 1.0f-s1;
        float[][] a0 = _at[is  ];
        float[][] a1 = _at[is+1];
        float[] a0a = a0[ia];
        float[] a0b = a0[ib];
        float[] a0c = a0[ic];
        float[] a1a = a1[ia];
        float[] a1b = a1[ib];
        float[] a1c = a1[ic];
        //
        if (i1==52 && i2==52 && i3==52) {
          trace("s0="+s0+" s1="+s1);
          trace("ia="+ia+" ib="+ib+" ic="+ic);
          trace("wa="+wa+" wb="+wb+" wc="+wc);
          dumpA(a1[ia]);
        }
        //
        int n = a0a.length;
        for (int j=0; j<n; ++j)
          a[j] = s0*(wa*a0a[j]+wb*a0b[j]+wc*a0c[j]) +
                 s1*(wa*a1a[j]+wb*a1b[j]+wc*a1c[j]);
      }
      private float _sigma;
      private float[][][] _at;
      private float[][][] _ds;
      private short[][][] _iw;
    }

    // Sampling of sigma for tabulated filter coefficients.
    private static final float _sigmaMin =  0.1f;
    private static final float _sigmaMax = 16.0f;
    private static final int _nsigma = 2;
    private static final float _fsigma = _sigmaMin;
    private static final float _dsigma = (_sigmaMax-_sigmaMin) /
                                         (float)(_nsigma-1);
    private static final float _ssigma = (1.0f-4.0f*FLT_EPSILON)/_dsigma;

    // Sampling of the unit sphere for tabulated filter coefficients.
    // This sampling must be coarser (using fewer bits) than the 16-bit
    // sampling used to encode unit vectors.
    private static final UnitSphereSampling _uss = new UnitSphereSampling(7);
    private static final int _nvec = _uss.getMaxIndex();

    // Tables of coefficients for both inline and normal filters.
    private static float[][][] _atableInline = new float[_nsigma][_nvec][];
    private static float[][][] _atableNormal = new float[_nsigma][_nvec][];

    // Precomputed weights for linear interpolation within triangles
    // of the coarse unit-sphere sampling. Because unit-vectors are
    // quantized to 16 bits, only these weights will be needed to
    // interpolate filter coefficients.
    private static int _nw = _uss16.getMaxIndex(); // number of weights
    private static int[] _ia = new int[1+_nw]; // _ia[0] unused
    private static int[] _ib = new int[1+_nw]; // _ib[0] unused
    private static int[] _ic = new int[1+_nw]; // _ic[0] unused
    private static float[] _wa = new float[1+_nw]; // _wa[0] unused
    private static float[] _wb = new float[1+_nw]; // _wb[0] unused
    private static float[] _wc = new float[1+_nw]; // _wc[0] unused

    // Lags for minimum-phase factors with the following stencil:
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
    private static int _mlag = 4;
    private static int _nlag = 2+3*(2+_mlag)+_mlag*(_mlag+1+_mlag); // = 56
    private static int[] _lag1 = new int[_nlag];
    private static int[] _lag2 = new int[_nlag];
    private static int[] _lag3 = new int[_nlag];
    private static LocalCausalFilter _lcf;

    // Initialization of static fields.
    static {
      //trace("nlag="+_nlag);
      for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
        int jlag2 = (ilag3==0)?0:-_mlag;
        int klag2 = 1;
        for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
          int jlag1 = (ilag3==0 && ilag2==0)?0:-_mlag;
          int klag1 = (ilag3==0 || ilag2>=0)?1: _mlag;
          for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
            _lag1[ilag] = ilag1;
            _lag2[ilag] = ilag2;
            _lag3[ilag] = ilag3;
            //trace("ilag="+ilag+" lag1="+ilag1+" ilag2="+ilag2+" ilag3="+ilag3);
          }
        }
      }
      for (int iw=1; iw<_nw; ++iw) {
        float[] wi = _uss16.getPoint(iw);
        int[] iabc = _uss.getTriangle(wi);
        float[] wabc = _uss.getWeights(wi,iabc);
        _ia[iw] = iabc[0]-1; // subtract 1 because of
        _ib[iw] = iabc[1]-1; // zero-based indexing
        _ic[iw] = iabc[2]-1; // in arrays of filters
        _wa[iw] = wabc[0];
        _wb[iw] = wabc[1];
        _wc[iw] = wabc[2];
      }
      _lcf = new LocalCausalFilter(_lag1,_lag2,_lag3);
    }

    private static int itype(Type type) {
      return 3*10+type.ordinal();
    }

    private static int nbyte(Type type) {
      return 4*_nlag*_nvec*_nsigma;
    }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

    static void dumpA(float[] a) {
      int n = _mlag+1+_mlag;
      float[][] b = new float[2*n][n];
      for (int ilag=0; ilag<_nlag; ++ilag) {
        int i1 = n-1-(_lag1[ilag]+_mlag);
        int i2 = n-1-(_lag2[ilag]+_mlag);
        int i3 = _lag3[ilag];
        b[i3*n+i2][i1] = a[ilag];
      }
      edu.mines.jtk.mosaic.SimplePlot.asPixels(b);
    }
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  public static void main(String[] args) {
    //testFactoredFilter3();
    test3();
  }

  private static void testFactoredFilter3() {
    String ffile = "filters.dat";
    //FactoredFilter2 ff2 = new FactoredFilter2(FactoredFilter2.Type.INLINE);
    //FactoredFilter3 ff3 = new FactoredFilter3(FactoredFilter3.Type.INLINE);
    FactoredFilter3 ff3 = new FactoredFilter3(
      FactoredFilter3.Type.INLINE,ffile);
    try {
      ArrayFile af = new ArrayFile(ffile,"rw");
      af.writeInt(_fileFormat);
      ff3.save(af);
      af.close();
    } catch (IOException ioe) {
      Check.state(false,"no exception "+ioe);
    }
  }

  private static void test3() {
    float sigma = 20.0f;
    String ffile = "filters.dat";
    /*
    LocalDiffusionFilterMp ldf = new LocalDiffusionFilterMp(sigma);
    ldf.save(ffile);
    */
    LocalDiffusionFilterMp ldf = new LocalDiffusionFilterMp(sigma,ffile);
    int n1 = 101;
    int n2 = 101;
    int n3 = 101;
    float[][][] x = Array.randfloat(n1,n2,n3);
    float[][][] y = Array.zerofloat(n1,n2,n3);
    short[][][] iw = Array.fillshort((short)1,n1,n2,n3);
    ldf.applyInlinePass(null,iw,x,y);
  }
}
