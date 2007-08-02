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

// FOR EXPERIMENTS ONLY!
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.sgl.test.*;

/**
 * Local anisotropic diffusion filter via minimum-phase factors.
 * Uses tables of minimum-phase filters that are applied recursively.
 * These filters require no iteration, but cannot be parallelized.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.08.02
 */
public class LocalDiffusionFilterMp extends LocalDiffusionFilter {

  /**
   * Sets the name of the file from which to read pre-computed filters.
   * This method checks the version number of the specified file (and 
   * ensures that the file is readable), but does not actually read any 
   * filters. Filters are read from the file only when needed.
   * <p>
   * This method is static and is typically called only once if at all. 
   * Tables of pre-computed filters are shared for all instances of this 
   * class. Therefore, this method forgets any previously read filters,
   * and filters will then be read again when needed.
   * @param ffile name of file containing filters.
   */
  public static void setFiltersFile(String fileName) {
    checkFileVersionNumber(fileName);
    _ffile = fileName;
    _ff2 = null;
    _ff3 = null;
  }

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   */
  public LocalDiffusionFilterMp(double sigma) {
    super(sigma);
    _sigma = (float)sigma;
  }

  /**
   * Saves filter coefficients to a file with specified name.
   * @param ffile name of file to contain filter coefficents.
   */
  public void save(String ffile) {
    loadFactoredFilter2();
    loadFactoredFilter3();
    try {
      ArrayFile af = new ArrayFile(ffile,"rw");
      af.writeInt(_ffileVersionNumber);
      _ff2.save(af);
      _ff3.save(af);
      af.close();
    } catch (IOException ioe) {
      Check.state(false,"no exception "+ioe);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void solveLinear(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    loadFactoredFilter2();
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] t = new float[n2][n1];
    _dlf.applyLinear(null,v1,x,t);
    _ff2.applyInverse(_sigma,ds,v1,t,y);
    Array.sub(x,y,y);
  }

  protected void solveLinear(
    byte[][][] is, short[][][] iw, float[][][] x, float[][][] y) 
  {
    loadFactoredFilter3();
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] t = new float[n3][n2][n1];
    _dlf.applyLinear(null,iw,x,t);
    _ff3.applyInverse(false,_sigma,is,iw,t,y);
    Array.sub(x,y,y);
  }

  protected void solvePlanar(
    byte[][][] is, short[][][] iu, float[][][] x, float[][][] y) 
  {
    loadFactoredFilter3();
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;

    // For all constant-i3 slices, apply 2-D filters for in-plane vectors v.
    float[][] v1 = new float[n2][n1];
    float[][] ds = (is!=null)?new float[n2][n1]:null;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        if (is!=null) {
          float[] ds32 = (is!=null)?ds[i2]:null;
          byte[] is32 = (is!=null)?is[i3][i2]:null;
          for (int i1=0; i1<n1; ++i1) {
            ds32[i1] = _sigma*0.01f*(float)is32[i1];
          }
        }
        short[] iu32 = iu[i3][i2];
        float[] v132 = v1[i2];
        for (int i1=0; i1<n1; ++i1) {
          int iv = _iv[abs(iu32[i1])];
          float[] vi = _uss16.getPoint(iv);
          v132[i1] = vi[2];
        }
      }
      solveLinear(ds,v1,x[i3],y[i3]);
    }

    // Apply 3-D filter for in-plane vectors w orthogonal to v.
    float[][][] t = new float[n3][n2][n1];
    _dlf.applyPlanar(null,iu,y,t);
    _ff3.applyInverse(true,_sigma,is,iu,t,y);
    Array.sub(x,y,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma;

  // All instances of this class use a single directional Laplacian filter.
  // The half-width sigma = sqrt(2) for this filter compensates for the 
  // scaling by (sigma*sigma)/2 when it is applied.
  private static DirectionalLaplacianFilter _dlf =
    new DirectionalLaplacianFilter(sqrt(2.0f));

  // Factored filters shared by all instances of this class. If the file 
  // name is not null, filters will be read from that file when needed.
  private static final int _ffileVersionNumber = 1;
  private static String _ffile;
  private static FactoredFilter2 _ff2;
  private static FactoredFilter3 _ff3;

  // Mapping of normal vectors u to orthogonal in-plane vectors v and w.
  private static UnitSphereSampling _uss16 = new UnitSphereSampling(16);
  private static int _nu = _uss16.getMaxIndex();
  private static int[] _iv = new int[1+_nu]; // _iv[0] unused
  private static int[] _iw = new int[1+_nu]; // _iw[0] unused
  static {
    for (int iu=1; iu<=_nu; ++iu) {
      float[] ui = _uss16.getPoint(iu);
      float u1 = ui[2];
      float u2 = ui[1];
      float u3 = ui[0];
      float u12 = sqrt(u1*u1+u2*u2);
      float v1 = (u12>0.0f)?-u2/u12:0.0f;
      float v2 = (u12>0.0f)? u1/u12:1.0f;
      float v3 = 0.0f;
      float w1 = -u3*v2;
      float w2 =  u3*v1;
      float w3 = u12;
      _iv[iu] = abs(_uss16.getIndex(v3,v2,v1));
      _iw[iu] = abs(_uss16.getIndex(w3,w2,w1));
    }
  }

  private static void checkFileVersionNumber(String ffile) {
    try {
      ArrayFile af = new ArrayFile(ffile,"r");
      if (_ffileVersionNumber!=af.readInt())
        throw new RuntimeException("bad version number in file "+ffile);
      af.close();
    } catch (IOException ioe) {
      throw new RuntimeException("cannot read file "+ffile+": "+ioe);
    }
  }

  private static void loadFactoredFilter2() {
    if (_ff2==null) {
      if (_ffile==null) {
        _ff2 = new FactoredFilter2();
      } else {
        try {
          ArrayFile af = new ArrayFile(_ffile,"r");
          af.readInt(); // skip file version number, already checked
          while (af.readInt()!=2)
            af.skipBytes(af.readInt());
          af.readInt(); // skip number of bytes
          _ff2 = new FactoredFilter2(af);
          af.close();
        } catch (IOException ioe) {
          throw new RuntimeException("cannot read 2-D filters: "+ioe);
        }
      }
    }
  }

  private static void loadFactoredFilter3() {
    if (_ff3==null) {
      if (_ffile==null) {
        _ff3 = new FactoredFilter3();
      } else {
        try {
          ArrayFile af = new ArrayFile(_ffile,"r");
          af.readInt(); // skip file version number, already checked
          while (af.readInt()!=3)
            af.skipBytes(af.readInt());
          af.readInt(); // skip number of bytes
          _ff3 = new FactoredFilter3(af);
          af.close();
        } catch (IOException ioe) {
          throw new RuntimeException("cannot read 3-D filters: "+ioe);
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  // A local linear filter approximated with minimum-phase factors.
  // Factors are tabulated as a function of sigma and theta.
  private static class FactoredFilter2 {

    FactoredFilter2() {
      trace("FactoredFilter2: constructing filters ...");

      // A causal filter with same lags as the local causal filter.
      CausalFilter cf = new CausalFilter(_lag1,_lag2);

      // Arrays used to get 3x3 auto-correlations to be factored.
      float[][] r = new float[3][3];
      float[][] t = new float[3][3];
      t[1][1] = 1.0f;

      // Tabulated factors for all angles theta and half-widths sigma.
      for (int isigma=0; isigma<_nsigma; ++isigma) {
        float sigma = _fsigma+isigma*_dsigma;
        float scale = 2.0f/(sigma*sigma);
        for (int itheta=0; itheta<_ntheta; ++itheta) {
          float theta = _ftheta+itheta*_dtheta;
          Array.mul(scale,t,r);
          float v1 = sin(theta);
          float v2 = cos(theta);
          _dlf.applyLinear(1.0f,v1,v2,t,r);
          cf.factorWilsonBurg(100,0.000001f,r);
          _atable[isigma][itheta] = cf.getA();
        }
      }

      trace("...  done.");
    }

    FactoredFilter2(ArrayFile af) {
      trace("FactoredFilter2: begin load ...");
      try {
        // Read arrays of filter coefficients.
        for (int isigma=0; isigma<_nsigma; ++isigma) {
          for (int itheta=0; itheta<_ntheta; ++itheta) {
            _atable[isigma][itheta] = new float[_nlag];
            af.readFloats(_atable[isigma][itheta]);
          }
        }
      } catch (IOException ioe) {
        Check.state(false,"loaded filters successfully: "+ioe);
      }
      trace("...  done.");
    }

    void save(ArrayFile af) {
      trace("FactoredFilter2: begin save ...");
      try {
        af.writeInt(2); // filter type
        af.writeInt(4*_nlag*_ntheta*_nsigma); // number of bytes
        for (int isigma=0; isigma<_nsigma; ++isigma) {
          for (int itheta=0; itheta<_ntheta; ++itheta) {
            af.writeFloats(_atable[isigma][itheta]);
          }
        }
      } catch (IOException ioe) {
        Check.state(false,"saved filters successfully: "+ioe);
      }
      trace("...  done.");
    }

    void applyInverse(
      float sigma, float[][] ds, float[][] v1, float[][] x, float[][] y) 
    {
      A2 a2 = new A2(_atable,sigma,ds,v1);
      _lcf.applyInverseTranspose(a2,x,y);
      _lcf.applyInverse(a2,y,y);
    }

    // Provides linear filter coefficients for each sample index via
    // bilinear interpolation of pre-computed filter coefficients.
    private static class A2 implements LocalCausalFilter.A2 {
      A2(float[][][] atable, float sigma, float[][] ds, float[][] v1) {
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
    private static final float _sigmaMin =  0.01f;
    private static final float _sigmaMax = 16.00f;
    private static final int _nsigma = 9;
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

    // Tables of filter coefficients.
    private static float[][][] _atable = new float[_nsigma][_ntheta][];

    // Lags for minimum-phase factors with the following stencil:
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

    FactoredFilter3() {
      trace("FactoredFilter3: begin make ...");

      // Causal filter with same lags as our local causal filter.
      CausalFilter cf = new CausalFilter(_lag1,_lag2,_lag3);

      // Arrays used to compute 3x3x3 auto-correlations to be factored.
      float[][][] r = new float[3][3][3];
      float[][][] t = new float[3][3][3];
      t[1][1][1] = 1.0f;

      // Filters for all half-widths sigma and unit-vector indices.
      trace("nsigma="+_nsigma+" nvec="+_nvec);
      //for (int isigma=0; isigma<_nsigma; ++isigma) {
      for (int isigma=_nsigma-1; isigma>=0; --isigma) {
        float sigma = _fsigma+isigma*_dsigma;
        float scale = 2.0f/(sigma*sigma);
        for (int ivec=0; ivec<_nvec; ++ivec) {
          trace("isigma="+isigma+" ivec="+ivec);
          Array.mul(scale,t,r);
          float[] w = _uss.getPoint(1+ivec); // {wx,wy,wz}
          float w1 = w[2]; // w1 = wz
          float w2 = w[1]; // w2 = wy
          float w3 = w[0]; // w3 = wx
          _dlf.applyLinear(1.0f,w1,w2,w3,t,r);
          trace("  sigma="+sigma+" w1="+w1+" w2="+w2+" w3="+w3);
          cf.factorWilsonBurg(100,0.000001f,r);
          _atable[isigma][ivec] = cf.getA();
        }
      }
      trace("...  done.");
    }

    FactoredFilter3(ArrayFile af) {
      trace("FactoredFilter3: begin load ...");
      try {
        // Read arrays of filter coefficients.
        for (int isigma=0; isigma<_nsigma; ++isigma) {
          for (int ivec=0; ivec<_nvec; ++ivec) {
            _atable[isigma][ivec] = new float[_nlag];
            af.readFloats(_atable[isigma][ivec]);
          }
        }
      } catch (IOException ioe) {
        Check.state(false,"loaded filters successfully: "+ioe);
      }
      trace("...  done.");
    }

    void save(ArrayFile af) {
      trace("FactoredFilter3: begin save ...");
      try {
        af.writeInt(3); // filter type
        af.writeInt(4*_nlag*_nvec*_nsigma); // number of bytes
        for (int isigma=0; isigma<_nsigma; ++isigma) {
          for (int ivec=0; ivec<_nvec; ++ivec) {
            af.writeFloats(_atable[isigma][ivec]);
          }
        }
      } catch (IOException ioe) {
        Check.state(false,"saved filters successfully: "+ioe);
      }
      trace("...  done.");
    }

    void applyInverse(
      boolean mapuw, float sigma, 
      byte[][][] is, short[][][] iw, 
      float[][][] x, float[][][] y) 
    {
      A3 a3 = new A3(mapuw,_atable,sigma,is,iw);
      _lcf.applyInverseTranspose(a3,x,y);
      _lcf.applyInverse(a3,y,y);
    }

    // Provides filter coefficients for each sample index via
    // linear interpolation of pre-computed filter coefficients.
    private static class A3 implements LocalCausalFilter.A3 {
      A3(boolean mapuw, float[][][] atable, 
         float sigma, byte[][][] is, short[][][] iw) 
      {
        _mapuw = mapuw;
        _sigma = sigma;
        _at = atable;
        _is = is;
        _iw = iw;
      }
      public void get(int i1, int i2, int i3, float[] a) {
        int iw = abs(_iw[i3][i2][i1]);
        if (_mapuw) 
          iw = LocalDiffusionFilterMp._iw[iw];
        int ia = _ia[iw];
        int ib = _ib[iw];
        int ic = _ic[iw];
        float wa = _wa[iw];
        float wb = _wb[iw];
        float wc = _wc[iw];
        float sigma = _sigma;
        if (_is!=null) 
          sigma *= 0.01f*(float)_is[i3][i2][i1];
        if (sigma<_sigmaMin) sigma = _sigmaMin;
        if (sigma>_sigmaMax) sigma = _sigmaMax;
        float s = (sigma-_fsigma)*_ssigma;
        int is = (int)s;
        float s1 = s-(float)is;
        float s0 = 1.0f-s1;
        float[][] as = (s0>s1)?_at[is  ]:_at[is+1];
        float[] asv = as[ia];
        if (wa<wb && wc<wb) asv = as[ib];
        if (wb<wc && wa<wc) asv = as[ic];
        Array.copy(asv,a);
        /*
        float[][] a0 = _at[is  ];
        float[][] a1 = _at[is+1];
        float[] a0a = a0[ia];
        float[] a0b = a0[ib];
        float[] a0c = a0[ic];
        float[] a1a = a1[ia];
        float[] a1b = a1[ib];
        float[] a1c = a1[ic];
        float[] ai = a0a;
        int n = a0a.length;
        for (int j=0; j<n; ++j)
          a[j] = s0*(wa*a0a[j]+wb*a0b[j]+wc*a0c[j]) +
                 s1*(wa*a1a[j]+wb*a1b[j]+wc*a1c[j]);
        */
      }
      private boolean _mapuw;
      private float _sigma;
      private float[][][] _at;
      private byte[][][] _is;
      private short[][][] _iw;
    }

    // Sampling of sigma for tabulated filter coefficients.
    private static final float _sigmaMin =  0.01f;
    private static final float _sigmaMax = 16.00f;
    private static final int _nsigma = 9;
    private static final float _fsigma = _sigmaMin;
    private static final float _dsigma = (_sigmaMax-_sigmaMin) /
                                         (float)(_nsigma-1);
    private static final float _ssigma = (1.0f-4.0f*FLT_EPSILON)/_dsigma;

    // Sampling of the unit sphere for tabulated filter coefficients.
    // This sampling must be coarser (using fewer bits) than the 16-bit
    // sampling used to encode unit vectors.
    private static final UnitSphereSampling _uss = new UnitSphereSampling(10);
    private static final int _nvec = _uss.getMaxIndex();

    // Table of filter coefficients.
    private static float[][][] _atable = new float[_nsigma][_nvec][];

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
    private static int _nlag = 2+3*(2+_mlag)+_mlag*(1+2*_mlag); // = 56
    private static int[] _lag1 = new int[_nlag];
    private static int[] _lag2 = new int[_nlag];
    private static int[] _lag3 = new int[_nlag];
    private static LocalCausalFilter _lcf;

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

    // Initialization of static fields.
    static {

      // Lags for filters.
      for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
        int jlag2 = (ilag3==0)?0:-_mlag;
        int klag2 = 1;
        for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
          int jlag1 = (ilag3==0 && ilag2==0)?0:-_mlag;
          int klag1 = (ilag3==1 && ilag2<0)?_mlag:1;
          for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
            _lag1[ilag] = ilag1;
            _lag2[ilag] = ilag2;
            _lag3[ilag] = ilag3;
          }
        }
      }
      _lcf = new LocalCausalFilter(_lag1,_lag2,_lag3);

      // Filter indices and weights for 16-bit unit sphere sample indices.
      for (int iw=1; iw<=_nw; ++iw) {
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
    }

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

    private static void checkA(CausalFilter cf, float[][][] r) {
      float[][][] t = new float[21][21][21];
      float[][][] t2 = new float[21][21][21];
      t[10][10][10] = 1.0f;
      cf.apply(t,t2);
      cf.applyTranspose(t2,t);
      float[][][] s = Array.copy(3,3,3,9,9,9,t);
      Array.dump(r);
      Array.dump(s);
    }
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // FOR EXPERIMENTS ONLY.
  // 2-D filters are efficient and accurate enough for practical use.
  // Current status for 3-D filters:
  // Approximations to planar filters require too many lags to be useful. 
  // Line filters seem to require fewer lags, although these too show
  // errors. For example the width of the notch is too large, and amplitudes
  // away from the notch are not constant. They exhibit oscillations.
  // Increasing the number of lags improves the accuracy of these filters,
  // but makes them more costly.
  // An alterna

  public static void main(String[] args) {
    testFactorizations();
    testFactorizations2();
  }

  // Test harness for experimenting with minimum-phase factors. 
  private static void testFactorizations() {
    int n = 105; // number of samples
    int k = n/2; // index of central sample
    float[][][] x = new float[n][n][n];
    float[][][] y = new float[n][n][n];
    float[][][] z = new float[n][n][n];
    x[k][k][k] = 1.0f; // input is unit impulse

    // Unit vector.
    float theta = 0.0f*FLT_PI/180.0f;
    float   phi = 0.0f*FLT_PI/180.0f;
    float v1 = cos(theta);
    float v2 = sin(phi)*sin(theta);
    float v3 = cos(phi)*sin(theta);

    // Diffusion coefficients for linear or planar filter.
    boolean linear = false;
    float d11,d22,d33,d12,d13,d23;
    if (linear) {
      d11 = v1*v1;
      d22 = v2*v2;
      d33 = v3*v3;
      d12 = v1*v2;
      d13 = v1*v3;
      d23 = v2*v3;
    } else {
      d11 = 1.0f-v1*v1;
      d22 = 1.0f-v2*v2;
      d33 = 1.0f-v3*v3;
      d12 =     -v1*v2;
      d13 =     -v1*v3;
      d23 =     -v2*v3;
    }
    
    // Numerator y = A'Ax.
    applyDiffusionFilter27(d11,d12,d13,d22,d23,d33,x,y);

    // Denominator z = (eps*I+A'A)x.
    float sigma = 16.0f;
    float eps = 2.0f/(sigma*sigma);
    Array.mul(eps,x,z);
    applyDiffusionFilter27(d11,d12,d13,d22,d23,d33,x,z);

    // Print eps*I+A'A.
    float[][][] r = Array.copy(3,3,3,k-1,k-1,k-1,z);
    Array.dump(r);

    // Plot ratio Ay/Az which equals the desired amplitude spectrum.
    float[][][] ay = amplitude(y);
    float[][][] az = amplitude(z);
    plot3d(Array.sub(1.0f,Array.div(ay,az)));

    // Minimum-phase causal filter.
    int[][] lags = (linear)?makeLagsLinear27():makeLagsPlanar27();
    int[] lag1 = lags[0], lag2 = lags[1], lag3 = lags[2];
    CausalFilter cf = new CausalFilter(lag1,lag2,lag3);
    cf.factorWilsonBurg(100,0.000001f,r);
    float[] a = cf.getA();
    int nlag = lag1.length;
    for (int ilag=0; ilag<nlag; ++ilag) {
      System.out.println(
        "l1="+lag1[ilag]+" l2="+lag2[ilag]+" l3="+lag3[ilag]+" a="+a[ilag]);
    }

    // Use minimum-phase factors to recompute denominator z.
    cf.apply(x,z);
    cf.applyTranspose(z,z);

    // Plot amplitude spectrum.
    az = amplitude(z);
    plot3d(Array.sub(1.0f,Array.div(ay,az)));

    // Print approximation to eps*I+A'A implied by causal filter.
    r = Array.copy(7,7,7,k-3,k-3,k-3,z);
    Array.dump(r);
  }

  // Another test harness for experimenting with minimum-phase factors. 
  // This one approximates a planar filter by cascading two linear filters.
  private static void testFactorizations2() {
    int n = 105; // number of samples
    int k = n/2; // index of central sample
    float[][][] x = new float[n][n][n];
    float[][][] y1 = new float[n][n][n];
    float[][][] y2 = new float[n][n][n];
    float[][][] za1 = new float[n][n][n];
    float[][][] za2 = new float[n][n][n];
    float[][][] zb1 = new float[n][n][n];
    float[][][] zb2 = new float[n][n][n];
    x[k][k][k] = 1.0f; // input is unit impulse

    // Filter epsilon.
    float sigma = 16.0f;
    float eps = 2.0f/(sigma*sigma);

    // Unit vector.
    float theta =  0.0f*FLT_PI/180.0f;
    float   phi =  0.0f*FLT_PI/180.0f;
    float u1 = cos(theta);
    float u2 = sin(phi)*sin(theta);
    float u3 = cos(phi)*sin(theta);
    float v1,v2,v3;
    float u12 = sqrt(u1*u1+u2*u2);
    if (u12>0.0f) {
      v1 = -u2/u12;
      v2 =  u1/u12;
    } else {
      v1 = 0.0f;
      v2 = 1.0f;
    }
    v3 = 0.0f;
    float w1 = -v2*u3;
    float w2 =  v1*u3;
    float w3 =  u12;
    trace("u1="+u1+" v2="+v2+" w3="+w3);

    // 1st (v) linear filter.
    float d11,d22,d33,d12,d13,d23;
    d11 = v1*v1;
    d22 = v2*v2;
    d33 = v3*v3; // = 0
    d12 = v1*v2;
    d13 = v1*v3; // = 0
    d23 = v2*v3; // = 0
    applyDiffusionFilter9(d11,d12,d22,x[k],y1[k]);
    Array.mul(eps,x,za1);
    applyDiffusionFilter9(d11,d12,d22,x[k],za1[k]);
    int[][] lags = makeLagsLinear9();
    int[] lag1,lag2,lag3;
    lag1 = lags[0];
    lag2 = lags[1];
    CausalFilter cf = new CausalFilter(lag1,lag2);
    float[][] r1 = Array.copy(3,3,k-1,k-1,za1[k]);
    Array.dump(r1);
    cf.factorWilsonBurg(100,0.000001f,r1);
    cf.apply(x[k],zb1[k]);
    cf.applyTranspose(zb1[k],zb1[k]);
    float[][][] a1 = Array.sub(1.0f,Array.div(amplitude(y1),amplitude(za1)));
    float[][][] b1 = Array.sub(1.0f,Array.div(amplitude(y1),amplitude(zb1)));

    // 2nd (w) linear filter.
    d11 = w1*w1;
    d22 = w2*w2;
    d33 = w3*w3;
    d12 = w1*w2;
    d13 = w1*w3;
    d23 = w2*w3;
    applyDiffusionFilter27(d11,d12,d13,d22,d23,d33,x,y2);
    Array.mul(eps,x,za2);
    applyDiffusionFilter27(d11,d12,d13,d22,d23,d33,x,za2);
    lags = makeLagsLinear27();
    lag1 = lags[0];
    lag2 = lags[1];
    lag3 = lags[2];
    cf = new CausalFilter(lag1,lag2,lag3);
    Array.mul(eps,x,zb2);
    applyDiffusionFilter27(d11,d12,d13,d22,d23,d33,x,zb2);
    float[][][] r2 = Array.copy(3,3,3,k-1,k-1,k-1,zb2);
    Array.dump(r2);
    cf.factorWilsonBurg(100,0.000001f,r2);
    cf.apply(x,zb2);
    cf.applyTranspose(zb2,zb2);
    float[][][] a2 = Array.sub(1.0f,Array.div(amplitude(y2),amplitude(za2)));
    float[][][] b2 = Array.sub(1.0f,Array.div(amplitude(y2),amplitude(zb2)));

    // Amplitude spectra.
    float[][][] aa = Array.mul(a1,a2);
    float[][][] bb = Array.mul(b1,b2);
    plot3d(aa);
    plot3d(bb);
    //plot3d(amplitude(za1));
    //plot3d(amplitude(zb1));
    //plot3d(amplitude(za2));
    //plot3d(amplitude(zb2));
  }

  // Computes y = y+G'DGx with 9-point (2-D) stencil.
  private static void applyDiffusionFilter9(
   float d11, float d12, float d22,
   float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.5f*(xa-xb);
        float x2 = 0.5f*(xa+xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Simple filter y += G'DGx with 27-point stencil.
  private static void applyDiffusionFilter27(
   float d11, float d12, float d13, float d22, float d23, float d33,
   float[][][] x, float[][][] y)
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float x000 = x[i3  ][i2  ][i1  ];
          float x001 = x[i3  ][i2  ][i1-1];
          float x010 = x[i3  ][i2-1][i1  ];
          float x100 = x[i3-1][i2  ][i1  ];
          float x011 = x[i3  ][i2-1][i1-1];
          float x101 = x[i3-1][i2  ][i1-1];
          float x110 = x[i3-1][i2-1][i1  ];
          float x111 = x[i3-1][i2-1][i1-1];
          //float x1 = 0.25f*(x000+x010+x100+x110-x001-x011-x101-x111);
          //float x2 = 0.25f*(x000+x001+x100+x101-x010-x011-x110-x111);
          //float x3 = 0.25f*(x000+x001+x010+x011-x100-x101-x110-x111);
          float xa = x000-x111;
          float xb = x001-x110;
          float xc = x010-x101;
          float xd = x100-x011;
          float x1 = 0.25f*(xa-xb+xc+xd);
          float x2 = 0.25f*(xa+xb-xc+xd);
          float x3 = 0.25f*(xa+xb+xc-xd);
          float y1 = d11*x1+d12*x2+d13*x3;
          float y2 = d12*x1+d22*x2+d23*x3;
          float y3 = d13*x1+d23*x2+d33*x3;
          float ya = 0.25f*(y1+y2+y3);
          float yb = 0.25f*(y1-y2+y3);
          float yc = 0.25f*(y1+y2-y3);
          float yd = 0.25f*(y1-y2-y3);
          y[i3  ][i2  ][i1  ] += ya;
          y[i3  ][i2  ][i1-1] -= yd;
          y[i3  ][i2-1][i1  ] += yb;
          y[i3-1][i2  ][i1  ] += yc;
          y[i3  ][i2-1][i1-1] -= yc;
          y[i3-1][i2  ][i1-1] -= yb;
          y[i3-1][i2-1][i1  ] += yd;
          y[i3-1][i2-1][i1-1] -= ya;
        }
      }
    }
  }

  // Experimental filter y += G'DGx with 19-point stencil.
  private static void applyDiffusionFilter19(
   float d11, float d12, float d13, float d22, float d23, float d33,
   float[][][] x, float[][][] y)
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float e11 = 0.5f*d11;
    float e22 = 0.5f*d22;
    float e33 = 0.5f*d33;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float x000 = x[i3  ][i2  ][i1  ];
          float x001 = x[i3  ][i2  ][i1-1];
          float x010 = x[i3  ][i2-1][i1  ];
          float x100 = x[i3-1][i2  ][i1  ];
          float x011 = x[i3  ][i2-1][i1-1];
          float x101 = x[i3-1][i2  ][i1-1];
          float x110 = x[i3-1][i2-1][i1  ];
          float x111 = x[i3-1][i2-1][i1-1];
          float x12 = 0.5f*(x000-x001+x010-x011);
          float x21 = 0.5f*(x000-x010+x001-x011);
          float x13 = 0.5f*(x000-x001+x100-x101);
          float x31 = 0.5f*(x000-x100+x001-x101);
          float x23 = 0.5f*(x000-x010+x100-x110);
          float x32 = 0.5f*(x000-x100+x010-x110);
          float y12 = e11*x12+d12*x21;
          float y21 = d12*x12+e22*x21;
          float y13 = e11*x13+d13*x31;
          float y31 = d13*x13+e33*x31;
          float y23 = e22*x23+d23*x32;
          float y32 = d23*x23+e33*x32;
          float y000 = 0.5f*( y12+y21+y13+y31+y23+y32);
          float y001 = 0.5f*(-y12+y21-y13+y31        );
          float y010 = 0.5f*( y12-y21        -y23+y32);
          float y100 = 0.5f*(         y13-y31+y23-y32);
          float y011 = 0.5f*(-y12-y21                );
          float y101 = 0.5f*(        -y13-y31        );
          float y110 = 0.5f*(                -y23-y32);
          y[i3  ][i2  ][i1  ] += y000;
          y[i3  ][i2  ][i1-1] += y001;
          y[i3  ][i2-1][i1  ] += y010;
          y[i3-1][i2  ][i1  ] += y100;
          y[i3  ][i2-1][i1-1] += y011;
          y[i3-1][i2  ][i1-1] += y101;
          y[i3-1][i2-1][i1  ] += y110;
        }
      }
    }
  }

  // Makes lags for 2-D linear filters with 9-point stencil.
  private static int[][] makeLagsLinear9() {
    int mlag = 4;
    int nlag = 4+mlag;
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    for (int ilag=0; ilag<nlag; ++ilag) {
      lag1[ilag] = (ilag<=1)?ilag:ilag-2-mlag;
      lag2[ilag] = (ilag<=1)?0:1;
    }
    return new int[][]{lag1,lag2};
  }

  // Makes lags for filters with 27-point stencil.
  private static int[][] makeLagsLinear27() {
    int mlag = 4;
    int nlag = 2+3*(2+mlag)+mlag*(1+2*mlag); // = 56
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    int[] lag3 = new int[nlag];
    for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
      int jlag2 = (ilag3==0)?0:-mlag;
      int klag2 = 1;
      for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
        int jlag1 = (ilag3==0 && ilag2==0)?0:-mlag;
        int klag1 = (ilag3==1 && ilag2<0)?mlag:1;
        for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
          lag1[ilag] = ilag1;
          lag2[ilag] = ilag2;
          lag3[ilag] = ilag3;
        }
      }
    }
    return new int[][]{lag1,lag2,lag3};
  }

  // Makes lags for planar filters with 27-point stencil.
  private static int[][] makeLagsPlanar27() {
    int mlag = 4;
    int nlag = 3+2*mlag+(1+2*mlag)*(1+2*mlag);
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    int[] lag3 = new int[nlag];
    for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
      int jlag2 = (ilag3==0)?0:-mlag;
      int klag2 = (ilag3==0)?mlag:1;
      for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
        int jlag1 = (ilag3==0 && ilag2==0)?0:-mlag;
        int klag1 = (ilag3==1 && ilag2>0)?1:mlag;
        for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
          lag1[ilag] = ilag1;
          lag2[ilag] = ilag2;
          lag3[ilag] = ilag3;
        }
      }
    }
    return new int[][]{lag1,lag2,lag3};
  }

  // Makes lags for planar filters with 19-point stencil.
  private static int[][] makeLagsPlanar19() {
    int mlag = 4;
    int nlag = 1+mlag+mlag*(mlag+1+mlag)+1+mlag+(1+mlag)*(mlag+1+mlag);
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    int[] lag3 = new int[nlag];
    for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
      int jlag2 = (ilag3==0)?0:-mlag;
      int klag2 = (ilag3==0)?mlag:1;
      for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
        int jlag1 = (ilag3==0 && ilag2==0)?0:-mlag;
        int klag1 = (ilag3==1 && ilag2==1)?0:mlag;
        for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
          lag1[ilag] = ilag1;
          lag2[ilag] = ilag2;
          lag3[ilag] = ilag3;
        }
      }
    }
    return new int[][]{lag1,lag2,lag3};
  }

  // Makes lags for entire NSHP stencil used to find a smaller stencil.
  private static int[][] makeLagsAll() {
    int mlag = 4;
    int nlag = 1+mlag+mlag*(mlag+1+mlag)+(mlag+1+mlag)*(mlag+1+mlag);
    int[] lag1 = new int[nlag];
    int[] lag2 = new int[nlag];
    int[] lag3 = new int[nlag];
    for (int ilag3=0,ilag=0; ilag3<2; ++ilag3) {
      int jlag2 = (ilag3==0)?0:-mlag;
      int klag2 = mlag;
      for (int ilag2=jlag2; ilag2<=klag2; ++ilag2) {
        int jlag1 = (ilag3==0 && ilag2==0)?0:-mlag;
        int klag1 = mlag;
        for (int ilag1=jlag1; ilag1<=klag1; ++ilag1,++ilag) {
          lag1[ilag] = ilag1;
          lag2[ilag] = ilag2;
          lag3[ilag] = ilag3;
        }
      }
    }
    return new int[][]{lag1,lag2,lag3};
  }

  // Computes 3-D Fourier amplitude spectrum.
  private static float[][][] amplitude(float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    n1 = FftComplex.nfftSmall(n1);
    n2 = FftComplex.nfftSmall(n2);
    n3 = FftComplex.nfftSmall(n3);
    float[][][] xr = Array.copy(n1,n2,n3,x);
    float[][][] xi = Array.zerofloat(n1,n2,n3);
    float[][][] cx = Array.cmplx(xr,xi);
    FftComplex fft1 = new FftComplex(n1);
    FftComplex fft2 = new FftComplex(n2);
    FftComplex fft3 = new FftComplex(n3);
    fft1.complexToComplex1(1,n2,n3,cx,cx);
    fft2.complexToComplex2(1,n1,n3,cx,cx);
    fft3.complexToComplex3(1,n1,n2,cx,cx);
    float[][][] ax = Array.cabs(cx);
    float[][][] a = Array.zerofloat(n1,n2,n3);
    int j1 = n1/2;
    int j2 = n2/2;
    int j3 = n3/2;
    Array.copy(n1-j1,n2-j2,n3-j3,0,0,0,ax,j1,j2,j3,a);
    Array.copy(j1,n2-j2,n3-j3,n1-j1,0,0,ax,0,j2,j3,a);
    Array.copy(n1-j1,j2,n3-j3,0,n2-j2,0,ax,j1,0,j3,a);
    Array.copy(n1-j1,n2-j2,j3,0,0,n3-j3,ax,j1,j2,0,a);
    Array.copy(j1,j2,n3-j3,n1-j1,n2-j2,0,ax,0,0,j3,a);
    Array.copy(n1-j1,j2,j3,0,n2-j2,n3-j3,ax,j1,0,0,a);
    Array.copy(j1,n2-j2,j3,n1-j1,0,n3-j3,ax,0,j2,0,a);
    Array.copy(j1,j2,j3,n1-j1,n2-j2,n3-j3,ax,0,0,0,a);
    return a;
  }

  // Plots a 3-D array.
  public static void plot3d(float[][][] x) {
    System.out.println("x min="+Array.min(x)+" max="+Array.max(x));
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    ImagePanelGroup ipg = new ImagePanelGroup(s1,s2,s3,new SimpleFloat3(x));
    ipg.setClips(0.0f,1.0f);
    ipg.setColorModel(ColorMap.JET);
    World world = new World();
    world.addChild(ipg);
    TestFrame frame = new TestFrame(world);
    frame.setVisible(true);
  }
}
