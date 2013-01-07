/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lss;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.dsp.LocalDiffusionKernel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes local semblance images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.03.07
 */
public class LocalSemblanceFilterX {

  /**
   * Types of smoothing filters used to compute semblance.
   */
  public enum Smoothing {
    NONE,
    BOXCAR,
    GAUSSIAN,
    LAPLACIAN
  }

  /**
   * 2D smoothing directions correspond to eigenvectors of tensors.
   */
  public enum Direction2 {
    U,V,UV
  }

  /**
   * 3D smoothing directions correspond to eigenvectors of tensors.
   */
  public enum Direction3 {
    U,V,W,UV,UW,VW,UVW
  }

  /**
   * Smoothing filter used to compute semblance.
   */
  public interface Smoother {
    public void apply(float[] f, float[] g);
    public void apply(
      Direction2 d, EigenTensors2 t, float[][] f, float[][] g);
    public void apply(
      Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g);
  }

  /**
   * Constructs a local semblance filter.
   * @param smoothing1 type of 1st smoothing filter.
   * @param halfWidth1 half-width of 1st smoothing filter.
   * @param smoothing2 type of 2nd smoothing filter.
   * @param halfWidth2 half-width of 2nd smoothing filter.
   */
  public LocalSemblanceFilterX(
    Smoothing smoothing1, int halfWidth1,
    Smoothing smoothing2, int halfWidth2) 
  {
    _smoother1 = makeSmoothingFilter(smoothing1,halfWidth1);
    _smoother2 = makeSmoothingFilter(smoothing2,halfWidth2);
  }

  public void smooth1(float[] f, float[] g) {
    _smoother1.apply(f,g);
  }
  public float[] smooth1(float[] f) {
    float[] g = shapeOf(f);
    smooth1(f,g);
    return g;
  }
  public void smooth1(
    Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
  {
    _smoother1.apply(d,t,f,g);
  }
  public float[][] smooth1(Direction2 d, EigenTensors2 t, float[][] f) {
    float[][] g = shapeOf(f);
    smooth1(d,t,f,g);
    return g;
  }
  public void smooth1(
    Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
  {
    _smoother1.apply(d,t,f,g);
  }
  public float[][][] smooth1(Direction3 d, EigenTensors3 t, float[][][] f) {
    float[][][] g = shapeOf(f);
    smooth1(d,t,f,g);
    return g;
  }

  public void smooth2(float[] f, float[] g) {
    _smoother2.apply(f,g);
  }
  public float[] smooth2(float[] f) {
    float[] g = shapeOf(f);
    smooth2(f,g);
    return g;
  }
  public void smooth2(
    Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
  {
    _smoother2.apply(d,t,f,g);
  }
  public float[][] smooth2(Direction2 d, EigenTensors2 t, float[][] f) {
    float[][] g = shapeOf(f);
    smooth2(orthogonal(d),t,f,g);
    return g;
  }
  public void smooth2(
    Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
  {
    _smoother2.apply(d,t,f,g);
  }
  public float[][][] smooth2(Direction3 d, EigenTensors3 t, float[][][] f) {
    float[][][] g = shapeOf(f);
    smooth2(orthogonal(d),t,f,g);
    return g;
  }

  public void semblance(float[] f, float[] s) {
    int n1 = f.length;
    float[] sn,sd;
    sn = smooth1(f);
    sn = mul(sn,sn);
    sn = smooth2(sn);
    sd = mul(f,f);
    sd = smooth1(sd);
    sd = smooth2(sd);
    for (int i1=0; i1<n1; ++i1) {
      float sni = sn[i1];
      float sdi = sd[i1];
      if (sdi<=0.0f)
        s[i1] = 0.0f;
      else if (sdi<sni)
        s[i1] = 1.0f;
      else
        s[i1] = sni/sdi;
    }
  }

  public float[] semblance(float[] f) {
    float[] s = shapeOf(f);
    semblance(f,s);
    return s;
  }

  public void semblance(
    Direction2 d, EigenTensors2 t, float[][] f, float[][] s) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] sn,sd;
    sn = smooth1(d,t,f);
    sn = mul(sn,sn);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(sn);
    sn = smooth2(d,t,sn);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(sn);
    sd = mul(f,f);
    sd = smooth1(d,t,sd);
    sd = smooth2(d,t,sd);
    //edu.mines.jtk.mosaic.SimplePlot.asPixels(sd);
    int count0 = 0;
    int count1 = 0;
    for (int i2=0; i2<n2; ++i2) {
      if (allZero(f[i2])) {
        for (int i1=0; i1<n1; ++i1) {
          s[i2][i1] = 0.0f;
        }
      } else {
        for (int i1=0; i1<n1; ++i1) {
          float sni = sn[i2][i1];
          float sdi = sd[i2][i1];
          if (sdi<=0.0f || sni<0.0f) {
            s[i2][i1] = 0.0f;
            ++count0;
          } if (sdi<sni) {
            s[i2][i1] = 1.0f;
            ++count1;
          } else {
            s[i2][i1] = sni/sdi;
          }
        }
      }
    }
    trace("semblance2: count0="+count0+" count1="+count1);
  }

  public float[][] semblance(Direction2 d, EigenTensors2 t, float[][] f) {
    float[][] s = shapeOf(f);
    semblance(d,t,f,s);
    return s;
  }

  public void semblance(
    Direction3 d, EigenTensors3 t, float[][][] f, float[][][] s) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] sn,sd;
    sn = smooth1(d,t,f);
    sn = mul(sn,sn);
    sn = smooth2(d,t,sn);
    sd = mul(f,f);
    sd = smooth1(d,t,sd);
    sd = smooth2(d,t,sd);
    int count0 = 0;
    int count1 = 0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        if (allZero(f[i3][i2])) {
          for (int i1=0; i1<n1; ++i1) {
            s[i3][i2][i1] = 0.0f;
          }
        } else {
          for (int i1=0; i1<n1; ++i1) {
            float sni = sn[i3][i2][i1];
            float sdi = sd[i3][i2][i1];
            if (sdi<=0.0f || sni<0.0f) {
              s[i3][i2][i1] = 0.0f;
              ++count0;
            } else if (sdi<sni) {
              s[i3][i2][i1] = 1.0f;
              ++count1;
            } else {
              s[i3][i2][i1] = sni/sdi;
            }
          }
        }
      }
    }
    trace("semblance3: count0="+count0+" count1="+count1);
  }

  public float[][][] semblance(Direction3 d, EigenTensors3 t, float[][][] f) {
    float[][][] s = shapeOf(f);
    semblance(d,t,f,s);
    return s;
  }

  /**
   * Computes semblances using slopes in a specified window.
   * This is a rather common method for computing local semblance
   * when image features tend to be oriented along the 2nd dimension;
   * that is, when f[i2][i1] varies most with the sample index i1.
   * Windows for both i1 and i2 are rectangular with uniform weights.
   * <p>
   * The only unusual aspect of this implementation is that the slopes 
   * of lines along which semblance is computed are determined by 
   * eigen-decompositions of specified tensors. This method does not 
   * scan over slopes to maximize semblance.
   * <p>
   * Specifically, in an (x1,x2) coordinate system, slopes of lines 
   * x1 = i1+p*(x2-i2) are computed from eigenvectors u of tensors via 
   * p = -u2/u1, where u1, u2, and p may vary with sample indices i1 and 
   * i2 according to the specified tensors. Semblance numerators and 
   * denominators are computed along these linear trajectories within 
   * windows of i2, and then smoothed in windows of i1 before computing 
   * semblance ratios.
   * @param pmax maximum slope such that -pmax &lt;= p &lt;= pmax.
   * @param hw1 half-width of window in 1st dimension.
   * @param hw2 half-width of window in 2nd dimension.
   * @param et eigen-decomposition of tensor field.
   * @param f input image.
   * @return output semblance image.
   */
  public static float[][] semblanceForSlopes(
    double pmax, int hw1, int hw2, EigenTensors2 et, float[][] f) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[] u = new float[2];
    float[] p = new float[n1];
    float[] x1 = new float[n1];
    float[] f1 = new float[n1];
    float[] snum = new float[n1];
    float[] sden = new float[n1];
    float[] tnum = x1;
    float[] tden = f1;
    float[][] s = new float[n2][n1];
    SincInterp si = new SincInterp();
    float sigma2 = sigmaGaussian(hw2);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma2);

    // For all samples i2, ...
    for (int i2=0; i2<n2; ++i2) {

      // Compute slopes and zero semblance numerators and denominators.
      for (int i1=0; i1<n1; ++i1) {
        et.getEigenvectorU(i1,i2,u);
        float u1 = u[0];
        float u2 = u[1];
        double pi = -u2/u1;
        if ( pi>pmax) pi =  pmax;
        if (-pi>pmax) pi = -pmax;
        p[i1] = (float)pi;
        snum[i1] = 0.0f;
        sden[i1] = 0.0f;
      }

      // Window of samples i2 in semblance computations.
      int j2lo = max(   0,i2-hw1);
      int j2hi = min(n2-1,i2+hw1);
      float scale2 = 1.0f/(j2hi-j2lo+1);

      // For all samples i2 in window, ...
      for (int j2=j2lo; j2<=j2hi; ++j2) {

        // Interpolate according to slopes of lines.
        float d2 = j2-i2;
        for (int i1=0; i1<n1; ++i1)
          x1[i1] = i1+p[i1]*d2;
        si.interpolate(n1,1.0,0.0,f[j2],n1,x1,f1);

        // Accumulate semblance numerators and denominators.
        for (int i1=0; i1<n1; ++i1) {
          float fi = f1[i1];
          snum[i1] += fi;
          sden[i1] += fi*fi;
        }
      }

      // Complete (square and scale) semblance numerators.
      for (int i1=0; i1<n1; ++i1)
        snum[i1] = snum[i1]*snum[i1]*scale2;

      // Compute semblance from smoothed numerators and denominators.
      rgf.apply0(snum,tnum);
      rgf.apply0(sden,tden);
      for (int i1=0; i1<n1; ++i1) {
        //s[i2][i1] = tnum[i1]/tden[i1];
        float sni = tnum[i1];
        float sdi = tden[i1];
        if (sdi<=0.0f || sni<0.0f) {
          s[i2][i1] = 0.0f;
        } else if (sdi<sni) {
          s[i2][i1] = 1.0f;
        } else {
          s[i2][i1] = sni/sdi;
        }
      }

      /*
      // Smooth semblance numerators and denominators over i1.
      float ssnum = 0.0f;
      float ssden = 0.0f;
      for (int i1=0; i1<hw2; ++i1) {
        ssnum += snum[i1];
        ssden += sden[i1];
      }
      for (int i1=0,i1m=i1-hw2-1,i1p=i1+hw2; i1<n1; ++i1,++i1m,++i1p) {
        if (i1p<n1) {
          ssnum += snum[i1p];
          ssden += sden[i1p];
        }
        if (i1m>=0) {
          ssnum -= snum[i1m];
          ssden -= sden[i1m];
        }
        s[i2][i1] = ssnum/ssden;
      }
      */
    }
    return s;
  }
  public static float[][][] semblanceForSlopes(
    double pmax, int hw1, int hw2, EigenTensors3 et, float[][][] f) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[] u = new float[3];
    float[] p = new float[n1];
    float[] q = new float[n1];
    float[] x1 = new float[n1];
    float[] f1 = new float[n1];
    float[] snum = new float[n1];
    float[] sden = new float[n1];
    float[] tnum = x1;
    float[] tden = f1;
    float[][][] s = new float[n3][n2][n1];
    SincInterp si = new SincInterp();
    float sigma2 = sigmaGaussian(hw2);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma2);

    // For all samples i2,i3, ...
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {

      // Compute slopes and zero semblance numerators and denominators.
      for (int i1=0; i1<n1; ++i1) {
        et.getEigenvectorU(i1,i2,i3,u);
        float u1 = u[0];
        float u2 = u[1];
        float u3 = u[2];
        double pi = -u2/u1;
        double qi = -u3/u1;
        if ( pi>pmax) pi =  pmax;
        if (-pi>pmax) pi = -pmax;
        if ( qi>pmax) qi =  pmax;
        if (-qi>pmax) qi = -pmax;
        p[i1] = (float)pi;
        q[i1] = (float)qi;
        snum[i1] = 0.0f;
        sden[i1] = 0.0f;
      }

      // Window of samples i2,i3 in semblance computations.
      int j2lo = max(   0,i2-hw1);
      int j2hi = min(n2-1,i2+hw1);
      int j3lo = max(   0,i3-hw1);
      int j3hi = min(n3-1,i3+hw1);

      // For all samples i2,i3 in window, ...
      for (int j3=j3lo; j3<=j3hi; ++j3) {
      for (int j2=j2lo; j2<=j2hi; ++j2) {

        // Interpolate according to slopes of lines.
        float d2 = j2-i2;
        float d3 = j3-i3;
        for (int i1=0; i1<n1; ++i1)
          x1[i1] = i1+p[i1]*d2+q[i1]*d3;
        si.interpolate(n1,1.0,0.0,f[j3][j2],n1,x1,f1);

        // Accumulate semblance numerators and denominators.
        for (int i1=0; i1<n1; ++i1) {
          float fi = f1[i1];
          snum[i1] += fi;
          sden[i1] += fi*fi;
        }
      }}

      // Complete (square and scale) semblance numerators.
      float scale23 = 1.0f/((j2hi-j2lo+1)*(j3hi-j3lo+1));
      for (int i1=0; i1<n1; ++i1)
        snum[i1] = snum[i1]*snum[i1]*scale23;

      // Compute semblance from smoothed numerators and denominators.
      rgf.apply0(snum,tnum);
      rgf.apply0(sden,tden);
      for (int i1=0; i1<n1; ++i1) {
        //s[i3][i2][i1] = tnum[i1]/tden[i1];
        float sni = tnum[i1];
        float sdi = tden[i1];
        if (sdi<=0.0f || sni<0.0f) {
          s[i3][i2][i1] = 0.0f;
        } else if (sdi<sni) {
          s[i3][i2][i1] = 1.0f;
        } else {
          s[i3][i2][i1] = sni/sdi;
        }
      }
    }}
    return s;
  }

  public static float[][][] applyForTensorTraces(
    double sigma1, double sigma2, float[][][] f)
  {
    LocalOrientFilter lof1 = new LocalOrientFilter(sigma1);
    LocalOrientFilter lof2 = new LocalOrientFilter(sigma2);
    EigenTensors3 t1 = lof1.applyForTensors(f);
    EigenTensors3 t2 = lof2.applyForTensors(f);
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] s = new float[n3][n2][n1];
    float[] d1 = new float[6];
    float[] d2 = new float[6];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t1.getTensor(i1,i2,i3,d1);
          t2.getTensor(i1,i2,i3,d2);
          float a11 = d1[0], a12 = d1[1], a13 = d1[2],
                             a22 = d1[3], a23 = d1[4],
                                          a33 = d1[5];
          float b11 = d2[0], b12 = d2[1], b13 = d2[2],
                             b22 = d2[3], b23 = d2[4],
                                          b33 = d2[5];
          float ta = a11+a22+a33;
          float tb = b11+b22+b33;
          float tab = a11*b11+a12*b12+a13*b13+
                      a12*b12+a22*b22+a23*b23+
                      a13*b13+a23*b23+a33*b33;
          float tatb = ta*tb;
          s[i3][i2][i1] = (tatb>0.0f)?tab/tatb:0.0f;
        }
      }
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static class CopySmoother implements Smoother {
    public void apply(float[] f, float[] g) {
      copy(f,g);
    }
    public void apply(
      Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
    {
      copy(f,g);
    }
    public void apply(
      Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
    {
      copy(f,g);
    }
  }

  private static class LaplacianSmoother implements Smoother {
    LaplacianSmoother(int halfWidth) {
      //_lsf = new LocalSmoothingFilter(0.01,1000);
      _lsf = new LocalSmoothingFilter(0.01,1000,
        new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D71));
      _scale = scaleLaplacian(halfWidth);
    }
    public void apply(float[] f, float[] g) {
      _lsf.apply(_scale,f,g);
    }
    public void apply(
      Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
    {
      setEigenvalues(d,t);
      _lsf.apply(t,_scale,f,g);
    }
    public void apply(
      Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
    {
      setEigenvalues(d,t);
      _lsf.apply(t,_scale,f,g);
    }
    private float _scale;
    private LocalSmoothingFilter _lsf;
  }

  private static class FirSmoother implements Smoother {
    FirSmoother(Smoothing smoothing, int halfWidth) {
      if (smoothing==Smoothing.BOXCAR) {
        _w1 = makeBoxcar1(halfWidth);
        _w2 = makeBoxcar2(halfWidth);
        _w3 = makeBoxcar3(halfWidth);
      } else if (smoothing==Smoothing.GAUSSIAN) {
        _w1 = makeGaussian1(halfWidth);
        _w2 = makeGaussian2(halfWidth);
        _w3 = makeGaussian3(halfWidth);
      }
    }
    public void apply(float[] f, float[] g) {
      applyAll(f,g);
    }
    public void apply(
      Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
    {
      if (d==Direction2.U || d==Direction2.V)
        applyLinear(d,t,f,g);
      else
        applyAll(f,g);
    }
    public void applyAll(float[] f, float[] g) {
      int n1 = f.length;
      float[] w = _w1;
      int n = w.length;
      int m = (n-1)/2;
      for (int i1=0; i1<n1; ++i1) {
        float gi = 0.0f;
        for (int k1=0; k1<n; ++k1) {
          int j1 = i1+k1-m; if (j1<0) j1 = 0; if (j1>=n1) j1 = n1-1;
          gi += w[k1]*f[j1];
        }
        g[i1] = gi;
      }
    }
    public void applyLinear(
      Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
    {
      int n1 = f[0].length;
      int n2 = f.length;
      float[] w = _w1;
      int n = w.length;
      int m = (n-1)/2;
      float[] r = new float[2];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (d==Direction2.U)
            t.getEigenvectorU(i1,i2,r);
          else
            t.getEigenvectorV(i1,i2,r);
          float r1 = r[0], r2 = r[1];
          float gi = 0.0f;
          for (int i=0,ir=-m; i<n; ++i,++ir) {
            float ir1 = ir*r1, ir2 = ir*r2;
            float x1 = i1+ir1, x2 = i2+ir2;
            gi += w[i]*interpolateBilinear(x1,x2,n1,n2,f);
          }
          g[i2][i1] = gi;
        }
      }
    }
    public void applyLinearSinc(
      Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
    {
      int n1 = f[0].length;
      int n2 = f.length;
      float[] w = _w1;
      int n = w.length;
      int m = (n-1)/2;
      float[] r = new float[2];
      SincInterp si = new SincInterp();
      si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (d==Direction2.U)
            t.getEigenvectorU(i1,i2,r);
          else
            t.getEigenvectorV(i1,i2,r);
          float r1 = r[0], r2 = r[1];
          float gi = 0.0f;
          for (int i=0,ir=-m; i<n; ++i,++ir) {
            float ir1 = ir*r1, ir2 = ir*r2;
            float x1 = i1+ir1, x2 = i2+ir2;
            gi += w[i]*si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1,x2);
          }
          g[i2][i1] = gi;
        }
      }
    }
    public void applyAll(float[][] f, float[][] g) {
      int n1 = f[0].length;
      int n2 = f.length;
      float[][] w = _w2;
      int n = w.length;
      int m = (n-1)/2;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float gi = 0.0f;
          for (int k2=0; k2<n; ++k2) {
            int j2 = i2+k2-m; if (j2<0) j2 = 0; if (j2>=n2) j2 = n2-1;
            for (int k1=0; k1<n; ++k1) {
              int j1 = i1+k1-m; if (j1<0) j1 = 0; if (j1>=n1) j1 = n1-1;
              gi += w[k2][k1]*f[j2][j1];
            }
          }
          g[i2][i1] = gi;
        }
      }
    }
    public void apply(
      Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
    {
      if (d==Direction3.U || d==Direction3.V || d==Direction3.W)
        applyLinear(d,t,f,g);
      else if (d==Direction3.UV || d==Direction3.UW || d==Direction3.VW)
        applyPlanar(d,t,f,g);
      else
        applyAll(f,g);
    }
    public void applyLinear(
      Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
    {
      int n1 = f[0][0].length;
      int n2 = f[0].length;
      int n3 = f.length;
      float[] w = _w1;
      int n = w.length;
      int m = (n-1)/2;
      float[] r = new float[3];
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            if (d==Direction3.U)
              t.getEigenvectorU(i1,i2,i3,r);
            else if (d==Direction3.V)
              t.getEigenvectorV(i1,i2,i3,r);
            else 
              t.getEigenvectorW(i1,i2,i3,r);
            float r1 = r[0], r2 = r[1], r3 = r[2];
            float gi = 0.0f;
            for (int i=0,ir=-m; i<n; ++i,++ir) {
              float ir1 = ir*r1, ir2 = ir*r2, ir3 = ir*r3;
              float x1 = i1+ir1, x2 = i2+ir2, x3 = i3+ir3;
              gi += w[i]*interpolateTrilinear(x1,x2,x3,n1,n2,n3,f);
            }
            g[i3][i2][i1] = gi;
          }
        }
      }
    }
    public void applyPlanar(
      Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
    {
      int n1 = f[0][0].length;
      int n2 = f[0].length;
      int n3 = f.length;
      float[][] w = _w2;
      int n = w.length;
      int m = (n-1)/2;
      float[] r = new float[3];
      float[] s = new float[3];
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            if (d==Direction3.UV) {
              t.getEigenvectorU(i1,i2,i3,r);
              t.getEigenvectorV(i1,i2,i3,s);
            } else if (d==Direction3.UW) {
              t.getEigenvectorU(i1,i2,i3,r);
              t.getEigenvectorW(i1,i2,i3,s);
            } else {
              t.getEigenvectorV(i1,i2,i3,r);
              t.getEigenvectorW(i1,i2,i3,s);
            }
            float r1 = r[0], r2 = r[1], r3 = r[2];
            float s1 = s[0], s2 = s[1], s3 = s[2];
            float gi = 0.0f;
            for (int i=0,ir=-m; i<n; ++i,++ir) {
              float ir1 = ir*r1, ir2 = ir*r2, ir3 = ir*r3;
              float x1 = i1+ir1, x2 = i2+ir2, x3 = i3+ir3;
              for (int j=0,js=-m; j<n; ++j,++js) {
                float js1 = js*s1, js2 = js*s2, js3 = js*s3;
                float y1 = x1+js1, y2 = x2+js2, y3 = x3+js3;
                gi += w[i][j]*interpolateTrilinear(y1,y2,y3,n1,n2,n3,f);
              }
            }
            g[i3][i2][i1] = gi;
          }
        }
      }
    }
    public void applyAll(float[][][] f, float[][][] g) {
      int n1 = f[0][0].length;
      int n2 = f[0].length;
      int n3 = f.length;
      float[][][] w = _w3;
      int n = w.length;
      int m = (n-1)/2;
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float gi = 0.0f;
            for (int k3=0; k3<n; ++k3) {
              int j3 = i3+k3-m; if (j3<0) j3 = 0; if (j3>=n3) j3 = n3-1;
              for (int k2=0; k2<n; ++k2) {
                int j2 = i2+k2-m; if (j2<0) j2 = 0; if (j2>=n2) j2 = n2-1;
                for (int k1=0; k1<n; ++k1) {
                  int j1 = i1+k1-m; if (j1<0) j1 = 0; if (j1>=n1) j1 = n1-1;
                  gi += w[k3][k2][k1]*f[j3][j2][j1];
                }
              }
            }
            g[i3][i2][i1] = gi;
          }
        }
      }
    }
    private float[] _w1;
    private float[][] _w2;
    private float[][][] _w3;
  }

  private static boolean allZero(float[] f) {
    int n1 = f.length;
    for (int i1=0; i1<n1; ++i1)
      if (f[i1]!=0.0f)
        return false;
    return true;
  }

  private static Smoother makeSmoothingFilter(Smoothing s, int hw) {
    if (s==Smoothing.BOXCAR || s==Smoothing.GAUSSIAN)
      return new FirSmoother(s,hw);
    else if (s==Smoothing.LAPLACIAN)
      return new LaplacianSmoother(hw);
    else
      return new CopySmoother();
  }

  private static void setEigenvalues(Direction2 d, EigenTensors2 t) {
    float au = 0.0f;
    float av = 0.0f;
    if (d==Direction2.U || d==Direction2.UV)
      au = 1.0f;
    if (d==Direction2.V || d==Direction2.UV)
      av = 1.0f;
    t.setEigenvalues(au,av);
  }
  private static void setEigenvalues(Direction3 d, EigenTensors3 t) {
    float au = 0.0f;
    float av = 0.0f;
    float aw = 0.0f;
    if (d==Direction3.U || 
        d==Direction3.UV || 
        d==Direction3.UW ||
        d==Direction3.UVW)
      au = 1.0f;
    if (d==Direction3.V || 
        d==Direction3.UV || 
        d==Direction3.VW ||
        d==Direction3.UVW)
      av = 1.0f;
    if (d==Direction3.W || 
        d==Direction3.UW || 
        d==Direction3.VW ||
        d==Direction3.UVW)
      aw = 1.0f;
    t.setEigenvalues(au,av,aw);
  }

  private Smoother _smoother1,_smoother2;

  private static float[] shapeOf(float[] f) {
    return new float[f.length];
  }
  private static float[][] shapeOf(float[][] f) {
    return new float[f.length][f[0].length];
  }
  private static float[][][] shapeOf(float[][][] f) {
    return new float[f.length][f[0].length][f[0][0].length];
  }

  private static float sigmaGaussian(int halfWidth) {
    return sqrt(halfWidth*(halfWidth+1)/3.0f);
  }
  private static float scaleLaplacian(int halfWidth) {
    return halfWidth*(halfWidth+1)/6.0f;
  }

  private static Direction2 orthogonal(Direction2 d) {
    if (d==Direction2.U)
      return Direction2.V;
    else
      return Direction2.U;
  }
  private static Direction3 orthogonal(Direction3 d) {
    if (d==Direction3.U)
      return Direction3.VW;
    else if (d==Direction3.V)
      return Direction3.UW;
    else if (d==Direction3.W)
      return Direction3.UV;
    else if (d==Direction3.UV)
      return Direction3.W;
    else if (d==Direction3.UW)
      return Direction3.V;
    else
      return Direction3.U;
  }
  private static float[] normalize(float[] w) {
    return mul(w,1.0f/sum(w));
  }
  private static float[][] normalize(float[][] w) {
    return mul(w,1.0f/sum(w));
  }
  private static float[][][] normalize(float[][][] w) {
    return mul(w,1.0f/sum(w));
  }
  private static float[] makeBoxcar1(int halfWidth) {
    return normalize(fillfloat(1.0f,2*halfWidth+1));
  }
  private static float[][] makeBoxcar2(int halfWidth) {
    int m = halfWidth;
    int n = 2*m+1;
    float[][] w = new float[n][n];
    float mm = m*m;
    for (int i=0,ii=-m; i<n; ++i,++ii)
      for (int j=0,jj=-m; j<n; ++j,++jj)
        if (ii*ii+jj*jj<=mm)
          w[i][j] = 1.0f;
    return normalize(w);
  }
  private static float[][][] makeBoxcar3(int halfWidth) {
    int m = halfWidth;
    int n = 2*m+1;
    float[][][] w = new float[n][n][n];
    float mm = m*m;
    for (int i=0,ii=-m; i<n; ++i,++ii)
      for (int j=0,jj=-m; j<n; ++j,++jj)
        for (int k=0,kk=-m; k<n; ++k,++kk)
          if (ii*ii+jj*jj+kk*kk<=mm)
            w[i][j][k] = 1.0f;
    return normalize(w);
  }
  private static float[] makeGaussian1(int halfWidth) {
    float sigma = sigmaGaussian(halfWidth);
    float s = 0.5f/(sigma*sigma);
    int m = (int)(1.0f+3.0f*sigma);
    int n = m+1+m;
    float[] w = new float[n];
    for (int i=0,ii=-m; i<n; ++i,++ii)
      w[i] = exp(-s*(ii*ii));
    return normalize(w);
  }
  private static float[][] makeGaussian2(int halfWidth) {
    float sigma = sigmaGaussian(halfWidth);
    float s = 0.5f/(sigma*sigma);
    int m = (int)(1.0f+3.0f*sigma);
    int n = m+1+m;
    float[][] w = new float[n][n];
    for (int i=0,ii=-m; i<n; ++i,++ii)
      for (int j=0,jj=-m; j<n; ++j,++jj)
        w[i][j] = exp(-s*(ii*ii+jj*jj));
    return normalize(w);
  }
  private static float[][][] makeGaussian3(int halfWidth) {
    float sigma = sigmaGaussian(halfWidth);
    float s = 0.5f/(sigma*sigma);
    int m = (int)(1.0f+3.0f*sigma);
    int n = m+1+m;
    float[][][] w = new float[n][n][n];
    for (int i=0,ii=-m; i<n; ++i,++ii)
      for (int j=0,jj=-m; j<n; ++j,++jj)
        for (int k=0,kk=-m; k<n; ++k,++kk)
          w[i][j][k] = exp(-s*(ii*ii+jj*jj+kk*kk));
    return normalize(w);
  }

  private static float interpolateBilinear(
    float x1, float x2, int n1, int n2, float[][] f) 
  {
    int j1 = (int)x1;
    int j2 = (int)x2;
    int k1 = j1+1;
    int k2 = j2+1;
    float b1 = x1-j1;
    float b2 = x2-j2;
    if (x1<0.0f) { j1 = 0; k1 = 1; b1 = 0.0f; }
    if (x2<0.0f) { j2 = 0; k2 = 1; b2 = 0.0f; }
    if (k1>=n1) { j1 = n1-2; k1 = n1-1; b1 = 1.0f; }
    if (k2>=n2) { j2 = n2-2; k2 = n2-1; b2 = 1.0f; }
    float a1 = 1.0f-b1;
    float a2 = 1.0f-b2;
    return a2*(a1*f[j2][j1]+b1*f[j2][k1]) +
           b2*(a1*f[k2][j1]+b1*f[k2][k1]);
  }

  private static float interpolateTrilinear(
    float x1, float x2, float x3, int n1, int n2, int n3, float[][][] f) 
  {
    int j1 = (int)x1;
    int j2 = (int)x2;
    int j3 = (int)x3;
    int k1 = j1+1;
    int k2 = j2+1;
    int k3 = j3+1;
    float b1 = x1-j1;
    float b2 = x2-j2;
    float b3 = x3-j3;
    if (x1<0.0f) { j1 = 0; k1 = 1; b1 = 0.0f; }
    if (x2<0.0f) { j2 = 0; k2 = 1; b2 = 0.0f; }
    if (x3<0.0f) { j3 = 0; k3 = 1; b3 = 0.0f; }
    if (k1>=n1) { j1 = n1-2; k1 = n1-1; b1 = 1.0f; }
    if (k2>=n2) { j2 = n2-2; k2 = n2-1; b2 = 1.0f; }
    if (k3>=n3) { j3 = n3-2; k3 = n3-1; b3 = 1.0f; }
    float a1 = 1.0f-b1;
    float a2 = 1.0f-b2;
    float a3 = 1.0f-b3;
    return a3*(a2*(a1*f[j3][j2][j1]+b1*f[j3][j2][k1]) +
               b2*(a1*f[j3][k2][j1]+b1*f[j3][k2][k1])) +
           b3*(a2*(a1*f[k3][j2][j1]+b1*f[k3][j2][k1]) +
               b2*(a1*f[k3][k2][j1]+b1*f[k3][k2][k1]));
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }
} 
