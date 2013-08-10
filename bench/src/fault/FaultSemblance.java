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

/**
 * Computes semblances that can be used to find faults.
 * Uses slopes of locally linear or planar features to align such
 * features before computing semblance numerators and denominators.
 * After alignment, each numerator is a squared average of image 
 * values, and each denominator is an average of squared values.
 * <p>
 * Before computing semblance, semblance numerators and denominators 
 * should be smoothed within fault lines or planes. This smoothing
 * must be performed by methods not in this class.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.08.13
 */
public class FaultSemblance {

  /**
   * Constructs a fault semblance computer.
   */
  public FaultSemblance() {
    _si = new SincInterp();
    _si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
  }

  /**
   * Returns slopes of locally linear features in the specified image.
   * @param f the image.
   * @return array of slopes.
   */
  public float[][] slopes(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] p = new float[n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(SIGMA1,SLOPE_MAX);
    lsf.findSlopes(f,p);
    return p;
  }

  /**
   * Returns slopes of locally planar features in the specified image.
   * @param f the image.
   * @return array {p2,p3} of slopes.
   */
  public float[][][][] slopes(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(SIGMA1,SLOPE_MAX);
    lsf.findSlopes(f,p2,p3,null);
    return new float[][][][]{p2,p3};
  }

  /**
   * Returns semblance numerators and denominators for the specified image.
   * @param p input array of slopes.
   * @param f input image.
   * @return array {snum,sden} of semblance numerators and denominators.
   */
  public float[][][] semblanceNumDen(float[][] p, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[] xm = new float[n1];
    float[] xp = new float[n1];
    float[] fm = new float[n1];
    float[] fp = new float[n1];
    float[][] sn = new float[n2][n1];
    float[][] sd = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float[] f2 = f[i2];
      float[] f2m = f[i2m];
      float[] f2p = f[i2p];
      float[] p2m = p[i2m];
      float[] p2p = p[i2p];
      float[] sn2 = sn[i2];
      float[] sd2 = sd[i2];
      for (int i1=0; i1<n1; ++i1) {
        xm[i1] = i1-p2m[i1];
        xp[i1] = i1+p2p[i1];
      }
      _si.interpolate(n1,1.0,0.0,f2m,n1,xm,fm);
      _si.interpolate(n1,1.0,0.0,f2p,n1,xp,fp);
      float[] gm = fm, g0 = f2, gp = fp;
      if (i2m==i2) gm = g0;
      if (i2p==i2) gp = g0;
      for (int i1=0; i1<n1; ++i1) {
        float gmi = gm[i1];
        float g0i = g0[i1];
        float gpi = gp[i1];
        float sumn = gmi+g0i+gpi;
        float sumd = gmi*gmi+g0i*g0i+gpi*gpi;
        sn2[i1] = sumn*sumn;
        sd2[i1] = 3.0f*sumd;
      }
    }
    return new float[][][]{sn,sd};
  }

  /**
   * Returns semblance numerators and denominators.
   * Each numerator is a squared average of image values, and each
   * denominator is an average of squared values. Specified slopes 
   * are used to align image samples before this averaging.
   * @param p2 array of slopes in 2nd dimension.
   * @param p3 array of slopes in 3rd dimension.
   * @param f input image.
   * @return array {snum,sden} of semblance numerators and denominators.
   */
  public float[][][][] semblanceNumDen(
    final float[][][] p2, final float[][][] p3, final float[][][] f) 
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] sn = new float[n3][n2][n1];
    final float[][][] sd = new float[n3][n2][n1];
    final SincInterp si = new SincInterp();
    si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      float[] xmm = new float[n1];
      float[] xm0 = new float[n1];
      float[] xmp = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xpm = new float[n1];
      float[] xp0 = new float[n1];
      float[] xpp = new float[n1];
      float[] gmm = new float[n1];
      float[] gm0 = new float[n1];
      float[] gmp = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gpm = new float[n1];
      float[] gp0 = new float[n1];
      float[] gpp = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p2mm = p2[i3m][i2m];
        float[] p2m0 = p2[i3m][i2 ];
        float[] p2mp = p2[i3m][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p200 = p2[i3 ][i2 ];
        float[] p20p = p2[i3 ][i2p];
        float[] p2pm = p2[i3p][i2m];
        float[] p2p0 = p2[i3p][i2 ];
        float[] p2pp = p2[i3p][i2p];
        float[] p3mm = p3[i3m][i2m];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3mp = p3[i3m][i2p];
        float[] p30m = p3[i3 ][i2m];
        float[] p300 = p3[i3 ][i2 ];
        float[] p30p = p3[i3 ][i2p];
        float[] p3pm = p3[i3p][i2m];
        float[] p3p0 = p3[i3p][i2 ];
        float[] p3pp = p3[i3p][i2p];
        float[] sn32 = sn[i3][i2];
        float[] sd32 = sd[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xmm[i1] = i1-p3mm[i1]-p2mm[i1];
          xm0[i1] = i1-p3m0[i1]         ;
          xmp[i1] = i1-p3mp[i1]+p2mp[i1];
          x0m[i1] = i1         -p20m[i1];
          x0p[i1] = i1         +p20p[i1];
          xpm[i1] = i1+p3pm[i1]-p2pm[i1];
          xp0[i1] = i1+p3p0[i1]         ;
          xpp[i1] = i1+p3pp[i1]+p2pp[i1];
        }
        si.interpolate(n1,1.0,0.0,fmm,n1,xmm,gmm);
        si.interpolate(n1,1.0,0.0,fm0,n1,xm0,gm0);
        si.interpolate(n1,1.0,0.0,fmp,n1,xmp,gmp);
        si.interpolate(n1,1.0,0.0,f0m,n1,x0m,g0m);
        si.interpolate(n1,1.0,0.0,f0p,n1,x0p,g0p);
        si.interpolate(n1,1.0,0.0,fpm,n1,xpm,gpm);
        si.interpolate(n1,1.0,0.0,fp0,n1,xp0,gp0);
        si.interpolate(n1,1.0,0.0,fpp,n1,xpp,gpp);
        float[] hmm = gmm, hm0 = gm0, hmp = gmp;
        float[] h0m = g0m, h00 = f00, h0p = g0p;
        float[] hpm = gpm, hp0 = gp0, hpp = gpp;
        if (            i2==0   ) h0m = h00;
        if (            i2==n2-1) h0p = h00;
        if (i3==0               ) hm0 = h00;
        if (i3==n3-1            ) hp0 = h00;
        if (i3==0    && i2==0   ) hmm = h00;
        if (i3==0    && i2==n2-1) hmp = h00;
        if (i3==n3-1 && i2==0   ) hpm = h00;
        if (i3==n3-1 && i2==n2-1) hpp = h00;
        for (int i1=0; i1<n1; ++i1) {
          float hmmi = hmm[i1];
          float hm0i = hm0[i1];
          float hmpi = hmp[i1];
          float h0mi = h0m[i1];
          float h00i = h00[i1];
          float h0pi = h0p[i1];
          float hpmi = hpm[i1];
          float hp0i = hp0[i1];
          float hppi = hpp[i1];
          float sumn = hmmi+hm0i+hmpi+
                       h0mi+h00i+h0pi+
                       hpmi+hp0i+hppi;
          float sumd = hmmi*hmmi+hm0i*hm0i+hmpi*hmpi+
                       h0mi*h0mi+h00i*h00i+h0pi*h0pi+
                       hpmi*hpmi+hp0i*hp0i+hppi*hppi;
          sn32[i1] = sumn*sumn;
          sd32[i1] = 9.0f*sumd;
        }
      }
    }});
    return new float[][][][]{sn,sd};
  }

  /**
   * Returns semblance numerators and denominators.
   * Each numerator is a squared average of image values, and each
   * denominator is an average of squared values. Specified slopes 
   * are used to align image samples before this averaging.
   * @param p array {p2,p3} of slopes in 2nd and 3rd dimensions.
   * @param f input image.
   * @return array {snum,sden} of semblance numerators and denominators.
   */
  public float[][][][] semblanceNumDen(
    final float[][][][] p, final float[][][] f) 
  {
    return semblanceNumDen(p[0],p[1],f);
  }

  /**
   * Returns semblance numerators and denominators for the specified image.
   * First, slopes of locally linear image features are computed. These 
   * slopes are then used to align image samples when computing semblance 
   * numerators and denominators.
   * @param f input image.
   * @return array {snum,sden} of semblance numerators and denominators.
   */
  public float[][][] semblanceNumDen(float[][] f) {
    float[][] p = slopes(f);
    return semblanceNumDen(p,f);
  }

  /**
   * Returns semblance numerators and denominators for the specified image.
   * First, slopes of locally planar image features are computed. These 
   * slopes are then used to align image samples when computing semblance 
   * numerators and denominators.
   * @param f input image.
   * @return array {snum,sden} of semblance numerators and denominators.
   */
  public float[][][][] semblanceNumDen(float[][][] f) {
    float[][][][] p = slopes(f);
    return semblanceNumDen(p,f);
  }

  /**
   * Returns semblance computed for specified numerators and denominators.
   * @param sn array of semblance numerators.
   * @param sd array of semblance denominators.
   * @return array of semblances.
   */
  public float[][] semblanceFromNumDen(float[][] sn, float[][] sd) {
    int n1 = sn[0].length;
    int n2 = sn.length;
    float[][] sr = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] sn2 = sn[i2];
      float[] sd2 = sd[i2];
      float[] sr2 = sr[i2];
      for (int i1=0; i1<n1; ++i1) {
        float sni = sn2[i1];
        float sdi = sd2[i1];
        if (sdi<=0.0f || sni<0.0f) {
          sr2[i1] = 0.0f;
        } else if (sdi<sni) {
          sr2[i1] = 1.0f;
        } else {
          sr2[i1] = sni/sdi;
        }
      }
    }
    return sr;
  }

  /**
   * Returns semblance computed for specified numerators and denominators.
   * @param sn array of semblance numerators.
   * @param sd array of semblance denominators.
   * @return array of semblances.
   */
  public float[][][] semblanceFromNumDen(float[][][] sn, float[][][] sd) {
    int n3 = sn.length;
    final float[][][] fsn = sn;
    final float[][][] fsd = sd;
    final float[][][] sr = new float[n3][][];
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      sr[i3] = semblanceFromNumDen(fsn[i3],fsd[i3]);
    }});
    return sr;
  }

  /**
   * Returns semblance computed for specified numerators and denominators.
   * @param snd array {snum,sden} of semblance numerators and denomiators.
   * @return array of semblances.
   */
  public float[][] semblanceFromNumDen(float[][][] snd) {
    return semblanceFromNumDen(snd[0],snd[1]);
  }

  /**
   * Returns semblance computed for specified numerators and denominators.
   * @param snd array {snum,sden} of semblance numerators and denomiators.
   * @return array of semblances.
   */
  public float[][][] semblanceFromNumDen(float[][][][] snd) {
    return semblanceFromNumDen(snd[0],snd[1]);
  }

  /**
   * Returns an image with specified samples taper to zero at edges.
   * Tapering enables simple zero-value boundary conditions to be
   * used in fault detection without artifacts caused by abrupt 
   * truncations of strong image features near edges.
   * @param m width of the tapered band of samples at each edge.
   * @param f input image.
   * @return the tapered image.
   */
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
    /* disable horizontal taper
    for (int i2=0,j2=n2-1; i2<m; ++i2,--j2) {
      float ti = t[i2];
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] *= ti;
        g[j2][i1] *= ti;
      }
    }
    */
    return g;
  }

  /**
   * Returns an image with specified samples taper to zero at edges.
   * Tapering enables simple zero-value boundary conditions to be
   * used in fault detection without artifacts caused by abrupt 
   * truncations of strong image features near edges.
   * @param m width of the tapered band of samples at each edge.
   * @param f input image.
   * @return the tapered image.
   */
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

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static float SIGMA1 = 8.0f;
  private static float SLOPE_MAX = 5.0f;

  private SincInterp _si;
}
