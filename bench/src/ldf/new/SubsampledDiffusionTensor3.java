/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Subsampled diffusion tensors for 3-D images. Tensors for each sample of
 * an image can consume significantly more memory than the image itself.
 * This class subsamples tensors by a factor of two in each image dimension.
 * Eight tensors are averaged to obtain each subsampled tensor. The resulting 
 * data structure requires less memory than the image.
 * <p>
 * Methods that return tensor elements, eigenvectors or eigenvalues return
 * values corresponding to the nearest average tensor. In other words, these 
 * return the same average tensor for up to eight different 3-D sample 
 * indices. These methods do not interpolate tensors.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.26
 */
public class SubsampledDiffusionTensor3 implements DiffusionTensor3 {

  public SubsampledDiffusionTensor3(double sigma, float[][][] x) {
    LocalOrientFilter lof = new LocalOrientFilter(sigma);
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] eu = new float[n3][n2][n1];
    float[][][] ev = new float[n3][n2][n1];
    float[][][] ew = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    lof.apply(x,null,null,
      null,u2,u3,
      null,null,null,
      null,w2,w3,
      eu,ev,ew,
      null,null);
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _eu = subsample(eu);
    _ev = subsample(ev);
    _ew = subsample(ew);
    _u2 = subsample(u2);
    _u3 = subsample(u3);
    _w2 = subsample(w2);
    _w3 = subsample(w3);
  }

  public SubsampledDiffusionTensor3(
    float[][][] d11, float[][][] d12, float[][][] d13,
    float[][][] d22, float[][][] d23, float[][][] d33)
  {
  }

  /**
   * Gets array of elements {D11, D12, D13, D22, D23, D33}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getElements(int i1, int i2, int i3) {
    float[] u = getEigenvectorU(i1,i2,i3);
    float[] v = getEigenvectorV(i1,i2,i3);
    float[] w = getEigenvectorW(i1,i2,i3);
    float[] e = getEigenvalues(i1,i2,i3);
    float u1 = u[0], u2 = u[1], u3 = u[2];
    float v1 = v[0], v2 = v[1], v3 = v[2];
    float w1 = w[0], w2 = w[1], w3 = w[2];
    float eu = e[0], ev = e[1], ew = e[2];
    float d11 = eu*u1*u1+ev*v1*v1+ew*w1*w1;
    float d12 = eu*u1*u2+ev*v1*v2+ew*w1*w2;
    float d13 = eu*u1*u3+ev*v1*v3+ew*w1*w3;
    float d22 = eu*u2*u2+ev*v2*v2+ew*w2*w2;
    float d23 = eu*u2*u3+ev*v2*v3+ew*w2*w3;
    float d33 = eu*u3*u3+ev*v3*v3+ew*w3*w3;
    return new float[]{d11,d12,d13,d22,d23,d33};
  }

  /**
   * Gets the eigenvector u corresponding to the largest eigenvalue.
   * The array contains the vector components {u1, u2, u3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvectorU(int i1, int i2, int i3) {
    int j1 = i1/2;
    int j2 = i2/2;
    int j3 = i3/2;
    float u2 = _u2[j3][j2][j1];
    float u3 = _u3[j3][j2][j1];
    float us = u2*u2+u3*u3;
    float u1 = (us<1.0f)?sqrt(1.0f-us):0.0f;
    return new float[]{u1,u2,u3};
  }

  /**
   * Gets the eigenvector v corresponding to the second largest eigenvalue.
   * The array contains the vector components {v1, v2, v3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvectorV(int i1, int i2, int i3) {
    int j1 = i1/2;
    int j2 = i2/2;
    int j3 = i3/2;
    float u2 = _u2[j3][j2][j1];
    float u3 = _u3[j3][j2][j1];
    float us = u2*u2+u3*u3;
    float u1 = (us<1.0f)?sqrt(1.0f-us):0.0f;
    float w2 = _w2[j3][j2][j1];
    float w3 = _w3[j3][j2][j1];
    float ws = w2*w2+w3*w3;
    float w1 = (ws<1.0f)?sqrt(1.0f-ws):0.0f;
    float v1 = u2*w3-u3*w2;
    float v2 = u3*w1-u1*w3;
    float v3 = u1*w2-u2*w1;
    return new float[]{v1,v2,v3};
  }

  /**
   * Gets the eigenvector u corresponding to the smallest eigenvalue.
   * The array contains the vector components {w1, w2, w3}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvectorW(int i1, int i2, int i3) {
    int j1 = i1/2;
    int j2 = i2/2;
    int j3 = i3/2;
    float w2 = _w2[j3][j2][j1];
    float w3 = _w3[j3][j2][j1];
    float ws = w2*w2+w3*w3;
    float w1 = (ws<1.0f)?sqrt(1.0f-ws):0.0f;
    return new float[]{w1,w2,w3};
  }

  /**
   * Gets the three eigenvalues {eu, ev, ew}.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   */
  public float[] getEigenvalues(int i1, int i2, int i3) {
    int j1 = i1/2;
    int j2 = i2/2;
    int j3 = i3/2;
    float eu = _eu[j3][j2][j1];
    float ev = _ev[j3][j2][j1];
    float ew = _ew[j3][j2][j1];
    return new float[]{eu,ev,ew};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _m1,_m2,_m3;
  private int _n1,_n2,_n3;
  private float[][][] _eu,_ev,_ew;
  private float[][][] _u2,_u3;
  private float[][][] _w2,_w3;

  /**
   * Downsample from [n3][n2][n1] to [(n3+1)/2][(n2+1)/2][(n1+1)/2] samples.
   */
  private static float[][][] subsample(float[][][] x) {
    // Example for n even:
    // 0   0   0   0         n = 4 before subsampling
    //   x       x           m = 2  after subsampling
    //
    // Example for n odd:
    // 0   0   0   0   0     n = 5 before subsampling
    //   x       x       x   m = 3  after subsampling
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int m1 = (n1+1)/2;
    int m2 = (n2+1)/2;
    int m3 = (n3+1)/2;
    float s8 = 1.0f/8.0f;
    float[][][] y =  new float[m3][m2][m1];
    int j1max = n1-1;
    int j2max = n2-1;
    int j3max = n3-1;
    for (int i3=0,j3=0; i3<m3; ++i3,j3+=2) {
      int j30 = j3;
      int j3p = (j3<j3max)?j3+1:j3max;
      for (int i2=0,j2=0; i2<m2; ++i2,j2+=2) {
        int j20 = j2;
        int j2p = (j2<j2max)?j2+1:j2max;
        float[] x00 = x[j30][j20];
        float[] x0p = x[j30][j2p];
        float[] xp0 = x[j3p][j20];
        float[] xpp = x[j3p][j2p];
        for (int i1=0,j1=0; i1<m1; ++i1,j1+=2) {
          int j10 = j1;
          int j1p = (j1<j1max)?j1+1:j1max;
          y[i3][i2][i1] = s8*(x00[j10]+x00[j1p] +
                              x0p[j10]+x0p[j1p] +
                              xp0[j10]+xp0[j1p] +
                              xpp[j10]+xpp[j1p]);
        }
      }
    }
    return y;
  }

  /**
   * Downsample from [n3][n2][n1] to [(n3+1)/2][(n2+1)/2][(n1+1)/2] samples.
   * The 27-sample stencil used to average samples is
   * <pre><code>
   *   [ 1/64 1/32 1/64 ] [ 1/32 1/16 1/32 ] [ 1/64 1/32 1/64 ]
   *   [ 1/32 1/16 1/32 ] [ 1/16 1/8  1/16 ] [ 1/32 1/16 1/32 ]
   *   [ 1/64 1/32 1/64 ] [ 1/32 1/16 1/32 ] [ 1/64 1/32 1/64 ]
   * </code></pre>
   */
  private static float[][][] xsubsample(float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int m1 = (n1+1)/2;
    int m2 = (n2+1)/2;
    int m3 = (n3+1)/2;
    float s8 = 1.0f/8.0f;
    float s16 = 1.0f/16.0f;
    float s32 = 1.0f/32.0f;
    float s64 = 1.0f/64.0f;
    float[][][] y =  new float[m3][m2][m1];
    int j1min = 0;
    int j2min = 0;
    int j3min = 0;
    int j1max = n1-1;
    int j2max = n2-1;
    int j3max = n3-1;
    for (int i3=0,j3=0; i3<m3; ++i3,j3+=2) {
      int j3m = (j3min<j3)?j3-1:j3min;
      int j30 = j3;
      int j3p = (j3<j3max)?j3+1:j3max;
      for (int i2=0,j2=0; i2<m2; ++i2,j2+=2) {
        int j2m = (j2min<j2)?j2-1:j2min;
        int j20 = j2;
        int j2p = (j2<j2max)?j2+1:j2max;
        float[] xmm = x[j3m][j2m];
        float[] xm0 = x[j3m][j20];
        float[] xmp = x[j3m][j2p];
        float[] x0m = x[j30][j2m];
        float[] x00 = x[j30][j20];
        float[] x0p = x[j30][j2p];
        float[] xpm = x[j3p][j2m];
        float[] xp0 = x[j3p][j20];
        float[] xpp = x[j3p][j2p];
        for (int i1=0,j1=0; i1<m1; ++i1,j1+=2) {
          int j1m = (j1min<j1)?j1-1:j1min;
          int j10 = j1;
          int j1p = (j1<j1max)?j1+1:j1max;
          y[i3][i2][i1] = s8*(x00[j10]) +
                         s16*(x00[j1m]+x00[j1p] +
                              x0m[j10]+x0p[j10] +
                              xm0[j10]+xp0[j10]) +
                         s32*(xm0[j1m]+xm0[j1p]+xmm[j10]+xmp[j10] +
                              x0m[j1m]+x0m[j1p]+x0p[j1m]+x0p[j1p] +
                              xp0[j1m]+xp0[j1p]+xpm[j10]+xpp[j10]) +
                         s64*(xmm[j1m]+xmm[j1p]+xmp[j1m]+xmp[j1p] +
                              xpm[j1m]+xpm[j1p]+xpp[j1m]+xpp[j1p]);
        }
      }
    }
    return y;
  }
}
