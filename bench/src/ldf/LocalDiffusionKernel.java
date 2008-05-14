/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * A local diffusion kernel for use in anisotropic diffusion filtering.
 * <p>
 * This kernel is a filter that computes y = y+G'DGx where G is the 
 * gradient operator, G' is its adjoint, and D is a local diffusion 
 * tensor that determines for each image sample the filter coefficients.
 * <p>
 * A local diffusion kernel is typically used in combinations with others.
 * For example, the filter implied by (I+G'DG)y = G'DGx acts as a notch
 * filter. It attenuates features for which G'DG is zero while preserving 
 * other features. The diffusion tensors D control the width, orientation,
 * and anisotropy of the spectral notch. Note that application of this filter 
 * requires solving a sparse symmetric positive-definite system of equations.
 * <p>
 * An even simpler example is the filter implied by (I+G'DG)y = x. This
 * filter smooths features in the directions implied by the tensors D.
 * Again, application of this filter requires solving a sparse symmetric 
 * positive-definite system of equations.
 * <p>
 * The accumulation of the kernel output in y = y+G'DGx is useful when
 * constructing such combination filters. Given y = 0, this kernel
 * computes y = G'DGx. Given y = x, it computes y = (I+G'DG)x.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.10
 */
public class LocalDiffusionKernel {

  /**
   * Constructs a local diffusion kernel.
   */
  public LocalDiffusionKernel() {
    this(1.0/12.0);
  }

  /**
   * Constructs a local diffusion kernel with an experimental factor.
   * @param rs the experimental factor for tuning gradient approximations.
   */
  public LocalDiffusionKernel(double rs) {
    _rs = (float)rs;
    _ers = _rs;
    _frs = 1.0f-2.0f*_rs;
    float t = 5.0f/12.0f-1.0f/sqrt(6.0f);
    //t = 0.25f; // TESTING!
    float r = (1.0f-sqrt(t))*(1.0f-sqrt(t));
    float s = sqrt(r*t);
    System.out.println("r="+r+" s="+s+" t="+t);
    _erst = s*s;
    _frst = s*(r+t);
    _grst = (r+t)*(r+t);
  }

  ///////////////////////////////////////////////////////////////////////////
  // 2-D

  /**
   * Computes y = y+G'DGx, for specified local 2-D diffusion tensors D.
   * @param ldt local diffusion tensors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void apply(LocalDiffusionTensors2 ldt, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[] d = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
        ldt.getTensor(i1p,i2p,d);
        float a = 0.5f*d[0];
        float b = 0.5f*d[1];
        float c = 0.5f*d[2];
        float t = 2.0f*_rs*(a+c);
        //float t = abs(b);
        float xpp = x[i2p][i1p];
        float xpm = x[i2p][i1m];
        float xmp = x[i2m][i1p];
        float xmm = x[i2m][i1m];
        float apppm = (a-t)*(xpp-xpm);
        float ampmm = (a-t)*(xmp-xmm);
        float bppmm = (b+t)*(xpp-xmm);
        float bpmmp = (b-t)*(xpm-xmp);
        float cppmp = (c-t)*(xpp-xmp);
        float cpmmm = (c-t)*(xpm-xmm);
        y[i2p][i1p] += apppm+bppmm+cppmp;
        y[i2p][i1m] -= apppm+bpmmp-cpmmm;
        y[i2m][i1p] += ampmm+bpmmp-cppmp;
        y[i2m][i1m] -= ampmm+bppmm+cpmmm;
        if (ZERO_SLOPE_BOUNDARIES) {
          if (i1m==0) {
            cpmmm = 2.0f*c*(xpm-xmm);
            y[i2m][i1m] -= cpmmm;
            y[i2p][i1m] += cpmmm;
          } else if (i1p==n1m) {
            cppmp = 2.0f*c*(xpp-xmp);
            y[i2m][i1p] -= cppmp;
            y[i2p][i1p] += cppmp;
          }
          if (i2m==0) {
            ampmm = 2.0f*a*(xmp-xmm);
            y[i2m][i1m] -= ampmm;
            y[i2m][i1p] += ampmm;
          } else if (i2p==n2m) {
            apppm = 2.0f*a*(xpp-xpm);
            y[i2p][i1m] -= apppm;
            y[i2p][i1p] += apppm;
          }
        }
      }
    }
  }

  /**
   * Gets filter coefficients for this kernel.
   * @param ldt local diffusion 2x2 tensors.
   * @return arrays[5][n2][n1] of coefficients.
   *  The five arrays are {s00,s0p,spm,sp0,spp}, where s00 corresponds to 
   *  lag[0][0], s0p to lag[0][1], spm to lag[1][-1], and so on. 
   *  (Here, m, 0, and p denote minus, zero, and plus, respectively.)
   * @return the filter.
   */
  public float[][][] getCoefficients(LocalDiffusionTensors2 ldt) {
    int n1 = ldt.getN1();
    int n2 = ldt.getN2();
    int n1m = n1-1;
    int n2m = n2-1;
    float[][][] s = new float[5][n2][n1];
    float[][] s00 = s[0];
    float[][] s0p = s[1];
    float[][] spm = s[2];
    float[][] sp0 = s[3];
    float[][] spp = s[4];
    float[] d = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
        ldt.getTensor(i1p,i2p,d);
        float a = 0.5f*d[0];
        float b = 0.5f*d[1];
        float c = 0.5f*d[2];
        float t = 2.0f*_rs*(a+c);
        //float t = abs(b);
        s00[i2p][i1p] += a+b+c-t;
        s00[i2p][i1m] += a-b+c-t;
        s00[i2m][i1p] += a-b+c-t;
        s00[i2m][i1m] += a+b+c-t;
        s0p[i2p][i1m] -= a-t;
        s0p[i2m][i1m] -= a-t;
        spp[i2m][i1m] -= b+t;
        spm[i2m][i1p] += b-t;
        sp0[i2m][i1p] -= c-t;
        sp0[i2m][i1m] -= c-t;
        if (ZERO_SLOPE_BOUNDARIES) {
          if (i1m==0) {
            s00[i2p][i1m] += 2.0f*c;
            s00[i2m][i1m] += 2.0f*c;
            sp0[i2m][i1m] -= 2.0f*c;
          } else if (i1p==n1m) {
            s00[i2p][i1p] += 2.0f*c;
            s00[i2m][i1p] += 2.0f*c;
            sp0[i2m][i1p] -= 2.0f*c;
          }
          if (i2m==0) {
            s00[i2m][i1p] += 2.0f*a;
            s00[i2m][i1m] += 2.0f*a;
            s0p[i2m][i1m] -= 2.0f*a;
          } else if (i2p==n2m) {
            s00[i2p][i1p] += 2.0f*a;
            s00[i2p][i1m] += 2.0f*a;
            s0p[i2p][i1m] -= 2.0f*a;
          }
        }
      }
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////
  // 3-D

  /**
   * Computes y = y+G'DGx, for specified local 3-D diffusion tensors D.
   * @param ldt local diffusion tensors.
   * @param x input image. Must be distinct from the array y.
   * @param y input/output image. Must be distinct from the array x.
   */
  public void apply(DiffusionTensors3 ldt, float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float e = _erst;
    float f = _frst;
    float g = _grst;
    float fp50 = f*0.50f;
    float gp25 = g*0.25f;
    float[] d = new float[6];
    for (int i3m=0,i3p=1; i3p<n3; ++i3m,++i3p) {
      for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
        for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
          ldt.getTensor(i1m,i2m,i3m,d);
          float d11 = d[0];
          float d12 = d[1];
          float d13 = d[2];
          float d22 = d[3];
          float d23 = d[4];
          float d33 = d[5];
          float pd12 = 0.125f*d12;
          float pd13 = 0.125f*d13;
          float pd23 = 0.125f*d23;
          float ed11 = e*d11;
          float ed22 = e*d22;
          float ed33 = e*d33;
          float eddd = ed11+ed22+ed33;
          float fd11 = fp50*(d22+d33);
          float fd22 = fp50*(d11+d33);
          float fd33 = fp50*(d11+d22);
          float td11 = fd11-ed11;
          float td22 = fd22-ed22;
          float td33 = fd33-ed33;
          float hd11 = gp25*d11-fd11;
          float hd22 = gp25*d22-fd22;
          float hd33 = gp25*d33-fd33;
          float hd12 = td33+pd12;
          float hd21 = td33-pd12;
          float hd13 = td22+pd13;
          float hd31 = td22-pd13;
          float hd23 = td11+pd23;
          float hd32 = td11-pd23;
          float hdpp = eddd+pd12+pd13+pd23;
          float hdpm = eddd-pd12+pd13-pd23;
          float hdmp = eddd+pd12-pd13-pd23;
          float hdmm = eddd-pd12-pd13+pd23;
          float xppp = x[i3p][i2p][i1p];
          float xppm = x[i3p][i2p][i1m];
          float xpmp = x[i3p][i2m][i1p];
          float xpmm = x[i3p][i2m][i1m];
          float xmpp = x[i3m][i2p][i1p];
          float xmpm = x[i3m][i2p][i1m];
          float xmmp = x[i3m][i2m][i1p];
          float xmmm = x[i3m][i2m][i1m];
          float ypppppm = hd11*(xppp-xppm); // hd11
          float ypmppmm = hd11*(xpmp-xpmm);
          float ymppmpm = hd11*(xmpp-xmpm);
          float ymmpmmm = hd11*(xmmp-xmmm);
          float yppppmp = hd22*(xppp-xpmp); // hd22
          float yppmpmm = hd22*(xppm-xpmm);
          float ymppmmp = hd22*(xmpp-xmmp);
          float ympmmmm = hd22*(xmpm-xmmm);
          float ypppmpp = hd33*(xppp-xmpp); // hd33
          float yppmmpm = hd33*(xppm-xmpm);
          float ypmpmmp = hd33*(xpmp-xmmp);
          float ypmmmmm = hd33*(xpmm-xmmm);
          float yppppmm = hd12*(xppp-xpmm); // hd12, hd21
          float ymppmmm = hd12*(xmpp-xmmm);
          float ypmpppm = hd21*(xpmp-xppm);
          float ymmpmpm = hd21*(xmmp-xmpm);
          float ypppmpm = hd13*(xppp-xmpm); // hd13, hd31
          float ypmpmmm = hd13*(xpmp-xmmm);
          float ymppppm = hd31*(xmpp-xppm);
          float ymmppmm = hd31*(xmmp-xpmm);
          float ypppmmp = hd23*(xppp-xmmp); // hd23, hd32
          float yppmmmm = hd23*(xppm-xmmm);
          float ympppmp = hd32*(xmpp-xpmp);
          float ympmpmm = hd32*(xmpm-xpmm);
          float ypppmmm = hdpp*(xppp-xmmm); // hdpp, hdpm, hdmp, hdmm
          float ypmpmpm = hdpm*(xpmp-xmpm);
          float ympppmm = hdmp*(xmpp-xpmm);
          float ymmpppm = hdmm*(xmmp-xppm);
          y[i3p][i2p][i1p] += 
            ypppppm+yppppmp+ypppmpp+yppppmm+ypppmpm+ypppmmp+ypppmmm;
          y[i3p][i2p][i1m] -=
            ypppppm-yppmpmm-yppmmpm+ypmpppm+ymppppm-yppmmmm+ymmpppm;
          y[i3p][i2m][i1p] +=
            ypmppmm-yppppmp+ypmpmmp+ypmpppm+ypmpmmm-ympppmp+ypmpmpm;
          y[i3p][i2m][i1m] -=
            ypmppmm+yppmpmm-ypmmmmm+yppppmm+ymmppmm+ympmpmm+ympppmm;
          y[i3m][i2p][i1p] +=
            ymppmpm+ymppmmp-ypppmpp+ymppmmm+ymppppm+ympppmp+ympppmm;
          y[i3m][i2p][i1m] -=
            ymppmpm-ympmmmm+yppmmpm+ymmpmpm+ypppmpm-ympmpmm+ypmpmpm;
          y[i3m][i2m][i1p] +=
            ymmpmmm-ymppmmp-ypmpmmp+ymmpmpm+ymmppmm-ypppmmp+ymmpppm;
          y[i3m][i2m][i1m] -=
            ymmpmmm+ympmmmm+ypmmmmm+ymppmmm+ypmpmmm+yppmmmm+ypppmmm;
          if (ZERO_SLOPE_BOUNDARIES) {
            if (i1m==0) {
              hd22 = 0.5f*d22;
              hd33 = 0.5f*d33;
              hd23 = 0.5f*d23;
              hd32 = -hd23;
              yppmpmm = hd22*(xppm-xpmm);
              ympmmmm = hd22*(xmpm-xmmm);
              yppmmpm = hd33*(xppm-xmpm);
              ypmmmmm = hd33*(xpmm-xmmm);
              yppmmmm = hd23*(xppm-xmmm);
              ympmpmm = hd32*(xmpm-xpmm);
              y[i3m][i2m][i1m] -= ypmmmmm+ympmmmm+yppmmmm;
              y[i3m][i2p][i1m] -= yppmmpm-ympmmmm-ympmpmm;
              y[i3p][i2m][i1m] += ypmmmmm-yppmpmm-ympmpmm;
              y[i3p][i2p][i1m] += yppmpmm+yppmmpm+yppmmmm;
              if (i2m==0) {
                ypmmmmm = d33*(xpmm-xmmm);
                y[i3m][i2m][i1m] -= ypmmmmm;
                y[i3p][i2m][i1m] += ypmmmmm;
              } else if (i2p==n2m) {
                yppmmpm = d33*(xppm-xmpm);
                y[i3m][i2p][i1m] -= yppmmpm;
                y[i3p][i2p][i1m] += yppmmpm;
              }
              if (i3m==0) {
                ympmmmm = d22*(xmpm-xmmm);
                y[i3m][i2m][i1m] -= ympmmmm;
                y[i3m][i2p][i1m] += ympmmmm;
              } else if (i3p==n3m) {
                yppmpmm = d22*(xppm-xpmm);
                y[i3p][i2m][i1m] -= yppmpmm;
                y[i3p][i2p][i1m] += yppmpmm;
              }
            } else if (i1p==n1m) {
              hd22 = 0.5f*d22;
              hd33 = 0.5f*d33;
              hd23 = 0.5f*d23;
              hd32 = -hd23;
              yppppmp = hd22*(xppp-xpmp);
              ymppmmp = hd22*(xmpp-xmmp);
              ypppmpp = hd33*(xppp-xmpp);
              ypmpmmp = hd33*(xpmp-xmmp);
              ypppmmp = hd23*(xppp-xmmp);
              ympppmp = hd32*(xmpp-xpmp);
              y[i3m][i2m][i1p] -= ypmpmmp+ymppmmp+ypppmmp;
              y[i3m][i2p][i1p] -= ypppmpp-ymppmmp-ympppmp;
              y[i3p][i2m][i1p] += ypmpmmp-yppppmp-ympppmp;
              y[i3p][i2p][i1p] += yppppmp+ypppmpp+ypppmmp;
              if (i2m==0) {
                ypmpmmp = d33*(xpmp-xmmp);
                y[i3m][i2m][i1p] -= ypmpmmp;
                y[i3p][i2m][i1p] += ypmpmmp;
              } else if (i2p==n2m) {
                ypppmpp = d33*(xppp-xmpp);
                y[i3m][i2p][i1p] -= ypppmpp;
                y[i3p][i2p][i1p] += ypppmpp;
              }
              if (i3m==0) {
                ymppmmp = d22*(xmpp-xmmp);
                y[i3m][i2m][i1p] -= ymppmmp;
                y[i3m][i2p][i1p] += ymppmmp;
              } else if (i3p==n3m) {
                yppppmp = d22*(xppp-xpmp);
                y[i3p][i2m][i1p] -= yppppmp;
                y[i3p][i2p][i1p] += yppppmp;
              }
            }
            if (i2m==0) {
              hd11 = 0.5f*d11;
              hd33 = 0.5f*d33;
              hd13 = 0.5f*d13;
              hd31 = -hd13;
              ypmppmm = hd11*(xpmp-xpmm);
              ymmpmmm = hd11*(xmmp-xmmm);
              ypmpmmp = hd33*(xpmp-xmmp);
              ypmmmmm = hd33*(xpmm-xmmm);
              ypmpmmm = hd13*(xpmp-xmmm);
              ymmppmm = hd31*(xmmp-xpmm);
              y[i3m][i2m][i1m] -= ypmmmmm+ymmpmmm+ypmpmmm;
              y[i3m][i2m][i1p] -= ypmpmmp-ymmpmmm-ymmppmm;
              y[i3p][i2m][i1m] += ypmmmmm-ypmppmm-ymmppmm;
              y[i3p][i2m][i1p] += ypmppmm+ypmpmmp+ypmpmmm;
              if (i1m==0) {
                ypmmmmm = d33*(xpmm-xmmm);
                y[i3m][i2m][i1m] -= ypmmmmm;
                y[i3p][i2m][i1m] += ypmmmmm;
              } else if (i1p==n1m) {
                ypmpmmp = d33*(xpmp-xmmp);
                y[i3m][i2m][i1p] -= ypmpmmp;
                y[i3p][i2m][i1p] += ypmpmmp;
              }
              if (i3m==0) {
                ymmpmmm = d11*(xmmp-xmmm);
                y[i3m][i2m][i1m] -= ymmpmmm;
                y[i3m][i2m][i1p] += ymmpmmm;
              } else if (i3p==n3m) {
                ypmppmm = d11*(xpmp-xpmm);
                y[i3p][i2m][i1m] -= ypmppmm;
                y[i3p][i2m][i1p] += ypmppmm;
              }
            } else if (i2p==n2m) {
              hd11 = 0.5f*d11;
              hd33 = 0.5f*d33;
              hd13 = 0.5f*d13;
              hd31 = -hd13;
              ypppppm = hd11*(xppp-xppm);
              ymppmpm = hd11*(xmpp-xmpm);
              ypppmpp = hd33*(xppp-xmpp);
              yppmmpm = hd33*(xppm-xmpm);
              ypppmpm = hd13*(xppp-xmpm);
              ymppppm = hd31*(xmpp-xppm);
              y[i3m][i2p][i1m] -= yppmmpm+ymppmpm+ypppmpm;
              y[i3m][i2p][i1p] -= ypppmpp-ymppmpm-ymppppm;
              y[i3p][i2p][i1m] += yppmmpm-ypppppm-ymppppm;
              y[i3p][i2p][i1p] += ypppppm+ypppmpp+ypppmpm;
              if (i1m==0) {
                yppmmpm = d33*(xppm-xmpm);
                y[i3m][i2p][i1m] -= yppmmpm;
                y[i3p][i2p][i1m] += yppmmpm;
              } else if (i1p==n1m) {
                ypppmpp = d33*(xppp-xmpp);
                y[i3m][i2p][i1p] -= ypppmpp;
                y[i3p][i2p][i1p] += ypppmpp;
              }
              if (i3m==0) {
                ymppmpm = d11*(xmpp-xmpm);
                y[i3m][i2p][i1m] -= ymppmpm;
                y[i3m][i2p][i1p] += ymppmpm;
              } else if (i3p==n3m) {
                ypppppm = d11*(xppp-xppm);
                y[i3p][i2p][i1m] -= ypppppm;
                y[i3p][i2p][i1p] += ypppppm;
              }
            }
            if (i3m==0) {
              hd11 = 0.5f*d11;
              hd22 = 0.5f*d22;
              hd12 = 0.5f*d12;
              hd21 = -hd12;
              ymppmpm = hd11*(xmpp-xmpm);
              ymmpmmm = hd11*(xmmp-xmmm);
              ymppmmp = hd22*(xmpp-xmmp);
              ympmmmm = hd22*(xmpm-xmmm);
              ymppmmm = hd12*(xmpp-xmmm);
              ymmpmpm = hd21*(xmmp-xmpm);
              y[i3m][i2m][i1m] -= ympmmmm+ymmpmmm+ymppmmm;
              y[i3m][i2m][i1p] -= ymppmmp-ymmpmmm-ymmpmpm;
              y[i3m][i2p][i1m] += ympmmmm-ymppmpm-ymmpmpm;
              y[i3m][i2p][i1p] += ymppmpm+ymppmmp+ymppmmm;
              if (i1m==0) {
                ympmmmm = d22*(xmpm-xmmm);
                y[i3m][i2m][i1m] -= ympmmmm;
                y[i3m][i2p][i1m] += ympmmmm;
              } else if (i1p==n1m) {
                ymppmmp = d22*(xmpp-xmmp);
                y[i3m][i2m][i1p] -= ymppmmp;
                y[i3m][i2p][i1p] += ymppmmp;
              }
              if (i2m==0) {
                ymmpmmm = d11*(xmmp-xmmm);
                y[i3m][i2m][i1m] -= ymmpmmm;
                y[i3m][i2m][i1p] += ymmpmmm;
              } else if (i2p==n2m) {
                ymppmpm = d11*(xmpp-xmpm);
                y[i3m][i2p][i1m] -= ymppmpm;
                y[i3m][i2p][i1p] += ymppmpm;
              }
            } else if (i3p==n3m) {
              hd11 = 0.5f*d11;
              hd22 = 0.5f*d22;
              hd12 = 0.5f*d12;
              hd21 = -hd12;
              ypppppm = hd11*(xppp-xppm);
              ypmppmm = hd11*(xpmp-xpmm);
              yppppmp = hd22*(xppp-xpmp);
              yppmpmm = hd22*(xppm-xpmm);
              yppppmm = hd12*(xppp-xpmm);
              ypmpppm = hd21*(xpmp-xppm);
              y[i3p][i2m][i1m] -= yppmpmm+ypmppmm+yppppmm;
              y[i3p][i2m][i1p] -= yppppmp-ypmppmm-ypmpppm;
              y[i3p][i2p][i1m] += yppmpmm-ypppppm-ypmpppm;
              y[i3p][i2p][i1p] += ypppppm+yppppmp+yppppmm;
              if (i1m==0) {
                yppmpmm = d22*(xppm-xpmm);
                y[i3p][i2m][i1m] -= yppmpmm;
                y[i3p][i2p][i1m] += yppmpmm;
              } else if (i1p==n1m) {
                yppppmp = d22*(xppp-xpmp);
                y[i3p][i2m][i1p] -= yppppmp;
                y[i3p][i2p][i1p] += yppppmp;
              }
              if (i2m==0) {
                ypmppmm = d11*(xpmp-xpmm);
                y[i3p][i2m][i1m] -= ypmppmm;
                y[i3p][i2m][i1p] += ypmppmm;
              } else if (i2p==n2m) {
                ypppppm = d11*(xppp-xppm);
                y[i3p][i2p][i1m] -= ypppppm;
                y[i3p][i2p][i1p] += ypppppm;
              }
            }
          }
        }
      }
    }
  }
  public float[][][][] getCoefficients(DiffusionTensors3 ldt) {
    int n1 = ldt.getN1();
    int n2 = ldt.getN2();
    int n3 = ldt.getN3();
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float e = _erst;
    float f = _frst;
    float g = _grst;
    float fp50 = f*0.50f;
    float gp25 = g*0.25f;
    float[][][][] s = new float[14][n3][n2][n1];
    float[][][] s000 = s[0];
    float[][][] s00p = s[1];
    float[][][] s0pm = s[2];
    float[][][] s0p0 = s[3];
    float[][][] s0pp = s[4];
    float[][][] spmm = s[5];
    float[][][] spm0 = s[6];
    float[][][] spmp = s[7];
    float[][][] sp0m = s[8];
    float[][][] sp00 = s[9];
    float[][][] sp0p = s[10];
    float[][][] sppm = s[11];
    float[][][] spp0 = s[12];
    float[][][] sppp = s[13];
    float[] d = new float[6];
    for (int i3m=0,i3p=1; i3p<n3; ++i3m,++i3p) {
      for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
        for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
          ldt.getTensor(i1m,i2m,i3m,d);
          float d11 = d[0];
          float d12 = d[1];
          float d13 = d[2];
          float d22 = d[3];
          float d23 = d[4];
          float d33 = d[5];

          // Coefficients of 2x2x2 stencils exactly as in method apply.
          float pd12 = 0.125f*d12;
          float pd13 = 0.125f*d13;
          float pd23 = 0.125f*d23;
          float ed11 = e*d11;
          float ed22 = e*d22;
          float ed33 = e*d33;
          float eddd = ed11+ed22+ed33;
          float fd11 = fp50*(d22+d33);
          float fd22 = fp50*(d11+d33);
          float fd33 = fp50*(d11+d22);
          float td11 = fd11-ed11;
          float td22 = fd22-ed22;
          float td33 = fd33-ed33;
          float hd11 = gp25*d11-fd11;
          float hd22 = gp25*d22-fd22;
          float hd33 = gp25*d33-fd33;
          float hd12 = td33+pd12;
          float hd21 = td33-pd12;
          float hd13 = td22+pd13;
          float hd31 = td22-pd13;
          float hd23 = td11+pd23;
          float hd32 = td11-pd23;
          float hdpp = eddd+pd12+pd13+pd23;
          float hdpm = eddd-pd12+pd13-pd23;
          float hdmp = eddd+pd12-pd13-pd23;
          float hdmm = eddd-pd12-pd13+pd23;

          // The comment at the end of each line below indicates the element 
          // of the 2x2x2 stencil by which the preceding coefficient(s) is 
          // multiplied in the method apply. For example, the comment // mmp 
          // in the line for s00p[i3m][i2m][i1m] below tells us to find the 
          // coefficient that is multiplied by x[i3m][i2m][i1p] to compute
          // y[i3m][i2m][i1m]. In the method apply, that coefficient is hd11.
          // Note that sxxx[i3m][i2m][i1m] has eight updates, and that
          // sxxx[i3m][i2m][i1p] has seven updates, and so on, so that
          // sxxx[i3p][i2p][i1p] has one update. 
          float hd112233 = hd11+hd22+hd33;
          s000[i3m][i2m][i1m] += hd112233+hd12+hd13+hd23+hdpp; // mmm
          s000[i3m][i2m][i1p] += hd112233+hd21+hd31+hd23+hdmm; // mmp
          s000[i3m][i2p][i1m] += hd112233+hd21+hd13+hd32+hdpm; // mpm
          s000[i3m][i2p][i1p] += hd112233+hd12+hd31+hd32+hdmp; // mpp
          s000[i3p][i2m][i1m] += hd112233+hd12+hd31+hd32+hdmp; // pmm
          s000[i3p][i2m][i1p] += hd112233+hd21+hd13+hd32+hdpm; // pmp
          s000[i3p][i2p][i1m] += hd112233+hd21+hd31+hd23+hdmm; // ppm
          s000[i3p][i2p][i1p] += hd112233+hd12+hd13+hd23+hdpp; // ppp
          s00p[i3m][i2m][i1m] -= hd11; // mmp
          s00p[i3m][i2p][i1m] -= hd11; // mpp
          s00p[i3p][i2p][i1m] -= hd11; // ppp
          s00p[i3p][i2m][i1m] -= hd11; // pmp
          s0p0[i3m][i2m][i1m] -= hd22; // mpm
          s0p0[i3m][i2m][i1p] -= hd22; // mpp
          s0p0[i3p][i2m][i1m] -= hd22; // ppm
          s0p0[i3p][i2m][i1p] -= hd22; // ppp
          sp00[i3m][i2m][i1m] -= hd33; // pmm
          sp00[i3m][i2m][i1p] -= hd33; // pmp
          sp00[i3m][i2p][i1m] -= hd33; // ppm
          sp00[i3m][i2p][i1p] -= hd33; // ppp
          s0pp[i3m][i2m][i1m] -= hd12; // mpp
          s0pp[i3p][i2m][i1m] -= hd12; // ppp
          s0pm[i3m][i2m][i1p] -= hd21; // mpm
          s0pm[i3p][i2m][i1p] -= hd21; // ppm
          sp0p[i3m][i2m][i1m] -= hd13; // pmp
          sp0p[i3m][i2p][i1m] -= hd13; // ppp
          sp0m[i3m][i2m][i1p] -= hd31; // pmm
          sp0m[i3m][i2p][i1p] -= hd31; // ppm
          spp0[i3m][i2m][i1m] -= hd23; // ppm
          spp0[i3m][i2m][i1p] -= hd23; // ppp
          spm0[i3m][i2p][i1m] -= hd32; // pmm
          spm0[i3m][i2p][i1p] -= hd32; // pmp
          sppp[i3m][i2m][i1m] -= hdpp; // ppp
          sppm[i3m][i2m][i1p] -= hdmm; // ppm
          spmp[i3m][i2p][i1m] -= hdpm; // pmp
          spmm[i3m][i2p][i1p] -= hdmp; // pmm

          // Again, boundary conditions must be consistent with method apply.
          if (ZERO_SLOPE_BOUNDARIES) {
            if (i1m==0) {
              hd22 = 0.5f*d22;
              hd33 = 0.5f*d33;
              hd23 = 0.5f*d23;
              hd32 = -hd23;
              s000[i3m][i2m][i1m] += hd22+hd33+hd23; // mmm
              s000[i3m][i2p][i1m] += hd22+hd33+hd32; // mpm
              s000[i3p][i2m][i1m] += hd22+hd33+hd32; // pmm
              s000[i3p][i2p][i1m] += hd22+hd33+hd23; // ppm
              s0p0[i3m][i2m][i1m] -= hd22; // mpm
              s0p0[i3p][i2m][i1m] -= hd22; // ppm
              sp00[i3m][i2m][i1m] -= hd33; // pmm
              sp00[i3m][i2p][i1m] -= hd33; // ppm
              spp0[i3m][i2m][i1m] -= hd23; // ppm
              spm0[i3m][i2p][i1m] -= hd32; // pmm
              if (i2m==0) {
                s000[i3m][i2m][i1m] += d33; // mmm
                s000[i3p][i2m][i1m] += d33; // pmm
                sp00[i3m][i2m][i1m] -= d33; // pmm
              } else if (i2p==n2m) {
                s000[i3m][i2p][i1m] += d33; // mpm
                s000[i3p][i2p][i1m] += d33; // ppm
                sp00[i3m][i2p][i1m] -= d33; // ppm
              }
              if (i3m==0) {
                s000[i3m][i2m][i1m] += d22; // mmm
                s000[i3m][i2p][i1m] += d22; // mpm
                s0p0[i3m][i2m][i1m] -= d22; // mpm
              } else if (i3p==n3m) {
                s000[i3p][i2m][i1m] += d22; // pmm
                s000[i3p][i2p][i1m] += d22; // ppm
                s0p0[i3p][i2m][i1m] -= d22; // ppm
              }
            } else if (i1p==n1m) {
              hd22 = 0.5f*d22;
              hd33 = 0.5f*d33;
              hd23 = 0.5f*d23;
              hd32 = -hd23;
              s000[i3m][i2m][i1p] += hd22+hd33+hd23; // mmp
              s000[i3m][i2p][i1p] += hd22+hd33+hd32; // mpp
              s000[i3p][i2m][i1p] += hd22+hd33+hd32; // pmp
              s000[i3p][i2p][i1p] += hd22+hd33+hd23; // ppp
              s0p0[i3m][i2m][i1p] -= hd22; // mpp
              s0p0[i3p][i2m][i1p] -= hd22; // ppp
              sp00[i3m][i2m][i1p] -= hd33; // pmp
              sp00[i3m][i2p][i1p] -= hd33; // ppp
              spp0[i3m][i2m][i1p] -= hd23; // ppp
              spm0[i3m][i2p][i1p] -= hd32; // pmp
              if (i2m==0) {
                s000[i3m][i2m][i1p] += d33; // mmp
                s000[i3p][i2m][i1p] += d33; // pmp
                sp00[i3m][i2m][i1p] -= d33; // pmp
              } else if (i2p==n2m) {
                s000[i3m][i2p][i1p] += d33; // mpp
                s000[i3p][i2p][i1p] += d33; // ppp
                sp00[i3m][i2p][i1p] -= d33; // ppp
              }
              if (i3m==0) {
                s000[i3m][i2m][i1p] += d22; // mmp
                s000[i3m][i2p][i1p] += d22; // mpp
                s0p0[i3m][i2m][i1p] -= d22; // mpp
              } else if (i3p==n3m) {
                s000[i3p][i2m][i1p] += d22; // pmp
                s000[i3p][i2p][i1p] += d22; // ppp
                s0p0[i3p][i2m][i1p] -= d22; // ppp
              }
            }
            if (i2m==0) {
              hd11 = 0.5f*d11;
              hd33 = 0.5f*d33;
              hd13 = 0.5f*d13;
              hd31 = -hd13;
              s000[i3m][i2m][i1m] += hd11+hd33+hd13; // mmm
              s000[i3m][i2m][i1p] += hd11+hd33+hd31; // mmp
              s000[i3p][i2m][i1m] += hd11+hd33+hd31; // pmm
              s000[i3p][i2m][i1p] += hd11+hd33+hd13; // pmp
              s00p[i3m][i2m][i1m] -= hd11; // mmp
              s00p[i3p][i2m][i1m] -= hd11; // pmp
              sp00[i3m][i2m][i1m] -= hd33; // pmm
              sp00[i3m][i2m][i1p] -= hd33; // pmp
              sp0p[i3m][i2m][i1m] -= hd13; // pmp
              sp0m[i3m][i2m][i1p] -= hd31; // pmm
              if (i1m==0) {
                s000[i3m][i2m][i1m] += d33; // mmm
                s000[i3p][i2m][i1m] += d33; // pmm
                sp00[i3m][i2m][i1m] -= d33; // pmm
              } else if (i1p==n1m) {
                s000[i3m][i2m][i1p] += d33; // mmp
                s000[i3p][i2m][i1p] += d33; // pmp
                sp00[i3m][i2m][i1p] -= d33; // pmp
              }
              if (i3m==0) {
                s000[i3m][i2m][i1m] += d11; // mmm
                s000[i3m][i2m][i1p] += d11; // mmp
                s00p[i3m][i2m][i1m] -= d11; // mmp
              } else if (i3p==n3m) {
                s000[i3p][i2m][i1m] += d11; // pmm
                s000[i3p][i2m][i1p] += d11; // pmp
                s00p[i3p][i2m][i1m] -= d11; // pmp
              }
            } else if (i2p==n2m) {
              hd11 = 0.5f*d11;
              hd33 = 0.5f*d33;
              hd13 = 0.5f*d13;
              hd31 = -hd13;
              s000[i3m][i2p][i1m] += hd11+hd33+hd13; // mpm
              s000[i3m][i2p][i1p] += hd11+hd33+hd31; // mpp
              s000[i3p][i2p][i1m] += hd11+hd33+hd31; // ppm
              s000[i3p][i2p][i1p] += hd11+hd33+hd13; // ppp
              s00p[i3m][i2p][i1m] -= hd11; // mpp
              s00p[i3p][i2p][i1m] -= hd11; // ppp
              sp00[i3m][i2p][i1m] -= hd33; // ppm
              sp00[i3m][i2p][i1p] -= hd33; // ppp
              sp0p[i3m][i2p][i1m] -= hd13; // ppp
              sp0m[i3m][i2p][i1p] -= hd31; // ppm
              if (i1m==0) {
                s000[i3m][i2p][i1m] += d33; // mpm
                s000[i3p][i2p][i1m] += d33; // ppm
                sp00[i3m][i2p][i1m] -= d33; // ppm
              } else if (i1p==n1m) {
                s000[i3m][i2p][i1p] += d33; // mpp
                s000[i3p][i2p][i1p] += d33; // ppp
                sp00[i3m][i2p][i1p] -= d33; // ppp
              }
              if (i3m==0) {
                s000[i3m][i2p][i1m] += d11; // mpm
                s000[i3m][i2p][i1p] += d11; // mpp
                s00p[i3m][i2p][i1m] -= d11; // mpp
              } else if (i3p==n3m) {
                s000[i3p][i2p][i1m] += d11; // ppm
                s000[i3p][i2p][i1p] += d11; // ppp
                s00p[i3p][i2p][i1m] -= d11; // ppp
              }
            }
            if (i3m==0) {
              hd11 = 0.5f*d11;
              hd22 = 0.5f*d22;
              hd12 = 0.5f*d12;
              hd21 = -hd12;
              s000[i3m][i2m][i1m] += hd11+hd22+hd12; // mmm
              s000[i3m][i2m][i1p] += hd11+hd22+hd21; // mmp
              s000[i3m][i2p][i1m] += hd11+hd22+hd21; // mpm
              s000[i3m][i2p][i1p] += hd11+hd22+hd12; // mpp
              s00p[i3m][i2m][i1m] -= hd11; // mmp
              s00p[i3m][i2p][i1m] -= hd11; // mpp
              s0p0[i3m][i2m][i1m] -= hd22; // mpm
              s0p0[i3m][i2m][i1p] -= hd22; // mpp
              s0pp[i3m][i2m][i1m] -= hd12; // mpp
              s0pm[i3m][i2m][i1p] -= hd21; // mpm
              if (i1m==0) {
                s000[i3m][i2m][i1m] += d22; // mmm
                s000[i3m][i2p][i1m] += d22; // mpm
                s0p0[i3m][i2m][i1m] -= d22; // mpm
              } else if (i1p==n1m) {
                s000[i3m][i2m][i1p] += d22; // mmp
                s000[i3m][i2p][i1p] += d22; // mpp
                s0p0[i3m][i2m][i1p] -= d22; // mpp
              }
              if (i2m==0) {
                s000[i3m][i2m][i1m] += d11; // mmm
                s000[i3m][i2m][i1p] += d11; // mmp
                s00p[i3m][i2m][i1m] -= d11; // mmp
              } else if (i2p==n2m) {
                s000[i3m][i2p][i1m] += d11; // mpm
                s000[i3m][i2p][i1p] += d11; // mpp
                s00p[i3m][i2p][i1m] -= d11; // mpp
              }
            } else if (i3p==n3m) {
              hd11 = 0.5f*d11;
              hd22 = 0.5f*d22;
              hd12 = 0.5f*d12;
              hd21 = -hd12;
              s000[i3p][i2m][i1m] += hd11+hd22+hd12; // pmm
              s000[i3p][i2m][i1p] += hd11+hd22+hd21; // pmp
              s000[i3p][i2p][i1m] += hd11+hd22+hd21; // ppm
              s000[i3p][i2p][i1p] += hd11+hd22+hd12; // ppp
              s00p[i3p][i2p][i1m] -= hd11; // ppp
              s00p[i3p][i2m][i1m] -= hd11; // pmp
              s0p0[i3p][i2m][i1m] -= hd22; // ppm
              s0p0[i3p][i2m][i1p] -= hd22; // ppp
              s0pp[i3p][i2m][i1m] -= hd12; // ppp
              s0pm[i3p][i2m][i1p] -= hd21; // ppm
              if (i1m==0) {
                s000[i3p][i2m][i1m] += d22; // pmm
                s000[i3p][i2p][i1m] += d22; // ppm
                s0p0[i3p][i2m][i1m] -= d22; // ppm
              } else if (i1p==n1m) {
                s000[i3p][i2m][i1p] += d22; // pmp
                s000[i3p][i2p][i1p] += d22; // ppp
                s0p0[i3p][i2m][i1p] -= d22; // ppp
              }
              if (i2m==0) {
                s000[i3p][i2m][i1m] += d11; // pmm
                s000[i3p][i2m][i1p] += d11; // pmp
                s00p[i3p][i2m][i1m] -= d11; // pmp
              } else if (i2p==n2m) {
                s000[i3p][i2p][i1m] += d11; // ppm
                s000[i3p][i2p][i1p] += d11; // ppp
                s00p[i3p][i2p][i1m] -= d11; // ppp
              }
            }
          }
        }
      }
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean ZERO_SLOPE_BOUNDARIES = true;
  private float _rs; // experimental parameter for gradient approximation
  private float _ers,_frs; // for 2D 2x2-sample stencils
  private float _erst,_frst,_grst; // for 3D 2x2x2-sample stencils

  // This version is more like the one used for 3D 2x2x2 stencils.
  // But the version above works, so this one is currently unused.
  private void applyX(LocalDiffusionTensors2 ldt, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float e = _ers;
    float f = _frs;
    float fp5 = f*0.5f;
    float[] d = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
        ldt.getTensor(i1p,i2p,d);
        float d11 = d[0];
        float d12 = d[1];
        float d22 = d[2];
        float edd = e*(d11+d22);
        float p12 = 0.50f*d12;
        float h11 = fp5*d11-e*d22;
        float h22 = fp5*d22-e*d11;
        float h12 = edd+p12;
        float h21 = edd-p12;
        float xpp = x[i2p][i1p];
        float xpm = x[i2p][i1m];
        float xmp = x[i2m][i1p];
        float xmm = x[i2m][i1m];
        float ypppm = h11*(xpp-xpm);
        float ympmm = h11*(xmp-xmm);
        float ypmmm = h22*(xpm-xmm);
        float yppmp = h22*(xpp-xmp);
        float yppmm = h12*(xpp-xmm);
        float ymppm = h21*(xmp-xpm);
        y[i2p][i1p] += ypppm+yppmp+yppmm;
        y[i2p][i1m] -= ypppm-ypmmm+ymppm;
        y[i2m][i1p] += ympmm-yppmp+ymppm;
        y[i2m][i1m] -= ympmm+ypmmm-yppmm;
        if (ZERO_SLOPE_BOUNDARIES) {
          if (i1m==0) {
            ypmmm = d22*(xpm-xmm);
            y[i2p][i1m] += ypmmm;
            y[i2m][i1m] -= ypmmm;
          } else if (i1p==n1m) {
            yppmp = d22*(xpp-xmp);
            y[i2p][i1p] += yppmp;
            y[i2m][i1p] -= yppmp;
          }
          if (i2m==0) {
            ympmm = d11*(xmp-xmm);
            y[i2m][i1p] += ympmm;
            y[i2m][i1m] -= ympmm;
          } else if (i2p==n2m) {
            ypppm = d11*(xpp-xpm);
            y[i2p][i1p] += ypppm;
            y[i2p][i1m] -= ypppm;
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void testGetCoefficients2() {
    int n1 = 5;
    int n2 = 5;
    float[][] x = Array.zerofloat(n1,n2); 
    //float[][] x = Array.randfloat(n1,n2); 
    //x[0][0] = x[n2-1][0] = x[0][n1-1] = x[n2-1][n1-1] = 1.0f;
    x[n2/2][n1/2] = 1.0f;
    float[][] y = Array.zerofloat(n1,n2);
    float[][] z = Array.zerofloat(n1,n2);
    float theta = FLT_PI*0.0f/8.0f;
    float[][] d0 = Array.fillfloat(1.0f,n1,n2);
    float[][] d1 = Array.fillfloat(1.0f,n1,n2);
    float[][] v1 = Array.fillfloat(sin(theta),n1,n2);
    LocalDiffusionTensors2 ldt = new LocalDiffusionTensors2(0.0,1.0,d0,d1,v1);
    float rs = 3.0f/12.0f;
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(rs);
    LocalSpd9Filter lsf = new LocalSpd9Filter(ldk.getCoefficients(ldt));
    ldk.apply(ldt,x,y);
    lsf.apply(x,z);
    Array.dump(x);
    Array.dump(y);
    Array.dump(z);
    System.out.println("error = "+Array.sum(Array.abs(Array.sub(y,z))));
  }

  private static void testGetCoefficients3() {
    int n1 = 5;
    int n2 = 4;
    int n3 = 3;
    //float[][][] x = Array.randfloat(n1,n2,n3); 
    float[][][] x = Array.zerofloat(n1,n2,n3); 
    x[n3/2][n2/2][n1/2] = 1.0f;
    //for (int i3=0; i3<n3; ++i3)
    //  x[i3][n2/2][n1/2] = 1.0f;
    float[][][] y = Array.zerofloat(n1,n2,n3); 
    float[][][] z = Array.zerofloat(n1,n2,n3); 
    //DiffusionTensors3 ldt = makeRandomDiffusionTensors3(n1,n2,n3);
    float theta = FLT_PI*2.0f/8.0f;
    float phi = FLT_PI*0.0f/8.0f;
    //DiffusionTensors3 ldt = makePlanarDiffusionTensors3(n1,n2,n3,theta,phi);
    DiffusionTensors3 ldt = makeLinearDiffusionTensors3(n1,n2,n3,theta,phi);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    LocalSpd27Filter lsf = new LocalSpd27Filter(ldk.getCoefficients(ldt));
    ldk.apply(ldt,x,y);
    lsf.apply(x,z);
    Array.dump(x);
    Array.dump(y);
    Array.dump(z);
    System.out.println("error = "+Array.sum(Array.abs(Array.sub(y,z))));
    System.out.println("sum y = "+Array.sum(y));
    System.out.println("sum z = "+Array.sum(z));
  }

  private static DiffusionTensors3 makeLinearDiffusionTensors3(
    int n1, int n2, int n3, float theta, float phi) 
  {
    DiffusionTensors3 ldt = new DiffusionTensors3(n1,n2,n3,1.0f,1.0f,1.0f);
    float d1 = 1.000f;
    float d2 = 0.000f;
    float d3 = 0.000f;
    float w1 = cos(theta);
    float w2 = sin(theta)*cos(phi);
    float w3 = sin(theta)*sin(phi);
    float[] d = {d1,d2,d3};
    float[] w = {w1,w2,w3};
    float[] u = makeOrthogonalVector(w);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          ldt.setCoefficients(i1,i2,i3,d);
          ldt.setEigenvectorU(i1,i2,i3,u);
          ldt.setEigenvectorW(i1,i2,i3,w);
        }
      }
    }
    return ldt;
  }
  private static DiffusionTensors3 makePlanarDiffusionTensors3(
    int n1, int n2, int n3, float theta, float phi) 
  {
    DiffusionTensors3 ldt = new DiffusionTensors3(n1,n2,n3,1.0f,1.0f,1.0f);
    float d1 = 0.000f;
    float d2 = 1.000f;
    float d3 = 0.000f;
    float u1 = cos(theta);
    float u2 = sin(theta)*cos(phi);
    float u3 = sin(theta)*sin(phi);
    float[] d = {d1,d2,d3};
    float[] u = {u1,u2,u3};
    float[] w = makeOrthogonalVector(u);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          ldt.setCoefficients(i1,i2,i3,d);
          ldt.setEigenvectorU(i1,i2,i3,u);
          ldt.setEigenvectorW(i1,i2,i3,w);
        }
      }
    }
    return ldt;
  }
  private static DiffusionTensors3 makeRandomDiffusionTensors3(
    int n1, int n2, int n3) 
  {
    DiffusionTensors3 ldt = new DiffusionTensors3(n1,n2,n3,1.0f,1.0f,1.0f);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] d = makeRandomCoefficients();
          float[] u = makeRandomVector();
          float[] w = makeOrthogonalVector(u);
          ldt.setCoefficients(i1,i2,i3,d);
          ldt.setEigenvectorU(i1,i2,i3,u);
          ldt.setEigenvectorW(i1,i2,i3,w);
        }
      }
    }
    return ldt;
  }
  private static java.util.Random r = new java.util.Random();
  private static float[] makeRandomCoefficients() {
    float d1 = r.nextFloat();
    float d2 = r.nextFloat();
    float d3 = r.nextFloat();
    float ds = 1.0f/(d1+d2+d3);
    return new float[]{d1*ds,d2*ds,d3*ds};
  }
  private static float[] makeRandomVector() {
    float a = r.nextFloat()-0.5f;
    float b = r.nextFloat()-0.5f;
    float c = r.nextFloat()-0.5f;
    float s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }
  private static float[] makeOrthogonalVector(float[] v1) {
    float a1 = v1[0];
    float b1 = v1[1];
    float c1 = v1[2];
    float a2 = r.nextFloat()-0.5f;
    float b2 = r.nextFloat()-0.5f;
    float c2 = r.nextFloat()-0.5f;
    float d11 = a1*a1+b1*b1+c1*c1;
    float d12 = a1*a2+b1*b2+c1*c2;
    float s = d12/d11;
    float a = a2-s*a1;
    float b = b2-s*b1;
    float c = c2-s*c1;
    s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }


  public static void main(String[] args) {
    //testGetCoefficients2();
    testGetCoefficients3();
  }
}
