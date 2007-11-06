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
 * Local symmetric positive-definite filter with a 3x3 (9-point) stencil.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.05
 */
public class LocalSpdFilter33 {

  public interface Stencils22 {

    /**
     * Gets stencils for four 2x2 filters, one for each quadrant.
     * @param i1 sample index in 1st dimension.
     * @param i2 sample index in 2nd dimension.
     * @param abcd array[4][4] of stencil coefficients.
     *  Array a[0] contains coefficients {a,b,c,d} for 0th quadrant;
     *  array a[1] contains coefficients {a,b,c,d} for 1st quadrant;
     *  and so on.
     */
    public void getStencils(int i1, int i2, float[][] abcd);
  }

  /**
   * Constructs a local symmetric positive-definite 3x3 filter.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   */
  public LocalSpdFilter33(int n1, int n2, Stencils22 s22) {
    _a = new float[5][n2][n1];
    float[][] a00 = _a[0];
    float[][] a0p = _a[1];
    float[][] apm = _a[2];
    float[][] ap0 = _a[3];
    float[][] app = _a[4];
    float[][] abcd = new float[4][4];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        s22.getStencils(i1,i2,abcd);
        float a0=abcd[0][0], b0=abcd[0][1], c0=abcd[0][2], d0=abcd[0][3];
        float a1=abcd[1][0], b1=abcd[1][1], c1=abcd[1][2], d1=abcd[1][3];
        float a2=abcd[2][0], b2=abcd[2][1], c2=abcd[2][2], d2=abcd[2][3];
        float a3=abcd[3][0], b3=abcd[3][1], c3=abcd[3][2], d3=abcd[3][3];
        a00[i2  ][i1  ] += a0*a0+a1*a1+a2*a2+a3*a3;
        a0p[i2  ][i1  ] += a0*b0+a2*b2;
        ap0[i2  ][i1  ] += a0*c0+a1*c1;
        apm[i2  ][i1  ] += a1*d1;
        app[i2  ][i1  ] += a0*d0;
        if (0<i1) {
          a00[i2  ][i1-1] += b1*b1+b3*b3;
          a0p[i2  ][i1-1] += a1*b1+a3*b3;
          ap0[i2  ][i1-1] += b1*d1;
          app[i2  ][i1-1] += b1*c1;
          if (0<i2) {
            a00[i2-1][i1-1] += d3*d3;
            a0p[i2-1][i1-1] += c3*d3;
            ap0[i2-1][i1-1] += b3*d3;
            app[i2-1][i1-1] += a3*d3;
          }
          if (i2<n2-1) {
            a00[i2+1][i1-1] += d1*d1;
            a0p[i2+1][i1-1] += c1*d1;
          }
        }
        if (i1<n1-1) {
          a00[i2  ][i1+1] += b0*b0+b2*b2;
          apm[i2  ][i1+1] += b0*c0;
          ap0[i2  ][i1+1] += b0*d0;
          if (0<i2) {
            a00[i2-1][i1+1] += d2*d2;
            apm[i2-1][i1+1] += a2*d2;
            ap0[i2-1][i1+1] += b2*d2;
          }
          if (i2<n2-1) {
            a00[i2+1][i1+1] += d0*d0;
          }
        }
        if (0<i2) {
          a00[i2-1][i1  ] += c2*c2+c3*c3;
          ap0[i2-1][i1  ] += a2*c2+a3*c3;
          a0p[i2-1][i1  ] += c2*d2;
          apm[i2-1][i1  ] += b3*c3;
          app[i2-1][i1  ] += b2*c2;
        }
        if (i2<n2-1) {
          a00[i2+1][i1  ] += c0*c0+c1*c1;
          a0p[i2+1][i1  ] += c0*d0;
        }

        /*
        a00[i2  ][i1  ] += a0*a0+a1*a1+a2*a2+a3*a3;
        a00[i2  ][i1+1] += b0*b0+b2*b2;
        a00[i2  ][i1-1] += b1*b1+b3*b3;
        a00[i2+1][i1  ] += c0*c0+c1*c1;
        a00[i2-1][i1  ] += c2*c2+c3*c3;
        a00[i2+1][i1+1] += d0*d0;
        a00[i2+1][i1-1] += d1*d1;
        a00[i2-1][i1+1] += d2*d2;
        a00[i2-1][i1-1] += d3*d3;

        a0p[i2  ][i1  ] += a0*b0+a2*b2;
        a0p[i2  ][i1-1] += a1*b1+a3*b3;
        a0p[i2+1][i1  ] += c0*d0;
        a0p[i2-1][i1  ] += c2*d2;
        a0p[i2+1][i1-1] += c1*d1;
        a0p[i2-1][i1-1] += c3*d3;

        ap0[i2  ][i1  ] += a0*c0+a1*c1;
        ap0[i2-1][i1  ] += a2*c2+a3*c3;
        ap0[i2  ][i1+1] += b0*d0;
        ap0[i2  ][i1-1] += b1*d1;
        ap0[i2-1][i1+1] += b2*d2;
        ap0[i2-1][i1-1] += b3*d3;

        apm[i2  ][i1  ] += a1*d1;
        apm[i2-1][i1+1] += a2*d2;
        apm[i2  ][i1+1] += b0*c0;
        apm[i2-1][i1  ] += b3*c3;

        app[i2  ][i1  ] += a0*d0;
        app[i2-1][i1-1] += a3*d3;
        app[i2  ][i1-1] += b1*c1;
        app[i2-1][i1  ] += b2*c2;
        */
      }
    }
  }

  public void apply(float[][] x, float[][] y) {
    int n1 = _a[0][0].length;
    int n2 = _a[0].length;
    float[][] a00 = _a[0];
    float[][] a0p = _a[1];
    float[][] apm = _a[2];
    float[][] ap0 = _a[3];
    float[][] app = _a[4];
    float t00,t0p,tpm,tp0,tpp;
    float x00,x0p,xpm,xp0,xpp;
    int i1 = 0;
    int i2 = 0;
    for (i2=0; i2<n2-1; ++i2) {
      t00 = a00[i2][i1];
      t0p = a0p[i2][i1];
      tp0 = ap0[i2][i1];
      tpp = app[i2][i1];
      x00 = x[i2  ][i1  ];
      x0p = x[i2  ][i1+1];
      xp0 = x[i2+1][i1  ];
      xpp = x[i2+1][i1+1];
      y[i2  ][i1  ] += t00*x00+t0p*x0p+tp0*xp0+tpp*xpp;
      y[i2  ][i1+1] += t0p*x00;
      y[i2+1][i1  ] += tp0*x00;
      y[i2+1][i1+1] += tpp*x00;
      for (i1=1; i1<n1-1; ++i1) {
        t00 = a00[i2][i1];
        t0p = a0p[i2][i1];
        tpm = apm[i2][i1];
        tp0 = ap0[i2][i1];
        tpp = app[i2][i1];
        x00 = x[i2  ][i1  ];
        x0p = x[i2  ][i1+1];
        xpm = x[i2+1][i1-1];
        xp0 = x[i2+1][i1  ];
        xpp = x[i2+1][i1+1];
        y[i2  ][i1  ] += t00*x00+t0p*x0p+tpm*xpm+tp0*xp0+tpp*xpp;
        y[i2  ][i1+1] += t0p*x00;
        y[i2+1][i1-1] += tpm*x00;
        y[i2+1][i1  ] += tp0*x00;
        y[i2+1][i1+1] += tpp*x00;
      }
      t00 = a00[i2][i1];
      tpm = apm[i2][i1];
      tp0 = ap0[i2][i1];
      x00 = x[i2  ][i1  ];
      xpm = x[i2+1][i1-1];
      xp0 = x[i2+1][i1  ];
      y[i2  ][i1  ] += t00*x00+tpm*xpm+tp0*xp0;
      y[i2+1][i1-1] += tpm*x00;
      y[i2+1][i1  ] += tp0*x00;
    }
    for (i1=0; i1<n1-1; ++i1) {
      t00 = a00[i2][i1];
      t0p = a0p[i2][i1];
      x00 = x[i2  ][i1  ];
      x0p = x[i2  ][i1+1];
      y[i2  ][i1  ] += t00*x00+t0p*x0p;
      y[i2  ][i1+1] += t0p*x00;
    }
    t00 = a00[i2][i1];
    x00 = x[i2  ][i1  ];
    y[i2  ][i1  ] += t00*x00;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Coefficients for finite-difference approximation.
  //private static final float R = 0.5f*(1.0f+sqrt(2.0f/3.0f));
  //private static final float S = 0.5f*(1.0f-sqrt(2.0f/3.0f));
  private static final float R = 0.5f;
  private static final float S = 0.5f;

  private float[][][] _a; // stencil coefficients

  /**
   * Solves L*D*L'*x = b, for elements of factors L, D, and L' stored in a.
   */
  private static void solveWithFactors(
    float[][][] a, float[][] b, float[][] x) 
  {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    float[][] d00 = a[0];
    int i1,i2;
    float yi,zi;

    // Solve L*z = b.
    float[][] z = x;
    Array.zero(z);
    for (i2=0; i2<n2-1; ++i2) {
      i1 = 0;
      zi = d00[i2][i1]*b[i2][i1];
      z[i2  ][i1  ] = zi;
      z[i2  ][i1+1] -= a0p[i2][i1]*zi;
      z[i2+1][i1  ] -= ap0[i2][i1]*zi;
      z[i2+1][i1+1] -= app[i2][i1]*zi;
      for (i1=1; i1<n1-1; ++i1) {
        zi = d00[i2][i1]*b[i2][i1];
        z[i2  ][i1  ] = zi;
        z[i2  ][i1+1] -= a0p[i2][i1]*zi;
        z[i2+1][i1-1] -= apm[i2][i1]*zi;
        z[i2+1][i1  ] -= ap0[i2][i1]*zi;
        z[i2+1][i1+1] -= app[i2][i1]*zi;
      }
      zi = d00[i2][i1]*b[i2][i1];
      z[i2  ][i1  ] = zi;
      z[i2+1][i1-1] -= apm[i2][i1]*zi;
      z[i2+1][i1  ] -= ap0[i2][i1]*zi;
    }
    i1 = 0;
    zi = d00[i2][i1]*b[i2][i1];
    z[i2  ][i1  ] = zi;
    z[i2  ][i1+1] -= a0p[i2][i1]*zi;
    for (i1=1; i1<n1-1; ++i1) {
      zi = d00[i2][i1]*b[i2][i1];
      z[i2  ][i1  ] = zi;
      z[i2  ][i1+1] -= a0p[i2][i1]*zi;
    }
    zi = d00[i2][i1]*b[i2][i1];
    z[i2  ][i1  ] = zi;

    // Solve D*y = z and L'*x = y.
    for (i2=n2-1; i2>=0; --i2) {
      for (i1=n1-1; i1>=0; --i1) {
        yi = z[i2][i1]/d00[i2][i1];
        yi -= app[i2][i1]*x[i2+1][i1+1];
        yi -= ap0[i2][i1]*x[i2+1][i1  ];
        yi -= apm[i2][i1]*x[i2+1][i1-1];
        yi -= a0p[i2][i1]*x[i2  ][i1+1];
        x[i2][i1] = d00[i2][i1]*yi;
      }
    }
  }


  /**
   * Makes factors of incomplete Cholesky decomposition IC(0) (zero fill-in).
   * Elements of the diagonal matrix d00 (inverse of a00) are stored in a[0].
   */
  private static void makeFactors(float[][][] a) {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] a00 = a[0];
    float[][] a0p = a[1];
    float[][] apm = a[2];
    float[][] ap0 = a[3];
    float[][] app = a[4];
    float[][] d00 = a[0];
    int i1 = 0;
    int i2 = 0;
    d00[i2][i1] = 1.0f/a00[i2][i1];
    for (i1=1; i1<n1; ++i1) {
      a00[i2][i1] -= d00[i2  ][i1-1]*a0p[i2  ][i1-1]*a0p[i2  ][i1-1];
      apm[i2][i1] -= d00[i2  ][i1-1]*ap0[i2  ][i1-1]*a0p[i2  ][i1-1];
      ap0[i2][i1] -= d00[i2  ][i1-1]*app[i2  ][i1-1]*a0p[i2  ][i1-1];
      d00[i2][i1] = 1.0f/a00[i2][i1];
    }
    for (i2=1; i2<n2; ++i2) {
      i1 = 0;
      a00[i2][i1] -= d00[i2-1][i1+1]*apm[i2-1][i1+1]*apm[i2-1][i1+1] +
                     d00[i2-1][i1  ]*ap0[i2-1][i1  ]*ap0[i2-1][i1  ];
      a0p[i2][i1] -= d00[i2-1][i1  ]*app[i2-1][i1  ]*ap0[i2-1][i1  ];
      d00[i2][i1] = 1.0f/a00[i2][i1];
      for (i1=1; i1<n1-1; ++i1) {
        a00[i2][i1] -= d00[i2  ][i1-1]*a0p[i2  ][i1-1]*a0p[i2  ][i1-1] +
                       d00[i2-1][i1+1]*apm[i2-1][i1+1]*apm[i2-1][i1+1] +
                       d00[i2-1][i1  ]*ap0[i2-1][i1  ]*ap0[i2-1][i1  ] +
                       d00[i2-1][i1-1]*app[i2-1][i1-1]*app[i2-1][i1-1];
        a0p[i2][i1] -= d00[i2-1][i1+1]*ap0[i2-1][i1+1]*apm[i2-1][i1+1] +
                       d00[i2-1][i1  ]*app[i2-1][i1  ]*ap0[i2-1][i1  ];
        apm[i2][i1] -= d00[i2  ][i1-1]*ap0[i2  ][i1-1]*a0p[i2  ][i1-1];
        ap0[i2][i1] -= d00[i2  ][i1-1]*app[i2  ][i1-1]*a0p[i2  ][i1-1];
        d00[i2][i1] = 1.0f/a00[i2][i1];
      }
      a00[i2][i1] -= d00[i2  ][i1-1]*a0p[i2  ][i1-1]*a0p[i2  ][i1-1] +
                     d00[i2-1][i1  ]*ap0[i2-1][i1  ]*ap0[i2-1][i1  ] +
                     d00[i2-1][i1-1]*app[i2-1][i1-1]*app[i2-1][i1-1];
      a0p[i2][i1] -= d00[i2-1][i1  ]*app[i2-1][i1  ]*ap0[i2-1][i1  ];
      apm[i2][i1] -= d00[i2  ][i1-1]*ap0[i2  ][i1-1]*a0p[i2  ][i1-1];
      ap0[i2][i1] -= d00[i2  ][i1-1]*app[i2  ][i1-1]*a0p[i2  ][i1-1];
      d00[i2][i1] = 1.0f/a00[i2][i1];
    }
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
  }
} 
