/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package het;

import static edu.mines.jtk.util.ArrayMath.*;

// for testing only
import edu.mines.jtk.mosaic.*;

/**
 * Local estimation of peak frequency.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.03.25
 */
public class LocalPeakFrequencyFilter {

  public LocalPeakFrequencyFilter(double sigma) {
    _sigma = (float)sigma;
  }

  public float[] findPeakFrequencies(float[] x) {
    int n = x.length;
    float[][] c = burg(_sigma,2,x);
    float[] c0 = c[0];
    float[] c1 = c[1];
    //plot("c0",c0);
    //plot("c1",c1);
    float[] pf = new float[n];
    float o2pi = 0.5f/FLT_PI;
    for (int i=0; i<n; ++i) {
      float c0i = c0[i];
      float c1i = c1[i];
      float cnum = -c0i*(1.0f-c1i)*(1.0f-c1i);
      float cden = 4.0f*c1i;
      float pw = (abs(cnum)<=abs(cden))?acos(cnum/cden):0.0f;
      pf[i] = pw*o2pi;
    }
    return pf;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float SIGMA_MIN = sqrt(2.0f)/2.0f;

  private float _sigma;

  /**
   * Two-sided exponential smoothing with zero-value or zero-slope ends.
   */
  private static void smooth(boolean zs, float sigma, float[] x, float[] y) {
    if (x==y) x = copy(x);
    if (sigma<SIGMA_MIN) {
      copy(x,y);
      return;
    }
    float a = (sigma-SIGMA_MIN)/(sigma+SIGMA_MIN);
    int n1 = x.length;
    float s = (1.0f-a)/(1.0f+a);
    y[0] = zs?x[0]*s/(1.0f-a):x[0]*s;
    for (int i1=1; i1<n1; ++i1)
      y[i1] = a*y[i1-1]+s*x[i1];
    float yip1 = zs?x[n1-1]*s*a/(1.0f-a):0.0f;
    y[n1-1] += yip1;
    for (int i1=n1-2; i1>=0; --i1) {
      float yi = a*(yip1+s*x[i1+1]);
      y[i1] += yi;
      yip1 = yi;
    }
  }

  /**
   * Returns arrays of local Burg reflection coefficients.
   */
  private static float[][] burg(float sigma, int m, float[] x) {
    int n = x.length;
    float[][] c = new float[m][n];
    float[] f = copy(x);
    float[] b = copy(x);
    float[] bf = new float[n];
    float[] bb = new float[n];
    float[] ff = new float[n];
    for (int k=0; k<m; ++k) {
      for (int i=n-1; i>k; --i) {
        float bi = b[i-1];
        float fi = f[i  ];
        bf[i] = bi*fi;
        bb[i] = bi*bi;
        ff[i] = fi*fi;
      }
      for (int i=k; i>=0; --i) {
        bf[i] = 0.0f;
        bb[i] = 0.0f;
        ff[i] = 0.0f;
      }
      smooth(false,sigma,bf,bf);
      smooth(false,sigma,bb,bb);
      smooth(false,sigma,ff,ff);
      for (int i=n-1; i>k; --i) {
        float ci = (2.0f*bf[i])/(bb[i]+ff[i]);
        //assert ci>=-1.0f:"ci>=-1.0f";
        //assert ci<= 1.0f:"ci<= 1.0f";
        float bi = b[i-1];
        float fi = f[i  ];
        f[i] = fi-ci*bi;
        b[i] = bi-ci*fi;
        c[k][i] = ci;
      }
      for (int i=k; i>=0; --i)
        c[k][i] = c[k][k+1];
    }
    return c;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static float[] sweep(int n, double f1, double f2) {
    float twopi = 2.0f*FLT_PI;
    float wa = (float)(twopi*f1);
    float wb = (float)(twopi*0.5*(f2-f1)/(n-1));
    float[] x = new float[n];
    for (int i=0; i<n; ++i) {
      float wi = wa+wb*i;
      x[i] = sin(wi*i);
    }
    return x;
  }

  public static void main(String[] args) {
    int n = 201;
    float w = 2.0f*FLT_PI*0.2f;
    float[] x = sweep(n,0.1,0.3);
    LocalPeakFrequencyFilter lpff = new LocalPeakFrequencyFilter(10);
    float[] f = lpff.findPeakFrequencies(x);
    plot("x",x);
    plot("f",f,0.0f,0.4f);
  }

  public static void plot(String title, float[] x) {
    plot(title,x,min(x),max(x));
  }

  public static void plot(String title, float[] x, float xmin, float xmax) {
    SimplePlot sp = SimplePlot.asPoints(x);
    sp.setTitle(title);
    sp.setVLimits(xmin,xmax);
  }
}
