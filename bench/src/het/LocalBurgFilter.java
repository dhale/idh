/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package het;

import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

// for testing only
import edu.mines.jtk.mosaic.*;

/**
 * Local Burg forward and inverse prediction error filtering.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.03.29
 */
public class LocalBurgFilter {

  public float[][] getCoefficients(double sigma, int m, float[] x) {
    int n = x.length;
    float[][] c = new float[m][n];
    float[] f = copy(x);
    float[] b = copy(x);
    float[] cn = new float[n];
    float[] cd = new float[n];
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    for (int k=0; k<m; ++k) {
      computeC1(k,f,b,cn,cd);
      ref.apply(cn,cn);
      ref.apply(cd,cd);
      computeC2(k,cn,cd,c[k],f,b);
    }
    return c;
  }

  public float[][][] getCoefficients(
    double sigma1, double sigma2, int m, float[][] x) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][][] c = new float[m][n2][n1];
    float[][] f = copy(x);
    float[][] b = copy(x);
    float[][] cn = new float[n2][n1];
    float[][] cd = new float[n2][n1];
    RecursiveExponentialFilter ref = new
      RecursiveExponentialFilter(sigma1,sigma2);
    for (int k=0; k<m; ++k) {
      for (int i2=0; i2<n2; ++i2) {
        computeC1(k,b[i2],f[i2],cn[i2],cd[i2]);
      }
      ref.apply(cn,cn);
      ref.apply(cd,cd);
      for (int i2=0; i2<n2; ++i2) {
        computeC2(k,cn[i2],cd[i2],c[k][i2],b[i2],f[i2]);
      }
    }
    return c;
  }

  public float[] getPeakFrequencies(float[][] c) {
    return getPeakFrequencies(c[0],c[1]);
  }

  public float[][] getPeakFrequencies(float[][][] c) {
    int n2 = c[0].length;
    float[][] pf = new float[n2][];
    for (int i2=0; i2<n2; ++i2) {
      pf[i2] = getPeakFrequencies(c[0][i2],c[1][i2]);
    }
    return pf;
  }

  public float[][] shiftPeakFrequencies(double df, float[][] c) {
    float dw = (float)(2.0*PI*df);
    float[] c1 = c[0];
    float[] c2 = c[1];
    c1 = copy(c1);
    int n = c1.length;
    for (int i=0; i<n; ++i) {
      float c1i = c1[i];
      c1[i] = c1i-sqrt(1.0f-c1i*c1i)*dw;
    }
    return new float[][]{c1,c2};
  }

  public float[][] zeroPeakFrequencies(float[][] c) {
    float[] pf = getPeakFrequencies(c);
    float sf = 2.0f*FLT_PI;
    int n = c[0].length;
    float[] c1 = new float[n];
    float[] c2 = new float[n];
    for (int i=0; i<n; ++i) {
      c1[i] = cos(sf*pf[i]);
      c2[i] = -1.0f;
    }
    return new float[][]{c1,c2};
  }

  public float[] applyForward(float[][] c, float[] x) {
    float[] y = new float[x.length];
    applyForward(c,x,y);
    return y;
  }

  public void applyForward(float[][] c, float[] x, float[] y) {
    int m = c.length;
    if (m==2) {
      applyForward(c[0],c[1],x,y);
    } else {
      assert false:"implemented for m!=2";
    }
  }

  public float[] applyInverse(float[][] c, float[] x) {
    float[] y = new float[x.length];
    applyInverse(c,x,y);
    return y;
  }

  public void applyInverse(float[][] c, float[] x, float[] y) {
    int m = c.length;
    if (m==2) {
      applyInverse(c[0],c[1],x,y);
    } else {
      assert false:"implemented for m!=2";
    }
  }

  public float[] applyNotch(double r, float[][] c, float[] x) {
    float[] y = new float[x.length];
    applyNotch(r,c,x,y);
    return y;
  }

  public void applyNotch(double r, float[][] c, float[] x, float[] y) {
    Check.argument(c.length==2,"c.length==2");
    int n = x.length;
    float s = (float)r;
    float ss = s*s;
    float[] t = new float[n];
    applyForward(c,x,t);
    c = copy(c);
    float[] c1 = c[0];
    float[] c2 = c[1];
    for (int i=0; i<n; ++i) {
      float c1i = c1[i];
      float c2i = c2[i];
      c2[i] = ss*c2i;
      c1[i] = s*c1i*(1.0f-c2i)/(1.0f-c2[i]);
    }
    applyInverse(c,t,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  /**
   * Compute coefficients c, 1st stage, before smoothing cn and cd.
   */
  private static void computeC1(
    int k, float[] b, float[] f, 
    float[] cn, float[] cd)
  {
    int n = b.length;
    for (int i=n-1; i>k; --i) {
      float bi = b[i-1];
      float fi = f[i  ];
      cn[i] = 2.0f*bi*fi;
      cd[i] = bi*bi+fi*fi;
    }
    for (int i=k; i>=0; --i) {
      cn[i] = 0.0f;
      cd[i] = 0.0f;
    }
  }

  /**
   * Compute coefficients c, 2nd stage, after smoothing cn and cd.
   */
  private static void computeC2(
    int k, float[] cn, float[] cd,
    float[] c, float[] b, float[] f)
  {
    int n = c.length;
    for (int i=n-1; i>k; --i) {
      float ci = cn[i]/cd[i];
      //assert ci>=-1.0f:"ci>=-1.0f";
      //assert ci<= 1.0f:"ci<= 1.0f";
      float bi = b[i-1];
      float fi = f[i  ];
      c[i] = ci;
      f[i] = fi-ci*bi;
      b[i] = bi-ci*fi;
    }
    for (int i=k; i>=0; --i)
      c[i] = c[k+1];
  }

  /**
   * Returns peak frequencies (in cycles/sample) from c1 and c2.
   */
  private static float[] getPeakFrequencies(float[] c1, float[] c2) {
    int n = c1.length;
    float[] pf = new float[n];
    float o2pi = 0.5f/FLT_PI;
    for (int i=0; i<n; ++i) {
      float c1i = c1[i];
      float c2i = c2[i];
      float cnum = -c1i*(1.0f-c2i)*(1.0f-c2i);
      float cden = 4.0f*c2i;
      float pw = (abs(cnum)<=abs(cden))?acos(cnum/cden):0.0f;
      pf[i] = pw*o2pi;
    }
    return pf;
  }

  /**
   * 2-coefficient forward filter with inner loop unrolled.
   */
  private static void applyForward(
    float[] c1, float[] c2, float[] x, float[] y) 
  {
    int n = x.length;
    float b0 = 0.0f, b1 = 0.0f, b2 = 0.0f;
    float f0 = 0.0f, f1 = 0.0f, f2 = 0.0f;
    for (int i=0; i<n; ++i) {
      float c1i = c1[i];
      float c2i = c2[i];
      f0 = x[i];
      f1 = f0-c1i*b0;
      f2 = f1-c2i*b1;
      b2 = b1-c2i*f1;
      b1 = b0-c1i*f0;
      b0 = f0;
      y[i] = f2;
    }
  }

  /**
   * 2-coefficient inverse filter with inner loop unrolled.
   */
  private void applyInverse(
    float[] c1, float[] c2, float[] x, float[] y) 
  {
    int n = x.length;
    float b0 = 0.0f, b1 = 0.0f, b2 = 0.0f;
    float f0 = 0.0f, f1 = 0.0f, f2 = 0.0f;
    for (int i=0; i<n; ++i) {
      float c1i = c1[i];
      float c2i = c2[i];
      f2 = x[i];
      f1 = f2+c2i*b1;
      f0 = f1+c1i*b0;
      b2 = b1-c2i*f1;
      b1 = b0-c1i*f0;
      b0 = f0;
      y[i] = f0;
    }
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
    test1();
  }

  public static void test1() {
    int n = 201;
    float w = 2.0f*FLT_PI*0.2f;
    float[] x = sweep(n,0.1,0.3);
    LocalBurgFilter lbf = new LocalBurgFilter();
    float[][] c = lbf.getCoefficients(10.0,2,x);
    plot("x",x);
    plot("c1",c[0]);
    plot("c2",c[1]);
    float[] pf = lbf.getPeakFrequencies(c);
    plot("pf",pf,0.0,0.4);
    float[] yf = copy(x);
    float[] yi = copy(x);
    float[] yn = copy(x);
    lbf.applyForward(c,x,yf);
    lbf.applyInverse(c,yf,yi);
    lbf.applyNotch(0.9,c,x,yn);
    plot("yf",yf);
    plot("yi",yi);
    plot("yn",yn);
  }

  public static void plot(String title, float[] x) {
    plot(title,x,min(x),max(x));
  }

  public static void plot(String title, float[] x, double xmin, double xmax) {
    SimplePlot sp = SimplePlot.asPoints(x);
    sp.setTitle(title);
    sp.setVLimits(xmin,xmax);
  }
}
