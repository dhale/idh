/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package het;

import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

// for testing only
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;

/**
 * Recursive symmetric exponential smoothing filter. For low frequencies 
 * this filter approximates a Gaussian with specified half-widths sigma.
 * <p>
 * This filter is faster than a recursive Gaussian filter. However, unlike 
 * the Gaussian filter, the symmetric exponential filter is not isotropic 
 * in multiple dimensions.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.04.03
 */
public class RecursiveExponentialFilter {

  /**
   * The method used to extrapolate input samples when filtering.
   * Input samples are defined implicitly by extrapolation outside the 
   * domain for which samples are specified explicitly, with either zero 
   * or constant values. If constant, the extrapolated values are the 
   * first and last uniform sample values. The default is extrapolation 
   * with zeros.
   */
  public enum Extrapolation {
    ZERO,
    CONSTANT
  };

  public RecursiveExponentialFilter(double sigma) {
    this(sigma,sigma,sigma);
  }

  public RecursiveExponentialFilter(double sigma1, double sigma23) {
    this(sigma1,sigma23,sigma23);
  }

  public RecursiveExponentialFilter(
    double sigma1, double sigma2, double sigma3)
  {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma3;
    _a1 = (_sigma1<SIGMA_MIN)?0.0f:(_sigma1-SIGMA_MIN)/(_sigma1+SIGMA_MIN);
    _a2 = (_sigma2<SIGMA_MIN)?0.0f:(_sigma2-SIGMA_MIN)/(_sigma2+SIGMA_MIN);
    _a3 = (_sigma3<SIGMA_MIN)?0.0f:(_sigma3-SIGMA_MIN)/(_sigma3+SIGMA_MIN);
    _s1 = (1.0f-_a1)/(1.0f+_a1);
    _s2 = (1.0f-_a2)/(1.0f+_a2);
    _s3 = (1.0f-_a3)/(1.0f+_a3);
  }

  /**
   * Sets the method used to implicitly extrapolate input values.
   * @param e the extrapolation method.
   */
  public void setExtrapolation(Extrapolation e) {
    _ec = e==Extrapolation.CONSTANT;
  }

  public void apply(float[] x, float[] y) {
    apply1(x,y);
  }

  public void apply(float[][] x, float[][] y) {
    apply2(x,y);
    apply1(y,y);
  }

  public void apply(float[][][] x, float[][][] y) {
    apply3(x,y);
    apply2(y,y);
    apply1(y,y);
  }

  public void apply1(float[] x, float[] y) {
    if (!distinct(x,y))
      x = copy(x);
    smooth1(_ec,_a1,_s1,x,y);
  }

  public void apply1(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = null;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] y2 = y[i2];
      if (!distinct(x2,y2)) {
        if (xt==null) 
          xt = new float[n1];
        copy(x2,xt);
        x2 = xt;
      }
      smooth1(_ec,_a1,_s1,x2,y2);
    }
  }

  public void apply2(float[][] x, float[][] y) {
    if (!distinct(x,y))
      x = copy(x);
    smooth2(_ec,_a2,_s2,x,y);
  }

  public void apply1(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] xtt = null;
        for (int i2=0; i2<n2; ++i2) {
          float[] x32 = x[i3][i2];
          float[] y32 = y[i3][i2];
          if (!distinct(x32,y32)) {
            if (xtt==null) 
              xtt = new float[n1];
            copy(x32,xtt);
            x32 = xtt;
          }
          smooth1(_ec,_a1,_s1,x32,y32);
        }
      }
    });
  }

  public void apply2(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[][] xt = null;
        float[][] x3 = x[i3];
        float[][] y3 = y[i3];
        if (!distinct(x3,y3)) {
          if (xt==null) 
            xt = new float[n2][n1];
          copy(x3,xt);
          x3 = xt;
        }
        smooth2(_ec,_a2,_s2,x3,y3);
      }
    });
  }

  public void apply3(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[][] x2 = new float[n3][];
        float[][] y2 = new float[n3][];
        float[][] xt = null;
        for (int i3=0; i3<n3; ++i3) {
          x2[i3] = x[i3][i2];
          y2[i3] = y[i3][i2];
        }
        if (!distinct(x2,y2)) {
          if (xt==null)
            xt = new float[n3][n1];
          copy(x2,xt);
          x2 = xt;
        }
        smooth2(_ec,_a3,_s3,x2,y2);
      }
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float SIGMA_MIN = sqrt(2.0f)/2.0f;

  private float _sigma1,_a1,_s1;
  private float _sigma2,_a2,_s2;
  private float _sigma3,_a3,_s3;
  private boolean _ec;

  private static boolean distinct(float[] x, float[] y) {
    return x!=y;
  }

  private static boolean distinct(float[][] x, float[][] y) {
    if (x==y)
      return false;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      if (!distinct(x[i2],y[i2]))
        return false;
    return true;
  }

  private static boolean distinct(float[][][] x, float[][][] y) {
    if (x==y)
      return false;
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3)
      if (!distinct(x[i3],y[i3]))
        return false;
    return true;
  }

  private static void smooth1(
    boolean ec, float a, float s, float[] x, float[] y)
  {
    if (a==0.0f) {
      copy(x,y);
    } else {
      int n1 = x.length;
      y[0] = ec?x[0]*s/(1.0f-a):x[0]*s;
      for (int i1=1; i1<n1; ++i1)
        y[i1] = a*y[i1-1]+s*x[i1];
      float yip1 = ec?x[n1-1]*s*a/(1.0f-a):0.0f;
      y[n1-1] += yip1;
      for (int i1=n1-2; i1>=0; --i1) {
        float yi = a*(yip1+s*x[i1+1]);
        y[i1] += yi;
        yip1 = yi;
      }
    }
  }

  private static void smooth2(
    boolean ec, float a, float s, float[][] x, float[][] y)
  {
    if (a==0.0f) {
      copy(x,y);
    } else {
      int n1 = x[0].length;
      int n2 = x.length;
      float[] xi,xp,yi,yp,ym;
      xi = x[0];
      yi = y[0];
      float s0 = ec?s/(1.0f-a):s;
      for (int i1=0; i1<n1; ++i1)
        yi[i1] = s0*xi[i1];
      for (int i2=1; i2<n2; ++i2) {
        xi = x[i2];
        yi = y[i2];
        ym = y[i2-1];
        for (int i1=0; i1<n1; ++i1)
          yi[i1] = a*ym[i1]+s*xi[i1];
      }
      yp = new float[n1];
      float sn = ec?s*a/(1.0f-a):0.0f;
      xp = x[n2-1];
      yi = y[n2-1];
      for (int i1=0; i1<n1; ++i1) {
        float yii = sn*xp[i1];
        yp[i1] = yii;
        yi[i1] += yii;
      }
      for (int i2=n2-2; i2>=0; --i2) {
        xp = x[i2+1];
        yi = y[i2];
        for (int i1=0; i1<n1; ++i1) {
          float yii = a*(yp[i1]+s*xp[i1]);
          yp[i1] = yii;
          yi[i1] += yii;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    test1();
    test2();
    test3();
  }

  public static void test1() {
    int n1 = 201;
    float[] x = new float[n1];
    float[] y = new float[n1];
    x[0] = x[n1/2] = x[n1-1] = 1.0f;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(20.0);
    copy(x,y);
    ref.apply(y,y);
    plot("y",y);
    fill(1.0f,x);
    ref.setExtrapolation(RecursiveExponentialFilter.Extrapolation.CONSTANT);
    ref.apply(x,y);
    System.out.println("x constant, y min ="+min(y)+" max ="+max(y));
    //plot("y",y);
  }

  public static void test2() {
    int n1 = 203;
    int n2 = 201;
    float[][] x = new float[n2][n1];
    float[][] y = new float[n2][n1];
    x[   0][   0] = 1.0f;
    x[   0][n1-1] = 1.0f;
    x[n2-1][   0] = 1.0f;
    x[n2-1][n1-1] = 1.0f;
    x[n2/2][n1/2] = 1.0f;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(30.0);
    copy(x,y);
    ref.apply(y,y);
    plot("y",y);
    fill(1.0f,x);
    ref.setExtrapolation(RecursiveExponentialFilter.Extrapolation.CONSTANT);
    ref.apply(x,y);
    System.out.println("x constant, y min ="+min(y)+" max ="+max(y));
    //plot("y",y);
  }

  public static void test3() {
    int n1 = 203;
    int n2 = 201;
    int n3 = 199;
    float[][][] x = new float[n3][n2][n1];
    float[][][] y = new float[n3][n2][n1];
    x[   0][   0][   0] = 1.0f;
    x[   0][   0][n1-1] = 1.0f;
    x[   0][n2-1][   0] = 1.0f;
    x[   0][n2-1][n1-1] = 1.0f;
    x[n3-1][   0][   0] = 1.0f;
    x[n3-1][   0][n1-1] = 1.0f;
    x[n3-1][n2-1][   0] = 1.0f;
    x[n3-1][n2-1][n1-1] = 1.0f;
    x[n3/2][n2/2][n1/2] = 1.0f;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(40.0);
    copy(x,y);
    ref.apply(y,y);
    plot("y",copy(y));
    fill(1.0f,x);
    ref.setExtrapolation(RecursiveExponentialFilter.Extrapolation.CONSTANT);
    ref.apply(x,y);
    System.out.println("x constant, y min ="+min(y)+" max ="+max(y));
    //plot("y",y);
  }

  public static void plot(String title, float[] x) {
    SimplePlot sp = SimplePlot.asPoints(x);
    sp.setTitle(title);
  }

  public static void plot(String title, float[][] x) {
    SimplePlot sp = new SimplePlot();
    PixelsView pv = sp.addPixels(x);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setColorModel(ColorMap.JET);
    sp.addColorBar();
    sp.setTitle(title);
  }

  public static void plot(String title, float[][][] x) {
    SimpleFrame sf = SimpleFrame.asImagePanels(x);
    sf.setSize(800,800);
  }
}
