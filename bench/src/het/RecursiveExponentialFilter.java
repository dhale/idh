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
 * Recursive symmetric exponential smoothing filter. Except perhaps near
 * the edges of input and output arrays, the impulse response of this 
 * two-sided filter is symmetric and decays exponentially from its peak 
 * value at zero lag.
 * <p>
 * Like the Gaussian filter, the impulse response of the exponential
 * filter is nowhere zero. The half-width sigma for the exponential 
 * filter is here defined so that, for low frequencies, the frequency 
 * response of the exponential filter approximates that for a Gaussian 
 * filter with the same specified half-width sigma. Specifically, the
 * value, slope and curvature of the frequency responses will be the
 * same for the exponential and Gaussian filters if the same half-widths
 * are specified.
 * <p>
 * This smoothing filter is faster than a recursive Gaussian filter. 
 * This filter also provides a variety of boundary conditions that
 * can be used to control the filtering of samples near the edges
 * of arrays. For most (but not all) of these boundary conditions, 
 * this filter is symmetric and positive-definite (SPD). This means, 
 * for example, that it can be used as a preconditioner in 
 * conjugate-gradient solutions of SPD systems of equations.
 * <p>
 * Multidimensional filters are applied as a cascade of one-dimensional
 * filters applied for each dimension of multidimensional arrays. In 
 * contrast to the Gaussian filter, this cascade for the exponential 
 * filter does not have isotropic impulse or frequency responses.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.08.11
 */
public class RecursiveExponentialFilter {

  /**
   * Boundary condition used at edges of either input or output samples.
   * All except the input-zero-slope condition yield a filter that is
   * symmetric positive-definite (SPD).
   * <p>
   * The default boundary condition is output-zero-slope.
   */
  public enum Edges {
    /**
     * Extrapolate input samples beyond edges with zero values.
     * A filter with this boundary condition is SPD.
     */
    INPUT_ZERO_VALUE,
    /**
     * Extrapolate input samples beyond edges with values at the edges.
     * Extrapolated values will be constant, so that the slope is zero. 
     * <em>A filter with this boundary condition is not SPD.</em>
     */
    INPUT_ZERO_SLOPE,
    /**
     * Constrain output values beyond edges to be zero.
     * Output samples near the edges will have nearly zero value.
     * A filter with this boundary condition is SPD.
     */
    OUTPUT_ZERO_VALUE,
    /**
     * Constrain output values beyond edges to be constant.
     * Output samples near the edges will have nearly zero slope.
     * A filter with this boundary condition (the default) is SPD.
     */
    OUTPUT_ZERO_SLOPE
  };

  /**
   * Constructs a filter with specified half-width.
   * The same half-width is used when applying the filter for all 
   * dimensions of multidimensional arrays.
   * @param sigma filter half-width.
   */
  public RecursiveExponentialFilter(double sigma) {
    this(sigma,sigma,sigma);
  }

  /**
   * Constructs a filter with specified half-widths.
   * @param sigma1 filter half-width for the 1st dimension.
   * @param sigma23 filter half-width for 2nd and 3rd dimensions. 
   */
  public RecursiveExponentialFilter(double sigma1, double sigma23) {
    this(sigma1,sigma23,sigma23);
  }

  /**
   * Constructs a filter with specified half-widths.
   * @param sigma1 filter half-width for the 1st dimension.
   * @param sigma2 filter half-width for the 2nd dimension.
   * @param sigma3 filter half-width for the 3rd dimension.
   */
  public RecursiveExponentialFilter(
    double sigma1, double sigma2, double sigma3)
  {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma3;
    _a1 = aFromSigma(sigma1);
    _a2 = aFromSigma(sigma2);
    _a3 = aFromSigma(sigma3);
  }

  /**
   * Sets the boundary condition used for samples beyond edges.
   * @param edges the boundary condition.
   */
  public void setEdges(Edges edges) {
    _ei = edges==Edges.INPUT_ZERO_VALUE || edges==Edges.INPUT_ZERO_SLOPE;
    _zs = edges==Edges.INPUT_ZERO_SLOPE || edges==Edges.OUTPUT_ZERO_SLOPE;
  }

  /**
   * Applies this filter.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float[] x, float[] y) {
    apply1(x,y);
  }

  /**
   * Applies this filter along all array dimensions.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float[][] x, float[][] y) {
    apply2(x,y);
    apply1(y,y);
  }

  /**
   * Applies this filter along all array dimensions.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float[][][] x, float[][][] y) {
    apply3(x,y);
    apply2(y,y);
    apply1(y,y);
  }

  /**
   * Applies this filter along the 1st (only) array dimension.
   * @param x input array.
   * @param y output array.
   */
  public void apply1(float[] x, float[] y) {
    if (!distinct(x,y))
      x = copy(x);
    smooth1(_ei,_zs,_a1,x,y);
  }

  /**
   * Applies this filter along the 1st array dimension.
   * @param x input array.
   * @param y output array.
   */
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
      smooth1(_ei,_zs,_a1,x2,y2);
    }
  }

  /**
   * Applies this filter along the 2nd array dimension.
   * @param x input array.
   * @param y output array.
   */
  public void apply2(float[][] x, float[][] y) {
    if (!distinct(x,y))
      x = copy(x);
    smooth2(_ei,_zs,_a2,x,y);
  }

  /**
   * Applies this filter along the 1st array dimension.
   * @param x input array.
   * @param y output array.
   */
  public void apply1(float[][][] x, float[][][] y) {
    final float[][][] xx = x;
    final float[][][] yy = y;
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] xtt = null;
        for (int i2=0; i2<n2; ++i2) {
          float[] x32 = xx[i3][i2];
          float[] y32 = yy[i3][i2];
          if (!distinct(x32,y32)) {
            if (xtt==null) 
              xtt = new float[n1];
            copy(x32,xtt);
            x32 = xtt;
          }
          smooth1(_ei,_zs,_a1,x32,y32);
        }
      }
    });
  }

  /**
   * Applies this filter along the 2nd array dimension.
   * @param x input array.
   * @param y output array.
   */
  public void apply2(float[][][] x, float[][][] y) {
    final float[][][] xx = x;
    final float[][][] yy = y;
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[][] xt = null;
        float[][] x3 = xx[i3];
        float[][] y3 = yy[i3];
        if (!distinct(x3,y3)) {
          if (xt==null) 
            xt = new float[n2][n1];
          copy(x3,xt);
          x3 = xt;
        }
        smooth2(_ei,_zs,_a2,x3,y3);
      }
    });
  }

  /**
   * Applies this filter along the 3rd array dimension.
   * @param x input array.
   * @param y output array.
   */
  public void apply3(float[][][] x, float[][][] y) {
    final float[][][] xx = x;
    final float[][][] yy = y;
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[][] x2 = new float[n3][];
        float[][] y2 = new float[n3][];
        float[][] xt = null;
        for (int i3=0; i3<n3; ++i3) {
          x2[i3] = xx[i3][i2];
          y2[i3] = yy[i3][i2];
        }
        if (!distinct(x2,y2)) {
          if (xt==null)
            xt = new float[n3][n1];
          copy(x2,xt);
          x2 = xt;
        }
        smooth2(_ei,_zs,_a3,x2,y2);
      }
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static float aFromSigma(double sigma) {
    if (sigma<=0.0f)
      return 0.0f;
    double ss = sigma*sigma;
    return (float)((1.0+ss-sqrt(1.0+2.0*ss))/ss);
  }

  private float _sigma1,_a1;
  private float _sigma2,_a2;
  private float _sigma3,_a3;
  private boolean _ei = false; // true, iff b.c. specified for input edges
  private boolean _zs = true; // true, iff zero slope boundary conditions

  private static void smooth1(
    boolean ei, boolean zs, float a, float[] x, float[] y) 
  {
    if (a==0.0f) {
      copy(x,y);
    } else if (ei) {
      smooth1Ei(zs,a,x,y);
    } else {
      smooth1Eo(zs,a,x,y);
    }
  }

  private static void smooth2(
    boolean ei, boolean zs, float a, float[][] x, float[][] y) 
  {
    if (a==0.0f) {
      copy(x,y);
    } else if (ei) {
      smooth2Ei(zs,a,x,y);
    } else {
      smooth2Eo(zs,a,x,y);
    }
  }

  private static void smooth1Ei(
    boolean zs, float a, float[] x, float[] y)
  {
    int n1 = x.length;
    float b = 1.0f-a;

    // forward
    y[0] = zs?x[0]:b*x[0];
    for (int i1=1; i1<n1; ++i1)
      y[i1] = a*y[i1-1]+b*x[i1];

    // reverse
    y[n1-1] = zs?(y[n1-1]+a*x[n1-1])/(1.0f+a):y[n1-1]/(1.0f+a);
    for (int i1=n1-2; i1>=0; --i1)
      y[i1] = a*y[i1+1]+b*y[i1];
  }

  private static void smooth2Ei(
    boolean zs, float a, float[][] x, float[][] y)
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float b = 1.0f-a;
    float[] xi,yi,yp,ym;
    float sx,sy;

    // forward
    xi = x[0];
    yi = y[0];
    sx = zs?1.0f:b;
    for (int i1=0; i1<n1; ++i1)
      yi[i1] = sx*xi[i1];
    for (int i2=1; i2<n2; ++i2) {
      xi = x[i2];
      yi = y[i2];
      ym = y[i2-1];
      for (int i1=0; i1<n1; ++i1)
        yi[i1] = a*ym[i1]+b*xi[i1];
    }

    // reverse
    sx = zs?a/(1.0f+a):0.0f;
    sy = 1.0f/(1.0f+a);
    xi = x[n2-1];
    yi = y[n2-1];
    for (int i1=0; i1<n1; ++i1)
      yi[i1] = sy*yi[i1]+sx*xi[i1];
    for (int i2=n2-2; i2>=0; --i2) {
      xi = x[i2];
      yi = y[i2];
      yp = y[i2+1];
      for (int i1=0; i1<n1; ++i1)
        yi[i1] = a*yp[i1]+b*yi[i1];
    }
  }

  private static void smooth1Eo(
    boolean zs, float a, float[] x, float[] y)
  {
    int n1 = x.length;
    float aa = a*a;
    float ss = zs?1.0f-a:1.0f;
    float gg = zs?aa-a:aa;
    float c = (1.0f-aa-ss)/ss;
    float d = 1.0f/(1.0f-aa+gg*(1.0f+c*pow(aa,n1-1)));
    float e = (1.0f-a)*(1.0f-a)*FLT_EPSILON/4.0f;

    // copy scaled input to output
    mul((1.0f-a)*(1.0f-a),x,y);

    // reversed triangular factorization
    int k1 = min((int)ceil(log(e)/log(a)),2*n1-2); // 1 <= k1 <= 2*n1-2
    float ynm1 = 0.0f;
    int m1 = k1-n1+1; // 2-n1 <= m1 <= n1-1
    for (int i1=m1; i1>0; --i1)
      ynm1 = a*ynm1+y[i1];
    ynm1 *= c;
    if (n1-k1<1)
      ynm1 = a*ynm1+(1.0f+c)*y[0];
    m1 = max(n1-k1,1); // 1 <= m1 <= n1-1
    for (int i1=m1; i1<n1; ++i1)
      ynm1 = a*ynm1+y[i1];
    ynm1 *= d;

    // reverse substitution
    y[n1-1] -= gg*ynm1;
    for (int i1=n1-2; i1>=0; --i1)
      y[i1] += a*y[i1+1];
    y[0] /= ss;

    // forward substitution
    for (int i1=1; i1<n1-1; ++i1)
      y[i1] += a*y[i1-1];
    y[n1-1] = ynm1;
  }

  private static void smooth2Eo(
    boolean zs, float a, float[][] x, float[][] y)
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float aa = a*a;
    float ss = zs?1.0f-a:1.0f;
    float gg = zs?aa-a:aa;
    float c = (1.0f-aa-ss)/ss;
    float d = 1.0f/(1.0f-aa+gg*(1.0f+c*pow(aa,n1-1)));
    float e = (1.0f-a)*(1.0f-a)*FLT_EPSILON/4.0f;

    // copy scaled input to output
    mul((1.0f-a)*(1.0f-a),x,y);

    // reversed triangular factorization
    int k2 = min((int)ceil(log(e)/log(a)),2*n2-2);
    float[] ynm1 = new float[n1];
    int m2 = k2-n2+1;
    for (int i2=m2; i2>0; --i2) {
      float[] yi = y[i2];
      for (int i1=0; i1<n1; ++i1)
        ynm1[i1] = a*ynm1[i1]+yi[i1];
    }
    for (int i1=0; i1<n1; ++i1)
      ynm1[i1] *= c;
    if (n2-k2<1) {
      for (int i1=0; i1<n1; ++i1)
        ynm1[i1] = a*ynm1[i1]+(1.0f+c)*y[0][i1];
    }
    m2 = max(n2-k2,1);
    for (int i2=m2; i2<n2; ++i2) {
      float[] yi = y[i2];
      for (int i1=0; i1<n1; ++i1)
        ynm1[i1] = a*ynm1[i1]+yi[i1];
    }
    for (int i1=0; i1<n1; ++i1)
      ynm1[i1] *= d;

    // reverse substitution
    for (int i1=0; i1<n1; ++i1)
      y[n2-1][i1] -= gg*ynm1[i1];
    for (int i2=n2-2; i2>=0; --i2) {
      float[] yi = y[i2];
      float[] yp = y[i2+1];
      for (int i1=0; i1<n1; ++i1)
        yi[i1] += a*yp[i1];
    }
    float oss = 1.0f/ss;
    for (int i1=0; i1<n1; ++i1)
      y[0][i1] *= oss;

    // forward substitution
    for (int i2=1; i2<n2-1; ++i2) {
      float[] yi = y[i2];
      float[] ym = y[i2-1];
      for (int i1=0; i1<n1; ++i1)
        yi[i1] += a*ym[i1];
    }
    for (int i1=0; i1<n1; ++i1)
      y[n2-1][i1] = ynm1[i1];
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    test1();
    //test2();
    //test3();
  }

  public static void test1() {
    int n1 = 201;
    float[] xi = makeImpulses(n1);
    float[] xs = makeSteps(n1);
    float[][] xAll = {xi,xs};
    RecursiveExponentialFilter.Edges[] edgesAll = {
      RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE,
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE,
      RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE,
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE
    };
    for (float[] x:xAll) {
      for (Edges edges:edgesAll) {
        test1ForEdges(edges,x);
      }
    }
  }
  private static float[] makeImpulses(int n1) {
    float[] x = new float[n1];
    x[0] = x[n1/2] = x[n1-1] = 1.0f;
    return x;
  }
  private static float[] makeSteps(int n1) {
    float[] x = fillfloat(1.0f,n1);
    for (int i1=0; i1<n1/2; ++i1)
      x[i1] = -1.0f;
    return x;
  }
  private static void test1ForEdges(Edges edges, float[] x) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(10.0);
    ref.setEdges(edges);
    float[] y = copy(x);
    ref.apply(y,y);
    System.out.println("x constant, y min="+min(y)+" max="+max(y));
    plot("y for "+edges,y);
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
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE);
    copy(x,y);
    ref.apply(y,y);
    plot("y",y);
    fill(1.0f,x);
    for (int i2=0; i2<n2; ++i2) {
      float s2 = (i2<n2/2)?-1.0f:1.0f;
      for (int i1=0; i1<n1; ++i1) {
        float s1 = (i1<n1/2)?-1.0f:1.0f;
          x[i2][i1] *= s1*s2;
      }
    }
    ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    ref.apply(x,y);
    System.out.println("x constant, y min="+min(y)+" max="+max(y));
    plot("y",y);
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
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE);
    copy(x,y);
    ref.apply(y,y);
    plot("y",copy(y));
    fill(1.0f,x);
    for (int i3=0; i3<n3; ++i3) {
      float s3 = (i3<n3/2)?-1.0f:1.0f;
      for (int i2=0; i2<n2; ++i2) {
        float s2 = (i2<n2/2)?-1.0f:1.0f;
        for (int i1=0; i1<n1; ++i1) {
          float s1 = (i1<n1/2)?-1.0f:1.0f;
            x[i3][i2][i1] *= s1*s2*s3;
        }
      }
    }
    ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    ref.apply(x,y);
    System.out.println("x constant, y min="+min(y)+" max="+max(y));
    plot("y",y);
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
