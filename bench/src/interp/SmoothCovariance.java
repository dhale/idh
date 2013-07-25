/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

//import static java.lang.Math.*;
import edu.mines.jtk.opt.BrentZeroFinder;
import edu.mines.jtk.util.Check;
import static interp.SpecialMath.*;

// Testing only.
import static edu.mines.jtk.util.ArrayMath.*;
import java.awt.*;
import javax.swing.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * A covariance function corresponding to efficient smoothing filters.
 * Parameters for these filters are computed to approximate the Matern
 * covariance function.
 * @author Dave Hale, Colorado School of Mines.
 * @version 2012.09.15
 */
public class SmoothCovariance implements Covariance {

  /**
   * Constructs the covariance for specified parameters.
   * The number of spatial dimensions is currently required to be 2.
   * @param sigma the sqrt(variance) for distance r = 0.
   * @param shape the shape; greater than 0, not greater than 3.
   * @param range the effective range; greater than 0.
   * @param ndim the number of spatial dimensions; must be 2.
   */
  public SmoothCovariance(
    double sigma, double shape, double range, int ndim) {
    Check.argument(shape>0.0,"shape > 0");
    Check.argument(shape<=3.0,"shape <= 3");
    Check.argument(range>0.0,"range > 0");
    Check.argument(ndim==2,"ndim = 2");
    _sigma = sigma;
    _shape = shape;
    _range = range;
    _ndim = ndim;
    if (ndim==2) {
      init2();
      //init2(1.0+_shape);
    }
    _ldk = new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D22);
    _lsf = new LocalSmoothingFilter(1.0e-6,10000,_ldk);
    //_lsf.setPreconditioner(true);
    _stf = new SymmetricTridiagonalFilter(0.76,0.52,0.76,0.24);
    //_stf = new SymmetricTridiagonalFilter(0.80,0.60,0.80,0.20);
  }

  /**
   * Gets the shape for this covariance function.
   * @return the shape.
   */
  public double getShape() {
    return _shape;
  }

  /**
   * Gets the sigma for this covariance function.
   * @return the sigma.
   */
  public double getSigma() {
    return _sigma;
  }

  /**
   * Gets the range for this covariance function.
   * @return the range.
   */
  public double getRange() {
    return _range;
  }

  /**
   * Gets the number of dimensions for this covariance function.
   * @return the number of dimensions.
   */
  public double getDimensions() {
    return _ndim;
  }

  public double evaluate(double r) {
    return _li.interpolate(r);
  }

  /**
   * Applies this covariance function via 2D tensor-guided smoothing.
   * <em>Sampling intervals delta for s1 and s2 must be equal.</em>
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param tensors the tensor field.
   * @param q the input/output image.
   */
  public void apply(Sampling s1, Sampling s2, Tensors2 tensors, float[][] q) {
    Check.argument(s1.getDelta()==s2.getDelta(),"s1 and s2 have same delta");
    int n1 = q[0].length;
    int n2 = q.length;
    float[][] s = q;
    float[][] t = new float[n2][n1];
    float dx = (float)s1.getDelta();
    float akk = _ascl*_kscl*_kscl/dx/dx;
    float bkk = _bscl*_kscl*_kscl/dx/dx;
    float ckk = _cscl*_kscl*_kscl/dx/dx;
    cscale(ckk,tensors,q);
    _stf.apply1(q,s);
    _stf.apply2(s,s);
    _lsf.apply(tensors,bkk,s,t);
    for (int ifac=0; ifac<_nfac; ++ifac) {
      float[][] st = s; s = t; t = st;
      _lsf.apply(tensors,akk,s,t);
    }
    _lsf.apply(tensors,akk,s,t);
    _stf.apply2(t,q);
    _stf.apply1(q,q);
    cscale(ckk,tensors,q);
  }

  public void applyInverse(
    Sampling s1, Sampling s2, Tensors2 tensors, float[][] q) 
  {
    Check.argument(s1.getDelta()==s2.getDelta(),"s1 and s2 have same delta");
    int n1 = q[0].length;
    int n2 = q.length;
    float[][] s = q;
    float[][] t = new float[n2][n1];
    float dx = (float)s1.getDelta();
    float akk = _ascl*_kscl*_kscl/dx/dx;
    float bkk = _bscl*_kscl*_kscl/dx/dx;
    float ckk = _cscl*_kscl*_kscl/dx/dx;
    cscaleInverse(ckk,tensors,q);
    _stf.applyInverse1(q,t);
    _stf.applyInverse2(t,t);
    for (int ifac=0; ifac<_nfac; ++ifac) {
      float[][] st = s; s = t; t = st;
      copy(s,t); _ldk.apply(tensors,akk,s,t);
    }
    copy(t,s); _ldk.apply(tensors,bkk,t,s);
    _stf.applyInverse2(s,q);
    _stf.applyInverse1(q,q);
    cscaleInverse(ckk,tensors,q);
  }

  public float[][] apply(
    Sampling s1, Sampling s2, Tensors2 tensors,
    float[] x1, float[] x2, float[] z) 
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n = z.length;
    float[][] q = new float[n2][n1];
    for (int i=0; i<n; ++i) {
      int i1 = s1.indexOfNearest(x1[i]);
      int i2 = s2.indexOfNearest(x2[i]);
      q[i2][i1] = z[i];
    }
    apply(s1,s2,tensors,q);
    return q;
  }

  /**
   * Applies half of this covariance function via 2D tensor-guided smoothing.
   * <em>This method can be used for only shape = 1.</em>
   * <em>Sampling intervals delta for s1 and s2 must be equal.</em>
   * If the input image contains independently and identically distributed
   * Gaussian random numbers sampled from a standard normal distribution, then
   * the output image will exhibit spatial correlation consistent with this
   * covariance function.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param tensors the tensor field.
   * @param q the input/output image.
   */
  public void applyHalf(
    Sampling s1, Sampling s2, Tensors2 tensors, float[][] q) 
  {
    Check.state(_shape==1.0,"shape = 1");
    Check.argument(s1.getDelta()==s2.getDelta(),"s1 and s2 have same delta");
    int n1 = q[0].length;
    int n2 = q.length;
    float[][] t = new float[n2][n1];
    float dx = (float)s1.getDelta();
    float akk = _ascl*_kscl*_kscl/dx/dx;
    float ckk = _cscl*_kscl*_kscl/dx/dx;
    cscale(ckk,tensors,q);
    _stf.apply1(q,t);
    _stf.apply2(t,t);
    _lsf.apply(tensors,akk,t,q);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _shape;
  private double _sigma;
  private double _range;
  private int _ndim;
  private int _nfac; // number of smoothings using ascl
  private float _ascl,_bscl,_cscl,_kscl; // scale factors
  private LinearInterpolator _li;
  private LocalSmoothingFilter _lsf;
  private LocalDiffusionKernel _ldk;
  private SymmetricTridiagonalFilter _stf;

  /**
   * Applies scaling by factors that depend on 2D tensor determinants.
   * Is applied twice, both before and after smoothing filters, which
   * ensures that this covariance operator is symmetric.
   */
  private void cscale(float ckk, Tensors2 tensors, float[][] t) {
    float[] d = new float[3];
    int n1 = t[0].length;
    int n2 = t.length;
    float rckk = sqrt(ckk);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        tensors.getTensor(i1,i2,d);
        float d11 = d[0], d12 = d[1], d22 = d[2];
        float det = d11*d22-d12*d12;
        t[i2][i1] *= rckk*sqrt(sqrt(det));
      }
    }
  }

  /**
   * Applies inverse scaling by factors that depend on 2D tensor determinants.
   */
  private void cscaleInverse(float ckk, Tensors2 tensors, float[][] t) {
    float[] d = new float[3];
    int n1 = t[0].length;
    int n2 = t.length;
    float rckk = 1.0f/sqrt(ckk);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        tensors.getTensor(i1,i2,d);
        float d11 = d[0], d12 = d[1], d22 = d[2];
        float det = 1.0f/(d11*d22-d12*d12);
        t[i2][i1] *= rckk*sqrt(sqrt(det));
      }
    }
  }

  /**
   * Initializes this covariance function for 2D.
   */
  private void init2() {
    MaternMatch2 mm2 = new MaternMatch2(_shape);
    _nfac = mm2.n;
    _ascl = (float)mm2.a;
    _bscl = (float)mm2.b;
    _cscl = (float)(mm2.c*_sigma*_sigma);
    _kscl = (float)(0.5*_range/sqrt(_shape));
    trace("nfac="+_nfac+" ascl="+_ascl+" bscl="+_bscl+" cscl="+_cscl);
    double rmax =  16.0*_range; // c(r)/c(0) < 0.001, for r > rmax
    int nr = 10001;
    double dr = rmax/(nr-1);
    float[] cr = new float[nr];
    for (int ir=0; ir<nr; ++ir) {
      double rn = ir*dr/_kscl; // normalized distance
      cr[ir] = (float)(_sigma*_sigma*eval2(_nfac,_ascl,_bscl,rn));
    }
    _li = new LinearInterpolator();
    _li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    _li.setUniform(nr,dr,0.0,cr);
  }

  /**
   * Finds parameters a and b for a 2D smooth covariance function that
   * approximates a Matern covariance function with specified shape.
   * Chooses the two parameters a and b to match the 90% and 10% values
   * of the Matern function. This choice is somewhat arbitrary. However,
   * the 90% value helps to match the Matern shape near the origin, and
   * the 10% value helps to match the Matern range far from the origin.
   */
  private static class MaternMatch2 {
    int n; // number of primary smoothings
    double a; // parameter for all primary smoothings
    double b; // parameter for one secondary smoothing
    double c; // scale parameter to ensure that c(0) = 1
    MaternMatch2(double shape) {
      Check.argument(0.0<shape && shape<=3.0,"0 < shape <= 3");
      double e = 1.0+shape;
      n = (int)e;
      a = 1.000000;
      b = 0.000000;
      c = 4.0*PI*shape;

      // If a secondary smoothing is required (so that b>0), ...
      if (e>n) {
        a = AMIN;
        b = BMIN;

        // Compute Matern covariance function and distances r to be 
        // matched. Scale the effective range in the Matern covariance 
        // function so that we can perform the matching without worrying 
        // about the scale factor. Matching two values of the Matern
        // covariance gives us two (nonlinear) equations to be solved 
        // for the two parameters a and b.
        MaternCovariance mc = new MaternCovariance(1.0,shape,2.0*sqrt(shape));
        double c1 = 0.90;
        double c2 = 0.10;
        double r1 = r(mc,c1);
        double r2 = r(mc,c2);

        // Perturbation of a and b used to estimate derivatives.
        double d = sqrt(DBL_EPSILON);

        // Now use Newton-Raphson iterations to find parameters 
        // a and b such that f(r1) = c1 and f(r2) = c2.
        double apb = a+b;
        while (apb>0.001) {
          double aold = a;
          double bold = b;

          // Evaluate smooth covariance values f(r1) and f(r2).
          double f1 = eval2(n,a,b,r1);
          double f2 = eval2(n,a,b,r2);

          // Numerical derivatives for Jacobian matrix.
          double f11 = (eval2(n,a+d,b,r1)-f1)/d;
          double f21 = (eval2(n,a+d,b,r2)-f2)/d;
          double f12 = (eval2(n,a,b+d,r1)-f1)/d;
          double f22 = (eval2(n,a,b+d,r2)-f2)/d;

          // Inverse of Jacobian matrix.
          double det = f11*f22-f21*f12;
          double a11 =  f22/det;
          double a21 = -f21/det;
          double a12 = -f12/det;
          double a22 =  f11/det;

          // Newton-Raphson update.
          a -= a11*(f1-c1)+a12*(f2-c2);
          b -= a21*(f1-c1)+a22*(f2-c2);

          // Ensure that a and b stay in bounds.
          a = max(AMIN,a);
          b = max(BMIN,min(BMAX,b));

          // Converged when this sum is small.
          apb = abs(a-aold)+abs(b-bold);
        }
        if (n==1) {
          c = 4.0*PI*(a-b)/log(a/b);
        } else if (n==2) {
          c = 4.0*PI*(a-b)*(a-b)/(a-b-b*log(a/b));
        } else if (n==3) {
          c = 8.0*PI*(a-b)*(a-b)*(a-b)/((a-3.0*b)*(a-b)+2.0*b*b*log(a/b));
        }
      }
    }
    private static double AMIN = 1.000000;
    private static double BMIN = 0.000001;
    private static double BMAX = 0.999999;
    private double r(final MaternCovariance mc, final double c) {
      BrentZeroFinder bzf = new BrentZeroFinder(
        new BrentZeroFinder.Function() {
          public double evaluate(double r) {
            return mc.evaluate(r)-c;
          }
        });
      return bzf.findZero(0.0,8.0,0.001);
    }
  }

  /**
   * Evaluates smooth covariance for specified shape parameters.
   * The range (not the so-called effective range) is one. The 
   * variance sigma is also one, so that the covariance is one
   * for distance zero.
   * This function is too costly to use in our evaluate method
   * above, but is useful for computing shape parameters and 
   * building a table for linear interpolation.
   * @param n number of primary smoothings.
   * @param a parameter for all primary smoothings.
   * @param b parameter for one secondary smoothing.
   */
  private static double eval2(int n, double a, double b, double r) {
    if (r==0.0) {
      return 1.0;
    } else if (a==1.0 && b==0.0) {
      if (n==2) {
        return r*k1(r);
      } else if (n==3) {
        return 0.5*r*r*kn(2,r);
      } else { // n==4
        return 0.125*r*r*r*kn(3,r);
      }
    } else if (n==1) {
      return 2.0*(k0(r/sqrt(a))-k0(r/sqrt(b)))/log(a/b);
    } else if (n==2) {
      return ((a-b)*r*k1(r/sqrt(a))/sqrt(a) -
              2.0*b*(k0(r/sqrt(a)) -
              k0(r/sqrt(b)))) /
             (a-b+b*log(b/a));
    } else { // n==3
      return (2.0*sqrt(a)*(a-3.0*b)*(a-b)*r*k1(r/sqrt(a)) -
              8.0*a*b*b*k0(r/sqrt(b)) +
              (8.0*a*b*b+(a-b)*(a-b)*r*r)*k0(r/sqrt(a))) /
             (2.0*a*((a-3.0*b)*(a-b)+2.0*b*b*log(a/b)));
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // test

  public void testSpd(int n1, int n2, Tensors2 t) {
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    float[][] x = sub(randfloat(n1,n2),0.5f);
    float[][] y = sub(randfloat(n1,n2),0.5f);
    float[][] ax = copy(x);
    float[][] ay = copy(y);
    //apply(s1,s2,t,ax);
    //apply(s1,s2,t,ay);
    _lsf.applySmoothS(x,ax);
    _lsf.applySmoothS(y,ay);
    float xay = sum(mul(x,ay));
    float yax = sum(mul(y,ax));
    float xax = sum(mul(x,ax));
    float yay = sum(mul(y,ay));
    System.out.println("xax="+xax+" yay="+yay);
    System.out.println("xay="+xay+" yax="+yax);
  }

  public static void tabulateMaternMatch2() {
    int nshape = 30;
    double dshape = 0.1;
    double fshape = dshape;
    Sampling sshape = new Sampling(nshape,dshape,fshape);
    for (int ishape=0; ishape<nshape; ++ishape) {
      double shape = sshape.getValue(ishape);
      if (shape>3.0) shape = 3.0;
      MaternMatch2 mm2 = new MaternMatch2(shape);
      int n = mm2.n;
      double a = mm2.a;
      double b = mm2.b;
      double c = mm2.c;
      System.out.printf("%3.1f %1d %8.6f %8.6f %9.6f%n",
        shape,mm2.n,mm2.a,mm2.b,mm2.c);
    }
  }

  public static void trace(String s) {
    System.out.println(s);
  }

  public static void plot(Sampling sx, float[] y, 
    double ymin, double ymax, String title) 
  {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    PointsView pv = sp.addPoints(sx,y);
    sp.setVLimits(ymin,ymax);
    sp.setSize(790,700);
    sp.setTitle(title);
  }

  public static void plot(Sampling sx, float[][] y, 
    double ymin, double ymax, String title) 
  {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    Color[] colors = {
      Color.LIGHT_GRAY,Color.GRAY,Color.DARK_GRAY,Color.BLACK
    };
    for (int iy=0; iy<y.length; ++iy) {
      PointsView pv = sp.addPoints(sx,y[iy]);
      pv.setLineColor(colors[iy%colors.length]);
    }
    sp.setVLimits(ymin,ymax);
    sp.setSize(790,700);
    sp.setTitle(title);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
    public void run() {
      //go();
      tabulateMaternMatch2();
    }});
  }
  public static void go() {
    int ndim = 2;
    double sigma = 1.0;
    double range = 1.0;
    double rmin = 0.0;
    double rmax = 4*range;
    double cmin = 0.0;
    double cmax = sigma*sigma+0.03;
    int nr = 1001;
    double dr = (rmax-rmin)/(nr-1);
    double fr = rmin;
    Sampling sr = new Sampling(nr,dr,fr);
    double[] shapes = {0.5,0.75,1.0,1.5};
    int nshape = shapes.length;
    float[][] cm = new float[nshape][nr];
    float[][] cs = new float[nshape][nr];
    for (int ishape=0; ishape<nshape; ++ishape) {
      double shape = shapes[ishape];
      System.out.println("shape="+shape);
      MaternCovariance mc = new MaternCovariance(sigma,shape,range);
      SmoothCovariance sc = new SmoothCovariance(sigma,shape,range,ndim);
      for (int ir=0; ir<nr; ++ir) {
        double r = sr.getValue(ir);
        cm[ishape][ir] = (float)mc.evaluate(r);
        cs[ishape][ir] = (float)sc.evaluate(r);
      }
      float[][] cms = new float[][]{cm[ishape],cs[ishape]};
      plot(sr,cms,cmin,cmax,"shape = "+shape);
    }
    plot(sr,cm,cmin,cmax,"Matern");
    plot(sr,cs,cmin,cmax,"Smooth");
  }
}
