/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

//import static java.lang.Math.*;
import edu.mines.jtk.util.Check;

// Testing only.
import static edu.mines.jtk.util.ArrayMath.*;
import java.awt.*;
import javax.swing.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * A covariance function corresponding to efficient smoothing filters.
 * @author Dave Hale, Colorado School of Mines.
 */
public class SmoothCovariance implements Covariance {

  /**
   * Constructs the covariance for specified parameters.
   * @param shape the shape; must be greater than zero.
   * @param sigma the sqrt(variance) for distance r = 0.
   * @param range the effective range.
   * @param ndim the number of spatial dimensions.
   */
  public SmoothCovariance(
    double shape, double sigma, double range, int ndim) {
    Check.argument(ndim==1 || ndim==2,"ndim = 1 or ndim = 2");
    _shape = shape;
    _sigma = sigma;
    _range = range;
    _ndim = ndim;
    _tensors = new IdentityTensors2();
    if (ndim==1) {
      makeInterpolator1();
    } else if (ndim==2) {
      makeInterpolator2();
    }
    _lsf = new LocalSmoothingFilter(1.0e-6,1000);
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

  /**
   * Gets the tensors used to compute distances and apply this covariance.
   * @return the tensors.
   */
  public Tensors2 getTensors() {
    return _tensors;
  }

  /**
   * Sets the tensors used to compute distances and apply this covariance.
   * The default is a constant isotropic tensor.
   * @param tensors the tensors; null for default tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = (tensors!=null)?tensors:new IdentityTensors2();
  }

  public double evaluate(double r) {
    if (r<0.0)
      r = -r;
    if (r>_rmax)
      return 0.0;
    return _li.interpolate(r);
  }

  public void apply(float[][] q) {
    int n1 = q[0].length;
    int n2 = q.length;
    float[] d = new float[3];
    float[][] t = new float[n2][n1];
    _lsf.applySmoothS(q,q);
    _lsf.apply(_tensors,_bscl*_ascl*_ascl,q,t);
    for (int ifac=0; ifac<_nfac; ++ifac) {
      float[][] qold = q; q = t; t = qold;
      _lsf.apply(_tensors,_ascl*_ascl,q,t);
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        _tensors.getTensor(i1,i2,d);
        float d11 = d[0], d12 = d[1], d22 = d[2];
        float det = d11*d22-d12*d12;
        t[i2][i1] *= _cscl*det;
      }
    }
    for (int ifac=0; ifac<_nfac; ++ifac) {
      float[][] qold = q; q = t; t = qold;
      _lsf.apply(_tensors,_ascl*_ascl,q,t);
    }
    _lsf.apply(_tensors,_bscl*_ascl*_ascl,t,q);
    _lsf.applySmoothS(q,q);
  }

  public float[][] apply(Sampling s1, Sampling s2, 
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
    apply(q);
    return q;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _shape;
  private double _sigma;
  private double _range;
  private int _ndim;
  private int _nfac;
  private float _ascl,_bscl,_cscl;
  private double _rmax;
  private Tensors2 _tensors;
  private LinearInterpolator _li;
  private LocalSmoothingFilter _lsf;

  private static class IdentityTensors2 implements Tensors2 {
    public void getTensor(int i1, int i2, float[] d) {
      d[0] = 1.0f;
      d[1] = 0.0f;
      d[2] = 1.0f;
    }
  }

  private void makeInterpolator1() {
    double cb = 4.0;
    double rmax = 20.0*_range;
    double e = (0.5+_shape)/2.0;
    _ascl = (float)(_range*0.5/sqrt(_shape));
    _nfac = (int)e;
    double kb = sqrt(pow(cb,1.0/e)-1.0)/_ascl;
    _bscl = (float)((cb/pow(1.0+_ascl*_ascl*kb*kb,_nfac)-1.0) /
                    (_ascl*_ascl*kb*kb));
    int nx = 1001;
    double dx = rmax/(nx-1);
    double fx = 0.0;
    Sampling sx = new Sampling(nx,dx,fx);
    Fft fft = new Fft(sx);
    Sampling sk = fft.getFrequencySampling1();
    int nk = sk.getCount();
    float[] g = new float[2*nk];
    for (int kk=0,kr=0; kk<nk; ++kk,kr+=2) {
      double k = 2.0*PI*sk.getValue(kk);
      double aks = _ascl*_ascl*k*k;
      double gkh = 1.0/((1.0+_bscl*aks)*pow(1.0+aks,_nfac));
      g[kr] = (float)(gkh*gkh);
    }
    float[] f = fft.applyInverse(g);
    _cscl = (float)(_sigma*_sigma)/f[0];
    mul(_cscl,f,f);
    nx = (nx-1)/2;
    float[] c = copy(nx,f);
    _li = new LinearInterpolator();
    _li.setUniform(nx,dx,fx,c);
    _rmax = fx+(nx-1)*dx;
  }

  private void makeInterpolator2() {
    double rmax = 20.0*_range;
    _ascl = (float)(_range*0.5/sqrt(_shape));
    double e = (1.0+_shape)/2.0;
    _nfac = (int)e;
    double cb = 16.0;
    double kb = sqrt(pow(cb,1.0/e)-1.0)/_ascl;
    _bscl = (float)((cb/pow(1.0+_ascl*_ascl*kb*kb,_nfac)-1.0) /
                    (_ascl*_ascl*kb*kb));
    System.out.println("e="+e+" nfac="+_nfac+" bscl="+_bscl);
    int nx = 1001;
    double dx = rmax/(nx-1);
    double fx = 0.0;
    Sampling sx = new Sampling(nx,dx,fx);
    Fft fft = new Fft(sx,sx);
    Sampling sk1 = fft.getFrequencySampling1();
    Sampling sk2 = fft.getFrequencySampling2();
    int nk1 = sk1.getCount();
    int nk2 = sk2.getCount();
    float[][] g = new float[nk2][2*nk1];
    for (int k2=0; k2<nk2; ++k2) {
      double k2k = 2.0*PI*sk2.getValue(k2);
      for (int k1=0,kr=0; k1<nk1; ++k1,kr+=2) {
        double k1k = 2.0*PI*sk1.getValue(k1);
        double aks = _ascl*_ascl*(k1k*k1k+k2k*k2k);
        double gkh = 1.0/((1.0+_bscl*aks)*pow(1.0+aks,_nfac));
        g[k2][kr] = (float)(gkh*gkh);
      }
    }
    float[][] f = fft.applyInverse(g);
    _cscl = (float)(_sigma*_sigma)/f[0][0];
    mul(_cscl,f[0],f[0]);
    nx = (nx-1)/2;
    float[] c = copy(nx,f[0]);
    _li = new LinearInterpolator();
    _li.setUniform(nx,dx,fx,c);
    _rmax = fx+(nx-1)*dx;
  }

  /**/
  ///////////////////////////////////////////////////////////////////////////
  // test

  public void testSpd(int n1, int n2) {
    float[][] x = sub(randfloat(n1,n2),0.5f);
    float[][] y = sub(randfloat(n1,n2),0.5f);
    float[][] ax = copy(x);
    float[][] ay = copy(y);
    apply(ax);
    apply(ay);
    float xay = sum(mul(x,ay));
    float yax = sum(mul(y,ax));
    float xax = sum(mul(x,ax));
    float yay = sum(mul(y,ay));
    System.out.println("xax="+xax+" yay="+yay);
    System.out.println("xay="+xay+" yax="+yax);
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
      Color.RED,Color.GREEN,Color.BLUE,
      Color.CYAN,Color.MAGENTA,Color.YELLOW,
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
      go();
    }});
  }

  public static void go() {
    double sigma = 2.0;
    double range = 2.0;
    double xmin = 0.0;
    double xmax = 2.0*range;
    double ymin = 0.0;
    double ymax = sigma*sigma+0.1;
    int nx = 1001;
    double dx = (xmax-xmin)/(nx-1);
    double fx = xmin;
    Sampling sx = new Sampling(nx,dx,fx);
    float[] y05 = new float[nx];
    float[] y10 = new float[nx];
    float[] y15 = new float[nx];
    float[] y20 = new float[nx];
    float[] ygg = new float[nx];
    float[] z05 = new float[nx];
    float[] z10 = new float[nx];
    float[] z15 = new float[nx];
    float[] z20 = new float[nx];
    SmoothCovariance c05 = new SmoothCovariance(0.5,sigma,range,2);
    SmoothCovariance c10 = new SmoothCovariance(1.0,sigma,range,2);
    SmoothCovariance c15 = new SmoothCovariance(1.5,sigma,range,2);
    SmoothCovariance c20 = new SmoothCovariance(2.0,sigma,range,2);
    MaternCovariance m05 = new MaternCovariance(0.5,sigma,range);
    MaternCovariance m10 = new MaternCovariance(1.0,sigma,range);
    MaternCovariance m15 = new MaternCovariance(1.5,sigma,range);
    MaternCovariance m20 = new MaternCovariance(2.0,sigma,range);
    for (int ix=0; ix<nx; ++ix) {
      float x = (float)sx.getValue(ix);
      y05[ix] = (float)c05.evaluate(x);
      y10[ix] = (float)c10.evaluate(x);
      y15[ix] = (float)c15.evaluate(x);
      y20[ix] = (float)c20.evaluate(x);
      ygg[ix] = (float)(sigma*sigma*exp(-x*x/(range*range)));
      z05[ix] = (float)m05.evaluate(x);
      z10[ix] = (float)m10.evaluate(x);
      z15[ix] = (float)m15.evaluate(x);
      z20[ix] = (float)m20.evaluate(x);
    }
    plot(sx,new float[][]{y05,z05},ymin,ymax,"c05(x) and m05(x)");
    plot(sx,new float[][]{y10,z10},ymin,ymax,"c10(x) and m10(x)");
    plot(sx,new float[][]{y15,z15},ymin,ymax,"c15(x) and m15(x)");
    plot(sx,new float[][]{y20,z20},ymin,ymax,"c20(x) and m20(x)");
    plot(sx,ygg,ymin,ymax,"Gaussian");
    plot(sx,new float[][]{y05,y10,y15,y20,ygg},ymin,ymax,"Smooth");
    plot(sx,new float[][]{z05,z10,z15,z20,ygg},ymin,ymax,"Matern");
  }
  /**/
}
