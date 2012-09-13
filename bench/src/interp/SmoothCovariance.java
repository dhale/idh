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
      init1();
    } else if (ndim==2) {
      init2();
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
    return _li.interpolate(r);
  }

  public void apply(float[][] q) {
    int n1 = q[0].length;
    int n2 = q.length;
    float[][] s = q;
    float[][] t = new float[n2][n1];
    cscale(q);
    _lsf.applySmoothS(q,s);
    _lsf.apply(_tensors,_bscl*_kscl*_kscl,s,t);
    for (int ifac=0; ifac<_nfac; ++ifac) {
      float[][] st = s; s = t; t = st;
      _lsf.apply(_tensors,_ascl*_kscl*_kscl,s,t);
    }
    _lsf.applySmoothS(t,q);
    cscale(q);
  }
  private void cscale(float[][] t) {
    float[] d = new float[3];
    int n1 = t[0].length;
    int n2 = t.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        _tensors.getTensor(i1,i2,d);
        float d11 = d[0], d12 = d[1], d22 = d[2];
        float det = d11*d22-d12*d12;
        t[i2][i1] *= _cscl*sqrt(sqrt(det));
      }
    }
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
  private float _ascl,_bscl,_cscl,_kscl;
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

  private void init1() {
    /*
    double cb = 16.0;
    double rmax = 20.0*_range;
    double e = 0.5+_shape;
    _kscl = (float)(_range*0.5/sqrt(_shape));
    _nfac = (int)e;
    double kb = sqrt(pow(cb,1.0/e)-1.0)/_kscl;
    _bscl = (float)((cb/pow(1.0+_kscl*_kscl*kb*kb,_nfac)-1.0) /
                    (_kscl*_kscl*kb*kb));
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
      double aks = _kscl*_kscl*k*k;
      double gkr = 1.0/((1.0+_bscl*aks)*pow(1.0+aks,_nfac));
      g[kr] = (float)gkr;
    }
    float[] f = fft.applyInverse(g);
    _cscl = (float)sqrt((_sigma*_sigma)/f[0]);
    mul(_cscl*_cscl,f,f);
    nx = (nx-1)/2;
    float[] c = copy(nx,f);
    _li = new LinearInterpolator();
    _li.setUniform(nx,dx,fx,c);
    */
  }

  /**
   * Initializes for 2D smooth covariance function. Computes various 
   * constants and constructs the linear interpolator used to evaluate 
   * this function.
   */
  private void init2() {

    // Factor by which to (approximately) scale wavenumbers k in 
    // our smoothing filter approximations of 1/(1+(kscl*k)^2).
    // The factor 0.5/sqrt(_shape) makes the specified range the
    // so-called "effective range."
    _kscl = (float)(0.5*_range/sqrt(_shape));

    // The exponent in the amplitude spectrum of the Matern 
    // covariance function, which is proportional to 1/(1+k^2)^e.
    double e = 1.0+_shape;

    // Integer number of factors 1/(1+k^2) in our approximation
    // 1/(1+b*k^2) * 1/(1+k^2)^nfac of the Matern covariance
    // function. Each factor will, in turn, be approximated by
    // application of a smoothing filter.
    _nfac = (int)e;

    // If the exponent e is not an integer, then we need one more
    // more smoothing filter to complete our approximation. This
    // filter has an amplitude spectrum that approximates 1/(1+b*k^2),
    // we choose b to approximate the missing Matern factor 
    // 1/(1+k^2)^(e-nfac). Specifically, we choose b to match a 
    // small Matern amplitude for some high wavenumber k. The
    // small amplitude was determined empirically by comparing the 
    // resulting covariance functions with Matern functions having 
    // the same shapes.
    //double akk = 0.001; // will match this Matern amplitude for some k^2
    //double kka = pow(akk,-1.0/e)-1.0; // k^2 at which to match amplitude
    //_bscl = (float)((pow(1.0+kka,e-_nfac)-1.0)/kka); // the factor b
    //System.out.println("e="+e+" nfac="+_nfac+" bscl="+_bscl);
    float[] ab = computeAB(e);
    _ascl = ab[0];
    _bscl = ab[1];

    // Choose spatial sampling of covariance so aliasing is negligible.
    // We use the same sampling for both x1 and x2 axes, because the
    // covariance function should have radial symmetry. Ours will not,
    // because of finite-difference approximations in the smoothing 
    // filters, but the approximations seem to be good enough.
    float rmax = 20.0f*_kscl;
    int nx = 4001;
    float dx = rmax/(nx-1);
    float fx = 0.0f;
    Sampling sx = new Sampling(nx,dx,fx);

    // Construct FFT and get samplings for wavenumber
    Fft fft = new Fft(sx,sx);
    Sampling sk1 = fft.getFrequencySampling1();
    Sampling sk2 = fft.getFrequencySampling2();
    int nk1 = sk1.getCount();
    int nk2 = sk2.getCount();

    // Compute the Fourier transform of our covariance function.
    // Account for finite-difference approximations in the smoothing
    // filters used to implement the factors 1/(1+k^2) and 1/(1+b*k^2)
    // described above. Note that the numerator is not one in this
    // Fourier transform (except where k1=k2=0) because the finite-
    // difference approximations include some small filters that
    // attenuate amplitudes for wavenumbers near Nyquist.
    float[][] g = new float[nk2][2*nk1];
    for (int ik2=0; ik2<nk2; ++ik2) {
      double k2 = 2.0*PI*sk2.getValue(ik2);
      for (int ik1=0; ik1<nk1; ++ik1) {
        double k1 = 2.0*PI*sk1.getValue(ik1);
        double cosk1 = cos(k1*dx);
        double cosk2 = cos(k2*dx);
        double kk = _kscl*_kscl*2.0*(1.0-cosk1*cosk2)/(dx*dx);
        double gnum = pow((0.5+0.5*cosk1)*(0.5+0.5*cosk2),2);
        double gden = pow(1.0+_ascl*kk,_nfac)*(1.0+_bscl*kk);
        g[ik2][2*ik1] = (float)(gnum/gden);
      }
    }

    // Compute the covariance function by inverse Fourier transform.
    // Use only one array in the returned array of arrays, because
    // the covariance function should have radial symmetry. Ours will
    // not, due to finite-difference approximations noted above, but 
    // here we ignore any anisotropy.
    float[][] f = fft.applyInverse(g);

    // Keep only the first half of the array obtained by inverse
    // Fourier transform. The second half contains the symmetric
    // left (negative r) part of the covariance function.
    nx = (nx-1)/2;
    float[] c = copy(nx,f[0]);

    // Scale covariance function c(r) so that c(0) = sigma^2.
    float cscl = (float)(_sigma*_sigma)/c[0];
    mul(cscl,c,c);

    // Scale factor used twice when applying this covariance function.
    _cscl = sqrt(cscl)*dx;

    // Construct the linear interpolator used to evaluate the
    // our uniformly sampled covariance function. Specify
    // constant extrapolation to be used when the range r is
    // beyond the sampled values.
    _li = new LinearInterpolator();
    _li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    _li.setUniform(nx,dx,fx,c);
  }
  
  private static float[] computeAB(double e) {
    int n = (int)e;
    double ak1 = 0.002; // Matern amplitude to match for k = k1
    double ak2 = 0.200; // Matern amplitude to match for k = k2
    double kk1 = pow(ak1,-1.0/e)-1.0; // k^2 for amplitude ak1
    double kk2 = pow(ak2,-1.0/e)-1.0; // k^2 for amplitude ak2
    // Solve for a and b such that
    // (1+a*kk1)^n (1+b*kk1) = 1/ak1
    // (1+a*kk2)^n (1+b*kk2) = 1/ak2
    double a = 1.0;
    double b = (1.0/(ak1*pow(1.0+a*kk1,n))-1.0)/kk1;
    double d = a+b;
    while (d>0.001) {
      double aold = a;
      double bold = b;
      a = (1.0/pow(ak2*(1.0+b*kk2),1.0/n)-1.0)/kk2;
      b = (1.0/(ak1*pow(1.0+a*kk1,n))-1.0)/kk1;
      if (b<0.0) b = 0.0;
      d = abs(a-aold)+abs(b-bold);
    }
    System.out.println("e="+e+" n="+n+" a="+a+" b="+b);
    return new float[]{(float)a,(float)b};
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
    double range = 80.0;
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
  private static float[][] VAB = {
    {0.00f, 1.00000f, 0.00000000f},
    {0.01f, 1.00831f, 0.00265717f},
    {0.02f, 1.01660f, 0.00535567f},
    {0.03f, 1.02484f, 0.00809606f},
    {0.04f, 1.03305f, 0.0108789f},
    {0.05f, 1.04123f, 0.0137049f},
    {0.06f, 1.04936f, 0.0165746f},
    {0.07f, 1.05746f, 0.0194887f},
    {0.08f, 1.06551f, 0.0224478f},
    {0.09f, 1.07353f, 0.0254526f},
    {0.10f, 1.08150f, 0.0285039f},
    {0.11f, 1.08942f, 0.0316023f},
    {0.12f, 1.09730f, 0.0347488f},
    {0.13f, 1.10513f, 0.0379439f},
    {0.14f, 1.11292f, 0.0411887f},
    {0.15f, 1.12065f, 0.0444838f},
    {0.16f, 1.12833f, 0.0478302f},
    {0.17f, 1.13596f, 0.0512288f},
    {0.18f, 1.14353f, 0.0546804f},
    {0.19f, 1.15105f, 0.0581861f},
    {0.20f, 1.15851f, 0.0617469f},
    {0.21f, 1.16591f, 0.0653637f},
    {0.22f, 1.17324f, 0.0690377f},
    {0.23f, 1.18051f, 0.0727698f},
    {0.24f, 1.18772f, 0.0765614f},
    {0.25f, 1.19486f, 0.0804134f},
    {0.26f, 1.20193f, 0.0843272f},
    {0.27f, 1.20893f, 0.0883040f},
    {0.28f, 1.21585f, 0.0923451f},
    {0.29f, 1.22270f, 0.0964518f},
    {0.30f, 1.22947f, 0.100626f},
    {0.31f, 1.23616f, 0.104868f},
    {0.32f, 1.24277f, 0.109180f},
    {0.33f, 1.24929f, 0.113564f},
    {0.34f, 1.25572f, 0.118020f},
    {0.35f, 1.26206f, 0.122552f},
    {0.36f, 1.26831f, 0.127160f},
    {0.37f, 1.27447f, 0.131846f},
    {0.38f, 1.28053f, 0.136613f},
    {0.39f, 1.28648f, 0.141462f},
    {0.40f, 1.29233f, 0.146395f},
    {0.41f, 1.29807f, 0.151414f},
    {0.42f, 1.30371f, 0.156522f},
    {0.43f, 1.30922f, 0.161720f},
    {0.44f, 1.31462f, 0.167012f},
    {0.45f, 1.31990f, 0.172399f},
    {0.46f, 1.32506f, 0.177885f},
    {0.47f, 1.33008f, 0.183471f},
    {0.48f, 1.33497f, 0.189161f},
    {0.49f, 1.33973f, 0.194958f},
    {0.50f, 1.34435f, 0.200865f},
    {0.51f, 1.34882f, 0.206885f},
    {0.52f, 1.35314f, 0.213022f},
    {0.53f, 1.35730f, 0.219278f},
    {0.54f, 1.36131f, 0.225659f},
    {0.55f, 1.36515f, 0.232167f},
    {0.56f, 1.36882f, 0.238808f},
    {0.57f, 1.37231f, 0.245585f},
    {0.58f, 1.37562f, 0.252503f},
    {0.59f, 1.37875f, 0.259567f},
    {0.60f, 1.38167f, 0.266782f},
    {0.61f, 1.38440f, 0.274154f},
    {0.62f, 1.38692f, 0.281688f},
    {0.63f, 1.38922f, 0.289392f},
    {0.64f, 1.39130f, 0.297271f},
    {0.65f, 1.39314f, 0.305332f},
    {0.66f, 1.39474f, 0.313583f},
    {0.67f, 1.39609f, 0.322032f},
    {0.68f, 1.39718f, 0.330688f},
    {0.69f, 1.39799f, 0.339560f},
    {0.70f, 1.39852f, 0.348658f},
    {0.71f, 1.39875f, 0.357992f},
    {0.72f, 1.39868f, 0.367575f},
    {0.73f, 1.39827f, 0.377418f},
    {0.74f, 1.39753f, 0.387537f},
    {0.75f, 1.39644f, 0.397945f},
    {0.76f, 1.39497f, 0.408660f},
    {0.77f, 1.39310f, 0.419699f},
    {0.78f, 1.39082f, 0.431082f},
    {0.79f, 1.38809f, 0.442833f},
    {0.80f, 1.38490f, 0.454974f},
    {0.81f, 1.38122f, 0.467534f},
    {0.82f, 1.37700f, 0.480545f},
    {0.83f, 1.37222f, 0.494041f},
    {0.84f, 1.36683f, 0.508063f},
    {0.85f, 1.36078f, 0.522658f},
    {0.86f, 1.35402f, 0.537880f},
    {0.87f, 1.34648f, 0.553792f},
    {0.88f, 1.33808f, 0.570470f},
    {0.89f, 1.32873f, 0.588005f},
    {0.90f, 1.31831f, 0.606509f},
    {0.91f, 1.30669f, 0.626121f},
    {0.92f, 1.29369f, 0.647020f},
    {0.93f, 1.27906f, 0.669439f},
    {0.94f, 1.26249f, 0.693697f},
    {0.95f, 1.24352f, 0.720252f},
    {0.96f, 1.22144f, 0.749798f},
    {0.97f, 1.19511f, 0.783496f},
    {0.98f, 1.16226f, 0.823593f},
    {0.99f, 1.11715f, 0.875836f},
    {1.00f, 1.00000f, 0.00000000f},
    {1.01f, 1.00416f, 0.00265408f},
    {1.02f, 1.00829f, 0.00534334f},
    {1.03f, 1.01241f, 0.00806835f},
    {1.04f, 1.01651f, 0.0108298f},
    {1.05f, 1.02059f, 0.0136282f},
    {1.06f, 1.02465f, 0.0164643f},
    {1.07f, 1.02869f, 0.0193387f},
    {1.08f, 1.03271f, 0.0222521f},
    {1.09f, 1.03671f, 0.0252053f},
    {1.10f, 1.04068f, 0.0281989f},
    {1.11f, 1.04463f, 0.0312337f},
    {1.12f, 1.04856f, 0.0343105f},
    {1.13f, 1.05246f, 0.0374300f},
    {1.14f, 1.05634f, 0.0405932f},
    {1.15f, 1.06020f, 0.0438007f},
    {1.16f, 1.06402f, 0.0470536f},
    {1.17f, 1.06783f, 0.0503526f},
    {1.18f, 1.07160f, 0.0536987f},
    {1.19f, 1.07534f, 0.0570929f},
    {1.20f, 1.07906f, 0.0605361f},
    {1.21f, 1.08275f, 0.0640294f},
    {1.22f, 1.08641f, 0.0675738f},
    {1.23f, 1.09003f, 0.0711703f},
    {1.24f, 1.09363f, 0.0748201f},
    {1.25f, 1.09719f, 0.0785243f},
    {1.26f, 1.10072f, 0.0822842f},
    {1.27f, 1.10422f, 0.0861009f},
    {1.28f, 1.10768f, 0.0899758f},
    {1.29f, 1.11110f, 0.0939101f},
    {1.30f, 1.11449f, 0.0979052f},
    {1.31f, 1.11784f, 0.101963f},
    {1.32f, 1.12115f, 0.106084f},
    {1.33f, 1.12442f, 0.110270f},
    {1.34f, 1.12765f, 0.114523f},
    {1.35f, 1.13084f, 0.118845f},
    {1.36f, 1.13398f, 0.123236f},
    {1.37f, 1.13708f, 0.127700f},
    {1.38f, 1.14014f, 0.132237f},
    {1.39f, 1.14315f, 0.136850f},
    {1.40f, 1.14610f, 0.141541f},
    {1.41f, 1.14901f, 0.146311f},
    {1.42f, 1.15187f, 0.151164f},
    {1.43f, 1.15468f, 0.156100f},
    {1.44f, 1.15743f, 0.161123f},
    {1.45f, 1.16013f, 0.166235f},
    {1.46f, 1.16277f, 0.171439f},
    {1.47f, 1.16535f, 0.176737f},
    {1.48f, 1.16787f, 0.182131f},
    {1.49f, 1.17033f, 0.187626f},
    {1.50f, 1.17273f, 0.193224f},
    {1.51f, 1.17505f, 0.198929f},
    {1.52f, 1.17731f, 0.204743f},
    {1.53f, 1.17950f, 0.21067f},
    {1.54f, 1.18162f, 0.216715f},
    {1.55f, 1.18366f, 0.222881f},
    {1.56f, 1.18562f, 0.229172f},
    {1.57f, 1.18750f, 0.235593f},
    {1.58f, 1.18930f, 0.242148f},
    {1.59f, 1.19101f, 0.248843f},
    {1.60f, 1.19264f, 0.255683f},
    {1.61f, 1.19417f, 0.262672f},
    {1.62f, 1.19560f, 0.269818f},
    {1.63f, 1.19694f, 0.277127f},
    {1.64f, 1.19817f, 0.284605f},
    {1.65f, 1.19929f, 0.292259f},
    {1.66f, 1.20030f, 0.300098f},
    {1.67f, 1.20119f, 0.308130f},
    {1.68f, 1.20197f, 0.316362f},
    {1.69f, 1.20261f, 0.324806f},
    {1.70f, 1.20312f, 0.333472f},
    {1.71f, 1.20349f, 0.342369f},
    {1.72f, 1.20372f, 0.351512f},
    {1.73f, 1.20379f, 0.360913f},
    {1.74f, 1.20370f, 0.370586f},
    {1.75f, 1.20344f, 0.380548f},
    {1.76f, 1.20300f, 0.390815f},
    {1.77f, 1.20237f, 0.401407f},
    {1.78f, 1.20154f, 0.412345f},
    {1.79f, 1.20049f, 0.423653f},
    {1.80f, 1.19922f, 0.435357f},
    {1.81f, 1.19770f, 0.447487f},
    {1.82f, 1.19592f, 0.460076f},
    {1.83f, 1.19386f, 0.473163f},
    {1.84f, 1.19150f, 0.486791f},
    {1.85f, 1.18880f, 0.501010f},
    {1.86f, 1.18574f, 0.515881f},
    {1.87f, 1.18228f, 0.531473f},
    {1.88f, 1.17839f, 0.547868f},
    {1.89f, 1.17400f, 0.565167f},
    {1.90f, 1.16905f, 0.583494f},
    {1.91f, 1.16347f, 0.603002f},
    {1.92f, 1.15715f, 0.623892f},
    {1.93f, 1.14996f, 0.646425f},
    {1.94f, 1.14172f, 0.670959f},
    {1.95f, 1.13216f, 0.698009f},
    {1.96f, 1.12090f, 0.728362f},
    {1.97f, 1.10726f, 0.763338f},
    {1.98f, 1.08996f, 0.805497f},
    {1.99f, 1.06570f, 0.861426f},
    {2.00f, 1.00000f, 0.00000000f},
    {2.01f, 1.00277f, 0.00265305f},
    {2.02f, 1.00553f, 0.00533919f},
    {2.03f, 1.00827f, 0.00805902f},
    {2.04f, 1.01101f, 0.0108131f},
    {2.05f, 1.01373f, 0.0136022f},
    {2.06f, 1.01643f, 0.0164267f},
    {2.07f, 1.01912f, 0.0192875f},
    {2.08f, 1.02180f, 0.0221851f},
    {2.09f, 1.02446f, 0.0251203f},
    {2.10f, 1.02711f, 0.0280938f},
    {2.11f, 1.02974f, 0.0311063f},
    {2.12f, 1.03236f, 0.0341585f},
    {2.13f, 1.03496f, 0.0372513f},
    {2.14f, 1.03754f, 0.0403854f},
    {2.15f, 1.04011f, 0.0435617f},
    {2.16f, 1.04266f, 0.0467810f},
    {2.17f, 1.04519f, 0.0500442f},
    {2.18f, 1.04770f, 0.0533522f},
    {2.19f, 1.05020f, 0.0567059f},
    {2.20f, 1.05267f, 0.0601063f},
    {2.21f, 1.05513f, 0.0635543f},
    {2.22f, 1.05757f, 0.0670511f},
    {2.23f, 1.05999f, 0.0705976f},
    {2.24f, 1.06238f, 0.0741949f},
    {2.25f, 1.06476f, 0.0778442f},
    {2.26f, 1.06711f, 0.0815467f},
    {2.27f, 1.06944f, 0.0853035f},
    {2.28f, 1.07175f, 0.0891160f},
    {2.29f, 1.07404f, 0.0929853f},
    {2.30f, 1.07630f, 0.0969129f},
    {2.31f, 1.07854f, 0.100900f},
    {2.32f, 1.08075f, 0.104948f},
    {2.33f, 1.08293f, 0.109059f},
    {2.34f, 1.08509f, 0.113234f},
    {2.35f, 1.08723f, 0.117475f},
    {2.36f, 1.08933f, 0.121783f},
    {2.37f, 1.09141f, 0.126160f},
    {2.38f, 1.09345f, 0.130609f},
    {2.39f, 1.09547f, 0.135130f},
    {2.40f, 1.09745f, 0.139726f},
    {2.41f, 1.09941f, 0.144399f},
    {2.42f, 1.10133f, 0.149151f},
    {2.43f, 1.10321f, 0.153984f},
    {2.44f, 1.10507f, 0.158900f},
    {2.45f, 1.10688f, 0.163903f},
    {2.46f, 1.10866f, 0.168994f},
    {2.47f, 1.11040f, 0.174177f},
    {2.48f, 1.11211f, 0.179453f},
    {2.49f, 1.11377f, 0.184827f},
    {2.50f, 1.11539f, 0.190300f},
    {2.51f, 1.11697f, 0.195877f},
    {2.52f, 1.11851f, 0.201561f},
    {2.53f, 1.12000f, 0.207354f},
    {2.54f, 1.12144f, 0.213262f},
    {2.55f, 1.12284f, 0.219288f},
    {2.56f, 1.12418f, 0.225435f},
    {2.57f, 1.12547f, 0.231710f},
    {2.58f, 1.12671f, 0.238115f},
    {2.59f, 1.12790f, 0.244657f},
    {2.60f, 1.12902f, 0.251340f},
    {2.61f, 1.13009f, 0.258171f},
    {2.62f, 1.13110f, 0.265155f},
    {2.63f, 1.13204f, 0.272298f},
    {2.64f, 1.13291f, 0.279607f},
    {2.65f, 1.13371f, 0.287089f},
    {2.66f, 1.13445f, 0.294753f},
    {2.67f, 1.13510f, 0.302607f},
    {2.68f, 1.13568f, 0.310659f},
    {2.69f, 1.13617f, 0.318920f},
    {2.70f, 1.13658f, 0.327399f},
    {2.71f, 1.13690f, 0.336109f},
    {2.72f, 1.13712f, 0.345061f},
    {2.73f, 1.13724f, 0.354269f},
    {2.74f, 1.13726f, 0.363747f},
    {2.75f, 1.13716f, 0.373513f},
    {2.76f, 1.13695f, 0.383582f},
    {2.77f, 1.13661f, 0.393976f},
    {2.78f, 1.13615f, 0.404715f},
    {2.79f, 1.13554f, 0.415824f},
    {2.80f, 1.13477f, 0.427330f},
    {2.81f, 1.13385f, 0.439263f},
    {2.82f, 1.13276f, 0.451658f},
    {2.83f, 1.13148f, 0.464553f},
    {2.84f, 1.12999f, 0.477996f},
    {2.85f, 1.12828f, 0.492036f},
    {2.86f, 1.12633f, 0.506736f},
    {2.87f, 1.12412f, 0.522167f},
    {2.88f, 1.12160f, 0.538416f},
    {2.89f, 1.11876f, 0.555586f},
    {2.90f, 1.11554f, 0.573807f},
    {2.91f, 1.11188f, 0.593238f},
    {2.92f, 1.10772f, 0.614089f},
    {2.93f, 1.10297f, 0.636631f},
    {2.94f, 1.09749f, 0.661242f},
    {2.95f, 1.09111f, 0.688459f},
    {2.96f, 1.08355f, 0.719111f},
    {2.97f, 1.07433f, 0.754585f},
    {2.98f, 1.06256f, 0.797583f},
    {2.99f, 1.04590f, 0.855063f},
    {3.00f, 1.00000f, 0.000000f}
  };
}
