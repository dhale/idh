package lcc;

import java.awt.*;
import java.util.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

import util.*;

public class Warp1 {

  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Warp1(args);
      }
    });
  }

  private int _fontSize = 24;
  private int _width = 900;
  private int _height = 800;

  private static final LocalCorrelationFilter.Type SIMPLE =
    LocalCorrelationFilter.Type.SIMPLE;
  private static final LocalCorrelationFilter.Type SYMMETRIC =
    LocalCorrelationFilter.Type.SYMMETRIC;

  private static final LocalCorrelationFilter.Window GAUSSIAN = 
    LocalCorrelationFilter.Window.GAUSSIAN;
  private static final LocalCorrelationFilter.Window RECTANGLE = 
    LocalCorrelationFilter.Window.RECTANGLE;

  private int _length = 315;
  private float _flo = 0.15f;
  private float _fhi = 0.25f;
  private float _dmax = 0.5f;
  private int _lmax = 1;
  private int _lmin = -_lmax;
  private LocalCorrelationFilter.Type _type = SIMPLE; 
  private LocalCorrelationFilter.Window _window = GAUSSIAN; 
  private float _sigma = 8.0f;
  private Displacement _disp = new SinusoidDisplacement(_dmax,_length);
  //private Displacement _disp = new ConstantDisplacement(_dmax,_length);
  private float[] _sequence = makeRandom(_length,_flo,_fhi);

  private Warp1(String[] args) {
    //testSimpleSymmetric();
    testGaussianRectangle();
  }

  private void testSimpleSymmetric() {
    _window = RECTANGLE;
    _type = SIMPLE;
    testWarp();
    _type = SYMMETRIC;
    testWarp();
  }

  private void testGaussianRectangle() {
    _window = GAUSSIAN;
    testWarp();
    _window = RECTANGLE;
    testWarp();
  }

  private void testWarp() {
    float sigma = _sigma;
    if (_window==RECTANGLE)
      sigma = (float)(0.5*(sqrt(2.0*PI)*sigma-1.0)); // area same as Gaussian
    int n = _length;
    int lmin = _lmin;
    int lmax = _lmax;
    int nlag = 1+lmax-lmin;
    float[] f = _sequence;
    float[] g = _disp.warp(f);
    float[] d = (_type==SIMPLE)?_disp.ux():_disp.um();
    float[][] c = new float[nlag][n];
    float[] u = new float[n];
    lcc(sigma,lmin,lmax,f,g,c,u);
    plot(f,g,c,u,d);
  }

  private static float[] applyGaussianWindow(double x, double s, float[] f) {
    int n = f.length;
    float[] g = new float[n];
    for (int i=0; i<n; ++i) {
      double e = (i-x)/s;
      g[i] = f[i]*(float)exp(-0.5*e*e);
    }
    return g;
  }

  private float[] makeRandom(int n, double fl, double fh) {

    // Random sequence in [-1,1].
    float[] f = Array.randfloat(new Random(314159),n);
    f = Array.sub(Array.mul(2.0f,f),1.0f);

    // Lowpass filter.
    ButterworthFilter lp = new ButterworthFilter(0.01,0.01,0.02,0.99);
    lp.applyForward(f,f);

    // Highpass filter.
    ButterworthFilter hp = new ButterworthFilter(fl,0.99,fh,0.01);
    hp.applyForward(f,f);

    return f;
  }

  private float[] makeSweep(int n, double fl, double fh) {
    float[] f = new float[n];
    for (int i=0; i<n; ++i) {
      double phase = 2.0*DBL_PI*(fl+0.5*i*(fh-fl)/(n-1))*i;
      f[i] = (float)cos(phase);
    }
    return f;
  }

  private void lcc(
    double sigma, int lmin, int lmax, float[] f, float[] g,
    float[][] c, float[] u)
  {
    int n = f.length;
    int nlag = 1+lmax-lmin;
    LocalCorrelationFilter.Type type = _type;
    LocalCorrelationFilter.Window window = _window;
    LocalCorrelationFilter lcf = 
      new LocalCorrelationFilter(type,window,sigma);
    lcf.setInputs(f,g);
    for (int lag=lmin; lag<=lmax; ++lag) {
      int k = lag-lmin;
      lcf.correlate(lag,c[k]);
      lcf.normalize(lag,c[k]);
    }
    byte[] l = new byte[n];
    lcf.findMaxLags(lmin,lmax,l);
    lcf.refineLags(l,u);
  }

  private void plot(float[] f, float[] g, float[][] c, float[] u, float[] d) {
    int n = f.length;
    Sampling s = new Sampling(n,1.0,0.0);
    int nlag = c.length;
    Sampling slag = new Sampling(nlag,1.0,-(nlag-1)/2);
    PlotPanel panel = new PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP);
    SequenceView fv = panel.addSequence(0,0,s,f);
    SequenceView gv = panel.addSequence(1,0,s,g);
    PixelsView cv = panel.addPixels(2,0,s,slag,c);
    cv.setClips(-1.0f,1.0f);
    cv.setColorModel(ColorMap.JET);
    PointsView dv = panel.addPoints(2,0,s,d);
    dv.setLineColor(Color.WHITE);
    dv.setLineStyle(PointsView.Line.DOT);
    PointsView uv = panel.addPoints(2,0,s,u);
    uv.setLineColor(Color.WHITE);
    panel.setHLabel("sample");
    panel.setVLabel(0,"f");
    panel.setVLabel(1,"g");
    panel.setVLabel(2,"lag");
    frame(panel);
  }

  private PlotFrame frame(PlotPanel panel) {
    PlotFrame frame = new PlotFrame(panel);
    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    frame.setFontSize(_fontSize);
    frame.setSize(_width,_height);
    frame.setVisible(true);
    return frame;
  }

  /**
   * Abstract base class for synthetic displacements.
   * The function u(x) is displacement. The function e(x) is strain.
   * A point x is displaced to a point y(x) = x+u(x).
   * <p>
   * Warping is the computation of the sequence g(y) = f(x(y)).
   * Unwarping is the computation of the sequence f(x) = g(y(x)).
   * <p>
   * For warping, we need the function x(y) = y-u(x(y)) = y-uy(y). We 
   * compute the displacement uy(y) by iteration so that uy(y) = u(x(y).
   * <p>
   * We also define a midpoint m(x) = (x+y(x))/2, and compute the 
   * displacement um(m) = u(x(m)) from u(x) by iteration so that 
   * um(m) = u(x(m)).
   */
  private static abstract class Displacement {
    public abstract double u(double x);
    public abstract double e(double x);
    public abstract double umax();
    public abstract double emax();
    public double ux(double x) {
      return u(x);
    }
    public double um(double m) {
      double um = 0.0;
      double up;
      do {
        up = um;
        um = u(m-0.5*um);
      } while (abs(um-up)>0.0001);
      return um;
    }
    public double uy(double y) {
      double uy = 0.0;
      double up;
      do {
        up = uy;
        uy = u(y-uy);
      } while (abs(uy-up)>0.0001);
      return uy;
    }
    public float[] ux() {
      float[] u = new float[_n];
      for (int i=0; i<_n; ++i) {
        double x = i;
        u[i] = (float)ux(x);
      }
      return u;
    }
    public float[] um() {
      float[] u = new float[_n];
      for (int i=0; i<_n; ++i) {
        double m = i;
        u[i] = (float)um(m);
      }
      return u;
    }
    public float[] uy() {
      float[] u = new float[_n];
      for (int i=0; i<_n; ++i) {
        double y = i;
        u[i] = (float)uy(y);
      }
      return u;
    }
    public float[] warp(float[] f) {
      SincInterpolator si = new SincInterpolator();
      si.setUniform(_n,1.0,0.0,f);
      float[] g = new float[_n];
      for (int i=0; i<_n; ++i) {
        double y = i;
        double x = y-uy(y);
        g[i] = si.interpolate(x);
      }
      return g;
    }
    public float[] unwarp(float[] g) {
      SincInterpolator si = new SincInterpolator();
      si.setUniform(_n,1.0,0.0,g);
      float[] f = new float[_n];
      for (int i=0; i<_n; ++i) {
        double x = i;
        double y = x+ux(x);
        f[i] = si.interpolate(y);
      }
      return f;
    }
    protected Displacement(int n) {
      _n = n;
    }
    private int _n;
  }

  /**
   * Constant (zero-strain) displacement.
   */
  private static class ConstantDisplacement extends Displacement {
    public ConstantDisplacement(double u, int n) {
      super(n);
      _u = u;
      System.out.println("ConstantDisplacement: max u="+umax());
      System.out.println("ConstantDisplacement: max e="+emax());
    }
    public double u(double x) {
      return _u;
    }
    public double umax() {
      return _u;
    }
    public double e(double x) {
      return 0.0;
    }
    public double emax() {
      return 0.0;
    }
    private double _u;
  }

  /**
   * Derivative-of-Gaussian displacement.
   */
  private static class GaussianDisplacement extends Displacement {
    public GaussianDisplacement(double umax, int n) {
      super(n);
      _a = (n-1)/2.0;
      _b = _a/3;
      _c = umax*exp(0.5)/_b;
      _umax = umax;
      _emax = _c;
      System.out.println("GaussianDisplacement: max u="+umax());
      System.out.println("GaussianDisplacement: max e="+emax());
    }
    public double u(double x) {
      double xa = x-_a;
      return -_c*xa*exp(-0.5*(xa*xa)/(_b*_b));
    }
    public double umax() {
      return _umax;
    }
    public double e(double x) {
      double xa = x-_a;
      return -_c*(1.0-(xa*xa)/(_b*_b))*exp(-0.5*(xa*xa)/(_b*_b));
    }
    public double emax() {
      return _emax;
    }
    private double _a;
    private double _b;
    private double _c;
    private double _umax;
    private double _emax;
  }

  /**
   * Sinusoid displacement.
   */
  private static class SinusoidDisplacement extends Displacement {
    public SinusoidDisplacement(double umax, int n) {
      super(n);
      double l = n-1;
      _a = umax;
      _b = 2.0*PI/l;
      _umax = umax;
      _emax = _a*_b;
      System.out.println("SinusoidDisplacement: max u="+umax());
      System.out.println("SinusoidDisplacement: max e="+emax());
    }
    public double u(double x) {
      return _a*sin(_b*x);
    }
    public double umax() {
      return _umax;
    }
    public double e(double x) {
      return _a*_b*cos(_b*x);
    }
    public double emax() {
      return _emax;
    }
    private double _a;
    private double _b;
    private double _umax;
    private double _emax;
  }
}
