package fault;

import edu.mines.jtk.dsp.SincInterp;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Synthetic warping of 1-D sequences.
 * The function u(x) describes the warping, in which a point x
 * is displaced to a point y(x) = x+u(x).
 * <p>
 * Warping is the computation of the sequence g(y) = f(x(y)).
 * Unwarping is the computation of the sequence f(x) = g(y(x)).
 * <p>
 * For warping, we need the function x(y) = y-u(x(y)) = y-uy(y). We 
 * compute the displacement uy(y) by iteration so that uy(y) = u(x(y)).
 * <p>
 * We also define a midpoint m(x) = (x+y(x))/2, and compute the 
 * displacement um(m) = u(x(m)) from u(x) by iteration so that 
 * um(m) = u(x(m)).
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.08.30
 */
public abstract class Warp1 {

  public static Warp1 constant(double u, int n) {
    return new ConstantWarp1(u,n);
  }

  public static Warp1 gaussian(double u, int n) {
    return new GaussianWarp1(u,n);
  }

  public static Warp1 sinusoid(double u, int n) {
    return new SinusoidWarp1(u,n);
  }

  public abstract double u(double x);

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
    SincInterp si = new SincInterp();
    float[] g = new float[_n];
    for (int i=0; i<_n; ++i) {
      double y = i;
      double x = y-uy(y);
      g[i] = si.interpolate(_n,1.0,0.0,f,x);
    }
    return g;
  }

  public float[] unwarp(float[] g) {
    SincInterp si = new SincInterp();
    float[] f = new float[_n];
    for (int i=0; i<_n; ++i) {
      double x = i;
      double y = x+ux(x);
      f[i] = si.interpolate(_n,1.0,0.0,g,y);
    }
    return f;
  }

  protected Warp1(int n) {
    _n = n;
  }

  private int _n;

  /**
   * Constant (zero-strain) displacement.
   */
  private static class ConstantWarp1 extends Warp1 {
    public ConstantWarp1(double u, int n) {
      super(n);
      _u = u;
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
  private static class GaussianWarp1 extends Warp1 {
    public GaussianWarp1(double umax, int n) {
      super(n);
      _a = (n-1)/2.0;
      _b = _a/3;
      _c = umax*exp(0.5)/_b;
      _umax = umax;
      _emax = _c;
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
  private static class SinusoidWarp1 extends Warp1 {
    public SinusoidWarp1(double umax, int n) {
      super(n);
      double l = n-1;
      _a = umax;
      _b = 2.0*PI/l;
      _umax = umax;
      _emax = _a*_b;
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
