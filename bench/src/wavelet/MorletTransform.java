package wavelet;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A continuous wavelet transform with the Morlet wavelet.
 * Transforms real-valued functions of time x(t) to complex-valued functions
 * of time and frequency y(t,f). The Morlet transform is equivalent to a
 * collection of Fourier transforms for seamlessly overlapping Gaussian
 * windows, where window width is inversely proportional to frequency.
 * <p>
 * In the Morlet transform, the logarithm of frequency is sampled uniformly,
 * and frequency is sampled non-uniformly. Specifically, the increment between
 * sampled frequencies increases linearly with frequency. This class provides
 * methods that simplify specification of this log-frequency sampling.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.12.31
 */
public class MorletTransform {

  /**
   * Returns a sampling of frequencies that is uniform in log-frequency.
   * Note that the returned sampling is of frequencies, not log-frequencies;
   * the returned sampling will not be uniform in frequency.
   * @param nf number of frequencies; must be greater than one.
   * @param fmin minimum frequency; must be greater than zero.
   * @param fmax maximum frequency; must be greater than fmin.
   */
  public static Sampling frequencySampling(int nf, double fmin, double fmax) {
    Check.argument(nf>1,"nf>1");
    Check.argument(fmin>0.0,"fmin>0.0");
    Check.argument(fmin<fmax,"fmin<fmax");
    double logfmin = log10(fmin);
    double logfmax = log10(fmax);
    double dlogf = (logfmax-logfmin)/(nf-1);
    double flogf = logfmin;
    double[] fs = new double[nf];
    for (int jf=0; jf<nf; ++jf)
      fs[jf] = pow(10.0f,flogf+jf*dlogf);
    return new Sampling(fs);
  }

  /**
   * Construct a Morlet transform for specified time and frequency sampling.
   * Frequencies will be sampled non-uniformly so that log-frequency is 
   * sampled uniformly.
   * @param st sampling of time.
   * @param nf number of frequencies; must be greater than 1.
   * @param fmin minimum frequency; must be positive.
   * @param fmax maximum frequency; must not be less than fmin.
   */
  public MorletTransform(Sampling st, int nf, double fmin, double fmax) {
    this(st,frequencySampling(nf,fmin,fmax));
  }

  /**
   * Construct a Morlet transform for specified time and frequency sampling.
   * With this constructor, an arbitrary sampling of frequencies may be
   * specified, subject to constraints that the minumum frequency must be
   * greater than zero, and the maximum frequency must not exceed the Nyquist
   * frequency.
   * @param st sampling of time; must be uniform.
   * @param sf sampling of frequency.
   */
  public MorletTransform(Sampling st, Sampling sf) {
    Check.argument(st.isUniform(),"st is uniform");
    Check.argument(sf.getFirst()>0.0,"minimum f > 0");
    Check.argument(sf.getLast()<=0.5/st.getDelta(),"maximum f <= Nyquist");
    // Here we set w0, the frequency (in radians/sample) of the base Morlet
    // wavelet with scale sigma = 1. Typically w0 > 5 so that the sum of the
    // Morlet wavelet coefficients is approximately zero. The scale sigma of
    // the Morlet wavelet is inversely proportional to its frequency w; that
    // is, sigma = w0/w = w0/(2*PI*f), where f denotes frequency in
    // cycles/sample. The choice w0 = 5.336 is typical. For this choice, the
    // value of the second-highest peak in the real (cosine) part of the
    // Morlet wavelet m(t) equals 1/2 the value of the highest peak m(t=0).
    double w0 = 5.336;
    double sigmaScale = w0/(2.0*PI); // sigma = w0/(2*PI*f)
    int nt = st.getCount();
    int nf = sf.getCount();
    double dt = st.getDelta();
    _ks = new Kernel[nf];
    for (int jf=0; jf<nf; ++jf) {
      double freqj = sf.getValue(jf);
      double freqn = freqj*dt; // normalized freq (cycles/sample)
      double sigma = (freqn>0.0)?sigmaScale/freqn:3.0*nt;
      _ks[jf] = new Kernel(nt,sigma,freqn);
    }
    _st = st;
    _sf = sf;
  }

  /**
   * Gets the sampling of time for this transform.
   * @return sampling of frequency.
   */
  public Sampling getTimeSampling() {
    return _st;
  }

  /**
   * Gets the sampling of frequency for this transform.
   * The returned sampling may or may not be uniform, depending on how this
   * transform was constructed.
   * @return sampling of frequency.
   */
  public Sampling getFrequencySampling() {
    return _sf;
  }

  /**
   * Gets the sampling of log-frequency for this transform.
   * Here, the logarithm is a power of 10 (not a power of e). The returned
   * sampling may or may not be uniform, depending on how this transform was
   * constructed.
   * @return sampling of log-frequency.
   */
  public Sampling getLogFrequencySampling() {
    int nf = _sf.getCount();
    double[] logf = new double[nf];
    for (int jf=0; jf<nf; ++jf)
      logf[jf] = log10(_sf.getValue(jf));
    Sampling slogf = new Sampling(logf);
    return slogf;
  }

  /**
   * Returns a transformed array {yr,yi} of complex-valued y(t,f).
   * The returned array yr contains the real parts of y(t,f);
   * the returned array yi contains the imaginary parts of y(t,f).
   * @param x input array[nt] of x(t) to be transformed.
   * @return array[2][nf][nt] {yr,yi} of y(t,f).
   */
  public float[][][] apply(float[] x) {
    int nf = _sf.getCount();
    int nt = _st.getCount();
    float[][][] y = new float[2][nf][nt];
    float[][] yr = y[0];
    float[][] yi = y[1];
    for (int jf=0; jf<nf; ++jf)
      _ks[jf].apply(x,yr[jf],yi[jf]);
    return y;
  }

  /**
   * Returns the magnitude (abs) of a complex-valued y(t,f).
   * @param array[2][nf][nt] {yr,yi} of y(t,f).
   * @return array[nf][nt] of magnitudes of y(t,f).
   */
  public static float[][] abs(float[][][] y) {
    int nt = y[0][0].length;
    int nf = y[0].length;
    float[][] yr = y[0];
    float[][] yi = y[1];
    float[][] ya = new float[nf][nt];
    for (int jf=0; jf<nf; ++jf) {
      for (int jt=0; jt<nt; ++jt) {
        float yrj = yr[jf][jt];
        float yij = yi[jf][jt];
        ya[jf][jt] = sqrt(yrj*yrj+yij*yij);
      }
    }
    return ya;
  }

  /**
   * Returns the phase (arg, in cycles) of a complex-valued y(t,f).
   * @param array[2][nf][nt] {yr,yi} of y(t,f).
   * @return array[nf][nt] of phases of y(t,f).
   */
  public static float[][] arg(float[][][] y) {
    int nt = y[0][0].length;
    int nf = y[0].length;
    float[][] yr = y[0];
    float[][] yi = y[1];
    float[][] ya = new float[nf][nt];
    float otwopi = (float)(0.5/PI);
    for (int jf=0; jf<nf; ++jf) {
      for (int jt=0; jt<nt; ++jt) {
        float yrj = yr[jf][jt];
        float yij = yi[jf][jt];
        ya[jf][jt] = atan2(yij,yrj)*otwopi;
      }
    }
    return ya;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  
  private Sampling _st; // time sampling
  private Sampling _sf; // frequency sampling
  private Kernel[] _ks; // array of kernels, one per frequency

  private static class Kernel {
    Kernel(int n, double sigma, double freqn) {
      _n = n;
      _rgf = new RecursiveGaussianFilter(sigma);
      // The recursive Gaussian filter scales by 1/(sqrt(2*PI)*sigma),
      // but Morlet wavelet should scale by 1/sqrt(sqrt(PI)*sigma), so
      // we need this scale factor to convert the former to the latter.
      _scale = (float)sqrt(sqrt(PI)*2.0*sigma);
      float w = (float)(2.0*PI*freqn);
      _c = new float[n];
      _s = new float[n];
      for (int i=0; i<n; ++i) {
        _c[i] = cos(w*i);
        _s[i] = sin(w*i);
      }
    }
    void apply(float[] x, float[] yr, float[] yi) {
      Check.argument(x.length==_n,"x.length == n");
      Check.argument(yr.length==_n,"yr.length == n");
      Check.argument(yi.length==_n,"yi.length == n");
      for (int i=0; i<_n; ++i) {
        float sx = _scale*x[i];
        yr[i] =  _c[i]*sx;
        yi[i] = -_s[i]*sx;
      }
      _rgf.apply0(yr,yr);
      _rgf.apply0(yi,yi);
      for (int i=0; i<_n; ++i) {
        float zr = yr[i];
        float zi = yi[i];
        yr[i] = _c[i]*zr-_s[i]*zi;
        yi[i] = _s[i]*zr+_c[i]*zi;
      }
    }
    private RecursiveGaussianFilter _rgf;
    private int _n; // number of samples input/output
    private float _scale; // scale factor (converts RGF to Morlet)
    private float[] _c; // cos(2*pi*f*i), for i = 0,1,...,n-1
    private float[] _s; // sin(2*pi*f*i), for i = 0,1,...,n-1
  }
}
