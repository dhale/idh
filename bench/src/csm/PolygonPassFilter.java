package csm;

import edu.mines.jtk.dsp.FftReal;
import edu.mines.jtk.dsp.FftComplex;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A 2D filter that passes a specified polygonal area in the frequency domain.
 * This filter is a generalization of a 1D band-pass filter. Frequencies to be
 * passed are those that lie inside a closed polygon with specified vertices.
 * Frequencies to be rejected lie outside the polgyon.
 * <p>
 * For efficiency, this filter is applied using fast Fourier transforms. It
 * has a frequency response that only approximates the ideal frequency
 * response, equal to one inside the polygon and zero outside. Approximation
 * enables the spatial extent of the impulse response to be limited by
 * specified half-widths. In effect, these half-widths control the sharpness
 * of the transition between passed and rejected frequencies. Smaller
 * half-widths imply a smoother transition, while larger half-widths imply a
 * sharper transition.
 * <p>
 * When this filter is applied, samples outside the specified input array are
 * assumed to be zero. This assumption is possible only because the impulse
 * response of this filter has finite spatial extent.
 * <p>
 * A polygon pass filter is thread safe. Once constructed, it may be applied
 * in parallel to multiple 2D arrays.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.10.13
 */
public class PolygonPassFilter {

  /**
   * Constructs a filter for a specified polygon in the frequency domain.
   * The polygon must be closed; the last vertex must equal the first vertex.
   * Polygon vertices (f1,f2) are frequencies in cycles/sample, which means
   * that the Nyquist frequency is 0.5. Typically, frequencies f1 are in the
   * range [0,0.5], and frequencies f2 are in the range [-0.5,0.5]. Polygon
   * vertices (f1,f2) outside of these ranges will yield a frequency response
   * that is aliased.
   * <p> 
   * The frequency response of the filter is symmetric about the origin.
   * Therefore, for each polygon specified, a corresponding polygon with 
   * negated vertices is assumed, and should not be specified.
   * @param n1 max number of samples in 1st dimension of input/output arrays.
   * @param n2 max number of samples in 2nd dimension of input/output arrays.
   * @param m1 half-width of filter in 1st dimension.
   * @param m2 half-width of filter in 2nd dimension.
   * @param f1 array of frequencies (cycles/sample) for 1st dimension.
   * @param f2 array of frequencies (cycles/sample) for 2nd dimension.
   */
  public PolygonPassFilter(
    int n1, int n2, int m1, int m2, float[] f1, float[] f2) 
  {
    Check.argument(f1.length==f2.length,"length of f1 == length of f2");
    Check.argument(f1[0]==f1[f1.length-1],"first f1 == last f1");
    Check.argument(f2[0]==f2[f2.length-1],"first f2 == last f2");
    _n1 = n1;
    _n2 = n2;
    _m1 = m1;
    _m2 = m2;
    _f1 = f1;
    _f2 = f2;

    // FFTs for 1st and 2nd dimensions.
    _n1fft = FftReal.nfftSmall(n1+m1);
    _n2fft = FftComplex.nfftSmall(n2+m2);
    _fft1 = new FftReal(_n1fft);
    _fft2 = new FftComplex(_n2fft);

    // Numbers of frequencies in 1st and 2nd dimensions.
    _n1f = 1+_n1fft/2;
    _n2f = _n2fft;

    // Compute windowed h(x1,x2) for integer x1 and x2 by accumulating
    // contributions for each edge (and its symmetric mate).
    double twopi = 2.0*FLT_PI;
    double[][] h = new double[_m2+1+_m2][_m1+1+_m1];
    for (int jf=1; jf<_f1.length; ++jf) {
      double a1 = twopi*_f1[jf-1];
      double a2 = twopi*_f2[jf-1];
      double b1 = twopi*_f1[jf];
      double b2 = twopi*_f2[jf];
      double c1 = 0.5*(b1+a1);
      double c2 = 0.5*(b2+a2);
      double d1 = 0.5*(b1-a1);
      double d2 = 0.5*(b2-a2);
      double hs = 0.5/(FLT_PI*FLT_PI);
      for (int k2=-_m2,j2=0; k2<=_m2; ++k2,++j2) {
        double x2 = k2;
        double w2 = hamming(_m2,x2);
        for (int k1=-_m1,j1=0; k1<=_m1; ++k1,++j1) {
          double x1 = k1;
          double w1 = hamming(_m1,x1);
          double ht;
          if (x1==0.0 && x2==0.0) {
            ht = c1*d2-c2*d1;
          } else if (x1==0.0) {
            ht = -2.0*c2*d1*sinc(c2*x2)*sinc(d2*x2);
          } else if (x2==0.0) {
            ht =  2.0*c1*d2*sinc(c1*x1)*sinc(d1*x1);
          } else {
            ht = (d2/x1-d1/x2)*sin(c1*x1+c2*x2)*sinc(d1*x1+d2*x2);
          }
          h[j2][j1] += w1*w2*hs*ht;
        }
      }
    }

    // H(k1,k2) is the 2D Fourier transform of h(x1,x2). Include the FFT scale
    // factor here so that we need no scaling when we later apply the filter.
    float[][] hfft = new float[_n2f][_n1f*2];
    for (int k2=-_m2,j2=_n2fft-_m2; k2<=_m2; ++k2,++j2) {
      if (j2==_n2fft) j2 = 0;
      for (int k1=-_m1,j1=_n1fft-_m1; k1<=_m1; ++k1,++j1) {
        if (j1==_n1fft) j1 = 0;
        hfft[j2][j1] = (float)(h[_m2+k2][_m1+k1]/_n1fft/_n2fft);
      }
    }
    _fft1.realToComplex1(-1,_n2f,hfft,hfft);
    _fft2.complexToComplex2(-1,_n1f,hfft,hfft);

    // Keep only the real part. The imaginary parts should be zero.
    _h = new float[_n2f][_n1f];
    for (int i2=0; i2<_n2f; ++i2) {
      for (int i1=0; i1<_n1f; ++i1) {
        _h[i2][i1] = hfft[i2][i1*2];
      }
    }
  }

  /**
   * Applies this filter.
   * May be applied in-place; input and output arrays may be the same array.
   * @param x input array.
   * @param y output array.
   */
  public void apply(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    Check.argument(n1<=_n1,"dimension n1 of x is small enough");
    Check.argument(n2<=_n2,"dimension n2 of x is small enough");

    // Pad input array with zeros and forward Fourier transform.
    float[][] xy = new float[_n2fft][_n1fft+2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xy[i2][i1] = x[i2][i1];
      }
      for (int i1=n1; i1<_n1fft; ++i1) {
        xy[i2][i1] = 0.0f;
      }
    }
    for (int i2=n2; i2<_n2fft; ++i2) {
      for (int i1=0; i1<_n1fft; ++i1) {
        xy[i2][i1] = 0.0f;
      }
    }
    _fft1.realToComplex1(-1,_n2,xy,xy);
    _fft2.complexToComplex2(-1,_n1f,xy,xy);

    // Multiply by real-valued frequency response of this filter.
    for (int i2=0; i2<_n2f; ++i2) {
      for (int i1=0,i1r=0,i1i=1; i1<_n1f; ++i1,i1r+=2,i1i+=2) {
        xy[i2][i1r] *= _h[i2][i1];
        xy[i2][i1i] *= _h[i2][i1];
      }
    }

    // Inverse Fourier transform and copy to output array.
    _fft2.complexToComplex2(1,_n1f,xy,xy);
    _fft1.complexToReal1(1,_n2,xy,xy);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        y[i2][i1] = xy[i2][i1];
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2,_m1,_m2;
  private float[] _f1,_f2;
  private int _n1fft,_n2fft,_n1f,_n2f;
  private FftReal _fft1;
  private FftComplex _fft2;
  private float[][] _h;

  private static double sinc(double x) {
    return (x!=0.0)?sin(x)/x:1.0;
  }

  // The Hamming window for odd window length m+1+m.
  private static double hamming(int m, double x) {
    return 0.54+0.46*cos(PI*x/m);
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    int n1 = 51;
    int n2 = 52;
    float[][] y = new float[n2][n1];
    float[][] x = new float[n2][n1]; 
    x[n2/2][n1/2] = 1.0f;
    //float[] f1 = {0.0f, 0.0f, 0.5f, 0.5f, 0.0f}; // all-pass
    //float[] f2 = {0.5f,-0.5f,-0.5f, 0.5f, 0.5f};
    //float[] f1 = {0.0f, 0.0f, 0.5f, 0.5f, 0.0f}; // quadrant
    //float[] f2 = {0.5f, 0.0f, 0.0f, 0.5f, 0.5f};
    float[] f1 = {0.0f, 0.4f, 0.5f, 0.5f, 0.0f}; // dip filter
    float[] f2 = {0.0f,-0.5f,-0.5f,-0.4f, 0.0f};
    int m1 = 10;
    int m2 = 10;
    PolygonPassFilter ppf = new PolygonPassFilter(n1,n2,m1,m2,f1,f2);
    ppf.apply(x,y);
    //dump(y);
    trace("y: min="+min(y)+" max="+max(y));
    edu.mines.jtk.mosaic.SimplePlot.asPixels(y);
  }
}
