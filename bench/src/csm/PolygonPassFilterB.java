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
public class PolygonPassFilterB {

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
   * @param f1 array of frequencies (cycles/sample) for 1st dimension.
   * @param f2 array of frequencies (cycles/sample) for 2nd dimension.
   */
  public PolygonPassFilterB(
    int n1, int n2, float[] f1, float[] f2) 
  {
    Check.argument(f1.length==f2.length,"length of f1 == length of f2");
    Check.argument(f1[0]==f1[f1.length-1],"first f1 == last f1");
    Check.argument(f2[0]==f2[f2.length-1],"first f2 == last f2");
    _n1 = n1;
    _n2 = n2;

    // FFTs for 1st and 2nd dimensions.
    _n1fft = FftReal.nfftSmall(n1);
    _n2fft = FftComplex.nfftSmall(n2);
    _fft1 = new FftReal(_n1fft);
    _fft2 = new FftComplex(_n2fft);

    // Sampling of frequencies in 1st and 2nd dimensions.
    _n1f = 1+_n1fft/2;
    _n2f = _n2fft;
    float d1f = 1.0f/_n1fft;
    float d2f = 1.0f/_n2fft;

    // Inflate polygon slightly.
    f1 = mul(1.0f+FLT_EPSILON,f1);
    f2 = mul(1.0f+FLT_EPSILON,f2);

    // Compute Fourier transform H(k1,k2), non-zero inside polygon.
    _h = new float[_n2f][_n1f];
    float hs = 1.0f/_n1fft/_n2fft;
    for (int i2=0; i2<_n2f; ++i2) {
      float k2 = i2*d2f;
      if (k2>=0.5f) k2 -= 1.0f;
      for (int i1=0; i1<_n1f; ++i1) {
        float k1 = i1*d1f;
        if (pointInsidePolygon(k1,k2,f1,f2)) {
          _h[i2][i1] = hs;
        }
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
    _fft1.realToComplex1(-1,n2,xy,xy);
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
    _fft1.complexToReal1(1,n2,xy,xy);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        y[i2][i1] = xy[i2][i1];
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2;
  private float[] _f1,_f2;
  private int _n1fft,_n2fft,_n1f,_n2f;
  private FftReal _fft1;
  private FftComplex _fft2;
  private float[][] _h;

  private static boolean pointInsidePolygon(
    float px, float py, float[] vx, float[] vy)
  {
    int n = vx.length;
    boolean inside = false;
    for (int i=0,j=n-1; i<n; j=i++) {
      if (((vy[i]>py)!=(vy[j]>py)) &&
          (px<(vx[j]-vx[i])*(py-vy[i])/(vy[j]-vy[i])+vx[i]))
        inside = !inside;
    }
    return inside;
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
    PolygonPassFilterB ppf = new PolygonPassFilterB(n1,n2,f1,f2);
    ppf.apply(x,y);
    //dump(y);
    trace("y: min="+min(y)+" max="+max(y));
    edu.mines.jtk.mosaic.SimplePlot.asPixels(y);
  }
}
