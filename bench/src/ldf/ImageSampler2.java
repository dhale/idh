package ldf;

import static edu.mines.jtk.util.MathPlus.*;

/**
 * A 2-D image sampler.
 * Here, an image is a sampled function of two (x1,x2) coordinates.
 * It is represented by an array[n2][n1] of floats. (Specifically,
 * the x1 dimension is fastest, and the x2 dimension is slowest.) Given 
 * such an image, an image sampler determines which samples correspond 
 * to various geometric primitives, such as points, lines, etc.
 * It provides various methods to get/set those image sample values.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 12/18/2001
 * @version 09/06/2007
 */
public class ImageSampler2 {

  /**
   * Constructs an image sampler for the specified image.
   * @param image the image.
   */
  public ImageSampler2(float[][] image) {
    _image = image;
    _n1 = _image[0].length;
    _n2 = _image.length;
  }

  /**
   * Determines whether the specified sample value equals the null value.
   * @param v the image sample value.
   * @return true, if equal to the null value; false, otherwise.
   */
  public final boolean isNull(float v) {
    return v==_vnull;
  }

  /**
   * Gets the null image sample value.
   * This sample value corresponds to samples that lie outside
   * the bounds of the image or some geometric primitive.
   * @return the null value.
   */
  public float getNullValue() {
    return _vnull;
  }

  /**
   * Sets the null image sample value.
   * This sample value corresponds to samples that lie outside
   * the bounds of the image or some geometric primitive.
   * @param vnull the null value.
   */
  public void setNullValue(float vnull) {
    _vnull = vnull;
  }

  /**
   * Returns the image sample value corresponding to a single sample.
   * @param i1 x1 index of sample.
   * @param i2 x2 index of sample.
   * @return the image sample value;
     *  null, if the sample lies outside the image bounds.
   */
  public final float get(int i1, int i2) {
    try {
      return _image[i2][i1];
    } catch (ArrayIndexOutOfBoundsException e) {
      return _vnull;
    }
  }

  /**
   * Returns the image sample value corresponding to a single point.
   * @param x1 x1 coordinate of point.
   * @param x2 x2 coordinate of point.
   * @return the image sample value;
     *  null, if the sample lies outside the image bounds.
   */
  public final float get(float x1, float x2) {
    int i1 = round(x1);
    int i2 = round(x2);
    return get(i1,i2);
  }

  /**
   * Sets one image sample value.
   * Does nothing if the sample lies outside the image bounds.
   * @param i1 x1 index of sample.
   * @param i2 x2 index of sample.
   * @param v image sample value.
   */
  public final void set(int i1, int i2, float v) {
    try {
      _image[i2][i1] = v;
    } catch (ArrayIndexOutOfBoundsException e) {
    }
  }

  /**
   * Sets one image sample value.
   * Does nothing if the sample lies outside the image bounds.
   * @param x1 x1 coordinate of point.
   * @param x2 x2 coordinate of point.
   * @param v image sample value.
   */
  public final void set(float x1, float x2, float v) {
    int i1 = round(x1);
    int i2 = round(x2);
    set(i1,i2,v);
  }

  /**
   * Returns the image samples corresponding to a line segment.
   * @param x1a x1 coordinate of point a.
   * @param x2a x2 coordinate of point a.
   * @param x1b x1 coordinate of point b.
   * @param x2b x2 coordinate of point b.
   */
  public final Line sampleLine(
    float x1a, float x2a, 
    float x1b, float x2b)
  {
    return new Line(x1a,x2a,x1b,x2b);
  }

  /**
   * Image samples corresponding to a line segment.
   * On a line, the image is a function of one (r) coordinate.
   * Image (x1,x2) coordinates are projected onto the r axis, so that
   * each integer value of r corresponds to exactly one image sample.
   */
  public class Line {

    /** 
     * Returns the r coordinate of point a.
     * @return the r coordinate of point a.
     */
    public float ra() {
      return _ra;
    }

    /** 
     * Returns the r coordinate of point b.
     * @return the r coordinate of point b.
     */
    public float rb() {
      return _rb;
    }

    /**
     * Returns the length of this line.
     * @return the length.
     */
    public float length() {
      return _length;
    }

    /**
     * Returns the length of this line's projection onto the r axis.
     * @return the length.
     */
    public float lengthProjected() {
      return _lengthProjected;
    }

    /**
     * Returns the number of image samples spanned by this line.
     * Use this number as a bound in loops over image samples.
     * @return the number of samples (possibly including nulls).
     */
    public int nr() {
      return _nr;
    }

    /**
     * Converts sample index to image coordinate x1.
     * @param ir sample index.
     * @return the image coordinate x1.
     */
    public float x1(int ir) {
      return x1((float)ir);
    }

    /**
     * Converts sample coordinate to image coordinate x1.
     * @param r sample coordinate.
     * @return the image coordinate x1.
     */
    public float x1(float r) {
      return _x1q+_x1r*r;
    }

    /**
     * Converts sample index to image coordinate x2.
     * @param ir sample index.
     * @return the image coordinate x2.
     */
    public float x2(int ir) {
      return x2((float)ir);
    }

    /**
     * Converts sample coordinate to image coordinate x2.
     * @param r sample coordinate.
     * @return the image coordinate x2.
     */
    public float x2(float r) {
      return _x2q+_x2r*r;
    }

    /**
     * Returns an image sample value.
     * @param ir sample index.
     * @return image sample value; 
     *  null, if the sample lies outside the line or image bounds.
     */
    public float get(int ir) {
      return get((float)ir);
    }

    /**
     * Returns an image sample value.
     * @param r sample coordinate.
     * @return image sample value.
     *  null, if the sample lies outside the line or image bounds.
     */
    public float get(float r) {
      if (_rmin<=r && r<=_rmax) {
        return ImageSampler2.this.get(x1(r),x2(r));
      } else {
        return _vnull;
      }
    }
    
    /**
     * Returns an array of image sample values.
     * @return array[nr] of image sample values;
     *  values are null for samples outside the line or image bounds.
     */
    public float[] get() {
      float[] v = new float[_nr];
      for (int ir=0; ir<_nr; ++ir) {
        v[ir] = get(ir);
      }
      return v;
    }

    /**
     * Sets one image sample value.
     * Does nothing if the sample lies outside the line or image bounds.
     * @param ir sample index.
     * @param v image sample value.
     */
    public void set(int ir, float v) {
      set((float)ir,v);
    }

    /**
     * Sets one image sample value.
     * Does nothing if the sample lies outside the line or image bounds.
     * @param r sample coordinate.
     * @param v image sample value.
     */
    public void set(float r, float v) {
      if (_rmin<=r && r<=_rmax) {
        ImageSampler2.this.set(x1(r),x2(r),v);
      }
    }

    /**
     * Sets an array of image sample values.
     * @param v array[nr] of image sample values.
     */
    public void set(float[] v) {
      for (int ir=0; ir<_nr; ++ir) {
        set(ir,v[ir]);
      }
    }

    private float _ra,_rb,_rmin,_rmax;
    private int _nr;
    private float _x1q,_x2q;
    private float _x1r,_x2r;
    private float _length,_lengthProjected;

    /** 
     * Constructs image samples corresponding to the line segment ab.
     * @param x1a x1 coordinate of point a.
     * @param x2a x2 coordinate of point a.
     * @param x1b x1 coordinate of point b.
     * @param x2b x2 coordinate of point b.
     */
    private Line(
      float x1a, float x2a,
      float x1b, float x2b)
    {
      float x1min = min(x1a,x1b);
      float x2min = min(x2a,x2b);
      float x1max = max(x1a,x1b);
      float x2max = max(x2a,x2b);
      float x1dif = x1max-x1min;
      float x2dif = x2max-x2min;
      if (x1dif>x2dif) {
        _ra = x1a-x1min;
        _rb = x1b-x1min;
        _nr = 1+(int)x1dif;
      } else {
        _ra = x2a-x2min;
        _rb = x2b-x2min;
        _nr = 1+(int)x2dif;
      }
      _rmin = min(_ra,_rb);
      _rmax = max(_ra,_rb);
      float rf = 1.0f/(_rb-_ra);
      _x1q = (_rb*x1a-_ra*x1b)*rf;
      _x2q = (_rb*x2a-_ra*x2b)*rf;
      _x1r = (x1b-x1a)*rf;
      _x2r = (x2b-x2a)*rf;
      _length = sqrt(x1dif*x1dif+x2dif*x2dif);
      _lengthProjected = _rmax-_rmin;
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2;
  private float[][] _image;
  private float _vnull = FLT_MAX;
}
