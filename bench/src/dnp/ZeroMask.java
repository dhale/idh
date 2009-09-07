package dnp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A mask for samples that are zero or near zero.
 * Values in the mask image are either true or false.
 * Samples for which the mask is false may be unreliable
 * in some applications.
 * For example, where the mask is false, we might set structure 
 * tensors to default tensors that represent horizontal layering.
 * <p> 
 * Algorithm:
 * (1) Compute global mean absolute amplitude (gabs) for the entire image.
 * (2) Compute local mean absolute amplitude (labs) in Gaussian windows.
 * Mask is zero for all samples where labs is less than small*gabs.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.04
 */
public class ZeroMask {

  /**
   * Constructs a zero mask for a 2D image.
   * @param small small value; zeros in mask where labs &lt; small*gabs.
   * @param sigma1 Gaussian window half-width for 1st dimension.
   * @param sigma2 Gaussian window half-width for 2nd dimension.
   * @param x array of image values from which mask is derived.
   */
  public ZeroMask(
    double small, double sigma1, double sigma2,
    float[][] x) 
  {
    _n1 = x[0].length;
    _n2 = x.length;
    float[][] t = abs(x);
    float a = ((sum(t)/_n1)/_n2); // global mean absolute amplitude
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(sigma1);
    RecursiveGaussianFilter rgf2 = new RecursiveGaussianFilter(sigma2);
    float[][] b = zerofloat(_n1,_n2);
    rgf1.apply0X(t,b);
    rgf2.applyX0(b,t);
    _mask2 = new boolean[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (t[i2][i1]<small*a) {
          _mask2[i2][i1] = false;
        } else {
          _mask2[i2][i1] = true;
        }
      }
    }
  }

  /**
   * Constructs a zero mask from specified array of floats.
   * Mask is true for all non-zero samples in the array.
   * @param x array of values from which mask is derived.
   */
  public ZeroMask(float[][] x) {
    _n1 = x[0].length;
    _n2 = x.length;
    _mask2 = new boolean[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (x[i2][i1]!=0.0f)
          _mask2[i2][i1] = true;
      }
    }
  }

  /**
   * Constructs a zero mask for a 3D image.
   * @param small small value; zeros in mask where labs &lt; small*gabs.
   * @param sigma1 Gaussian window half-width for 1st dimension.
   * @param sigma2 Gaussian window half-width for 2nd dimension.
   * @param sigma3 Gaussian window half-width for 3rd dimension.
   * @param x array of image values from which mask is derived.
   */
  public ZeroMask(
    double small, double sigma1, double sigma2, double sigma3,
    float[][][] x) 
  {
    _n1 = x[0][0].length;
    _n2 = x[0].length;
    _n3 = x.length;
    float[][][] t = abs(x);
    float a = ((sum(t)/_n1)/_n2)/_n3; // global mean absolute amplitude
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(sigma1);
    RecursiveGaussianFilter rgf2 = new RecursiveGaussianFilter(sigma2);
    RecursiveGaussianFilter rgf3 = new RecursiveGaussianFilter(sigma3);
    float[][][] b = zerofloat(_n1,_n2,_n3);
    rgf1.apply0XX(t,b);
    rgf2.applyX0X(b,t);
    rgf3.applyXX0(t,b); // local mean absolute amplitude
    _mask3 = new boolean[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (b[i3][i2][i1]<small*a) {
            _mask3[i3][i2][i1] = false;
          } else {
            _mask3[i3][i2][i1] = true;
          }
        }
      }
    }
  }

  /**
   * Constructs a zero mask from specified array of floats.
   * Mask is true for all non-zero samples in the array.
   * @param x array of values from which mask is derived.
   */
  public ZeroMask(float[][][] x) {
    _n1 = x[0][0].length;
    _n2 = x[0].length;
    _n3 = x.length;
    _mask3 = new boolean[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (x[i3][i2][i1]!=0.0f)
            _mask3[i3][i2][i1] = true;
        }
      }
    }
  }

  /**
   * Fills an array of floats representing this mask.
   * The returned array has values 0.0f and 1.0f.
   * @param f array of floats.
   */
  public void getAsFloats(float[][] f) {
    Check.state(_mask2!=null,"mask constructed for a 2D image");
    for (int i2=0; i2<_n2; ++i2)
      for (int i1=0; i1<_n1; ++i1)
        f[i2][i1] = (_mask2[i2][i1])?1.0f:0.0f;
  }

  /**
   * Fills an array of floats representing this mask.
   * The returned array has values 0.0f and 1.0f.
   * @param f array of floats.
   */
  public void getAsFloats(float[][][] f) {
    Check.state(_mask3!=null,"mask constructed for a 3D image");
    for (int i3=0; i3<_n3; ++i3)
      for (int i2=0; i2<_n2; ++i2)
        for (int i1=0; i1<_n1; ++i1)
          f[i3][i2][i1] = (_mask3[i3][i2][i1])?1.0f:0.0f;
  }

  /**
   * Applies this mask to a specified array of values.
   * @param vfalse value to use where mask is false.
   * @param v array of values to be masked.
   */
  public void apply(float vfalse, float[][] v) {
    Check.state(_mask2!=null,"mask constructed for a 2D image");
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (!_mask2[i2][i1])
          v[i2][i1] = vfalse;
      }
    }
  }

  /**
   * Applies this mask to a specified array of values.
   * @param vfalse value to use where mask is false.
   * @param v array of values to be masked.
   */
  public void apply(float vfalse, float[][][] v) {
    Check.state(_mask3!=null,"mask constructed for a 3D image");
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (!_mask3[i3][i2][i1])
            v[i3][i2][i1] = vfalse;
        }
      }
    }
  }

  /**
   * Applies this mask to a specified eigentensor field.
   * @param efalse eigentensor to use where mask is false.
   * @param e eigentensors to be masked.
   */
  public void apply(float[] efalse, EigenTensors2 e) {
    Check.state(_mask2!=null,"mask constructed for a 2D image");
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (!_mask2[i2][i1])
          e.setTensor(i1,i2,efalse);
      }
    }
  }

  /**
   * Applies this mask to a specified eigentensor field.
   * @param efalse eigentensor to use where mask is false.
   * @param e eigentensors to be masked.
   */
  public void apply(float[] efalse, EigenTensors3 e) {
    Check.state(_mask3!=null,"mask constructed for a 3D image");
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (!_mask3[i3][i2][i1])
            e.setTensor(i1,i2,i3,efalse);
        }
      }
    }
  }

  private int _n1,_n2,_n3;
  private boolean[][] _mask2;
  private boolean[][][] _mask3;
}
