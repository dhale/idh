/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Fomel's plane-wave destruction filter.
 * A plane-wave destruction (PWD) filter (Fomel, 2002) uses lateral
 * derivatives to attenuate locally linear or planar features in images. The
 * slopes of such features are the filter parameters, and these slopes can be
 * estimated by minimizing a sum of squared PWD output values. Indeed, one of
 * the most useful applications of PWD is that of estimating slopes of locally
 * linear or planar image features.
 * <p>
 * We must typically constrain the slope estimated for one image sample to be
 * similar to slopes estimated for adjacent samples. Therefore, while
 * minimizing PWD output values, we also seek to minimize sample-to-sample
 * variations in slopes. The balance between these two goals is controlled by
 * smoothness parameters that may be specified for each image dimension.
 * <p>
 * As described by Fomel (2002), the requirement that slopes vary smoothly,
 * and the fact that the output of a plane-wave destruction filter is a
 * non-linear function of those slopes, leads to the use of an iterative
 * conjugate-gradient method for estimating slopes. One difference between
 * this implementation and that described by Fomel is that this class uses the
 * smoothness parameters to control directly the smoothness of the final
 * estimated slopes p, whereas Fomel describes using smoothness parameters to
 * control the smoothness of incremental updates dp in the iterative method
 * used to estimate slopes. In practice, I find that repeatedly accumulating
 * smooth updates dp need not yield final slopes p that are smooth.
 * <p>
 * Fomel's (2002) description of PWD uses either backward or forward
 * finite-difference approximations to lateral derivatives, which means that
 * slope estimates are obtained at locations halfway between image columns.
 * This implementation provides a third centered-difference option for which
 * slope estimates are obtained precisely on image samples, and this option is
 * the default.
 * <p>
 * In all methods for this class, the units of slope are samples per sample.
 * For example, the value 0.5 would represent a slope of one-half vertical
 * sample per horizontal sample.
 * <p>
 * Reference: Fomel, S., 2002, Applications of plane-wave destruction
 * filters: Geophysics, 67, 1946–196.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.01.14
 */
public class PlaneWaveDestructor {

  /**
   * Constructs a slope finder with specified bounds.
   * @param pmin minimum slope.
   * @param pmax maximum slope.
   */
  public PlaneWaveDestructor(double pmin, double pmax) {
    _pmin = (float)pmin;
    _pmax = (float)pmax;
  }

  /**
   * Sets smoothness of updates to slopes.
   * The default smoothness is 1.0 for all dimensions, which may
   * be too small for many applications.
   * @param smooth1 smoothness in 1st dimension.
   * @param smooth2 smoothness in 2nd dimension.
   */
  public void setSmoothness(double smooth1, double smooth2) {
    _eps1 = (float)smooth1;
    _eps2 = (float)smooth2;
  }

  /**
   * Sets the lateral bias in the location of PWD slope estimates.
   * The only bias values permitted are -0.5, 0.0, and 0.5. For example,
   * a value of 0.5 means that slopes for lateral sample index i are actually
   * estimated at locations i+0.5. The default is 0.0, for no lateral bias.
   * <p>
   * In effect, this parameter determines which finite-difference
   * approximation is used for horizontal derivatives in PWD. The values -0.5,
   * 0.0, and 0.5 correspond to backward, centered, and forward differences,
   * respectively.
   */
  public void setLateralBias(double bias) {
    Check.argument(bias==-0.5 || bias==0.0 || bias==0.5,"valid bias");
    _bias = (float)bias;
  }

  /**
   * Sets the number of outer iterations used to update slope estimates.
   * These iterations are necessary because each PWD output value is a
   * non-linear function of slope. The default number is 5.
   * @param nouter the number of outer iterations.
   */
  public void setOuterIterations(int nouter) {
    _nouter = nouter;
  }

  /**
   * Sets the number of inner iterations used to update slope estimates.
   * These iterations are necessary because we seek to estimate smoothly
   * varying slopes for all image samples. The default number is 20.
   * @param ninner the number of inner iterations.
   */
  public void setInnerIterations(int ninner) {
    _ninner = ninner;
  }

  /**
   * Applies the PWD filter for specified slopes.
   * @param p array of slopes.
   * @param f input array of image samples.
   * @return array of filtered image samples.
   */
  public float[][] applyFilter(float[][] p, float[][] f) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[][] g = new float[n2][n1];
    applyFilter(p,f,g);
    return g;
  }

  /**
   * Returns an array of slopes estimated for a specified image.
   * Assumes initial slope estimates are zero.
   * @param f input array of image samples.
   * @return array of estimated slopes.
   */
  public float[][] findSlopes(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] p = new float[n2][n1];
    updateSlopes(f,p);
    return p;
  }

  /**
   * Updates slope estimates to minimize PWD filter output.
   * This method is useful if a good initial estimate of slopes is available.
   * @param f input array of image samples.
   * @param p input/output array of slopes to be updated.
   */
  public void updateSlopes(float[][] f, float[][] p) {
    int n1 = f[0].length;
    int n2 = f.length;

    // Normalize input so that the local rms is one, so that the smoothness of
    // slope updates will be independent of the local amplitudes of image
    // features. Without this normalization, slopes in low-amplitude parts of
    // an image would tend to be smoother than in high-amplitude parts.
    f = divideByRms(16.0,16.0,f);

    // Preconditioner for conjugate-gradient iterations.
    M2 m2 = new M2(_eps1,_eps2);

    // Work arrays.
    float[][] g = new float[n2][n1]; // PWD output array
    float[][] p0 = new float[n2][n1]; // initial slopes
    float[][] dp = new float[n2][n1]; // increments for slopes
    VecArrayFloat2 vdp = new VecArrayFloat2(dp); // in a vector for CG

    // Outer loop over linear approximations.
    for (int iouter=0; iouter<_nouter; ++iouter) {

      // The current slopes p become the initial slopes p0.
      copy(p,p0);

      // Norm of PWD output computed for initial slopes p0.
      double g0norm = normPwd(p0,f,g);

      // SPD linear operator A in the system A dp = b.
      A2 a2 = new A2(_eps1,_eps2,p0,f);

      // Right-hand side vector b in the system A dp = b.
      VecArrayFloat2 vb = a2.makeRhs(p0,f);

      // Solve A dp = b using preconditioned conjugate-gradient method.
      // We use only a small number of CG iterations because we will
      CgSolver cs = new CgSolver(0.01,_ninner);
      cs.solve(a2,m2,vb,vdp);

      // Because PWD output is not a linear function of slopes, the update dp
      // might actually increase the norm of PWD output. Therefore, here we
      // loop over decreasing step sizes until we decrease that norm. Note
      // that we give up after a fixed number of step sizes. The maximum of 8
      // steps here is consistent with Fomel's implementation.
      float step = 1.0f;
      double gnorm = 2.0*g0norm;
      for (int istep=0; istep<8 && g0norm<=gnorm; ++istep,step*=0.5f) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float pi = p0[i2][i1]+step*dp[i2][i1];
            if (pi<_pmin) pi = _pmin;
            if (pi>_pmax) pi = _pmax;
            p[i2][i1] = pi;
          }
        }
        gnorm = normPwd(p,f,g);
      }
    }
  }

  /**
   * Applies the PWD filter for specified slopes.
   * @param p array of slopes.
   * @param f input array of image samples.
   * @param g output array of filtered image samples.
   */
  public void applyFilter(float[][] p, float[][] f, float[][] g) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[] fm = new float[n1];
    float[] tm = new float[n1];
    float[] fp = new float[n1];
    float[] tp = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = (_bias== 0.5f || i2==0   )?i2:i2-1;
      int i2p = (_bias==-0.5f || i2==n2-1)?i2:i2+1;
      if (i2m<i2p) {
        float dx2 = 0.5f*(i2p-i2m);
        for (int i1=0; i1<n1; ++i1) {
          tm[i1] = i1-p[i2][i1]*dx2;
          tp[i1] = i1+p[i2][i1]*dx2;
        }
        _si.interpolate(n1,1.0,0.0,f[i2m],n1,tm,fm);
        _si.interpolate(n1,1.0,0.0,f[i2p],n1,tp,fp);
        float scale = 0.5f/dx2;
        for (int i1=0; i1<n1; ++i1)
          g[i2][i1] = scale*(fp[i1]-fm[i1]);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _pmin;
  private float _pmax;
  private float _bias = 0.0f;
  private float _eps1 = 1.0f;
  private float _eps2 = 1.0f;
  private int _nouter = 5;
  private int _ninner = 20;
  private SincInterpolator _si = new SincInterpolator();
  private RecursiveGaussianFilter _rgf = new RecursiveGaussianFilter(1.0);

  /**
   * SPD linear operator A in equations to be solved for slope updates dp.
   */
  private class A2 implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      int n1 = x[0].length;
      int n2 = x.length;
      zero(y);
      for (int i2=1; i2<n2; ++i2) { // for most image samples, ...
        for (int i1=1; i1<n1; ++i1) {
          float ai = _a[i2][i1];
          float x0 = x[i2][i1];
          float x1 = x[i2][i1-1];
          float x2 = x[i2-1][i1];
          float y0 = ai*x0;                 // gather from inputs
          float y1 = _e1*(x0-x1);
          float y2 = _e2*(x0-x2);
          y[i2][i1] += ai*y0+_e1*y1+_e2*y2; // scatter to outputs
          y[i2][i1-1] -= _e1*y1;
          y[i2-1][i1] -= _e2*y2;
        }
      }
      for (int i1=0,i2=1; i2<n2; ++i2) { // special case: i1 = 0
        float ai = _a[i2][i1];
        float x0 = x[i2][i1];
        float x2 = x[i2-1][i1];
        float y0 = ai*x0;
        float y2 = _e2*(x0-x2);
        y[i2][i1] += ai*y0+_e2*y2;
        y[i2-1][i1] -= _e2*y2;
      }
      for (int i1=1,i2=0; i1<n1; ++i1) { // special case: i2 = 0
        float ai = _a[i2][i1];
        float x0 = x[i2][i1];
        float x1 = x[i2][i1-1];
        float y0 = ai*x0;
        float y1 = _e1*(x0-x1);
        y[i2][i1] += ai*y0+_e1*y1;
        y[i2][i1-1] -= _e1*y1;
      }
    }
    A2(float e1, float e2, float[][] p, float[][] f) {
      int n1 = f[0].length;
      int n2 = f.length;
      float[][] a = new float[n2][n1];
      float[] tm = new float[n1];
      float[] tp = new float[n1];
      float[] fm = new float[n1];
      float[] fp = new float[n1];
      float[] am = new float[n1];
      float[] ap = new float[n1];
      for (int i2=0; i2<n2; ++i2) {
        int i2m = (_bias== 0.5f || i2==0   )?i2:i2-1;
        int i2p = (_bias==-0.5f || i2==n2-1)?i2:i2+1;
        if (i2m<i2p) {
          _rgf.apply1(f[i2m],fm);
          _rgf.apply1(f[i2p],fp);
          float dx2 = 0.5f*(i2p-i2m);
          for (int i1=0; i1<n1; ++i1) {
            tm[i1] = i1-p[i2][i1]*dx2;
            tp[i1] = i1+p[i2][i1]*dx2;
          }
          _si.interpolate(n1,1.0,0.0,fm,n1,tm,am);
          _si.interpolate(n1,1.0,0.0,fp,n1,tp,ap);
          for (int i1=0; i1<n1; ++i1)
            a[i2][i1] = 0.5f*(ap[i1]+am[i1]);
        }
      }
      _a = a;
      _e1 = e1;
      _e2 = e2;
    }
    VecArrayFloat2 makeRhs(float[][] p, float[][] f) {
      int n1 = p[0].length;
      int n2 = p.length;
      float[][] e = new float[n2][n1];
      _rgf.apply0X(f,e);
      e = applyFilter(p,e);
      float[][] y = new float[n2][n1];
      for (int i2=1; i2<n2; ++i2) { // for most image samples, ...
        for (int i1=1; i1<n1; ++i1) {
          float ai = _a[i2][i1];
          float p0 = p[i2][i1];
          float p1 = p[i2][i1-1];
          float p2 = p[i2-1][i1];
          float y0 = e[i2][i1];                 // gather from inputs
          float y1 = _e1*(p0-p1);
          float y2 = _e2*(p0-p2);
          y[i2][i1] -= ai*y0+_e1*y1+_e2*y2; // scatter to outputs
          y[i2][i1-1] += _e1*y1;
          y[i2-1][i1] += _e2*y2;
        }
      }
      for (int i1=0,i2=1; i2<n2; ++i2) { // special case: i1 = 0
        float ai = _a[i2][i1];
        float p0 = p[i2][i1];
        float p2 = p[i2-1][i1];
        float y0 = e[i2][i1];
        float y2 = _e2*(p0-p2);
        y[i2][i1] -= ai*y0+_e2*y2;
        y[i2-1][i1] += _e2*y2;
      }
      for (int i1=1,i2=0; i1<n1; ++i1) { // special case: i2 = 0
        float ai = _a[i2][i1];
        float p0 = p[i2][i1];
        float p1 = p[i2][i1-1];
        float y0 = e[i2][i1];
        float y1 = _e1*(p0-p1);
        y[i2][i1] -= ai*y0+_e1*y1;
        y[i2][i1-1] += _e1*y1;
      }
      return new VecArrayFloat2(y);
    }
    private float[][] _a;
    private float _e1,_e2;
  }

  /**
   * SPD preconditioner implemented with smoothing filters.
   */
  private class M2 implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      _ref.apply1(x,y);
      _ref.apply2(y,y);
      _ref.apply2(y,y);
      _ref.apply1(y,y);
    }
    M2(float eps1, float eps2) {
      float sigma1 = sqrt(eps1);
      float sigma2 = sqrt(eps2);
      _ref = new RecursiveExponentialFilter(sigma1,sigma2);
    }
    private RecursiveExponentialFilter _ref;
  }

  /**
   * Returns the sum of squared outputs g for the PWD filter.
   */
  private double normPwd(float[][] p, float[][] f, float[][] g) {
    int n1 = p[0].length;
    int n2 = p.length;
    applyFilter(p,f,g);
    double gnorm = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        gnorm += g[i2][i1]*g[i2][i1];
      }
    }
    return gnorm;
  }

  private float[][] divideByRms(double sigma1, double sigma2, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = mul(f,f);
    RecursiveExponentialFilter ref = 
      new RecursiveExponentialFilter(16.0,16.0);
    ref.apply(g,g);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (g[i2][i1]>0.0f) {
          g[i2][i1] = f[i2][i1]/sqrt(g[i2][i1]);
        } else {
          g[i2][i1] = 0.0f;
        }
      }
    }
    return g;
  }
}
