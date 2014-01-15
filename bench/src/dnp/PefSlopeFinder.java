/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Finds slopes of image features using lateral prediction-error filters.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.01.14
 */
public class PefSlopeFinder {

  /**
   * Constructs a slope finder with specified bounds.
   * @param pmin minimum slope.
   * @param pmax maximum slope.
   */
  public PefSlopeFinder(double pmin, double pmax) {
    _pmin = (float)pmin;
    _pmax = (float)pmax;
  }

  /**
   * Sets smoothness of updates to slopes.
   * The default smoothness is 1.0 for all dimensions.
   * @param smooth1 smoothness in 1st dimension.
   * @param smooth2 smoothness in 2nd dimension.
   */
  public void setSmoothness(double smooth1, double smooth2) {
    _eps1 = (float)smooth1;
    _eps2 = (float)smooth2;
  }

  /**
   * Applies the prediction-error filter for specified slopes.
   * @param p array of slopes.
   * @param f array input to prediction-error filter.
   * @return array of prediction errors.
   */
  public float[][] applyPef(float[][] p, float[][] f) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[][] e = new float[n2][n1];
    float[] fm = new float[n1];
    float[] tm = new float[n1];
    float[] fp = new float[n1];
    float[] tp = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      if (i2>0) { // if possible, use previous to predict current
        for (int i1=0; i1<n1; ++i1)
          tm[i1] = i1-p[i2][i1];
        _si.interpolate(n1,1.0,0.0,f[i2-1],n1,tm,fm);
      } else {
        copy(f[i2],fm);
      }
      if (i2<n2-1) { // if possible, use next to predict current
        for (int i1=0; i1<n1; ++i1)
          tp[i1] = i1+p[i2][i1];
        _si.interpolate(n1,1.0,0.0,f[i2+1],n1,tp,fp);
      } else {
        copy(f[i2],fp);
      }
      for (int i1=0; i1<n1; ++i1) // compute prediction error
        e[i2][i1] = f[i2][i1]-0.5f*(fm[i1]+fp[i1]);
    }
    return e;
  }

  /**
   * Returns an array of slopes that minimize prediction errors.
   * @param f array input to prediction-error filter.
   * @return array of slopes.
   */
  public float[][] findSlopes(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] p = new float[n2][n1];
    updateSlopes(f,p);
    return p;
  }

  /**
   * Updates specified slopes to minimize prediction errors.
   * @param f array input to prediction-error filter.
   * @param p input/output array of slopes to be updated.
   */
  public void updateSlopes(float[][] f, float[][] p) {
    int n1 = f[0].length;
    int n2 = f.length;
    M2 m2 = new M2(_eps1,_eps2);
    A2 a2 = new A2(_eps1,_eps2,p,f);
    VecArrayFloat2 vb = a2.makeRhs(p,f);
    VecArrayFloat2 vdp = new VecArrayFloat2(n1,n2);
    CgSolver cs = new CgSolver(0.01,100);
    cs.solve(a2,m2,vb,vdp);
    float[][] dp = vdp.getArray();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float pi = p[i2][i1]+dp[i2][i1];
        if (pi<_pmin) pi = _pmin;
        if (pi>_pmax) pi = _pmax;
        p[i2][i1] = pi;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _pmin,_pmax;
  private float _eps1 = 1.0f;
  private float _eps2 = 1.0f;
  private SincInterp _si = new SincInterp();
  private RecursiveGaussianFilter _rgf = new RecursiveGaussianFilter(1.0);

  /**
   * SPD linear operator A in equations to be solved for slope updates dp.
   */
  private class A2 implements CgSolver.A {
    public A2(float e1, float e2, float[][] p, float[][] f) {
      int n1 = f[0].length;
      int n2 = f.length;
      float[][] a = new float[n2][n1];
      float[] tm = new float[n1];
      float[] fm = new float[n1];
      float[] am = new float[n1];
      float[] tp = new float[n1];
      float[] fp = new float[n1];
      float[] ap = new float[n1];
      for (int i2=0; i2<n2; ++i2) {
        if (i2>0) {
          _rgf.apply1(f[i2-1],fm);
          for (int i1=0; i1<n1; ++i1)
            tm[i1] = i1-p[i2][i1];
          _si.interpolate(n1,1.0,0.0,fm,n1,tm,am);
        } else {
          zero(am);
        }
        if (i2<n2-1) {
          _rgf.apply1(f[i2+1],fp);
          for (int i1=0; i1<n1; ++i1)
            tp[i1] = i1+p[i2][i1];
          _si.interpolate(n1,1.0,0.0,fp,n1,tp,ap);
        } else {
          zero(ap);
        }
        for (int i1=0; i1<n1; ++i1)
          a[i2][i1] = 0.5f*(ap[i1]-am[i1]);
      }
      _a = a;
      _e1 = e1;
      _e2 = e2;
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      int n1 = x[0].length;
      int n2 = x.length;
      zero(y);
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ai = _a[i2][i1];
          float x0 = x[i2][i1];
          float x1 = x[i2][i1-1];
          float x2 = x[i2-1][i1];
          float y0 = ai*x0;
          float y1 = _e1*(x0-x1);
          float y2 = _e2*(x0-x2);
          y[i2][i1] += ai*y0+_e1*y1+_e2*y2;
          y[i2][i1-1] -= _e1*y1;
          y[i2-1][i1] -= _e2*y2;
        }
      }
      for (int i1=0,i2=1; i2<n2; ++i2) {
        float ai = _a[i2][i1];
        float x0 = x[i2][i1];
        float x2 = x[i2-1][i1];
        float y0 = ai*x0;
        float y2 = _e2*(x0-x2);
        y[i2][i1] += ai*y0+_e2*y2;
        y[i2-1][i1] -= _e2*y2;
      }
      for (int i1=1,i2=0; i1<n1; ++i1) {
        float ai = _a[i2][i1];
        float x0 = x[i2][i1];
        float x1 = x[i2][i1-1];
        float y0 = ai*x0;
        float y1 = _e1*(x0-x1);
        y[i2][i1] += ai*y0+_e1*y1;
        y[i2][i1-1] -= _e1*y1;
      }
    }
    public VecArrayFloat2 makeRhs(float[][] p, float[][] f) {
      int n1 = p[0].length;
      int n2 = p.length;
      float[][] e = new float[n2][n1];
      _rgf.apply0X(f,e);
      e = applyPef(p,e);
      float[][] y = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float ai = _a[i2][i1];
          float y0 = e[i2][i1];
          y[i2][i1] += _a[i2][i1]*e[i2][i1];
        }
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
    M2(float eps1, float eps2) {
      float sigma1 = sqrt(eps1);
      float sigma2 = sqrt(eps2);
      _ref = new RecursiveExponentialFilter(sigma1,sigma2);
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      _ref.apply1(x,y);
      _ref.apply2(y,y);
      _ref.apply2(y,y);
      _ref.apply1(y,y);
    }
    private RecursiveExponentialFilter _ref;
  }
}
