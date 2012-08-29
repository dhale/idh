/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;

/**
 * Gridding in 2D by tensor-guided simple kriging.
 * Gridding is interpolation onto a uniformly sampled grid of a set 
 * of specified scattered sample values. Here that interpolation is 
 * performed by simple kriging.
 * <p>
 * More precisely, let d denote a data vector with elements equal to
 * the scattered known sample values, and let m denote a model vector
 * computed by gridding according to the following equation:
 * <pre>
 * m = m0 + Cm K' inv(K Cm K' + Cd) (d - K m0). 
 * </pre>
 * Here, m0 denotes an a priori model, K is a matrix that samples the
 * gridded model space at locations where scattered data are specified, 
 * and Cm and Cd are model and covariance matrices.
 * <p>
 * If specified, the tensor field that guides the interpolation is
 * used to implement multiplication by the model covariance matrix
 * Cm. In practice, the elements of this matrix are never computed
 * explicitly, because the gridded model space is assumed to be much
 * larger than the scattered data space. Instead, multiplication by
 * Cm is implemented by tensor-guided smoothing filters.
 * <p>
 * The data covariance matrix Cd is diagonal, with variances that
 * may be specified by one constant standard deviation or by an 
 * array of standard deviations, one for each sample of data.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.03.05
 */
public class KrigingGridder2 implements Gridder2 {

  /**
   * Constructs a gridder for default isotropic constant tensors.
   */
  public KrigingGridder2() {
    this(null);
  }

  /**
   * Constructs a gridder for specified samples.
   * The specified arrays are referenced; not copied.
   * @param f array of sample values f(x1,x2).
   * @param x1 array of sample x1 coordinates.
   * @param x2 array of sample x2 coordinates.
   */
  public KrigingGridder2(float[] f, float[] x1, float[] x2) {
    setTensor(1.0,0.0,1.0);
    setScattered(f,x1,x2);
  }

  /**
   * Constructs a gridder for the specified tensors.
   * @param tensors the tensors.
   */
  public KrigingGridder2(Tensors2 tensors) {
    setTensors(tensors);
  }

  /**
   * Constructs a gridder for the specified tensors and samples.
   * The specified arrays are referenced; not copied.
   * @param tensors the tensors.
   * @param f array of sample values f(x1,x2).
   * @param x1 array of sample x1 coordinates.
   * @param x2 array of sample x2 coordinates.
   */
  public KrigingGridder2(
    Tensors2 tensors, float[] f, float[] x1, float[] x2) 
  {
    setTensors(tensors);
    setScattered(f,x1,x2);
  }

  /**
   * Sets the constant tensor D used by this gridder.
   * @param d11 tensor element D(1,1). 
   * @param d12 tensor element D(1,2). 
   * @param d22 tensor element D(2,2). 
   */
  public void setTensor(double d11, double d12, double d22) {
    Check.argument(d11*d22>d12*d12,"tensor must be positive-definite");
    _d11 = d11;
    _d12 = d12;
    _d22 = d22;
    _tensors = null;
  }

  /**
   * Sets the variable tensor field used to compute distances.
   * The default is a constant isotropic tensor.
   * @param tensors the tensors; null for default tensors.
   */
  public void setTensors(Tensors2 tensors) {
    _tensors = tensors;
  }

  /**
   * Enables use of Paciorek's (2003) approximation in tensor-guided kriging.
   * This approximation is relevant only if tensors are specified.
   * @param p true, to use Paciorek's approximation; false, otherwise.
   */
  public void setPaciorek(boolean p) {
    _paciorek = p;
  }

  /**
   * For testing, only.
   */
  public void setIdentityTensors() {
    _tensors = new Tensors2() {
      public void getTensor(int i1, int i2, float[] d) {
        d[0] = 1.0f;
        d[1] = 0.0f;
        d[2] = 1.0f;
      }
    };
  }

  /**
   * Sets the model covariance function of distance.
   * The default is a smooth covariance function with 
   * parameters shape = sigma = range = 1.
   * @param cm the model covariance function.
   */
  public void setModelCovariance(Covariance model) {
    _cm = model;
  }

  /**
   * Sets a constant standard deviation for all data errors.
   * Data errors are assumed to be uncorrelated.
   * @param sd the standard deviation = sqrt(variance).
   */
  public void setDataError(double sd) {
    _sdConstant = sd;
  }

  /**
   * Sets standard deviations for data errors.
   * Data errors are assumed to be uncorrelated.
   * @param sd array of standard deviations = sqrt(variance).
   */
  public void setDataError(float[] sd) {
    _sd = copy(sd);
  }

  /**
   * Sets the order of the polynomial trend to be fit to sample values.
   * This trend is removed before kriging and restored after kriging.
   * The default order is -1, so that no trend is removed.
   * @param order the order of the polynomial fit; must be -1, 0, 1, or 2.
   */
  public void setPolyTrend(int order) {
    Check.argument(-1<=order,"-1<=order");
    Check.argument(order<=2,"order<=2");
    if (_trend!=null) {
      _trend.restore(_f,_x1,_x2);
      _trend = null;
    }
    if (order!=-1) {
      _trend = new PolyTrend2(order,_f,_x1,_x2);
      _trend.detrend(_f,_x1,_x2);
      //dump(_f);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // interface Gridder2

  public void setScattered(float[] f, float[] x1, float[] x2) {
    _f = copy(f);
    _x1 = copy(x1);
    _x2 = copy(x2);
  }

  public float[][] grid(Sampling s1, Sampling s2) {
    checkSamplings(s1,s2);
    checkScattered();
    ensureDataErrors();
    float[][] q = null;
    if (_tensors!=null) {
      if (_paciorek) {
        q = gridForVariableTensorsP(s1,s2);
      } else {
        q = gridForVariableTensors(s1,s2);
      }
    } else {
      q = gridForConstantTensors(s1,s2);
    }
    if (_trend!=null)
      _trend.restore(q,s1,s2);
    return q;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Tensors2 _tensors; // variable tensors D; null for constant D
  private boolean _paciorek; // true, if using Paciorek's approximation
  private double _d11,_d12,_d22; // elements of constant tensor D
  private float[] _f,_x1,_x2; // scattered data
  private float[] _sd; // array of std devs for data errors
  private double _sdConstant = 0.0;  // std dev for data errors, if constant
  private Covariance _cm = new SmoothCovariance(1.0,1.0,1.0,2);
  private PolyTrend2 _trend; // polynomial trend; null, if none
  private int _order = -1; // order of poly trend; -1, if none

  private void ensureDataErrors() {
    if (_sd==null || _sd.length!=_f.length)
      _sd = fillfloat((float)_sdConstant,_f.length);
  }

  /**
   * Simple kriging for a constant tensor field, which may be anisotropic.
   */
  private float[][] gridForConstantTensors(Sampling s1, Sampling s2) {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n = _f.length;
    DMatrix cf = new DMatrix(n,1);
    DMatrix cm = new DMatrix(n,n);
    double det = _d11*_d22-_d12*_d12;
    double t11 =  _d22/det;
    double t12 = -_d12/det;
    double t22 =  _d11/det;
    for (int i=0; i<n; ++i) {
      double x1i = _x1[i];
      double x2i = _x2[i];
      for (int j=0; j<n; ++j) {
        double x1j = _x1[j];
        double x2j = _x2[j];
        double r = distance(t11,t12,t22,x1i,x2i,x1j,x2j);
        cm.set(i,j,_cm.evaluate(r));
      }
      cm.set(i,i,cm.get(i,i)+_sd[i]*_sd[i]);
      cf.set(i,0,_f[i]);
    }
    DMatrix cw = cm.solve(cf);
    float[][] q = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      double x2i = s2.getValue(i2);
      for (int i1=0; i1<n1; ++i1) {
        double x1i = s1.getValue(i1);
        double qi = 0.0f;
        for (int j=0; j<n; ++j) {
          double x1j = _x1[j];
          double x2j = _x2[j];
          double r = distance(t11,t12,t22,x1i,x2i,x1j,x2j);
          qi += cw.get(j,0)*_cm.evaluate(r);
        }
        q[i2][i1] = (float)qi;
      }
    }
    return q;
  }

  private static double distance(
    double t11, double t12, double t22,
    double x1a, double x2a, 
    double x1b, double x2b) 
  {
    double dx1 = x1a-x1b;
    double dx2 = x2a-x2b;
    return sqrt(dx1*(t11*dx1+t12*dx2)+dx2*(t12*dx1+t22*dx2));
  }

  private static double evaluatePaciorek(
    Sampling s1, Sampling s2, Tensors2 tensors, Covariance cm,
    double d11i, double d12i, double d22i, double deti, float[] d,
    double x1i, double x2i, double x1j, double x2j)
  {
    int j1 = s1.indexOfNearest(x1j);
    int j2 = s2.indexOfNearest(x2j);
    tensors.getTensor(j1,j2,d);
    double d11j = d[0], d12j = d[1], d22j = d[2];
    double detj = d11j*d22j-d12j*d12j;
    double d11 = 0.5*(d11i+d11j);
    double d12 = 0.5*(d12i+d12j);
    double d22 = 0.5*(d22i+d22j);
    double det = d11*d22-d12*d12;
    double t11 =  d22/det;
    double t12 = -d12/det;
    double t22 =  d11/det;
    double s = sqrt(sqrt(deti*detj)/det);
    double r = distance(t11,t12,t22,x1i,x2i,x1j,x2j);
    return s*cm.evaluate(r);
  }

  /**
   * Paciorek's approximation to tensor-guided kriging.
   */
  private float[][] gridForVariableTensorsP(Sampling s1, Sampling s2) {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n = _f.length;
    DMatrix cf = new DMatrix(n,1);
    DMatrix cm = new DMatrix(n,n);
    float[] d = new float[3];
    for (int i=0; i<n; ++i) {
      double x1i = _x1[i];
      double x2i = _x2[i];
      int i1 = s1.indexOfNearest(x1i);
      int i2 = s2.indexOfNearest(x2i);
      _tensors.getTensor(i1,i2,d);
      double d11i = d[0], d12i = d[1], d22i = d[2];
      double deti = d11i*d22i-d12i*d12i;
      for (int j=0; j<n; ++j) {
        double x1j = _x1[j];
        double x2j = _x2[j];
        double cij = evaluatePaciorek(s1,s2,_tensors,_cm,
                                      d11i,d12i,d22i,deti,d,
                                      x1i,x2i,x1j,x2j);
        cm.set(i,j,cij);
      }
      cm.set(i,i,cm.get(i,i)+_sd[i]*_sd[i]);
      cf.set(i,0,_f[i]);
    }
    DMatrix cw = cm.solve(cf);
    float[][] q = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      double x2i = s2.getValue(i2);
      for (int i1=0; i1<n1; ++i1) {
        double x1i = s1.getValue(i1);
        _tensors.getTensor(i1,i2,d);
        double d11i = d[0], d12i = d[1], d22i = d[2];
        double deti = d11i*d22i-d12i*d12i;
        double qi = 0.0f;
        for (int j=0; j<n; ++j) {
          double x1j = _x1[j];
          double x2j = _x2[j];
          double cij = evaluatePaciorek(s1,s2,_tensors,_cm,
                                        d11i,d12i,d22i,deti,d,
                                        x1i,x2i,x1j,x2j);
          qi += cij*cw.get(j,0);
        }
        q[i2][i1] = (float)qi;
      }
    }
    return q;
  }

  /**
   * Tensor-guided kriging for only smooth model covariance.
   * Uses Paciorek's approximation as a preconditioner in
   * an iterative CG computation of the kriging weights.
   */
  private float[][] gridForVariableTensors(Sampling s1, Sampling s2) {
    Check.state(_cm instanceof SmoothCovariance,
      "model covariance is a SmoothCovariance");
    SmoothCovariance scm = (SmoothCovariance)_cm;
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float[][] fxs = SimpleGridder2.samplesOnGrid(s1,s2,_f,_x1,_x2,_sd);
    float[] f = fxs[0];
    float[] x1 = fxs[1];
    float[] x2 = fxs[2];
    float[] sd = fxs[3];
    int n = f.length;
    float[] af = f;
    float[] az = new float[n];
    VecArrayFloat1 vf = new VecArrayFloat1(af);
    VecArrayFloat1 vz = new VecArrayFloat1(az);
    SmoothA a = new SmoothA(x1,x2,s1,s2,scm,sd);
    SmoothM m = new SmoothM(x1,x2,s1,s2,scm,sd);
    CgSolver cgs = new CgSolver(0.001,100);
    CgSolver.Info info = cgs.solve(a,m,vf,vz);
    //CgSolver.Info info = cgs.solve(a,vf,vz);
    cgs = null;
    float[][] q = scm.apply(s1,s2,x1,x2,az);
    return q;
  }

  private static class SmoothA implements CgSolver.A {

    SmoothA(float[] x1, float[] x2, Sampling s1, Sampling s2,
      SmoothCovariance cm, float[] sd) 
    {
      int n = x1.length;
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      _k1 = new int[n];
      _k2 = new int[n];
      for (int i=0; i<n; ++i) {
        _k1[i] = s1.indexOfNearest(x1[i]);
        _k2[i] = s2.indexOfNearest(x2[i]);
      }
      _bx = new float[n2][n1];
      _cm = cm;
      _sd = sd;
    }

    public void apply(Vec x, Vec y) {
      float[] ax = ((VecArrayFloat1)x).getArray();
      float[] ay = ((VecArrayFloat1)y).getArray();
      int n = ax.length;
      zero(_bx);
      for (int i=0; i<n; ++i)
        _bx[_k2[i]][_k1[i]] = ax[i];
      _cm.apply(_bx);
      for (int i=0; i<n; ++i)
        ay[i] = _bx[_k2[i]][_k1[i]]+_sd[i]*_sd[i]*ax[i];
    }

    private int[] _k1,_k2;
    private float[][] _bx;
    private float[] _sd;
    private SmoothCovariance _cm;
  }

  private static class SmoothM implements CgSolver.A {

    SmoothM(float[] x1, float[] x2, Sampling s1, Sampling s2,
      SmoothCovariance cm, float[] sd) 
    {
      Tensors2 tensors = cm.getTensors();
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n = x1.length;
      DMatrix am = new DMatrix(n,n);
      float[] d = new float[3];
      for (int i=0; i<n; ++i) {
        double x1i = x1[i];
        double x2i = x2[i];
        int i1 = s1.indexOfNearest(x1i);
        int i2 = s2.indexOfNearest(x2i);
        tensors.getTensor(i1,i2,d);
        double d11i = d[0], d12i = d[1], d22i = d[2];
        double deti = d11i*d22i-d12i*d12i;
        for (int j=0; j<n; ++j) {
          double x1j = x1[j];
          double x2j = x2[j];
          double cij = evaluatePaciorek(s1,s2,tensors,cm,
                                        d11i,d12i,d22i,deti,d,
                                        x1i,x2i,x1j,x2j);
          am.set(i,j,cij);
        }
        am.set(i,i,am.get(i,i)+sd[i]*sd[i]);
      }
      am = am.inverse();
      _am = new float[n][n];
      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          _am[i][j] = (float)am.get(i,j);
        }
      }
    }

    public void apply(Vec x, Vec y) {
      float[] ax = ((VecArrayFloat1)x).getArray();
      float[] ay = ((VecArrayFloat1)y).getArray();
      int n = ax.length;
      for (int i=0; i<n; ++i) {
        float ayi = 0.0f;
        for (int j=0; j<n; ++j)
          ayi += _am[i][j]*ax[j];
        ay[i] = ayi;
      }
    }

    private float[][] _am;
  }

  private void checkSamplings(Sampling s1, Sampling s2) {
    Check.argument(s1.isUniform(),"s1 is uniform");
    Check.argument(s2.isUniform(),"s2 is uniform");
  }

  private void checkScattered() {
    Check.state(_f!=null,"scattered samples have been set");
    Check.state(_x1!=null,"scattered samples have been set"); 
    Check.state(_x2!=null,"scattered samples have been set");
  }
}
