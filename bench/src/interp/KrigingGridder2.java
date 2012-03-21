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
   * The default is a Matern covariance function with 
   * parameters shape = sigma = range = 1.
   * @param cm the model covariance function.
   */
  public void setModelCovariance(Covariance model) {
    _modelC = model;
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
  private Covariance _modelC = new Matern(1.0,1.0,1.0); // model covariance
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
        cm.set(i,j,_modelC.evaluate(r));
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
          qi += cw.get(j,0)*_modelC.evaluate(r);
        }
        q[i2][i1] = (float)qi;
      }
    }
    return q;
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
        int j1 = s1.indexOfNearest(x1j);
        int j2 = s2.indexOfNearest(x2j);
        _tensors.getTensor(j1,j2,d);
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
        cm.set(i,j,s*_modelC.evaluate(r));
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
          int j1 = s1.indexOfNearest(x1j);
          int j2 = s2.indexOfNearest(x2j);
          _tensors.getTensor(j1,j2,d);
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
          qi += cw.get(j,0)*s*_modelC.evaluate(r);
        }
        q[i2][i1] = (float)qi;
      }
    }
    return q;
  }

  /**
   * Tensor-guided kriging for Matern model covariance with shape = 1. 
   * This is a tensor-guided extension of Whittle's (1954) covariance
   * for 2D models. Uses Paciorek's approximation as a preconditioner
   * in an iterative CG computation of the kriging weights.
   */
  private float[][] gridForVariableTensors(Sampling s1, Sampling s2) {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    WhittleSmoother ws = new WhittleSmoother(n1,n2,_tensors,_modelC);
    //ws.testSpd(n1,n2); if (ws!=null) return null;
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
    A a = new A(x1,x2,s1,s2,ws,sd);
    M m = new M(x1,x2,s1,s2,ws,sd);
    CgSolver cgs = new CgSolver(0.001,100);
    CgSolver.Info info = cgs.solve(a,m,vf,vz);
    //CgSolver.Info info = cgs.solve(a,vf,vz);
    cgs = null;
    float[][] q = ws.apply(s1,s2,x1,x2,az);
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

  private static class WhittleSmoother {
    WhittleSmoother(int n1, int n2, Tensors2 tensors, Covariance modelC) {
      Check.argument((modelC instanceof Matern),"Covariance is Matern");
      double shapeM = ((Matern)modelC).getShape();
      double sigmaM = ((Matern)modelC).getSigma();
      double rangeM = ((Matern)modelC).getRange();
      Check.argument(shapeM==1.0,"Matern shape = 1.0");
      _tensors = tensors;
      _modelC = (Matern)modelC;
      _cscale = new float[n2][n1];
      float[] d = new float[3];
      float cscale = (float)(PI*sigmaM*sigmaM*rangeM*rangeM);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          _tensors.getTensor(i1,i2,d);
          float d11 = d[0], d12 = d[1], d22 = d[2];
          float det = d11*d22-d12*d12;
          _cscale[i2][i1] = cscale*det;
        }
      }
      _dscale = (float)(rangeM*rangeM/4.0);
      _lsf = new LocalSmoothingFilter(1.0e-6,1000);
    }
    void apply(float[][] q) {
      int n1 = q[0].length;
      int n2 = q.length;
      float[][] t = new float[n2][n1];
      _lsf.applySmoothS(q,q);
      _lsf.apply(_tensors,_dscale,q,t);
      mul(_cscale,t,t);
      _lsf.apply(_tensors,_dscale,t,q);
      _lsf.applySmoothS(q,q);
    }
    float[][] apply(Sampling s1, Sampling s2, 
      float[] x1, float[] x2, float[] z) 
    {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n = z.length;
      float[][] q = new float[n2][n1];
      for (int i=0; i<n; ++i) {
        int i1 = s1.indexOfNearest(x1[i]);
        int i2 = s2.indexOfNearest(x2[i]);
        q[i2][i1] = z[i];
      }
      apply(q);
      return q;
    }
    Tensors2 getTensors() {
      return _tensors;
    }
    Covariance getModelC() {
      return _modelC;
    }
    void testSpd(int n1, int n2) {
      float[][] x = sub(randfloat(n1,n2),0.5f);
      float[][] y = sub(randfloat(n1,n2),0.5f);
      float[][] ax = copy(x);
      float[][] ay = copy(y);
      apply(ax);
      apply(ay);
      float xay = sum(mul(x,ay));
      float yax = sum(mul(y,ax));
      float xax = sum(mul(x,ax));
      float yay = sum(mul(y,ay));
      System.out.println("xax="+xax+" yay="+yay);
      System.out.println("xay="+xay+" yax="+yax);
    }
    private Tensors2 _tensors;
    private Matern _modelC;
    private float[][] _cscale;
    private float _dscale;
    private LocalSmoothingFilter _lsf;
  }

  private static class A implements CgSolver.A {

    A(float[] x1, float[] x2, Sampling s1, Sampling s2,
      WhittleSmoother ws, float[] sd) 
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
      _sd = sd;
      _ws = ws;
    }

    public void apply(Vec x, Vec y) {
      float[] ax = ((VecArrayFloat1)x).getArray();
      float[] ay = ((VecArrayFloat1)y).getArray();
      int n = ax.length;
      zero(_bx);
      for (int i=0; i<n; ++i)
        _bx[_k2[i]][_k1[i]] = ax[i];
      _ws.apply(_bx);
      for (int i=0; i<n; ++i)
        ay[i] = _bx[_k2[i]][_k1[i]]+_sd[i]*_sd[i]*ax[i];
    }

    private int[] _k1,_k2;
    private float[][] _bx;
    private float[] _sd;
    private WhittleSmoother _ws;
  }

  private static class M implements CgSolver.A {

    M(float[] x1, float[] x2, Sampling s1, Sampling s2,
      WhittleSmoother ws, float[] sd) 
    {
      Tensors2 tensors = ws.getTensors();
      Covariance modelC = ws.getModelC();
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n = x1.length;
      DMatrix cm = new DMatrix(n,n);
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
          cm.set(i,j,s*modelC.evaluate(r));
        }
        cm.set(i,i,cm.get(i,i)+sd[i]*sd[i]);
      }
      cm = cm.inverse();
      _am = new float[n][n];
      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          _am[i][j] = (float)cm.get(i,j);
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
