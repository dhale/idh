/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import edu.mines.jtk.util.*;

/**
 * Local anisotropic diffusion filter via conjugate gradient iterations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.07.05
 */
public class LocalDiffusionFilterCg extends LocalDiffusionFilter {

  /**
   * Constructs a local diffusion filter.
   * @param sigma the nominal half-width for this filter.
   * @param small stop when sum of residuals squared decreases by this factor.
   * @param niter stop when number of iterations exceeds this number.
   */
  public LocalDiffusionFilterCg(double sigma, double small, int niter) {
    super(sigma);
    _dlf = new DirectionalLaplacianFilter(sigma);
    _sigma = (float)sigma;
    _small = (float)small;
    _niter = niter;
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void solveInline(
    float[][] ds, float[][] v1, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    final float[][] dsf = ds;
    final float[][] v1f = v1;
    Operator op = new Operator() {
      public void apply(float[][] x, float[][] y) {
        Array.copy(x,y);
        _dlf.applyInline(dsf,v1f,x,y);
      }
    };
    solveCg(op,x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma; // maximum filter half-width
  private float _small; // stop iterations when rr decreases by this factor
  private int _niter; // number of iterations
  private DirectionalLaplacianFilter _dlf;

  /**
   * Symmetric positive-definite operator to be inverted.
   */
  private static interface Operator {
    public void apply(float[][] x, float[][] y);
  }

  /**
   * Solves Ay = x via conjugate gradient iterations.
   */
  private void solveCg(Operator a, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    float[][] t = new float[n2][n1];
    a.apply(y,t);
    double rr = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      float[] x2 = x[i2];
      float[] r2 = r[i2];
      float[] s2 = s[i2];
      float[] t2 = t[i2];
      for (int i1=0; i1<n1; ++i1) {
        float ri = x2[i1]-t2[i1];
        r2[i1] = ri;
        s2[i1] = ri;
        rr += ri*ri;
      }
    }
    trace("solveCg: r="+rr);
    int miter;
    double rrsmall = rr*_small;
    for (miter=0; miter<_niter && rr>rrsmall; ++miter) {
      a.apply(s,t);
      double st = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] s2 = s[i2];
        float[] t2 = t[i2];
        for (int i1=0; i1<n1; ++i1)
          st += s2[i1]*t2[i1];
      }
      float alpha = (float)(rr/st);
      double rrold = rr;
      rr = 0.0;
      for (int i2=0; i2<n2; ++i2) {
        float[] y2 = y[i2];
        float[] r2 = r[i2];
        float[] s2 = s[i2];
        float[] t2 = t[i2];
        for (int i1=0; i1<n1; ++i1) {
          y2[i1] += alpha*s2[i1];
          r2[i1] -= alpha*t2[i1];
          rr += r2[i1]*r2[i1];
        }
      }
      if (rr<=rrsmall)
        break;
      float beta = (float)(rr/rrold);
      for (int i2=0; i2<n2; ++i2) {
        float[] r2 = r[i2];
        float[] s2 = s[i2];
        for (int i1=0; i1<n1; ++i1)
          s2[i1] = r2[i1]+beta*s2[i1];
      }
    }
    trace("  miter="+miter+" rr="+rr);
  }

  private static final boolean TRACE = true;
  private static void trace(String s) {
    if (TRACE)
      System.out.println(s);
  }
} 
