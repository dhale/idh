/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

// FOR TESTING ONLY
import javax.swing.*;
import edu.mines.jtk.mosaic.*;

/**
 * Wraps RecursiveGaussianFilter for linearly extrapolated end conditions.
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.03.23
 */
public class RgSmoother {

  public RgSmoother(double sigma) {
    _sigma = sigma;
    _rgf = new RecursiveGaussianFilter(sigma);
  }

  public void apply0(float[] x, float[] y) {
    int n = x.length;
    float[] z = extrap(x);
    _rgf.apply0(z,z);
    dxtrap(z,y);
  }

  public void apply1(float[] x, float[] y) {
    int n = x.length;
    float[] z = extrap(x);
    _rgf.apply1(z,z);
    dxtrap(z,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _sigma;
  private RecursiveGaussianFilter _rgf;

  private int nextrap() {
    return 10+(int)(6*_sigma);
  }
  
  private float slopeB(float[] x) {
    int n = x.length;
    int m = min(n-1,nextrap());
    float a = 1.0f;
    float b = 1.0f/m;
    float s = 0.0f;
    float t = 0.0f;
    for (int i=1; i<=m; ++i) {
      s += a*(x[i]-x[i-1]); // weighted sum of differences
      t += a; // sum of weights
      a -= b; // weight decreases linearly
    }
    return s/t;
  }

  private float slopeE(float[] x) {
    int n = x.length;
    int m = min(n-1,nextrap());
    float a = 1.0f;
    float b = 1.0f/m;
    float s = 0.0f;
    float t = 0.0f;
    for (int i=n-1; i>=n-m; --i) {
      s += a*(x[i]-x[i-1]); // weighted sum of differences
      t += a; // sum of weights
      a -= b; // weight decreases linearly
    }
    return s/t;
  }

  // Returns x padded with extrapolated values.
  private float[] extrap(float[] x) {
    int m = nextrap();
    int n = x.length;
    int np = m+n+m;
    float sb = slopeB(x);
    float se = slopeE(x);
    float xb = x[0];
    float xe = x[n-1];
    float[] xp = new float[np];
    for (int i=0; i<m; ++i) {
      xp[m-1-i] = xb-(i+1)*sb;
      xp[m+n+i] = xe+(i+1)*se;
    }
    copy(n,0,x,m,xp);
    return xp;
  }

  // Copies, discarding the extrapolated part.
  private void dxtrap(float[] xp, float[] x) {
    int n = x.length;
    int m = nextrap();
    copy(n,m,xp,0,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Testing

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testSmoother();
      }
    });
  }
  public static void testSmoother() {
    int n = 501;
    float sigma = 50.0f;
    float[] x = rampfloat(100.0f,0.1f,n);
    float[] y = new float[n];
    float[] z = new float[n];
    RgSmoother rs = new RgSmoother(sigma);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rs.apply1(x,y);
    rgf.apply1(x,z);
    //dump(y);
    SimplePlot.asPoints(y);
    SimplePlot.asPoints(z);
  }
}
