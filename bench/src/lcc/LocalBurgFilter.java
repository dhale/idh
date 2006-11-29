/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Local Burg prediction error filtering for 2-D and 3-D images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.11.28
 */
public class LocalBurgFilter {

  public LocalBurgFilter(double sigma) {
    _rgf = new RecursiveGaussianFilter(sigma);
  }

  public float[][][] applyQ1(int m, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] f = y;
    float[][] b = new float[n2][n1];
    float[][] t = new float[n2][n1];
    float[][] fb = new float[n2][n1];
    float[][] ff = new float[n2][n1];
    float[][] bb = new float[n2][n1];
    float[][][] c = new float[m][n2][n1];
    Array.copy(x,f);
    Array.copy(x,b);
    for (int k=0; k<m; ++k) {
      int k1 = (k+2)/2;
      int k2 = (k+1)/2;
      int j1 = (k+1)%2;
      int j2 = (k+0)%2;
      for (int i2=k2; i2<n2; ++i2) {
        for (int i1=k1; i1<n1; ++i1) {
          float fk = f[i2][i1];
          float bk = b[i2-j2][i1-j1];
          fb[i2][i1] = fk*bk;
          ff[i2][i1] = fk*fk;
          bb[i2][i1] = bk*bk;
        }
      }
      _rgf.apply0X(fb,t);  _rgf.applyX0(t,fb);
      _rgf.apply0X(ff,t);  _rgf.applyX0(t,ff);
      _rgf.apply0X(bb,t);  _rgf.applyX0(t,bb);
      for (int i2=n2-1; i2>=k2; --i2) {
        for (int i1=n1-1; i1>=k1; --i1) {
          float cn = 2.0f*fb[i2][i1];
          float cd = ff[i2][i1]+bb[i2][i1];
          float ck = (cd!=0.0f)?cn/cd:0.0f;
          float fk = f[i2][i1];
          float bk = b[i2-j2][i1-j1];
          f[i2][i1] = fk-ck*bk;
          b[i2][i1] = bk-ck*fk;
          c[k][i2][i1] = ck;
        }
      }
    }
    return c;
  }

  public float[][][] applyQ4(int m, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] f = y;
    float[][] b = new float[n2][n1];
    float[][] t = new float[n2][n1];
    float[][] fb = new float[n2][n1];
    float[][] ff = new float[n2][n1];
    float[][] bb = new float[n2][n1];
    float[][][] c = new float[m][n2][n1];
    Array.copy(x,f);
    Array.copy(x,b);
    for (int k=0; k<m; ++k) {
      int k2 = (k+2)/2;
      int k1 = (k+1)/2;
      int j2 = (k+1)%2;
      int j1 = (k+0)%2;
      for (int i2=0; i2<n2-k2; ++i2) {
        for (int i1=k1; i1<n1; ++i1) {
          float fk = f[i2][i1];
          float bk = b[i2+j2][i1-j1];
          fb[i2][i1] = fk*bk;
          ff[i2][i1] = fk*fk;
          bb[i2][i1] = bk*bk;
        }
      }
      _rgf.apply0X(fb,t);  _rgf.applyX0(t,fb);
      _rgf.apply0X(ff,t);  _rgf.applyX0(t,ff);
      _rgf.apply0X(bb,t);  _rgf.applyX0(t,bb);
      for (int i2=0; i2<n2-k2; ++i2) {
        for (int i1=n1-1; i1>=k1; --i1) {
          float cn = 2.0f*fb[i2][i1];
          float cd = ff[i2][i1]+bb[i2][i1];
          float ck = (cd!=0.0f)?cn/cd:0.0f;
          float fk = f[i2][i1];
          float bk = b[i2+j2][i1-j1];
          f[i2][i1] = fk-ck*bk;
          b[i2][i1] = bk-ck*fk;
          c[k][i2][i1] = ck;
        }
      }
    }
    return c;
  }

  public float[][][] applyQ3(int m, float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] f = y;
    float[][] b = new float[n2][n1];
    float[][] t = new float[n2][n1];
    float[][] fb = new float[n2][n1];
    float[][] ff = new float[n2][n1];
    float[][] bb = new float[n2][n1];
    float[][][] c = new float[m][n2][n1];
    Array.copy(x,f);
    Array.copy(x,b);
    for (int k=0; k<m; ++k) {
      int k1 = (k+2)/2;
      int k2 = (k+1)/2;
      int j1 = (k+1)%2;
      int j2 = (k+0)%2;
      for (int i2=0; i2<n2-k2; ++i2) {
        for (int i1=0; i1<n1-k1; ++i1) {
          float fk = f[i2][i1];
          float bk = b[i2+j2][i1+j1];
          fb[i2][i1] = fk*bk;
          ff[i2][i1] = fk*fk;
          bb[i2][i1] = bk*bk;
        }
      }
      _rgf.apply0X(fb,t);  _rgf.applyX0(t,fb);
      _rgf.apply0X(ff,t);  _rgf.applyX0(t,ff);
      _rgf.apply0X(bb,t);  _rgf.applyX0(t,bb);
      for (int i2=0; i2<n2-k2; ++i2) {
        for (int i1=0; i1<n1-k1; ++i1) {
          float cn = 2.0f*fb[i2][i1];
          float cd = ff[i2][i1]+bb[i2][i1];
          float ck = (cd!=0.0f)?cn/cd:0.0f;
          float fk = f[i2][i1];
          float bk = b[i2+j2][i1+j1];
          f[i2][i1] = fk-ck*bk;
          b[i2][i1] = bk-ck*fk;
          c[k][i2][i1] = ck;
        }
      }
    }
    return c;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private RecursiveGaussianFilter _rgf;
}
