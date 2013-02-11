/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import static java.lang.Math.*;

/**
 * Simplified dynamic warping to find shifts between two sequences.
 * This class implements the pseudocode in Hale and Compton (2013).
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.02.04
 */
public class DynamicWarpingS {

  public static float[][] computeErrors(int nl, float[] f, float[] g) {
    int ni = f.length;
    float[][] e = new float[ni][nl];
    for (int i=0; i<ni; ++i) {
      for (int l=0; l<nl; ++l) {
        float fmg = f[i]-g[i+l];
        e[i][l] = fmg*fmg;
      }
    }
    return e;
  }

  public static float[] findShifts(double rl, double ru, float[][] e) {
    int ni = e.length;
    int nl = e[0].length;
    float[][] d = new float[ni][nl];
    int[][] m = new int[ni][nl];
    float[] u = new float[ni];
    for (int l=0; l<nl; ++l) // initialize
      d[0][l] = e[0][l];
    for (int i=1; i<ni; ++i) { // accumulate
      for (int l=0; l<nl; ++l) {
        float dl = Float.MAX_VALUE;
        int ml = -1;
        int ql = max((int)ceil(rl),l-nl+1);
        int qu = min((int)floor(ru),l);
        for (int q=ql; q<=qu; ++q) {
          float dq = d[i-1][l-q]+e[i][l];
          if (dq<dl) {
            dl = dq;
            ml = q;
          }
        }
        d[i][l] = dl;
        m[i][l] = ml;
      }
    }
    int i = ni-1; // minimize d
    float di = Float.MAX_VALUE;
    for (int l=0; l<nl; ++l) {
      if (d[i][l]<di) {
        di = d[i][l];
        u[i] = l;
      }
    }
    while (i>0) { // backtrack
      u[i-1] = u[i]-m[i][(int)u[i]];
      i = i-1;
    }
    return u;
  }

  public static float[] findShiftsI(
    double rl, double ru, float[][] e, int[] i) 
  {
    int ni = e.length;
    int nl = e[0].length;
    int nj = i.length;
    float[][] d = new float[nj][nl];
    int[][] m = new int[nj][nl];
    float[] ui = new float[nj];
    for (int l=0; l<nl; ++l) // initialize
      d[0][l] = e[0][l];
    for (int j=1; j<nj; ++j) { // accumulate
      int h = i[j]-i[j-1];
      for (int l=0; l<nl; ++l) {
        float dl = Float.MAX_VALUE;
        int ml = -1;
        int ql = max((int)ceil(h*rl),l-nl+1);
        int qu = min((int)floor(h*ru),l);
        for (int q=ql; q<=qu; ++q) {
          float dq = d[j-1][l-q];
          for (int p=0; p<h; ++p) {
            int s = (int)(l-(float)p*q/h+0.5);
            dq += e[i[j]-p][s];
          }
          if (dq<dl) {
            dl = dq;
            ml = q;
          }
        }
        d[j][l] = dl;
        m[j][l] = ml;
      }
    }
    int j = nj-1; // minimize d
    float dj = Float.MAX_VALUE;
    for (int l=0; l<nl; ++l) {
      if (d[j][l]<dj) {
        dj = d[j][l];
        ui[j] = l;
      }
    }
    while (j>0) { // backtrack
      ui[j-1] = ui[j]-m[j][(int)ui[j]];
      j = j-1;
    }
    return ui;
  }
}
