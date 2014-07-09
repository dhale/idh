/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;
import static fah.FaultGeometry.*;

/**
 * Uses image samples alongside fault skins to estimate fault dip slips.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */
public class FaultSlipper {

  public FaultSlipper(FaultSkin[] skins) {
    _skinList = new ArrayList<FaultSkin>(skins.length);
    for (FaultSkin skin:skins)
      _skinList.add(skin);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private ArrayList<FaultSkin> _skinList;

  private void computeErrorsAndInitShifts(
      int lmax, float d, float[][][] f, float[][][] p2, float[][][] p3) {
    for (FaultSkin skin:_skinList) {
    }
  }

  /**
   * Computes alignment errors and initializes shifts for this cell. Computes
   * both minus-plus (emp) and plus-minus (epm) errors. The minus-plus errors
   * correspond to differences between the sample value on the minus side of
   * this cell and those for the plus sides of cells up and down dip from this
   * cell. The plus-minus errors are similarly defined. 
   * <p> 
   * This method uses specified slopes to initialize both minus-plus and
   * plus-minus shifts to compensate for the fact that shifts are estimated
   * using image samples located a horizontal distance d away from this cell. 
   * <p> 
   * For lags where image sample values are unavailable, say, near surface
   * boundaries, errors are extrapolated from other lags, but are negated, so
   * that extrapolated errors can be detected and modified later, after errors
   * for all relevant cells have been computed.
   */
  void computeErrorsAndInitShifts(
      FaultCell cell, int lmax, float d, 
      float[][][] f, float[][][] p2, float[][][] p3) {
    int n1 = f[0][0].length;
    float[] y = new float[3];

    // New arrays for alignment errors.
    cell.emp = new float[lmax+1+lmax];
    cell.epm = new float[lmax+1+lmax];

    // Errors for lag zero.
    float d2 =  d*cell.v3;
    float d3 = -d*cell.v2;
    float y1 = cell.x1, y2 = cell.x2, y3 = cell.x3;
    float fm = imageValueAt(y1,y2-d2,y3-d3,f);
    float fp = imageValueAt(y1,y2+d2,y3+d3,f);
    float gm = fm;
    float gp = fp;
    float p2m = imageValueAt(y1,y2-d2,y3-d3,p2);
    float p2p = imageValueAt(y1,y2+d2,y3+d3,p2);
    float p3m = imageValueAt(y1,y2-d2,y3-d3,p3);
    float p3p = imageValueAt(y1,y2+d2,y3+d3,p3);
    float empl = cell.emp[lmax] = alignmentError(fm,gp);
    float epml = cell.epm[lmax] = alignmentError(fp,gm);

    // Initial shifts compensate for horizontal distance d.
    float s23 = d*((p2m+p2p)*cell.w2+(p3m+p3p)*cell.w3);
    cell.smp = -s23;
    cell.spm =  s23;

    // Errors for samples above; make any extrapolated errors negative.
    FaultCell ca = cell;
    int nlaga = min(lmax,ca.i1);
    y1 = cell.x1; y2 = cell.x2; y3 = cell.x3;
    for (int ilag=1; ilag<=lmax; ++ilag) {
      if (ilag<=nlaga) {
        y[0] = y1; y[1] = y2; y[2] = y3;
        ca = ca.walkUpDipFrom(y);
        y1 = y[0]; y2 = y[1]; y3 = y[2];
        d2 =  d*ca.v3;
        d3 = -d*ca.v2;
        gm = imageValueAt(y1,y2-d2,y3-d3,f);
        gp = imageValueAt(y1,y2+d2,y3+d3,f);
        empl = ca.emp[lmax-ilag] = alignmentError(fm,gp);
        epml = ca.epm[lmax-ilag] = alignmentError(fp,gm);
      } else {
        ca.emp[lmax-ilag] = -empl;
        ca.epm[lmax-ilag] = -epml;
      }
    }

    // Errors for samples below; make any extrapolated errors negative.
    FaultCell cb = cell;
    int nlagb = min(lmax,n1-1-cb.i1);
    y1 = cell.x1; y2 = cell.x2; y3 = cell.x3;
    for (int ilag=1; ilag<=lmax; ++ilag) {
      if (ilag<=nlagb) {
        y[0] = y1; y[1] = y2; y[2] = y3;
        cb = cb.walkDownDipFrom(y);
        y1 = y[0]; y2 = y[1]; y3 = y[2];
        d2 =  d*cb.v3;
        d3 = -d*cb.v2;
        gm = imageValueAt(y1,y2-d2,y3-d3,f);
        gp = imageValueAt(y1,y2+d2,y3+d3,f);
        empl = cb.emp[lmax+ilag] = alignmentError(fm,gp);
        epml = cb.epm[lmax+ilag] = alignmentError(fp,gm);
      } else {
        cb.emp[lmax+ilag] = -empl;
        cb.epm[lmax+ilag] = -epml;
      }
    }
  }

  /**
   * Extrapolates alignment errors emp and epm where not computed. Errors that
   * could not be computed are negative, and are copies of errors for smaller
   * lags that could be computed. (Errors for zero lag can always be
   * computed.) 
   * <p>
   * For each lag with a negative error, this method first attempts to
   * extrapolate using other errors for the same lag stored in cell nabors
   * above or below. This first extrapolation works best when shifts vary
   * slowly with depth.
   * <p> 
   * If this first extrapolation is impossible, because the number of above
   * and below nabors for some lag is too small, then errors are extrapolated
   * using the errors already computed for other lags. Those errors are
   * already stored in the cells, but are negative, so in this second
   * extrapolation we simply change their sign.
   */
  /*
  private static void extrapolateErrors(FaultCell[][] cab) {
    int mns = qns.length;

    // For all arrays of quads linked north-south, ...
    for (int ins=0; ins<mns; ++ins) {
      int lns = qns[ins].length;
      float[][] emp = new float[lns][];
      float[][] epm = new float[lns][];

      // Get arrays of errors for all quads in the array.
      for (int jns=0; jns<lns; ++jns) {
        Vert v = qns[ins][jns].getVert();
        emp[jns] = v.emp;
        epm[jns] = v.epm;
      }
      int lmax = (emp[0].length-1)/2;

      // For each array of errors, ...
      for (int jns=0; jns<lns; ++jns) {

        // For all lags, negative and positive, ...
        for (int ilag=1; ilag<=lmax; ++ilag) {
          int ilagm = lmax-ilag;
          int ilagp = lmax+ilag;

          // Extrapolate for negative lags.
          float empim = emp[jns][ilagm];
          float epmim = epm[jns][ilagm];
          for (int kns=jns; kns<lns && empim<0.0f; ++kns)
            if (emp[kns][ilagm]>=0.0f)
              empim = emp[kns][ilagm];
          for (int kns=jns; kns<lns && epmim<0.0f; ++kns)
            if (epm[kns][ilagm]>=0.0f)
              epmim = epm[kns][ilagm];
          if (empim<0.0f) empim = -empim;
          if (epmim<0.0f) epmim = -epmim;
          emp[jns][ilagm] = empim;
          epm[jns][ilagm] = epmim;

          // Extrapolate for positive lags.
          float empip = emp[jns][ilagp];
          float epmip = epm[jns][ilagp];
          for (int kns=jns; kns>=0 && empip<0.0f; --kns)
            if (emp[kns][ilagp]>=0.0f)
              empip = emp[kns][ilagp];
          for (int kns=jns; kns>=0 && epmip<0.0f; --kns)
            if (epm[kns][ilagp]>=0.0f)
              epmip = epm[kns][ilagp];
          if (empip<0.0f) empip = -empip;
          if (epmip<0.0f) epmip = -epmip;
          emp[jns][ilagp] = empip;
          epm[jns][ilagp] = epmip;
        }
      }
    }
  }
  */

  private static float alignmentError(float f, float g) {
    float fmg = f-g;
    return fmg*fmg;
  }

  private static float imageValueAt(
    float p1, float p2, float p3, float[][][]f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    int i1 = max(0,min(n1-1,round(p1)));
    int i2 = max(0,min(n2-1,round(p2)));
    int i3 = max(0,min(n3-1,round(p3)));
    return f[i3][i2][i1];
  }
}
