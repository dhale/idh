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
 * Uses seismic image samples alongside fault skins to estimate fault dip
 * slips.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */
public class FaultSlipper {

    /**
     * Computes alignment errors and initializes shifts. Computes both
     * minus-plus (emp) and plus-minus (epm) errors. The minus-plus errors
     * correspond to differences between the sample value on the minus side of
     * this cell and those for the plus sides of cells up and down dip from
     * this cell. The plus-minus errors are defined similarly.
     * <p>
     * Uses the specified slopes to initialize both minus-plus and plus-minus
     * shifts to compensate for the fact that shifts are estimated using image
     * samples located a horizontal distance d away from this cell.
     * <p>
     * For lags where image sample values are unavailable, say, near surface
     * boundaries, errors are extrapolated from other lags, but are negated,
     * so that extrapolated errors can be detected and modified later after
     * errors for all cells in the surface have been computed.
     */
    /*
    void computeErrorsAndInitShifts(
      int lmax, float d, float[][][] f, float[][][] p2, float[][][] p3) 
    {
      int n1 = f[0][0].length;
      int nlag;
      float y1,y2,y3,gm,gp;
      Cell ci;

      // Errors for lag zero.
      emp = new float[lmax+1+lmax];
      epm = new float[lmax+1+lmax];
      float us = q.us23, u1 = us*q.u1, u2 = us*q.u2, u3 = us*q.u3;
      float ut = (i1-q.c1)*u1;
      float x1 = i1, x2 = q.c2-ut*u2, x3 = q.c3-ut*u3;

      float x2m = x2-d*w2, x3m = x3-d*w3;
      float x2p = x2+d*w2, x3p = x3+d*w3;
      float fm = valueAt(x1,x2m,x3m,f);
      float fp = valueAt(x1,x2p,x3p,f);
      float p2m = valueAt(x1,x2m,x3m,p2);
      float p2p = valueAt(x1,x2p,x3p,p2);
      float p3m = valueAt(x1,x2m,x3m,p3);
      float p3p = valueAt(x1,x2p,x3p,p3);
      float empl = emp[lmax] = error(fm,fp);
      float epml = epm[lmax] = error(fp,fm);

      // Initial shifts compensate for horizontal distance d.
      float s23 = d*((p2m+p2p)*w2+(p3m+p3p)*w3);
      smp = -s23;
      spm =  s23;

      // Errors for samples south; any extrapolated errors are negative.
      us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
      y1 = x1+1.0f; y2 = x2-u1*u2; y3 = x3-u1*u3;
      qi = quadSouth(q,y1,y2,y3);
      nlag = min(lmax,n1-1-i1);
      for (int ilag=1; ilag<=lmax; ++ilag) {
        if (qi!=null && ilag<=nlag) {
          us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
          ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
          y2 -= ut*u2; y3 -= ut*u3;
          gm = valueAt(y1,y2-d*u2,y3-d*u3,f);
          gp = valueAt(y1,y2+d*u2,y3+d*u3,f);
          empl = emp[lmax+ilag] = error(fm,gp);
          epml = epm[lmax+ilag] = error(fp,gm);
          y1 += 1.0f; y2 -= u1*u2; y3 -= u1*u3;
          qi = quadSouth(qi,y1,y2,y3);
        } else {
          emp[lmax+ilag] = -empl;
          epm[lmax+ilag] = -epml;
        }
      }

      // Errors for samples north; any extrapolated errors are negative.
      us = q.us23; u1 = us*q.u1; u2 = us*q.u2; u3 = us*q.u3;
      y1 = x1-1.0f; y2 = x2+u1*u2; y3 = x3+u1*u3;
      qi = quadNorth(q,y1,y2,y3);
      nlag = min(lmax,i1);
      for (int ilag=1; ilag<=lmax; ++ilag) {
        if (qi!=null && ilag<=nlag) {
          us = qi.us23; u1 = us*qi.u1; u2 = us*qi.u2; u3 = us*qi.u3;
          ut = (y1-qi.c1)*u1+(y2-qi.c2)*u2+(y3-qi.c3)*u3;
          y2 -= ut*u2; y3 -= ut*u3;
          gm = valueAt(y1,y2-d*u2,y3-d*u3,f);
          gp = valueAt(y1,y2+d*u2,y3+d*u3,f);
          empl = emp[lmax-ilag] = error(fm,gp);
          epml = epm[lmax-ilag] = error(fp,gm);
          y1 -= 1.0f; y2 += u1*u2; y3 += u1*u3;
          qi = quadNorth(qi,y1,y2,y3);
        } else {
          emp[lmax-ilag] = -empl;
          epm[lmax-ilag] = -epml;
        }
      }
    }
  */
}
