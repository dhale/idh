/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

/**
 * Computes a sampled function x(y) for a sampled increasing function y(x).
 * The function y(x) is specified by an input sampling of x and an input 
 * array of strictly (monotonically) increasing y values. The computed 
 * function x(y) is represented by a specified output sampling of y and 
 * computed values x in an output array that also increase monotonically.
 * 
 * Currently, this class uses only inverse linear interpolation (and
 * extrapolation) to compute the sampled values of x(y) from the specified 
 * sampled values of y(x).
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.08
 */
public class InverseInterpolator {

  public InverseInterpolator(int ni, int no) {
    this(new Sampling(ni),new Sampling(no));
  }

  public InverseInterpolator(Sampling si, Sampling so) {
    Check.argument(si.getCount()>1,"at least two input samples");
    _si = si;
    _so = so;
  }

  public void invert(float[] y, float[] x) {
    int nxi = _si.getCount();
    int nyo = _so.getCount();
    Check.argument(y.length==nxi,"y.length equals number of input samples");
    Check.argument(x.length==nyo,"x.length equals number of output samples");
    int nxim1 = nxi-1;
    int jxi1 = 0;
    int jxi2 = 1;
    float xi1 = (float)_si.getValue(jxi1);
    float xi2 = (float)_si.getValue(jxi2);
    float yi1 = y[jxi1];
    float yi2 = y[jxi2];
    Check.argument(yi1<yi2,"y values strictly increasing");
    float dxody = (xi2-xi1)/(yi2-yi1);
    int jyo = 0;
    float yo = (float)_so.getValue(jyo);
    while (jyo<nyo) {
      if (yo<=yi2 || jxi2==nxim1) {
        x[jyo++] = xi1+(yo-yi1)*dxody;
        if (jyo<nyo) 
          yo = (float)_so.getValue(jyo);
      } else if (jxi2<nxim1) {
        xi1 = (float)_si.getValue(++jxi1);
        xi2 = (float)_si.getValue(++jxi2);
        yi1 = y[jxi1];
        yi2 = y[jxi2];
        Check.argument(yi1<yi2,"y values strictly increasing");
        dxody = (xi2-xi1)/(yi2-yi1);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _si,_so;
}
