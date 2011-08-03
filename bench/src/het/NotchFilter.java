/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package het;

import edu.mines.jtk.dsp.RecursiveParallelFilter;
import edu.mines.jtk.util.Cdouble;
import static edu.mines.jtk.util.ArrayMath.*;

// for testing only
import edu.mines.jtk.mosaic.*;

/**
 * A symmetric notch filter.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.03.25
 */
public class NotchFilter extends RecursiveParallelFilter {

  public NotchFilter(double fn, double rn) {
    double wn = 2.0*PI*fn;
    Cdouble pn = Cdouble.polar(rn,wn);
    Cdouble zn = Cdouble.polar(1.0,wn);
    Cdouble[] poles = {pn,pn.conj()};
    Cdouble[] zeros = {zn,zn.conj()};
    init(poles,zeros,1.0);
  }
}
