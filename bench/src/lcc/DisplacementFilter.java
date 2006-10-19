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
 * Filters displacements to minimize strain.
 * Minimizes J(u) = u* Q u + s u* D* D u.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.10.19
 */
public class DisplacementFilter {

  /**
   * Construct a displacement filter with specified parameters.
   * @param smooth the amount of smoothing.
   */
  public DisplacementFilter(double smooth) {
  }

  /**
   * Applies this filter for the specified quadratic fit.
   * @param q coefficients for the quadratic fit.
   *  q[0] contains a measure of quality between 0 and 1; and
   *  q[1], q[2], and q[3] contain coefficients of the matrix Q:
   *  <pre>
   *    Q = | q[1]  q[2] |
   *        | q[2]  q[3] |
   *  </pre>
   * @param u the displacements to be filtered, in-place.
   *  The components (u1,u2) are (u[0],u[1]).
   */
  public void apply(float[][][] q, float[][][] u) {
    
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
}
