/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package interp;

import static java.lang.Math.*;

/**
 * A covariance function of distance.
 * @author Dave Hale, Colorado School of Mines.
 */
public interface Covariance {

  /**
   * Returns the covariance for the specified distance.
   * @param r the distance.
   * @return the covariance.
   */
  public double evaluate(double r);
}
