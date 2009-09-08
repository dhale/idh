/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp.test;

import java.util.Random;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import dnp.*;
import edu.mines.jtk.dsp.Sampling;

/**
 * Tests {@link dnp.InverseInterpolator}.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.08
 */
public class InverseInterpolatorTest extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(InverseInterpolatorTest.class);
    junit.textui.TestRunner.run(suite);
  }

  public void testLinear() {
    int n = 1000, nx = n, ny = n;
    Sampling sx = makeRandomSampling(nx);
    Sampling sy = new Sampling(ny);
    float[] y = new float[nx];
    for (int ix=0; ix<nx; ++ix)
      y[ix] = (float)sx.getValue(ix); // simple line y(x) = x
    float[] x = new float[ny];
    InverseInterpolator ii = new InverseInterpolator(sx,sy);
    ii.invert(y,x);
    float tiny = 0.001f;
    for (int iy=0; iy<ny; ++iy) {
      if (Math.abs(x[iy]-iy)>tiny) {
        System.out.println("iy="+iy+" x[iy]="+x[iy]);
        System.out.println(" last x="+sx.getValue(nx-2)+","+sx.getLast());
        System.out.println(" last y="+y[nx-2]+","+y[nx-1]);
      }
      assertEquals(iy,x[iy],tiny); // should have x(y) = y
    }
  }

  private static Sampling makeRandomSampling(int n) {
    double[] x = new double[n];
    Random r = new Random();
    x[0] = 3.14159;
    for (int i=1; i<n; ++i)
      x[i] = x[i-1]+0.001+1.999*r.nextFloat();
    return new Sampling(x);
  }
}
