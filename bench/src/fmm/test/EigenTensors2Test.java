/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm.test;

import junit.framework.*;

import fmm.EigenTensors2;

/**
 * Tests {@link fmm.EigenTensors2}.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.09
 */
public class EigenTensors2Test extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(EigenTensors2Test.class);
    junit.textui.TestRunner.run(suite);
  }

  public static void testRandom() {
    testRandom(0.1,1.0e-6);
  }

  private static void testRandom(double errorAngle, double errorCoeff) {
    int n1 = 13, n2 = 14;
    EigenTensors2 et = new EigenTensors2(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float[] a = makeRandomCoefficients();
        float[] u = makeRandomVector();
        et.setCoefficients(i1,i2,a);
        et.setEigenvectorU(i1,i2,u);
        float[] c;
        c = et.getEigenvectorU(i1,i2); checkVectors(u,c,errorAngle);
        c = et.getCoefficients(i1,i2); checkCoefficients(c,a,errorCoeff);
        et.setTensor(i1,i2,et.getTensor(i1,i2));
        c = et.getEigenvectorU(i1,i2); checkVectors(u,c,errorAngle);
        c = et.getCoefficients(i1,i2); checkCoefficients(c,a,errorCoeff);
      }
    }
  }

  private static void checkCoefficients(float[] c, float[] a, double e) {
    assertEquals(c[0],a[0],e);
    assertEquals(c[1],a[1],e);
  }

  private static void checkVectors(float[] u, float[] v, double e) {
    float uv = u[0]*v[0]+u[1]*v[1];
    double ca = Math.max(-1.0,Math.min(uv,1.0f));
    double a = Math.toDegrees(Math.acos(ca));
    if (a>90.0f)
      a -= 180.0f;
    assertEquals(0.0,a,e);
  }

  private static java.util.Random r = new java.util.Random();

  // Random coefficients.
  private static float[] makeRandomCoefficients() {
    float a1 = r.nextFloat();
    float a2 = r.nextFloat();
    return new float[]{a1,a2};
  }

  // Random unit vector with non-negative 3rd component.
  private static float[] makeRandomVector() {
    float a = r.nextFloat()-0.5f;
    float b = r.nextFloat()-0.5f;
    float s = 1.0f/(float)Math.sqrt(a*a+b*b);
    return new float[]{a*s,b*s};
  }
}
