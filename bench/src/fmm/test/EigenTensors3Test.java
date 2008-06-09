/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm.test;

import junit.framework.*;

import fmm.EigenTensors3;

/**
 * Tests {@link fmm.EigenTensors3}.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.09
 */
public class EigenTensors3Test extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(EigenTensors3Test.class);
    junit.textui.TestRunner.run(suite);
  }

  public static void testRandom() {
    testRandom(false,0.1,1.0e-6);
    testRandom(true,1.0,1.0e-3);
  }

  private static void testRandom(
    boolean compressed, double errorAngle, double errorCoeff) 
  {
    int n1 = 3, n2 = 4, n3 = 5;
    EigenTensors3 et = new EigenTensors3(n1,n2,n3,false);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] a = makeRandomCoefficients();
          float[] u = makeRandomVector();
          float[] w = makeOrthogonalVector(u);
          et.setCoefficients(i1,i2,i3,a);
          et.setEigenvectorU(i1,i2,i3,u);
          et.setEigenvectorW(i1,i2,i3,w);
          float[] c;
          c = et.getEigenvectorU(i1,i2,i3); checkVectors(u,c,errorAngle);
          c = et.getEigenvectorW(i1,i2,i3); checkVectors(w,c,errorAngle);
          c = et.getCoefficients(i1,i2,i3); checkCoefficients(c,a,errorCoeff);
          et.setTensor(i1,i2,i3,et.getTensor(i1,i2,i3));
          c = et.getEigenvectorU(i1,i2,i3); checkVectors(u,c,errorAngle);
          c = et.getEigenvectorW(i1,i2,i3); checkVectors(w,c,errorAngle);
          c = et.getCoefficients(i1,i2,i3); checkCoefficients(c,a,errorCoeff);
        }
      }
    }
  }

  private static void checkCoefficients(float[] c, float[] a, double e) {
    assertEquals(c[0],a[0],e);
    assertEquals(c[1],a[1],e);
    assertEquals(c[2],a[2],e);
  }

  private static void checkVectors(float[] u, float[] v, double e) {
    float uv = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
    float uu = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
    float vv = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    double ca = Math.min(uv,1.0f);
    double a = Math.toDegrees(Math.acos(ca));
    assertEquals(0.0,a,e);
  }

  private static java.util.Random r = new java.util.Random();

  // Random coefficients.
  private static float[] makeRandomCoefficients() {
    float a1 = r.nextFloat();
    float a2 = r.nextFloat();
    float a3 = r.nextFloat();
    return new float[]{a1,a2,a3};
  }

  // Random unit vector with non-negative 3rd component.
  private static float[] makeRandomVector() {
    float a = r.nextFloat()-0.5f;
    float b = r.nextFloat()-0.5f;
    float c = r.nextFloat()-0.5f;
    if (c<0.0f) {
      a = -a;
      b = -b;
      c = -c;
    }
    float s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }

  // Random unit vector orthogonal to specified vector.
  private static float[] makeOrthogonalVector(float[] v1) {
    float a1 = v1[0];
    float b1 = v1[1];
    float c1 = v1[2];
    float a2 = r.nextFloat()-0.5f;
    float b2 = r.nextFloat()-0.5f;
    float c2 = r.nextFloat()-0.5f;
    float d11 = a1*a1+b1*b1+c1*c1;
    float d12 = a1*a2+b1*b2+c1*c2;
    float s = d12/d11;
    float a = a2-s*a1;
    float b = b2-s*b1;
    float c = c2-s*c1;
    if (c<0.0f) {
      a = -a;
      b = -b;
      c = -c;
    }
    s = 1.0f/(float)Math.sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }
}
