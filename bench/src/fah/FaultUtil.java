/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Utilities for working with faults.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.03
 */
public class FaultUtil {

  /**
   * Returns fault normal vector for specified fault strike and dip.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return u output array {u1,u2,u3} with components of normal vector.
   */
  public static float[] faultNormalFromStrikeAndDip(float phi, float theta) {
    float[] u = new float[3];
    faultNormalFromStrikeAndDip(phi,theta,u);
    return u;
  }

  /**
   * Computes fault normal vector for specified fault strike and dip.
   * The normal vector will point upwards.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @param u output array {u1,u2,u3} with components of normal vector.
   */
  public static void faultNormalFromStrikeAndDip(
      float phi, float theta, float[] u) {
    float u1 = 0.0f;
    float u2 = 0.0f;
    float u3 = 0.0f;
    if (phi==90.0f && theta==0.0f) {
      u2 = -1.0f;
    } else if (phi==-90.0f && theta==0.0f) {
      u2 = 1.0f;
    } else {
      float p = toRadians(phi);
      float t = toRadians(theta);
      float cp = cos(p);
      float sp = sin(p);
      float ct = cos(t);
      float st = sin(t);
      u1 = -st;
      u2 = -sp*ct;
      u3 =  cp*ct;
      if (u1>0.0f) {
        u1 = -u1;
        u2 = -u2;
        u3 = -u3;
      }
    }
    u[0] = u1;
    u[1] = u2;
    u[2] = u3;
  }

  /**
   * Returns fault strike angle for specified fault normal vector.
   * @param u1 1st component of fault normal vector.
   * @param u2 2nd component of fault normal vector.
   * @param u3 3rd component of fault normal vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromNormal(float u1, float u2, float u3) {
    return toDegrees(atan(-u2/u3));
  }

  /**
   * Returns fault strike angle for specified fault normal vector.
   * @param u array {u1,u2,u3} of components of fault normal vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromNormal(float[] u) {
    return faultStrikeFromNormal(u[0],u[1],u[2]);
  }

  /**
   * Returns fault dip angle for specified fault normal vector.
   * @param u1 1st component of fault normal vector.
   * @param u2 2nd component of fault normal vector.
   * @param u3 3rd component of fault normal vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromNormal(float u1, float u2, float u3) {
    return u3>=0.0f ? toDegrees(asin(-u1)) : toDegrees(asin(u1));
  }

  /**
   * Returns fault dip angle for specified fault normal vector.
   * @param u array {u1,u2,u3} of components of fault normal vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromNormal(float[] u) {
    return faultDipFromNormal(u[0],u[1],u[2]);
  }

  /**
   * Returns fault strike vector v for the specified fault normal vector u.
   * The fault strike vector lies in a horizontal plane and is parallel
   * to fault strike. The cross product of the fault normal and strike
   * vectors equals the fault dip vector, which points down the fault.
   * @param u input array {u1,u2,u3} of components of fault normal vector.
   * @return array {v1,v2,v3} of components of fault strike vector.
   */
  public static float[] faultStrikeVectorFromNormal(float[] u) {
    float u1 = u[0], u2 = u[1], u3 = u[2];
    float v1 = 0.0f;
    float v2,v3;
    float u2s = u2*u2;
    float u3s = u3*u3;
    if (u2s>u3s) {
      v3 = 1.0f/sqrt(1.0f+u3s/u2s);
      v2 = -v3*u3/u2;
    } else {
      v2 = 1.0f/sqrt(1.0f+u2s/u3s);
      v3 = -v2*u2/u3;
    }
    if (u3*v2<u2*v3) {
      v2 = -v2;
      v3 = -v3;
    }
    return new float[]{v1,v2,v3};
  }

  /**
   * Returns fault dip vector w for the specified fault normal vector u.
   * The fault dip vector w points down the fault; it equals the cross
   * product of the fault normal vector u and fault strike vector v.
   */
  public static float[] faultDipVectorFromNormal(float[] u) {
    float[] v = faultStrikeVectorFromNormal(u);
    float u1 = u[0], u2 = u[1], u3 = u[2];
    float v1 = v[0], v2 = v[1], v3 = v[2];
    float w1 = u2*v3-u3*v2;
    float w2 = u3*v1-u1*v3;
    float w3 = u1*v2-u2*v1;
    return new float[]{w1,w2,w3};
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    float[] fp = {  0.0f, 90.0f,-90.0f,
                    0.0f, 90.0f,-90.0f,
                    0.0f, 90.0f,-90.0f};
    float[] ft = {  0.0f,  0.0f,  0.0f,
                    1.0f,  1.0f,  1.0f,
                   -1.0f, -1.0f, -1.0f};
    for (int i=0; i<fp.length; ++i)
      test(fp[i],ft[i]);
  }
  public static void test(float phik, float thetak) {
    float[] uk = faultNormalFromStrikeAndDip(phik,thetak);
    float phit = faultStrikeFromNormal(uk);
    float thetat = faultDipFromNormal(uk);
    float[] ut = faultNormalFromStrikeAndDip(phit,thetat);
    //trace("phik = "+phik+" thetak = "+thetak);
    //trace("phit = "+phit+" thetat = "+thetat);
    //trace("uk = ("+uk[0]+","+uk[1]+","+uk[2]+")");
    //trace("ut = ("+ut[0]+","+ut[1]+","+ut[2]+")");
    assertEqual(uk,ut);
  }
  public static void assertEqual(float x, float y) {
    assert abs(x-y)<0.01f;
  }
  public static void assertEqual(float[] x, float[] y) {
    for (int i=0; i<x.length; ++i)
      assertEqual(x[i],y[i]);
  }
  public static void trace(String s) {
    System.out.println(s);
  }
}
