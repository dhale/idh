/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gph;

import static java.lang.Math.*;

import edu.mines.jtk.util.Check;

/**
 * An oblate spheroidal model.
 * Enables conversion from (x,y,z) Cartesian coordinates 
 * to (latitude,longitude,altitude) coordinates and vice-versa. 
 * <p>
 * Units of latitude are degrees in [-90,90].
 * Units of longitude are degrees in [-180,180].
 * Units of (x,y,z) and altitude are arbitrary.
 *
 * @author Dave Hale, Colorado School of Mines.
 * @version 2010.12.29
 */
public class OblateModel {

  // Adapted from lla2ecef.m and ecef2lla.m, by Michael Kleder, 2005.

  /**
   * Returns a WGS-84 model of the earth.
   * Units of (x,y,z) and altitude are kilometers.
   */
  public static OblateModel forWGS84() {
    return new OblateModel(6378.1370000,
                           6356.7523142);
  }

  /**
   * Constructs an spherical model with specified radius.
   * @param r the radius.
   */
  public OblateModel(double r) {
    this(r,r);
  }

  /**
   * Constructs a model with specified semi-major and semi-minor axes.
   * @param a the semi-major axis (equitorial radius).
   * @param b the semi-minor axis (polar distance).
   */
  public OblateModel(double a, double b) {
    Check.argument(b<=a,"b<=a");
    _a = a;
    _b = b;
    _one_minus_ee = (b*b)/(a*a);
    _ee = 1.0-_one_minus_ee;
    _tinyxy = a*1.0e-4;
  }

  /**
   * Returns (x,y,z) coordinates for specified (lat,lon,alt).
   * @param lat latitude.
   * @param lon longitude.
   * @param alt altitude.
   * @return array of {x,y,z}.
   */
  public double[] xyz(double lat, double lon, double alt) {
    lat *= D2R; lon *= D2R;
    double slat = sin(lat);
    double clat = cos(lat);
    double slon = sin(lon);
    double clon = cos(lon);
    double a = _a;
    double b = a/sqrt(1.0-_ee*slat*slat);
    double balt = b+alt;
    double x = balt*clat*clon;
    double y = balt*clat*slon;
    double z = (_one_minus_ee*b+alt)*slat;
    return new double[]{x,y,z};
  }

  /**
   * Returns (lat,lon,alt) coordinates for specified (x,y,z).
   * @param x the x coordinate.
   * @param y the y coordinate.
   * @param z the z coordinate.
   * @return array of {lat,lon,alt}.
   */
  public double[] lla(double x, double y, double z) {
    double a = _a;
    double b = sqrt(a*a*_one_minus_ee);
    double c = sqrt((a*a-b*b)/(b*b));
    double p = sqrt(x*x+y*y);
    double t = atan2(a*z,b*p);
    double st = sin(t);
    double ct = cos(t);
    double lon = atan2(y,x);
    double lat = atan2((z+c*c*b*st*st*st),(p-_ee*a*ct*ct*ct));
    double slat = sin(lat);
    double clat = cos(lat);
    double d = a/sqrt(1.0-_ee*slat*slat);
    double alt = p/clat-d;
    if (abs(x)<_tinyxy && abs(y)<_tinyxy)
      alt = abs(z)-b;
    return new double[]{lat*R2D,lon*R2D,alt};
  }

  /**
   * Gets the equatorial radius of this ellipsoid.
   * @return the equatorial radius (= semi-major axis).
   */
  public double getRadius() {
    return _a;
  }

  /**
   * Gets the semi-major axis of this ellipsoid.
   * @return the semi-major axis (= equitorial radius).
   */
  public double getA() {
    return _a;
  }

  /**
   * Gets the semi-minor axis of this ellipsoid.
   * @return the semi-minor axis (= polar distance).
   */
  public double getB() {
    return _b;
  }

  /////////////////////////////////////////////////////////////////////////////
  // private

  private static final double D2R = PI/180.0; // degrees to radians
  private static final double R2D = 180.0/PI; // radians to degrees

  private double _a; // semi-major axis = equitorial radius
  private double _b; // semi-minor axis = polar distance
  private double _ee; // eccentricity squared = 1-(b*b)/(a*a)
  private double _one_minus_ee; // (b*b)/(a*a);
  private double _tinyxy; // improves accuracy near poles
}
