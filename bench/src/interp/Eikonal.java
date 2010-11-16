package interp;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Analytical solutions for linear slowness-squared and linear velocity
 * are those published by Fomel, Luo and Zhao (20XX).
 */
public class Eikonal {

  /**
   * Returns traveltime for linear velocity.
   * For a specified source point (p1,p2) where t(p1,p2) = 0, computes 
   * time t(x1,x2) for velocity v(x1,x2) = v0 + v1*x1 + v2*x2, 
   * where v0, v1 and v2 are specified constants.
   * @param v0 velocity v(0,0) at the origin.
   * @param v1 1st component of gradient of v(x1,x2).
   * @param v2 2nd component of gradient of v(x1,x2).
   * @param p1 1st coordinate of point where time is zero.
   * @param p2 2nd coordinate of point where time is zero.
   * @param x1 1st coordinate of point for which to compute time.
   * @param x2 2nd coordinate of point for which to compute time.
   * @return the traveltime.
   */
  public static double timeForLinearVelocity(
    double v0, double v1, double v2, 
    double p1, double p2, 
    double x1, double x2)
  {
    v0 += v1*p1+v2*p2;
    x1 -= p1;
    x2 -= p2;
    double xs = x1*x1+x2*x2;
    double s0 = 1.0/v0;
    double sx = 1.0/(v0+v1*x1+v2*x2);
    double grads = v1*v1+v2*v2;
    double grada = sqrt(grads);
    return (float)(acosh(1.0+0.5*sx*s0*grads*xs)/grada);
  }

  /**
   * Returns sampled traveltimes for linear velocity.
   * For a specified source point (p1,p2) where t(p1,p2) = 0, computes 
   * times t(x1,x2) for velocity v(x1,x2) = v0 + v1*x1 + v2*x2, 
   * where v0, v1 and v2 are specified constants.
   * @param v0 velocity v(0,0) at the origin.
   * @param v1 1st component of gradient of v(x1,x2).
   * @param v2 2nd component of gradient of v(x1,x2).
   * @param p1 1st coordinate of point where time is zero.
   * @param p2 2nd coordinate of point where time is zero.
   * @param s1 sampling of coordinate x1.
   * @param s2 sampling of coordinate x2.
   * @return array[n2][n1] of traveltimes.
   */
  public static float[][] timesForLinearVelocity(
    double v0, double v1, double v2, double p1, double p2,
    Sampling s1, Sampling s2)
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float[][] t = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      double x2 = s2.getValue(i2);
      for (int i1=0; i1<n1; ++i1) {
        double x1 = s1.getValue(i1);
        t[i2][i1] = (float)timeForLinearVelocity(v0,v1,v2,p1,p2,x1,x2);
      }
    }
    return t;
  }

  /**
   * Returns traveltime for linear slowness-squared.
   * For a specified source point (p1,p2) where t(p1,p2) = 0, computes 
   * time t(x1,x2) for slowness-squared s(x1,x2) = ss0 + ss1*x1 + ss2*x2, 
   * where ss0, ss1 and ss2 are specified constants.
   * @param ss0 slowness-squared s(0,0) at the origin.
   * @param ss1 1st component of gradient of s(x1,x2).
   * @param ss2 2nd component of gradient of s(x1,x2).
   * @param p1 1st coordinate of point where time is zero.
   * @param p2 2nd coordinate of point where time is zero.
   * @param x1 1st coordinate of point for which to compute time.
   * @param x2 2nd coordinate of point for which to compute time.
   * @return the traveltime.
   */
  public static double timeForLinearSloth(
    double ss0, double ss1, double ss2, 
    double p1, double p2, 
    double x1, double x2)
  {
    ss0 += ss1*p1+ss2*p2;
    ss1 *= 0.5;
    ss2 *= 0.5;
    x1 -= p1;
    x2 -= p2;
    double xs = x1*x1+x2*x2;
    double grads = ss1*ss1+ss2*ss2;
    double ssavg = ss0+ss1*x1+ss2*x2;
    double sigma = sqrt((2.0*xs)/(ssavg+sqrt(ssavg*ssavg-grads*xs)));
    return (float)((ssavg-SIXTH*grads*sigma*sigma)*sigma);
  }

  /**
   * Returns sampled traveltimes for linear slowness-squared.
   * For a specified source point (p1,p2) where t(p1,p2) = 0, computes 
   * times t(x1,x2) for slowness-squared s(x1,x2) = ss0 + ss1*x1 + ss2*x2, 
   * where ss0, ss1 and ss2 are specified constants.
   * @param ss0 slowness-squared s(0,0) at the origin.
   * @param ss1 1st component of gradient of s(x1,x2).
   * @param ss2 2nd component of gradient of s(x1,x2).
   * @param p1 1st coordinate of point where time is zero.
   * @param p2 2nd coordinate of point where time is zero.
   * @param s1 sampling of coordinate x1.
   * @param s2 sampling of coordinate x2.
   * @return array[n2][n1] of traveltimes.
   */
  public static float[][] timesForLinearSloth(
    double ss0, double ss1, double ss2, double p1, double p2,
    Sampling s1, Sampling s2)
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float[][] t = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      double x2 = s2.getValue(i2);
      for (int i1=0; i1<n1; ++i1) {
        double x1 = s1.getValue(i1);
        t[i2][i1] = (float)timeForLinearSloth(ss0,ss1,ss2,p1,p2,x1,x2);
      }
    }
    return t;
  }

  public static float covExp(float sigma, float delta, float d) {
    return sigma*sigma*exp(-d/delta);
  }

  public static float covGauss(float sigma, float delta, float d) {
    return sigma*sigma*exp(-0.5f*d*d/(delta*delta));
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static double acosh(double z) {
    return log(z+sqrt(z*z-1.0));
  }

  private static final double SIXTH = 1.0/6.0;
}
