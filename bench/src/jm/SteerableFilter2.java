package jm;

import static java.lang.Math.*;
import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

public class SteerableFilter2 {

  private static final double THETA0 = 0.0*PI/3.0;
  private static final double THETA1 = 1.0*PI/3.0;
  private static final double THETA2 = 2.0*PI/3.0;
  private static final double ONE_THIRD = 1.0/3.0;
  private static final double FOUR_THIRDS = 4.0/3.0;
  private static final double EIGHT_THIRDS = 8.0/3.0;
  private static final double SQRT_THREE = sqrt(3.0);

  // Returns the specified angle modulo PI.
  // The returned angle is in the range [0,PI].
  private static double modPi(double theta) {
    return theta-floor(theta/PI)*PI;
  }

  // Steered output (zeroth derivative).
  private static double eval0(double f0, double f1, double f2, double theta) {
    double k0 = ONE_THIRD*(1.0+2.0*cos(2.0*(theta-THETA0)));
    double k1 = ONE_THIRD*(1.0+2.0*cos(2.0*(theta-THETA1)));
    double k2 = ONE_THIRD*(1.0+2.0*cos(2.0*(theta-THETA2)));
    return k0*f0+k1*f1+k2*f2;
  }

  // First derivative of steered output.
  private static double eval1(double f0, double f1, double f2, double theta) {
    double k0 = -FOUR_THIRDS*sin(2.0*(theta-THETA0));
    double k1 = -FOUR_THIRDS*sin(2.0*(theta-THETA1));
    double k2 = -FOUR_THIRDS*sin(2.0*(theta-THETA2));
    return k0*f0+k1*f1+k2*f2;
  }

  // Second derivative of steered output.
  private static double eval2(double f0, double f1, double f2, double theta) {
    double k0 = -EIGHT_THIRDS*cos(2.0*(theta-THETA0));
    double k1 = -EIGHT_THIRDS*cos(2.0*(theta-THETA1));
    double k2 = -EIGHT_THIRDS*cos(2.0*(theta-THETA2));
    return k0*f0+k1*f1+k2*f2;
  }

  /**
   * Finds extrema of output values for a 2nd-order steerable filter. 
   * The filter is comprised of three steering filters for angles 0,
   * PI/3 and 2*PI/3.
   * <p>
   * Given three values output from the three steering filters, this method
   * computes two angles and corresponding steered output values such that 
   * the derivative of steered output with respect to angle is zero. 
   * Returned angles are in the range [0,PI].
   * @param f0 the value for steering filter with angle 0*PI/3.
   * @param f1 the value for steering filter with angle 1*PI/3.
   * @param f2 the value for steering filter with angle 2*PI/3.
   * @return array {{theta1,value1},{theta2,value2}}. The angle theta1
   *  corresponds to value1 with maximum absolute value; the angle theta2 
   *  corresponds to value2 with minimum absolute value. Both theta1 and
   *  theta2 are in the range [0,PI] and the absolute difference in these
   *  angles is PI/2.
   */
  public static double[][] findExtrema(double f0, double f1, double f2) {
    double theta1 = 0.5*(PI+atan2(SQRT_THREE*(f1-f2),(2.0*f0-f1-f2)));
    double value1 = eval0(f0,f1,f2,theta1);
    double theta2 = modPi(theta1+0.5*PI);
    double value2 = eval0(f0,f1,f2,theta2);
    if (abs(value1)>=abs(value2)) {
      return new double[][]{{theta1,value1},{theta2,value2}};
    } else {
      return new double[][]{{theta2,value2},{theta1,value1}};
    }
  }

  /**
   * Finds extrema of output values for a 2nd-order steerable filter. 
   * The filter is comprised of three steering filters for angles 0,
   * PI/3 and 2*PI/3.
   * <p>
   * Given three values output from the three steering filters, this method
   * computes two angles and corresponding steered output values such that 
   * the derivative of steered output with respect to angle is zero. 
   * Returned angles are in the range [0,PI].
   * @param f0 the value for steering filter with angle 0*PI/3.
   * @param f1 the value for steering filter with angle 1*PI/3.
   * @param f2 the value for steering filter with angle 2*PI/3.
   * @param stheta a small angle. The search for extrema ends when the 
   *  absolute change in angle theta is less than this small angle.
   * @return array {{theta1,value1},{theta2,value2}}. The angle theta1
   *  corresponds to value1 with maximum absolute value; the angle theta2 
   *  corresponds to value2 with minimum absolute value. Both theta1 and
   *  theta2 are in the range [0,PI] and the absolute difference in these
   *  angles is PI/2.
   */
  private static double[][] findExtremaNewtonsMethod(
    double f0, double f1, double f2, double stheta) 
  {
    // Newton's method with safeguard to ensure convergence.
    double theta = PI/3.0; // initial angle
    double btheta = theta; // bound for change in angle theta
    double dtheta = eval1(f0,f1,f2,theta)/eval2(f0,f1,f2,theta);
    while (abs(dtheta)>stheta) {
      //System.out.println("findExtrema: theta="+theta+" dtheta="+dtheta);
      if (abs(dtheta)>btheta)
        dtheta *= btheta/dtheta;
      theta -= dtheta;
      dtheta = eval1(f0,f1,f2,theta)/eval2(f0,f1,f2,theta);
    }
    double theta1 = modPi(theta);
    double value1 = eval0(f0,f1,f2,theta1);
    double theta2 = modPi(theta1+0.5*PI);
    double value2 = eval0(f0,f1,f2,theta2);
    if (abs(value1)>=abs(value2)) {
      return new double[][]{{theta1,value1},{theta2,value2}};
    } else {
      return new double[][]{{theta2,value2},{theta1,value1}};
    }
  }

  /**
   * Like the other methods above but uses a slow linear search.
   */
  private static double[][] findExtremaSlow(
    double f0, double f1, double f2, double dtheta) 
  {
    double theta1 = 0.0;
    double value1 = 0.0;
    for (double theta=0.0; theta<=PI; theta+=dtheta) {
      double value = eval0(f0,f1,f2,theta);
      if (abs(value)>abs(value1)) {
        value1 = value;
        theta1 = theta;
      }
    }
    double theta2 = modPi(theta1+0.5*PI);
    double value2 = eval0(f0,f1,f2,theta2);
    if (abs(value1)>=abs(value2)) {
      return new double[][]{{theta1,value1},{theta2,value2}};
    } else {
      return new double[][]{{theta2,value2},{theta1,value1}};
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Testing.

  private static void plotSteeredOutput(
    double f0, double f1, double f2, double dtheta) 
  {
    int ntheta = 1+(int)(PI/dtheta);
    float[] g = new float[ntheta];
    for (int itheta=0; itheta<ntheta; ++itheta) {
      double theta= itheta*dtheta;
      g[itheta] = (float)eval0(f0,f1,f2,theta);
    }
    Sampling s = new Sampling(ntheta,dtheta,0.0);
    SimplePlot.asPoints(s,g);
  }

  public static void main(String[] args) {
    double stheta = 0.1*PI/180.0;
    Random r = new Random();
    int ntest = 1;
    double[][] f = ArrayMath.sub(ArrayMath.randdouble(3,ntest),0.5);
    /*
    f[0] = new double[]{1.0,1.0,1.0};
    f[1] = new double[]{1.0,0.0,0.0};
    f[2] = new double[]{0.0,1.0,0.0};
    f[3] = new double[]{0.0,0.0,1.0};
    */
    for (int itest=0; itest<ntest; ++itest) {
      double f0 = f[itest][0];
      double f1 = f[itest][1];
      double f2 = f[itest][2];
      System.out.println("f0="+f0+" f1="+f1+" f2="+f2);
      double[][] tvf = findExtrema(f0,f1,f2);
      System.out.println("fast: ("+tvf[0][0]+","+tvf[0][1]+")");
      System.out.println("      ("+tvf[1][0]+","+tvf[1][1]+")");
      double[][] tvm = findExtremaNewtonsMethod(f0,f1,f2,stheta);
      System.out.println("newt: ("+tvm[0][0]+","+tvm[0][1]+")");
      System.out.println("      ("+tvm[1][0]+","+tvm[1][1]+")");
      double[][] tvs = findExtremaSlow(f0,f1,f2,stheta);
      System.out.println("slow: ("+tvs[0][0]+","+tvs[0][1]+")");
      System.out.println("      ("+tvs[1][0]+","+tvs[1][1]+")");
      plotSteeredOutput(f0,f1,f2,stheta);
    }
  }
}
