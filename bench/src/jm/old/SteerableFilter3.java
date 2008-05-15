package jm;

import static java.lang.Math.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

public class SteerableFilter3 {

  /**
   * Finds a critical point of the function:
   * <pre>
   * f(a,b,c) = ((a+b)*(a+b)-c*c)*f[0] + ((a-b)*(a-b)-c*c)*f[1] +
   *            ((a+c)*(a+c)-b*b)*f[2] + ((a-c)*(a-c)-b*b)*f[3] +
   *            ((b+c)*(b+c)-a*a)*f[4] + ((b-c)*(b-c)-a*a)*f[5]
   * </pre>
   * where (a,b,c) is a point on the unit-sphere.
   */
  private static void findCriticalAbc(double[] f, double[] abcf) {

    // Smallest significant change to a, b, or c (for convergence test).
    double abcsmall = 1.0e-2;

    // Smallest significant determinant (denominator in Newton's method).
    double detsmall = 1.0e-40;

    // Coefficients for evaluating the function f' = (f-fsum)/2.
    // The function f' can be evaluated more efficiently and has
    // the same critical points (a,b,c) as the specified function f.
    double fsum = f[0]+f[1]+f[2]+f[3]+f[4]+f[5];
    double fab = f[0]-f[1];
    double fac = f[2]-f[3];
    double fbc = f[4]-f[5];
    double faa = f[4]+f[5];
    double fbb = f[2]+f[3];
    double fcc = f[0]+f[1];

    // If necessary, initialize (a,b,c).
    double a = abcf[0], b = abcf[1], c = abcf[2];
    double one = a*a+b*b+c*c;
    if (one<0.999 || 1.001<one) {
      double haa,hbb,hcc,hab,hac,hbc;

      // a = 1.0, b = c = 0.0
      hbb = 2.0*(fbb-faa);
      hcc = 2.0*(fcc-faa);
      hbc = -fbc;
      double da = hbb*hcc-hbc*hbc;
      double dda = da*da;

      // b = 1.0, a = c = 0.0
      haa = 2.0*(faa-fbb);
      hcc = 2.0*(fcc-fbb);
      hac = -fac;
      double db = haa*hcc-hac*hac;
      double ddb = db*db;

      // c = 1.0, a = b = 0.0
      haa = 2.0*(faa-fcc);
      hbb = 2.0*(fbb-fcc);
      hab = -fab;
      double dc = haa*hbb-hab*hab;
      double ddc = dc*dc;

      // Choose the point with largest curvature (determinant-squared).
      if (dda>=ddb && dda>=ddc) {
        a = 1.0; b = 0.0; c = 0.0;
      } else if (ddb>=ddc) {
        a = 0.0; b = 1.0; c = 0.0;
      } else {
        a = 0.0; b = 0.0; c = 1.0;
      }
    }

    // Search for a critical point using Newton's method.
    for (int niter=0; niter<50; ++niter) {

      // We have three cases, depending on which coordinate of the 
      // point (a,b,c) has largest magnitude. For example, if that 
      // coordinate is c, we express c as c(a,b) so that f is a 
      // function of f(a,b). We then compute a Newton update (da,db), 
      // constrained so that after the update a*a+b*b <= 1.0. This 
      // last condition is necessary for a*a+b*b+c*c = 1.0, after updating
      // (a,b,c). We eliminate the component with largest magnitude 
      // (in this case, c), because that choice permits the largest 
      // updates (da,db).
      double aa = a*a;
      double bb = b*b;
      double cc = c*c;
      double da,db,dc;

      // If the coordinate c has largest magnitude, ...
      if (aa<=cc && bb<=cc) {
        double aoc = a/c;
        double boc = b/c;
        double aocs = aoc*aoc;
        double bocs = boc*boc;
        double ga = 2.0*a*(faa-fcc)+a*boc*fbc-c*(1.0-aocs)*fac-b*fab;
        double gb = 2.0*b*(fbb-fcc)+b*aoc*fac-c*(1.0-bocs)*fbc-a*fab;
        double haa = 2.0*(faa-fcc)+boc*(1.0+aocs)*fbc+aoc*(3.0+aocs)*fac;
        double hbb = 2.0*(fbb-fcc)+aoc*(1.0+bocs)*fac+boc*(3.0+bocs)*fbc;
        double hab = boc*(1.0+aocs)*fac+aoc*(1.0+bocs)*fbc-fab;
        double det = haa*hbb-hab*hab;
        if (det<=detsmall && -det<=detsmall) det = detsmall;
        double odet = 1.0/det;
        da = odet*(hbb*ga-hab*gb);
        db = odet*(haa*gb-hab*ga);
        dc = 0.0;
        for (double at=a-da,bt=b-db; at*at+bt*bt>=1.0; at=a-da,bt=b-db) {
          da *= 0.5;
          db *= 0.5;
        }
        a -= da;
        b -= db;
        c = (c>=0.0)?sqrt(1.0-a*a-b*b):-sqrt(1.0-a*a-b*b);
      }

      // Else if the coordinate b has largest magnitude, ...
      else if (aa<=bb) {
        double aob = a/b;
        double cob = c/b;
        double aobs = aob*aob;
        double cobs = cob*cob;
        double ga = 2.0*a*(faa-fbb)+a*cob*fbc-b*(1.0-aobs)*fab-c*fac;
        double gc = 2.0*c*(fcc-fbb)+c*aob*fab-b*(1.0-cobs)*fbc-a*fac;
        double haa = 2.0*(faa-fbb)+cob*(1.0+aobs)*fbc+aob*(3.0+aobs)*fab;
        double hcc = 2.0*(fcc-fbb)+aob*(1.0+cobs)*fab+cob*(3.0+cobs)*fbc;
        double hac = cob*(1.0+aobs)*fab+aob*(1.0+cobs)*fbc-fac;
        double det = haa*hcc-hac*hac;
        if (det<=detsmall && -det<=detsmall) det = detsmall;
        double odet = 1.0/det;
        da = odet*(hcc*ga-hac*gc);
        db = 0.0;
        dc = odet*(haa*gc-hac*ga);
        for (double at=a-da,ct=c-dc; at*at+ct*ct>=1.0; at=a-da,ct=c-dc) {
          da *= 0.5;
          dc *= 0.5;
        }
        a -= da;
        c -= dc;
        b = (b>=0.0)?sqrt(1.0-a*a-c*c):-sqrt(1.0-a*a-c*c);
      }

      // Else if the coordinate a has largest magnitude, ...
      else {
        double boa = b/a;
        double coa = c/a;
        double boas = boa*boa;
        double coas = coa*coa;
        double gb = 2.0*b*(fbb-faa)+b*coa*fac-a*(1.0-boas)*fab-c*fbc;
        double gc = 2.0*c*(fcc-faa)+c*boa*fab-a*(1.0-coas)*fac-b*fbc;
        double hbb = 2.0*(fbb-faa)+coa*(1.0+boas)*fac+boa*(3.0+boas)*fab;
        double hcc = 2.0*(fcc-faa)+boa*(1.0+coas)*fab+coa*(3.0+coas)*fac;
        double hbc = coa*(1.0+boas)*fab+boa*(1.0+coas)*fac-fbc;
        double det = hbb*hcc-hbc*hbc;
        if (det<=detsmall && -det<=detsmall) det = detsmall;
        double odet = 1.0/det;
        da = 0.0;
        db = odet*(hcc*gb-hbc*gc);
        dc = odet*(hbb*gc-hbc*gb);
        for (double bt=b-db,ct=c-dc; bt*bt+ct*ct>=1.0; bt=b-db,ct=c-dc) {
          db *= 0.5;
          dc *= 0.5;
        }
        b -= db;
        c -= dc;
        a = (a>=0.0)?sqrt(1.0-b*b-c*c):-sqrt(1.0-b*b-c*c);
      }

      // Test for convergence.
      if (da<abcsmall && -da<abcsmall &&
          db<abcsmall && -db<abcsmall &&
          dc<abcsmall && -dc<abcsmall)
        break;
    }
    double fabc = fsum+2.0*(fab*a*b+fac*a*c+fbc*b*c-faa*a*a-fbb*b*b-fcc*c*c);
    abcf[0] = a; abcf[1] = b; abcf[2] = c; abcf[3] = fabc;
  }

  /**
   * Finds the minimum of the function:
   * <pre>
   * f(a,b,c) = ((a+b)*(a+b)-c*c)*f[0] + ((a-b)*(a-b)-c*c)*f[1] +
   *            ((a+c)*(a+c)-b*b)*f[2] + ((a-c)*(a-c)-b*b)*f[3] +
   *            ((b+c)*(b+c)-a*a)*f[4] + ((b-c)*(b-c)-a*a)*f[5]
   * </pre>
   * where (a,b,c) is a point on the unit-sphere.
   * Uses the Nelder-Mead method that requires no derivatives.
   * <em>Currently broken without a method to initialize (a,b,c).</em>
   */
  private static void findMinimumAbcNelderMead(double[] f, double[] abcf) {

    // Coefficients for evaluating the function f' = (f-fsum)/2.
    // The function f' can be evaluated more efficiently and has
    // the same minimizing (a,b,c) as the specified function f.
    double fsum = f[0]+f[1]+f[2]+f[3]+f[4]+f[5];
    double faa = f[4]+f[5];
    double fbb = f[2]+f[3];
    double fcc = f[0]+f[1];
    double fbc = f[4]-f[5];
    double fac = f[2]-f[3];
    double fab = f[0]-f[1];

    // Initial the three (a,b,c) and corresponding function values.
    //initMinimumAbc(f,abcf);
    double delta = 0.4;
    double a0 = abcf[0], b0 = abcf[1], c0 = abcf[2];
    double as = a0*a0, bs = b0*b0, cs = c0*c0;
    double a1,b1,c1,a2,b2,c2;
    if (cs>=as && cs>=bs) {
      a1 = a0+delta; b1 = b0; c1 = c0; 
      a2 = a0; b2 = b0+delta; c2 = c0; 
    } else if (bs>=as) {
      a1 = a0+delta; b1 = b0; c1 = c0; 
      a2 = a0; b2 = b0; c2 = c0+delta; 
    } else {
      a1 = a0; b1 = b0+delta; c1 = c0; 
      a2 = a0; b2 = b0; c2 = c0+delta; 
    }
    double r1 = 1.0/sqrt(a1*a1+b1*b1+c1*c1);
    double r2 = 1.0/sqrt(a2*a2+b2*b2+c2*c2);
    a1 *= r1; b1 *= r1; c1 *= r1;
    a2 *= r2; b2 *= r2; c2 *= r2;
    //double a0 = 1.0, b0 = 0.0, c0 = 0.0;
    //double a1 = 0.0, b1 = 1.0, c1 = 0.0;
    //double a2 = 0.0, b2 = 0.0, c2 = 1.0;
    double f0 = fab*a0*b0+fac*a0*c0+fbc*b0*c0 -
                faa*a0*a0-fbb*b0*b0-fcc*c0*c0;
    double f1 = fab*a1*b1+fac*a1*c1+fbc*b1*c1 -
                faa*a1*a1-fbb*b1*b1-fcc*c1*c1;
    double f2 = fab*a2*b2+fac*a2*c2+fbc*b2*c2 -
                faa*a2*a2-fbb*b2*b2-fcc*c2*c2;

    // Iterate until convergence.
    double ptol = 0.999;
    int niter;
    for (niter=0; niter<100; ++niter) {

      // Order three points such that f0 <= f1 <= f2.
      if (f0>f1) {
        double at = a0; a0 = a1; a1 = at;
        double bt = b0; b0 = b1; b1 = bt;
        double ct = c0; c0 = c1; c1 = ct;
        double ft = f0; f0 = f1; f1 = ft;
      }
      if (f1>f2) {
        double at = a1; a1 = a2; a2 = at;
        double bt = b1; b1 = b2; b2 = bt;
        double ct = c1; c1 = c2; c2 = ct;
        double ft = f1; f1 = f2; f2 = ft;
      }
      if (f0>f1) {
        double at = a0; a0 = a1; a1 = at;
        double bt = b0; b0 = b1; b1 = bt;
        double ct = c0; c0 = c1; c1 = ct;
        double ft = f0; f0 = f1; f1 = ft;
      }
      double p01 = a0*a1+b0*b1+c0*c1;
      double p02 = a0*a2+b0*b2+c0*c2;
      if (p01>ptol && p02>ptol)
        break;

      // Middle and reflected points, both normalized to lie on sphere.
      double am = 0.5*(a0+a1);
      double bm = 0.5*(b0+b1);
      double cm = 0.5*(c0+c1);
      double sm = 1.0/sqrt(am*am+bm*bm+cm*cm); 
      am *= sm; bm *= sm; cm *= sm;
      double ss = am*a2+bm*b2+cm*c2;
      double ar = 2.0*ss*am-a2;
      double br = 2.0*ss*bm-b2;
      double cr = 2.0*ss*cm-c2;

      // If not min, but less than max, replace the max.
      double fr = fab*ar*br+fac*ar*cr+fbc*br*cr -
                  faa*ar*ar-fbb*br*br-fcc*cr*cr;
      if (f0<=fr && fr<f1) {
        a2 = ar; b2 = br; c2 = cr; f2 = fr;
        System.out.println("reflect");
      }

      // Else if new min, either expand or use reflected point to replace max.
      else if (fr<f0) {
        double ae = 2.0*ss*ar-am;
        double be = 2.0*ss*br-bm;
        double ce = 2.0*ss*cr-cm;
        double fe = fab*ae*be+fac*ae*ce+fbc*be*ce -
                    faa*ae*ae-fbb*be*be-fcc*ce*ce;
        if (fe<fr) {
          a2 = ae; b2 = be; c2 = ce; f2 = fe;
          System.out.println("expand");
        } else {
          a2 = ar; b2 = br; c2 = cr; f2 = fr;
          System.out.println("reflect");
        }
      }

      // Else if contraction (outside or inside) or shrink, ...
      else {
        boolean shrink = false;

        // If outside, either contract or must shrink.
        if (fr<f2) {
          double ac = 0.5*(ar+am);
          double bc = 0.5*(br+bm);
          double cc = 0.5*(cr+cm);
          double sc = 1.0/sqrt(ac*ac+bc*bc+cc*cc);
          ac *= sc; bc *= sc; cc *= sc;
          double fc = fab*ac*bc+fac*ac*cc+fbc*bc*cc -
                      faa*ac*ac-fbb*bc*bc-fcc*cc*cc;
          if (fc<=fr) {
            a2 = ac; b2 = bc; c2 = cc; f2 = fc;
            System.out.println("contract outside");
          } else {
            shrink = true;
          }
        } 
        
        // Else if inside, either contract or must shrink.
        else {
          double ac = 0.5*(a2+am);
          double bc = 0.5*(b2+bm);
          double cc = 0.5*(c2+cm);
          double sc = 1.0/sqrt(ac*ac+bc*bc+cc*cc);
          ac *= sc; bc *= sc; cc *= sc;
          double fc = fab*ac*bc+fac*ac*cc+fbc*bc*cc -
                      faa*ac*ac-fbb*bc*bc-fcc*cc*cc;
          if (fc<f2) {
            a2 = ac; b2 = bc; c2 = cc; f2 = fc;
            System.out.println("contract inside");
          } else {
            shrink = true;
          }
        }

        // If must shrink, ...
        if (shrink) {
          a1 = 0.5*(a0+a1); b1 = 0.5*(b0+b1); c1 = 0.5*(c0+c1);
          a2 = 0.5*(a0+a2); b2 = 0.5*(b0+b2); c2 = 0.5*(c0+c2);
          double s1 = 1.0/sqrt(a1*a1+b1*b1+c1*c1);
          double s2 = 1.0/sqrt(a2*a2+b2*b2+c2*c2);
          a1 *= s1; b1 *= s1; c1 *= s1;
          a2 *= s2; b2 *= s2; c2 *= s2;
          f1 = fab*a1*b1+fac*a1*c1+fbc*b1*c1 -
               faa*a1*a1-fbb*b1*b1-fcc*c1*c1;
          f2 = fab*a2*b2+fac*a2*c2+fbc*b2*c2 -
               faa*a2*a2-fbb*b2*b2-fcc*c2*c2;
          System.out.println("shrink");
        }
      }
    }
    f0 = fsum+2.0*f0;
    System.out.println("niter="+niter+" f0="+f0);
    abcf[0] = a0; abcf[1] = b0; abcf[2] = c0; abcf[3] = f0;
  }

  public static double[] initCriticalPhiThetaXX(double[] f) {
    double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
    double a = f0+f1+f2+f3;
    double b = f0+f1+f4+f5;
    double c = f2+f3+f4+f5;
    double s = 1.0/sqrt(a*a+b*b+c*c);
    a *= s;
    b *= s;
    c *= s;
    double p = modPi2(atan2(b,a));
    double t = acos(c);
    return new double[]{p,t};
  }
  private static double[] initCriticalPhiThetaX(double[] f) {
    double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
    double f0m1 = f0-f1, f2m3 = f2-f3, f4m5 = f4-f5;
    double f0p1 = f0+f1, f2p3 = f2+f3, f4p5 = f4+f5;
    double dab = 4.0*f0*f1-(f2p3-f4p5)*(f2p3-f4p5);
    double dac = 4.0*f2*f3-(f0p1-f4p5)*(f0p1-f4p5);
    double dbc = 4.0*f4*f5-(f0p1-f2p3)*(f0p1-f2p3);
    double sab = dab*dab;
    double sac = dac*dac;
    double sbc = dbc*dbc;
    double a,b,c,s;
    if (sab>=sac && sab>=sbc) {
      s = 1.0/dab;
      a = s*(f2m3*(f2p3-f0p1-f4p5)-f0m1*f4m5);
      b = s*(f4m5*(f4p5-f0p1-f2p3)-f0m1*f2m3);
      c = 1.0;
    } else if (sac>=sab & sbc>=sab) {
      s = 1.0/dac;
      a = s*(f0m1*(f0p1-f2p3-f4p5)-f2m3*f4m5);
      b = 1.0;
      c = s*(f4m5*(f4p5-f2p3-f0p1)-f2m3*f0m1);
    } else {
      s = 1.0/dbc;
      a = 1.0;
      b = s*(f0m1*(f0p1-f2p3-f4p5)-f2m3*f4m5);
      c = s*(f2m3*(f2p3-f0p1-f4p5)-f0m1*f4m5);
    }
    s = 1.0/sqrt(a*a+b*b+c*c);
    a *= s;
    b *= s;
    c *= s;
    double p = modPi2(atan2(b,a));
    double t = acos(c);
    return new double[]{p,t};
  }

  /**
   * Finds extrema of output values for a 2nd-order steerable filter. 
   * Requires six basis filter outputs (s = sqrt(2)/2):
   * <dl>
   * <dt> f0 <dd> corresponds to basis vector (s, s, 0).
   * <dt> f1 <dd> corresponds to basis vector (s,-s, 0).
   * <dt> f2 <dd> corresponds to basis vector (s, 0, s).
   * <dt> f3 <dd> corresponds to basis vector (s, 0,-s).
   * <dt> f4 <dd> corresponds to basis vector (0, s, s).
   * <dt> f5 <dd> corresponds to basis vector (0, s,-s).
   * </dl>
   * @param array {f0,f1,f2,f3,f4,f5} of basis filter outputs.
   * @return array of critical points {{pa,ta,fa},{pb,tb,fb},{pc,tc,fc}}.
   */
  public static double[][] findCriticalPhiTheta(double[] f) {

    // Initial phi and theta (away from poles).
    double[] pt = initCriticalPhiTheta(f);
    //if (pt!=null) {
    //  double[] ptf = {pt[0],pt[1]};
    //  return new double[][]{ptf,ptf,ptf};
    //}
    double p = pt[0]; // phi
    double t = pt[1]; // theta
    //double p = PIO3; // phi
    //double t = PIO3; // theta

    // Sums and differences of basis filter outputs.
    double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
    double f0p1 = f0+f1, f2p3 = f2+f3, f4p5 = f4+f5;
    double f0m1 = f0-f1, f2m3 = f2-f3, f4m5 = f4-f5;

    // Min and max values, used to determine threshold for Hessian.
    double fmin = f0;
    if (fmin>f1) fmin = f1;
    if (fmin>f2) fmin = f2;
    if (fmin>f3) fmin = f3;
    if (fmin>f4) fmin = f4;
    if (fmin>f5) fmin = f5;
    double fmax = f0;
    if (fmax<f1) fmax = f1;
    if (fmax<f2) fmax = f2;
    if (fmax<f3) fmax = f3;
    if (fmax<f4) fmax = f4;
    if (fmax<f5) fmax = f5;

    // Lower limit on Hessian
    double hsmall = 0.0001*(fmax-fmin);
    double hhsmall = hsmall*hsmall;
    if (hhsmall<1.0e-30) hhsmall = 1.0e-30;
    System.out.println("hhsmall="+hhsmall);

    // Constants for computation of gradient and Hessian.
    double haa = f0p1+f2p3-f4p5;
    double hab = f0m1;
    double hac = f2m3;
    double hba = hab;
    double hbb = f0p1-f2p3+f4p5;
    double hbc = f4m5;
    double hca = hac;
    double hcb = hbc;
    double hcc = f2p3+f4p5-f0p1;

    // Newton's method until converged or too many iterations.
    int niter = 0;
    int maxiter = 50;
    double small = 1.0e-3*PI;
    boolean converged = false;
    while (niter<maxiter && !converged) {

      // Unit vector (a,b,c) and partial derivatives wrt phi and theta.
      double cosp = cos(p), sinp = sin(p);
      double cost = cos(t), sint = sin(t);
      //double a   =  cosp*sint, b   =  sinp*sint, c   =  cost;
      //double at  =  cosp*cost, bt  =  sinp*cost, ct  = -sint;
      //double ap  = -sinp*sint, bp  =  cosp*sint, cp  =   0.0;
      //double att = -cosp*sint, btt = -sinp*sint, ctt = -cost;
      //double atp = -sinp*cost, btp =  cosp*cost, ctp =   0.0;
      //double app = -cosp*sint, bpp = -sinp*sint, cpp =   0.0;
      double a   =  cosp*sint, b   =  sinp*sint, c   =  cost;
      double at  =  cosp*cost, bt  =  sinp*cost, ct  = -sint;
      double ap  =         -b, bp  =          a;
      double att = -bp,        btt = ap,         ctt = -cost;
      double atp = -bt,        btp = at;
      double app = -bp,        bpp = ap;

      // Gradient (via chain rule).
      double ga = haa*a+hab*b+hac*c;
      double gb = hba*a+hbb*b+hbc*c;
      double gc = hca*a+hcb*b+hcc*c;
      double gt = ga*at+gb*bt+gc*ct;
      double gp = ga*ap+gb*bp;

      // Hessian (via chain rule).
      double hat = haa*at+hab*bt+hac*ct;
      double hbt = hba*at+hbb*bt+hbc*ct;
      double hct = hca*at+hcb*bt+hcc*ct;
      double hap = haa*ap+hab*bp;
      double hbp = hba*ap+hbb*bp;
      double htt = hat*at+hbt*bt+hct*ct+ga*att+gb*btt+gc*ctt;
      double htp = hat*ap+hbt*bp       +ga*atp+gb*btp;
      double hpp = hap*ap+hbp*bp       +ga*app+gb*bpp;

      // Update phi and theta.
      double det = hpp*htt-htp*htp; // determinant
      if (det<0.0) det = -det;
      if (det<=hhsmall) det = hhsmall;
      double odet = 1.0/det;
      double dp = odet*(htt*gp-htp*gt);
      double dt = odet*(hpp*gt-htp*gp);
      p -= dp;
      t -= dt;
      converged = abs(dp)*sint<=small && abs(dt)<=small;
      ++niter;
      //System.out.println("niter="+niter+" p="+p+" t="+t);
      //System.out.println("gp="+gp+" gt="+gt);
      //System.out.println("hpp="+hpp+" htp="+htp+" htt="+htt);
    }
    p = modPi2(p);
    t = modPi(t);
    System.out.println("niter="+niter+" p="+p+" t="+t);

    // One critical point for unit vector (a0,b0,c0).
    double p0 = p;
    double t0 = t;
    double cp0 = cos(p0), sp0 = sin(p0);
    double ct0 = cos(t0), st0 = sin(t0);
    double a0 = st0*cp0;
    double b0 = st0*sp0;
    double c0 = ct0;

    // Three unit vectors orthogonal to (a0,b0,c0).
    double pa = p0;
    double ta = modPi(t0+0.5*PI);
    double cpa = cos(pa), spa = sin(pa);
    double cta = cos(ta), sta = sin(ta);
    double xa = sta*cpa;
    double ya = sta*spa;
    double za = cta;
    double[] xyzb = rotate(PIO3,a0,b0,c0,xa,ya,za);
    double xb = xyzb[0], yb = xyzb[1], zb = xyzb[2];
    double[] xyzc = rotate(TWO_PIO3,a0,b0,c0,xa,ya,za);
    double xc = xyzc[0], yc = xyzc[1], zc = xyzc[2];

    // Two critical points with unit vectors (a1,b1,c1) and (a2,b2,c2).
    double fa = eval(f,xa,ya,za);
    double fb = eval(f,xb,yb,zb);
    double fc = eval(f,xc,yc,zc);
    double[][] e = SteerableFilter2.findExtrema(fa,fb,fc);
    double alpha1 = e[0][0], alpha2 = e[1][0];
    double[] xyz1 = rotate(alpha1,a0,b0,c0,xa,ya,za);
    double[] xyz2 = rotate(alpha2,a0,b0,c0,xa,ya,za);
    double a1 = xyz1[0], b1 = xyz1[1], c1 = xyz1[2];
    double a2 = xyz2[0], b2 = xyz2[1], c2 = xyz2[2];
    double p1 = atan2(b1,a1), t1 = acos(c1);
    double p2 = atan2(b2,a2), t2 = acos(c2);
    p1 = modPi2(p1); t1 = modPi(t1);
    p2 = modPi2(p2); t2 = modPi(t2);

    return new double[][]{{p0,t0},{p1,t1},{p2,t2}};
  }

  public static double[][] findCriticalAbcX(double[] f) {

    // Initial components (a,b,c) and Lagrange multiplier l.
    double a = 1.0, b = 0.0, c = 0.0, l = 1.0;

    // Sums and differences of basis filter outputs.
    double f0p1 = f[0]+f[1], f2p3 = f[2]+f[3], f4p5 = f[4]+f[5];
    double f0m1 = f[0]-f[1], f2m3 = f[2]-f[3], f4m5 = f[4]-f[5];

    // Constants for computation of gradient and Hessian.
    double faa = f0p1+f2p3-f4p5;
    double fab = f0m1;
    double fac = f2m3;
    double fba = fab;
    double fbb = f0p1-f2p3+f4p5;
    double fbc = f4m5;
    double fca = fac;
    double fcb = fbc;
    double fcc = f2p3+f4p5-f0p1;

    // Newton's method until converged or too many iterations.
    int niter = 0;
    int maxiter = 10;
    double small = 1.0e-3;
    boolean converged = false;
    double[] h = new double[16];
    while (niter<maxiter && !converged) {

      // Gradient 
      double ga = (faa+l)*a+fab*b+fac*c;
      double gb = fba*a+(fbb+l)*b+fbc*c;
      double gc = fca*a+fcb*b+(fcc+l)*c;
      double gl = 0.5*(a*a+b*b+c*c-1.0);

      // Hessian.
      h[ 0] = l+faa; h[ 1] =   fab; h[ 2] =   fac; h[ 3] =   a;
      h[ 4] =   fab; h[ 5] = l+fbb; h[ 6] =   fbc; h[ 7] =   b;
      h[ 8] =   fac; h[ 9] =   fbc; h[10] = l+fcc; h[11] =   c;
      h[12] =     a; h[13] =     b; h[14] =     c; h[15] = 0.0;
      invert(h,h);

      // Update (a,b,c) and l.
      double da = h[ 0]*ga+h[ 1]*gb+h[ 2]*gc+h[ 3]*gl;
      double db = h[ 4]*ga+h[ 5]*gb+h[ 6]*gc+h[ 7]*gl;
      double dc = h[ 8]*ga+h[ 9]*gb+h[10]*gc+h[11]*gl;
      double dl = h[12]*ga+h[13]*gb+h[14]*gc+h[15]*gl;
      a -= da;
      b -= db;
      c -= dc;
      l -= dl;
      ++niter;
      //converged = abs(da)<=small && abs(db)<=small && abs(dc)<=small;
      converged = -small<=da && da<=small &&
                  -small<=db && db<=small &&
                  -small<=dc && dc<=small;
     
    }
    System.out.println("niter="+niter+" a="+a+" b="+b+" c="+c);

    // One critical point at (a0,b0,c0).
    double a0 = a;
    double b0 = b;
    double c0 = c;

    // Three unit vectors r, s, and t orthogonal to (a0,b0,c0).
    // Angles between vectors are PI/3.
    double aa = a*a;
    double bb = b*b;
    double cc = c*c;
    double ar,br,cr;
    if (aa>=cc && bb>=cc) {
      double sr = 1.0/sqrt(aa+bb);
      ar =  b*sr;
      br = -a*sr;
      cr = 0.0;
    } else if (aa>=bb && cc>=bb) {
      double sr = 1.0/sqrt(aa+cc);
      ar =  c*sr;
      cr = -a*sr;
      br = 0.0;
    } else {
      double sr = 1.0/sqrt(bb+cc);
      br =  c*sr;
      cr = -b*sr;
      ar = 0.0;
    }
    double[] abcs = rotate(PIO3,a0,b0,c0,ar,br,cr);
    double as = abcs[0], bs = abcs[1], cs = abcs[2];
    double[] abct = rotate(TWO_PIO3,a0,b0,c0,ar,br,cr);
    double at = abct[0], bt = abct[1], ct = abct[2];

    // Steered outputs for unit vectors r, s, and t.
    double fr = eval(f,ar,br,cr);
    double fs = eval(f,as,bs,cs);
    double ft = eval(f,at,bt,ct);

    // Two critical points with unit vectors (a1,b1,c1) and (a2,b2,c2).
    double[][] e = SteerableFilter2.findExtrema(fr,fs,ft);
    double alpha1 = e[0][0], alpha2 = e[1][0];
    double[] abc1 = rotate(alpha1,a0,b0,c0,ar,br,cr);
    double[] abc2 = rotate(alpha2,a0,b0,c0,ar,br,cr);
    double a1 = abc1[0], b1 = abc1[1], c1 = abc1[2];
    double a2 = abc2[0], b2 = abc2[1], c2 = abc2[2];
    return new double[][]{{a0,b0,c0},{a1,b1,c1},{a2,b2,c2}};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final double PIO2 = PI/2.0;
  private static final double PIO3 = PI/3.0;
  private static final double PIO4 = PI/4.0;
  private static final double PIO6 = PI/6.0;
  private static final double PIO8 = PI/8.0;
  private static final double PIO32 = PI/32.0;
  private static final double TWO_PI = 2.0*PI;
  private static final double TWO_PIO3 = 2.0*PIO3;

  // Rotates (px,py,pz) by angle a around an axis with unit vector (ux,uy,uz).
  private static double[] rotate(double a, 
    double ux, double uy, double uz,
    double px, double py, double pz) 
  {
    double c = cos(a), r = 1.0-c, s = sin(a);
    double uxs = ux*s, uys = uy*s, uzs = uz*s;
    double uxr = ux*r, uyr = uy*r, uzr = uz*r;
    double uxyr = ux*uyr, uyzr = uy*uzr, uzxr = uz*uxr;
    double qx = (ux*uxr+c)*px + (uxyr-uzs)*py + (uzxr+uys)*pz;
    double qy = (uxyr+uzs)*px + (uy*uyr+c)*py + (uyzr-uxs)*pz;
    double qz = (uzxr-uys)*px + (uyzr+uxs)*py + (uz*uzr+c)*pz;
    return new double[]{qx,qy,qz};
  }

  private static double eval(double[] f, double p, double t) {
    double cp = cos(p), sp = sin(p);
    double ct = cos(t), st = sin(t);
    double a = cp*st, b = sp*st, c = ct;
    return eval(f,a,b,c);
  }

  private static double eval(double[] f, double a, double b, double c) {
    double apb = a+b;
    double amb = a-b;
    double apc = a+c;
    double amc = a-c;
    double bpc = b+c;
    double bmc = b-c;
    double w0 = apb*apb-c*c;
    double w1 = amb*amb-c*c;
    double w2 = apc*apc-b*b;
    double w3 = amc*amc-b*b;
    double w4 = bpc*bpc-a*a;
    double w5 = bmc*bmc-a*a;
    return w0*f[0]+w1*f[1]+w2*f[2]+w3*f[3]+w4*f[4]+w5*f[5];
  }

  // Returns the specified angle modulo PI.
  // The returned angle is in the range [0,PI].
  private static double modPi(double theta) {
    return theta-floor(theta/PI)*PI;
  }
  private static double modPi2(double theta) {
    return theta-floor(theta/(2.0*PI))*2.0*PI;
  }

  private static double[] invert(double[] a, double[] b) {

    // transpose
    double t00 = a[ 0];
    double t01 = a[ 4];
    double t02 = a[ 8];
    double t03 = a[12];
    double t04 = a[ 1];
    double t05 = a[ 5];
    double t06 = a[ 9];
    double t07 = a[13];
    double t08 = a[ 2];
    double t09 = a[ 6];
    double t10 = a[10];
    double t11 = a[14];
    double t12 = a[ 3];
    double t13 = a[ 7];
    double t14 = a[11];
    double t15 = a[15];

    // pairs for first 8 elements (cofactors)
    double u00 = t10*t15;
    double u01 = t11*t14;
    double u02 = t09*t15;
    double u03 = t11*t13;
    double u04 = t09*t14;
    double u05 = t10*t13;
    double u06 = t08*t15;
    double u07 = t11*t12;
    double u08 = t08*t14;
    double u09 = t10*t12;
    double u10 = t08*t13;
    double u11 = t09*t12;

    // first 8 elements (cofactors)
    b[ 0] = u00*t05+u03*t06+u04*t07-u01*t05-u02*t06-u05*t07;
    b[ 1] = u01*t04+u06*t06+u09*t07-u00*t04-u07*t06-u08*t07;
    b[ 2] = u02*t04+u07*t05+u10*t07-u03*t04-u06*t05-u11*t07;
    b[ 3] = u05*t04+u08*t05+u11*t06-u04*t04-u09*t05-u10*t06;
    b[ 4] = u01*t01+u02*t02+u05*t03-u00*t01-u03*t02-u04*t03;
    b[ 5] = u00*t00+u07*t02+u08*t03-u01*t00-u06*t02-u09*t03;
    b[ 6] = u03*t00+u06*t01+u11*t03-u02*t00-u07*t01-u10*t03;
    b[ 7] = u04*t00+u09*t01+u10*t02-u05*t00-u08*t01-u11*t02;

    // pairs for second 8 elements (cofactors)
    u00 = t02*t07;
    u01 = t03*t06;
    u02 = t01*t07;
    u03 = t03*t05;
    u04 = t01*t06;
    u05 = t02*t05;
    u06 = t00*t07;
    u07 = t03*t04;
    u08 = t00*t06;
    u09 = t02*t04;
    u10 = t00*t05;
    u11 = t01*t04;

    // second 8 elements (cofactors)
    b[ 8] = u00*t13+u03*t14+u04*t15-u01*t13-u02*t14-u05*t15;
    b[ 9] = u01*t12+u06*t14+u09*t15-u00*t12-u07*t14-u08*t15;
    b[10] = u02*t12+u07*t13+u10*t15-u03*t12-u06*t13-u11*t15;
    b[11] = u05*t12+u08*t13+u11*t14-u04*t12-u09*t13-u10*t14;
    b[12] = u02*t10+u05*t11+u01*t09-u04*t11-u00*t09-u03*t10;
    b[13] = u08*t11+u00*t08+u07*t10-u06*t10-u09*t11-u01*t08;
    b[14] = u06*t09+u11*t11+u03*t08-u10*t11-u02*t08-u07*t09;
    b[15] = u10*t10+u04*t08+u09*t09-u08*t09-u11*t10-u05*t08;

    // determinant
    double d = t00*b[0]+t01*b[1]+t02*b[2]+t03*b[3];
    
    // inverse
    d = 1.0/d;
    b[ 0] *= d;  b[ 1] *= d;  b[ 2] *= d;  b[ 3] *= d;  
    b[ 4] *= d;  b[ 5] *= d;  b[ 6] *= d;  b[ 7] *= d;  
    b[ 8] *= d;  b[ 9] *= d;  b[10] *= d;  b[11] *= d;  
    b[12] *= d;  b[13] *= d;  b[14] *= d;  b[15] *= d;  
    return b;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Testing.

  private static void plot(double[] f, double p, double t) {
    int np = 201;
    int nt = 101;
    double dp = 2.0*PI/(np-1);
    double dt = PI/(nt-1);
    double fp = 0.0;
    double ft = 0.0;
    float[][] fs = new float[np][nt];
    for (int ip=0; ip<np; ++ip) {
      double pi = fp+ip*dp;
      for (int it=0; it<nt; ++it) {
        double ti = ft+it*dt;
        fs[ip][it] = (float)eval(f,pi,ti);
      }
    }
    Sampling sp = new Sampling(np,dp,fp);
    Sampling st = new Sampling(nt,dt,ft);
    SimplePlot plot = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    plot.setSize(800,420);
    PixelsView pv = plot.addPixels(st,sp,fs);
    pv.setColorModel(ColorMap.JET);
    double[] pa = {p};
    double[] ta = {t};
    PointsView mv = plot.addPoints(ta,pa);
    mv.setMarkColor(Color.WHITE);
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
  }

  private static void plot(double[] f, double[] abc) {
    double a = abc[0], b = abc[1], c = abc[2];
    double p = modPi2(atan2(b,a)), t = acos(c);
    plot(f,p,t);
  }

  private static void plot(double[] f,
    double p0, double t0, double p1, double t1, double p2, double t2) 
  {
    int np = 201;
    int nt = 101;
    double dp = 2.0*PI/(np-1);
    double dt = PI/(nt-1);
    double fp = 0.0;
    double ft = 0.0;
    float[][] fs = new float[np][nt];
    for (int ip=0; ip<np; ++ip) {
      double pi = fp+ip*dp;
      for (int it=0; it<nt; ++it) {
        double ti = ft+it*dt;
        fs[ip][it] = (float)eval(f,pi,ti);
      }
    }
    Sampling sp = new Sampling(np,dp,fp);
    Sampling st = new Sampling(nt,dt,ft);
    SimplePlot plot = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    plot.setSize(800,420);
    PixelsView pv = plot.addPixels(st,sp,fs);
    pv.setColorModel(ColorMap.JET);
    double[] pa = {p0,p1,p2};
    double[] ta = {t0,t1,t2};
    PointsView mv = plot.addPoints(ta,pa);
    mv.setLineStyle(PointsView.Line.NONE);
    mv.setMarkColor(Color.WHITE);
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
  }

  private static void plot(double[] f, double[][] abc) {
    double a0 = abc[0][0], b0 = abc[0][1], c0 = abc[0][2];
    double a1 = abc[1][0], b1 = abc[1][1], c1 = abc[1][2];
    double a2 = abc[2][0], b2 = abc[2][1], c2 = abc[2][2];
    double p0 = modPi2(atan2(b0,a0)), t0 = acos(c0);
    double p1 = modPi2(atan2(b1,a1)), t1 = acos(c1);
    double p2 = modPi2(atan2(b2,a2)), t2 = acos(c2);
    plot(f,p0,t0,p1,t1,p2,t2);
  }

  private static void testPt() {
    Random r = new Random();
    int ntest = 1;
    double[][] ftests = Array.sub(Array.randdouble(6,ntest),0.5);
    /*
    ftests[0] = new double[]{1.0,0.0,0.0,0.0,0.0,0.0};
    ftests[1] = new double[]{0.0,1.0,0.0,0.0,0.0,0.0};
    ftests[2] = new double[]{0.0,0.0,1.0,0.0,0.0,0.0};
    ftests[3] = new double[]{0.0,0.0,0.0,1.0,0.0,0.0};
    ftests[4] = new double[]{0.0,0.0,0.0,0.0,1.0,0.0};
    ftests[5] = new double[]{0.0,0.0,0.0,0.0,0.0,1.0};
    ftests[6] = new double[]{0.0,0.0,1.0,1.0,1.0,1.0};
    ftests[7] = new double[]{0.8,0.0,0.8,0.0,1.0,1.0};
    ftests[8] = new double[]{0.1174,-0.1017,0.0535,0.0488,-0.0273,0.3474};
    ftests[9] = new double[]{0.1126,0.0602,0.4663,-0.4853,0.4832,0.1081};
    ftests[10] = new double[]{-0.1910,0.0813,-0.2290,0.3274,-0.0368,-0.3965};
    ftests[11] = new double[]{0.0582,-0.379,0.2959,-0.4895,-0.4966,0.4590};
    ftests[12] = new double[]{0.1,0.1,0.5,-0.5,0.5,0.1};
    ftests[13] = new double[]{0.2334,0.1102,-0.0485,0.2451,0.2162,-0.0368};
    ftests[0] = new double[]{0.0,0.0,1.0,1.0,1.0,1.0};
    */
    for (int itest=0; itest<ntest; ++itest) {
      double[] f = ftests[itest];
      double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
      System.out.println("f0="+f0+" f1="+f1+" f2="+f2);
      System.out.println("f3="+f3+" f4="+f4+" f5="+f5);
      //double[][] pt = findCriticalPhiThetaNm(f);
      //double[][] pt = findCriticalPhiTheta(f);
      double[][] pt = findCriticalAbcX(f);
      double p0 = pt[0][0], t0 = pt[0][1];
      double p1 = pt[1][0], t1 = pt[1][1];
      double p2 = pt[2][0], t2 = pt[2][1];
      plot(f,p0,t0,p1,t1,p2,t2);
    }
  }

  private static void testMinimumAbc() {
    Random r = new Random();
    int ntest = 5;
    double[][] ftests = Array.sub(Array.randdouble(6,ntest),0.5);
    /*
    ftests[0] = new double[]{1.0,0.0,0.0,0.0,0.0,0.0};
    ftests[1] = new double[]{0.0,1.0,0.0,0.0,0.0,0.0};
    ftests[2] = new double[]{0.0,0.0,1.0,0.0,0.0,0.0};
    ftests[3] = new double[]{0.0,0.0,0.0,1.0,0.0,0.0};
    ftests[4] = new double[]{0.0,0.0,0.0,0.0,1.0,0.0};
    ftests[5] = new double[]{0.0,0.0,0.0,0.0,0.0,1.0};
    ftests[6] = new double[]{0.0,0.0,1.0,1.0,1.0,1.0};
    */
    double[] abcf = new double[4];
    for (int itest=0; itest<ntest; ++itest) {
      double[] f = ftests[itest];
      double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
      System.out.println("f0="+f0+" f1="+f1+" f2="+f2);
      System.out.println("f3="+f3+" f4="+f4+" f5="+f5);
      findCriticalAbc(f,abcf);
      plot(f,abcf);
    }
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //testPt();
        testMinimumAbc();
      }
    });
  }
}
