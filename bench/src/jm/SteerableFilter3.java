package jm;

import static java.lang.Math.*;
import java.awt.*;
import java.util.*;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

public class SteerableFilter3 {

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Smallest significant change to a, b, or c (for convergence test).
  private static final double ABC_SMALL = 1.0e-2;

  // Smallest significant determinant (denominator in Newton's method).
  private static final double DET_SMALL = 1.0e-40;

  // Other constants.
  private static final double COS_PIO3 = cos(PI/3.0);
  private static final double SIN_PIO3 = sin(PI/3.0);
  private static final double SQRT3 = sqrt(3.0);

  /**
   * Finds three critical points of the 3D steering function.
   * The function is
   * <pre>
   * f(a,b,c) = ((a+b)*(a+b)-c*c)*f[0] + ((a-b)*(a-b)-c*c)*f[1] +
   *            ((a+c)*(a+c)-b*b)*f[2] + ((a-c)*(a-c)-b*b)*f[3] +
   *            ((b+c)*(b+c)-a*a)*f[4] + ((b-c)*(b-c)-a*a)*f[5]
   * </pre>
   * where (a,b,c) are components of a unit-vector.
   * @param f array[6] of function values.
   * @param abcf array of critical points; if null, a new array is returned.
   * @return array{{a0,b0,c0,f0},{a1,b1,c1,f1},{a2,b2,c2,f2}} of critical 
   *  points. If a non-null array is specified, that array will be returned
   *  with possibly modified values.
   */
  private static double[][] findCriticalPoints(double[] f, double[][] abcf) {

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

    // Initialize (a,b,c).
    double a,b,c;
    {
      double haa,hbb,hcc,hab,hac,hbc;

      // Determinant for a = 1.0, b = c = 0.0.
      hbb = 2.0*(fbb-faa);
      hcc = 2.0*(fcc-faa);
      hbc = -fbc;
      double da = hbb*hcc-hbc*hbc;
      double dda = da*da;

      // Determinant for b = 1.0, a = c = 0.0.
      haa = 2.0*(faa-fbb);
      hcc = 2.0*(fcc-fbb);
      hac = -fac;
      double db = haa*hcc-hac*hac;
      double ddb = db*db;

      // Determinant for c = 1.0, a = b = 0.0.
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

      // We have three cases, depending on which component of the 
      // vector (a,b,c) has largest magnitude. For example, if that 
      // component is c, we express c as c(a,b) so that f is a 
      // function of f(a,b). We then compute a Newton update (da,db), 
      // constrained so that after the update a*a+b*b <= 1.0. This 
      // last condition is necessary so that a*a+b*b+c*c = 1.0 after 
      // updating (a,b,c). We eliminate the component with largest 
      // magnitude (in this case, c), because that choice permits the 
      // largest updates (da,db).
      double da,db,dc;
      double aa = a*a;
      double bb = b*b;
      double cc = c*c;

      // If the component c has largest magnitude, ...
      if (aa<=cc && bb<=cc) {
        double aoc = a/c;
        double boc = b/c;
        double aocs = aoc*aoc;
        double bocs = boc*boc;
        double faamcc2 = (faa-fcc)*2.0;
        double fbbmcc2 = (fbb-fcc)*2.0;
        double aocsp1 = aocs+1.0;
        double bocsp1 = bocs+1.0;
        double ga = a*(faamcc2+boc*fbc)-c*(1.0-aocs)*fac-b*fab;
        double gb = b*(fbbmcc2+aoc*fac)-c*(1.0-bocs)*fbc-a*fab;
        double haa = faamcc2+boc*aocsp1*fbc+aoc*(3.0+aocs)*fac;
        double hbb = fbbmcc2+aoc*bocsp1*fac+boc*(3.0+bocs)*fbc;
        double hab = boc*aocsp1*fac+aoc*bocsp1*fbc-fab;
        double det = haa*hbb-hab*hab;
        if (det<=DET_SMALL && -det<=DET_SMALL) det = DET_SMALL;
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

      // Else if the component b has largest magnitude, ...
      else if (aa<=bb) {
        double aob = a/b;
        double cob = c/b;
        double aobs = aob*aob;
        double cobs = cob*cob;
        double faambb2 = (faa-fbb)*2.0;
        double fccmbb2 = (fcc-fbb)*2.0;
        double aobsp1 = aobs+1.0;
        double cobsp1 = cobs+1.0;
        double ga = a*(faambb2+cob*fbc)-b*(1.0-aobs)*fab-c*fac;
        double gc = c*(fccmbb2+aob*fab)-b*(1.0-cobs)*fbc-a*fac;
        double haa = faambb2+cob*aobsp1*fbc+aob*(3.0+aobs)*fab;
        double hcc = fccmbb2+aob*cobsp1*fab+cob*(3.0+cobs)*fbc;
        double hac = cob*aobsp1*fab+aob*cobsp1*fbc-fac;
        double det = haa*hcc-hac*hac;
        if (det<=DET_SMALL && -det<=DET_SMALL) det = DET_SMALL;
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

      // Else if the component a has largest magnitude, ...
      else {
        double boa = b/a;
        double coa = c/a;
        double boas = boa*boa;
        double coas = coa*coa;
        double fbbmaa2 = (fbb-faa)*2.0;
        double fccmaa2 = (fcc-faa)*2.0;
        double boasp1 = boas+1.0;
        double coasp1 = coas+1.0;
        double gb = b*(fbbmaa2+coa*fac)-a*(1.0-boas)*fab-c*fbc;
        double gc = c*(fccmaa2+boa*fab)-a*(1.0-coas)*fac-b*fbc;
        double hbb = fbbmaa2+coa*boasp1*fac+boa*(3.0+boas)*fab;
        double hcc = fccmaa2+boa*coasp1*fab+coa*(3.0+coas)*fac;
        double hbc = coa*boasp1*fab+boa*coasp1*fac-fbc;
        double det = hbb*hcc-hbc*hbc;
        if (det<=DET_SMALL && -det<=DET_SMALL) det = DET_SMALL;
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
      if (da<ABC_SMALL && -da<ABC_SMALL &&
          db<ABC_SMALL && -db<ABC_SMALL &&
          dc<ABC_SMALL && -dc<ABC_SMALL)
        break;
    }

    // Assuming convergence, we now have one critical point (a0,b0,c0).
    // We find the other two critical points by searching in the plane
    // orthogonal to the vector p0 = (a0,b0,c0).
    double a0 = a, b0 = b, c0 = c;

    // A second unit vector pr orthogonal to p0.
    double ar,br,cr;
    double aa = a*a, bb = b*b, cc = c*c;
    if (aa<=bb && aa<=cc) {
      ar = 0.0; br = c0; cr = -b0;
    } else if (bb<=cc) {
      ar = c0; br = 0.0; cr = -a0;
    } else {
      ar = b0; br = -a0; cr = 0.0;
    }
    double sr0 = ar*a0+br*b0+cr*c0; // dot product
    ar -= sr0*a0; br -= sr0*b0; cr -= sr0*c0; // Gram-Schmidt
    double sr = 1.0/sqrt(ar*ar+br*br+cr*cr);
    ar *= sr; br*= sr; cr *= sr; // unit vector

    // A third unit vector ps orthogonal to both pr and p0.
    double as = br*c0-b0*cr;
    double bs = cr*a0-c0*ar;
    double cs = ar*b0-a0*br;

    // Three points p1, p2, and p3 in plane of pr and ps,
    // with angles 0, PI/3 and -PI/3 as measured from p1.
    // These are the points used in a 2D steering function.
    double a1 = ar;
    double b1 = br;
    double c1 = cr;
    double a2 = COS_PIO3*ar+SIN_PIO3*as;
    double b2 = COS_PIO3*br+SIN_PIO3*bs;
    double c2 = COS_PIO3*cr+SIN_PIO3*cs;
    double a3 = COS_PIO3*ar-SIN_PIO3*as;
    double b3 = COS_PIO3*br-SIN_PIO3*bs;
    double c3 = COS_PIO3*cr-SIN_PIO3*cs;

    // Function values for points p0, p1, p2, and p3.
    double f0 = fab*a0*b0+fac*a0*c0+fbc*b0*c0-faa*a0*a0-fbb*b0*b0-fcc*c0*c0;
    double f1 = fab*a1*b1+fac*a1*c1+fbc*b1*c1-faa*a1*a1-fbb*b1*b1-fcc*c1*c1;
    double f2 = fab*a2*b2+fac*a2*c2+fbc*b2*c2-faa*a2*a2-fbb*b2*b2-fcc*c2*c2;
    double f3 = fab*a3*b3+fac*a3*c3+fbc*b3*c3-faa*a3*a3-fbb*b3*b3-fcc*c3*c3;

    // Critical points p1 and p2 such that p0, p1, and p2 are orthogonal.
    // We use the analytic solution available for 2D steering functions.
    double denom = 2.0*f1-f2-f3;
    if (denom<DET_SMALL && -denom<DET_SMALL)
      denom = DET_SMALL;
    double theta = 0.5*atan(SQRT3*(f2-f3)/denom);
    double ctheta = cos(theta);
    double stheta = sin(theta);
    a1 = ctheta*ar+stheta*as;
    b1 = ctheta*br+stheta*bs;
    c1 = ctheta*cr+stheta*cs;
    a2 = stheta*ar-ctheta*as;
    b2 = stheta*br-ctheta*bs;
    c2 = stheta*cr-ctheta*cs;

    // Function values for p1 and p2.
    f1 = fab*a1*b1+fac*a1*c1+fbc*b1*c1-faa*a1*a1-fbb*b1*b1-fcc*c1*c1;
    f2 = fab*a2*b2+fac*a2*c2+fbc*b2*c2-faa*a2*a2-fbb*b2*b2-fcc*c2*c2;

    // Convert to the equivalent steering function.
    f0 = fsum+2.0*f0;
    f1 = fsum+2.0*f1;
    f2 = fsum+2.0*f2;

    // Order the three points such that f0 <= f1 <= f2.
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

    // Results.
    if (abcf==null) abcf = new double[3][4];
    abcf[0][0] = a0; abcf[0][1] = b0; abcf[0][2] = c0; abcf[0][3] = f0;
    abcf[1][0] = a1; abcf[1][1] = b1; abcf[1][2] = c1; abcf[1][3] = f1;
    abcf[2][0] = a2; abcf[2][1] = b2; abcf[2][2] = c2; abcf[2][3] = f2;
    return abcf;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  /**
   * Uses the Nelder-Mead simplex method that requires no derivatives.
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

  // Returns the specified angle modulo 2*PI.
  // The returned angle is in the range [0,2PI].
  private static double modPi2(double theta) {
    return theta-floor(theta/(2.0*PI))*2.0*PI;
  }

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

  private static void plot(double[] f, double[][] abcf) {
    double a0 = abcf[0][0], b0 = abcf[0][1], c0 = abcf[0][2];
    double a1 = abcf[1][0], b1 = abcf[1][1], c1 = abcf[1][2];
    double a2 = abcf[2][0], b2 = abcf[2][1], c2 = abcf[2][2];
    double p0 = modPi2(atan2(b0,a0)), t0 = acos(c0);
    double p1 = modPi2(atan2(b1,a1)), t1 = acos(c1);
    double p2 = modPi2(atan2(b2,a2)), t2 = acos(c2);
    plot(f,p0,t0,p1,t1,p2,t2);
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

  private static void testCriticalAbc() {
    Random r = new Random();
    int ntest = 7;
    double[][] ftests = ArrayMath.sub(ArrayMath.randdouble(6,ntest),0.5);
    /*
    ftests[0] = new double[]{1.0,0.0,0.0,0.0,0.0,0.0};
    ftests[1] = new double[]{0.0,1.0,0.0,0.0,0.0,0.0};
    ftests[2] = new double[]{0.0,0.0,1.0,0.0,0.0,0.0};
    ftests[3] = new double[]{0.0,0.0,0.0,1.0,0.0,0.0};
    ftests[4] = new double[]{0.0,0.0,0.0,0.0,1.0,0.0};
    ftests[5] = new double[]{0.0,0.0,0.0,0.0,0.0,1.0};
    ftests[6] = new double[]{0.0,0.0,1.0,1.0,1.0,1.0};
    */
    ftests[0] = new double[]{
      -0.1985327030717713,
       0.09742635963656265,
      -0.3801624013886773,
       0.06162326013941677,
      -0.233385358971406,
      -0.26045753013799466
    };
    ftests[1] = new double[]{
      0.3230484517043609,
      0.4263561654944583,
     -0.4130341047874344,
     -0.004895721621528071,
      0.2025989210726613,
      0.0795547452177835
    };
    for (int itest=0; itest<ntest; ++itest) {
      double[] f = ftests[itest];
      double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
      System.out.println("f0="+f0+" f1="+f1+" f2="+f2);
      System.out.println("f3="+f3+" f4="+f4+" f5="+f5);
      double[][] abcf = findCriticalPoints(f,null);
      plot(f,abcf);
    }
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testCriticalAbc();
      }
    });
  }
}
