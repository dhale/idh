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
    double p = PIO2;
    double t = PIO2;

    // Sums and differences of basis filter outputs.
    double f0p1 = f[0]+f[1], f2p3 = f[2]+f[3], f4p5 = f[4]+f[5];
    double f0m1 = f[0]-f[1], f2m3 = f[2]-f[3], f4m5 = f[4]-f[5];

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
    int maxiter = 20;
    double small = 0.001*PI;
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
      //System.out.println("ap="+ap+" bp="+bp);
      //System.out.println("at="+at+" bt="+bt+" ct="+ct);
      //System.out.println("ga="+ga+" gb="+gb+" gc="+gc);

      // Ensure Hessian is positive-definite or negative-definite.
      double det = hpp*htt-htp*htp; // determinant
      System.out.println("niter="+niter+" det="+det+" p="+p+" t="+t);
      System.out.println("gp="+gp+" gt="+gt);
      System.out.println("hpp="+hpp+" htp="+htp+" htt="+htt);
      /*
      double dss = 0.001*(hpp*hpp+htt*htt);
      if (det*det<dss*dss) {
        if (det<0.0)
          det -= dss;
        else
          det += dss;
      }
      */

      // Update phi and theta.
      double odet = 1.0/det;
      double dp = odet*(htt*gp-htp*gt);
      double dt = odet*(hpp*gt-htp*gp);
      p -= dp;
      t -= dt;
      ++niter;
      converged = abs(dp*sint)<=small && abs(dt)<=small;
    }
    System.out.println("niter = "+niter);

    // One critical point for unit vector (a0,b0,c0).
    double p0 = modPi2(p);
    double t0 = modPi(t);
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

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final double PIO2 = PI/2.0;
  private static final double PIO3 = PI/3.0;
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
    ftests[0] = new double[]{0.0,0.0,1.0,1.0,1.0,1.0};
    ftests[0] = new double[]{0.8,0.0,0.8,0.0,1.0,1.0};
    ftests[0] = new double[]{0.1174,-0.1017,0.0535,0.0488,-0.0273,0.3474};
    */
    ftests[0] = new double[]{0.1126,0.0602,0.4663,-0.4853,0.4832,0.1081};
    for (int itest=0; itest<ntest; ++itest) {
      double[] f = ftests[itest];
      double f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5];
      System.out.println("f0="+f0+" f1="+f1+" f2="+f2);
      System.out.println("f3="+f3+" f4="+f4+" f5="+f5);
      double[][] pt = findCriticalPhiTheta(f);
      double p0 = pt[0][0], t0 = pt[0][1];
      double p1 = pt[1][0], t1 = pt[1][1];
      double p2 = pt[2][0], t2 = pt[2][1];
      plot(f,p0,t0,p1,t1,p2,t2);
    }
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testPt();
      }
    });
  }
}
