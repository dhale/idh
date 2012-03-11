package interp;

import java.util.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.mesh.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;

public class KrigingStudy {
  public static void compare12(
    int n, int nx, int ny,
    double v, double sigmaM, double rangeM, double sigmaD) 
  {
    Sampling sx = new Sampling(nx);
    Sampling sy = new Sampling(ny);
    float[][] fxy = fxy(n,sx,sy);
    float[] f = fxy[0];
    float[] x = fxy[1];
    float[] y = fxy[2];

    DMatrix k = matK(x,y,sx,sy);
    DMatrix kt = k.transpose();
    DMatrix cm = covM(v,sigmaM,rangeM,sx,sy);
    DMatrix cd = covD(n,sigmaD);
    DMatrix cmi = cm.inverse();
    //DMatrix cdi = cd.inverse();
    //DMatrix a1 = kt.times(cdi).times(k).plus(cmi);
    DMatrix a2 = k.times(cm).times(kt).plus(cd);
    //DMatrix b1 = a1.inverse().times(kt).times(cdi);
    DMatrix b2 = cm.times(kt).times(a2.inverse());
    DMatrix z = pack(f);

    //DMatrixSvd svd1 = new DMatrixSvd(a1);
    DMatrixSvd svd2 = new DMatrixSvd(a2);
    //double c1 = svd1.cond();
    double c2 = svd2.cond();
    //trace("c1/c2="+(c1/c2));
    //trace("c1="+c1);
    trace("c2="+c2);

    // CG
    //DMatrix m2 = null;
    //DMatrix m2 = a2.inverse(); // converge in 1 iteration
    DMatrix m2 = k.times(cmi.times(kt));
    DMatrix t2 = solveCg(a2,m2,z);

    // NC
    //DMatrix e = new DMatrix(n,1,1.0);
    //DMatrix s = cm.plus(kt.times(cd.times(k)));
    //double[] dw = s.times(kt.times(e)).getArray();
    //DMatrix w = DMatrix.diagonal(div(1.0,dw));
    //DMatrix mn = w.times(cm.times(kt.times(z)));

    //plot(k,"K");
    //plot(kt,"K'");
    //plot(cm,"Cm");
    //plot(cd,"Cd");
    //plot(cmi,"inv(Cm)");
    //plot(cdi,"inv(Cd)");
    //plot(a1,String.format("A1: cond(A1) = %.4g",c1));
    plot(a2,String.format("A2: cond(A2) = %.4g",c2));
    plot(a2.inverse(),"inv(A2)");
    if (m2!=null)
      plot(m2,"M2");

    float[][] ff = unpack(nx,ny,kt.times(z));
    //float[][] q1 = unpack(nx,ny,b1.times(z));
    //float[][] q1 = new BlendedGridder2(f,y,x).grid(sy,sx);
    float[][] q2 = unpack(nx,ny,b2.times(z));
    float[][] q3 = unpack(nx,ny,cm.times(kt.times(t2)));
    //float[][] qn = unpack(nx,ny,mn);
    float fmin = -0.5f; //min(ff);
    float fmax =  0.5f; //max(ff);
    //plot(q1,fmin,fmax,"q1");
    plot(q2,fmin,fmax,"q2");
    plot(q3,fmin,fmax,"q3");
    //plot(qn,fmin,fmax,"qn");
    plot(ff,fmin,fmax,"f");
  }

  private static DMatrix solveCg(
    final DMatrix ma, final DMatrix mm, DMatrix mb) {
    double[] db = mb.getArray();
    final int n = db.length;
    VecArrayDouble1 vb = new VecArrayDouble1(db);
    VecArrayDouble1 vx = new VecArrayDouble1(n);
    CgSolver cs = new CgSolver(0.001,1000);
    CgSolver.A opa = new CgSolver.A() {
      public void apply(Vec x, Vec y) {
        VecArrayDouble1 vx = (VecArrayDouble1)x;
        VecArrayDouble1 vy = (VecArrayDouble1)y;
        DMatrix mx = new DMatrix(n,1,vx.getArray());
        DMatrix my = ma.times(mx);
        copy(my.getArray(),vy.getArray());
      }
    };
    CgSolver.A opm = null;
    if (mm!=null) {
      opm = new CgSolver.A() {
        public void apply(Vec x, Vec y) {
          VecArrayDouble1 vx = (VecArrayDouble1)x;
          VecArrayDouble1 vy = (VecArrayDouble1)y;
          DMatrix mx = new DMatrix(n,1,vx.getArray());
          DMatrix my = mm.times(mx);
          copy(my.getArray(),vy.getArray());
        }
      };
    }
    DMatrix mx = new DMatrix(n,1,vx.getArray());
    CgSolver.Info info = (mm!=null) ?
      cs.solve(opa,opm,vb,vx) :
      cs.solve(opa,vb,vx);
    trace("solveCg: niter="+info.niter);
    return mx;
  }

  public static DMatrix pack(float[] f) {
    int n = f.length;
    DMatrix z = new DMatrix(n,1);
    for (int i=0; i<n; ++i)
      z.set(i,0,f[i]);
    return z;
  }
  public static float[][] unpack(int nx, int ny, DMatrix a) {
    float[][] q = new float[nx][ny];
    int m = a.getM();
    for (int i=0; i<m; ++i) {
      int ix = i%nx;
      int iy = i/ny;
      q[ix][iy] = (float)a.get(i,0);
    }
    return q;
  }

  public static float[][] fxy(int n, Sampling sx, Sampling sy) {
    int nx = sx.getCount();
    int ny = sy.getCount();
    Random r = new Random(314159);
    float[] f = new float[n];
    float[] x = new float[n];
    float[] y = new float[n];
    boolean[][] filled = new boolean[nx][ny];
    for (int i=0; i<n; ++i) {
      int ix = r.nextInt(nx);
      int iy = r.nextInt(ny);
      while (filled[ix][iy]) {
        ix = r.nextInt(nx);
        iy = r.nextInt(ny);
      }
      filled[ix][iy] = true;
      f[i] = r.nextFloat()-0.5f;
      x[i] = (float)sx.getValue(ix);
      y[i] = (float)sy.getValue(iy);
    }
    return new float[][]{f,x,y};
  }

  public static DMatrix matK(
    float[] x, float[] y, Sampling sx, Sampling sy)
  {
    int nx = sx.getCount();
    int ny = sy.getCount();
    int m = nx*ny;
    int n = x.length;
    DMatrix matk = new DMatrix(n,m);
    for (int j=0; j<n; ++j) {
      int ix = sx.indexOfNearest(x[j]);
      int iy = sy.indexOfNearest(y[j]);
      int i = ix+iy*nx;
      matk.set(j,i,1.0);
    }
    return matk;
  }

  public static final boolean SUBTREND = true;

  public static DMatrix covD(int n, double sigma) {
    return DMatrix.identity(n).times(sigma*sigma);
  }
  public static DMatrix covM(
    double v, double sigma, double range, Sampling sx, Sampling sy)
  {
    Matern cm = new Matern(v,sigma,range);
    int nx = sx.getCount();
    int ny = sy.getCount();
    int m = nx*ny;
    DMatrix covm = new DMatrix(m,m);
    for (int i=0; i<m; ++i) {
      int ix = i%nx;
      int iy = i/nx;
      double xi = sx.getValue(ix);
      double yi = sy.getValue(iy);
      for (int j=0; j<m; ++j) {
        int jx = j%nx;
        int jy = j/nx;
        double xj = sx.getValue(jx);
        double yj = sy.getValue(jy);
        covm.set(i,j,cm.evaluate(distance(xi,yi,xj,yj)));
      }
    }
    return covm;
  }

  public static float[][] gridKriging(
    float sigma, float range,
    float[] f, float[] x, float[] y, Sampling sx, Sampling sy)
  {
    return gridKriging(1.0f,sigma,range,f,x,y,sx,sy);
  }

  public static float[][] gridKriging(
    float v, float sigma, float range,
    float[] f, float[] x, float[] y, Sampling sx, Sampling sy)
  {
    Matern cm = new Matern(v,sigma,range);
    int n = f.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    DMatrix a = new DMatrix(n,n);
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = y[i];
      for (int j=0; j<n; ++j) {
        float xj = x[j];
        float yj = y[j];
        a.set(i,j,cm.evaluate(distance(xi,yi,xj,yj)));
      }
    }
    DMatrixLud lud = new DMatrixLud(a);
    DMatrix identity = DMatrix.identity(n,n);
    DMatrix ai = lud.solve(identity);
    DMatrix r = new DMatrix(n,1);
    float[] ft = null;
    float favg = 0.0f;
    if (SUBTREND) {
      ft = fitTrend(f,x,y);
      f = subTrend(ft,f,x,y);
    } else {
      favg = sum(f)/n; 
      f = sub(f,favg);
    }
    float[][] g = new float[ny][nx];
    for (int jy=0; jy<ny; ++jy) {
      float yj = (float)sy.getValue(jy);
      for (int jx=0; jx<nx; ++jx) {
        float xj = (float)sx.getValue(jx);
        for (int i=0; i<n; ++i) {
          r.set(i,0,cm.evaluate(distance(xj,yj,x[i],y[i])));
        }
        DMatrix w = ai.times(r);
        for (int i=0; i<n; ++i) {
          float wi = (float)w.get(i,0);
          g[jy][jx] += wi*f[i];
        }
        if (SUBTREND) {
          g[jy][jx] += ft[0]+ft[1]*xj+ft[2]*yj;
        } else {
          g[jy][jx] += favg;
        }
      }
    }
    return g;
  }

  public static double distance(double xa, double ya, double xb, double yb) {
    double dx = xb-xa;
    double dy = yb-ya;
    return sqrt(dx*dx+dy*dy);
  }

  public static float[] fitTrend(float[] f, float[] x1, float[] x2) {
    int n = f.length;
    float u1 = sum(x1)/n;
    float u2 = sum(x2)/n;
    x1 = sub(x1,u1);
    x2 = sub(x2,u2);
    float x00 = n;
    float x01 = sum(x1);
    float x02 = sum(x2);
    float x11 = sum(mul(x1,x1));
    float x12 = sum(mul(x1,x2));
    float x22 = sum(mul(x2,x2));
    float r0 = sum(f);
    float r1 = sum(mul(f,x1));
    float r2 = sum(mul(f,x2));
    double[][] aa = { {x00,x01,x02}, {x01,x11,x12}, {x02,x12,x22} };
    double[][] ra = { {r0},{r1},{r2} };
    DMatrix a = new DMatrix(aa);
    DMatrix r = new DMatrix(ra);
    DMatrixLud lud = new DMatrixLud(a);
    DMatrix s = lud.solve(r);
    float s0 = (float)s.get(0,0);
    float s1 = (float)s.get(1,0);
    float s2 = (float)s.get(2,0);
    s0 -= u1*s1+u2*s2;
    return new float[]{s0,s1,s2};
  }

  public static float[] subTrend(
    float[] s, float[] f, float[] x1, float[] x2) 
  {
    float s0 = s[0], s1 = s[1], s2 = s[2];
    return sub(f,add(s0,add(mul(s1,x1),mul(s2,x2))));
  }

  public static float[] addTrend(
    float[] s, float[] f, float[] x1, float[] x2) 
  {
    float s0 = s[0], s1 = s[1], s2 = s[2];
    return add(f,add(s0,add(mul(s1,x1),mul(s2,x2))));
  }

  public static float[][] addTrend(
    float[] t, float[][] f, Sampling s1, Sampling s2)
  {
    float t0 = t[0], t1 = t[1], t2 = t[2];
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float[][] g = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float x2 = (float)s2.getValue(i2);
      for (int i1=0; i1<n1; ++i1) {
        float x1 = (float)s1.getValue(i1);
        g[i2][i1] = f[i2][i1]+t0+t1*x1+t2*x2;
      }
    }
    return g;
  }

  public static void trace(String s) {
    System.out.println(s);
  }

  public static void plot(float[][] a, String title) {
    plot(a,0.0f,0.0f,title);
  }
  public static void plot(float[][] a, float amin, float amax, String title) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    PixelsView pv = sp.addPixels(a);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    if (amin<amax) 
      pv.setClips(amin,amax);
    //pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    pv.setColorModel(ColorMap.JET);
    sp.setSize(790,700);
    sp.setTitle(title);
    sp.addColorBar();
    sp.getPlotPanel().setColorBarWidthMinimum(100);
  }

  public static void plot(DMatrix a, String title) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    PixelsView pv = sp.addPixels(a.transpose().get());
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    //pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    pv.setColorModel(ColorMap.JET);
    sp.setSize(790,700);
    sp.setTitle(title);
    sp.addColorBar();
    sp.getPlotPanel().setColorBarWidthMinimum(100);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
    public void run() {
      go();
    }});
  }

  public static void go() {
    int nx = 25;
    int ny = 25;
    int n = 20;
    double sigmaD = 0.001;
    double vM = 1.0;
    double sigmaM = 1.000;
    double rangeM = 0.500*min(nx,ny);
    compare12(n,nx,ny,vM,sigmaM,rangeM,sigmaD);
  }
}
