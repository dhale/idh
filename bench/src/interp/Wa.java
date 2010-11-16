package interp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.mesh.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Wa {
  public static float[][] gridShepard(
    float[] z, float[] x, float[] y, Sampling sx, Sampling sy)
  {
    int n = x.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    for (int iy=0; iy<ny; ++iy) {
      double yi = sy.getValue(iy);
      for (int ix=0; ix<nx; ++ix) {
        double xi = sx.getValue(ix);
        double num = 0.0;
        double den = 0.0;
        for (int j=0; j<n; ++j) {
          double dx = x[j]-xi;
          double dy = y[j]-yi;
          double ds = dx*dx+dy*dy;
          double wj = 1.0/(ds*ds);
          num += wj*z[j];
          den += wj;
        }
        zi[iy][ix] = (float)(num/den);
      }
    }
    return zi;
  }

  public static final boolean GAUSSIAN = false;
  public static final boolean SUBTREND = true;

  public static float[][] gridKriging(
    float sigma, float delta,
    float[] f, float[] x, float[] y, Sampling sx, Sampling sy)
  {
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
        if (GAUSSIAN) {
          a.set(i,j,covGauss(sigma,delta,xi,yi,xj,yj));
        } else {
          a.set(i,j,covExp(sigma,delta,xi,yi,xj,yj));
        }
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
          if (GAUSSIAN) {
            r.set(i,0,covGauss(sigma,delta,xj,yj,x[i],y[i]));
          } else {
            r.set(i,0,covExp(sigma,delta,xj,yj,x[i],y[i]));
          }
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

  public static float[][] gridKriging(
    float sigma, float delta,
    float[] f, float[] x, float[] y, 
    Sampling sx, Sampling sy, float[][][] tm)
  {
    int n = f.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float dx = (float)sx.getDelta();
    float dy = (float)sy.getDelta();
    DMatrix a = new DMatrix(n,n);
    for (int i=0; i<n; ++i) {
      for (int j=0; j<n; ++j) {
        int jx = sx.indexOfNearest(x[j]);
        int jy = sy.indexOfNearest(y[j]);
        if (GAUSSIAN) {
          a.set(i,j,covGauss(sigma,delta,tm[i][jy][jx]));
        } else {
          a.set(i,j,covExp(sigma,delta,tm[i][jy][jx]));
        }
        //if (i==j) 
        //  a.set(i,j,0.000001*sigma*sigma+a.get(i,j));
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
          if (GAUSSIAN) {
            r.set(i,0,covGauss(sigma,delta,tm[i][jy][jx]));
          } else {
            r.set(i,0,covExp(sigma,delta,tm[i][jy][jx]));
          }
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
  public static float distance(float xa, float ya, float xb, float yb) {
    float dx = xb-xa;
    float dy = yb-ya;
    return sqrt(dx*dx+dy*dy);
  }
  public static float covExp(float sigma, float delta, float d) {
    return sigma*sigma*exp(-d/delta);
  }
  public static float covExp(
    float sigma, float delta, 
    float xa, float ya, float xb, float yb)
  {
    return covExp(sigma,delta,distance(xa,ya,xb,yb));
  }
  public static float covGauss(float sigma, float delta, float d) {
    return sigma*sigma*exp(-0.5f*d*d/(delta*delta));
  }
  public static float covGauss(
    float sigma, float delta, 
    float xa, float ya, float xb, float yb)
  {
    return covGauss(sigma,delta,distance(xa,ya,xb,yb));
  }

  public static float[][] gridPlanar(
    float[] z, float[] x, float[] y, Sampling sx, Sampling sy)
  {
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    TriMesh mesh = makeMesh(x,y,z,sx,sy,true);
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    TriMesh.Tri triLast = null;
    float xa = 0.0f, ya = 0.0f, za = 0.0f, ax = 0.0f, ay = 0.0f;
    float xb,yb,zb,xc,yc,zc;
    float znull = 0.5f*(min(z)+max(z));
    for (int iy=0; iy<ny; ++iy) {
      float yi = (float)sy.getValue(iy);
      for (int ix=0; ix<nx; ++ix) {
        float xi = (float)sx.getValue(ix);
        TriMesh.PointLocation pl = mesh.locatePoint(xi,yi);
        TriMesh.Tri tri = pl.tri();
        if (pl.isInside()) {
          if (tri!=triLast && tri!=null) {
            triLast = tri;
            TriMesh.Node na = tri.nodeA();
            TriMesh.Node nb = tri.nodeB();
            TriMesh.Node nc = tri.nodeC();
            xa = na.x(); ya = na.y(); za = (Float)zmap.get(na);
            xb = nb.x(); yb = nb.y(); zb = (Float)zmap.get(nb);
            xc = nc.x(); yc = nc.y(); zc = (Float)zmap.get(nc);
            xb -= xa; yb -= ya; zb -= za;
            xc -= xa; yc -= ya; zc -= za;
            float odet = 1.0f/(xb*yc-xc*yb);
            ax = (yc*zb-yb*zc)*odet;
            ay = (xb*zc-xc*zb)*odet;
          }
          zi[iy][ix] = za+ax*(xi-xa)+ay*(yi-ya);
        } else {
          zi[iy][ix] = znull;
        }
      }
    }
    return zi;
  }

  private static TriMesh makeMesh(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy, boolean extrap)
  {
    int n = x.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    TriMesh mesh = new TriMesh();
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    for (int i=0; i<n; ++i) {
      TriMesh.Node node = new TriMesh.Node(x[i],y[i]);
      mesh.addNode(node);
      zmap.put(node,new Float(z[i]));
    }
    extrap = true;
    if (extrap) {
      for (int iy=0; iy<ny; iy+=ny-1) {
        float yi = (float)sy.getValue(iy);
        yi += (iy==0)?-0.005f:0.005f;
        for (int ix=0; ix<nx; ix+=nx-1) {
          float xi = (float)sx.getValue(ix);
          xi += (ix==0)?-0.005f:0.005f;
          TriMesh.Node near = mesh.findNodeNearest(xi,yi);
          TriMesh.Node node = new TriMesh.Node(xi,yi);
          mesh.addNode(node);
          zmap.put(node,(Float)zmap.get(near));
        }
      }
    }
    return mesh;
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
}
