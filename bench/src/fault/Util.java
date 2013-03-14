package fault;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Utilities.
 */
public class Util {

  public static class QuadInsideCone implements FaultSurfer3.QuadFilter {
    public QuadInsideCone(
      double x1, double x2, double x3, double h, double r)
    {
      _x1 = x1;
      _x2 = x2;
      _x3 = x3;
      _xh = x1+h;
      _roh = r/h;
    }
    public boolean good(FaultSurfer3.Quad quad) {
      double c1 = quad.c1;
      boolean good = false;
      if (_x1<=c1 && c1<=_xh) {
        double d2 = quad.c2-_x2;
        double d3 = quad.c3-_x3;
        double rc = (c1-_x1)*_roh;
        if (d2*d2+d3*d3<=rc*rc)
          good = true;
      }
      return good;
    }
    private double _x1,_x2,_x3,_xh,_roh;
  }

  public static float[][][][] fakeSpheresFpt(int n1, int n2, int n3) {
    float ra=1.8f*n2, sa=0.0f, c1a= 0.0f*n1, c2a=-0.7f*n2, c3a=-0.7f*n3;
    float rb=1.9f*n2, sb=0.3f, c1b=-0.3f*n1, c2b=-0.9f*n2, c3b=-0.5f*n3;
    float[][][][] fpt = sphereFpt(n1,n2,n3,1.0f,ra,sa,c1a,c2a,c3a);
    fpt = mergeFpt2(fpt,sphereFpt(n1,n2,n3,1.0f,rb,sb,c1b,c2b,c3b));
    return fpt;
  }
  public static float[][][][] fakePlanesFpt(int n1, int n2, int n3) {
    float[][][][] fpt = new float[3][n3][n2][n1];
    fpt = mergeFpt(fpt,planeFpt(n1,n2,n3,1.0f, 40.0f,10.0f,n1/2,n2/2,n3/3));
    fpt = mergeFpt(fpt,planeFpt(n1,n2,n3,0.8f,-15.0f,10.0f,n1/2,n2/2,n3/3));
    return fpt;
  }
  public static float[][][][] planeFpt(
    int n1, int n2, int n3,
    float fl, float fp, float ft,
    float x1, float x2, float x3)
  {
    float pr = toRadians(fp);
    float tr = toRadians(ft);
    float cp = cos(pr); 
    float sp = sin(pr);
    float ct = cos(tr);
    float st = sin(tr);
    float u1 = -st;
    float u2 = -sp*ct;
    float u3 =  cp*ct;
    float[][][] f = new float[n3][n2][n1];
    float[][][] p = new float[n3][n2][n1];
    float[][][] t = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float d = (i1-x1)*u1+(i2-x2)*u2+(i3-x3)*u3+0.001f;
          f[i3][i2][i1] = fl*exp(-0.125f*d*d);
          p[i3][i2][i1] = fp;
          t[i3][i2][i1] = ft;
        }
      }
    }
    return new float[][][][]{f,p,t};
  }
  public static float[][][][] sphereFpt(
    int n1, int n2, int n3,
    float fl, float r, float s,
    float c1, float c2, float c3)
  {
    c1 += 0.001f;
    c2 += 0.001f;
    c3 += 0.001f;
    float[][][] f = new float[n3][n2][n1];
    float[][][] p = new float[n3][n2][n1];
    float[][][] t = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      float d3 = i3-c3;
      for (int i2=0; i2<n2; ++i2) {
        float d2 = i2-c2;
        for (int i1=0; i1<n1; ++i1) {
          float d1 = i1-c1;
          float di = sqrt(d1*d1+d2*d2+d3*d3);
          float pi = -atan(d2/d3);
          float ti = -asin(d1/di);
          di += s*cos(toDegrees(pi));
          f[i3][i2][i1] = fl*exp(-0.005f*(di-r)*(di-r));
          p[i3][i2][i1] = toDegrees(pi);
          t[i3][i2][i1] = toDegrees(ti);
        }
      }
    }
    return new float[][][][]{f,p,t};
  }
  public static float[][][][] mergeFpt(
    float[][][][] fpta, float[][][][] fptb) 
  {
    int n1 = fpta[0][0][0].length;
    int n2 = fpta[0][0].length;
    int n3 = fpta[0].length;
    float[][][] fa = fpta[0], pa = fpta[1], ta = fpta[2];
    float[][][] fb = fptb[0], pb = fptb[1], tb = fptb[2];
    float[][][] fc = new float[n3][n2][n1];
    float[][][] pc = new float[n3][n2][n1];
    float[][][] tc = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fai = fa[i3][i2][i1];
          float fbi = fb[i3][i2][i1];
          if (fa[i3][i2][i1]>=fb[i3][i2][i1]) {
            fc[i3][i2][i1] = fa[i3][i2][i1];
            pc[i3][i2][i1] = pa[i3][i2][i1];
            tc[i3][i2][i1] = ta[i3][i2][i1];
          } else {
            fc[i3][i2][i1] = fb[i3][i2][i1];
            pc[i3][i2][i1] = pb[i3][i2][i1];
            tc[i3][i2][i1] = tb[i3][i2][i1];
          }
        }
      }
    }
    return new float[][][][]{fc,pc,tc};
  }
  public static float[][][][] mergeFpt2(
    float[][][][] fpta, float[][][][] fptb) 
  {
    int n1 = fpta[0][0][0].length;
    int n2 = fpta[0][0].length;
    int n3 = fpta[0].length;
    float[][][] fa = fpta[0], pa = fpta[1], ta = fpta[2];
    float[][][] fb = fptb[0], pb = fptb[1], tb = fptb[2];
    float[][][] fc = new float[n3][n2][n1];
    float[][][] pc = new float[n3][n2][n1];
    float[][][] tc = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fai = fa[i3][i2][i1];
          float fbi = fb[i3][i2][i1];
          if (fai>fbi) {
            fc[i3][i2][i1] = fa[i3][i2][i1];
            pc[i3][i2][i1] = pa[i3][i2][i1];
            tc[i3][i2][i1] = ta[i3][i2][i1];
          } else {
            fc[i3][i2][i1] = fb[i3][i2][i1];
            pc[i3][i2][i1] = pb[i3][i2][i1];
            tc[i3][i2][i1] = tb[i3][i2][i1];
          }
        }
      }
    }
    return new float[][][][]{fc,pc,tc};
  }
}
