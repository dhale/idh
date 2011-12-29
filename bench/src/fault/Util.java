package fault;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Utilities.
 */
public class Util {

  public static float[][][][] fakeFpt(int n1, int n2, int n3) {
    float[][][][] fpt = new float[3][n3][n2][n1];
    fpt = mergeFpt(fpt,planeFpt(n1,n2,n3,1.0f, 20.0f,10.0f,n1/2,n2/2,n3/3));
    fpt = mergeFpt(fpt,planeFpt(n1,n2,n3,0.8f,-40.0f,10.0f,n1/2,n2/2,n3/3));
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
          float d = (i1-x1)*u1+(i2-x2)*u2+(i3-x3)*u3;
          f[i3][i2][i1] = fl*exp(-0.125f*d*d);
          p[i3][i2][i1] = fp;
          t[i3][i2][i1] = ft;
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
}
