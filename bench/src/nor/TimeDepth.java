package nor;

import static java.lang.Math.*;
import java.util.*;
import java.io.FileInputStream;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.io.ArrayOutputStream;

/**
 * Time-depth conversion for Norne dataset.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.31
 */
public class TimeDepth {
  public static void main(String[] args) {
    makeTimeDepth(
      Coordinates.SX,
      Coordinates.SY,
      Coordinates.SZ,
      Coordinates.ST);
  }
  public static void makeTimeDepth(
    Sampling sx, Sampling sy, Sampling sz, Sampling st)
  {
    String fileName = "/data/seis/nor/norne/timedepth/";
    fileName += "St0103_Norne_2003_depthconversion.avf.txt";
    int nx = sx.getCount();
    int ny = sy.getCount();
    int nz = sz.getCount();
    int nt = st.getCount();
    float[][][] z = fillfloat(-0.001f,nt,ny,nx);
    float[][][] t = fillfloat(-0.001f,nz,ny,nx);
    try {
      FileInputStream fis = new FileInputStream(fileName);
      Scanner s = new Scanner(fis);
      for (int i=0; i<6; ++i)
        s.nextLine();
      ArrayList<float[]> tlist = new ArrayList<float[]>();
      ArrayList<float[]> zlist = new ArrayList<float[]>();
      while (s.hasNextLine()) {
        FloatList tl = new FloatList();
        FloatList zl = new FloatList();
        double xe = 0.0;
        double yn = 0.0;
        while (s.hasNext() && s.next().startsWith("Function")) {
          xe = s.nextDouble();
          yn = s.nextDouble();
          double time = 0.001*s.nextDouble();
          double vavg = 0.001*s.nextDouble();
          double depth = 0.5*vavg*time;
          tl.add((float)time);
          zl.add((float)depth);
        }
        float[] ts = tl.trim();
        float[] zs = zl.trim();
        if (ts!=null && zs!=null) {
          Coordinates.Norne cn = new Coordinates.Norne(xe,yn,0.0);
          Coordinates.Csm cc = new Coordinates.Csm(cn);
          int jx = sx.indexOfNearest(cc.x3);
          int jy = sy.indexOfNearest(cc.x2);
          if (0<jx && jx<nx-1 && 0<jy && jy<ny-1) {
            System.out.println("jx="+jx+" jy="+jy);
            CubicInterpolator ciz = 
              new CubicInterpolator(CubicInterpolator.Method.LINEAR,ts,zs);
            CubicInterpolator cit = 
              new CubicInterpolator(CubicInterpolator.Method.LINEAR,zs,ts);
            for (int jt=0; jt<nt; ++jt) {
              float tj = (float)st.getValue(jt);
              float zi = ciz.interpolate(tj);
              z[jx][jy][jt] = zi;
            }
            for (int jz=0; jz<nz; ++jz) {
              float zj = (float)sz.getValue(jz);
              float ti = ciz.interpolate(zj);
              t[jx][jy][jz] = ti;
            }
          }
        }
        s.nextLine();
      }
      System.out.println("writing zoft and tofz ...");
      String zfileName = "/data/seis/nor/csm/dat/zoft.dat";
      String tfileName = "/data/seis/nor/csm/dat/tofz.dat";
      ArrayOutputStream aosz = new ArrayOutputStream(zfileName);
      aosz.writeFloats(z);
      aosz.close();
      ArrayOutputStream aost = new ArrayOutputStream(tfileName);
      aost.writeFloats(t);
      aost.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }
}
