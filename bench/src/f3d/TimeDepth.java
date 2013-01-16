package f3d;

import java.io.*;
import java.util.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * Time-depth conversion for F3.
 * <em>Not yet working, because the interval velocity file (includes
 *     time-depth pairs) provided with the F3 dataset is wrong. </em>
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.12.29
 */
public class TimeDepth {

  public TimeDepth(String fileName) {
    _xtz = readOdt(fileName);
    float[][] ts = _xtz[1];
    float[][] zs = _xtz[2];
    trace("ts: min="+min(ts)+" max="+max(ts));
    trace("zs: min="+min(zs)+" max="+max(zs));
    float[] x2s = _xtz[0][0];
    float[] x3s = _xtz[0][1];
    SimplePlot sp = new SimplePlot();
    PointsView pv = sp.addPoints(x2s,x3s);
    pv.setLineStyle(PointsView.Line.NONE);
    pv.setMarkStyle(PointsView.Mark.POINT);
  }
  public float[][][] getT(Sampling s1, Sampling s2, Sampling s3) {
    return getTZ(false,s1,s2,s3);
  }
  public float[][][] getZ(Sampling s1, Sampling s2, Sampling s3) {
    return getTZ(true,s1,s2,s3);
  }
  private float[][][] getTZ(
    boolean getz, Sampling s1, Sampling s2, Sampling s3) 
  {
    float[] x2s = _xtz[0][0];
    float[] x3s = _xtz[0][1];
    float[][] ts = getz?_xtz[1]:_xtz[2];
    float[][] zs = getz?_xtz[2]:_xtz[1];
    int nxs = x2s.length;
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    float d1 = (float)s1.getDelta();
    float f1 = (float)s1.getFirst();
    float[][] zi = new float[nxs][n1];
    float[] ti = rampfloat(f1,d1,n1);
    for (int ixs=0; ixs<nxs; ++ixs) {
      quickSort(ts[ixs]);
      quickSort(zs[ixs]);
      CubicInterpolator ci = new CubicInterpolator(
        CubicInterpolator.Method.LINEAR,ts[ixs].length,ts[ixs],zs[ixs]);
      ci.interpolate(ti,zi[ixs]);
    }
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    PixelsView pv = sp.addPixels(zi);
    pv.setColorModel(ColorMap.PRISM);
    /*
    float[] zss = new float[nxs];
    for (int ixs=0; ixs<nxs; ++ixs)
      zss[ixs] = zi[ixs][3*n1/4];
    NearestGridder2 sg = new NearestGridder2(zss,x2s,x3s);
    float[][] zg = sg.grid(s2,s3);
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    PixelsView pv = sp.addPixels(zg);
    pv.setColorModel(ColorMap.JET);
    */
    float[][][] z = new float[n3][n2][n1];
    if (z!=null) return z;
    SibsonInterpolator2 si = new SibsonInterpolator2(x2s,x3s);
    si.setGradientPower(1.0);
    for (int i3=0; i3<n3; ++i3) {
      float x3i = (float)s3.getValue(i3);
      for (int i2=0; i2<n2; ++i2) {
        float x2i = (float)s2.getValue(i2);
        SibsonInterpolator2.IndexWeight[] iw = 
          si.getIndexWeights(x2i,x3i);
        if (iw!=null) {
          int nw = iw.length;
          for (int jw=0; jw<nw; ++jw) {
            int jxs = iw[jw].index;
            float wxs = iw[jw].weight;
            float[] zij = zi[jxs];
            for (int i1=0; i1<n1; ++i1) {
              z[i3][i2][i1] += zij[i1]*wxs;
            }
          }
        }
      }
    }
    extrapolate23(z);
    return z;
  }
  private void extrapolate23(float[][][] v) {
    int n1 = v[0][0].length;
    int n2 = v[0].length;
    int n3 = v.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        if (v[i3][i2][n1-1]==0.0f) {
          int j2lo = i2-1;
          int j2hi = i2+1;
          while (j2lo>=0 && v[i3][j2lo][n1-1]==0.0f) --j2lo;
          while (j2hi<n2 && v[i3][j2hi][n1-1]==0.0f) ++j2hi;
          int j2 = -1;
          if (0<=j2lo && j2hi<n2) {
            j2 = (abs(i2-j2lo)<abs(i2-j2hi))?j2lo:j2hi;
          } else if (0<=j2lo) {
            j2 = j2lo;
          } else {
            j2 = j2hi;
          }
          if (0<=j2 && j2<n2)
            copy(v[i3][j2],v[i3][i2]);
        }
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
        if (v[i3][i2][n1-1]==0.0f) {
          int j3lo = i3-1;
          int j3hi = i3+1;
          while (j3lo>=0 && v[j3lo][i2][n1-1]==0.0f) --j3lo;
          while (j3hi<n3 && v[j3hi][i2][n1-1]==0.0f) ++j3hi;
          int j3 = -1;
          if (0<=j3lo && j3hi<n3) {
            j3 = (abs(i3-j3lo)<abs(i3-j3hi))?j3lo:j3hi;
          } else if (0<=j3lo) {
            j3 = j3lo;
          } else {
            j3 = j3hi;
          }
          if (0<=j3 && j3<n3)
            copy(v[j3][i2],v[i3][i2]);
        }
      }
    }
  }
  
  private float[][][] _xtz;
  
  /**
   * Returns arrays {{x2,x3},t,z}.
   */
  private static float[][][] readOdt(String fileName) {
    //trace("readOdt");
    FloatList x2s = new FloatList();
    FloatList x3s = new FloatList();
    ArrayList<float[]> ts = new ArrayList<float[]>();
    ArrayList<float[]> zs = new ArrayList<float[]>();
    try {
      FileInputStream fis = new FileInputStream(fileName);
      Scanner s = new Scanner(fis);
      s.nextLine(); // skip first line with headings
      double xeOld = -1.0;
      double ynOld = -1.0;
      FloatList tss = null;
      FloatList zss = null;
      while (s.hasNextLine()) {
        Scanner line = new Scanner(s.nextLine());
        double xe = line.nextDouble(); // easting
        double yn = line.nextDouble(); // northing
        if (xe!=xeOld || yn!=ynOld) { // if new coordinates
          Coordinates.Odt odt = new Coordinates.Odt(xe,yn);
          Coordinates.Csm csm = new Coordinates.Csm(odt);
          x2s.add((float)csm.x2);
          x3s.add((float)csm.x3);
          if (tss!=null) { // if we have t and z for old coordinates
            ts.add(clean(tss.trim()));
            zs.add(clean(zss.trim()));
            if (1280<x2s.n && x2s.n<1290) {
              dump(clean(tss.trim()));
              dump(clean(zss.trim()));
              trace("xe="+xe+" yn="+yn+" x2="+csm.x2+" x3="+csm.x3);
            }
          }
          tss = new FloatList();
          zss = new FloatList();
          xeOld = xe;
          ynOld = yn;
        }
        float t = 0.001f*line.nextFloat(); // time (s)
        line.nextFloat(); // Vrms
        line.nextFloat(); // Vint
        line.nextFloat(); // Vavg
        float z = 0.001f*line.nextFloat(); // depth (km)
        tss.add(t);
        zss.add(z);
      }
      s.close();
      ts.add(clean(tss.trim()));
      zs.add(clean(zss.trim()));
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    float[] x2a = x2s.trim();
    float[] x3a = x3s.trim();
    float[][] ta = ts.toArray(new float[0][0]);
    float[][] za = zs.toArray(new float[0][0]);
    //trace("nx="+x2a.length);
    return new float[][][]{{x2a,x3a},ta,za};
  }
  private static float[] clean(float[] a) {
    return (a.length==5)?new float[]{a[0],a[1],a[3],a[4]}:a;
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
