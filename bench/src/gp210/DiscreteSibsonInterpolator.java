package gp210;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import javax.swing.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.sgl.test.TestFrame;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A simple discrete approximation of Sibson's natural neighbor interpolation.
 * Given a set of known samples of a function f(x1,x2) for scattered points 
 * x1 and x2, this approximation computes an interpolating function g(x1,x2) 
 * with specified uniform samplings of x1 and x2. The function g(x1,x2)
 * approximates the natural neighbor interpolant proposed by Sibson (1981).
 * <p>
 * In this approximation, all scattered points x1 and x2 are rounded to the 
 * nearest uniformly sampled points. If two or more known samples fall into 
 * a single uniform sample bin, their values f(x1,x2) are averaged to obtain 
 * the value g(x1,x2) for that bin. Values of g(x1,x2) for all other bins are 
 * computed to approximate natural neighbor interpolation.
 * <p>
 * The primary goal of this implementation is simplicity. Like the method 
 * of Park et al. (2006), it requires no Delaunay triangulation or Voronoi 
 * tesselation, and its cost decreases as the number of known samples 
 * increases. Moreover, unlike the method of Park et al., this method uses 
 * no auxilary data structure such as a k-d tree to find nearest known 
 * samples. Computational complexity of this method is within a constant 
 * factor of that of Park et al.
 * <p>
 * Discrete implementations of Sibson's interpolation can produce artifacts
 * (small axis-aligned ridges or valleys) caused by sampling circles on a 
 * rectangular grid. To attenuate these artifacts, this method applies some 
 * number of Gauss-Seidel iterations of bi-Laplacian smoothing to the 
 * interpolated samples, without modifying the known samples.
 * <pre>
 * References: 
 * Park, S.W., L. Linsen, O. Kreylos, J.D. Owens, B. Hamann, 2006,
 * Discrete Sibson interpolation: IEEE Transactions on Visualization 
 * and Computer Graphics,
 * v. 12, 243-253.
 * Sibson, R., 1981, A brief description of natural neighbor interpolation,
 * in V. Barnett, ed., Interpreting Multivariate Data: John Wiley and Sons,
 * 21-36.
 * </pre>
 * @author Dave Hale, Colorado School of Mines 
 * @version 2009.01.29
 */
public class DiscreteSibsonInterpolator {

  /**
   * Constructs an interpolator for the specified samplings.
   * During interpolation, scattered points (x1,x2) that lie outside the 
   * uniformly sampled domain are ignored, while those that lie inside are
   * are rounded to the nearest uniformly sampled coordinates.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   */
  public DiscreteSibsonInterpolator(Sampling s1, Sampling s2) {
    _n1 = s1.getCount();
    _n2 = s2.getCount();
    _d1 = s1.getDelta();
    _d2 = s2.getDelta();
    _f1 = s1.getFirst();
    _f2 = s2.getFirst();
    _n = _n1*_n2;

    // Sample offsets, sorted by increasing distance. These offsets are
    // used in expanding-circle searches for nearest known samples.
    // We tabulate offsets for only one quadrant of a circle, because 
    // offsets for the other three quadrants can be found by symmetry.
    int[] kk = new int[_n];
    short[] k1 = new short[_n];
    short[] k2 = new short[_n];
    float[] ds = new float[_n];
    for (int m2=0,k=0; m2<_n2; ++m2) {
      double x2 = m2*_d2;
      for (int m1=0; m1<_n1; ++m1,++k) {
        double x1 = m1*_d1;
        kk[k] = k;
        k1[k] = (short)m1;
        k2[k] = (short)m2;
        ds[k] = (float)(x1*x1+x2*x2);
      }
    }
    quickIndexSort(ds,kk); // <- the only significant external code!
    _k1 = new short[_n];
    _k2 = new short[_n];
    for (int k=0; k<_n; ++k) {
      int kkk = kk[k];
      _k1[k] = k1[kkk];
      _k2[k] = k2[kkk];
    }
  }

  /**
   * Sets the number of bi-Laplacian smoothing iterations.
   * The default number is 100.
   * @param niter the number of smoothing iterations.
   */
  public void setSmoothingIterations(int niter) {
    _niter = niter;
  }

  /**
   * Applies this interpolator for specified f(x1,x2). Computes a uniformly
   * sampled g(x1,x2) that interpolates the scattered values f(x1,x2).
   * @param x1 array of x1 coordinates for which f(x1,x2) is provided.
   * @param x2 array of x2 coordinates for which f(x1,x2) is provided.
   * @param f array of scattered values f(x1,x2) to be interpolated.
   * @return array of uniformly sampled interpolated values g(x1,x2).
   */
  public float[][] apply(float[] x1, float[] x2, float[] f) {
    float fx1 = (float)_f1;
    float fx2 = (float)_f2;
    float lx1 = (float)(_f1+(_n1-1)*_d1);
    float lx2 = (float)(_f2+(_n2-1)*_d2);
    float od1 = 1.0f/(float)_d1;
    float od2 = 1.0f/(float)_d2;

    // Accumulate known samples into bins, counting the number in each bin.
    int n = x1.length;
    float[][] g = new float[_n2][_n1];
    float[][] c = new float[_n2][_n1];
    for (int i=0; i<n; ++i) {
      float x1i = x1[i];
      float x2i = x2[i];
      if (x1i<fx1 || x1i>lx1) continue; // skip scattered values
      if (x2i<fx2 || x2i>lx2) continue; // that fall out of bounds
      int i1 = (int)(0.5f+(x1i-fx1)*od1);
      int i2 = (int)(0.5f+(x2i-fx2)*od2);
      c[i2][i1] += 1.0f; // count known values accumulated
      g[i2][i1] += f[i]; // accumulate known values
    }

    // Average where more than one known sample per bin.
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (c[i2][i1]>0.0f) {
          g[i2][i1] /= c[i2][i1]; // normalize sum of known sample values
          c[i2][i1] = -c[i2][i1]; // negative counts denote known samples
        }
      }
    }

    // For all uniform sample bins (centers of scattering circles), ...
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {

        // Which bins are inside circle extending to nearest known sample?
        // Determine the function value for that nearest known sample.
        int kn = -1;
        float fn = 0.0f;
        for (int k=0; k<_n && kn<0; ++k) {
          int k1 = _k1[k];
          int k2 = _k2[k];
          for (int m2=0,j2=i2-k2; m2<2; ++m2,j2=i2+k2) {
            if (j2<0 || j2>=_n2) continue;
            for (int m1=0,j1=i1-k1; m1<2; ++m1,j1=i1+k1) {
              if (j1<0 || j1>=_n1) continue;
              if (c[j2][j1]<0.0f) { // if sample is known, ...
                kn = k;
                fn = g[j2][j1];
              }
            }
          }
        }

        // By symmetry, each pair of offset indices (k1,k2) actually 
        // represents four uniformly sampled points. In a rectangular 
        // sampling grid, either four or eight such points lie equidistant 
        // from the origin. If eight offsets, then four of them may 
        // immediately follow the four that we found in the search 
        // above. Here we include these next four, if necessary.
        if (kn<_n-1 && _k1[kn]==_k2[kn+1] && _k2[kn]==_k1[kn+1])
          ++kn;

        // Scatter the nearest function value into all bins inside the circle.
        for (int k=0; k<=kn; ++k) {
          int k1 = _k1[k];
          int k2 = _k2[k];
          for (int m2=0,j2=i2-k2; m2<2; ++m2,j2=i2+k2) {
            if (j2<0 || j2>=_n2) continue;
            for (int m1=0,j1=i1-k1; m1<2; ++m1,j1=i1+k1) {
              if (j1<0 || j1>=_n1) continue;
              if (c[j2][j1]>=0.0f) { // if sample is unknown, ...
                g[j2][j1] += fn;
                c[j2][j1] += 1.0f;
              }
            }
          }
        }
      }
    }

    // Normalize accumulated values by the number scattered into each bin.
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (c[i2][i1]>0.0f)
          g[i2][i1] /= c[i2][i1];
      }
    }

    // Some Gauss-Seidel iterations of bi-Laplacian smoothing to
    // attenuate artifacts in discrete Sibson interpolation.
    int n1m = _n1-1;
    int n2m = _n2-1;
    float a1 =  8.0f/20.0f;
    float a2 = -2.0f/20.0f;
    float a3 = -1.0f/20.0f;
    for (int jiter=0; jiter<_niter; ++jiter) {
      for (int i2=0; i2<_n2; ++i2) {
        int i2m = (i2==0  )?i2:i2-1;
        int i2p = (i2==n2m)?i2:i2+1;
        int i2mm = (i2m==0  )?i2m:i2m-1;
        int i2pp = (i2p==n2m)?i2p:i2p+1;
        for (int i1=0; i1<_n1; ++i1) {
          int i1m = (i1==0  )?i1:i1-1;
          int i1p = (i1==n1m)?i1:i1+1;
          int i1mm = (i1m==0  )?i1m:i1m-1;
          int i1pp = (i1p==n1m)?i1p:i1p+1;
          if (c[i2][i1]>0.0f) {
            float g1 = a1*(g[i2 ][i1m]+g[i2 ][i1p]+g[i2m][i1 ]+g[i2p][i1 ]);
            float g2 = a2*(g[i2m][i1m]+g[i2m][i1p]+g[i2p][i1m]+g[i2p][i1p]);
            float g3 = a3*(g[i2][i1mm]+g[i2][i1pp]+g[i2mm][i1]+g[i2pp][i1]);
            g[i2][i1] = g1+g2+g3;
          }
        }
      }
    }

    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2; // numbers of samples
  private double _d1,_d2; // sampling intervals
  private double _f1,_f2; // first-sample values
  private short[] _k1,_k2; // sample offsets for circular search
  private int _n; // total number of samples = _n1*_n2
  private int _niter = 100; // number of smoothing iterations

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static float[][] readData() {

    // Scan data into one big list of floats.
    ArrayList<Float> al = new ArrayList<Float>();
    try {
      Scanner s = new Scanner(new File("src/gp210/"+DATA_NAME+".txt"));
      while (s.hasNextLine()) {
        al.add(s.nextFloat());
        al.add(s.nextFloat());
        al.add(s.nextFloat());
        s.nextLine();
      }
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }

    // Split data into separate arrays of floats.
    int n = al.size()/3;
    float[] x = new float[n];
    float[] y = new float[n];
    float[] z = new float[n];
    for (int i=0,j=0; i<n; ++i,j+=3) {
      x[i] = al.get(j+0);
      y[i] = al.get(j+1);
      z[i] = al.get(j+2);
    }
    float[][] data = new float[][]{x,y,z};
    return data;
  }

  private static Sampling[] makeSamplings(float[] x, float[] y) {
    float xmin = min(x);
    float xmax = max(x);
    float ymin = min(y);
    float ymax = max(y);
    int nx = 315;
    int ny = 315;
    double fx = 0.0;
    double fy = 0.0;
    double dx = 1.0/(nx-1);
    double dy = 1.0/(ny-1);
    Sampling sx = new Sampling(nx,dx,fx);
    Sampling sy = new Sampling(ny,dy,fy);
    return new Sampling[]{sx,sy};
  }

  private static void plot(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz,
    String title, String png) 
  {
    SimplePlot sp = new SimplePlot();
    sp.setSize(PLOT_WIDTH,PLOT_HEIGHT);
    sp.setTitle(title);
    sp.setHLabel("x (m)");
    sp.setVLabel("y (m)");
    sp.addColorBar();
    PlotPanel pp = sp.getPlotPanel();
    pp.setColorBarWidthMinimum(80);
    pp.setLimits(-0.02,-0.02,1.02,1.02);
    PixelsView iv = sp.addPixels(sx,sy,sz);
    iv.setInterpolation(PixelsView.Interpolation.LINEAR);
    iv.setColorModel(ColorMap.PRISM);
    PointsView dv = pp.addPoints(x,y,z);
    dv.setLineStyle(PointsView.Line.NONE);
    dv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    dv.setMarkSize(2.0f);
    dv.setTextFormat("%2.0f");
    sp.paintToPng(300,6,png+"_"+DATA_NAME+".png");
  }

  private static void plot3d(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz)
  {
    sz = mul(0.01f,sz);
    z = mul(0.01f,z);
    PointGroup pg = makePointGroup(x,y,z);
    TriangleGroup tg = makeTriangleGroup(sx,sy,sz);
    World world = new World();
    world.addChild(pg);
    world.addChild(tg);
    TestFrame frame = new TestFrame(world);
    OrbitView view = frame.getOrbitView();
    view.setScale(2.0f);
    view.setAxesOrientation(AxesOrientation.XOUT_YRIGHT_ZUP);
    frame.setSize(new Dimension(1200,800));
    frame.setVisible(true);
  }
  private static PointGroup makePointGroup(float[] x, float[] y, float[] z) {
    int n = x.length;
    float[] xyz = new float[3*n];
    copy(n,0,1,x,0,3,xyz);
    copy(n,0,1,y,1,3,xyz);
    copy(n,0,1,z,2,3,xyz);
    float size = 0.01f;
    PointGroup pg = new PointGroup(size,xyz);
    StateSet states = new StateSet();
    ColorState cs = new ColorState();
    cs.setColor(Color.RED);
    states.add(cs);
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    states.add(lms);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setSpecular(Color.WHITE);
    ms.setShininess(100.0f);
    states.add(ms);
    pg.setStates(states);
    return pg;
  }
  private static TriangleGroup makeTriangleGroup(
    Sampling sx, Sampling sy, float[][] sz) 
  {
    sz = transpose(sz);
    TriangleGroup tg = new TriangleGroup(true,sx,sy,sz);
    StateSet states = new StateSet();
    ColorState cs = new ColorState();
    cs.setColor(Color.LIGHT_GRAY);
    states.add(cs);
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    states.add(lms);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setSpecular(Color.WHITE);
    ms.setShininess(100.0f);
    states.add(ms);
    tg.setStates(states);
    return tg;
  }

  private static void go() {
    float[][] data = readData();
    float[] x = data[0], y = data[1], z = data[2];
    Sampling[] s = makeSamplings(x,y);
    Sampling sx = s[0], sy = s[1];
    int nx = sx.getCount(), ny = sy.getCount();
    DiscreteSibsonInterpolator dsi = new DiscreteSibsonInterpolator(sx,sy);
    dsi.setSmoothingIterations(100);
    float[][] zi = dsi.apply(x,y,z);
    plot(x,y,z,sx,sy,zi,"Discrete Sibson","ds");
    plot3d(x,y,z,sx,sy,zi);
  }

  //public static final String DATA_NAME = "data1";
  public static final String DATA_NAME = "data2";
  //public static final String DATA_NAME = "data5";
  public static final int PLOT_HEIGHT = 785;
  public static final int PLOT_WIDTH = 875;
  public static final int PLOT_MESH_WIDTH = PLOT_WIDTH-132;

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        go();
      }
    });
  }
}
