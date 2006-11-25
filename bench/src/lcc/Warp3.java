
package lcc;

import java.awt.*;
import java.awt.image.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

import util.*;

public class Warp3 {

  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Warp3(args);
      }
    });
  }

  private static final LocalCorrelationFilter.Type SIMPLE =
    LocalCorrelationFilter.Type.SIMPLE;
  private static final LocalCorrelationFilter.Type SYMMETRIC =
    LocalCorrelationFilter.Type.SYMMETRIC;

  private static final LocalCorrelationFilter.Window GAUSSIAN = 
    LocalCorrelationFilter.Window.GAUSSIAN;
  private static final LocalCorrelationFilter.Window RECTANGLE = 
    LocalCorrelationFilter.Window.RECTANGLE;

  private int _fontSize = 24;
  private int _width = 640;
  private int _height = 505;
  private int _widthColorBar = 80;

  private String _pngDir = System.getProperty("png.dir");
  private String _dataDir = "/data";

  int _n1 = 101;
  int _n2 = 103;
  int _n3 = 105;
  private float _d1max = 1.00f;
  private float _d2max = 1.00f;
  private float _d3max = 1.00f;
  private int _lmax = 2;
  private int _lmin = -_lmax;
  private LocalCorrelationFilter.Type _type = SIMPLE; 
  private LocalCorrelationFilter.Window _window = GAUSSIAN; 
  private float _sigma = 8.0f;
  private LocalCorrelationFilter _lcf = 
    new LocalCorrelationFilter(_type,_window,_sigma);
  private Displacement _disp = 
    new GaussianDisplacement(_d1max,_d2max,_d3max,_n1,_n2,_n3);
  private float[][][] _f,_g;

  private Warp3(String[] args) {
    if (_pngDir==null)
      _pngDir = ".";
    initImages();
    doIntro();
    doShift();
  }

  private float[][][] makeImage() {
    float[][][] f = Array.sub(Array.randfloat(_n1,_n2,_n3),0.5f); 
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.5);
    rgf.apply0XX(f,f);
    rgf.applyX0X(f,f);
    rgf.applyXX0(f,f);
    System.out.println("f min="+Array.min(f)+" max="+Array.max(f));
    return f;
  }

  private void initImages() {
    _f = makeImage();
    _g = _disp.warp(_f);
  }

  private void doIntro() {
    World world = new World();
    addSlices(world,_f,0.0f,ColorMap.GRAY);
    addSlices(world,_g,0.0f,ColorMap.GRAY);
    frame(world);
  }

  private void doShift() {
    int n1 = _n1;
    int n2 = _n2;
    int n3 = _n3;
    int l1 = _lmax;
    int l2 = _lmax;
    int l3 = _lmax;
    ShiftFinder _sf = new ShiftFinder(_sigma);
    float[][][] f = Array.copy(_f);
    float[][][] g = Array.copy(_g);
    {
      _sf.whiten(_f,f);
      _sf.whiten(_g,g);
      World world = new World();
      addSlices(world,f,0.0f,ColorMap.GRAY);
      addSlices(world,_f,0.0f,ColorMap.GRAY);
      frame(world);
    }
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] du = new float[n3][n2][n1];
    float[][][] ha = Array.copy(g);
    float[][][] hb = Array.copy(g);

    for (int iter=0; iter<4; ++iter) {
      float[][][] ht;
      _sf.find1(-l1,l1,f,ha,du);
      System.out.println("1: du min="+Array.min(du)+" max="+Array.max(du));
      _sf.shift1(du,ha,hb,u1,u2,u3);
      System.out.println("1: u1 min="+Array.min(u1)+" max="+Array.max(u1));
      ht = ha; ha = hb; hb = ht;
      _sf.find2(-l2,l2,f,ha,du);
      System.out.println("2: du min="+Array.min(du)+" max="+Array.max(du));
      _sf.shift2(du,ha,hb,u1,u2,u3);
      System.out.println("2: u2 min="+Array.min(u2)+" max="+Array.max(u2));
      ht = ha; ha = hb; hb = ht;
      _sf.find3(-l3,l3,f,ha,du);
      System.out.println("3: du min="+Array.min(du)+" max="+Array.max(du));
      _sf.shift3(du,ha,hb,u1,u2,u3);
      System.out.println("3: u3 min="+Array.min(u3)+" max="+Array.max(u3));
      ht = ha; ha = hb; hb = ht;
      World world = new World();
      addSlices(world,u1,_d1max,ColorMap.JET);
      addSlices(world,u2,_d2max,ColorMap.JET);
      addSlices(world,u3,_d3max,ColorMap.JET);
      frame(world);
    }
  }

  private void plot(float[][] f, float clip, String png) {
    PlotPanel panel = panel();
    int n1 = f[0].length;
    int n2 = f.length;
    Sampling s1 = new Sampling(n1,1,0);
    Sampling s2 = new Sampling(n2,1,0);
    if (n1<50 && n2<50) {
      s1 = new Sampling(n1,1,-(n1-1)/2);
      s2 = new Sampling(n2,1,-(n2-1)/2);
    }
    PixelsView pv = panel.addPixels(s1,s2,f);
    if (clip!=0.0f) {
      pv.setClips(-clip,clip);
    } else {
      pv.setPercentiles(0.0f,100.0f);
    }
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    frame(panel,png);
  }

  private void plotu(float[][] u, float clip, String png) {
    PlotPanel panel = panel();
    PixelsView pv = panel.addPixels(u);
    pv.setColorModel(ColorMap.JET);
    if (clip!=0.0f) {
      pv.setClips(-clip,clip);
    } else {
      pv.setPercentiles(1.0f,99.9f);
    }
    frame(panel,png);
  }

  private PlotPanel panel() {
    PlotPanel panel = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
    panel.addColorBar();
    panel.setColorBarWidthMinimum(_widthColorBar);
    return panel;
  }

  private PlotFrame frame(PlotPanel panel, String png) {
    PlotFrame frame = new PlotFrame(panel);
    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    frame.setFontSize(_fontSize);
    frame.setSize(_width,_height);
    frame.setVisible(true);
    //if (png!=null)
    //  frame.paintToPng(200,6,_pngDir+"/"+png+".png");
    return frame;
  }

  private JFrame frame(World world) {
    OrbitView view = new OrbitView(world);
    view.setAxesOrientation(View.AxesOrientation.XRIGHT_YOUT_ZDOWN);
    ViewCanvas canvas = new ViewCanvas(view);
    canvas.setView(view);
    ModeManager mm = new ModeManager();
    mm.add(canvas);
    OrbitViewMode ovm = new OrbitViewMode(mm);
    SelectDragMode sdm = new SelectDragMode(mm);
    JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
    toolBar.setRollover(true);
    toolBar.add(new ModeToggleButton(ovm));
    toolBar.add(new ModeToggleButton(sdm));
    JMenu modeMenu = new JMenu("Mode");
    modeMenu.add(new ModeMenuItem(ovm));
    modeMenu.add(new ModeMenuItem(sdm));
    JMenuBar menuBar = new JMenuBar();
    menuBar.add(modeMenu);
    ovm.setActive(true);
    JFrame frame = new JFrame();
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.setSize(new Dimension(600,600));
    frame.add(canvas,BorderLayout.CENTER);
    frame.add(toolBar,BorderLayout.WEST);
    frame.setJMenuBar(menuBar);
    frame.setVisible(true);
    return frame;
  }
  private void addSlices(
    World world, float[][][] f, float clip, IndexColorModel icm) 
  {
    SimpleFloat3 sf3 = new SimpleFloat3(f);
    int nx = _n3;
    int ny = _n2;
    int nz = _n1;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    double lx = nx-1;
    double ly = ny-1;
    double lz = nz-1;
    Sampling sx = new Sampling(nx);
    Sampling sy = new Sampling(ny);
    Sampling sz = new Sampling(nz);
    Point3 qmin = new Point3(fx,fy,fz);
    Point3 qmax = new Point3(lx,ly,lz);
    Axis[] axes = new Axis[]{Axis.X,Axis.Y,Axis.Z};
    for (int iaxis=0; iaxis<axes.length; ++iaxis) {
      AxisAlignedQuad aaq = new AxisAlignedQuad(axes[iaxis],qmin,qmax);
      AxisAlignedFrame aaf = aaq.getFrame();
      ImagePanel ip = new ImagePanel(sx,sy,sz,sf3);
      if (clip!=0.0f)
        ip.setClips(-clip,clip);
      ip.setColorModel(icm);
      aaf.addChild(ip);
      world.addChild(aaq);
    }
  }

  /**
   * Abstract base class for synthetic displacements.
   * The function u(x) is displacement. The function e(x) is strain.
   * A point x is displaced to a point y(x) = x+u(x).
   * <p>
   * Warping is the computation of the sequence g(y) = f(x(y)).
   * Unwarping is the computation of the sequence f(x) = g(y(x)).
   * <p>
   * For warping, we need the function x(y) = y-u(x(y)) = y-uy(y). We 
   * compute the displacement uy(y) by iteration so that uy(y) = u(x(y).
   * <p>
   * We also define a midpoint m(x) = (x+y(x))/2, and compute the 
   * displacement um(m) = u(x(m)) from u(x) by iteration so that 
   * um(m) = u(x(m)).
   */
  private static abstract class Displacement {
    public abstract double u1(double x1, double x2, double x3);
    public abstract double u2(double x1, double x2, double x3);
    public abstract double u3(double x1, double x2, double x3);
    public double u1x(double x1, double x2, double x3) {
      return u1(x1,x2,x3);
    }
    public double u2x(double x1, double x2, double x3) {
      return u2(x1,x2,x3);
    }
    public double u3x(double x1, double x2, double x3) {
      return u3(x1,x2,x3);
    }
    public double u1m(double m1, double m2, double m3) {
      double u1p;
      double u1m = 0.0;
      double u2m = 0.0;
      double u3m = 0.0;
      do {
        u1p = u1m;
        u1m = u1(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
        u2m = u2(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
        u3m = u3(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
      } while (abs(u1m-u1p)>0.0001);
      return u1m;
    }
    public double u2m(double m1, double m2, double m3) {
      double u2p;
      double u1m = 0.0;
      double u2m = 0.0;
      double u3m = 0.0;
      do {
        u2p = u2m;
        u1m = u1(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
        u2m = u2(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
        u3m = u3(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
      } while (abs(u2m-u2p)>0.0001);
      return u2m;
    }
    public double u3m(double m1, double m2, double m3) {
      double u3p;
      double u1m = 0.0;
      double u2m = 0.0;
      double u3m = 0.0;
      do {
        u3p = u3m;
        u1m = u1(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
        u2m = u2(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
        u3m = u3(m1-0.5*u1m,m2-0.5*u2m,m3-0.5*u3m);
      } while (abs(u3m-u3p)>0.0001);
      return u3m;
    }
    public double u1y(double y1, double y2, double y3) {
      double u1p;
      double u1y = 0.0;
      double u2y = 0.0;
      double u3y = 0.0;
      do {
        u1p = u1y;
        u1y = u1(y1-u1y,y2-u2y,y3-u3y);
        u2y = u2(y1-u1y,y2-u2y,y3-u3y);
        u3y = u3(y1-u1y,y2-u2y,y3-u3y);
      } while (abs(u1y-u1p)>0.0001);
      return u1y;
    }
    public double u2y(double y1, double y2, double y3) {
      double u2p;
      double u1y = 0.0;
      double u2y = 0.0;
      double u3y = 0.0;
      do {
        u2p = u2y;
        u1y = u1(y1-u1y,y2-u2y,y3-u3y);
        u2y = u2(y1-u1y,y2-u2y,y3-u3y);
        u3y = u3(y1-u1y,y2-u2y,y3-u3y);
      } while (abs(u2y-u2p)>0.0001);
      return u2y;
    }
    public double u3y(double y1, double y2, double y3) {
      double u3p;
      double u1y = 0.0;
      double u2y = 0.0;
      double u3y = 0.0;
      do {
        u3p = u3y;
        u1y = u1(y1-u1y,y2-u2y,y3-u3y);
        u2y = u2(y1-u1y,y2-u2y,y3-u3y);
        u3y = u3(y1-u1y,y2-u2y,y3-u3y);
      } while (abs(u3y-u3p)>0.0001);
      return u3y;
    }
    public float[][][] u1x() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double x3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double x2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double x1 = i1;
            u[i3][i2][i1] = (float)u1x(x1,x2,x3);
          }
        }
      }
      return u;
    }
    public float[][][] u2x() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double x3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double x2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double x1 = i1;
            u[i3][i2][i1] = (float)u2x(x1,x2,x3);
          }
        }
      }
      return u;
    }
    public float[][][] u3x() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double x3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double x2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double x1 = i1;
            u[i3][i2][i1] = (float)u3x(x1,x2,x3);
          }
        }
      }
      return u;
    }
    public float[][][] u1m() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double m3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double m2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double m1 = i1;
            u[i3][i2][i1] = (float)u1m(m1,m2,m3);
          }
        }
      }
      return u;
    }
    public float[][][] u2m() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double m3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double m2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double m1 = i1;
            u[i3][i2][i1] = (float)u2m(m1,m2,m3);
          }
        }
      }
      return u;
    }
    public float[][][] u3m() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double m3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double m2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double m1 = i1;
            u[i3][i2][i1] = (float)u3m(m1,m2,m3);
          }
        }
      }
      return u;
    }
    public float[][][] u1y() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double y3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double y2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double y1 = i1;
            u[i3][i2][i1] = (float)u1y(y1,y2,y3);
          }
        }
      }
      return u;
    }
    public float[][][] u2y() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double y3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double y2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double y1 = i1;
            u[i3][i2][i1] = (float)u2y(y1,y2,y3);
          }
        }
      }
      return u;
    }
    public float[][][] u3y() {
      float[][][] u = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double y3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double y2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double y1 = i1;
            u[i3][i2][i1] = (float)u3y(y1,y2,y3);
          }
        }
      }
      return u;
    }
    public float[][][] warp(float[][][] f) {
      SincInterpolator si = new SincInterpolator();
      si.setUniform(_n1,1.0,0.0,_n2,1.0,0.0,_n3,1.0,0.0,f);
      float[][][] g = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double y3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double y2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double y1 = i1;
            double x1 = y1-u1y(y1,y2,y3);
            double x2 = y2-u2y(y1,y2,y3);
            double x3 = y3-u3y(y1,y2,y3);
            g[i3][i2][i1] = si.interpolate(x1,x2,x3);
          }
        }
      }
      return g;
    }
    public float[][][] unwarp(float[][][] g) {
      SincInterpolator si = new SincInterpolator();
      si.setUniform(_n1,1.0,0.0,_n2,1.0,0.0,_n3,1.0,0.0,g);
      float[][][] f = new float[_n3][_n2][_n1];
      for (int i3=0; i3<_n3; ++i3) {
        double x3 = i3;
        for (int i2=0; i2<_n2; ++i2) {
          double x2 = i2;
          for (int i1=0; i1<_n1; ++i1) {
            double x1 = i1;
            double y1 = x1+u1x(x1,x2,x3);
            double y2 = x2+u2x(x1,x2,x3);
            double y3 = x3+u3x(x1,x2,x3);
            f[i3][i2][i1] = si.interpolate(y1,y2,y3);
          }
        }
      }
      return f;
    }
    protected Displacement(int n1, int n2, int n3) {
      _n1 = n1;
      _n2 = n2;
      _n3 = n3;
    }
    private int _n1,_n2,_n3;
  }

  /**
   * Constant (zero-strain) displacement.
   */
  private static class ConstantDisplacement extends Displacement {
    public ConstantDisplacement(
      double u1, double u2, double u3, 
      int n1, int n2, int n3) 
    {
      super(n1,n2,n3);
      _u1 = u1;
      _u2 = u2;
      _u3 = u3;
    }
    public double u1(double x1, double x2, double x3) {
      return _u1;
    }
    public double u2(double x1, double x2, double x3) {
      return _u2;
    }
    public double u3(double x1, double x2, double x3) {
      return _u3;
    }
    private double _u1,_u2,_u3;
  }

  /**
   * Derivative-of-Gaussian displacement.
   */
  private static class GaussianDisplacement extends Displacement {
    public GaussianDisplacement(
      double u1max, double u2max, double u3max,
      int n1, int n2, int n3) 
    {
      super(n1,n2,n3);
      _a1 = (n1-1)/2.0;
      _a2 = (n2-1)/2.0;
      _a3 = (n3-1)/2.0;
      _b1 = _a1/3.0;
      _b2 = _a2/3.0;
      _b3 = _a3/3.0;
      _c1 = u1max*exp(0.5)/_b1;
      _c2 = u2max*exp(0.5)/_b2;
      _c3 = u2max*exp(0.5)/_b3;
    }
    public double u1(double x1, double x2, double x3) {
      double xa1 = x1-_a1;
      double xa2 = x2-_a2;
      double xa3 = x3-_a3;
      return -_c1*xa1*exp(-0.5*((xa1*xa1)/(_b1*_b1) +
                                (xa2*xa2)/(_b2*_b2) +
                                (xa3*xa3)/(_b3*_b3)));
    }
    public double u2(double x1, double x2, double x3) {
      double xa1 = x1-_a1;
      double xa2 = x2-_a2;
      double xa3 = x3-_a3;
      return -_c2*xa2*exp(-0.5*((xa1*xa1)/(_b1*_b1) +
                                (xa2*xa2)/(_b2*_b2) +
                                (xa3*xa3)/(_b3*_b3)));
    }
    public double u3(double x1, double x2, double x3) {
      double xa1 = x1-_a1;
      double xa2 = x2-_a2;
      double xa3 = x3-_a3;
      return -_c3*xa3*exp(-0.5*((xa1*xa1)/(_b1*_b1) +
                                (xa2*xa2)/(_b2*_b2) +
                                (xa3*xa3)/(_b3*_b3)));
    }
    private double _a1,_a2,_a3;
    private double _b1,_b2,_b3;
    private double _c1,_c2,_c3;
  }

  /**
   * Sinusoid displacement.
   */
  private static class SinusoidDisplacement extends Displacement {
    public SinusoidDisplacement(
      double u1max, double u2max, double u3max,
      int n1, int n2, int n3) 
    {
      super(n1,n2,n3);
      double l1 = n1-1;
      double l2 = n2-1;
      double l3 = n3-1;
      _a1 = u1max;
      _a2 = u2max;
      _a3 = u3max;
      _b1 = 2.0*PI/l1;
      _b2 = 2.0*PI/l2;
      _b3 = 2.0*PI/l3;
    }
    public double u1(double x1, double x2, double x3) {
      return _a1*sin(_b1*x1)*sin(0.5*_b2*x2)*sin(0.5*_b3*x3);
    }
    public double u2(double x1, double x2, double x3) {
      return _a2*sin(0.5*_b1*x1)*sin(_b2*x2)*sin(0.5*_b3*x3);
    }
    public double u3(double x1, double x2, double x3) {
      return _a3*sin(0.5*_b1*x1)*sin(0.5*_b2*x2)*sin(_b3*x3);
    }
    private double _a1,_a2,_a3;
    private double _b1,_b2,_b3;
  }
}
