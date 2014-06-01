package lcc;

import javax.swing.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.Floats;

public class Warp2 {

  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Warp2(args);
      }
    });
  }

  private static final LocalCorrelationFilter.Type SIMPLE =
    LocalCorrelationFilter.Type.SIMPLE;
  private static final LocalCorrelationFilter.Window GAUSSIAN = 
    LocalCorrelationFilter.Window.GAUSSIAN;
  private int _fontSize = 24;
  private int _width = 640;
  private int _height = 505;
  private int _widthColorBar = 80;

  private String _pngDir = System.getProperty("png.dir");
  private String _dataDir = "/data";

  int _n1 = 315;
  int _n2 = 315;
  private float _d1max = 3.00f;
  private float _d2max = 3.00f;
  private int _lmax = 7;
  private int _lmin = -_lmax;
  private LocalCorrelationFilter.Type _type = SIMPLE; 
  private LocalCorrelationFilter.Window _window = GAUSSIAN; 
  private float _sigma = 12.0f;
  private boolean _whiten = true;
  private boolean _smoothAfterWhiten = true;
  private boolean _randomTestImage = false;
  private LocalCorrelationFilter _lcf = 
    new LocalCorrelationFilter(_type,_window,_sigma);
  private Displacement _disp = 
    new GaussianDisplacement(_d1max,_d2max,_n1,_n2);
  //  new SinusoidDisplacement(_d1max,_d2max,_n1,_n2);
  private float[][] _f,_g;

  private Warp2(String[] args) {
    if (_pngDir==null)
      _pngDir = ".";
    initImages();
    doIntro();
    //doLcc();
    if (_whiten) {
      doLpf();
      //doLcc();
    }
    //doEstimateDisplacements();
    doSequentialShifts();
  }

  private float[][] readImage() {
    int n1 = _n1;
    int n2 = _n2;
    String fileName = _dataDir+"/seis/vg/junks.dat";
    return Floats.readLittleEndian(fileName,n1,n2);
  }

  private float[][] makeImage() {
    float[][] f = mul(sub(randfloat(_n1,_n2),0.5f),20.0f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.5);
    rgf.apply0X(f,f);
    rgf.applyX0(f,f);
    System.out.println("makeImage: f min="+ min(f)+" max="+max(f));
    return f;
  }

  private void initImages() {
    if (_randomTestImage) {
      _f = makeImage();
    } else {
      _f = readImage();
    }
    _g = _disp.warp(_f);
  }

  private void doIntro() {
    plot(_f,6.0f,"imagef");
    plot(_g,6.0f,"imageg");
  }

  private void doSequentialShifts() {
    int n1 = _n1;
    int n2 = _n2;
    int l1 = _lmax;
    int l2 = _lmax;
    ShiftFinder _sf = new ShiftFinder(_sigma);

    float[][] f = copy(_f);
    float[][] g = copy(_g);
    /*
    if (_whiten) {
      _sf.whiten(_f,f);
      _sf.whiten(_g,g);
    }
    */
    plot(f,0.0f,null);
    plot(g,0.0f,null);

    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] du = new float[n2][n1];
    float[][] h = copy(g);

    plot(sub(g,f),0.0f,null);

    for (int iter=0; iter<4; ++iter) {
      _sf.find1(-l1,l1,f,h,du);
      System.out.println("1: du min="+min(du)+" max="+max(du));
      _sf.shift1(du,u1,u2,h);
      System.out.println("1: u1 min="+min(u1)+" max="+max(u1));
      _sf.find2(-l2,l2,f,h,du);
      System.out.println("2: du min="+min(du)+" max="+max(du));
      _sf.shift2(du,u1,u2,h);
      System.out.println("2: u2 min="+min(u2)+" max="+max(u2));
      plotu(u1,_d1max,"u1");
      plotu(u2,_d2max,"u2");
    }

    float[][] e1 = _disp.u1x();
    float[][] e2 = _disp.u2x();
    plotu(e1,_d1max,"e1");
    plotu(e2,_d2max,"e2");
    plot(sub(h,f),0.0f,null);
  }

  private void doLcc() {
    int n1 = _n1;
    int n2 = _n2;
    int l1 = _lmax;
    int l2 = _lmax;
    int m1 = 1+2*l1;
    int m2 = 1+2*l2;
    _lcf.setInputs(_f,_g);
    float[][] c = new float[n2][n1];
    float[][] t = new float[n2][n1];
    int[] k1 = { 38, 98,113,150,98,172};
    int[] k2 = {172,216,230,105,132,227};
    int nk = k1.length;
    float[][][] ck = new float[nk][m2][m1];
    for (int lag2=-l2; lag2<=l2; ++lag2) {
      for (int lag1=-l1; lag1<=l1; ++lag1) {
        _lcf.correlate(lag1,lag2,t);
        _lcf.normalize(lag1,lag2,t);
        copy(n1/m1,n2/m2,l1,l2,m1,m2,t,l1+lag1,l2+lag2,m1,m2,c);
        for (int k=0; k<nk; ++k) {
          ck[k][l2+lag2][l1+lag1] = t[k2[k]][k1[k]];
        }
      }
    }
    plot(c,0.0f,"lcc");
    for (int k=0; k<nk; ++k)
      plot(ck[k],0.0f,"lcc"+k);
  }

  private void doLpf() {
    int n1 = _n1;
    int n2 = _n2;
    LocalPredictionFilter lpf = new LocalPredictionFilter(_sigma);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sqrt(1.0));
    float[][] t = new float[n2][n1];
    /*
    int[] lag1 = {       1,
                  -1, 0, 1};
    int[] lag2 = {       0,
                   1, 1, 1};
    */
    int[] lag1 = {1, 0};
    int[] lag2 = {0, 1};
    lpf.applyPef(lag1,lag2,_f,t);
    copy(t,_f);
    lpf.applyPef(lag1,lag2,_g,t);
    copy(t,_g);
    if (_smoothAfterWhiten) { 
      rgf.apply0X(_f,_f);
      rgf.applyX0(_f,_f);
      rgf.apply0X(_g,_g);
      rgf.applyX0(_g,_g);
    }
    plot(_f,0.6f,"lpeff");
    plot(_g,0.6f,"lpefg");
  }

  private void doEstimateDisplacements() {
    _lcf.setInputs(_f,_g);
    int min1 = _lmin;
    int max1 = _lmax;
    int min2 = _lmin;
    int max2 = _lmax;
    byte[][] l1 = new byte[_n2][_n1];
    byte[][] l2 = new byte[_n2][_n1];
    float[][] u1 = new float[_n2][_n1];
    float[][] u2 = new float[_n2][_n1];
    _lcf.findMaxLags(min1,max1,min2,max2,l1,l2);

    float[][][] q = new float[4][_n2][_n1];
    _lcf.refineLags(l1,l2,u1,u2,q);
    plotu(u1,_d1max,"u1");
    plotu(u2,_d2max,"u2");

    float[][] e1 = (_type==SIMPLE)?_disp.u1x():_disp.u1m();
    float[][] e2 = (_type==SIMPLE)?_disp.u2x():_disp.u2m();
    plotu(e1,_d1max,"e1");
    plotu(e2,_d2max,"e2");
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
    public abstract double u1(double x1, double x2);
    public abstract double u2(double x1, double x2);
    public double u1x(double x1, double x2) {
      return u1(x1,x2);
    }
    public double u2x(double x1, double x2) {
      return u2(x1,x2);
    }
    public double u1m(double m1, double m2) {
      double u1p;
      double u1m = 0.0;
      double u2m = 0.0;
      do {
        u1p = u1m;
        u1m = u1(m1-0.5*u1m,m2-0.5*u2m);
        u2m = u2(m1-0.5*u1m,m2-0.5*u2m);
      } while (abs(u1m-u1p)>0.0001);
      return u1m;
    }
    public double u2m(double m1, double m2) {
      double u2p;
      double u1m = 0.0;
      double u2m = 0.0;
      do {
        u2p = u2m;
        u1m = u1(m1-0.5*u1m,m2-0.5*u2m);
        u2m = u2(m1-0.5*u1m,m2-0.5*u2m);
      } while (abs(u2m-u2p)>0.0001);
      return u2m;
    }
    public double u1y(double y1, double y2) {
      double u1p;
      double u1y = 0.0;
      double u2y = 0.0;
      do {
        u1p = u1y;
        u1y = u1(y1-u1y,y2-u2y);
        u2y = u2(y1-u1y,y2-u2y);
      } while (abs(u1y-u1p)>0.0001);
      return u1y;
    }
    public double u2y(double y1, double y2) {
      double u2p;
      double u1y = 0.0;
      double u2y = 0.0;
      do {
        u2p = u2y;
        u1y = u1(y1-u1y,y2-u2y);
        u2y = u2(y1-u1y,y2-u2y);
      } while (abs(u2y-u2p)>0.0001);
      return u2y;
    }
    public float[][] u1x() {
      float[][] u = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double x2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double x1 = i1;
          u[i2][i1] = (float)u1x(x1,x2);
        }
      }
      return u;
    }
    public float[][] u2x() {
      float[][] u = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double x2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double x1 = i1;
          u[i2][i1] = (float)u2x(x1,x2);
        }
      }
      return u;
    }
    public float[][] u1m() {
      float[][] u = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double m2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double m1 = i1;
          u[i2][i1] = (float)u1m(m1,m2);
        }
      }
      return u;
    }
    public float[][] u2m() {
      float[][] u = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double m2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double m1 = i1;
          u[i2][i1] = (float)u2m(m1,m2);
        }
      }
      return u;
    }
    public float[][] u1y() {
      float[][] u = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double y2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double y1 = i1;
          u[i2][i1] = (float)u1y(y1,y2);
        }
      }
      return u;
    }
    public float[][] u2y() {
      float[][] u = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double y2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double y1 = i1;
          u[i2][i1] = (float)u2y(y1,y2);
        }
      }
      return u;
    }
    public float[][] warp(float[][] f) {
      SincInterpolator si = new SincInterpolator();
      float[][] g = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double y2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double y1 = i1;
          double x1 = y1-u1y(y1,y2);
          double x2 = y2-u2y(y1,y2);
          g[i2][i1] = si.interpolate(_n1,1.0,0.0,_n2,1.0,0.0,f,x1,x2);
        }
      }
      return g;
    }
    public float[][] unwarp(float[][] g) {
      SincInterpolator si = new SincInterpolator();
      float[][] f = new float[_n2][_n1];
      for (int i2=0; i2<_n2; ++i2) {
        double x2 = i2;
        for (int i1=0; i1<_n1; ++i1) {
          double x1 = i1;
          double y1 = x1+u1x(x1,x2);
          double y2 = x2+u2x(x1,x2);
          f[i2][i1] = si.interpolate(_n1,1.0,0.0,_n2,1.0,0.0,g,y1,y2);
        }
      }
      return f;
    }
    protected Displacement(int n1, int n2) {
      _n1 = n1;
      _n2 = n2;
    }
    private int _n1,_n2;
  }

  /**
   * Constant (zero-strain) displacement.
   */
  private static class ConstantDisplacement extends Displacement {
    public ConstantDisplacement(double u1, double u2, int n1, int n2) {
      super(n1,n2);
      _u1 = u1;
      _u2 = u2;
    }
    public double u1(double x1, double x2) {
      return _u1;
    }
    public double u2(double x1, double x2) {
      return _u2;
    }
    private double _u1,_u2;
  }

  /**
   * Derivative-of-Gaussian displacement.
   */
  private static class GaussianDisplacement extends Displacement {
    public GaussianDisplacement(double u1max, double u2max, int n1, int n2) {
      super(n1,n2);
      _a1 = (n1-1)/2.0;
      _a2 = (n2-1)/2.0;
      _b1 = _a1/3.0;
      _b2 = _a2/3.0;
      _c1 = u1max*exp(0.5)/_b1;
      _c2 = u2max*exp(0.5)/_b2;
    }
    public double u1(double x1, double x2) {
      double xa1 = x1-_a1;
      double xa2 = x2-_a2;
      return -_c1*xa1*exp(-0.5*((xa1*xa1)/(_b1*_b1)+(xa2*xa2)/(_b2*_b2)));
    }
    public double u2(double x1, double x2) {
      double xa1 = x1-_a1;
      double xa2 = x2-_a2;
      return -_c2*xa2*exp(-0.5*((xa1*xa1)/(_b1*_b1)+(xa2*xa2)/(_b2*_b2)));
    }
    private double _a1,_a2;
    private double _b1,_b2;
    private double _c1,_c2;
  }

  /**
   * Sinusoid displacement.
   */
  private static class SinusoidDisplacement extends Displacement {
    public SinusoidDisplacement(double u1max, double u2max, int n1, int n2) {
      super(n1,n2);
      double l1 = n1-1;
      double l2 = n2-1;
      _a1 = u1max;
      _a2 = u2max;
      _b1 = 2.0*PI/l1;
      _b2 = 2.0*PI/l2;
    }
    public double u1(double x1, double x2) {
      return _a1*sin(_b1*x1)*sin(0.5*_b2*x2);
    }
    public double u2(double x1, double x2) {
      return _a2*sin(_b2*x2)*sin(0.5*_b1*x1);
    }
    private double _a1,_a2;
    private double _b1,_b2;
  }
}
