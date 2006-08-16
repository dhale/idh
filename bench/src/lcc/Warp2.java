package lcc;

import java.awt.*;
import javax.swing.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

import util.*;

public class Warp2 {

  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Warp2(args);
      }
    });
  }

  private float[][] _f,_g;
  private int _fontSize = 24;
  private int _width = 650;
  private int _height = 500;
  private String _pngDir = System.getProperty("png.dir");
  private String _dataDir = "/data";

  private Warp2(String[] args) {
    if (_pngDir==null)
      _pngDir = ".";
    initImages();
    doIntro();
  }

  private float[][] readImage() {
    int n1 = 315;
    int n2 = 315;
    String fileName = _dataDir+"/seis/vg/junks.dat";
    return Floats.readLittleEndian(fileName,n1,n2);
  }

  private void initImages() {
    _f = readImage();
    int n1 = _f[0].length;
    int n2 = _f.length;
    float[][][][] xye = computeWarpFunctions(n1,n2);
    float[][] x1 = xye[0][0];
    float[][] x2 = xye[0][1];
    float[][] y1 = xye[1][0];
    float[][] y2 = xye[1][1];
    float[][] e1 = xye[2][0];
    float[][] e2 = xye[2][1];
    _g = warp(x1,x2,_f);
  }

  private void doIntro() {
    float clip = 6.0f;
    plot(_f,clip,"imagef");
    plot(_g,clip,"imageg");
  }

  private void plot(float[][] f, float clip, String png) {
    PlotPanel panel = panel();
    PixelsView pv = panel.addPixels(f);
    if (clip!=0.0f) {
      pv.setClips(-clip,clip);
    } else {
      pv.setPercentiles(1.0f,99.0f);
    }
    frame(panel,png);
  }

  private PlotPanel panel() {
    PlotPanel panel = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
    panel.addColorBar();
    return panel;
  }

  private PlotFrame frame(PlotPanel panel, String png) {
    PlotFrame frame = new PlotFrame(panel);
    frame.setFontSize(_fontSize);
    frame.setSize(_width,_height);
    frame.setVisible(true);
    if (png!=null)
      frame.paintToPng(200,6,_pngDir+"/"+png+".png");
    return frame;
  }

  /**
   * Computes functions that define warping, unwarping and displacement.
   * Warping is q(y) = p(x(y)). Unwarping is p(x) = q(y(x)). The displacement
   * implied by warping is e = y(x)-x. (The displacement implied by unwarping
   * is x(y)-y, but is not computed.) All functions are uniformly-sampled.
   */
  private static float[][][][] computeWarpFunctions(int n1, int n2) {
    float a1 = 0.20f;
    float a2 = 0.10f;
    float b1 = (float)((n1-1)/2);
    float b2 = (float)((n2-1)/2);
    float s1 = 50.0f;
    float s2 = 50.0f;
    float[][][] x = warpGauss(a1,a2,b1,b2,s1,s2,n1,n2);
    float[][] x1 = x[0];
    float[][] x2 = x[1];
    float[][][] y = unwarpGauss(a1,a2,b1,b2,s1,s2,n1,n2);
    float[][] y1 = y[0];
    float[][] y2 = y[1];
    float[][] z1 = Array.rampfloat(0.0f,1.0f,0.0f,n1,n2);
    float[][] z2 = Array.rampfloat(0.0f,0.0f,1.0f,n1,n2);
    float[][] e1 = Array.sub(y1,z1);
    float[][] e2 = Array.sub(y2,z2);
    return new float[][][][]{{x1,x2},{y1,y2},{e1,e2}};
  }

  /**
   * Uses 2-D sinc interpolation to compute q(i1,i2) = p(x1(i1,i2),x2[i1,i2]).
   * Here, functions f(i1,i2) are specified and returned as arrays f[i2][i1].
   */
  private static float[][] warp(float[][] x1, float[][] x2, float[][] p) {
    int n1 = p[0].length;
    int n2 = p.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniform(n1,1.0,0.0,n2,1.0,0.0,p);
    float[][] q = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        q[i2][i1] = si.interpolate(x1[i2][i1],x2[i2][i1]);
      }
    }
    return q;
  }

  /**
   * Computes a uniformly-sampled warping function based on a gaussian. 
   * The warping function is x = (x1,x2), where x1 = x1(y1,y2) and 
   * x2 = x2(y1,y2) for sampled values of (y1,y2). The function x is defined 
   * so that x equals y when y equals b, x is less than y for y less than b, 
   * and x is greater than y for y greater than b. Thus the parameters b1 and 
   * b2 control the location of displacment. The parameters a control the 
   * maximum displacement, and the parameters s control the spatial extent of 
   * displacement.
   */
  private static float[][][] warpGauss(
    float a1, float a2, 
    float b1, float b2, 
    float s1, float s2, 
    int n1, int n2)
  {
    float[][][] x = new float[2][n2][n1];
    float[][] x1 = x[0];
    float[][] x2 = x[1];
    for (int i2=0; i2<n2; ++i2) {
      float y2 = (float)i2-b2;
      for (int i1=0; i1<n1; ++i1) {
        float y1 = (float)i1-b1;
        x1[i2][i1] = b1+y1*(1.0f+a1*gauss(s1,s2,y1,y2));
        x2[i2][i1] = b2+y2*(1.0f+a2*gauss(s1,s2,y1,y2));
      }
    }
    return x;
  }

  /**
   * Computes a uniformly-sampled unwarping function based on a gaussian. 
   * The unwarping function is y = (y1,y2), where y1 = y1(x1,x2) and 
   * y2 = y2(x1,x2). In other words, the unwarping function y(x) is the 
   * inverse of the warping function x(y). The inverse y(x) is computed 
   * by iteration, because x(y) and y(x) are transcendental functions.
   */
  private static float[][][] unwarpGauss(
    float a1, float a2, 
    float b1, float b2, 
    float s1, float s2, 
    int n1, int n2)
  {
    float[][][] y = new float[2][n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float x2 = (float)i2-b2;
      for (int i1=0; i1<n1; ++i1) {
        float x1 = (float)i1-b1;
        float y1 = x1;
        float y2 = x2;
        float y1p,y2p;
        do {
          y1p = y1;
          y2p = y2;
          y1 = x1/(1.0f+a1*gauss(s1,s2,y1p,y2p));
          y2 = x2/(1.0f+a2*gauss(s1,s2,y1p,y2p));
        } while (abs(y1-y1p)>0.001f || abs(y2-y2p)>0.001f);
        y[0][i2][i1] = y1+b1;
        y[1][i2][i1] = y2+b2;
      }
    }
    return y;
  }

  /**
   * 2-D gaussian function with specified widths (sigmas) s1 and s2.
   */
  private static float gauss(float s1, float s2, float x1, float x2) {
    float e1 = x1/s1;
    float e2 = x2/s2;
    return exp(-0.5f*(e1*e1+e2*e2));
  }
}
