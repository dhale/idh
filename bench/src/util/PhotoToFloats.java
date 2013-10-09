package util;

// Standard Java packages.
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*;
import javax.swing.*;

// Mines Java Toolkit (http://www.mines.edu/~dhale/jtk).
import edu.mines.jtk.awt.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Converts a color photograph to a 2D array of floats.
 */
public class PhotoToFloats {

  /**
   * Reads, processes, and plots one or more images.
   * @param args array of file names specified on the command line.
   */
  public static void main(String[] args) {

    // If no file names, print a helpful message and exit.
    if (args.length<2) {
      System.out.println("Must specify names of input and output files!");
      System.exit(1);
    }
    String fileInput = args[0];
    String fileOutput = args[1];

    // Read and plot red, green, and blue components.
    float[][][] rgb = rgbFromFile(fileInput);
    float[][] r = rgb[0];
    float[][] g = rgb[1];
    float[][] b = rgb[2];
    plot(transpose(r));
    plot(transpose(g));
    plot(transpose(b));

    // Convert to gray image, transpose, and scale.
    float[][] f = transpose(gray(r,g,b));
    mul(1.0f/255.0f,f,f);

    // Write to simple flat file of floats.
    try {
      System.out.println("n1="+f[0].length+" n2="+f.length);
      ArrayOutputStream aos = new ArrayOutputStream(fileOutput);
      aos.writeFloats(f);
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }

    // Plot the image.
    plot(f);
  }

  /**
   * Converts red, green, and blue components to grays. The transform used 
   * here (and in television) is gray = 0.30*red + 0.59*green + 0.11*blue.
   * @param r array of red components.
   * @param g array of green components.
   * @param b array of blue components.
   * @return output gray image.
   */
  public static float[][] gray(float[][] r, float[][] g, float[][] b) {
    //return copy(b); // TODO: fix this!
    int n1 = r[0].length;
    int n2 = r.length;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[i2][i1] = 0.30f*r[i2][i1]+0.59f*g[i2][i1]+0.11f*b[i2][i1];
      }
    }
    return f;
  }

  /**
   * Returns the transpose of the specified image. In the transpose,
   * rows become columns and columns become rows. In other words,
   * output g[i][j] = f[j][i], for all image sample indices (i,j).
   * @param f input image.
   * @return output image.
   */
  public static float[][] transpose(float[][] f) {
    //return copy(f); // TODO: fix this!
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i1][i2] = f[i2][i1];
      }
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private (You need not change anything below this line.)

  /**
   * Returns an array of 256 grays from black to white.
   */
  private static Color[] makeGrayColors() {
    Color[] colors = new Color[256];
    for (int i=0; i<256; ++i) {
      int r = i;
      int g = i;
      int b = i;
      colors[i] = new Color(r,g,b);
    }
    return colors;
  }

  /**
   * Plots a 2-D regular array as a gray-scale image.
   * @param f the array to plot.
   */
  private static void plot(float[][] f)  {
    plot(f,makeGrayColors());
  }

  /**
   * Plots a 2-D regular array using specified array of 256 colors.
   * @param f the array to plot.
   * @param colors array[256] of colors.
   */
  private static void plot(float[][] f, Color[] colors)  {
    if (colors.length!=256)
      throw new RuntimeException("Need exactly 256 colors!");

    // Width and height of frame to preserve aspect ratio, approximately.
    int n1 = f[0].length;
    int n2 = f.length;
    int size = 600;
    int w,h;
    if (n1<n2) {
      w = size;
      h = w*n1/n2;
    } else {
      h = size;
      w = h*n2/n1;
    }
    int widthAxis = 100;
    int heightAxis = 50;
    int widthColorBar = 70;
    w += widthAxis+widthColorBar;
    h += heightAxis;

    // Plot panel and frame using the package edu.mines.jtk.mosaic.
    PlotPanel panel = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = panel.addPixels(f);
    //pv.setClips(0.0f,255.0f);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setColorModel(ColorMap.makeIndexColorModel(colors));
    panel.addColorBar();
    panel.setColorBarWidthMinimum(widthColorBar);
    PlotFrame frame = new PlotFrame(panel);
    frame.setVisible(true);
    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    frame.setFontSize(24);
    frame.setSize(w,h);
  }

  /**
   * Returns RGB (red, green, blue) components of an image read from a file.
   * @param fileName the name of the file containing a JPEG image.
   * @return array of RGB components; the components (red,green,blue) are 
   * stored in {rgb[0],rgb[1],rgb[2]}.
   */
  private static float[][][] rgbFromFile(String fileName) {
    BufferedImage bi = null;
    try {
      File file = new File(fileName);
      bi = ImageIO.read(file);
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }
    int w = bi.getWidth();
    int h = bi.getHeight();
    int[] pixels = bi.getRGB(0,0,w,h,null,0,w);
    float[][][] rgb = new float[3][h][w];
    for (int i=0,k=0; i<h; ++i) {
      for (int j=0; j<w; ++j,++k) {
        int p = pixels[k];
        int r = (p>>16)&0xff;
        int g = (p>> 8)&0xff;
        int b = (p>> 0)&0xff;
        rgb[0][i][j] = (float)r;
        rgb[1][i][j] = (float)g;
        rgb[2][i][j] = (float)b;
      }
    }
    return rgb;
  }
}
