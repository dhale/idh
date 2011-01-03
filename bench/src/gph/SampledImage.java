/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gph;

import java.awt.image.*;
import java.io.*;
import java.net.URL;
import java.nio.IntBuffer;
import javax.imageio.ImageIO;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/** 
 * An image with uniformly sampled x and y axes.
 * The x axis is horizontal and corresponds to the image width.
 * The y axis is vertical and corresponds to the image height.
 * The image pixel with indices (0,0) corresponds to the sample
 * with indices (0,0), the sample with minimum x and y coordinates.
 * Samplings for horizontal and vertical image dimensions associate
 * x and y coordinates with each image pixel. 
 * <p>
 * This class includes methods that facilitate rendering with either
 * Java2D or OpenGL. For example, the samplings associated with an
 * image can be used to compute texture coordinates for OpenGL.
 * <p>
 * The image is stored by rows, and the first row with index iy = 0
 * corresponds to the minimum y coordinate. This convention is fine
 * for OpenGL, but is inconsistent with the way pixels are usually 
 * stored in Java2D buffered images or in files, such as PNG files.
 * When converting to or from such formats, the image is flipped
 * automatically to account for this inconsistency.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2010.12.29
 */
public class SampledImage {

  /**
   * Returns a sampled image constructed from the specified image file.
   * @param sx sampling of x, the horizontal coordinate.
   * @param sy sampling of y, the vertical coordinate.
   * @param fileName name of image file (e.g., a PNG file).
   * @return the sampled image.
   */
  public static SampledImage fromFile(
    Sampling sx, Sampling sy, String fileName) 
  {
    BufferedImage bi = null;
    try {
      bi = ImageIO.read(new File(fileName));
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }
    return fromImage(sx,sy,bi);
  }

  /**
   * Returns a sampled image constructed from the specified resource file.
   * Looks for the specified file in a subdirectory named resources of 
   * the directory that contains the specified class.
   * @param sx sampling of x, the horizontal coordinate.
   * @param sy sampling of y, the vertical coordinate.
   * @param klass class with which to find resource file.
   * @param fileName name of resource file (e.g., a PNG file).
   * @return the sampled image.
   */
  public static SampledImage fromResource(
    Sampling sx, Sampling sy, Class klass, String fileName) 
  {
    BufferedImage bi = null;
    try {
      URL url = klass.getResource("resources/"+fileName);
      bi = ImageIO.read(url);
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }
    return fromImage(sx,sy,bi);
  }

  /**
   * Returns a sampled image constructed from the specified buffered image.
   * @param sx sampling of x, the horizontal coordinate.
   * @param sy sampling of y, the vertical coordinate.
   * @param bi the buffered image.
   * @return the sampled image.
   */
  public static SampledImage fromImage(
    Sampling sx, Sampling sy, BufferedImage bi) 
  {
    int w = bi.getWidth();
    int h = bi.getHeight();
    byte[][] r = new byte[h][w];
    byte[][] g = new byte[h][w];
    byte[][] b = new byte[h][w];
    byte[][] a = new byte[h][w];
    int[] pixels = bi.getRGB(0,0,w,h,null,0,w);
    for (int i=h-1,k=0; i>=0; --i) { // flip image vertically so that
      for (int j=0; j<w; ++j,++k) {  // first row <=> min yitude
        int p = pixels[k];
        a[i][j] = (byte)((p>>24)&0xff);
        r[i][j] = (byte)((p>>16)&0xff);
        g[i][j] = (byte)((p>> 8)&0xff);
        b[i][j] = (byte)((p>> 0)&0xff);
      }
    }
    return new SampledImage(sx,sy,r,g,b,a);
  }

  /**
   * Returns a sampled image constructed from the specified floats.
   * @param sx sampling of x, the horizontal coordinate.
   * @param sy sampling of y, the vertical coordinate.
   * @param f array of floats to be mapped to image pixels.
   * @param cm color map used to convert floats to colors.
   * @return the sampled image.
   */
  public static SampledImage fromFloats(
    Sampling sx, Sampling sy,
    float[][] f, ColorMap cm) 
  {
    IndexColorModel icm = cm.getColorModel();
    float fmin = (float)cm.getMinValue();
    float fmax = (float)cm.getMaxValue();
    float fscale = (fmax>fmin)?1.0f/(fmax-fmin):1.0f;
    float fshift = fmin;
    int m = f.length;
    int n = f[0].length;
    byte[][] rb = new byte[m][n];
    byte[][] gb = new byte[m][n];
    byte[][] bb = new byte[m][n];
    byte[][] ab = new byte[m][n];
    for (int i=0; i<m; ++i) {
      for (int j=0; j<n; ++j) {
        int kf = toInt((f[i][j]-fshift)*fscale);
        rb[i][j] = (byte)icm.getRed(kf);
        gb[i][j] = (byte)icm.getGreen(kf);
        bb[i][j] = (byte)icm.getBlue(kf);
        ab[i][j] = BYTE255;
      }
    }
    return new SampledImage(sx,sy,rb,gb,bb,ab);
  }

  /**
   * Constructs a sampled image with specified samplings and colors.
   * Alphas are assumed to equal 255, for no transparency.
   * @param sx sampling of x, the horizontal coordinate.
   * @param sy sampling of y, the vertical coordinate.
   * @param r array of reds; by reference, not by copy.
   * @param g array of greens; by reference, not by copy.
   * @param b array of blues; by reference, not by copy.
   */
  public SampledImage(
    Sampling sx, Sampling sy,
    byte[][] r, byte[][] g, byte[][] b) 
  {
    this(sx,sy,r,g,b,fillbyte(BYTE255,r[0].length,r.length));
  }

  /**
   * Constructs a sampled image with specified samplings and colors.
   * @param sx sampling of x, the horizontal coordinate.
   * @param sy sampling of y, the vertical coordinate.
   * @param r array of reds; by reference, not by copy.
   * @param g array of greens; by reference, not by copy.
   * @param b array of blues; by reference, not by copy.
   * @param a array of alphas; by reference, not by copy.
   */
  public SampledImage(
    Sampling sx, Sampling sy,
    byte[][] r, byte[][] g, byte[][] b, byte[][] a) 
  {
    Check.argument(r[0].length==sx.getCount(),
      "sampling sx is consistent with image width");
    Check.argument(r.length==sy.getCount(),
      "sampling sy is consistent with image height");
    _sx = sx;
    _sy = sy;
    _nx = sx.getCount();
    _ny = sy.getCount();
    for (_nx2=8; _nx2<_nx; _nx2*=2); // power of 2
    for (_ny2=8; _ny2<_ny; _ny2*=2); // power of 2
    _r = r;
    _g = g;
    _b = b;
    _a = a;
  }

  /**
   * Sets the alpha for all pixels in this image.
   * @param alpha the alpha value, in [0,1].
   */
  public void setAlpha(float alpha) {
    fill(toByte(alpha),_a);
  }

  /**
   * Sets the alphas for all pixels in this image.
   * @param a array of alpha values, in [0,1].
   */
  public void setAlpha(float[][] a) {
    for (int iy=0; iy<_ny; ++iy) {
      for (int ix=0; ix<_nx; ++ix) {
        _a[iy][ix] = toByte(a[iy][ix]);
      }
    }
  }

  /**
   * Gets the sampling of x, the horizontal coordinate.
   * @return the x sampling.
   */
  public Sampling getSamplingX() {
    return _sx;
  }

  /**
   * Gets the sampling of y, the vertical coordinate.
   * @return the y sampling.
   */
  public Sampling getSamplingY() {
    return _sy;
  }

  /**
   * Gets the byte array of reds in this image.
   * @return byte array of reds; by reference, not by copy.
   */
  public byte[][] getRed() {
    return _r;
  }

  /**
   * Gets the byte array of greens in this image.
   * @return byte array of greens; by reference, not by copy.
   */
  public byte[][] getGreen() {
    return _g;
  }

  /**
   * Gets the byte array of blues in this image.
   * @return byte array of blues; by reference, not by copy.
   */
  public byte[][] getBlue() {
    return _b;
  }

  /**
   * Gets the byte array of alphas in this image.
   * @return byte array of alphas; by reference, not by copy.
   */
  public byte[][] getAlpha() {
    return _a;
  }

  /**
   * Gets the image width.
   * @return the image width.
   */
  public int getWidth() {
    return _nx;
  }

  /**
   * Gets the image height.
   * @return the image height.
   */
  public int getHeight() {
    return _ny;
  }

  /**
   * Gets the image width padded to a power of 2.
   * This width is used for OpenGL textures.
   * @return the padded image width.
   */
  public int getWidth2() {
    return _nx2;
  }

  /**
   * Gets the image height padded to a power of 2.
   * This height is used for OpenGL textures.
   * @return the padded image height.
   */
  public int getHeight2() {
    return _ny2;
  }

  /**
   * Gets the texture coordinate s for the specified x coordinate.
   * Accounts for any power-of-2 padding of the texture width.
   * @param x the x coordinate.
   * @return texture coordinate s.
   */
  public float getTextureS(float x) {
    float fx = (float)_sx.getFirst();
    float dx = (float)_sx.getDelta();
    return (0.5f+(x-fx)/dx)/_nx2;
  }

  /**
   * Gets the texture coordinate t for the specified y coordinate.
   * Accounts for any power-of-2 padding of the texture height.
   * @param y the y coordinate.
   * @return texture coordinate t.
   */
  public float getTextureT(float y) {
    float fy = (float)_sy.getFirst();
    float dy = (float)_sy.getDelta();
    return (0.5f+(y-fy)/dy)/_ny2;
  }

  /**
   * Gets array {r,g,b,a} of arrays of image colors.
   * All returned color components are in the range [0,1].
   * @param flipVertical true, to flip image vertically.
   * @return array {r,g,b,a} of arrays of image colors.
   */
  public float[][][] getFloatsRGBA(boolean flipVertical) {
    float scale = 1.0f/255.0f;
    float[][] r = new float[_ny][_nx];
    float[][] g = new float[_ny][_nx];
    float[][] b = new float[_ny][_nx];
    float[][] a = new float[_ny][_nx];
    for (int iy=0; iy<_ny; ++iy) {
      int jy = flipVertical?_ny-1-iy:iy;
      for (int ix=0; ix<_nx; ++ix) {
        int ri = _r[jy][ix]; if (ri<0) ri += 256;
        int gi = _g[jy][ix]; if (gi<0) gi += 256;
        int bi = _b[jy][ix]; if (bi<0) bi += 256;
        int ai = _a[jy][ix]; if (ai<0) ai += 256;
        r[iy][ix] = ri*scale;
        g[iy][ix] = gi*scale;
        b[iy][ix] = bi*scale;
        a[iy][ix] = ai*scale;
      }
    }
    return new float[][][]{r,g,b,a};
  }

  /**
   * Gets array of image gray values, ignoring all alphas.
   * Returned gray values are in the range [0,1].
   * @param flipVertical true, to flip image vertically.
   * @return array {r,g,b,a} of arrays of image colors.
   */
  public float[][] getFloatsGray(boolean flipVertical) {
    float[][][] rgba = getFloatsRGBA(flipVertical);
    float[][] r = rgba[0], g = rgba[1], b = rgba[2];
    float[][] gray = new float[_ny][_nx];
    for (int iy=0; iy<_ny; ++iy)
      for (int ix=0; ix<_nx; ++ix)
        gray[iy][ix] = 0.30f*r[iy][ix]+0.59f*g[iy][ix]+0.11f*b[iy][ix];
    return gray;
  }

  /**
   * Gets this image as a buffered image, suitable for Java2D.
   * Image pixels are in the default Java2D ARGB format.
   * Image pixels are flipped, so that the sample with indices (0,ny-1) 
   * is the first pixel in the returned buffered image.
   * @return the buffered image.
   */
  public BufferedImage getBufferedImage() {
    BufferedImage bi = new BufferedImage(_nx,_ny,BufferedImage.TYPE_INT_ARGB);
    int[] pixels = new int[_nx*_ny];
    for (int iy=_ny-1,k=0; iy>=0; --iy) { // flip image vertically
      for (int ix=0; ix<_nx; ++ix,++k) {
        int a = _a[iy][ix];
        int r = _r[iy][ix];
        int g = _g[iy][ix];
        int b = _b[iy][ix];
        int p = ((a&0xff)<<24)|((r&0xFF)<<16)|((g&0xFF)<<8)|((b&0xff));
        pixels[k] = p;
      }
    }
    bi.setRGB(0,0,_nx,_ny,pixels,0,_nx);
    return bi;
  }

  /**
   * Gets this image as an int buffer, suitable for an OpenGL texture.
   * Image pixels are in the OpenGL RGBA format, and image dimensions
   * are padded to powers of 2. The sample with indices (0,0) is the
   * first pixel in the returned buffer, with texture coordinates
   * (0.5/nx2,0.5/ny2). The sample with indices (nx-1,ny-1) has texture
   * coordinates ((0.5+(nx-1))/nx2,(0.5+(ny-1))/ny2), where (nx,ny) and
   * (nx2,ny2) are image dimensions before and after padding. Any padded
   * pixels are computed by assuming the image is periodic.
   * @return the int buffer.
   */
  public IntBuffer getIntBuffer() {
    IntBuffer ib = Direct.newIntBuffer(_nx2*_ny2);
    for (int iy=0,k=0; iy<_ny2; ++iy) {
      int jy = iy%_ny;
      for (int ix=0; ix<_nx2; ++ix,++k) {
        int jx = ix%_nx;
        int a = _a[jy][jx];
        int r = _r[jy][jx];
        int g = _g[jy][jx];
        int b = _b[jy][jx];
        int p = ((a&0xff)<<24)|((b&0xff)<<16)|((g&0xff)<<8)|((r&0xff)<<0);
        ib.put(k,p);
      }
    }
    /*
    for (int iy=0; iy<_ny; ++iy) {
      for (int ix=0,k=ix+iy*_nx2; ix<_nx; ++ix,++k) {
        int a = _a[iy][ix];
        int r = _r[iy][ix];
        int g = _g[iy][ix];
        int b = _b[iy][ix];
        int p = ((a&0xff)<<24)|((b&0xff)<<16)|((g&0xff)<<8)|((r&0xff)<<0);
        ib.put(k,p);
      }
      for (int ix=_nx,k=ix+iy*_nx2; ix<_nx2; ++ix,++k)
        ib.put(k,0);
    }
    for (int iy=_ny; iy<_ny2; ++iy)
      for (int ix=0,k=iy*_nx2; ix<_nx2; ++ix,++k)
        ib.put(k,0);
    */
    return ib;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static byte BYTE255 = (byte)255;

  private Sampling _sx,_sy; // x and y sampling
  private int _nx,_ny; // width and height
  private int _nx2,_ny2; // width and height padded to power of 2
  private byte[][] _r,_g,_b,_a; // red, green, blue, alpha

  private static byte toByte(float v) {
    if (v<0.0f) v = 0.0f;
    if (v>1.0f) v = 1.0f;
    return (byte)(v*255.0f+0.5f);
  }

  private static int toInt(float v) {
    if (v<0.0f) v = 0.0f;
    if (v>1.0f) v = 1.0f;
    return (int)(v*255.0f+0.5f);
  }
}
