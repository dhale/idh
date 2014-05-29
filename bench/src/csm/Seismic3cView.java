/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package edu.mines.jtk.mosaic;

import java.awt.*;
import java.awt.image.*;
import java.util.ArrayList;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.Check;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.05.20
 */
public class Seismic3cView extends TiledView implements ColorMapped {

  public enum Orientation {
    X1RIGHT_X2UP,
    X1DOWN_X2RIGHT
  }

  public Seismic3cView(float[] x1, float[] x2, float[][][] f) {
    set(x1,x2,f);
  }

  public void set(float[] x1, float[] x2, float[][][] f) {
    _nx = x1.length;
    _ns = f[0][0].length;
    _x1 = copy(x1);
    _x2 = copy(x2);
    _f = copy(f);
    updateBestProjectors();
    repaint();
  }
  
  public void setSample(int ks) {
    _ks = ks;
    repaint();
  }

  /**
   * Sets the orientation of (x1,x2) axes.
   * @param orientation the orientation.
   */
  public void setOrientation(Orientation orientation) {
    if (_orientation!=orientation) {
      _orientation = orientation;
      updateBestProjectors();
      repaint();
    }
  } 

  public void setColorModel(IndexColorModel colorModel) {
    _colorMap.setColorModel(colorModel);
    repaint();
  }

  public ColorMap getColorMap() {
    return _colorMap;
  }

  public void paint(Graphics2D g2d) {
    g2d.setRenderingHint(
      RenderingHints.KEY_ANTIALIASING,
      RenderingHints.VALUE_ANTIALIAS_ON);

    Projector hp = getHorizontalProjector();
    Projector vp = getVerticalProjector();
    Transcaler ts = getTranscaler();

    // Font size and line width from graphics context.
    float fontSize = g2d.getFont().getSize2D();
    float lineWidth = 1.0f;
    Stroke stroke = g2d.getStroke();
    if (stroke instanceof BasicStroke) {
      BasicStroke bs = (BasicStroke)stroke;
      lineWidth = bs.getLineWidth();
    }

    // Graphics context for marks.
    int markSize = round(fontSize/2.0f);
    Graphics2D gmark = (Graphics2D)g2d.create();
    if (_markSize>=0.0f)
      markSize = round(_markSize*lineWidth);
    float width = lineWidth;
    if (_lineWidth!=0.0f)
      width *= _lineWidth;
    BasicStroke bs = new BasicStroke(width);
    gmark.setStroke(bs);

    // Arrays for (x,y) coordinates.
    int[] x = new int[_nx];
    int[] y = new int[_nx];

    // Compute (x,y) coordinates.
    computeXY(hp,vp,ts,_x1,_x2,_f,_ks,x,y);

    // Draw marks at points.
    paintFilledCircles(gmark,markSize,x,y,_f[2]);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  int _ns; // number of samples
  int _nx; // number of stations with (x1,x2) coordinates
  float[] _x1; // x1 coordinates of stations
  float[] _x2; // x2 coordinates of stations
  float[][][] _f; // arrays {vert,inline,crossline} of 3C amplitudes
  int _ks; // index of displayed samples
  private Orientation _orientation = Orientation.X1RIGHT_X2UP;
  private float _markSize = 12.0f;
  private float _lineWidth = 1.0f;
  private ColorMap _colorMap = new ColorMap(ColorMap.BLUE_WHITE_RED);

  /**
   * Called when we might new realignment.
   */
  private void updateBestProjectors() {

    // Min and max (x1,x2) values.
    float x1min =  FLT_MAX;
    float x2min =  FLT_MAX;
    float x1max = -FLT_MAX;
    float x2max = -FLT_MAX;
    float dxsum = 0.0f;
    for (int ix=0; ix<_nx; ++ix) {
      float x1i = _x1[ix];
      float x2i = _x2[ix];
      if (ix>0) {
        float dx1 = x1i-_x1[ix-1];
        float dx2 = x2i-_x2[ix-1];
        dxsum += sqrt(dx1*dx1+dx2*dx2);
      }
      x1min = min(x1min,x1i);
      x2min = min(x2min,x2i);
      x1max = max(x1max,x1i);
      x2max = max(x2max,x2i);
    }
    float dxavg = dxsum/(_nx-1);
    x1min -= dxavg;
    x2min -= dxavg;
    x1max += dxavg;
    x2max += dxavg;

    // Ensure x1min<x1max and x2min<x2max.
    if (x1min==x1max) {
      x1min -= ulp(x1min);
      x1max += ulp(x1max);
    }
    if (x2min==x2max) {
      x2min -= ulp(x2min);
      x2max += ulp(x2max);
    }

    // Assume mark sizes and line widths less than 2% of plot dimensions.
    // The goal is to avoid clipping big marks and wide lines. The problem
    // is that mark sizes and line widths are specified in screen pixels
    // (or points), but margins u0 and u1 are specified in normalized 
    // coordinates, fractions of our tile's width and height. Here, we do 
    // not know those dimensions.
    double u0 = 0.01;
    double u1 = 0.99;

    // Best projectors.
    Projector bhp = null;
    Projector bvp = null;
    if (_orientation==Orientation.X1RIGHT_X2UP) {
      bhp = (x1min<x1max)?new Projector(x1min,x1max,u0,u1):null;
      bvp = (x2min<x2max)?new Projector(x2max,x2min,u0,u1):null;
    } else if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      bhp = (x2min<x2max)?new Projector(x2min,x2max,u0,u1):null;
      bvp = (x1min<x1max)?new Projector(x1min,x1max,u0,u1):null;
    }
    setBestProjectors(bhp,bvp);
  }

  private void computeXY(
    Projector hp, Projector vp, Transcaler ts,
    float[] x1, float[] x2, float[][][] f, int k, int[] x, int[] y) 
  {
    ts = ts.combineWith(hp,vp);
    float[] xv = null;
    float[] yv = null;
    float[] f1 = new float[_nx];
    float[] f2 = new float[_nx];
    for (int ix=0; ix<_nx; ++ix) {
      f1[ix] = f[1][ix][_ks];
      f2[ix] = f[0][ix][_ks];
    }
    float[] fx,fy;
    if (_orientation==Orientation.X1RIGHT_X2UP) {
      xv = x1;
      yv = x2;
      fx = f1;
      fy = f2;
    } else {
      xv = x2;
      yv = x1;
      fx = f2;
      fy = f1;
    }
    for (int i=0; i<_nx; ++i) {
      x[i] = ts.x(xv[i]+fx[i]);
      y[i] = ts.y(yv[i]+fy[i]);
    }
  }
  
  private void paintFilledCircles(
    Graphics2D g2d, int s, int[] x, int[] y, float[][] f) 
  {
    Graphics2D gmark = (Graphics2D)g2d.create();
    IndexColorModel cm = _colorMap.getColorModel();
    float clipMin = -4.0f;
    float clipMax =  4.0f;
    int n = x.length;
    int wh = 1+2*(s/2);
    int xy = wh/2;
    for (int i=0; i<n; ++i) {
      float fi = f[i][_ks];
      int index;
      if (fi<clipMin) {
        index = 0;
      } else if (fi>clipMax) {
        index = 255;
      } else {
        index = (int)((fi-clipMin)/(clipMax-clipMin)*255f);
      }
      gmark.setColor(new Color(cm.getRGB(index)));
      gmark.fillOval(x[i]-xy,y[i]-xy,wh,wh);
      gmark.setColor(Color.BLACK);
      gmark.drawOval(x[i]-xy,y[i]-xy,wh,wh);
    }
  }
}
