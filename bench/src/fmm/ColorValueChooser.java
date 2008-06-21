/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.util.EnumSet;
import java.text.*;
import javax.swing.*;
import javax.swing.text.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * A color bar with associated values and an interface for choosing a value.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.19
 */
public class ColorValueChooser extends JPanel implements ColorMapListener {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a new color value chooser.
   */
  public ColorValueChooser() {
    this(null);
  }

  /**
   * Constructs a new color value chooser with specified label.
   * @param label the label; null, if none.
   */
  public ColorValueChooser(String label) {
    super();
    _mosaic = new Mosaic(1,1,EnumSet.of(Mosaic.AxesPlacement.RIGHT));
    if (label!=null)
      _mosaic.getTileAxisRight(0).setLabel(label);
    _mosaic.setWidthMinimum(0,25);
    _mosaic.setWidthElastic(0,0);
    _tile = _mosaic.getTile(0,0);
    JPanel vpanel = new JPanel();
    _swatch = new JPanel();
    _swatch.setBorder(BorderFactory.createLineBorder(Color.black));
    _swatch.setPreferredSize(new Dimension(27,0));
    _swatch.setBackground(Color.RED);
    _vfield = new NumberTextField(0,100,true);
    _vfield.setValue(50);
    vpanel.setLayout(new GridBagLayout());
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.gridx = 0;
    gbc.weightx = 0;
    gbc.insets.top = 3;
    gbc.insets.bottom = 3;
    vpanel.add(_swatch,gbc);
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx = 1;
    gbc.weightx = 100;
    gbc.insets.top = 0;
    gbc.insets.bottom = 0;
    vpanel.add(_vfield,gbc);
    this.setLayout(new BorderLayout());
    this.add(_mosaic,BorderLayout.CENTER);
    this.add(vpanel,BorderLayout.SOUTH);
    this.add(new JTextField(),BorderLayout.NORTH);
  }

  /**
   * Sets the range.
   * @param vb begin value of range.
   * @param ve end value of range.
   */
   /*
  public void setRange(double vb, double ve) {
    if (!_hasRange || _vb!=vb || _ve!=ve) {
      _vb = (float)vb;
      _ve = (float)ve;
      _hasRange = true;
      updatePointsView();
    }
  }

  public float getRangeBegin() {
    return _vb;
  }

  public float getRangeEnd() {
    return _ve;
  }

  public boolean hasRange() {
    return _hasRange;
  }
  */

  public void colorMapChanged(ColorMap cm) {
    float vmin = (float)cm.getMinValue();
    float vmax = (float)cm.getMaxValue();
    if (vmin==vmax) {
      vmin -= Math.ulp(vmin);
      vmax += Math.ulp(vmax);
    }
    int nv = 256;
    double dv = (vmax-vmin)/(nv-1);
    double fv = vmin;
    Sampling vs = new Sampling(nv,dv,fv);
    float[][] va = new float[nv][1];
    Color[] ca = new Color[nv];
    for (int iv=0; iv<nv; ++iv) {
      float vi = (float)vs.getValue(iv);
      va[iv][0] = vi;
      ca[iv] = cm.getColor(vi);
    }
    if (_pixels==null) {
      _pixels = new PixelsView(va);
      _pixels.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
      _pixels.setInterpolation(PixelsView.Interpolation.LINEAR);
      _tile.addTiledView(_pixels);
    }
    IndexColorModel icm = ColorMap.makeIndexColorModel(ca);
    _pixels.setClips(vmin,vmax);
    _pixels.setColorModel(icm);
    Sampling s1 = new Sampling(1);
    Sampling s2 = vs;
    _pixels.set(s1,s2,va);
  }

  // Override base class implementation.
  public void setFont(Font font) {
    super.setFont(font);
    if (_mosaic!=null)
      _mosaic.setFont(font);
    revalidate();
  }

  // Override base class implementation.
  public void setForeground(Color color) {
    super.setForeground(color);
    if (_mosaic!=null)
      _mosaic.setForeground(color);
  }

  // Override base class implementation.
  public void setBackground(Color color) {
    super.setBackground(color);
    if (_mosaic!=null)
      _mosaic.setBackground(color);
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
    //paintToRect((Graphics2D)g,0,0,getWidth(),getHeight());
  }

  protected Mosaic getMosaic() {
    return _mosaic;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Mosaic _mosaic;
  private Tile _tile;
  private PixelsView _pixels;
  private JPanel _swatch;
  private JFormattedTextField _vfield;

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    ColorValueChooser cvc = new ColorValueChooser("Velocity (km/s)");
    ColorMap cmap = new ColorMap(0.0,100.0,ColorMap.getJet());
    cmap.addListener(cvc);
    JFrame frame = new JFrame();
    frame.setSize(100,500);
    frame.add(cvc);
    frame.setVisible(true);
  }
}

/*
  ///////////////////////////////////////////////////////////////////////////
  // private

  boolean _hasRange = false; // true only after range has been set
  float _vb,_ve; // begin and end of range (begin may exceed end)
  private PointsView _points; // points view that shows the range

  private void updatePointsView() {
    if (!_hasRange)
      return;
    Mosaic mosaic = super.getMosaic();
    Tile tile = mosaic.getTile(0,0);
    tile.addMouseListener(new RangeAdjuster());
    Projector hp = tile.getHorizontalProjector();
    Projector vp = tile.getVerticalProjector();
    float h0 = (float)hp.v0();
    float h1 = (float)hp.v1();
    float hm = 0.5f*(h0+h1);
    float[][] x1 = new float[3][2];
    float[][] x2 = new float[3][2];
    x1[0][0] = h0;  x2[0][0] = _vb;
    x1[0][1] = h1;  x2[0][1] = _vb;
    x1[1][0] = h0;  x2[1][0] = _ve;
    x1[1][1] = h1;  x2[1][1] = _ve;
    x1[2][0] = hm;  x2[2][0] = _vb;
    x1[2][1] = hm;  x2[2][1] = _ve;
    if (_points==null) {
      _points = new PointsView(x1,x2);
      tile.addTiledView(_points);
    } else {
      _points.set(x1,x2);
    }
  }

  private class RangeAdjuster extends MouseAdapter {
    public void mousePressed(MouseEvent e) {
      Tile tile = (Tile)e.getSource();
      _adjusting = true;
      _vbegin = _vmouse = getValue(e);
      setRange(_vbegin,_vmouse);
      tile.addMouseMotionListener(_mml);
    }
    public void mouseReleased(MouseEvent e) {
      if (_adjusting) {
        Tile tile = (Tile)e.getSource();
        tile.removeMouseMotionListener(_mml);
        _adjusting = false;
      }
    }
    private boolean _adjusting;
    private float _vbegin,_vmouse;
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        _vmouse = getValue(e);
        setRange(_vbegin,_vmouse);
      }
    };
    private float getValue(MouseEvent e) {
      Tile tile = (Tile)e.getSource();
      Transcaler ts = tile.getTranscaler();
      Projector vp = tile.getVerticalProjector();
      double y = ts.y(e.getY());
      double v = vp.v(y);
      return (float)v;
    }
  }
*/
