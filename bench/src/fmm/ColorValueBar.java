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
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;

/**
 * A color bar with an interface for choosing a color-mapped value.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.19
 */
public class ColorValueBar extends JPanel implements ColorMapListener {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a new color value chooser for the specified color map.
   * @param cmap the color map.
   */
  public ColorValueBar(ColorMap cmap) {
    this(cmap,null);
  }

  /**
   * Constructs a new color value chooser with specified label.
   * @param cmap the color map.
   * @param label the label; null, if none.
   */
  public ColorValueBar(ColorMap cmap, String label) {
    super();

    // Values from the specified color map.
    _cmap = cmap;
    double vmin = cmap.getMinValue();
    double vmax = cmap.getMaxValue();
    double vavg = 0.5*(vmin+vmax);

    // The mosaic contains one tile that contains a pixels view.
    // The pixels view is constructed later.
    _mosaic = new Mosaic(1,1,EnumSet.of(Mosaic.AxesPlacement.RIGHT));
    if (label!=null)
      _mosaic.getTileAxisRight(0).setLabel(label);
    _tile = _mosaic.getTile(0,0);

    // Color bar (not including the axis) will have width 25 pixels.
    int cbWidth = 25;
    _mosaic.setWidthMinimum(0,cbWidth);
    _mosaic.setWidthElastic(0,0);

    // Swatch and text field display the current color and value.
    _swatch = new JPanel();
    _swatch.setBorder(BorderFactory.createLineBorder(Color.black));
    _swatch.setPreferredSize(new Dimension(cbWidth+2,0));
    _swatch.setBackground(cmap.getColor(vavg));
    _vfield = new NumberTextField(vmin,vmax);
    _vfield.setValue(vavg);
    _vfield.setFormat("%1.3g");
    JPanel vpanel = new JPanel();
    vpanel.setLayout(new GridBagLayout());
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.gridx = 0;
    gbc.weightx = 0;
    gbc.insets.top = 3; // some fine tuning to align 
    gbc.insets.bottom = 3; // with the text field
    vpanel.add(_swatch,gbc);
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx = 1;
    gbc.weightx = 100;
    gbc.insets.top = 0;
    gbc.insets.bottom = 0;
    vpanel.add(_vfield,gbc);

    // Layout for everything.
    this.setLayout(new BorderLayout());
    this.add(_mosaic,BorderLayout.CENTER);
    this.add(vpanel,BorderLayout.SOUTH);

    // Listen for value entered in text field.
    _vfield.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        NumberTextField ntf = (NumberTextField)e.getSource();
        setValue(ntf.getDouble());
      }
    });

    // Listen for mouse events in color bar.
    _tile.addMouseListener(new ValueAdjuster());

    // Listen for changes to the color map.
    _cmap.addListener(this);
  }

  /**
   * Sets the currently chosen value.
   * @param v the value.
   */
  public void setValue(double v) {
    _vfield.setDouble(v);
    _swatch.setBackground(_cmap.getColor(v));
  }

  /**
   * Gets the currently chosen value.
   * @return the value.
   */
  public double getValue() {
    return _vfield.getDouble();
  }

  public void colorMapChanged(ColorMap cm) {
    float vmin = (float)cm.getMinValue();
    float vmax = (float)cm.getMaxValue();
    if (vmin==vmax) {
      vmin -= Math.ulp(vmin);
      vmax += Math.ulp(vmax);
    }
    setValueRange(vmin,vmax);
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
  // private

  private ColorMap _cmap;
  private Mosaic _mosaic;
  private Tile _tile;
  private PixelsView _pixels;
  private JPanel _swatch;
  private NumberTextField _vfield;

  private class ValueAdjuster extends MouseAdapter {
    public void mousePressed(MouseEvent e) {
      Tile tile = (Tile)e.getSource();
      _adjusting = true;
      setValue(getValue(e));
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
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        setValue(getValue(e));
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

  private void setValueRange(double vmin, double vmax) {
    _vfield.setValueRange(vmin,vmax);
    double v = getValue();
    if (v<vmin || v>vmax) {
      if (v<vmin) v = vmin;
      if (v>vmax) v = vmax;
      setValue(v);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    ColorMap cmap = new ColorMap(0.0,100.0,ColorMap.getJet());
    ColorValueBar cvb = new ColorValueBar(cmap,"Velocity (km/s)");
    JFrame frame = new JFrame();
    frame.setSize(100,500);
    frame.add(cvb);
    frame.add(new JTextField(),BorderLayout.NORTH);
    frame.setVisible(true);
  }
}
