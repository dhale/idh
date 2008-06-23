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
import java.beans.*;
import java.util.EnumSet;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.StringUtil;

/**
 * An interface for choosing a color-mapped value.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.19
 */
public class ColorValueChooser extends JPanel {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a new color value chooser.
   * @param colorMap the color map.
   */
  public ColorValueChooser(ColorMap colorMap) {
    super();

    // Values from the specified color map.
    _colorMap = colorMap;
    double vmin = colorMap.getMinValue();
    double vmax = colorMap.getMaxValue();
    double vnow = 0.5*(vmin+vmax);

    // Swatch of color corresponding to the current value.
    _colorSwatch = new JPanel();
    _colorSwatch.setBorder(BorderFactory.createLineBorder(Color.black));
    _colorSwatch.setPreferredSize(new Dimension(25,0));

    // Number text field for current value.
    _vnowField = new NumberTextField(vmin,vmax);
    _vnowField.setFormat(_format);
    _vnowField.setValue(vnow);

    // Labels for min and max values.
    _vminLabel = new JLabel(formatValue(vmin));
    _vmaxLabel = new JLabel(formatValue(vmax));

    // Layout.
    this.setLayout(new GridBagLayout());
    GridBagConstraints gbc;
    gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.gridx = 0;
    gbc.gridy = 0;
    gbc.insets.left = 3; // a narrow margin
    gbc.weightx = 0;
    this.add(new JLabel("max:"),gbc);
    gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx = 1;
    gbc.gridy = 0;
    gbc.insets.left = 5; // align with text in number text field
    gbc.weightx = 100;
    this.add(_vmaxLabel,gbc);
    gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.gridx = 0;
    gbc.gridy = 1;
    gbc.weightx = 0;
    gbc.insets.top = 3; // align with number text field
    gbc.insets.bottom = 3; // align with number text field
    this.add(_colorSwatch,gbc);
    gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx = 1;
    gbc.gridy = 1;
    gbc.weightx = 100;
    gbc.insets.top = 0;
    gbc.insets.bottom = 0;
    this.add(_vnowField,gbc);
    gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.gridx = 0;
    gbc.gridy = 2;
    gbc.insets.left = 3; // a narrow margin
    gbc.weightx = 0;
    this.add(new JLabel("min:"),gbc);
    gbc = new GridBagConstraints();
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx = 1;
    gbc.gridy = 2;
    gbc.insets.left = 5; // align with text in number text field
    gbc.weightx = 100;
    this.add(_vminLabel,gbc);

    // Listen for changes to color map.
    _colorMap.addListener(new ColorMapListener() {
      public void colorMapChanged(ColorMap cm) {
        double vmin = cm.getMinValue();
        double vmax = cm.getMaxValue();
        _vminLabel.setText(formatValue(vmin));
        _vmaxLabel.setText(formatValue(vmax));
        if (getValue()<vmin)
          setValue(vmin);
        if (getValue()>vmax)
          setValue(vmax);
        updateColorSwatch();
      }
    });

    // Listen for changes to value in number text field.
    _vnowField.addPropertyChangeListener("value",new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
        if (_vnowField==e.getSource())
          updateColorSwatch();
      }
    });

    // Select all in the number text field to make changes easier.
    // Note that invokeLater is necessary to workaround some internal
    // handling of the selection by JFormattedTextField.
    _vnowField.addFocusListener(new FocusAdapter() {
      public void focusGained(FocusEvent e) {
        if (_vnowField==e.getSource()) {
          SwingUtilities.invokeLater(new Runnable() {
            public void run(){
              _vnowField.selectAll();
            }
          });
        }
      }
    });
    _vnowField.addMouseListener(new MouseAdapter() {
      public void mouseReleased(MouseEvent e) {
        if (_vnowField==e.getSource()) {
          SwingUtilities.invokeLater(new Runnable() {
            public void run(){
              _vnowField.selectAll();
            }
          });
        }
      }
    });
  }

  /**
   * Sets the color bar that may be used to choose a color.
   * @param colorBar the color bar; null, if none.
   */
  public void setColorBar(ColorBar colorBar) {
    if (_colorBar!=null) {
      _colorBar.getTile().removeMouseListener(_valueAdjuster);
    }
    _colorBar = colorBar;
    if (_colorBar!=null) {
      _valueAdjuster = new ValueAdjuster();
      _colorBar.getTile().addMouseListener(_valueAdjuster);
    }
  }

  /**
   * Gets the current chosen value.
   * @return the value.
   */
  public double getValue() {
    return _vnowField.getDouble();
  }

  /**
   * Sets the current chosen value. If necessary, the specified value
   * will be adjusted to be within current minimum and maximum values.
   * @param v the value.
   */
  public void setValue(double v) {
    _vnowField.setDouble(v);
    updateColorSwatch();
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final String _format = "%1.3g";

  private ColorMap _colorMap;
  private ColorBar _colorBar;
  private ValueAdjuster _valueAdjuster;
  private JPanel _colorSwatch;
  private JLabel _vminLabel;
  private JLabel _vmaxLabel;
  private NumberTextField _vnowField;

  private static String formatValue(double v) {
    String s = String.format(_format,v);
    s = StringUtil.removeTrailingZeros(s);
    return s;
  }

  private void updateColorSwatch() {
    _colorSwatch.setBackground(_colorMap.getColor(getValue()));
  }

  // Listens for mouse events in tile of color bar.
  private class ValueAdjuster extends MouseAdapter {
    public void mousePressed(MouseEvent e) {
      _adjusting = true;
      setValue(getValue(e));
      Tile tile = (Tile)e.getSource();
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
      return (float)tile.pixelToWorldVertical(e.getY());
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    ColorMap cm = new ColorMap(0.0,100.0,ColorMap.getJet());
    ColorBar cb = new ColorBar("Velocity (km/s)");
    cm.addListener(cb);
    ColorValueChooser cvc = new ColorValueChooser(cm);
    cvc.setColorBar(cb);
    JFrame frame = new JFrame();
    frame.setSize(200,500);
    frame.add(cvc);
    frame.add(cb,BorderLayout.EAST);
    frame.add(new JTextField(),BorderLayout.NORTH);
    frame.setVisible(true);
  }
}
