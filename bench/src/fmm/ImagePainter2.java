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
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * 2D image painting using fast marching methods.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.23
 */
public class ImagePainter2 {

  public ImagePainter2(float[][] image) {
    _n1 = image[0].length;
    _n2 = image.length;
    _nv = 1;
    _image = image;
    _painting = new Painting2(_n1,_n2,_nv);
    _painting.paintAt(0,0,0.5f);
    _painting.extrapolate();

    // Plot panel.
    int fontSize = 24;
    int width = 900;
    int height = 750;
    int widthColorBar = 80;
    PlotPanel.Orientation ppo = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel.AxesPlacement ppap = PlotPanel.AxesPlacement.LEFT_TOP;
    _panel = new PlotPanel(1,1,ppo,ppap);
    _panel.setColorBarWidthMinimum(widthColorBar);
    _colorBar = _panel.addColorBar();

    // Image view.
    _imageView = _panel.addPixels(_image);
    _imageView.setInterpolation(PixelsView.Interpolation.NEAREST);

    // Paint view on top of image view.
    _paintView = _panel.addPixels(_painting.getValues());
    _paintView.setInterpolation(PixelsView.Interpolation.NEAREST);
    _paintView.setColorModel(ColorMap.JET);

    // Color map.
    _colorMap = _paintView.getColorMap();
    System.out.println("_colorMap="+_colorMap);

    // Make paint control only after we have the color map and color bar.
    _paintControl = new PaintControl();
    _paintControl.setVisible(true);

    // Plot frame.
    _frame = new PlotFrame(_panel);
    _frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    _frame.setFontSize(fontSize);
    _frame.setSize(width,height);
    _frame.setVisible(true);
  }

  /**
   * Sets the range of values that can be painted.
   * @param vmin the minimum value.
   * @param vmax the maximum value.
   */
  public void setValueRange(double vmin, double vmax) {
    _paintView.setClips((float)vmin,(float)vmax);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2,_nv;
  private float[][] _image;
  private Painting2 _painting;

  private PlotPanel _panel;
  private PlotFrame _frame;
  private PixelsView _imageView;
  private PixelsView _paintView;
  private ColorMap _colorMap;
  private ColorBar _colorBar;
  private PaintControl _paintControl;

  private void updatePaintColorModel(float alpha) {
    IndexColorModel icm = _paintView.getColorModel();
    int n = 256;
    byte[] r = new byte[n];
    byte[] g = new byte[n];
    byte[] b = new byte[n];
    byte[] a = new byte[n];
    icm.getReds(r);
    icm.getGreens(g);
    icm.getBlues(b);
    byte ba = (byte)(255.0*alpha);
    a[0] = 0; // transparent for byte index 0
    for (int i=1; i<n; ++i)
      a[i] = ba;
    icm = new IndexColorModel(8,256,r,g,b,a);
    _paintView.setColorModel(icm);
  }

  private class PaintControl extends JFrame {
    public PaintControl() {
      JPanel panel = new JPanel();
      panel.setLayout(new BorderLayout());
      _cvc = new ColorValueChooser(_colorMap);
      _cvc.setColorBar(_colorBar);
      _ac = new AlphaChooser(this);
      setAlpha(0.5f);
      panel.add(_cvc,BorderLayout.NORTH);
      panel.add(_ac,BorderLayout.SOUTH);
      this.add(panel);
      this.pack();
    }
    public float getAlpha() {
      return _alpha;
    }
    void setAlpha(float alpha) {
      _alpha = alpha;
      updatePaintColorModel(alpha);
    }
    private ColorValueChooser _cvc;
    private AlphaChooser _ac;
    private float _alpha = 0.5f;
  }

  private class AlphaChooser extends JPanel {
    public AlphaChooser(PaintControl pc) {
      _pc = pc;
      _slider = new JSlider(0,100);
      _slider.setValue((int)(100.0f*_pc.getAlpha()));
      add(new JLabel("opacity:"));
      add(_slider);
      _slider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent e) {
          if (_slider==e.getSource()) {
            int value = _slider.getValue();
            _pc.setAlpha(0.01f*(float)value);
          }
        }
      });
    }
    private PaintControl _pc;
    private JSlider _slider;
  }

  private static void testImagePainter() {
    int n1 = 101;
    int n2 = 101;
    float[][] image = Array.randfloat(n1,n2);
    ImagePainter2 ip = new ImagePainter2(image);
    ip.setValueRange(0.0,1.0);
  }

  private static void trace(String s) {
    //System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testImagePainter();
      }
    });
  }
}
