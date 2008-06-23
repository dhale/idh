/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import java.util.*;
import javax.swing.*;

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

    int fontSize = 24;
    int width = 800;
    int height = 750;
    int widthColorBar = 80;
    PlotPanel.Orientation ppo = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel.AxesPlacement ppap = PlotPanel.AxesPlacement.LEFT_TOP;
    _panel = new PlotPanel(1,1,ppo,ppap);
    _panel.setColorBarWidthMinimum(widthColorBar);

    _imageView = _panel.addPixels(_image);
    _paintView = _panel.addPixels(_painting.getValues());
    _imageView.setInterpolation(PixelsView.Interpolation.NEAREST);
    _paintView.setInterpolation(PixelsView.Interpolation.NEAREST);
    _paintView.setColorModel(ColorMap.getJet(0.7f));

    _frame = new PlotFrame(_panel);
    _frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    _frame.setFontSize(fontSize);
    _frame.setSize(width,height);
    _frame.setVisible(true);
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

  private static void testImagePainter() {
    int n1 = 101;
    int n2 = 101;
    float[][] image = Array.randfloat(n1,n2);
    ImagePainter2 ip = new ImagePainter2(image);
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
