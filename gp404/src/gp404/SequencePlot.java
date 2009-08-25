/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gp404;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * A plot of one or more sequences.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.08.19
 */
public class SequencePlot {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a plot for one sequence.
   * @param l1 a sequence label.
   * @param s1 a sequence.
   */
  public SequencePlot(String l1, Sequence s1) {
    this(new String[]{l1},
         new Sequence[]{s1});
  }

  /**
   * Constructs a plot for two sequences.
   * @param l1 a sequence label.
   * @param s1 a sequence.
   * @param l2 a sequence label.
   * @param s2 a sequence.
   */
  public SequencePlot(String l1, Sequence s1,
                      String l2, Sequence s2) {
    this(new String[]{l1,l2},
         new Sequence[]{s1,s2});
  }

  /**
   * Constructs a plot for three sequences.
   * @param l1 a sequence label.
   * @param s1 a sequence.
   * @param l2 a sequence label.
   * @param s2 a sequence.
   * @param l3 a sequence label.
   * @param s3 a sequence.
   */
  public SequencePlot(String l1, Sequence s1,
                      String l2, Sequence s2,
                      String l3, Sequence s3) {
    this(new String[]{l1,l2,l3},
         new Sequence[]{s1,s2,s3});
  }

  /**
   * Constructs a plot for multiple sequences.
   * @param al array of sequence labels.
   * @param as array of sequences.
   */
  public SequencePlot(String[] al, Sequence[] as) {
    makeFrame(al,as);
  }

  /**
   * Sets the visibility of function value zero in this view.
   * The default visibility is SequenceView.Zero.ALWAYS.
   * @param zero the visibility of function value zero.
   */
  public void setZero(final SequenceView.Zero zero) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        for (int is=0; is<_sviews.length; ++is)
          if (_sviews[is]!=null)
            _sviews[is].setZero(zero);
      }
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private PlotFrame _frame;
  private PlotPanel _panel;
  private SequenceView[] _sviews;
  private PointsView[] _pviews;

  private void makeFrame(final String[] al, final Sequence[] as) {
    Check.argument(al.length==as.length,"al.length==as.length");
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        int ns = as.length;
        _sviews = new SequenceView[ns];
        _pviews = new PointsView[ns];
        _panel = new PlotPanel(ns,1);
        _panel.setHLabel("time");
        _panel.setHFormat("%1.6f");
        //_panel.setVLimits(-1.0,1.0);
        for (int is=0; is<ns; ++is) {
          Sequence s = as[is];
          Sampling st = s.getSampling();
          int nt = st.getCount();
          float[] values = s.getValues();
          String l = al[is];
          _panel.setVLabel(is,l);
          if (nt<=1001) {
            _sviews[is] = _panel.addSequence(is,0,st,values);
          } else {
            _pviews[is] = _panel.addPoints(is,0,st,values);
          }
        }
        _frame = new PlotFrame(_panel);
        addButtons();
        _frame.setFontSize(12);
        _frame.setSize(950,240*as.length);
        _frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
        _frame.setVisible(true);
      }
    });
  }

  private void addButtons() {
    Action saveToPngAction = new AbstractAction("Save to PNG") {
      private static final long serialVersionUID = 1L;
      public void actionPerformed(ActionEvent event) {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.showSaveDialog(_frame);
        File file = fc.getSelectedFile();
        if (file!=null) {
          String filename = file.getAbsolutePath();
          _frame.paintToPng(300,6,filename);
        }
      }
    };
    JButton saveToPngButton = new JButton(saveToPngAction);
    JToolBar toolBar = new JToolBar();
    toolBar.add(saveToPngButton);
    _frame.add(toolBar,BorderLayout.NORTH);
  }
}
