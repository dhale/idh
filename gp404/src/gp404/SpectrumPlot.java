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
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A plot of a sequence with its amplitude (and optional phase) spectrum.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.10.29
 */
public class SpectrumPlot {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs an amplitude spectrum plot for the specified sequence.
   * @param label sequence label.
   * @param sequence a sequence.
   */
  public SpectrumPlot(String label, Sequence sequence) {
    this(label,sequence,false,false);
  }

  /**
   * Constructs an amplitude spectrum plot for the specified sequence.
   * @param label sequence label.
   * @param sequence a sequence.
   * @param db true, for amplitude in decibels (dB); false, otherwise.
   */
  public SpectrumPlot(String label, Sequence sequence, boolean db) {
    this(label,sequence,db,false);
  }

  /**
   * Constructs a spectrum plot for the specified sequence.
   * @param label sequence label.
   * @param sequence a sequence.
   * @param db true, for amplitude in decibels (dB); false, otherwise.
   * @param plotPhase true, for phase spectrum; false, for amplitude only.
   */
  public SpectrumPlot(
    String label, Sequence sequence, boolean db, boolean plotPhase) 
  {
    init(label,sequence,db,plotPhase);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean MENU = false;

  private PlotPanel _plotPanelS;
  private PlotPanel _plotPanelAP;
  private PlotFrame _plotFrame;

  // Computes the amplitude and phase spectra for the specified sequence.
  private static Sequence[] computeAmplitudeAndPhase(
    Sequence sx, boolean db) 
  {

    // Time sampling.
    Sampling ts = sx.getSampling();
    int nt = ts.getCount();
    double dt = ts.getDelta();
    double ft = ts.getFirst();

    // Frequency sampling.
    int nfft = FftReal.nfftSmall(2*nt);
    int nf = nfft/2+1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling fs = new Sampling(nf,df,ff);

    // Real-to-complex fast Fourier transform.
    float[] x = sx.getValues();
    FftReal fft = new FftReal(nfft);
    float[] cf = new float[2*nf];
    copy(nt,x,cf);
    fft.realToComplex(-1,cf,cf);

    // Adjust phase for possibly non-zero time of first sample.
    float[] wft = rampfloat(0.0f,-2.0f*FLT_PI*(float)(df*ft),nf);
    cf = cmul(cf,cmplx(cos(wft),sin(wft)));

    // Amplitude spectrum, normalized.
    float[] af = cabs(cf);
    float amax = max(max(af),FLT_EPSILON);
    af = mul(1.0f/amax,af);
    if (db) {
      af = log10(af);
      af = mul(20.0f,af);
    }
    Sequence a = new Sequence(fs,af);

    // Phase spectrum, in cycles.
    float[] pf = carg(cf);
    pf = mul(0.5f/FLT_PI,pf);
    Sequence p = new Sequence(fs,pf);

    return new Sequence[]{a,p};
  }

  private void init(
    final String label, final Sequence sequence, 
    final boolean db, final boolean plotPhase) 
  {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        makePlot(label,sequence,db,plotPhase);
      }
    });
  }

  private void makePlot(
    String label, Sequence sequence, boolean db, boolean plotPhase) 
  {
    // Amplitude (and optional phase) spectrum.
    Sequence[] ap = computeAmplitudeAndPhase(sequence,db);
    Sequence a = ap[0];
    Sequence p = ap[1];

    // One plot panel for the sequence.
    _plotPanelS = new PlotPanel();
    _plotPanelS.setHLabel("time");
    _plotPanelS.setHFormat("%1.6f");
    _plotPanelS.setVLabel(label);

    // Another plot panel for the amplitude and phase responses.
    // Amplitude and phase are both functions of frequency, so they
    // can share a horizontal axis.
    // The amplitude response is in tile (0,0) (on top).
    // The phase response is in tile (0,1) (on bottom).
    if (plotPhase) {
      _plotPanelAP = new PlotPanel(2,1);
    } else {
      _plotPanelAP = new PlotPanel(1,1);
    }
    if (db) {
      _plotPanelAP.setVLimits(0,-100.0,0.0);
      _plotPanelAP.setVLabel(0,"amplitude (dB)");
    } else {
      _plotPanelAP.setVLimits(0,0.0,1.0);
      _plotPanelAP.setVLabel(0,"amplitude");
    }
    if (plotPhase) {
      _plotPanelAP.setVLimits(1,-0.5,0.5);
      _plotPanelAP.setVLabel(1,"phase (cycles)");
    }
    _plotPanelAP.setHLabel("frequency");

    // Plots in the panel in the frame.
    SequenceView sv = 
      _plotPanelS.addSequence(sequence.getSampling(),sequence.getValues());
    sv.setZero(SequenceView.Zero.NORMAL);
    _plotPanelAP.addPoints(0,0,a.getSampling(),a.getValues());
    if (plotPhase)
      _plotPanelAP.addPoints(1,0,p.getSampling(),p.getValues());

    // Both plot panels share a common frame, with the sequence on top
    // and the spectra on the bottom.
    _plotFrame = new PlotFrame(
      _plotPanelS,_plotPanelAP,PlotFrame.Split.VERTICAL);

    // The menu bar or simple buttons.
    if (MENU) {
      JMenu fileMenu = new JMenu("File");
      fileMenu.setMnemonic('F');
      fileMenu.add(new SaveAsPngAction(_plotFrame)).setMnemonic('a');
      fileMenu.add(new ExitAction()).setMnemonic('x');
      JMenuBar menuBar = new JMenuBar();
      menuBar.add(fileMenu);
      _plotFrame.setJMenuBar(menuBar);
    } else {
      JButton saveToPngButton = new JButton(new SaveAsPngAction(_plotFrame));
      JToolBar toolBar = new JToolBar();
      toolBar.add(saveToPngButton);
      _plotFrame.add(toolBar,BorderLayout.NORTH);
    }

    // Make the plot frame visible.
    _plotFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    _plotFrame.setLocation(200,50);
    _plotFrame.setSize(800,700);
    _plotFrame.setFontSize(12);
    _plotFrame.setVisible(true);
  }

  private class ExitAction extends AbstractAction {
    private ExitAction() {
      super("Exit");
    }
    public void actionPerformed(ActionEvent event) {
      System.exit(0);
    }
  }
  private class SaveAsPngAction extends AbstractAction {
    private PlotFrame _plotFrame;
    private SaveAsPngAction(PlotFrame plotFrame) {
      super("Save as PNG");
      _plotFrame = plotFrame;
    }
    public void actionPerformed(ActionEvent event) {
      JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
      fc.showSaveDialog(_plotFrame);
      File file = fc.getSelectedFile();
      if (file!=null) {
        String filename = file.getAbsolutePath();
        _plotFrame.paintToPng(300,6,filename);
      }
    }
  }

  private void addButtons() {
    Action saveToPngAction = new AbstractAction("Save to PNG") {
      private static final long serialVersionUID = 1L;
      public void actionPerformed(ActionEvent event) {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.showSaveDialog(_plotFrame);
        File file = fc.getSelectedFile();
        if (file!=null) {
          String filename = file.getAbsolutePath();
          _plotFrame.paintToPng(300,6,filename);
        }
      }
    };
    JButton saveToPngButton = new JButton(saveToPngAction);
    JToolBar toolBar = new JToolBar();
    toolBar.add(saveToPngButton);
    _plotFrame.add(toolBar,BorderLayout.NORTH);
  }
}
