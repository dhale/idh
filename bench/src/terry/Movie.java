/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package terry;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import static java.lang.Math.max;
import static java.lang.Math.min;
import java.util.Timer;
import java.util.TimerTask;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;

/**
 * Abstract base class for movies.
 * This class handles most of the graphical user interface for movies.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.01.05
 */
abstract class Movie {

  Movie(double dt) {
    _dt = dt;
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected abstract void initSampledFunction();
  protected abstract void stepSampledFunction(double t, double dt);
  protected abstract PlotPanel initPanel();
  protected abstract void updateView();

  protected void initFrame(final int width, final int height) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        _panel = initPanel();
        init();
        _frame = new PlotFrame(_panel);
        _frame.setFontSize(18);
        addButtons();
        _frame.setSize(width,height);
        _frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
        _frame.setVisible(true);
      }
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final double FPS_MIN = 1.0;
  private static final double FPS_MAX = 100.0;

  private double _dt;
  private double _time;
  private PlotFrame _frame;
  private PlotPanel _panel;
  private Timer _timer;
  private double _fps = 10.0;

  private void setFrameRate(double fps) {
    _fps = max(FPS_MIN,min(FPS_MAX,fps));
    if (_timer!=null) {
      stop();
      play();
    }
  }

  private void updateDisplay() {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        updateView();
        _panel.setTitle(String.format("time: %1.6G",_time));
      }
    });
  }

  private void init() {
    _time = 0.0;
    initSampledFunction();
    updateDisplay();
  }

  private void step() {
    stepSampledFunction(_time,_dt);
    _time += _dt;
    updateDisplay();
  }

  private void play() {
    if (_timer==null) {
      TimerTask task = new TimerTask() {
        public void run() {
          step();
        }
      };
      _timer = new Timer();
      int delay = max(1,(int)(1000.0/max(0.001,_fps))); // in ms
      _timer.schedule(task,0,delay);
    }
  }

  private void stop() {
    if (_timer!=null) {
      _timer.cancel();
      _timer = null;
    }
  }

  private void addButtons() {
    JButton saveToPngButton = new JButton(new AbstractAction("Save to PNG") {
      public void actionPerformed(ActionEvent event) {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.showSaveDialog(_frame);
        File file = fc.getSelectedFile();
        if (file!=null) {
          String filename = file.getAbsolutePath();
          _frame.paintToPng(300,6,filename);
        }
      }
    });
    JButton rewindButton = new JButton(
      new AbstractAction("",loadIcon("Rewind16.gif")) {
      public void actionPerformed(ActionEvent event) {
        init();
      }
    });
    JButton stepButton = new JButton(
      new AbstractAction("",loadIcon("StepForward16.gif")) {
      public void actionPerformed(ActionEvent event) {
        step();
      }
    });
    JButton playButton = new JButton(
      new AbstractAction("",loadIcon("Play16.gif")) {
      public void actionPerformed(ActionEvent event) {
        play();
      }
    });
    JButton stopButton = new JButton(
      new AbstractAction("",loadIcon("Stop16.gif")) {
      public void actionPerformed(ActionEvent event) {
        stop();
      }
    });
    JToolBar toolBar = new JToolBar();
    toolBar.add(rewindButton);
    toolBar.add(stepButton);
    toolBar.add(playButton);
    toolBar.add(stopButton);
    JSlider fpsSlider = new JSlider(1,100,(int)_fps);
    fpsSlider.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
        JSlider source = (JSlider)e.getSource();
        double fps = source.getValue();
        setFrameRate(fps);
      }
    });
    toolBar.add(new JSeparator(SwingConstants.VERTICAL));
    toolBar.add(fpsSlider);
    toolBar.add(new JSeparator(SwingConstants.VERTICAL));
    toolBar.add(saveToPngButton);
    _frame.add(toolBar,BorderLayout.NORTH);
  }

  private Icon loadIcon(String iconFile) {
    java.net.URL url = getClass().getResource("resources/"+iconFile);
    return (url!=null)?new ImageIcon(url):null;
  }
}
