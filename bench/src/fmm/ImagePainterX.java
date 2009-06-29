/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.IndexColorModel;
import java.io.File;
import java.io.IOException;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.EigenTensors2;
import edu.mines.jtk.dsp.LocalOrientFilter;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.ArrayMath;
import static edu.mines.jtk.util.MathPlus.*;
import edu.mines.jtk.util.Quantiler;

/**
 * 2D image painting using fast marching methods.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.23
 */
public class ImagePainterX {

  public ImagePainterX(float[][] image) {
    _n1 = image[0].length;
    _n2 = image.length;
    _nv = 1;
    _image = image;
    _st = new StructureTensors(SIGMA,2.0f,1.0f,1.0f,_image);
    _painting = new PaintingX(_n1,_n2,_nv,_st);
    _painting.setDefaultValue(0.0f);

    int fontSize = 24;
    int width = 1330;
    int height = 860;
    int widthColorBar = 80;

    // Plot panel.
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

    // Tensors view, if visible, on top of paint view.
    float[][][] x12 = getTensorEllipses(_n1,_n2,10,_st);
    float[][] x1 = x12[0];
    float[][] x2 = x12[1];
    _tensorsView = new PointsView(x1,x2);
    _tensorsView.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    _tensorsView.setLineColor(Color.YELLOW);

    // Color map.
    _colorMap = _paintView.getColorMap();

    // Plot frame.
    _frame = new PlotFrame(_panel);
    _frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    _frame.setFontSize(fontSize);
    _frame.setSize(width,height);
    _frame.setVisible(true);

    // Make paint control only after we have the color map and color bar.
    _paintControl = new PaintControl(_frame);

    makeModesMenusAndToolBar();
  }

  /**
   * Sets the range of values that can be painted.
   * @param vmin the minimum value.
   * @param vmax the maximum value.
   */
  public void setValueRange(double vmin, double vmax) {
    _valueMin = (float)vmin;
    _valueMax = (float)vmax;
    _paintView.setClips(_valueMin,_valueMax);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static float SIGMA = 3.0f;

  private int _n1,_n2,_nv;
  private float[][] _image;
  private float _valueMin,_valueMax;
  private PaintingX _painting;
  private StructureTensors _st;

  private PlotPanel _panel;
  private PlotFrame _frame;
  private PixelsView _imageView;
  private PixelsView _paintView;
  private PointsView _tensorsView;
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
    for (int i=0; i<n; ++i)
      a[i] = ba;
    icm = new IndexColorModel(8,n,r,g,b,a);
    _paintView.setColorModel(icm);
  }

  // Paint controls are visible when paint mode is active.
  private class PaintControl extends JDialog {
    public PaintControl(JFrame parent) {
      super(parent,false); // not modal
      setAlwaysOnTop(true);
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
      this.setLocation(new Point(_frame.getWidth(),0));
    }
    public float getValue() {
      return (float)_cvc.getValue();
    }
    public void setValue(float value) {
      _cvc.setValue(value);
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
      _slider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent e) {
          if (_slider==e.getSource()) {
            int value = _slider.getValue();
            _pc.setAlpha(0.01f*(float)value);
          }
        }
      });
      Box box = Box.createHorizontalBox();
      box.add(new JLabel("opacity:"));
      box.add(_slider);
      this.add(box);
    }
    private PaintControl _pc;
    private JSlider _slider;
  }

  private class PaintMode extends Mode {
    public PaintMode(ModeManager modeManager) {
      super(modeManager);
      setName("Paint");
      setIcon(loadIcon(PaintMode.class,"resources/PaintIcon16.png"));
      setMnemonicKey(KeyEvent.VK_P);
      setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_P,0));
      setShortDescription("Paint samples");
    }
    protected void setActive(Component component, boolean active) {
      if (component instanceof Tile) {
        if (active) {
          component.addMouseListener(_ml);
          _paintControl.setVisible(true);
        } else {
          component.removeMouseListener(_ml);
          _paintControl.setVisible(false);
        }
      }
    }
    private boolean _down; // true, if mouse is down (painting or erasing)
    private boolean _erasing; // true, if erasing instead of painting
    private float _value; // the paint value
    private int _i1Paint,_i2Paint; // indices of last sample painted
    private Tile _tile; // tile in which painting began
    private MouseListener _ml = new MouseAdapter() {;
      public void mousePressed(MouseEvent e) {
        _erasing = e.isControlDown() || e.isAltDown();
        if (beginPaint(e)) {
          _down = true;
          _tile.addMouseMotionListener(_mml);
        }
      }
      public void mouseReleased(MouseEvent e) {
        if (_down) {
          _tile.removeMouseMotionListener(_mml);
          endPaint(e);
          _down = false;
        }
      }
    };
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        if (_down)
          duringPaint(e);
      }
    };
    private int getIndex1(MouseEvent e) {
      _tile = (Tile)e.getSource();
      double x1 = _tile.pixelToWorldVertical(e.getY());
      int i1 = (int)(x1+0.5);
      return (0<=i1 && i1<_n1)?i1:-1;
    }
    private int getIndex2(MouseEvent e) {
      _tile = (Tile)e.getSource();
      double x2 = _tile.pixelToWorldHorizontal(e.getX());
      int i2 = (int)(x2+0.5);
      return (0<=i2 && i2<_n2)?i2:-1;
    }
    private boolean beginPaint(MouseEvent e) {
      int i1 = getIndex1(e);
      int i2 = getIndex2(e);
      if (e.isShiftDown()) {
        _paintControl.setValue(valueAt(i1,i2));
        return false;
      } else {
        _value = _paintControl.getValue();
        return paintAt(i1,i2);
      }
    }
    private void duringPaint(MouseEvent e) {
      int i1 = getIndex1(e);
      int i2 = getIndex2(e);
      paintAt(i1,i2);
    }
    private void endPaint(MouseEvent e) {
      duringPaint(e);
      _i1Paint = -1;
      _i2Paint = -1;
    }
    private float valueAt(int i1, int i2) {
      if (0<=i1 && i1<_n1 &&
          0<=i2 && i2<_n2) {
        return _painting.getValue(i1,i2);
      }
      return _paintControl.getValue();
    }
    private boolean paintAt(int i1, int i2) {
      if ((i1!=_i1Paint || i2!=_i2Paint) &&
          0<=i1 && i1<_n1 &&
          0<=i2 && i2<_n2) {
        _i1Paint = i1;
        _i2Paint = i2;
        float vi = _paintControl.getValue();
        if (_erasing) { // eraser size is 3x3 pixels
          for (int j2=i2-1; j2<=i2+1; ++j2) {
            for (int j1=i1-1; j1<=i1+1; ++j1) {
               if (0<i1 && i1<_n1-1 &&
                   0<i2 && i2<_n2-1) 
                 _painting.eraseFixedAt(j1,j2);
            }
          }
        } else {
          _painting.paintAt(i1,i2,_value);
        }
        _paintView.set(_painting.getValues());
        return true;
      }
      return false;
    }
  }

  private void makeModesMenusAndToolBar() {

    // Modes.
    ModeManager mm = _frame.getModeManager();
    TileZoomMode tzm = _frame.getTileZoomMode();
    PaintMode pm = new PaintMode(mm);

    // Menus.
    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic('F');
    fileMenu.add(new SaveAsPngAction(_frame)).setMnemonic('a');
    fileMenu.add(new ExitAction()).setMnemonic('x');
    JMenu modeMenu = new JMenu("Mode");
    modeMenu.setMnemonic('M');
    modeMenu.add(new ModeMenuItem(tzm));
    modeMenu.add(new ModeMenuItem(pm));
    JMenu viewMenu = new JMenu("View");
    viewMenu.add(new JCheckBoxMenuItem(new ShowTensorsAction()));
    JMenu structureMenu = new JMenu("Structure");
    JMenuItem isotropicItem = new JRadioButtonMenuItem(
      new AbstractAction("Isotropic") {
        public void actionPerformed(ActionEvent e) {
          updateStructureTensors(0.0f,0.0f,0.0f);
        }
      });
    JMenuItem linearItem = new JRadioButtonMenuItem(
      new AbstractAction("Linear") {
        public void actionPerformed(ActionEvent e) {
          updateStructureTensors(0.0f,100.0f,1.0f);
        }
      });
    JMenuItem layersItem = new JRadioButtonMenuItem(
      new AbstractAction("Layers") {
        public void actionPerformed(ActionEvent e) {
          updateStructureTensors(1.0f,1.0f,1.0f);
        }
      });
    JMenuItem interfacesItem = new JRadioButtonMenuItem(
      new AbstractAction("Interfaces") {
        public void actionPerformed(ActionEvent e) {
          updateStructureTensors(0.0f,1.0f,1.0f);
        }
      });
    structureMenu.add(isotropicItem);
    structureMenu.add(linearItem);
    structureMenu.add(layersItem);
    structureMenu.add(interfacesItem);
    ButtonGroup structureGroup = new ButtonGroup();
    structureGroup.add(isotropicItem);
    structureGroup.add(linearItem);
    structureGroup.add(layersItem);
    structureGroup.add(interfacesItem);
    layersItem.setSelected(true);
    JMenuBar menuBar = new JMenuBar();
    menuBar.add(fileMenu);
    menuBar.add(modeMenu);
    menuBar.add(viewMenu);
    menuBar.add(structureMenu);
    _frame.setJMenuBar(menuBar);

    // Tool bar.
    JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
    toolBar.setRollover(true);
    toolBar.add(new ModeToggleButton(tzm));
    toolBar.add(new ModeToggleButton(pm));
    toolBar.add(new JButton(new AbstractAction("C") {
      public void actionPerformed(ActionEvent e) {
        _painting.clearAll();
        _paintView.set(_painting.getValues());
        _paintView.setClips(_valueMin,_valueMax);
      }
    }));
    toolBar.add(new JButton(new AbstractAction("F") {
      public void actionPerformed(ActionEvent e) {
        _painting.clearNotFixed();
        _paintView.set(_painting.getValues());
        _paintView.setClips(_valueMin,_valueMax);
      }
    }));
    toolBar.add(new JButton(new AbstractAction("E") {
      public void actionPerformed(ActionEvent e) {
        _painting.extrapolate();
        _paintView.set(_painting.getValues());
        _paintView.setClips(_valueMin,_valueMax);
      }
    }));
    toolBar.add(new JButton(new AbstractAction("I") {
      public void actionPerformed(ActionEvent e) {
        _painting.interpolate();
        _paintView.set(_painting.getValues());
        _paintView.setClips(_valueMin,_valueMax);
      }
    }));
    toolBar.add(new JButton(new AbstractAction("T") {
      public void actionPerformed(ActionEvent e) {
        _painting.extrapolate();
        _paintView.set(_painting.getTimes());
        _paintView.setPercentiles(0.0f,100.0f);
      }
    }));
    _frame.add(toolBar,BorderLayout.WEST);

    // Initially activate paint mode.
    pm.setActive(true);
  }

  public void showTensors(boolean show) {
    if (show) {
      _panel.getTile(0,0).addTiledView(_tensorsView);
    } else {
      _panel.getTile(0,0).removeTiledView(_tensorsView);
    }
  }

  // Actions.
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
  private class ShowTensorsAction extends AbstractAction {
    private ShowTensorsAction() {
      super("Tensors");
    }
    public void actionPerformed(ActionEvent event) {
      _show = !_show;
      showTensors(_show);
    }
    boolean _show = false;
  }
  private void updateStructureTensors(float alpha, float beta, float gamma) {
    _st = new StructureTensors(SIGMA,alpha,beta,gamma,_image);
    _painting.setTensors(_st);
    float[][][] x12 = getTensorEllipses(_n1,_n2,10,_st);
    float[][] x1 = x12[0];
    float[][] x2 = x12[1];
    _tensorsView.set(x1,x2);
  }

  private static float[][][] getTensorEllipses(
    int n1, int n2, int ns, EigenTensors2 et) 
  {
    int nt = 51;
    int m1 = 1+(n1-1)/ns;
    int m2 = 1+(n2-1)/ns;
    int j1 = (n1-1-(m1-1)*ns)/2;
    int j2 = (n2-1-(m2-1)*ns)/2;
    int nm = m1*m2;
    //double r = 0.45*ns;
    float[][] sm = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float[] s = et.getEigenvalues(i1,i2);
        sm[i2][i1] = s[1];
      }
    }
    sm = ArrayMath.copy(m1,m2,j1,j2,ns,ns,sm);
    float smq = Quantiler.estimate(0.05f,sm);
    double r = 0.45*ns*sqrt(smq);
    float[][] x1 = new float[nm][nt];
    float[][] x2 = new float[nm][nt];
    double dt = 2.0*PI/(nt-1);
    double ft = 0.0f;
    for (int i2=j2,im=0; i2<n2; i2+=ns) {
      double y2 = i2+r;
      for (int i1=j1; i1<n1; i1+=ns,++im) {
        float[] u = et.getEigenvectorU(i1,i2);
        float[] s = et.getEigenvalues(i1,i2);
        double u1 = u[0];
        double u2 = u[1];
        double v1 = -u2;
        double v2 =  u1;
        double su = s[0];
        double sv = s[1];
        su = max(su,smq);
        sv = max(sv,smq);
        double a = r/sqrt(sv);
        double b = r/sqrt(su);
        for (int it=0; it<nt; ++it) {
          double t = ft+it*dt;
          double cost = cos(t);
          double sint = sin(t);
          x1[im][it] = (float)(i1+b*cost*u1-a*sint*u2);
          x2[im][it] = (float)(i2+a*sint*u1+b*cost*u2);
        }
      }
    }
    return new float[][][]{x1,x2};
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static float[][] readImage(int n1, int n2, String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  private static float[][] gain(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        float ai = abs(xi);
        float si = (xi>=0.0f)?1.0f:-1.0f;
        y[i2][i1] = si*pow(ai,1.00f);
      }
    }
    return y;
  }

  private static class StructureTensors 
    extends EigenTensors2
    implements PaintingX.Tensors 
  {
    StructureTensors(float sigma, float[][] x) {
      this(sigma,-1.0f,1.0f,1.0f,x);
    }
    StructureTensors(
      float sigma, float alpha, float beta, float gamma, float[][] x) 
    {
      super(x[0].length,x.length);
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] u1 = new float[n2][n1];
      float[][] u2 = new float[n2][n1];
      float[][] su = new float[n2][n1];
      float[][] sv = new float[n2][n1];
      LocalOrientFilter lof = new LocalOrientFilter(sigma);
      lof.apply(x,null,u1,u2,null,null,su,sv,null);
      trace("su min="+ ArrayMath.min(su)+" max="+ ArrayMath.max(su));
      trace("sv min="+ ArrayMath.min(sv)+" max="+ ArrayMath.max(sv));
      float[][] sa = ArrayMath.pow(su,alpha);
      float[][] sb = ArrayMath.pow(ArrayMath.div(sv,su),beta);
      float[][] sc = ArrayMath.pow(ArrayMath.sub(1.0f,coherence(sigma,x)),gamma);
      su = ArrayMath.mul(sa,sc);
      sv = ArrayMath.mul(sb,su);
      trace("su min="+ ArrayMath.min(su)+" max="+ ArrayMath.max(su));
      trace("sv min="+ ArrayMath.min(sv)+" max="+ ArrayMath.max(sv));
      //SimplePlot.asPixels(su);
      //SimplePlot.asPixels(sv);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          setEigenvalues(i1,i2,su[i2][i1],sv[i2][i1]);
          setEigenvectorU(i1,i2,u1[i2][i1],u2[i2][i1]);
        }
      }
    }
  }

  private static float[][] coherence(double sigma, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    LocalOrientFilter lof1 = new LocalOrientFilter(sigma);
    LocalOrientFilter lof2 = new LocalOrientFilter(sigma*4);
    float[][] u11 = new float[n2][n1];
    float[][] u21 = new float[n2][n1];
    float[][] su1 = new float[n2][n1];
    float[][] sv1 = new float[n2][n1];
    float[][] u12 = new float[n2][n1];
    float[][] u22 = new float[n2][n1];
    float[][] su2 = new float[n2][n1];
    float[][] sv2 = new float[n2][n1];
    lof1.apply(x,null,u11,u21,null,null,su1,sv1,null);
    lof2.apply(x,null,u12,u22,null,null,su2,sv2,null);
    float[][] c = u11;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u11i = u11[i2][i1];
        float u21i = u21[i2][i1];
        float su1i = su1[i2][i1];
        float sv1i = sv1[i2][i1];
        float u12i = u12[i2][i1];
        float u22i = u22[i2][i1];
        float su2i = su2[i2][i1];
        float sv2i = sv2[i2][i1];
        float s111 = (su1i-sv1i)*u11i*u11i+sv1i;
        float s121 = (su1i-sv1i)*u11i*u21i     ;
        float s221 = (su1i-sv1i)*u21i*u21i+sv1i;
        float s112 = (su2i-sv2i)*u12i*u12i+sv2i;
        float s122 = (su2i-sv2i)*u12i*u22i     ;
        float s222 = (su2i-sv2i)*u22i*u22i+sv2i;
        float s113 = s111*s112+s121*s122;
        float s223 = s121*s122+s221*s222;
        float t1 = s111+s221;
        float t2 = s112+s222;
        float t3 = s113+s223;
        float t12 = t1*t2;
        c[i2][i1] = (t12>0.0f)?t3/t12:0.0f;
      }
    }
    return c;
  }

  private static void testImagePainterA() {
    int n1 = 251;
    int n2 = 357;
    float[][] image = readImage(n1,n2,"/data/seis/tp/tp73.dat");
    image = gain(image);
    ImagePainterX ip = new ImagePainterX(image);
    ip.setValueRange(0.0,1.0);
  }

  private static void testImagePainterB() {
    int n1 = 500;
    int n2 = 500;
    float[][] image = readImage(n1,n2,"/data/seis/atw/atwj1s.dat");
    image = gain(image);
    ImagePainterX ip = new ImagePainterX(image);
    ip.setValueRange(0.0,1.0);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testImagePainterA();
        testImagePainterB();
      }
    });
  }
}
