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
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import static edu.mines.jtk.ogl.Gl.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * 3D image painting using an anistropic brush shape.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.08.21
 */
public class ImagePainter3 {

  public ImagePainter3(float[][][] image, EigenTensors3 tensors) {
    _n1 = image[0][0].length;
    _n2 = image[0].length;
    _n3 = image.length;
    _image = image;
    _pt = tensors;
    _pb = new PaintBrush3(_n1,_n2,_n3,_pt);
    _pb.setSize(30);
    ImagePanelGroup ipg = new ImagePanelGroup(_image);
    ipg.setColorModel(ColorMap.GRAY);

    _world = new World();
    _world.addChild(ipg);
    makeFrame();
  }

  public static EigenTensors3 readTensors(String fileName) {
    try {
      FileInputStream fis = new FileInputStream(fileName);
      ObjectInputStream ois = new ObjectInputStream(fis);
      EigenTensors3 et = (EigenTensors3)ois.readObject();
      fis.close();
      return et;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static void writeTensors(EigenTensors3 et, String fileName) {
    try {
      FileOutputStream fos = new FileOutputStream(fileName);
      ObjectOutputStream oos = new ObjectOutputStream(fos);
      oos.writeObject(et);
      fos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static float SIGMA = 3.0f;

  private static final int SIZE = 800;

  private int _n1,_n2,_n3;
  private float[][][] _image;
  private EigenTensors3 _pt;
  private JFrame _frame;
  private ViewCanvas _canvas;
  private World _world;
  private PaintBrush3 _pb;
  private TriangleGroup _pbtg;

  private void makeFrame() {
    OrbitView view = new OrbitView(_world);
    view.setAxesOrientation(View.AxesOrientation.XRIGHT_YOUT_ZDOWN);
    _canvas = new ViewCanvas(view);
    _canvas.setView(view);

    // Modes.
    ModeManager mm = new ModeManager();
    mm.add(_canvas);
    OrbitViewMode ovm = new OrbitViewMode(mm);
    SelectDragMode sdm = new SelectDragMode(mm);
    PaintBrushMode pbm = new PaintBrushMode(mm);

    JPopupMenu.setDefaultLightWeightPopupEnabled(false);
    ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);

    // File menu.
    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic('F');
    Action exitAction = new AbstractAction("Exit") {
      private static final long serialVersionUID = 1L;
      public void actionPerformed(ActionEvent event) {
        System.exit(0);
      }
    };
    JMenuItem exitItem = fileMenu.add(exitAction);
    exitItem.setMnemonic('x');

    // Mode menu.
    JMenu modeMenu = new JMenu("Mode");
    modeMenu.setMnemonic('M');
    JMenuItem ovmItem = new ModeMenuItem(ovm);
    modeMenu.add(ovmItem);
    JMenuItem sdmItem = new ModeMenuItem(sdm);
    modeMenu.add(sdmItem);
    JMenuItem pbmItem = new ModeMenuItem(pbm);
    modeMenu.add(pbmItem);

    // Paint menu.
    JMenu paintMenu = new JMenu("Paint");
    JMenuItem isotropicItem = new JRadioButtonMenuItem(
      new AbstractAction("Isotropic") {
        public void actionPerformed(ActionEvent e) {
          updatePaintTensors(0.0f,0.0f,0.0f);
        }
      });
    JMenuItem layersItem = new JRadioButtonMenuItem(
      new AbstractAction("Layers") {
        public void actionPerformed(ActionEvent e) {
          updatePaintTensors(1.0f,1.0f,1.0f);
        }
      });
    paintMenu.add(isotropicItem);
    paintMenu.add(layersItem);
    ButtonGroup paintGroup = new ButtonGroup();
    paintGroup.add(isotropicItem);
    paintGroup.add(layersItem);
    layersItem.setSelected(true);

    // Menu bar.
    JMenuBar menuBar = new JMenuBar();
    menuBar.add(fileMenu);
    menuBar.add(modeMenu);
    menuBar.add(paintMenu);

    // Tool bar.
    JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
    toolBar.setRollover(true);
    toolBar.add(new ModeToggleButton(ovm));
    toolBar.add(new ModeToggleButton(sdm));
    toolBar.add(new ModeToggleButton(pbm));
    ovm.setActive(true);

    // Frame.
    _frame = new JFrame();
    _frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    _frame.setSize(new Dimension(SIZE,SIZE));
    _frame.add(_canvas,BorderLayout.CENTER);
    _frame.add(toolBar,BorderLayout.WEST);
    _frame.setJMenuBar(menuBar);
    _frame.setVisible(true);
  }

  private PickResult pick(MouseEvent event) {
    ViewCanvas canvas = (ViewCanvas)event.getSource();
    View view = canvas.getView();
    if (view==null)
      return null;
    World world = view.getWorld();
    if (world==null)
      return null;
    PickContext pc = new PickContext(event);
    world.pick(pc);
    PickResult pickResult = pc.getClosest();
    if (pickResult!=null) {
      Point3 pointLocal = pickResult.getPointLocal();
      Point3 pointWorld = pickResult.getPointWorld();
      System.out.println("Painting at:");
      System.out.println("  local="+pointLocal);
      System.out.println("  world="+pointWorld);
    } else {
      System.out.println("Painting nothing");
    }
    return pickResult;
  }

  private class PaintBrushMode extends Mode {
    public PaintBrushMode(ModeManager modeManager) {
      super(modeManager);
      setName("Paint");
      Class<PaintBrushMode> cls = PaintBrushMode.class;
      setIcon(loadIcon(cls,"resources/PaintIcon16.png"));
      setMnemonicKey(KeyEvent.VK_B);
      setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_B,0));
      setShortDescription("Paint brush");
    }
    protected void setActive(Component component, boolean active) {
      if (component instanceof ViewCanvas) {
        if (active) {
          component.addMouseListener(_ml);
        } else {
          component.removeMouseListener(_ml);
        }
      }
    }
    private boolean _painting;
    private MouseConstrained _mouseConstrained;
    private MouseListener _ml = new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        PickResult pr = pick(e);
        if (pr!=null) {
          Node node = pr.getNode(AxisAlignedQuad.class);
          if (node!=null) {
            _painting = true;
            _canvas.addMouseMotionListener(_mml);
            AxisAlignedQuad quad = (AxisAlignedQuad)node;
            AxisAlignedFrame frame = quad.getFrame();
            Axis axis = frame.getAxis();
            Point3 origin = pr.getPointWorld();
            Vector3 normal = null;
            if (axis==Axis.X) {
              normal = new Vector3(1.0,0.0,0.0);
            } else if (axis==Axis.Y) {
              normal = new Vector3(0.0,1.0,0.0);
            } else if (axis==Axis.Z) {
              normal = new Vector3(0.0,0.0,1.0);
            }
            Plane plane = new Plane(origin,normal);
            Matrix44 worldToPixel = pr.getWorldToPixel();
            _mouseConstrained = new MouseOnPlane(e,origin,plane,worldToPixel);
            paintAt(origin);
          }
        }
      }
      public void mouseReleased(MouseEvent e) {
        if (_painting) {
          _mouseConstrained = null;
          _canvas.removeMouseMotionListener(_mml);
          _painting = false;
        }
      }
    };
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        Point3 point = _mouseConstrained.getPoint(e);
        paintAt(point);
      }
    };

    private void paintAt(Point3 point) {
      int i1 = max(0,min(_n1-1,(int)(point.z+0.5)));
      int i2 = max(0,min(_n2-1,(int)(point.y+0.5)));
      int i3 = max(0,min(_n3-1,(int)(point.x+0.5)));
      //trace("paintAt: i1="+i1+" i2="+i2+" i3="+i3);
      _pb.setLocation(i1,i2,i3);
      if (_pbtg!=null)
        _world.removeChild(_pbtg);
      PaintBrush3.Contour contour = _pb.getContour();
      _pbtg = new TriangleGroup(contour.i,contour.x,contour.u);
      StateSet states = new StateSet();
      ColorState cs = new ColorState();
      cs.setColor(Color.CYAN);
      states.add(cs);
      LightModelState lms = new LightModelState();
      lms.setTwoSide(true);
      states.add(lms);
      MaterialState ms = new MaterialState();
      ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
      ms.setSpecular(Color.WHITE);
      ms.setShininess(100.0f);
      states.add(ms);
      _pbtg.setStates(states);
      _world.addChild(_pbtg);
    }
  }

  private static float[][][] readImage(
    int n1, int n2, int n3, String fileName) 
  {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[][][] x = new float[n3][n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  private void updatePaintTensors(float alpha, float beta, float gamma) {
    _pt = readTensors("/data/seis/tp/pt3s211.dat");
    _pb = new PaintBrush3(_n1,_n2,_n3,_pt);
  }

  private static void testImagePainter() {
    int n1 = 251;
    int n2 = 161;
    int n3 = 357;
    float[][][] image = readImage(n1,n2,n3,"/data/seis/tp/tp3s.dat");
    image = Array.add(image,Array.mul(0.001f,Array.randfloat(n1,n2,n3)));
    EigenTensors3 tensors = readTensors("/data/seis/tp/et3s211.dat");
    ImagePainter3 ip = new ImagePainter3(image,tensors);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testImagePainter();
      }
    });
  }
}
