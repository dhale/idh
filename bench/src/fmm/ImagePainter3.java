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

  public ImagePainter3(float[][][] image) {
    _n1 = image[0][0].length;
    _n2 = image[0].length;
    _n3 = image.length;
    _image = image;
    //_pt = new PaintTensors(SIGMA,2.0f,1.0f,1.0f,_image);
    //writeTensors(_pt,"/data/seis/tp/tpet.dat");
    _pt = readTensors("/data/seis/tp/tpet.dat");
    _pb = new PaintBrush3(_n1,_n2,_n3,_pt);

    Sampling s1 = new Sampling(_n1);
    Sampling s2 = new Sampling(_n2);
    Sampling s3 = new Sampling(_n3);
    ImagePanelGroup ipg = new ImagePanelGroup(s3,s2,s1,
                                              new SimpleFloat3(_image));
    ipg.setColorModel(ColorMap.GRAY);

    _world = new World();
    _world.addChild(ipg);
    makeFrame();
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static float SIGMA = 3.0f;

  private static final int SIZE = 800;

  private int _n1,_n2,_n3;
  private float[][][] _image;
  private PaintTensors _pt;
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
    PaintMode pm = new PaintMode(mm);

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
    toolBar.add(new ModeToggleButton(pm));
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

  private class PaintMode extends Mode {
    public PaintMode(ModeManager modeManager) {
      super(modeManager);
      setName("Paint");
      setIcon(loadIcon(PaintMode.class,"resources/PaintIcon16.png"));
      setMnemonicKey(KeyEvent.VK_P);
      setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_P,0));
      setShortDescription("Paint");
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
    private MouseListener _ml = new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        PickResult pr = pick(e);
        if (pr!=null) {
          Node node = pr.getNode(AxisAlignedPanel.class);
          //if (node!=null) {
            _canvas.addMouseMotionListener(_mml);
            _painting = true;
            paintAt(pr);
          //}
        }
      }
      public void mouseReleased(MouseEvent e) {
        if (_painting) {
          _canvas.removeMouseMotionListener(_mml);
          _painting = false;
        }
      }
    };
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        //PickResult pickResult = pick(e);
      }
    };
    private void paintAt(PickResult pr) {
      Point3 pw = pr.getPointWorld();
      int i1 = max(0,min(_n1-1,(int)(pw.z+0.5)));
      int i2 = max(0,min(_n2-1,(int)(pw.y+0.5)));
      int i3 = max(0,min(_n3-1,(int)(pw.x+0.5)));
      trace("paintAt: i1="+i1+" i2="+i2+" i3="+i3);
      _pb.setLocation(i1,i2,i3);
      if (_pbtg!=null)
        _world.removeChild(_pbtg);
      PaintBrush3.Contour contour = _pb.getContour();
      Array.dump(contour.x);
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
    _pt = new PaintTensors(SIGMA,alpha,beta,gamma,_image);
    _pb = new PaintBrush3(_n1,_n2,_n3,_pt);
  }

  private static class PaintTensors extends EigenTensors3 {
    PaintTensors(int n1, int n2, int n3) {
      super(n1,n2,n3,false);
    }
    PaintTensors(float sigma, float[][][] x) {
      this(sigma,1.0f,1.0f,1.0f,x);
    }
    PaintTensors(
      float sigma, float alpha, float beta, float gamma, float[][][] x) 
    {
      super(x[0][0].length,x[0].length,x.length,false);
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      float[][][] u1 = new float[n3][n2][n1];
      float[][][] u2 = new float[n3][n2][n1];
      float[][][] u3 = new float[n3][n2][n1];
      float[][][] w1 = new float[n3][n2][n1];
      float[][][] w2 = new float[n3][n2][n1];
      float[][][] w3 = new float[n3][n2][n1];
      float[][][] su = new float[n3][n2][n1];
      float[][][] sv = new float[n3][n2][n1];
      float[][][] sw = new float[n3][n2][n1];
      LocalOrientFilter lof = new LocalOrientFilter(sigma);
      lof.apply(x,null,null,
                u1,u2,u3,
                null,null,null,
                w1,w2,w3,
                su,sv,sw,
                null,null);
      float[][][] dc = Array.pow(Array.sub(1.0f,coherence(sigma,x)),-gamma);
      float[][][] du = Array.mul(dc,Array.pow(su,-alpha));
      float[][][] dv = Array.mul(du,Array.pow(Array.div(sv,su),-beta));
      float[][][] dw = Array.mul(du,Array.pow(Array.div(sw,su),-beta));
      float ds = 1.0f/Array.max(dw);
      du = Array.mul(ds,du);
      dv = Array.mul(ds,dv);
      dw = Array.mul(ds,dw);
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            setEigenvectorU(i1,i2,i3,
                            u1[i3][i2][i1],u2[i3][i2][i1],u3[i3][i2][i1]);
            setEigenvectorW(i1,i2,i3,
                            w1[i3][i2][i1],w2[i3][i2][i1],w3[i3][i2][i1]);
            setEigenvalues(i1,i2,i3,
                           du[i3][i2][i1],dv[i3][i2][i1],dw[i3][i2][i1]);
          }
        }
      }
    }
  }
  private static PaintTensors readTensors(String fileName) {
    PaintTensors pt = null;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      int n1 = ais.readInt();
      int n2 = ais.readInt();
      int n3 = ais.readInt();
      pt = new PaintTensors(n1,n2,n3);
      float[] u = new float[n1*3];
      float[] w = new float[n1*3];
      float[] e = new float[n1*3];
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          ais.readFloats(u);
          ais.readFloats(w);
          ais.readFloats(e);
          for (int i1=0,i=0; i1<n1; ++i1,i+=3) {
            pt.setEigenvectorU(i1,i2,i3,u[i],u[i+1],u[i+2]);
            pt.setEigenvectorW(i1,i2,i3,w[i],w[i+1],w[i+2]);
            pt.setEigenvalues(i1,i2,i3,e[i],e[i+1],e[i+2]);
          }
        }
      }
      ais.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    return pt;
  }
  private static void writeTensors(EigenTensors3 et, String fileName) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      int n1 = et.getN1();
      int n2 = et.getN2();
      int n3 = et.getN3();
      aos.writeInt(n1);
      aos.writeInt(n2);
      aos.writeInt(n3);
      float[] u = new float[3];
      float[] w = new float[3];
      float[] e = new float[3];
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            et.getEigenvectorU(i1,i2,i3,u);
            et.getEigenvectorW(i1,i2,i3,w);
            et.getEigenvalues(i1,i2,i3,e);
            aos.writeFloats(u);
            aos.writeFloats(w);
            aos.writeFloats(e);
          }
        }
      }
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  private static float[][][] coherence(double sigma, float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    return Array.fillfloat(1.0f,n1,n2,n3);
    /*
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
    */
  }

  private static void testImagePainter() {
    int n1 = 251;
    int n2 = 161;
    int n3 = 357;
    float[][][] image = readImage(n1,n2,n3,"/data/seis/tp/tp3s.dat");
    ImagePainter3 ip = new ImagePainter3(image);
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
