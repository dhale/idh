package sw;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.ogl.Gl.*;

/**
 * A frame for viewing sw data.
 * @author Dave Hale
 * @version 2006.09.18
 */
public class SwView extends JFrame {
  private static final long serialVersionUID = 1L;

  // Inline sampling
  private static final int NX = 367;
  private static final float DX = 0.025f;
  private static final float FX = 0.0f;
  //private static final float DX = 1.0f;
  //private static final float FX = 1325.0f;

  // Crossline sampling
  private static final int NY = 623;
  private static final float DY = 0.025f;
  private static final float FY = 0.0f;
  //private static final float DY = 2.0f;
  //private static final float FY = 1316.0f;

  // Time sampling
  private static final int NZ = 301;
  private static final float DZ = 0.004f;
  private static final float FZ = 3.6f;

  private static final int SIZE = 600;
  private static final String DATA_DIR = "/data/seis/sw/";

  public SwView(World world) {
    OrbitView view = (world!=null)?new OrbitView(world):new OrbitView();
    view.setAxesOrientation(View.AxesOrientation.XRIGHT_YOUT_ZDOWN);
    ViewCanvas canvas = new ViewCanvas(view);
    canvas.setView(view);

    ModeManager mm = new ModeManager();
    mm.add(canvas);
    OrbitViewMode ovm = new OrbitViewMode(mm);
    SelectDragMode sdm = new SelectDragMode(mm);

    JPopupMenu.setDefaultLightWeightPopupEnabled(false);
    ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);

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

    JMenu modeMenu = new JMenu("Mode");
    modeMenu.setMnemonic('M');
    JMenuItem ovmItem = new ModeMenuItem(ovm);
    modeMenu.add(ovmItem);
    JMenuItem sdmItem = new ModeMenuItem(sdm);
    modeMenu.add(sdmItem);

    JMenuBar menuBar = new JMenuBar();
    menuBar.add(fileMenu);
    menuBar.add(modeMenu);

    JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
    toolBar.setRollover(true);
    JToggleButton ovmButton = new ModeToggleButton(ovm);
    toolBar.add(ovmButton);
    JToggleButton sdmButton = new ModeToggleButton(sdm);
    toolBar.add(sdmButton);

    ovm.setActive(true);

    this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    this.setSize(new Dimension(SIZE,SIZE));
    this.add(canvas,BorderLayout.CENTER);
    this.add(toolBar,BorderLayout.WEST);
    this.setJMenuBar(menuBar);
  }

 
  private static SimpleFloat3 loadData(String fileName) {
    return loadData(fileName,false);
  }
  private static SimpleFloat3 loadData(String fileName, boolean bytes) {
    System.out.print("loading ... ");
    int nx = NX;
    int ny = NY;
    int nz = NZ;
    SimpleFloat3 sf3;
    try {
      ArrayFile af = new ArrayFile(DATA_DIR+fileName,"r");
      float[][][] f = Array.zerofloat(nz,ny,nx);
      if (bytes) {
        byte[][][] b = Array.zerobyte(nz,ny,nx);
        af.readBytes(b);
        f = Array.tofloat(b);
      } else {
        af.readFloats(f);
      }
      sf3 = new SimpleFloat3(f);
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }
    System.out.println("done");
    return sf3;
  }

  private static void addHorizons(World world, float[][][] h) {
    Color[] colors = {Color.YELLOW,Color.CYAN};
    int nh = h.length;
    for (int ih=0; ih<nh; ++ih) {
      float[] xyz = HorizonReader.xyz(h[ih]);
      TriangleGroup tg = new TriangleGroup(true,xyz);
      StateSet states = new StateSet();
      ColorState cs = new ColorState();
      cs.setColor(colors[ih]);
      states.add(cs);
      LightModelState lms = new LightModelState();
      lms.setTwoSide(true);
      states.add(lms);
      MaterialState ms = new MaterialState();
      ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
      ms.setSpecular(Color.white);
      ms.setShininess(100.0f);
      states.add(ms);
      tg.setStates(states);
      world.addChild(tg);
    }
  }

  private static void addSlices(
    World world, Float3 f, 
    float clip, IndexColorModel icm)
  {
    int nx = f.getN3();
    int ny = f.getN2();
    int nz = f.getN1();
    double dx = DX;
    double dy = DY;
    double dz = DZ;
    double fx = FX;
    double fy = FY;
    double fz = FZ;
    double lx = fx+(nx-1)*dx;
    double ly = fy+(ny-1)*dy;
    double lz = fz+(nz-1)*dz;
    Sampling sx = new Sampling(nx,dx,fx);
    Sampling sy = new Sampling(ny,dy,fy);
    Sampling sz = new Sampling(nz,dz,fz);
    Point3 qmin = new Point3(fx,fy,fz);
    Point3 qmax = new Point3(lx,ly,lz);
    Axis[] axes = new Axis[]{Axis.X,Axis.Y,Axis.Z};
    for (int iaxis=0; iaxis<axes.length; ++iaxis) {
      AxisAlignedQuad aaq = new AxisAlignedQuad(axes[iaxis],qmin,qmax);
      AxisAlignedFrame aaf = aaq.getFrame();
      ImagePanel ip = new ImagePanel(sx,sy,sz,f);
      ip.setClips(-clip,clip);
      ip.setColorModel(icm);
      System.out.println("clip min="+ip.getClipMin()+" max="+ip.getClipMax());
      aaf.addChild(ip);
      world.addChild(aaq);
    }
  }

  public static void main(String[] args) {
    World world = new World();
    //float[][][] h = HorizonReader.readTopBase();
    //addHorizons(world,h);
    addSlices(world,loadData("sw02a.dat"),5000.0f,ColorMap.GRAY);
    //addSlices(world,loadData("sw04a.dat"),5000.0f,ColorMap.GRAY);
    //addSlices(world,loadData("sw04a13.dat"),5000.0f,ColorMap.GRAY);
    //addSlices(world,loadData("r1.dat"),5000.0f,ColorMap.GRAY);
    //addSlices(world,loadData("r13.dat"),5000.0f,ColorMap.GRAY);
    addSlices(world,loadData("u1.dat"),1.0f,ColorMap.RED_WHITE_BLUE);
    addSlices(world,loadData("u3.dat"),0.5f,ColorMap.RED_WHITE_BLUE);
    SwView swView = new SwView(world);
    swView.setVisible(true);
  }
}
