package ivus;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;

/**
 * A frame for viewing IVUS data.
 * @author Dave Hale
 * @version 2006.06.28
 */
public class IvusFrame extends JFrame {
  private static final long serialVersionUID = 1L;

  public IvusFrame(World world) {
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

  public static void main(String[] args) {
    Float3 f = loadFC();
    World world = new World();
    addSlices(world,f);
    IvusFrame ivusFrame = new IvusFrame(world);
    ivusFrame.setVisible(true);
  }

  private static final int SIZE = 600;
  private static final String DATA_DIR = "/data/ivus/";

  private static SimpleFloat3 loadFC() {
    System.out.print("loading ... ");
    int nx = 256;
    int ny = 256;
    //int nz = 1000;
    int nz = 1000;
    int kz = 2000;
    SimpleFloat3 sf3;
    try {
      ArrayFile af = new ArrayFile(DATA_DIR+"CS_patient/CS_matrix.imgs","r");
      //ArrayFile af = 
      //  new ArrayFile(DATA_DIR+"FC_patient/FC_2000_2999.imgs","r");
      af.seek(175+kz*nx*ny);
      byte[][][] b = ArrayMath.zerobyte(nx,ny,nz);
      af.readBytes(b);
      float[][][] f = ArrayMath.zerofloat(nz,ny,nx);
      for (int ix=0; ix<nx; ++ix) {
        for (int iy=0; iy<ny; ++iy) {
          for (int iz=0; iz<nz; ++iz) {
            int i = b[iz][iy][ix];
            if (i<0)
              i += 256;
            f[ix][iy][iz] = (float)i;
          }
        }
      }
      sf3 = new SimpleFloat3(f);
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }
    System.out.println("done");
    return sf3;
  }

  private static void addSlices(World world, Float3 f) {
    int nx = f.getN3();
    int ny = f.getN2();
    int nz = f.getN1();
    double dx = 1.0;
    double dy = 1.0;
    double dz = 1.0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
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
      ImagePanel ip = new ImagePanel(sz,sy,sx,f);
      ip.setColorModel(ColorMap.GRAY);
      System.out.println("clip min="+ip.getClipMin()+" max="+ip.getClipMax());
      aaf.addChild(ip);
      world.addChild(aaq);
    }
  }
}
