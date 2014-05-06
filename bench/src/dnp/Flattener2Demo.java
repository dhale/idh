package dnp;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Interactive demo of flattening with constraints.
 */
public class Flattener2Demo {

  /**
   * Runs the program.
   * @param args arguments (ignored).
   */
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Flattener2Demo();
      }
    });
  }

  private Flattener2Demo() {
    int n1 = 251;
    int n2 = 357;
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    String fileName = "/data/seis/tpd/csm/oldslices/tp73.dat";
    float[][] f = new float[n2][n1];
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    float[][] g = copy(f);
    float[][] t = rampfloat(0.0f,1.0f,0.0f,n1,n2);
    new Plot(s1,s2,f,g,t);
  }

  // Control points are pairs of coordinates (x1,x2).
  private static class ControlPoint {
    double x1,x2;
    ControlPoint(double x1, double x2) {
      this.x1 = x1;
      this.x2 = x2;
    }
    double distanceTo(double x1, double x2) {
      double d1 = this.x1-x1;
      double d2 = this.x2-x2;
      return sqrt(d1*d1+d2*d2);
    }
    double distanceTo(ControlPoint cp) {
      return distanceTo(cp.x1,cp.x2);
    }
    public boolean equals(Object o) {
      ControlPoint that = (ControlPoint)o;
      return this.x1==that.x1 && this.x2==that.x2;
    }
    public String toString() {
      return "("+x1+","+x2+")";
    }
  }
  // A named and colored list of control points with a view.
  private static class ControlSet {
    String name;
    Color color;
    int key;
    ArrayList<ControlPoint> cps;
    PointsView pv;
    ControlSet(String name, Color color, int key) {
      this.name = name;
      this.color = color;
      this.key = key;
      this.cps = new ArrayList<ControlPoint>();
    }
    void add(ControlPoint cp) {
      if (!cps.contains(cp)) {
        cps.add(cp);
        updatePointsView();
      }
    }
    void remove(ControlPoint cp) {
      if (cps.contains(cp)) {
        cps.remove(cp);
        updatePointsView();
      }
    }
    void move(ControlPoint cpOld, ControlPoint cpNew) {
      if (cps.contains(cpOld)) {
        cps.remove(cpOld);
        cps.add(cpNew);
        updatePointsView();
      }
    }
    PointsView updatePointsView() {
      int ncp = cps.size();
      float[] x1 = new float[ncp];
      float[] x2 = new float[ncp];
      for (int icp=0; icp<ncp; ++icp) {
        ControlPoint cp = cps.get(icp);
        x1[icp] = (float)cp.x1;
        x2[icp] = (float)cp.x2;
      }
      if (pv!=null) {
        pv.set(x1,x2);
      } else {
        pv = new PointsView(x1,x2);
        pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
        pv.setLineStyle(PointsView.Line.NONE);
        pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE);
        pv.setMarkSize(18.0f);
        pv.setMarkColor(color);
      }
      return pv;
    }
  }
  // Six predefined control sets, initially empty.
  private ControlSet[] _css = {
    new ControlSet("R",Color.RED,KeyEvent.VK_R),
    new ControlSet("G",Color.GREEN,KeyEvent.VK_G),
    new ControlSet("B",Color.BLUE,KeyEvent.VK_B),
    new ControlSet("C",Color.CYAN,KeyEvent.VK_C),
    new ControlSet("M",Color.MAGENTA,KeyEvent.VK_M),
    new ControlSet("Y",Color.YELLOW,KeyEvent.VK_Y),
  };
  Tile _cssTile; // tile in which control sets are plotted
  // Currently active control set, if not null.
  private ControlSet _csa = null;

  private ControlPoint findNearestControlPoint(ControlPoint cp) {
    ControlPoint pmin = null;
    double dmin = 0.0;
    for (ControlPoint p: _csa.cps) {
      double d = cp.distanceTo(p);
      if (pmin==null || d<dmin) {
        pmin = p;
        dmin = d;
      }
    }
    return pmin;
  }

  ///////////////////////////////////////////////////////////////////////////

  // Location and size of plot.
  private static final int PLOT_X = 50;
  private static final int PLOT_Y = 50;
  private static final int PLOT_WIDTH = 950;
  private static final int PLOT_HEIGHT = 650;

  // One plot for everything.
  private class Plot {

    private Plot(Sampling s1, Sampling s2,
                 float[][] f, float[][] g, float[][] t)
    {
      PlotPanel pp = new PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT);
      pp.setHLabel(0,"Inline (samples)");
      pp.setHLabel(1,"Inline (samples)");
      pp.setVLabel("Time (samples)");
      pp.addPixels(0,0,s1,s2,f);
      pp.addPixels(0,1,s1,s2,g);
      ContoursView cv = pp.addContours(0,0,s1,s2,t);
      cv.setColorModel(ColorMap.JET);
      pp.addColorBar(cv,"Relative geologic time");
      for (ControlSet cs:_css) {
        pp.addTiledView(0,0,cs.updatePointsView());
      }
      _cssTile = _css[0].pv.getTile();

      // A plot frame has a mode for zooming in tiles or tile axes.
      PlotFrame pf = new PlotFrame(pp);
      TileZoomMode tzm = pf.getTileZoomMode();
      tzm.setActive(true);

      // Add modes for editing control sets.
      ModeManager modeManager = pf.getModeManager();
      int ncpm = _css.length;
      ControlPointMode[] cpms = new ControlPointMode[ncpm];
      for (int icpm=0; icpm<ncpm; ++icpm)
        cpms[icpm] = new ControlPointMode(modeManager,_css[icpm]);

      // The menu bar includes a mode menu for selecting a mode.
      JMenu fileMenu = new JMenu("File");
      fileMenu.setMnemonic('F');
      fileMenu.add(new ExitAction()).setMnemonic('x');
      JMenu modeMenu = new JMenu("Mode");
      modeMenu.setMnemonic('M');
      modeMenu.add(new ModeMenuItem(tzm));
      for (int icpm=0; icpm<ncpm; ++icpm)
        modeMenu.add(new ModeMenuItem(cpms[icpm]));
      JMenuBar menuBar = new JMenuBar();
      menuBar.add(fileMenu);
      menuBar.add(modeMenu);
      pf.setJMenuBar(menuBar);

      // The tool bar includes toggle buttons for selecting a mode.
      JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
      toolBar.setRollover(true);
      toolBar.add(new ModeToggleButton(tzm));
      for (int icpm=0; icpm<ncpm; ++icpm)
        toolBar.add(new ModeToggleButton(cpms[icpm]));
      pf.add(toolBar, BorderLayout.WEST);

      // Make the plot frame visible.
      pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      pf.setLocation(PLOT_X,PLOT_Y);
      pf.setSize(PLOT_WIDTH,PLOT_HEIGHT);
      pf.setFontSizeForSlide(1.0,1.0);
      pf.setVisible(true);
    }
  }

  ///////////////////////////////////////////////////////////////////////////

  // A mode for adding, removing, or moving control points.
  private class ControlPointMode extends Mode {

    public ControlPointMode(ModeManager modeManager, ControlSet cs) {
      super(modeManager);
      _cs = cs;
      setName(cs.name);
      setIcon(loadIcon(Flattener2Demo.class,"resources/Flat16.png"));
      setMnemonicKey(cs.key);
      setAcceleratorKey(KeyStroke.getKeyStroke(cs.key, 0));
      setShortDescription("Add (Shift), remove (Ctrl), or drag points");
    }

    protected void setActive(Component component, boolean active) {
      if (component instanceof Tile) {
        Tile tile = (Tile)component;
        if (tile==_cssTile) {
          if (active) {
            _csa = _cs; // make our control set the active one
            component.addMouseListener(_ml);
          } else {
            if (_csa==_cs) // if our control set is the active one,
              _csa = null; // then no control set will be active
            component.removeMouseListener(_ml);
          }
        }
      }
    }

    private ControlSet _cs; // control set for this mode
    private ControlPoint _cpEdit; // if editing, last point
    private boolean _editing; // true, if currently editing
    private Tile _tile; // tile in which editing began

    // Handles mouse pressed and released events.
    private MouseListener _ml = new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        if (e.isShiftDown()) {
          add(e);
        } else if (e.isControlDown()) {
          remove(e);
        } else {
          if (beginEdit(e)) {
            _editing = true;
            _tile.addMouseMotionListener(_mml);
          }
        }
      }
      public void mouseReleased(MouseEvent e) {
        if (_editing) {
          _tile.removeMouseMotionListener(_mml);
          endEdit(e);
          _editing = false;
        }
      }
    };

    // Handles mouse dragged events.
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        if (_editing)
          duringEdit(e);
      }
    };

    // Converts an point (x,y) in pixels to a complex number z.
    private ControlPoint pointToControlPoint(int x, int y) {
      Transcaler ts = _tile.getTranscaler();
      Projector hp = _tile.getHorizontalProjector();
      Projector vp = _tile.getVerticalProjector();
      double xu = ts.x(x);
      double yu = ts.y(y);
      double xv = hp.v(xu);
      double yv = vp.v(yu);
      return new ControlPoint(yv,xv);
    }

    // Converts control point p to a point (x,y) in pixels.
    private Point controlPointToPoint(ControlPoint cp) {
      Transcaler ts = _tile.getTranscaler();
      Projector hp = _tile.getHorizontalProjector();
      Projector vp = _tile.getVerticalProjector();
      double xu = hp.u(cp.x2);
      double yu = vp.u(cp.x1);
      int xp = ts.x(xu);
      int yp = ts.y(yu);
      return new Point(xp,yp);
    }

    // Determines whether a specified point (x,y) is within a small
    // fixed number of pixels to the specified control point cp.
    private boolean closeEnough(int x, int y, ControlPoint cp) {
      Point p = controlPointToPoint(cp);
      return abs(p.x-x)<12 && abs(p.y-y)<12;
    }

    // Adds a control point at mouse coordinates (x,y).
    private void add(MouseEvent e) {
      _tile = (Tile)e.getSource();
      int x = e.getX();
      int y = e.getY();
      ControlPoint cp = pointToControlPoint(x,y);
      _cs.add(cp);
    }

    // Removes a control point, if mouse (x,y) is close enough to one.
    private void remove(MouseEvent e) {
      _tile = (Tile)e.getSource();
      int x = e.getX();
      int y = e.getY();
      ControlPoint cp = pointToControlPoint(x,y);
      ControlPoint cq = findNearestControlPoint(cp);
      if (cq!=null && closeEnough(x,y,cq))
        _cs.remove(cq);
    }

    // Begins editing of an existing control point, if close enough.
    // Returns true, if close enough so that we have begun editing;
    // false, otherwise.
    private boolean beginEdit(MouseEvent e) {
      _tile = (Tile)e.getSource();
      int x = e.getX();
      int y = e.getY();
      ControlPoint cp = pointToControlPoint(x,y);
      ControlPoint cq = findNearestControlPoint(cp);
      if (cq!=null && closeEnough(x,y,cq)) {
        _cs.move(cq,cp);
        _cpEdit = cp;
        return true;
      }
      return false;
    }

    // Called while a pole or zero is being dragged during edited.
    private void duringEdit(MouseEvent e) {
      int x = e.getX();
      int y = e.getY();
      ControlPoint cp = pointToControlPoint(x,y);
      _cs.move(_cpEdit,cp);
      _cpEdit = cp;
    }

    // Called when done editing a pole or zero.
    private void endEdit(MouseEvent e) {
      duringEdit(e);
      _editing = false;
    }
  }

  private class ExitAction extends AbstractAction {
    private ExitAction() {
      super("Exit");
    }
    public void actionPerformed(ActionEvent event) {
      System.exit(0);
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
