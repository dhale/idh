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

import fault.*;

/**
 * Interactive demo of flattening with constraints.
 */
public class Flattener2Demo {

  private Flattener2C.Mappings _mappings; // mapping for time/depth <-> RGT
  private ConstraintSet[] _css = { // 6 empty sets of constraint points
    new ConstraintSet("R",Color.RED,KeyEvent.VK_R),
    new ConstraintSet("G",Color.GREEN,KeyEvent.VK_G),
    new ConstraintSet("B",Color.BLUE,KeyEvent.VK_B),
    new ConstraintSet("C",Color.CYAN,KeyEvent.VK_C),
    new ConstraintSet("M",Color.MAGENTA,KeyEvent.VK_M),
    new ConstraintSet("Y",Color.YELLOW,KeyEvent.VK_Y),
  };
  private ConstraintSet _csa = null; // currently active constraint set
  private Tile _cssTile; // tile in which constraint sets are plotted
  private Flattener2C _flattener; // for updating mappings with constraints
  private MappingsUpdater _mappingsUpdater; // runs in Swing worker thread
  private ContoursView _contoursView; // contours used for horizons
  private Sampling _s1,_s2; // samplings for image
  private float[][] _p2,_el; // slopes and linearities

  ///////////////////////////////////////////////////////////////////////////

  /**
   * Constraint points are pairs of coordinates (x1,x2).
   */
  private static class ConstraintPoint {
    double x1,x2;
    ConstraintPoint(double x1, double x2) {
      this.x1 = x1;
      this.x2 = x2;
    }
    double distanceTo(double x1, double x2) {
      double d1 = this.x1-x1;
      double d2 = this.x2-x2;
      return sqrt(d1*d1+d2*d2);
    }
    double distanceTo(ConstraintPoint cp) {
      return distanceTo(cp.x1,cp.x2);
    }
    public boolean equals(Object o) {
      ConstraintPoint that = (ConstraintPoint)o;
      return this.x1==that.x1 && this.x2==that.x2;
    }
    public String toString() {
      return "("+x1+","+x2+")";
    }
  }

  /**
   * A named and colored list of constraint points with a view.
   */
  private class ConstraintSet {
    String name;
    Color color;
    int key;
    ArrayList<ConstraintPoint> cps;
    PointsView pv;
    ConstraintSet(String name, Color color, int key) {
      this.name = name;
      this.color = color;
      this.key = key;
      this.cps = new ArrayList<ConstraintPoint>();
    }
    void add(ConstraintPoint cp) {
      if (!cps.contains(cp)) {
        cps.add(cp);
        updatePointsView();
        updateMappings();
      }
    }
    void remove(ConstraintPoint cp) {
      if (cps.contains(cp)) {
        cps.remove(cp);
        updatePointsView();
        updateMappings();
      }
    }
    void move(ConstraintPoint cpOld, ConstraintPoint cpNew) {
      if (cps.contains(cpOld)) {
        cps.remove(cpOld);
        cps.add(cpNew);
        updatePointsView();
        updateMappings();
      }
    }
    PointsView updatePointsView() {
      int ncp = cps.size();
      float[] x1 = new float[ncp];
      float[] x2 = new float[ncp];
      for (int icp=0; icp<ncp; ++icp) {
        ConstraintPoint cp = cps.get(icp);
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

  private float[][][] getConstraints() {
    int ncs = _css.length;
    float[][][] cs = new float[2][ncs][];
    for (int ics=0; ics<ncs; ++ics) {
      ArrayList<ConstraintPoint> cps = _css[ics].cps;
      int ncp = cps.size();
      float[] x1 = cs[0][ics] = new float[ncp];
      float[] x2 = cs[1][ics] = new float[ncp];
      for (int icp=0; icp<ncp; ++icp) {
        ConstraintPoint cp = cps.get(icp);
        x1[icp] = (float)cp.x1;
        x2[icp] = (float)cp.x2;
      }
    }
    return cs;
  }

  private void updateMappings() {
    if (_mappingsUpdater!=null) {
      trace("updateMappings: cancelling");
      _mappingsUpdater.cancel(false);
    }
    _mappingsUpdater = new MappingsUpdater();
    _mappingsUpdater.execute();
  }

  /**
   * Iteratively updates mappings in a Swing worker thread. Stops updating
   * when cancelled, say, because constraints have changed. A new mappings
   * updater must be constructed for every update, typically after any of
   * the parameters used to compute mappings have been modified.
   */
  private class MappingsUpdater
    extends SwingWorker<Flattener2C.Mappings,Void>
  {
    public Flattener2C.Mappings doInBackground() { // in Swing worker thread
      trace("MappingsUpdater.doInBackground");
      Flattener2C.Stopper stopper = new Flattener2C.Stopper() {
        public boolean stop() {
          return isCancelled();
        }
      };
      float[][][] cs = getConstraints();
      _mappings = _flattener.updateMappingsFromSlopes(
        _s1,_s2,_p2,_el,cs,_mappings,stopper);
      trace("returning mappings="+_mappings);
      return _mappings;
    }
    public void done() { // in Swing event dispatch thread
      try {
        _mappings = get();
      } catch (Exception e) {
        //throw new RuntimeException(e);
      }
      trace("MappingsUpdater.done: _mappings="+_mappings);
      if (_mappings!=null)
        _contoursView.set(_mappings.u1);
    }
  }

  private ConstraintPoint findNearestConstraintPoint(ConstraintPoint cp) {
    ConstraintPoint pmin = null;
    double dmin = 0.0;
    for (ConstraintPoint p : _csa.cps) {
      double d = cp.distanceTo(p);
      if (pmin == null || d < dmin) {
        pmin = p;
        dmin = d;
      }
    }
    return pmin;
  }

  // One plot for everything.
  private class Plot {

    private Plot(Sampling s1, Sampling s2, float[][] f, float[][] t)
    {
      PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
      pp.setHLabel(0,"Inline (samples)");
      pp.setVLabel("Time (samples)");
      pp.addPixels(0,0,s1,s2,f);
      _contoursView = pp.addContours(0,0,s1,s2,t);
      _contoursView.setContours(40);
      _contoursView.setLineColor(Color.YELLOW);
      //pp.addColorBar(_contoursView,"Relative geologic time");
      for (ConstraintSet cs:_css) {
        pp.addTiledView(0,0,cs.updatePointsView());
      }
      _cssTile = _css[0].pv.getTile();

      // A plot frame has a mode for zooming in tiles or tile axes.
      PlotFrame pf = new PlotFrame(pp);
      TileZoomMode tzm = pf.getTileZoomMode();
      tzm.setActive(true);

      // Add modes for editing constraint sets.
      ModeManager modeManager = pf.getModeManager();
      int ncpm = _css.length;
      ConstraintPointMode[] cpms = new ConstraintPointMode[ncpm];
      for (int icpm=0; icpm<ncpm; ++icpm)
        cpms[icpm] = new ConstraintPointMode(modeManager,_css[icpm]);

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
      pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
      pf.setLocation(50,50);
      pf.setSize(700,800);
      pf.setFontSizeForSlide(1.0,1.0);
      pf.setVisible(true);
    }
  }

  ///////////////////////////////////////////////////////////////////////////

  /**
   * A mode for adding, removing, or moving constraint points.
   */
  private class ConstraintPointMode extends Mode {

    public ConstraintPointMode(ModeManager modeManager, ConstraintSet cs) {
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
            _csa = _cs; // make our constraint set the active one
            component.addMouseListener(_ml);
          } else {
            if (_csa==_cs) // if our constraint set is the active one,
              _csa = null; // then no constraint set will be active
            component.removeMouseListener(_ml);
          }
        }
      }
    }

    private ConstraintSet _cs; // constraint set for this mode
    private ConstraintPoint _cpEdit; // if editing, last point
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
    private ConstraintPoint pointToConstraintPoint(int x, int y) {
      Transcaler ts = _tile.getTranscaler();
      Projector hp = _tile.getHorizontalProjector();
      Projector vp = _tile.getVerticalProjector();
      double xu = ts.x(x);
      double yu = ts.y(y);
      double xv = hp.v(xu);
      double yv = vp.v(yu);
      return new ConstraintPoint(yv,xv);
    }

    // Converts constraint point p to a point (x,y) in pixels.
    private Point constraintPointToPoint(ConstraintPoint cp) {
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
    // fixed number of pixels to the specified constraint point cp.
    private boolean closeEnough(int x, int y, ConstraintPoint cp) {
      Point p = constraintPointToPoint(cp);
      return abs(p.x-x)<12 && abs(p.y-y)<12;
    }

    // Adds a constraint point at mouse coordinates (x,y).
    private void add(MouseEvent e) {
      _tile = (Tile)e.getSource();
      int x = e.getX();
      int y = e.getY();
      ConstraintPoint cp = pointToConstraintPoint(x,y);
      _cs.add(cp);
    }

    // Removes a constraint point, if mouse (x,y) is close enough to one.
    private void remove(MouseEvent e) {
      _tile = (Tile)e.getSource();
      int x = e.getX();
      int y = e.getY();
      ConstraintPoint cp = pointToConstraintPoint(x,y);
      ConstraintPoint cq = findNearestConstraintPoint(cp);
      if (cq!=null && closeEnough(x,y,cq))
        _cs.remove(cq);
    }

    // Begins editing of an existing constraint point, if close enough.
    // Returns true, if close enough so that we have begun editing;
    // false, otherwise.
    private boolean beginEdit(MouseEvent e) {
      _tile = (Tile)e.getSource();
      int x = e.getX();
      int y = e.getY();
      ConstraintPoint cp = pointToConstraintPoint(x,y);
      ConstraintPoint cq = findNearestConstraintPoint(cp);
      if (cq!=null && closeEnough(x,y,cq)) {
        _cs.move(cq,cp);
        _cpEdit = cp;
        return true;
      }
      return false;
    }

    // Called while a constraint point is being dragged during editing.
    private void duringEdit(MouseEvent e) {
      int x = e.getX();
      int y = e.getY();
      ConstraintPoint cp = pointToConstraintPoint(x,y);
      _cs.move(_cpEdit,cp);
      _cpEdit = cp;
    }

    // Called when done editing a constraint point.
    private void endEdit(MouseEvent e) {
      duringEdit(e);
      _editing = false;
    }
  }

  /**
   * Handles program exit.
   */
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

  private static float[][] findFaults(float[][] f) {
    FaultSemblance fse = new FaultSemblance();
    float[][] p2 = fse.slopes(f);
    float[][][] snd = fse.semblanceNumDen(p2,f);
    FaultScanner2 fsc = new FaultScanner2(20,snd,FaultScanner2.Smoother.SHEAR);
    float[][][] flft = fsc.scan(-20,20);
    return flft[0];
  }

  ///////////////////////////////////////////////////////////////////////////

  private Flattener2Demo() {
    int n1 = 251;
    int n2 = 357;
    double d1 = 1.0;
    double d2 = 1.0;
    _s1 = new Sampling(n1,d1,0.0);
    _s2 = new Sampling(n2,d2,0.0);
    String fileName = "/data/seis/tpd/csm/oldslices/tp73.dat";
    float[][] f = new float[n2][n1];
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    double sigma = 8.0; // good for Teapot Dome image tp73
    double pmax = 10.0;
    _p2 = new float[n2][n1];
    _el = new float[n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma,pmax);
    lsf.findSlopes(f,_p2,_el);
    _p2 = mul((float)(d1/d2),_p2);
    _el = findFaults(f);
    _el = pow(sub(1.0f,_el),8);
    _flattener = new Flattener2C();
    _flattener.setWeight1(0.02);
    _flattener.setIterations(0.01,1000);
    _flattener.setSmoothings(4.0,8.0);
    //updateMappings();
    float[][][] cs = getConstraints();
    trace("number of constraint sets = "+cs[0].length);
    _mappings = _flattener.updateMappingsFromSlopes(
      _s1,_s2,_p2,_el,cs,null,null);
    new Plot(_s1,_s2,f,_mappings.u1);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Flattener2Demo();
      }
    });
  }
}
