/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

/**
 * A mode for editing points displayed in a PointsView.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.06
 */
public class PointsEditMode extends Mode {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a points edit mode with specified manager.
   * @param modeManager the mode manager for this mode.
   */
  public PointsEditMode(
    ModeManager modeManager, PointsView pv, float[][] x1, float[][] x2) 
  {
    super(modeManager);
    setName("PointsEdit");
    //setIcon(loadIcon(PointsEditMode.class,"resources/PointsEdit16.gif"));
    setMnemonicKey(KeyEvent.VK_E);
    setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_E,0));
    setShortDescription("Edit points");
    _pv = pv;
    _x1 = x1;
    _x2 = x2;
    _ns = x1.length;
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected
  
  protected void setActive(Component component, boolean active) {
    if (component instanceof Tile) {
      Tile tile = (Tile)component;
      InputMap im = tile.getInputMap();
      ActionMap am = tile.getActionMap();
      if (active) {
        tile.addMouseListener(_ml);
        im.put(KS_BACK_SPACE,"backspace");
        am.put("backspace",_bsa);
      } else {
        tile.removeMouseListener(_ml);
        im.remove(KS_BACK_SPACE);
        am.remove("backspace");
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static KeyStroke KS_BACK_SPACE = 
    KeyStroke.getKeyStroke(KeyEvent.VK_BACK_SPACE,0);

  private Tile _tile; // tile in which editing began; null, if not editing
  private PointsView _pv; // points view to update as points are edited
  private int _isSelected = -1; // index of selected segment
  private int _ipSelected = -1; // index of selected point
  private boolean _creating; // true if in process of creating a segment
  private int _ns; // number of segments.
  float[][] _x1; // x1 coordinates of points in segments
  float[][] _x2; // x2 coordinates of points in segments
  private int _xdown; // x coordinate where mouse down
  private int _ydown; // y coordinate where mouse down
  private int _xedit; // x coordinate of current edit
  private int _yedit; // y coordinate of current edit
  private boolean _hasMotionListener; // true if handling mouse drag

  // Event handlers.
  private Action _bsa = new AbstractAction() {
    public void actionPerformed(ActionEvent e) {
      onBackSpace();
    }
  };
  private MouseListener _ml = new MouseAdapter() {;
    public void mousePressed(MouseEvent e) {
      onMouseDown(e);
    }
    public void mouseReleased(MouseEvent e) {
      onMouseUp(e);
    }
  };
  private MouseMotionListener _mml = new MouseMotionAdapter() {
    public void mouseDragged(MouseEvent e) {
      onMouseDrag(e);
    }
  };

  private void onBackSpace() {
    trace("onBackSpace:");
    int is = _isSelected;
    int ip = _ipSelected;
    if (is>=0 && ip>=0) {
      trace("  removing point: is="+is+" ip="+ip);
      removePoint(is,ip);
      deselect();
      updatePointsView();
    }
  }

  private void onMouseDown(MouseEvent e) {
    if (!hasPointsView(getTile(e)))
      return;
    _tile = getTile(e);

    // Remember mouse down location.
    int x = _xdown = e.getX();
    int y = _ydown = e.getY();

    // What segment and point, if any, was clicked?
    int[] i = {0,0}, p = {0,0};
    int icode = getSegmentAndPoint(x,y,i,p);
    int is = i[0], ip = i[1];
    int xp = p[0], yp = p[1];
    trace("onMouseDown: is="+is+" ip="+ip);

    // If currently in the process of creating a new segment, ...
    if (_creating) {

      // If mouse is on previous (selected) point in the created segment,
      // then end segment creation.
      if (icode==1 && isSelected(is,ip)) {
        _creating = false;
      } 
      
      // Else append a new point to the created segment.
      else {
        is = _isSelected;
        ip = appendPoint(x,y);
        selectPoint(is,ip);
        addMotionListener();
      }
    }

    // Else if not currently in the process of creating a new segment, ...
    else {

      // If mouse on existing point, select the point
      if (icode==1) {
        selectPoint(is,ip);
        addMotionListener();
      }

      // Else if mouse on existing segment, ...
      else if (icode==2) {

        // If shift key down, add a new point to existing segment.
        if (e.isShiftDown()) {
          insertPoint(is,ip,x,y);
          selectPoint(is,ip);
          addMotionListener();
        }
      }

      // Else if mouse down away from all segments, ...
      else {

        // If shift key down, create a new segment.
        if (e.isShiftDown()) {
          is = appendSegment(x,y);
          ip = 0;
          selectPoint(is,ip);
          addMotionListener();
          _creating = true;
        }

        // Else deselect all points.
        else {
          deselect();
        }
      }
    }
    updatePointsView();
  }

  private void onMouseDrag(MouseEvent e) {
    changeSelectedPoint(e.getX(),e.getY());
    updatePointsView();
  }

  private void onMouseUp(MouseEvent e) {
    if (hasMotionListener()) {
      changeSelectedPoint(e.getX(),e.getY());
      removeMotionListener();
    }
    _tile = null;
    updatePointsView();
  }

  private void addMotionListener() {
    _tile.addMouseMotionListener(_mml);
    _hasMotionListener = true;
  }

  private void removeMotionListener() {
    _tile.removeMouseMotionListener(_mml);
    _hasMotionListener = false;
  }

  private boolean hasMotionListener() {
    return _hasMotionListener;
  }

  private Tile getTile(MouseEvent e) {
    return (Tile)e.getSource();
  }

  private boolean hasPointsView(Tile tile) {
    int ntv = tile.countTiledViews();
    for (int itv=ntv-1; itv>=0; --itv) {
      TiledView tv = tile.getTiledView(itv);
      if (tv instanceof PointsView) {
        PointsView pv = (PointsView)tv;
        if (_pv==pv)
          return true;
      }
    }
    return false;
  }

  private void updatePointsView() {
    _pv.set(_x1,_x2);
  }

  private int x(float x1, float x2) {
    Transcaler ts = _tile.getTranscaler();
    PointsView.Orientation pvo = _pv.getOrientation();
    if (pvo==PointsView.Orientation.X1RIGHT_X2UP) {
      Projector p = _tile.getHorizontalProjector();
      return ts.x(p.u(x1));
    } else {
      Projector p = _tile.getVerticalProjector();
      return ts.x(p.u(x2));
    }
  }

  private int y(float x1, float x2) {
    Transcaler ts = _tile.getTranscaler();
    PointsView.Orientation pvo = _pv.getOrientation();
    if (pvo==PointsView.Orientation.X1RIGHT_X2UP) {
      Projector p = _tile.getHorizontalProjector();
      return ts.y(p.u(x2));
    } else {
      Projector p = _tile.getVerticalProjector();
      return ts.y(p.u(x1));
    }
  }

  private float x1(int x, int y) {
    Transcaler ts = _tile.getTranscaler();
    PointsView.Orientation pvo = _pv.getOrientation();
    if (pvo==PointsView.Orientation.X1RIGHT_X2UP) {
      Projector p = _tile.getHorizontalProjector();
      return (float)p.v(ts.x(x));
    } else {
      Projector p = _tile.getVerticalProjector();
      return (float)p.v(ts.y(y));
    }
  }

  private float x2(int x, int y) {
    Transcaler ts = _tile.getTranscaler();
    PointsView.Orientation pvo = _pv.getOrientation();
    if (pvo==PointsView.Orientation.X1RIGHT_X2UP) {
      Projector p = _tile.getHorizontalProjector();
      return (float)p.v(ts.y(y));
    } else {
      Projector p = _tile.getVerticalProjector();
      return (float)p.v(ts.x(x));
    }
  }

  private static double pointToPointSquared(
    double xp, double yp,
    double xa, double ya)
  {
    double xq = xp-xa;
    double yq = yp-ya;
    return xq*xq+yq*yq;
  }

  private static double pointToSegmentSquared(
    double xp, double yp,
    double xa, double ya,
    double xb, double yb,
    double[] q)
  {
    double xq = xp-xa;
    double yq = yp-ya;
    double xr = xb-xa;
    double yr = yb-ya;
    double rn = xq*xr+yq*yr;
    if (rn<=0.0) {
      xq = xa;
      yq = ya;
    } else {
      double rd = xr*xr+yr*yr;
      if (rd<=rn) {
        xq = xb;
        yq = yb;
      } else {
        double r = rn/rd;
        xq = xa+r*xr;
        yq = ya+r*yr;
      }
    }
    if (q!=null) {
      q[0] = xq;
      q[1] = yq;
    }
    xq -= xp;
    yq -= yp;
    return xq*xq+yq*yq;
  }

  private void changeSelectedPoint(int x, int y) {
    changePoint(_isSelected,_ipSelected,x,y);
  }

  private void changePoint(int is, int ip, int x, int y) {
    trace("changePoint: is="+is+" ip="+ip);
    _x1[is][ip] = x1(x,y);
    _x2[is][ip] = x2(x,y);
  }

  private int appendPoint(int is, int x, int y) {
    int ip = _x1[is].length;
    insertPoint(is,ip,x,y);
    return ip;
  }

  private int appendPoint(int x, int y) {
    return appendPoint(_ns-1,_xdown,_ydown);
  }

  private int appendSegment(int x, int y) {
    insertPoint(_ns,0,x,y);
    return _ns-1;
  }

  private void insertPoint(int is, int ip, int x, int y) {
    Check.argument(is<=_ns,"is in bounds");
    Check.argument(
      is==_ns && ip==0 || is<_ns && ip<=_x1[is].length,
      "is in bounds");
    if (is==_ns) {
      int ns = is+1;
      float[][] x1 = new float[ns][0];
      float[][] x2 = new float[ns][0];
      for (int js=0; js<_ns; ++js) {
        x1[js] = _x1[js];
        x2[js] = _x2[js];
      }
      _x1 = x1;
      _x2 = x2;
      _ns = ns;
    }
    int np = _x1[is].length;
    float[] x1p = _x1[is];
    float[] x2p = _x2[is];
    float[] x1t = new float[np+1];
    float[] x2t = new float[np+1];
    for (int jp=0; jp<ip; ++jp) {
      x1t[jp] = x1p[jp];
      x2t[jp] = x2p[jp];
    }
    x1t[ip] = x1(x,y);
    x2t[ip] = x2(x,y);
    for (int jp=ip; jp<np; ++jp) {
      x1t[jp+1] = x1p[jp];
      x2t[jp+1] = x2p[jp];
    }
    _x1[is] = x1t;
    _x2[is] = x2t;
  }

  private void removePoint(int is, int ip) {
    int np = _x1[is].length;
    if (np==1) {
      float[][] x1t = new float[_ns-1][];
      float[][] x2t = new float[_ns-1][];
      for (int js=0; js<is; ++js) {
        x1t[js] = _x1[js];
        x2t[js] = _x2[js];
      }
      for (int js=is+1; js<_ns; ++js) {
        x1t[js-1] = _x1[js];
        x2t[js-1] = _x2[js];
      }
      --_ns;
      _x1 = x1t;
      _x2 = x2t;
    } else {
      float[] x1p = _x1[is];
      float[] x2p = _x2[is];
      float[] x1t = new float[np-1];
      float[] x2t = new float[np-1];
      for (int jp=0; jp<ip; ++jp) {
        x1t[jp] = x1p[jp];
        x2t[jp] = x2p[jp];
      }
      for (int jp=ip+1; jp<np; ++jp) {
        x1t[jp-1] = x1p[jp];
        x2t[jp-1] = x2p[jp];
      }
      _x1[is] = x1t;
      _x2[is] = x2t;
    }
  }

  private boolean isSelected(int is, int ip) {
    return _isSelected==is && _ipSelected==ip;
  }

  private void selectPoint(int is, int ip) {
    _isSelected = is;
    _ipSelected = ip;
  }

  private void deselect() {
    _isSelected = -1;
    _ipSelected = -1;
  }

  /**
   * Gets indices of a nearby segment and point, if close enough.
   * If an existing segment is close enough, both segment and point indices
   * {is,ip} are returned. If an existing point on the segment is close
   * enough, then ip is the index of that existing point. Otherwise, if a
   * line between two existing points is close enough, then ip is the 
   * larger index of those two existing points. In both cases, the pixel 
   * coordinates {x,y} of the closest point are returned. If no existing 
   * segment is close enough, then {is,ip} = {-1,-1}.
   * @param x the x coordinate in pixels.
   * @param y the y coordinate in pixels.
   * @param i array {is,ip} of indices of segment and point; {-1,-1} if none.
   * @param p array {xp,yp} coordinates of close point. This point could be an
   *          existing point or a point on a line between two existing points.
   * @return 0, if no existing segment (or point) is close enough; 
   *         1, if the closest point already exists in a segment; 
   *         2, if the closest point lies on a segment between two 
   *            existing points in that segment.
   */
  private int getSegmentAndPoint(int x, int y, int[] i, int[] p) {
    double[] pt = {0.0,0.0};
    double dpmin = 25.0; // distance must be less than 5 pixels
    double dsmin = 25.0; // 25 is this distance squared
    int icode = 0; // initially assume no segment is close enough
    int ipmin = -1;
    int ismin = -1;
    int xtmin = -1;
    int ytmin = -1;
    for (int is=0; is<_ns; ++is) {
      double xa = 0.0;
      double ya = 0.0;
      double xb = 0.0;
      double yb = 0.0;
      int np = _x1[is].length;
      for (int ip=0; ip<np; ++ip) {
        xb = xa;
        yb = ya;
        xa = x(_x1[is][ip],_x2[is][ip]);
        ya = y(_x1[is][ip],_x2[is][ip]);
        double dp = pointToPointSquared(x,y,xa,ya);
        if (dp<dpmin) {
          ipmin = ip;
          ismin = is;
          dpmin = dp;
          xtmin = (int)(xa+0.5);
          ytmin = (int)(ya+0.5);
          icode = 1;
        } else if (ip>0 && icode!=1) {
          double ds = pointToSegmentSquared(x,y,xa,ya,xb,yb,pt);
          if (ds<dsmin) {
            ipmin = ip;
            ismin = is;
            dsmin = ds;
            xtmin = (int)(pt[0]+0.5);
            ytmin = (int)(pt[1]+0.5);
            icode = 2;
          }
        }
      }
    }
    i[0] = ismin;
    i[1] = ipmin;
    p[0] = xtmin;
    p[1] = ytmin;
    return icode;
  }

  private void drawLine(JComponent c, int x, int y) {
    _xedit = x;
    _yedit = y;
    Graphics g = c.getGraphics();
    g.setColor(Color.BLACK);
    g.setXORMode(c.getBackground());
    g.drawLine(_xdown,_ydown,_xedit,_yedit);
    g.dispose();
  }

  /*
  To create a new segment, 
    Begin by shift-clicking away from any existing segment
    Continue clicking away to add points to segment
    Click again on last point added to end segment
    While new segment not yet ended,
      the last point added is selected
    When new segment has been ended,
      no point is selected
  To select a point in any segment,
    Click the point
  To unselect all points,
    Click away from any points
  To modify an existing segment,
    Click-drag a point to move it
    Shift-click to add a new point
    Select an existing point and press backspace key

  onMouseDown
    if creating
      if mouse on selected point (last point created)
        creating = false
      else
        append and select point
        add motion listener
    else 
      locate mouse
      if mouse on existing point
        select existing point
        add motion listener
        dragging = true
      else if mouse on existing segment
        if shiftDown
          insert and select new point in existing segment
          add motion listener
      else // if mouse is away from any segments
        if shiftDown
          append and select point in a new segment
          add motion listener
          creating = true
        else
          deselect all points
  onMouseDrag (should be called only when mouse is down and dragged)
    update coordinates of selected point
  onMouseUp
    if motion listener
      update coordinates of selected point
      remove motion listener
  */

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        go();
      }
    });
  }
  private static void go() {
    int n1 = 101;
    int n2 = 101;
    float[][] f = Array.sin(Array.rampfloat(0.0f,0.1f,0.1f,n1,n2));
    float[][] x1 = new float[0][0];
    float[][] x2 = new float[0][0];

    PlotPanel.Orientation orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel panel = new PlotPanel(orientation);

    PixelsView pxv = panel.addPixels(f);
    pxv.setColorModel(ColorMap.JET);

    PointsView ptv = panel.addPoints(x1,x2);
    ptv.setStyle("k-o");

    panel.addColorBar("velocity");
    panel.setVLabel("depth (km)");
    panel.setHLabel("distance (km)");

    PlotFrame frame = new PlotFrame(panel);
    frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
    frame.setSize(800,700);
    frame.setFontSize(24);
    frame.setVisible(true);

    ModeManager mm = frame.getModeManager();
    PointsEditMode pem = new PointsEditMode(mm,ptv,x1,x2);
    pem.setActive(true);
  } 

  private void trace(String s) {
    System.out.println(s);
  }
}

