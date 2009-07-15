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
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A mode for editing an image by drawing curves.
 * <em>
 * This class is a hack for testing diffusion/interpolation ideas.
 * </em>
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.06
 */
public class ImageEditMode extends Mode {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs an image edit mode with specified manager.
   * @param modeManager the mode manager for this mode.
   */
  public ImageEditMode(
    ModeManager modeManager, PixelsView pixels, float vnull, float[][] v)
  {
    super(modeManager);
    setName("Edit");
    //setIcon(loadIcon(ImageEditMode.class,"resources/ImageEdit16.gif"));
    //setIcon(loadIcon(MouseTrackMode.class,"resources/Track24.gif"));
    setMnemonicKey(KeyEvent.VK_E);
    setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_E,0));
    setShortDescription("Edit points");
    _tile = pixels.getTile();
    _pixels = pixels;
    fill(vnull,v);
    _vnull = vnull;
    _n1 = v[0].length;
    _n2 = v.length;
    _v = v;
    _is2 = new ImageSampler2(_v);
    _ns = 0;
    _x1 = new float[0][0];
    _x2 = new float[0][0];
    _vx = new float[0][0];
    _points = new PointsView(_x1,_x2);
    if (pixels.getOrientation()==PixelsView.Orientation.X1RIGHT_X2UP) {
      _points.setOrientation(PointsView.Orientation.X1RIGHT_X2UP);
    } else {
      _points.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    }
    _points.setStyle("w-o");
  }


  ///////////////////////////////////////////////////////////////////////////
  // protected
  
  protected void setActive(Component component, boolean active) {
    if (component instanceof Tile) {
      Tile tile = (Tile)component;
      if (active) {
        tile.addTiledView(_points);
        tile.addMouseListener(_ml);
        tile.addMouseWheelListener(_mwl);
        InputMap im = tile.getInputMap();
        ActionMap am = tile.getActionMap();
        im.put(KS_BACK_SPACE,"backspace");
        im.put(KS_UP,"up");
        im.put(KS_DOWN,"down");
        am.put("backspace",_bsa);
        am.put("up",_uaa);
        am.put("down",_daa);
      } else {
        tile.removeTiledView(_points);
        tile.removeMouseListener(_ml);
        tile.removeMouseWheelListener(_mwl);
        InputMap im = tile.getInputMap();
        ActionMap am = tile.getActionMap();
        im.remove(KS_BACK_SPACE);
        im.remove(KS_UP);
        im.remove(KS_DOWN);
        am.remove("backspace");
        am.remove("up");
        am.remove("down");
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Tile _tile; // tile containing the pixels view
  private PixelsView _pixels; // pixels view of image to edit
  private PointsView _points; // points view of editing curves
  private int _isSelected = -1; // index of selected segment
  private int _ipSelected = -1; // index of selected point
  private boolean _creating; // true if in process of creating a segment
  float _vnull; // the null (unknown) image value
  private int _n1; // number of samples in 1st dimension of image
  private int _n2; // number of samples in 2nd dimension of image
  float[][] _v; // array[n2][n1] of image values to edit
  private int _ns; // number of segments
  float[][] _x1; // x1 coordinates of points in segments
  float[][] _x2; // x2 coordinates of points in segments
  float[][] _vx; // values v(x1,x2)
  private int _xdown; // x coordinate where mouse down
  private int _ydown; // y coordinate where mouse down
  private int _xedit; // x coordinate of current edit
  private int _yedit; // y coordinate of current edit
  private boolean _hasMotionListener; // true if handling mouse drag
  private ImageSampler2 _is2; // used to get/set image samples

  // Event handlers.
  private static KeyStroke KS_BACK_SPACE = 
    KeyStroke.getKeyStroke(KeyEvent.VK_BACK_SPACE,0);
  private static KeyStroke KS_DOWN = 
    KeyStroke.getKeyStroke(KeyEvent.VK_DOWN,0);
  private static KeyStroke KS_UP = 
    KeyStroke.getKeyStroke(KeyEvent.VK_UP,0);
  private Action _bsa = new AbstractAction() {
    public void actionPerformed(ActionEvent e) {
      onBackSpace();
    }
  };
  private Action _uaa = new AbstractAction() {
    public void actionPerformed(ActionEvent e) {
      onArrow(1.0f);
    }
  };
  private Action _daa = new AbstractAction() {
    public void actionPerformed(ActionEvent e) {
      onArrow(-1.0f);
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
  private MouseWheelListener _mwl = new MouseWheelListener() {;
    public void mouseWheelMoved(MouseWheelEvent e) {
      onMouseWheel(e);
    }
  };
  private MouseMotionListener _mml = new MouseMotionAdapter() {
    public void mouseDragged(MouseEvent e) {
      onMouseDrag(e);
    }
  };

  private void onBackSpace() {
    removeSelectedPoint();
    updateAll();
  }

  private void removeSelectedPoint() {
    int is = _isSelected;
    int ip = _ipSelected;
    if (is>=0 && ip>=0) {
      removePoint(is,ip);
      deselect();
    }
  }

  private void onArrow(float dv) {
    int is = _isSelected;
    int ip = _ipSelected;
    if (is>=0 && ip>=0) {
      _vx[is][ip] += dv;
      trace("value="+_vx[is][ip]);
      updateImage();
      updatePixelsView();
    }
  }

  private void onMouseDown(MouseEvent e) {
    if (_tile!=getTile(e))
      return;

    // Remember mouse down location.
    int x = _xdown = e.getX();
    int y = _ydown = e.getY();

    // What segment and point, if any, was clicked?
    int[] i = {0,0}, p = {0,0};
    int icode = getSegmentAndPoint(x,y,i,p);
    int is = i[0], ip = i[1];
    int xp = p[0], yp = p[1];
    //trace("onMouseDown: is="+is+" ip="+ip);

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

      // If mouse on existing point, ...
      if (icode==1) {

        // If alt key down, delete the point.
        if (e.isAltDown()) {
          selectPoint(is,ip);
          removeSelectedPoint();
        } 
        
        // Else, select the point
        else {
          selectPoint(is,ip);
          addMotionListener();
        }
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
    updateAll();
  }

  private void onMouseDrag(MouseEvent e) {
    changeSelectedPoint(e.getX(),e.getY());
    updateAll();
  }

  private void onMouseUp(MouseEvent e) {
    if (hasMotionListener()) {
      changeSelectedPoint(e.getX(),e.getY());
      removeMotionListener();
    }
    updateAll();
  }

  private void onMouseWheel(MouseWheelEvent e) {
    float dv = (float)(-e.getWheelRotation());
    //trace("dv="+dv);
    onArrow(dv);
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

  private void updateAll() {
    updateImage();
    updatePixelsView();
    updatePointsView();
  }

  private void updateImage() {
    fill(_vnull,_v);
    for (int is=0; is<_ns; ++is) {
      float[] x1 = _x1[is];
      float[] x2 = _x2[is];
      float[] vx = _vx[is];
      int np = _x1[is].length;
      if (np==1) {
        _is2.set(x1[0],x2[0],vx[0]);
      } else if (np>1) {
        for (int ip=1; ip<np; ++ip) {
          float x1a = x1[ip-1];
          float x2a = x2[ip-1];
          float vxa = vx[ip-1];
          float x1b = x1[ip];
          float x2b = x2[ip];
          float vxb = vx[ip];
          ImageSampler2.Line line = _is2.sampleLine(x1a,x2a,x1b,x2b);
          int nr = line.nr();
          float ra = line.ra();
          float rb = line.rb();
          float dvx = (nr>1)?(vxb-vxa)/(float)(nr-1):vxa;
          for (int ir=0; ir<nr; ++ir) {
            float vxi = vxa+(float)ir*dvx;
            line.set(ir,vxi);
          }
        }
      }
    }
  }

  private void updatePixelsView() {
    _pixels.set(_v);
  }

  private void updatePointsView() {
    _points.set(_x1,_x2);
  }

  private int x(float x1, float x2) {
    Transcaler ts = _tile.getTranscaler();
    PointsView.Orientation pvo = _points.getOrientation();
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
    PointsView.Orientation pvo = _points.getOrientation();
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
    PointsView.Orientation pvo = _points.getOrientation();
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
    PointsView.Orientation pvo = _points.getOrientation();
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
    //trace("changePoint: is="+is+" ip="+ip);
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
      float[][] vx = new float[ns][0];
      for (int js=0; js<_ns; ++js) {
        x1[js] = _x1[js];
        x2[js] = _x2[js];
        vx[js] = _vx[js];
      }
      _x1 = x1;
      _x2 = x2;
      _vx = vx;
      _ns = ns;
    }
    int np = _x1[is].length;
    float[] x1p = _x1[is];
    float[] x2p = _x2[is];
    float[] vxp = _vx[is];
    float[] x1t = new float[np+1];
    float[] x2t = new float[np+1];
    float[] vxt = new float[np+1];
    for (int jp=0; jp<ip; ++jp) {
      x1t[jp] = x1p[jp];
      x2t[jp] = x2p[jp];
      vxt[jp] = vxp[jp];
    }
    x1t[ip] = x1(x,y);
    x2t[ip] = x2(x,y);
    vxt[ip] = _is2.get(x1t[ip],x2t[ip]);
    for (int jp=ip; jp<np; ++jp) {
      x1t[jp+1] = x1p[jp];
      x2t[jp+1] = x2p[jp];
      vxt[jp+1] = vxp[jp];
    }
    _x1[is] = x1t;
    _x2[is] = x2t;
    _vx[is] = vxt;
  }

  private void removePoint(int is, int ip) {
    int np = _x1[is].length;
    if (np==1) {
      float[][] x1t = new float[_ns-1][];
      float[][] x2t = new float[_ns-1][];
      float[][] vxt = new float[_ns-1][];
      for (int js=0; js<is; ++js) {
        x1t[js] = _x1[js];
        x2t[js] = _x2[js];
        vxt[js] = _vx[js];
      }
      for (int js=is+1; js<_ns; ++js) {
        x1t[js-1] = _x1[js];
        x2t[js-1] = _x2[js];
        vxt[js-1] = _vx[js];
      }
      --_ns;
      _x1 = x1t;
      _x2 = x2t;
      _vx = vxt;
    } else {
      float[] x1p = _x1[is];
      float[] x2p = _x2[is];
      float[] vxp = _vx[is];
      float[] x1t = new float[np-1];
      float[] x2t = new float[np-1];
      float[] vxt = new float[np-1];
      for (int jp=0; jp<ip; ++jp) {
        x1t[jp] = x1p[jp];
        x2t[jp] = x2p[jp];
        vxt[jp] = vxp[jp];
      }
      for (int jp=ip+1; jp<np; ++jp) {
        x1t[jp-1] = x1p[jp];
        x2t[jp-1] = x2p[jp];
        vxt[jp-1] = vxp[jp];
      }
      _x1[is] = x1t;
      _x2[is] = x2t;
      _vx[is] = vxt;
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
    //float[][] f = sin(rampfloat(0.0f,0.1f,0.1f,n1,n2));
    float[][] f = zerofloat(n1,n2);

    PlotPanel.Orientation orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel panel = new PlotPanel(orientation);

    PixelsView pv = panel.addPixels(f);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setColorModel(ColorMap.JET);

    panel.addColorBar("time");
    panel.setVLabel("depth (km)");
    panel.setHLabel("distance (km)");

    PlotFrame frame = new PlotFrame(panel);
    frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
    frame.setSize(800,700);
    frame.setFontSize(24);
    frame.setVisible(true);

    ModeManager mm = frame.getModeManager();
    ImageEditMode iem = new ImageEditMode(mm,pv,0.0f,f);
    iem.setActive(true);
  } 

  private void trace(String s) {
    System.out.println(s);
  }
}

