/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.awt.*;
import java.awt.event.*;

import javax.swing.JComponent;
import javax.swing.KeyStroke;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * A mode for painting a 2-D array of floats.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.06
 */
public class PaintMode extends Mode {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a paint mode with specified manager.
   * @param modeManager the mode manager for this mode.
   */
  public PaintMode(
    ModeManager modeManager, PaintBar pb, PixelsView pv, float[][] f) 
  {
    super(modeManager);
    setName("Paint");
    //setIcon(loadIcon(PaintMode.class,"resources/Paint16.gif"));
    setMnemonicKey(KeyEvent.VK_P);
    setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_P,0));
    setShortDescription("Paint");
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected
  
  protected void setActive(Component component, boolean active) {
    if (component instanceof Tile) {
      if (active) {
        component.addMouseListener(_ml);
      } else {
        component.removeMouseListener(_ml);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Tile _tile; // tile in which painting began
  private int _xbegin; // x coordinate where paint began
  private int _ybegin; // y coordinate where paint began
  private int _xdraw; // x coordinate to which line was last drawn
  private int _ydraw; // y coordinate to which line was last drawn

  private MouseListener _ml = new MouseAdapter() {;
    public void mousePressed(MouseEvent e) {
      beginPaint(e);
    }
    public void mouseReleased(MouseEvent e) {
      endPaint(e);
    }
  };

  private MouseMotionListener _mml = new MouseMotionAdapter() {
    public void mouseDragged(MouseEvent e) {
      duringPaint(e);
    }
  };

  private void beginPaint(MouseEvent e) {
    _xbegin = e.getX();
    _ybegin = e.getY();
    Object source = e.getSource();
    Tile tile = _tile = (Tile)source;
    drawPaint(tile,_xbegin,_ybegin);
    tile.addMouseMotionListener(_mml);
  }

  private void duringPaint(MouseEvent e) {
    int xdraw = e.getX();
    int ydraw = e.getY();
    drawPaint(_tile,_xdraw,_ydraw);
    drawPaint(_tile, xdraw, ydraw);
  }

  private void endPaint(MouseEvent e) {

    drawPaint(_tile,_xdraw,_ydraw);
    _tile.removeMouseMotionListener(_mml);
    Transcaler ts = tile.getTranscaler();
    Projector hp = tile.getHorizontalProjector();
    Projector vp = tile.getVerticalProjector();
    float xa = hp.v(ts.x(_xbegin));
    float ya = vp.v(ts.y(_ybegin));
    float xb = hp.v(ts.x(_xdraw));
    float yb = vp.v(ts.y(_ydraw));

    PixelsView.Orientation pvo = _pixels.getOrientation();
    float x1a,x2a,x1b,x2b;
    if (pvo==PixelsView.Orientation.X1RIGHT_X2UP) {
      x1a = xa;
      x1b = xb;
      x2a = ya;
      x2b = yb;
    } else {
      x1a = ya;
      x1b = yb;
      x2a = xa;
      x2b = xb;
    }
    ImageSampler2 is = new ImageSampler2(

    // If painting or unpainting x, ...
    if (zx) {
      vr.x = (xmin<xmax)?ts.x(xmin):0.0;
      vr.width = (xmin<xmax)?ts.x(xmax)-vr.x:1.0;
    }

    // If painting or unpainting y, ...
    if (zy) {
      vr.y = (ymin<ymax)?ts.y(ymin):0.0;
      vr.height = (ymin<ymax)?ts.y(ymax)-vr.y:1.0;
    }

    // Arbitrarily limit paint to a factor of 10000.
    double tiny = 0.0001;
    if (vr.width<tiny) {
      vr.x -= (tiny-vr.width)/2;
      vr.width = tiny;
    }
    if (vr.height<tiny) {
      vr.y -= (tiny-vr.height)/2;
      vr.height = tiny;
    }
    vr.x = min(1.0-vr.width,vr.x);
    vr.y = min(1.0-vr.height,vr.y);

    // Set view rectangle of one tile, and mosaic will set the others.
    tile.setViewRectangle(vr);
    _tile = null;
  }

  private void drawPaint(Tile tile, int x, int y) {
    if (tile==null)
      return;

    // Clip paint to tile bounds.
    x = max(0,min(tile.getWidth()-1,x));
    y = max(0,min(tile.getHeight()-1,y));

    // Draw line in this tile.
    drawLine(tile,x,y);
  }

  private void drawLine(JComponent c, int x, int y) {
    _xdraw = x;
    _ydraw = y;
    Graphics g = c.getGraphics();
    g.setColor(Color.BLACK);
    g.setXORMode(c.getBackground());
    g.drawLine(_xbegin,_ybegin,_xdraw,_ydraw);
    g.dispose();
  }
}

