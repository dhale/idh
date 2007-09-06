/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

import java.awt.*;
import java.util.EnumSet;

import edu.mines.jtk.mosaic.*;

/**
 * A color bar with a selected range of colors that can be used for painting.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.09.05
 */
public class PaintBar extends ColorBar {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a new paint bar with no label.
   */
  public PaintBar() {
    this(null);
  }

  /**
   * Constructs a new paint bar with specified label.
   * @param label the label; null, if none.
   */
  public PaintBar(String label) {
    super(label);
  }

  /**
   * Sets the range.
   * @param vb begin value of range.
   * @param ve end value of range.
   */
  public void setRange(double vb, double ve) {
    if (!_hasRange || _vb!=vb || _ve!=ve) {
      _vb = (float)vb;
      _ve = (float)ve;
      _hasRange = true;
      updatePointsView();
    }
  }

  public float getRangeBegin() {
    return _vb;
  }

  public float getRangeEnd() {
    return _ve;
  }

  public boolean hasRange() {
    return _hasRange;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  boolean _hasRange = false; // true only after range has been set
  float _vb,_ve; // begin and end of range (begin may exceed end)
  private PointsView _points; // points view that shows the range

  private void updatePointsView() {
    if (!_hasRange)
      return;
    Mosaic mosaic = super.getMosaic();
    Tile tile = mosaic.getTile(0,0);
    Projector hp = tile.getHorizontalProjector();
    Projector vp = tile.getVerticalProjector();
    float h0 = (float)hp.v0();
    float h1 = (float)hp.v1();
    float hm = 0.5f*(h0+h1);
    float[][] x1 = new float[3][2];
    float[][] x2 = new float[3][2];
    x1[0][0] = h0;  x2[0][0] = _vb;
    x1[0][1] = h1;  x2[0][1] = _vb;
    x1[1][0] = h0;  x2[1][0] = _ve;
    x1[1][1] = h1;  x2[1][1] = _ve;
    x1[2][0] = hm;  x2[2][0] = _vb;
    x1[2][1] = hm;  x2[2][1] = _ve;
    if (_points==null) {
      _points = new PointsView(x1,x2);
      tile.addTiledView(_points);
    } else {
      _points.set(x1,x2);
    }
  }
}
