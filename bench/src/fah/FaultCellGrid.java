/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;
import static fah.FaultGeometry.*;

/**
 * Fault cells in a 3D sampling grid. Each grid sample indexed by (i1,i2,i3)
 * contains either one fault cell or null. The grid facilitates searches for
 * cell nabors in skins and fast iterations along fault traces tangent to
 * fault strike and fault curves tangent to fault dip.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.06
 */
public class FaultCellGrid {

  /**
   * Constructs an empty fault grid with specified dimensions.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   */
  public FaultCellGrid(int n1, int n2, int n3) {
    this(n1,n2,n3,null);
  }

  /**
   * Constructs a fault grid with specified dimensions and cells.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param cells array of cells to be included in the grid.
   */
  public FaultCellGrid(int n1, int n2, int n3, FaultCell[] cells) {
    _n1 = n1; _n2 = n2; _n3 = n3;
    _cells = new FaultCell[n3][n2][n1];
    if (cells!=null) {
      for (FaultCell cell:cells) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        _cells[i3][i2][i1] = cell;
      }
    }
  }

  /**
   * Gets the fault cell with specified indices, if any.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @return the fault cell; null, if none or if indices are out of bounds.
   */
  public FaultCell get(int i1, int i2, int i3) {
    if (0<=i1 && i1<_n1 && 0<=i2 && i2<_n2 && 0<=i3 && i3<_n3) {
      return _cells[i3][i2][i1];
    } else {
      return null;
    }
  }

  /**
   * Sets the fault cell with specified indices.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @param cell the fault cell.
   */
  public void set(int i1, int i2, int i3, FaultCell cell) {
    _cells[i3][i2][i1] = cell;
  }

  /**
   * Finds a fault cell above the specified cell. Searches for a cell above
   * that lies nearest to the line containing the specified cell and its dip
   * vector. If the specified cell is already linked to a nabor cell above,
   * this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell above.
   * @return the cell above; null, if none.
   */
  public FaultCell findCellAbove(FaultCell cell) {
    if (cell==null) return null;
    if (cell.ca!=null) return cell.ca;
    return findCellAboveBelow(true,cell);
  }

  /**
   * Finds a fault cell below the specified cell. Searches for a cell below
   * that lies nearest to the line containing the specified cell and its dip
   * vector. If the specified cell is already linked to a nabor cell below,
   * this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell below.
   * @return the cell below; null, if none.
   */
  public FaultCell findCellBelow(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cb!=null) return cell.cb;
    return findCellAboveBelow(false,cell);
  }

  /**
   * Finds a fault cell left of the specified cell. Searches for a cell left
   * that lies nearest to the line containing the specified cell and its
   * strike vector. If the specified cell is already linked to a nabor cell
   * left, this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell left.
   * @return the cell left; null, if none.
   */
  public FaultCell findCellLeft(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cl!=null) return cell.cl;
    return findCellLeftRight(true,cell);
  }

  /**
   * Finds a fault cell right of the specified cell. Searches for a cell right
   * that lies nearest to the line containing the specified cell and its
   * strike vector. If the specified cell is already linked to a nabor cell
   * right, this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell right.
   * @return the cell right; null, if none.
   */
  public FaultCell findCellRight(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cr!=null) return cell.cr;
    return findCellLeftRight(false,cell);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _n1,_n2,_n3;
  private FaultCell[][][] _cells;

  private FaultCell findCellAboveBelow(boolean above, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float u1 = cell.u1;
    float u2 = cell.u2;
    float u3 = cell.u3;
    int k1 = 1;
    if (above) {
      k1 = -k1;
      u1 = -u1;
      u2 = -u2;
      u3 = -u3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        FaultCell c = get(i1+k1,i2+k2,i3+k3);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float d3 = c.x3-x3;
          float du = d1*u1+d2*u2+d3*u3;
          if (du>0.0f) {
            d1 -= du*u1;
            d2 -= du*u2;
            d3 -= du*u3;
            float d = d1*d1+d2*d2+d3*d3; // squared distance to dip line
            if (d<dmin) {
              cmin = c;
              dmin = d;
            }
          }
        }
      }
    }
    return cmin;
  }

  // The search for a cell left or right is not so straightforward as for a
  // cell above or below. The specified cell has eight adjacent samples. We
  // want a cell that is both nearby and located in the strike direction (if
  // right) or opposite direction (if left). We therefore first look for the
  // best cell among the N, E, S, and W adjacent samples, because they are
  // likely to be nearest, and we do not want to skip over them. If and only
  // if we do not find any candidate cells located in the specified direction,
  // we then look among the NE, SE, SW, and NW samples.
  private static final int[] K2LR = { 0, 1, 0,-1, 1, 1,-1,-1};
  private static final int[] K3LR = { 1, 0,-1, 0, 1,-1,-1, 1};
  private FaultCell findCellLeftRight(boolean left, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float v1 = cell.v1;
    float v2 = cell.v2;
    float v3 = cell.v3;
    if (left) {
      v1 = -v1;
      v2 = -v2;
      v3 = -v3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int ik=0; ik<8; ++ik) {
      if (ik==4 && cmin!=null)
        break;
      int k2 = K2LR[ik];
      int k3 = K3LR[ik];
      FaultCell c = get(i1,i2+k2,i3+k3);
      if (c!=null) {
        float d1 = c.x1-x1;
        float d2 = c.x2-x2;
        float d3 = c.x3-x3;
        float dv = d1*v1+d2*v2+d3*v3;
        if (dv>0.0f) {
          d1 -= dv*v1;
          d2 -= dv*v2;
          d3 -= dv*v3;
          float d = d1*d1+d2*d2+d3*d3; // squared distance to strike line
          if (d<dmin) {
            cmin = c;
            dmin = d;
          }
        }
      }
    }
    return cmin;
  }
}
