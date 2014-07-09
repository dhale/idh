
/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.util.*;

import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A linked list of fault cells that may be used to analyze faults.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */

public class FaultSkin {

  /**
   * Gets the cell that was the seed used to grow this skin.
   * @return the seed cell.
   */
  public FaultCell getSeed() {
    return _seed;
  }

  /**
   * Returns the number of cells in this skin.
   */
  public int size() {
    return _cellList.size();
  }

  /**
   * Gets an array of cells in this skin.
   * @return array of cells.
   */
  public FaultCell[] getCells() {
    return _cellList.toArray(new FaultCell[0]);
  }

  /**
   * Gets a list of cells in this skin. 
   * @return list of cells; by reference, not by copy.
   */
  public List<FaultCell> getCellList() {
    return _cellList;
  }

  /**
   * Returns array of arrays of cells linked above and below.
   * @return array of arrays of linked cells; by reference, not by copy.
   */
  public FaultCell[][] getCellsAB() {
    if (_cellsAB!=null)
      return _cellsAB;
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();

    // For all cells in this skin, ...
    for (FaultCell cell:_cellList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search above for the top cell.
        FaultCell c = cell;
        for (FaultCell ca=c.ca; ca!=null; ca=c.ca)
          c = ca;

        // Add the top cell and all cells below it.
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        for (; c!=null; c=c.cb) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the big list.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }

    // Convert the big list to the array to be returned.
    _cellsAB = cellsList.toArray(new FaultCell[0][]);
    return _cellsAB;
  }

  /**
   * Returns array of arrays of cells linked left and right.
   * @return array of arrays of linked cells.
   */
  public FaultCell[][] getCellsLR() {
    if (_cellsLR!=null)
      return _cellsLR;
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();

    // For all cells in this skin, ...
    for (FaultCell cell:_cellList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search left until we have no left nabor or until we return
        // to the cell with which we began, as for a conical fault.
        FaultCell c = cell;
        for (FaultCell cl=c.cl; cl!=null && c!=cell; cl=c.cl)
          c = cl;

        // Remember the leftmost cell found and add it to a new list.
        FaultCell cLeft = c;
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        cList.add(cLeft);
        cellSet.add(cLeft);

        // Add cells to the right. Again beware of cycles.
        for (c=c.cr; c!=null && c!=cLeft; c=c.cr) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the big list.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }

    // Convert the big list to the array to be returned.
    _cellsLR = cellsList.toArray(new FaultCell[0][]);
    return _cellsLR;
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors.
   * In these arrays, cells in this skin are represented by quads with
   * specified size.
   * @param size the size (in samples) of the quads representing cells.
   */
  public float[][] getCellXyzUvwRgb(float size) {
    return FaultCell.getXyzUvwRgb(size,getCells());
  }

  /**
   * Gets arrays of packed cell coordinates for cell links.
   * Each returned array contains packed (x,y,z) cell coordinates for
   * exactly one above-below or left-right linked list of cells.
   * @return array of arrays of packed xyz cell coordinates.
   */
  public float[][] getCellLinksXyz() {
    FaultCell[][] cellsAB = getCellsAB();
    FaultCell[][] cellsLR = getCellsLR();
    int nsAB = cellsAB.length; // number of segments for AB links
    int nsLR = cellsLR.length; // number of segments for LR links
    int ns = nsAB+nsLR; // total number of segments
    float[][] xyz = new float[ns][];
    float[][] rgb = new float[ns][];
    for (int is=0; is<ns; ++is) { // for all segments, ...
      FaultCell[] cells = (is<nsAB)?cellsAB[is]:cellsLR[is-nsAB]; // the cells
      int np = cells.length; // number of points in this segment
      float[] xyzi = new float[3*np]; // xyz for this segment
      for (int ip=0,ic=0; ip<np; ++ip) {
        FaultCell cell = cells[ip];
        xyzi[ic++] = cell.x3;
        xyzi[ic++] = cell.x2;
        xyzi[ic++] = cell.x1;
      }
      xyz[is] = xyzi;
    }
    return xyz;
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  private FaultCell _seed; // cell in this skin with highest fl; null, if empty
  private ArrayList<FaultCell> _cellList; // list of cells in this skin
  private FaultCell[][] _cellsAB; // arrays of cells from above to below
  private FaultCell[][] _cellsLR; // arrays of cells from left to right

  /**
   * Constructs an empty skin.
   */
  FaultSkin() {
    _cellList = new ArrayList<FaultCell>();
  }

  /**
   * Adds the specified cell to this skin.
   * @param cell the cell to be added.
   */
  void add(FaultCell cell) {
    cell.skin = this;
    if (_seed==null)
      _seed = cell;
    _cellList.add(cell);
    _cellsAB = null;
    _cellsLR = null;
  }

  /**
   * Gets a cell nearest the centroid of this skin.
   * In illustrations, this cell is often a good representative.
   * @return the cell nearest the centroid.
   */
  public FaultCell getCellNearestCentroid() {
    float c1 = 0.0f;
    float c2 = 0.0f;
    float c3 = 0.0f;
    float cs = 0.0f;
    for (FaultCell c:_cellList) {
      c1 += c.fl*c.x1;
      c2 += c.fl*c.x2;
      c3 += c.fl*c.x3;
      cs += c.fl;
    }
    c1 /= cs;
    c2 /= cs;
    c3 /= cs;
    float dmin = Float.MAX_VALUE;
    FaultCell cmin = null;
    for (FaultCell c:_cellList) {
      float d = c.distanceSquaredTo(c1,c2,c3);
      if (d<dmin) {
        cmin = c;
        dmin = d;
      }
    }
    return cmin;
  }
}
