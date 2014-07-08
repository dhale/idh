
/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.util.ArrayList;
import java.util.HashSet;

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
   * Gets an array of cells in this skin.
   * @return array of cells.
   */
  public FaultCell[] getCells() {
    return cellList.toArray(new FaultCell[0]);
  }

  /**
   * Gets the cell that was the seed used to grow this skin.
   * @return the seed cell.
   */
  public FaultCell getSeed() {
    return seed;
  }

  /**
   * Returns the number of cells in this skin.
   */
  public int size() {
    return cellList.size();
  }

  /**
   * Returns array of arrays of cells linked above and below.
   * @return array of arrays of linked cells.
   */
  public FaultCell[][] getCellsAB() {
    if (cellsAB!=null)
      return cellsAB;
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();
    for (FaultCell cell:cellList) {
      if (!cellSet.contains(cell)) {
        FaultCell c = cell;
        for (FaultCell ca=c.ca; ca!=null; ca=c.ca)
          c = ca;
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        for (; c!=null; c=c.cb) {
          cList.add(c);
          cellSet.add(c);
        }
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    cellsAB = cellsList.toArray(new FaultCell[0][]);
    return cellsAB;
  }

  /**
   * Returns array of arrays of cells linked left and right.
   * @return array of arrays of linked cells.
   */
  public FaultCell[][] getCellsLR() {
    if (cellsLR!=null)
      return cellsLR;
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();
    for (FaultCell cell:cellList) {
      if (!cellSet.contains(cell)) {
        FaultCell c = cell;
        for (FaultCell cl=c.cl; cl!=null; cl=c.cl)
          c = cl;
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        for (; c!=null; c=c.cr) {
          cList.add(c);
          cellSet.add(c);
        }
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    cellsLR = cellsList.toArray(new FaultCell[0][]);
    return cellsLR;
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

  FaultCell seed; // the cell in this skin with highest fl; null, if empty
  ArrayList<FaultCell> cellList; // list of cells in this skin
  FaultCell[][] cellsAB; // arrays of arrays of cells from above to below
  FaultCell[][] cellsLR; // arrays of arrays of cells from left to right

  /**
   * Constructs an empty skin.
   */
  FaultSkin() {
    cellList = new ArrayList<FaultCell>();
  }

  /**
   * Adds the specified cell to this skin.
   * @param cell the cell to be added.
   */
  void add(FaultCell cell) {
    if (this.seed==null)
      this.seed = cell;
    cellList.add(cell);
    cell.skin = this;
    cellsAB = null;
    cellsLR = null;
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
    for (FaultCell c:cellList) {
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
    for (FaultCell c:cellList) {
      float d = c.distanceSquaredTo(c1,c2,c3);
      if (d<dmin) {
        cmin = c;
        dmin = d;
      }
    }
    return cmin;
  }
}
