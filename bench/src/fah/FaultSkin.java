
/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import java.io.*;
import java.util.*;

import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A linked list of fault cells that may be used to analyze faults.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */

public class FaultSkin implements Iterable<FaultCell>,Serializable {
  private static final long serialVersionUID = 1L;

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
   * Returns an iterator for the cells in this skin. 
   * @return cell iterator.
   */
  public Iterator<FaultCell> iterator() {
    return _cellList.iterator();
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

        // Add the top cell and all cells below it to a new list.
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        for (; c!=null; c=c.cb) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    assert _cellList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _cellsAB = cellsList.toArray(new FaultCell[0][]);
    return _cellsAB;
  }

  /**
   * Returns array of arrays of cells linked left and right.
   * @return array of arrays of linked cells; by reference, not by copy.
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

        // Search left until we find no left nabor or until that left
        // nabor is the cell with which we began the search.
        FaultCell c = cell;
        for (FaultCell cl=c.cl; cl!=null && cl!=cell; cl=c.cl)
          c = cl;

        // Remember the leftmost cell found and add it to a new list.
        FaultCell cLeft = c;
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        cList.add(c);
        cellSet.add(c);

        // Add cells found to the right. Again beware of cycles.
        for (c=c.cr; c!=null && c!=cLeft; c=c.cr) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    assert _cellList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _cellsLR = cellsList.toArray(new FaultCell[0][]);
    checkCellArrays(_cellsLR);
    return _cellsLR;
  }

  /**
   * Smooths the normal vectors of cells in this skin.
   * @param nsmooth the number of smoothings.
   */
  public void smoothCellNormals(int nsmooth) {
    FaultCell.GetN getter = new FaultCell.GetN() {
      public float[] get(FaultCell cell) {
        return new float[]{cell.w1,cell.w2,cell.w3};
      }
    };
    FaultCell.SetN setter = new FaultCell.SetN() {
      public void set(FaultCell cell, float[] w) {
        float w1 = w[0]; 
        float w2 = w[1]; 
        float w3 = w[2];
        float ws = 1.0f/sqrt(w1*w1+w2*w2+w3*w3);
        w1 *= ws;
        w2 *= ws;
        w3 *= ws;
        cell.setNormalVector(w1,w2,w3);
      }
    };
    for (int ismooth=0; ismooth<nsmooth; ++ismooth)
      smoothN(getter,setter);
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

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, cells in this skin are represented by quads with specified
   * size, and colors corresponding to fault likelihoods.
   * @param size size (in samples) of the quads.
   * @param cmap colormap used to compute rgb colors from cell properties.
   */
  public float[][] getCellXyzUvwRgbForLikelihood(float size, ColorMap cmap) {
    return FaultCell.getXyzUvwRgbForLikelihood(size,cmap,getCells());
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, cells in this skin are represented by quads with specified
   * size, and colors corresponding to fault throws.
   * @param size size (in samples) of the quads.
   * @param cmap colormap used to compute rgb colors from cell properties.
   */
  public float[][] getCellXyzUvwRgbForThrow(float size, ColorMap cmap) {
    return FaultCell.getXyzUvwRgbForThrow(size,cmap,getCells());
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
      int ncell = cells.length; // number of cells in this segment
      int np = ncell; // number of points in this segment
      if (is>=nsAB && cells[0].cl==cells[ncell-1]) // if a LR cycle, ...
        ++np; // then add one more so we end with the starting point
      float[] xyzi = new float[3*np]; // xyz for this segment
      for (int ip=0,ic=0; ip<np; ++ip) {
        FaultCell cell = cells[ip%ncell]; // % to handle any LR cycle
        xyzi[ic++] = cell.x3;
        xyzi[ic++] = cell.x2;
        xyzi[ic++] = cell.x1;
      }
      xyz[is] = xyzi;
    }
    return xyz;
  }

  public static FaultSkin readFromFile(String fileName) {
    try {
      FileInputStream fis = new FileInputStream(fileName);
      ObjectInputStream ois = new ObjectInputStream(fis);
      FaultSkin skin = (FaultSkin)ois.readObject();
      ois.close();
      return skin;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static void writeToFile(String fileName, FaultSkin skin) {
    try {
      FileOutputStream fos = new FileOutputStream(fileName);
      ObjectOutputStream oos = new ObjectOutputStream(fos);
      oos.writeObject(skin);
      oos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  /**
   * Constructs an empty skin.
   */
  FaultSkin() {
    _cellList = new ArrayList<FaultCell>();
  }

  /**
   * Adds the specified skinless cell to this skin.
   * @param cell the cell to be added.
   */
  void add(FaultCell cell) {
    assert cell.skin==null;
    cell.skin = this;
    if (_seed==null)
      _seed = cell;
    _cellList.add(cell);
    _cellsAB = null;
    _cellsLR = null;
  }

  /**
   * Smooths one value stored in the cells of this skin. The value smoothed is
   * that accessed by the specified getter and setter. Each smoothed value is
   * an average of the values in a cell and its cell nabors. 
   */
  void smooth1(FaultCell.Get1 getter, FaultCell.Set1 setter) {
    int ncell = size();
    float[] vals = new float[ncell];
    float[] cnts = new float[ncell];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float valCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float valNabor = getter.get(cellNabor);
          vals[icell] += valCell+valNabor;
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float cnti = cnts[icell];
      float vali = vals[icell]/(cnti>0.0f?cnti:1.0f);
      setter.set(cell,vali);
    }
  }

  /**
   * Smooths multiple values stored in the cells of this skin. The values
   * smoothed are those accessed by the specified getter and setter. Each
   * smoothed value is an average of the values in a cell and its cell nabors. 
   */
  void smoothN(FaultCell.GetN getter, FaultCell.SetN setter) {
    int ncell = size();
    int nval = getter.get(_seed).length;
    float[][] vals = new float[ncell][nval];
    float[] cnts = new float[ncell];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float[] valsCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float[] valsNabor = getter.get(cellNabor);
          for (int ival=0; ival<nval; ++ival)
            vals[icell][ival] += valsCell[ival]+valsNabor[ival];
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float cnti = cnts[icell];
      float scli = 1.0f/(cnti>0.0f?cnti:1.0f);
      for (int ival=0; ival<nval; ++ival)
        vals[icell][ival] *= scli;
      setter.set(cell,vals[icell]);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  private FaultCell _seed; // cell in this skin with highest fl; null, if empty
  private ArrayList<FaultCell> _cellList; // list of cells in this skin
  private FaultCell[][] _cellsAB; // arrays of cells from above to below
  private FaultCell[][] _cellsLR; // arrays of cells from left to right

  private void checkCellArrays() {
    if (_cellsAB!=null)
      checkCellArrays(_cellsAB);
    if (_cellsLR!=null)
      checkCellArrays(_cellsLR);
  }

  private static void checkCellArrays(FaultCell[][] cells) {
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>();
    int ncell = cells.length;
    for (int icell=0; icell<ncell; ++icell) {
      int mcell = cells[icell].length;
      for (int jcell=0; jcell<mcell; ++jcell) {
        FaultCell c = cells[icell][jcell];
        assert cellSet.add(c);
      }
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
