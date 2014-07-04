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
 * Computes fault skins from images of fault likelihoods, strikes and dips. A
 * fault skin is a linked list of fault cells. Each fault cell is an oriented
 * point located on a ridge in an image of fault likelihood. After a fault
 * cell has been constructed, it can easily be found by simple indexing, as
 * for a 3D image. Each image sample corresponds to either no cell or one
 * cell.
 * <p>
 * A cell has up to four neighbors ("nabors") that lie above, below, left and
 * right of the cell when viewed from above the fault, that is, when looking
 * from the hanging wall toward the footwall. Links to nabors enables cells to
 * form a skin of connected cells, which represents a fault.
 * <p>
 * Links to left and right cell nabors can be used to iterate over all cells
 * along a fault trace, a path of constant depth that is everywhere tangent to
 * fault strike. Likewise, links to cell nabors above and below a cell can be
 * used to iterate up or down a fault. However, this simple up or down
 * iteration need not coincide with a fault curve that is everywhere tangent
 * to fault dip.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.03
 */
public class FaultSkinner {

  /**
   * An oriented point. Fault skins are comprised of cells.
   */
  public static class Cell {

    int i1,i2,i3; // fault cell indices
    float x1,x2,x3; // fault cell coordinates
    float fl,fp,ft; // fault likelihood, strike (phi) and dip (theta)
    float u1,u2,u3; // fault dip vector
    float v1,v2,v3; // fault strike vector
    float w1,w2,w3; // fault normal vector
    Cell ca,cb,cl,cr; // cell nabors above, below, left and right
    Skin skin; // if not null, the skin to which this cell belongs

    Cell(float x1, float x2, float x3, float fl, float fp, float ft) {
      this.i1 = round(x1);
      this.i2 = round(x2);
      this.i3 = round(x3);
      this.x1 = x1; 
      this.x2 = x2; 
      this.x3 = x3;
      this.fl = fl; 
      this.fp = fp; 
      this.ft = ft;
      float[] u = faultDipVectorFromStrikeAndDip(fp,ft);
      float[] v = faultStrikeVectorFromStrikeAndDip(fp,ft);
      float[] w = faultNormalVectorFromStrikeAndDip(fp,ft);
      this.u1 = u[0]; this.u2 = u[1]; this.u3 = u[2];
      this.v1 = v[0]; this.v2 = v[1]; this.v3 = v[2];
      this.w1 = w[0]; this.w2 = w[1]; this.w3 = w[2];
    }

    /*
    // TODO: are these necessary?
    @Override
    public boolean equals(Object object) {
      if (object==this)
        return true;
      if (object!=null && object.getClass()==this.getClass()) {
        Cell that = (Cell)object;
        return this.i1==that.i1 && this.i2==that.i2 && this.i3==that.i3;
      }
      return false;
    }
    @Override
    public int hashCode() {
      return i1^i2^i3;
    }
    */
  }

  /**
   * A linked list of fault cells.
   */
  public static class Skin {
    Cell seed; // the cell in this skin with highest fl; or null, if empty
    ArrayList<Cell> cellList; // list of cells in this skin; or null
    Cell[] cellArray; // array of cells in this skin; or null
    Cell[][] cellsAB; // arrays of arrays of cells ordered from above to below
    Cell[][] cellsLR; // arrays of arrays of cells ordered from left to right

    /**
     * Gets an array of cells in this skin.
     * @return array of cells.
     */
    public Cell[] getCells() {
      if (cellArray==null) {
        cellArray = cellList.toArray(new Cell[0]);
        cellList = null;
      }
      return cellArray;
    }

    /**
     * Adds the specified cell to this skin.
     * @param cell the cell to be added.
     */
    void add(Cell cell) {
      if (this.seed==null)
        this.seed = cell;
      cellList.add(cell);
      cell.skin = this;
    }

    /**
     * Constructs an empty skin.
     */
    Skin() {
      cellList = new ArrayList<Cell>();
    }
  }

  /**
   * Constructs a fault skinner for specified likelihoods and orientations.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   */
  public FaultSkinner(float[][][][] flpt) {
    _fl = flpt[0];
    _fp = flpt[1];
    _ft = flpt[2];
    _n1 = _fl[0][0].length;
    _n2 = _fl[0].length;
    _n3 = _fl.length;
    _fllo = 0.1f;
    _flhi = 0.5f;
    _dflmax = 0.1f;
    _dfpmax = 20.0f;
    _dftmax =  4.0f;
    _dabmax = 0.5f;
    _dcwmin = cos(toRadians(30.0f));
  }

  /**
   * Sets thresholds for fault likelihoods used in skinning. All cells in a
   * skin will have fault likelihoods not less than the lower threshold fllo.
   * At least one cell in a skin will have a fault likelihood not less than
   * the upper threshold flhi.
   * @param fllo the lower threshold.
   * @param flhi the upper threshold.
   */
  public void setLikelihoodThresholds(double fllo, double flhi) {
    _fllo = (float)_fllo;
    _flhi = (float)_flhi;
  }

  /**
   * Sets upper bounds for angle differences used in skinning. Cells can be
   * nabors only if the differences in their strikes and dips are less than
   * these thresholds.
   */
  public void setOrientationThresholds(double pmax, double tmax) {
    _dfpmax = (float)pmax;
    _dftmax = (float)tmax;
  }

  /**
   * Returns array of cells in ridge surfaces of fault likelihood.
   * @return array of cells.
   */
  public Cell[] findCells() {
    int n1 = _n1, n2 = _n2, n3 = _n3;
    float[][][] f = _fl, p = _fp, t = _ft;

    // Smooth fault likelihoods in 2nd and 3rd dimensions. This helps to
    // eliminate spurious ridges, and improves the accuracy of 2nd-order
    // finite-difference approximations (parabolic interpolation) used to
    // locate ridges.
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.applyX0X(f,fs);
    rgf.applyXX0(fs,fs);
    f = fs;

    // Vertical image boundaries are discontinuities that may look like
    // faults. If a fault appears to be near and nearly parallel to image
    // boundaries, then assume it is a boundary artifact and not truly a
    // fault.
    int imax = 5; // max number of samples considered to be near boundary
    float wwmax = 0.75f; // cosine of 30 degrees, squared

    // Loop over all samples. Construct cells for samples nearest to ridges.
    ArrayList<Cell> cells = new ArrayList<Cell>();
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmi = f[i3m][i2 ];
        float[] fim = f[i3 ][i2m];
        float[] fip = f[i3 ][i2p];
        float[] fpi = f[i3p][i2 ];
        float[] fmm = f[i3m][i2m];
        float[] fpp = f[i3p][i2p];
        float[] fmp = f[i3m][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fii = f[i3 ][i2 ];
        float[] pii = p[i3 ][i2 ];
        float[] tii = t[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          float fmii = fmi[i1 ];
          float fimi = fim[i1 ];
          float fipi = fip[i1 ];
          float fpii = fpi[i1 ];
          float fmmi = fmm[i1 ];
          float fppi = fpp[i1 ];
          float fmpi = fmp[i1 ];
          float fpmi = fpm[i1 ];
          float fiii = fii[i1 ];
          float piii = pii[i1 ];
          float tiii = tii[i1 ];

          // Most image samples will not have a fault cell.
          Cell cell = null;

          // If S-N ridge, ...
          if ((fipi<fiii && fimi<fiii) &&
              ((337.5f<=piii || piii<= 22.5f) || 
               (157.5f<=piii && piii<=202.5f))) {
            float f1 = 0.5f*(fipi-fimi); // 1st derivative
            float f2 = fipi-2.0f*fiii+fimi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fl = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fl>=_fllo && (cell==null || fl>cell.fl)) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w1 = w[0], w2 = w[1], w3 = w[2];
              if (imax<=i2 && i2<_n2-imax || w2*w2<=wwmax) {
                float x2 = i2+dr;
                cell = new Cell(i1,x2,i3,fl,piii,tiii);
              }
            }
          }

          // If SW-NE ridge, ...
          if ((fmpi<fiii && fpmi<fiii) &&
              (( 22.5f<=piii && piii<= 67.5f) || 
               (202.5f<=piii && piii<=247.5f))) {
            float f1 = 0.5f*(fmpi-fpmi); // 1st derivative
            float f2 = fmpi-2.0f*fiii+fpmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fl = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fl>=_fllo && (cell==null || fl>cell.fl)) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w1 = w[0], w2 = w[1], w3 = w[2];
              if ((imax<=i2 && i2<_n2-imax || w2*w2<=wwmax) &&
                  (imax<=i3 && i3<_n3-imax || w3*w3<=wwmax)) {
                float x2 = i2+dr;
                float x3 = i3-dr;
                cell = new Cell(i1,x2,x3,fl,piii,tiii);
              }
            }
          }

          // If W-E ridge, ...
          if ((fpii<fiii && fmii<fiii) &&
              ((67.5f<=piii && piii<=112.5f) ||
               (247.5f<=piii && piii<=292.5f))) {
            float f1 = 0.5f*(fpii-fmii); // 1st derivative
            float f2 = fmii-2.0f*fiii+fpii; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fl = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fl>=_fllo && (cell==null || fl>cell.fl)) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w1 = w[0], w2 = w[1], w3 = w[2];
              if (imax<=i3 && i3<_n3-imax || w3*w3<=wwmax) {
                float x3 = i3+dr;
                cell = new Cell(i1,i2,x3,fl,piii,tiii);
              }
            }
          }

          // If NW-SE ridge, ...
          if ((fppi<fiii && fmmi<fiii) &&
              ((112.5f<=piii && piii<=157.5f) || 
               (292.5f<=piii && piii<=337.5f))) {
            float f1 = 0.5f*(fppi-fmmi); // 1st derivative
            float f2 = fppi-2.0f*fiii+fmmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fl = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fl>=_fllo && (cell==null || fl>cell.fl)) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w1 = w[0], w2 = w[1], w3 = w[2];
              if ((imax<=i2 && i2<_n2-imax || w2*w2<=wwmax) &&
                  (imax<=i3 && i3<_n3-imax || w3*w3<=wwmax)) {
                float x2 = i2+dr;
                float x3 = i3+dr;
                cell = new Cell(i1,x2,x3,fl,piii,tiii);
              }
            }
          }

          // If we constructed a cell, add it to the list.
          if (cell!=null)
            cells.add(cell);
        }
      }
    }
    return cells.toArray(new Cell[0]);
  }

  /**
   * Returns an array of skins comprised of specified cells. Skins will
   * include only those cells with fault likelihoods not less than the higher
   * skinning threshold.
   * @param cells array of cells; will be sorted by fault likelihood.
   * @return array of skins.
   */
  public Skin[] findSkins(Cell[] cells) {
    return skin(cells);
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of quad coordinates, normals and colors.
   * Each quad represents one cell.
   */
  public static float[][] getXyzUvwRgb(Cell[] cells) {
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList fcl = new FloatList();
    float size = 0.3f;
    float[] qa = {0.0f,-size,-size};
    float[] qb = {0.0f, size,-size};
    float[] qc = {0.0f, size, size};
    float[] qd = {0.0f,-size, size};
    for (Cell cell:cells) {
      float x1 = cell.x1;
      float x2 = cell.x2;
      float x3 = cell.x3;
      float w1 = cell.w1;
      float w2 = cell.w2;
      float w3 = cell.w3;
      float fl = cell.fl;
      float fp = toRadians(cell.fp);
      float ft = toRadians(cell.ft);
      float cp = cos(fp);
      float sp = sin(fp);
      float ct = cos(ft);
      float st = sin(ft);
      float[] ra = rotatePoint(cp,sp,ct,st,qa);
      float[] rb = rotatePoint(cp,sp,ct,st,qb);
      float[] rc = rotatePoint(cp,sp,ct,st,qc);
      float[] rd = rotatePoint(cp,sp,ct,st,qd);
      float a1 = x1+ra[0], a2 = x2+ra[1], a3 = x3+ra[2];
      float b1 = x1+rb[0], b2 = x2+rb[1], b3 = x3+rb[2];
      float c1 = x1+rc[0], c2 = x2+rc[1], c3 = x3+rc[2];
      float d1 = x1+rd[0], d2 = x2+rd[1], d3 = x3+rd[2];
      xyz.add(a3); xyz.add(a2); xyz.add(a1); fcl.add(fl);
      xyz.add(b3); xyz.add(b2); xyz.add(b1); fcl.add(fl);
      xyz.add(c3); xyz.add(c2); xyz.add(c1); fcl.add(fl);
      xyz.add(d3); xyz.add(d2); xyz.add(d1); fcl.add(fl);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
    }
    float[] fc = fcl.trim();
    float fcmin = 0.0f;
    float fcmax = 1.0f;
    ColorMap cmap = new ColorMap(fcmin,fcmax,ColorMap.JET);
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Used to quickly search for potential cell nabors.
  private static class CellArray {
    CellArray(int n1, int n2, int n3, Cell[] cells) {
      _n1 = n1; _n2 = n2; _n3 = n3;
      _cells = new Cell[n3][n2][n1];
      if (cells!=null) {
        for (Cell cell:cells) {
          int i1 = cell.i1;
          int i2 = cell.i2;
          int i3 = cell.i3;
          _cells[i3][i2][i1] = cell;
        }
      }
    }
    Cell get(int i1, int i2, int i3) {
      if (0<=i1 && i1<_n1 && 0<=i2 && i2<_n2 && 0<=i3 && i3<_n3) {
        return _cells[i3][i2][i1];
      } else {
        return null;
      }
    }
    void set(int i1, int i2, int i3, Cell cell) {
      _cells[i3][i2][i1] = cell;
    }
    private int _n1,_n2,_n3;
    private Cell[][][] _cells;
  }

  private int _n1,_n2,_n3; // dimensions of fault images fl, fp and ft
  private float[][][] _fl,_fp,_ft; // fault likelihoods, strikes and dips
  private float _fllo; // lower threshold on fault likelihoods
  private float _flhi; // higher threshold on fault likelihoods
  private float _dflmax; // max difference between likelihoods of nabors
  private float _dfpmax; // max difference between strikes of nabors
  private float _dftmax; // max difference between dips of nabors
  private float _dabmax; // max distance to planes of nabors
  private float _dcwmin; // min dot product of normals of nabors
  private static Comparator<Cell> _cellComparator = new Comparator<Cell>() {
    public int compare(Cell c1, Cell c2) {
      if (c1.fl<c2.fl)
        return -1;
      else if (c1.fl>c2.fl)
        return 1;
      else
        return 0;
    }
  };

  private Skin[] skin(Cell[] cells) {

    // Array of cells used to quickly find cell nabors.
    CellArray cellArray = new CellArray(_n1,_n2,_n3,cells);

    // Empty list of skins.
    ArrayList<Skin> skinList = new ArrayList<Skin>();

    // Sort array of cells by increasing fault likelihoods.
    Arrays.sort(cells,_cellComparator);

    // While cells with high fault likelihood remain, ...
    int kseed = cells.length-1;
    while (kseed>=0 && cells[kseed].fl>_flhi) {

      // Look for a skinless cell with high fault likelihood.
      while (kseed>=0 && cells[kseed].skin!=null)
        --kseed;

      // If we found a cell with which to construct a new skin, ...
      if (kseed>=0 && cells[kseed].fl>_flhi) {
        Cell seed = cells[kseed];

        // Make a new empty skin and add it to the skin list.
        Skin skin = new Skin();
        skinList.add(skin);

        // Make a new sorted set of cells used to grow the skin from the seed.
        TreeSet<Cell> growSet = new TreeSet<Cell>(_cellComparator);
        growSet.add(seed);

        // While the grow set is not empty, ...
        while (!growSet.isEmpty()) {

          // Get the cell with highest fault likelihood, while removing it
          // from the grow set, and add it to the skin.
          Cell cell = growSet.pollLast();
          skin.add(cell);

          // Link with all good nabors, and add them to the grow set.
          Cell ca,cb,cl,cr;
          ca = findNaborAbove(cellArray,cell);
          cb = findNaborBelow(cellArray,ca);
          if (ca!=null && ca.skin==null && cb==cell) {
            linkAboveBelow(ca,cb);
            growSet.add(ca);
          }
          cb = findNaborBelow(cellArray,cell);
          ca = findNaborAbove(cellArray,cb);
          if (cb!=null && cb.skin==null && ca==cell) {
            linkAboveBelow(ca,cb);
            growSet.add(cb);
          }
          cl = findNaborLeft(cellArray,cell);
          cr = findNaborRight(cellArray,cl);
          if (cl!=null && cl.skin==null && cr==cell) {
            linkLeftRight(cl,cr);
            growSet.add(cl);
          }
          cr = findNaborRight(cellArray,cell);
          cl = findNaborLeft(cellArray,cr);
          if (cr!=null && cr.skin==null && cl==cell) {
            linkLeftRight(cl,cr);
            growSet.add(cr);
          }
        }
      }
    }
    return skinList.toArray(new Skin[0]);
  }

  // Returns true if the specified cells are nabors. This method assumes that
  // all links are mutual. For example, if c1 is the nabor above c2, then c2
  // must be the nabor below c1.
  private static boolean areNabors(Cell c1, Cell c2) {
    return c1.ca==c2 || c1.cb==c2 || c1.cl==c2 || c1.cr==c2;
  }

  // Methods to link mutually best nabors.
  private void linkLeftRight(Cell cl, Cell cr) {
    cr.cl = cl;
    cl.cr = cr;
  }
  private void linkAboveBelow(Cell ca, Cell cb) {
    cb.ca = ca;
    ca.cb = cb;
  }

  // Methods to return good nabors of a specified cell. Return null if no
  // nabor is good enough, based on various thresholds.
  private Cell findNaborAbove(CellArray cells, Cell cell) {
    if (cell==null) return null;
    if (cell.ca!=null) return cell.ca;
    return findNaborAboveBelow(true,cells,cell);
  }
  private Cell findNaborBelow(CellArray cells, Cell cell) {
    if (cell==null) return null;
    if (cell.cb!=null) return cell.cb;
    return findNaborAboveBelow(false,cells,cell);
  }
  private Cell findNaborLeft(CellArray cells, Cell cell) {
    if (cell==null) return null;
    if (cell.cl!=null) return cell.cl;
    return findNaborLeftRight(true,cells,cell);
  }
  private Cell findNaborRight(CellArray cells, Cell cell) {
    if (cell==null) return null;
    if (cell.cr!=null) return cell.cr;
    return findNaborLeftRight(false,cells,cell);
  }
  private Cell findNaborAboveBelow(boolean above, CellArray cells, Cell cell) {
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
      k1 = -1;
      u1 = -u1;
      u2 = -u2;
      u3 = -u3;
    }
    Cell cmax = null;
    float dmax = 0.0f;
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        Cell c = cells.get(i1+k1,i2+k2,i3+k3);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float d3 = c.x3-x3;
          float d = (d1*u1+d2*u2+d3*u3)/sqrt(d1*d1+d2*d2+d3*d3);
          if (d>dmax) {
            cmax = c;
            dmax = d;
          }
        }
      }
    }
    return canBeNabors(cell,cmax) ? cmax : null;
  }
  private Cell findNaborLeftRight(boolean left, CellArray cells, Cell cell) {
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
    Cell cmax = null;
    float dmax = 0.0f;
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        if (k2==0 && k3==0)
          continue;
        Cell c = cells.get(i1,i2+k2,i3+k3);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float d3 = c.x3-x3;
          float d = (d1*v1+d2*v2+d3*v3)/sqrt(d1*d1+d2*d2+d3*d3);
          if (d>dmax) {
            cmax = c;
            dmax = d;
          }
        }
      }
    }
    return canBeNabors(cell,cmax) ? cmax : null;
  }

  // Returns true if two specified cells can be nabors. The two cells are
  // assumed to be within one sample of each other. This method uses other
  // attributes of the cells to determine whether or not they can be nabors.
  private boolean canBeNabors(Cell ca, Cell cb) {
    if (ca==null || cb==null)
      return false;
    if (minFl(ca,cb)<_fllo)
      return false;
    if (absDeltaFl(ca,cb)>_dflmax)
      return false;
    if (absDeltaFp(ca,cb)>_dfpmax)
      return false;
    if (absDeltaFt(ca,cb)>_dftmax)
      return false;
    if (dotNormalVectors(ca,cb)<_dcwmin)
      return false;
    if (maxDistanceToPlane(ca,cb)>_dabmax)
      return false;
    return true;
  }
  private static float minFl(Cell ca, Cell cb) {
    return min(ca.fl,cb.fl);
  }
  private static float absDeltaFl(Cell ca, Cell cb) {
    return abs(ca.fl-cb.fl);
  }
  private static float absDeltaFp(Cell ca, Cell cb) {
    float del = ca.fp-cb.fp;
    return min(abs(del),abs(del+360.0f),abs(del-360.0f));
  }
  private static float absDeltaFt(Cell ca, Cell cb) {
    return abs(ca.ft-cb.ft);
  }
  private static float dotNormalVectors(Cell ca, Cell cb) {
    float aw1 = ca.w1, aw2 = ca.w2, aw3 = ca.w3;
    float bw1 = cb.w1, bw2 = cb.w2, bw3 = cb.w3;
    return aw1*bw1+aw2*bw2+aw3*bw3;
  }
  private static float maxDistanceToPlane(Cell ca, Cell cb) {
    float aw1 = ca.w1, aw2 = ca.w2, aw3 = ca.w3;
    float ax1 = ca.x1, ax2 = ca.x2, ax3 = ca.x3;
    float bw1 = cb.w1, bw2 = cb.w2, bw3 = cb.w3;
    float bx1 = cb.x1, bx2 = cb.x2, bx3 = cb.x3;
    float dx1 = ax1-bx1;
    float dx2 = ax2-bx2;
    float dx3 = ax3-bx3;
    float dab = aw1*dx1+aw2*dx2+aw3*dx3;
    float dba = bw1*dx1+bw2*dx2+bw3*dx3;
    return max(dab,dba);
  }

  // Rotates a specified point by strike (phi) and dip (theta) angles,
  // given specified cosines (cp and ct) and sines (sp and st) of those 
  // angles. The order of transformation is
  // (1) rotate around axis x3 by dip angle
  // (2) rotate around axis x1 by strike angle
  // Returns the coordinates of the rotated point.
  private static float[] rotatePoint(
      float cp, float sp, float ct, float st, float[] x) {
    float x1 = x[0], x2 = x[1], x3 = x[2];
    float y1 =     ct*x1+   st*x2;
    float y2 = -cp*st*x1+cp*ct*x2+sp*x3;
    float y3 =  sp*st*x1-sp*ct*x2+cp*x3;
    return new float[]{y1,y2,y3};
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static class FloatList {
    public int n;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
}
