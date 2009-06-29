/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.Tensors3;
import static edu.mines.jtk.util.MathPlus.max;
import static edu.mines.jtk.util.MathPlus.min;

/**
 * A paint brush for 3D images.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.08.21
 */
public class PaintBrush3 {

  /**
   * A contour surface represented by lists of vertices and triangles.
   */
  public static class Contour {

    /**
     * Vertex array with packed coordinates (x3,x2,x1). (Note x3 is first.)
     * The number of vertices equals x.length/3.
     */
    public float[] x;

    /**
     * Optional normal vector array with packed components (u3,u2,u1).
     * If not null, the number of normal vectors equals u.length/3.
     * This number also equals the number of vertices.
     */
    public float[] u;

    /**
     * Triangle array of packed vertex indices (i1,i2,i3). When multiplied by 
     * 3, each index references the coordinate x3 of a vertex (x3,x2,x1) 
     * stored in the packed vertex array. The number of triangles equals 
     * i.length/3. A vertex may be referenced by more than one triangle.
     */
    public int[] i;
  }

  /**
   * Constructs a paint brush for specified tensors.
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param n3 number of samples in 3rd dimension.
   * @param pt painting tensors.
   */
  public PaintBrush3(int n1, int n2, int n3, Tensors3 pt) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _pt = pt;
    setSize(10);
  }

  /**
   * Gets the location (k1,k2,k3) of this brush.
   * @return array {k1,k2,k3} with location.
   */
  public int[] getLocation() {
    int[] k = new int[3];
    k[0] = _k1;
    k[1] = _k2;
    k[2] = _k3;
    return k;
  }

  /**
   * Sets the location (k1,k2,k3) of this brush.
   * The default location is (0,0,0).
   * @param k1 sample index in 1st dimension.
   * @param k2 sample index in 2nd dimension.
   * @param k3 sample index in 3rd dimension.
   */
  public void setLocation(int k1, int k2, int k3) {
    if (_k1!=k1 || _k2!=k2 || _k3!=k3) {
      _k1 = k1;
      _k2 = k2;
      _k3 = k3;
      _dirty = true;
    }
  }

  /**
   * Gets the size of this paint brush.
   * @return the size.
   */
  public int getSize() {
    return _size;
  }

  /**
   * Sets the size of this paint brush. For an identity painting tensor
   * (for which time equals distance), the brush size equals its radius.
   * A brush with size zero covers only the one sample where it is located.
   * The default size is ten samples.
   * @param size the size.
   */
  public void setSize(int size) {
    _size = max(0,size);

    // Maximum time.
    _tmax = (float)(1+size);

    // Half of number of samples in array of brush times; modify as necessary.
    int nh = 1;
    while (nh<=size)
      nh += nh;
    
    // If number of samples in array of brush times has increased,
    // construct new brush tensors, time solver and marching cubes.
    if (_nh!=nh) {
      int nb = 1+2*nh;
      double db = 1.0;
      double fb = -nh;
      Sampling sb = new Sampling(nb,db,fb);
      _bt = new BrushTensors3();
      _ts = new TimeSolver3(nb,nb,nb,_bt);
      _mc = new MarchingCubes(sb,sb,sb,_ts.getTimes());
      _mc.setSwap13(true);
      _nb = nb;
      _nh = nh;
      trace("setSize: nh="+nh);
    }

    // Set maximum time for time solver, and note solution not valid.
    _ts.setMaxTime(2.0f*_tmax);
    _dirty = true;
  }

  public boolean[][][] getMask() {
    float[][][] ts = _ts.getTimes();
    int nm = 1+2*_size;
    boolean[][][] mask = new boolean[nm][nm][nm];
    int nmask = 0;
    for (int i3=0,j3=_nh-_size; i3<nm; ++i3,++j3) {
      for (int i2=0,j2=_nh-_size; i2<nm; ++i2,++j2) {
        for (int i1=0,j1=_nh-_size; i1<nm; ++i1,++j1) {
          float tj = ts[j3][j2][j1];
          if (tj<=_tmax) {
            mask[i3][i2][i1] = true;
            ++nmask;
          }
        }
      }
    }
    return mask;
  }

  /** 
   * Gets the contour for this paint brush.
   * @return the contour.
   */
  public Contour getContour() {
    if (_dirty) {
      _ts.reset();
      _ts.zeroAt(_nh,_nh,_nh);
      Sampling s1 = new Sampling(_nb,1.0,_k1-_nh);
      Sampling s2 = new Sampling(_nb,1.0,_k2-_nh);
      Sampling s3 = new Sampling(_nb,1.0,_k3-_nh);
      _mc.setSampling(s1,s2,s3);
      MarchingCubes.Contour contour = _mc.getContour((float)_size);
      _contour = new Contour();
      _contour.x = contour.x;
      _contour.u = contour.u;
      _contour.i = contour.i;
      _dirty = false;
    }
    return _contour;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Brush tensors are a subset of the painting tensors.
  private class BrushTensors3 implements Tensors3 {
    public void getTensor(int i1, int i2, int i3, float[] a) {
      i1 = max(0,min(_n1-1,i1+_k1-_nh));
      i2 = max(0,min(_n2-1,i2+_k2-_nh));
      i3 = max(0,min(_n3-1,i3+_k3-_nh));
      _pt.getTensor(i1,i2,i3,a);
    }
  }

  private int _n1,_n2,_n3;
  private Tensors3 _pt;
  private int _k1,_k2,_k3;
  private int _size,_nb,_nh;
  private float _tmax;
  private MarchingCubes _mc;
  private BrushTensors3 _bt;
  private TimeSolver3 _ts;
  private boolean _dirty;
  private Contour _contour;

  private static void trace(String s) {
    System.out.println(s);
  }
}
