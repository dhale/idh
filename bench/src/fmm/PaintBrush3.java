/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

// for testing
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;

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
     * Vertex array with packed coordinates (x1,x2,x3).
     * The number of vertices equals x.length/3.
     */
    public float[] x;

    /**
     * Optional normal vector array with packed components (u1,u2,u3).
     * If not null, the number of normal vectors equals u.length/3.
     * This number also equals the number of vertices.
     */
    public float[] u;

    /**
     * Triangle array of packed vertex indices (i1,i2,i3). When multiplied by 
     * 3, each index references the first coordinate x1 of a vertex (x1,x2,x3) 
     * stored in the packed vertex array. The number of triangles equals 
     * i.length/3. A vertex may be referenced by more than one triangle.
     */
    public int[] i;
  }

  /**
   * Constructs a paint brush for specified tensors.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param s3 sampling of 3rd dimension.
   * @param pt painting tensors.
   */
  public PaintBrush3(Sampling s1, Sampling s2, Sampling s3, Tensors3 pt) {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _pt = pt;
  }

  /**
   * Sets the location (sample indices) of this brush.
   * @param k1 sample index in 1st dimension.
   * @param k2 sample index in 2nd dimension.
   * @param k3 sample index in 3rd dimension.
   */
  public void setLocation(int k1, int k2, int k3) {
    _k1 = k1;
    _k2 = k2;
    _k3 = k3;
  }

  public Contour getContour() {
    return null;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _s1,_s2,_s3;
  private Tensors3 _pt;
  private int _k1,_k2,_k3;

  private static void trace(String s) {
    System.out.println(s);
  }
}
