/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import static fah.FaultGeometry.*;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A node on a fault surface, referenced by one or more quads.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.29
 */
public class FaultNode {

  float fl,fp,ft; // fault likelihood, strike and dip
  float x1,x2,x3; // point on fault
  float u1,u2,u3; // dip vector
  float v1,v2,v3; // strike vector
  float w1,w2,w3; // normal vector
  int nac; // number of values accumulated into this node

  /**
   * Accumulates specified values into this node. Nodes are constructed when
   * needed by quads. When first constructed, all node fields are zero. Values
   * needed to set those fields are accumulated as quads are constructed that
   * reference this node. The accumulated values are later used to complete
   * the computation of this node's fields.
   * @param fl fault likelihood.
   * @param x1 1st coordinate of node location.
   * @param x2 2nd coordinate of node location.
   * @param x3 3rd coordinate of node location.
   * @param w1 1st component of node normal vector.
   * @param w2 2nd component of node normal vector.
   * @param w3 3rd component of node normal vector.
   */
  void accumulate(
      float fl,
      float x1, float x2, float x3,
      float w1, float w2, float w3) {
    this.fl += fl;
    this.x1 += x1;
    this.x2 += x2;
    this.x3 += x3;
    this.w1 += w1;
    this.w2 += w2;
    this.w3 += w3;
    this.nac += 1;
  }

  /**
   * Completes the computation of fields for this node.
   * Called only after all field values have been accumulated.
   */
  void complete() {
    assert nac>0;
    float scale = 1.0f/nac;
    fl *= scale;
    x1 *= scale; x2 *= scale; x3 *= scale;
    scale = 1.0f/sqrt(w1*w1+w2*w2+w3*w3);
    w1 *= scale; w2 *= scale; w3 *= scale;
    fp = faultStrikeFromNormalVector(w1,w2,w3);
    ft = faultDipFromNormalVector(w1,w2,w3);
    float[] v = faultStrikeVectorFromStrikeAndDip(fp,ft);
    float[] u = faultDipVectorFromStrikeAndDip(fp,ft);
    v1 = v[0]; v2 = v[1]; v3 = v[2];
    u1 = u[0]; u2 = u[1]; u3 = u[2];
  }

  void blocky() {
    x1 = (int)x1+0.5f;
    x2 = (int)x2+0.5f;
    x3 = (int)x3+0.5f;
  }

  public String toString() {
    return "("+x1+","+x2+","+x3+")";
    //return "("+x1+","+x2+","+x3+"):("+fl+","+fp+","+ft+")";
  }
}
