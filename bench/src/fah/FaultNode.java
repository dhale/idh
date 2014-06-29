/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fah;

import static fah.FaultUtil.*;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A node on a fault surface, referenced by one or more quads.
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.03
 */
public class FaultNode {

  float fl,fp,ft; // fault likelihood, strike and dip
  float x1,x2,x3; // point on fault
  float u1,u2,u3; // normal vector
  float v1,v2,v3; // strike vector
  float w1,w2,w3; // dip vector
  int na; // number of values accumulated into this node

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
   * @param u1 1st component of node normal vector.
   * @param u2 2nd component of node normal vector.
   * @param u3 3rd component of node normal vector.
   */
  void accumulate(
      float fl,
      float x1, float x2, float x3,
      float u1, float u2, float u3) {

    // If not valid, do nothing.
    if (!isValid())
      return;

    // If normal vectors are inconsistent, this node is not valid.
    if (u1*this.u1+u2*this.u2+u3*this.u3<0.0) {
      this.na = -1;
      this.fl = 0.0f;
    }

    // If node not yet marked invalid, accumulate values.
    // Weight values by non-negative fault likelihood.
    if (na>=0) {
      this.fl += fl;
      this.x1 += fl*x1;
      this.x2 += fl*x2;
      this.x3 += fl*x3;
      this.u1 += fl*u1;
      this.u2 += fl*u2;
      this.u3 += fl*u3;
      this.na += 1;
    }
  }

  /**
   * Completes the computation of fields for this node.
   */
  void complete() {
    assert na>0;
    float scale = 1.0f/fl;
    x1 *= scale; 
    x2 *= scale; 
    x3 *= scale;
    u1 *= scale; 
    u2 *= scale; 
    u3 *= scale;
    float us = 1.0f/sqrt(u1*u1+u2*u2+u3*u3);
    u1 *= us; 
    u2 *= us; 
    u3 *= us;
    float[] u = {u1,u2,u3};
    float[] v = faultStrikeVectorFromNormal(u);
    float[] w = faultDipVectorFromNormal(u);
    v1 = v[0];
    v2 = v[1];
    v3 = v[2];
    w1 = w[0];
    w2 = w[1];
    w3 = w[2];
    fl /= na;
    fp = faultStrikeFromNormal(u);
    ft = faultDipFromNormal(u);
  }

  /**
   * Determines whether or not this node is valid.
   * A node is valid if it has accumulated consistent values.
   * Any inconsistent values, such as normal vectors pointing
   * in opposite directions, will forever invalidate this node.
   */
  boolean isValid() {
    return na>0;
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
