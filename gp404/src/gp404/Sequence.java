/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gp404;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A sequence is a real-valued sampled function of one variable.
 * A {@link Sequence} combines a {@link edu.mines.jtk.dsp.Sampling} with 
 * a reference to an an array of function values. Because array values 
 * are referenced (not copied), the cost of wrapping any array with a 
 * {@link Sequence} is small.
 * <p> 
 * One consequence of referencing the array of function values is that 
 * changes to elements in such an array are reflected in <em>all</em> 
 * {@link Sequences}s that reference that array. If this behavior is not
 * desired, the copy constructor {@link Sequence#Sequence(Sequence)} 
 * creates a new array copy of function values.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.08.19
 */
public class Sequence {

  /**
   * Constructs a sequence with specified sampling and values zero.
   * @param s the sampling.
   */
  public Sequence(Sampling s) {
    this(s,new float[s.getCount()]);
    _s = s;
    _v = new float[s.getCount()];
  }

  /**
   * Constructs a sequence with specified values and default sampling.
   * The default sampling has number (count) of samples = v.length, 
   * sampling interval (delta) = 1.0 and first sample value (first) = 0.0.
   * @param v array of sequence values; referenced, not copied.
   */
  public Sequence(float[] v) {
    _s = new Sampling(v.length,1.0,0.0);
    _v = v;
  }

  /**
   * Constructs a sequence with specified sampling and values.
   * @param s the sampling.
   * @param v array of sequence values; referenced, not copied.
   *  The array length v.length must equal the number of samples in s.
   */
  public Sequence(Sampling s, float[] v) {
    Check.argument(s.getCount()==v.length,
      "v.length equals the number of samples in s");
    _s = s;
    _v = v;
  }

  /**
   * Constructs a sequence with specified sampling and values zero.
   * @param n the number (count) of samples.
   * @param d the sampling interval (delta).
   * @param f the value of the first sample.
   */
  public Sequence(int n, double d, double f) {
    this(new Sampling(n,d,f));
  }

  /**
   * Constructs a sequence with specified sampling and values.
   * @param n the number (count) of time samples.
   * @param d the sampling interval (delta).
   * @param f the value of the first sample.
   * @param v array of sequence values; referenced, not copied. 
   *  The array length v.length must equal n.
   */
  public Sequence(int n, double d, double f, float[] v) {
    this(new Sampling(n,d,f),v);
  }

  /**
   * Constructs a copy of the specified sampled sequence. This constructor 
   * <em>copies</em> (does not reference) the array of function values from 
   * the specified sequence.
   * @param s the sequence to copy.
   */
  public Sequence(Sequence s) {
    this(s._s,copy(s._v));
  }

  /**
   * Gets the sampling for this sequence.
   * @return the sampling.
   */
  public Sampling getSampling() {
    return _s;
  }

  /**
   * Gets the array of values for this sequence.
   * @return the array of values; by reference, not by copy.
   */
  public float[] getValues() {
    return _v;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  Sampling _s;
  float[] _v;
}
