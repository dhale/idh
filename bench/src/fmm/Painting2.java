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
 * A 2D array of painted values, where most values are painted automatically.
 * Except for a relatively small number of fixed samples painted explicitly, 
 * most samples are painted by extrapolation and interpolation guided by
 * structure tensors. Intuitively, paint flows slowly through locations with
 * relatively high structure. Structure tensors may also cause paint to
 * flow anisotropically, in directions corresponding to relatively low
 * structure.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.07.28
 */
public class Painting2 {

  /**
   * Constructs a painting for the specified structure tensor field.
   * All painted values are initially clear (null).
   * @param n1 number of samples in 1st dimension.
   * @param n2 number of samples in 2nd dimension.
   * @param nv number of values painted for each sample
   * @param pt painting tensors.
   */
  public Painting2(int n1, int n2, int nv, Tensors2 pt) {
    _n1 = n1;
    _n2 = n2;
    _nv = nv;
    _dv = new float[nv];
    _k1 = new int[n2][n1];
    _k2 = new int[n2][n1];
    _tk = new float[n2][n1];
    _vk = new float[n2][n1][];
    _type = new byte[n2][n1];
    _heap = new TimeHeap2(TimeHeap2.Type.MAX,n1,n2);
    _tsol = new TimeSolver2(_tk,pt);
    _tsol.setConcurrency(TimeSolver2.Concurrency.PARALLEL);
    clearAll();
  }

  /**
   * Sets the tensors used in this painting.
   * @param pt painting tensors.
   */
  public void setTensors(Tensors2 pt) {
    _tsol.setTensors(pt);
  }

  /**
   * Sets the default value for value index zero.
   * Default values are used for all clear (not painted) samples.
   * @param value the default value.
   */
  public void setDefaultValue(float value) {
    _dv[0] = value;
  }

  /**
   * Sets all default values.
   * Default values are used for all clear (not painted) samples.
   * @param value array[nv] of default values.
   */
  public void setDefaultValues(float[] values) {
    for (int kv=0; kv<_nv; ++kv)
      _dv[kv] = values[kv];
  }

  /**
   * Clears all painted values, including all fixed values.
   */
  public void clearAll() {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        clear(i1,i2);
      }
    }
  }

  /**
   * Clears painted values that are not fixed.
   */
  public void clearNotFixed() {
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_type[i2][i1]!=FIXED) {
          clear(i1,i2);
        }
      }
    }
  }

  /**
   * Erases values for any fixed sample with specified indices.
   * If the specified sample is not fixed, this method does nothing.
   * @param i1 index in 1st dimension of sample to erase.
   * @param i2 index in 2nd dimension of sample to erase.
   */
  public void eraseFixedAt(int i1, int i2) {
    if (_type[i2][i1]==FIXED)
      clear(i1,i2);
  }
  private void clear(int i1, int i2) {
    _type[i2][i1] = CLEAR;
    _k1[i2][i1] = -1;
    _k2[i2][i1] = -1;
    _tk[i2][i1] = INFINITY;
    _vk[i2][i1] = null;
  }

  /**
   * Paints the specified sample with one specified value at index zero.
   * Paints default values for indices greater than zero.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param v value at index zero for the painted sample.
   */
  public void paintAt(int i1, int i2, float v) {
    _type[i2][i1] = FIXED;
    _k1[i2][i1] = i1;
    _k2[i2][i1] = i2;
    _tk[i2][i1] = INFINITY;
    _vk[i2][i1] = new float[_nv];
    _vk[i2][i1][0] = v;
    for (int iv=1; iv<_nv; ++iv)
      _vk[i2][i1][iv] = _dv[iv];
  }

  /**
   * Paints the specified sample with specified values.
   * After painting, the specified sample is fixed; its values will
   * not change in any subsequent extrapolation or interpolation.
   * @param i1 index in 1st dimension of sample to paint.
   * @param i2 index in 2nd dimension of sample to paint.
   * @param v array of values for painted sample; by copy, not by reference.
   */
  public void paintAt(int i1, int i2, float[] v) {
    _type[i2][i1] = FIXED;
    _k1[i2][i1] = i1;
    _k2[i2][i1] = i2;
    _tk[i2][i1] = INFINITY;
    _vk[i2][i1] = Array.copy(v);
  }

  /**
   * Extrapolates values from all fixed (explicitly painted) samples.
   * After extrapolation, all samples are either fixed or extrapolated.
   */
  public void extrapolate() {

    // Extrapolators listens for times decreased by time solver.
    Extrapolator ex = new Extrapolator();
    _tsol.addListener(ex);

    // Clear all samples that are not fixed, and insert all fixed 
    // samples into the max-heap with huge (invalid) times that can
    // only get smaller. After the max-heap is built, any one of the
    // fixed samples could be at the top of the heap.
    _heap.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        _tk[i2][i1] = INFINITY;
        if (_type[i2][i1]==FIXED) {
          _heap.insert(i1,i2,INFINITY);
        } else {
          _type[i2][i1] = CLEAR;
        }
      }
    }

    // Extrapolate from all fixed samples, one at a time, in order of
    // decreasing time. The fixed sample with the largest time is at the
    // top of the max-heap. As we extrapolate from this fixed sample, times 
    // for other fixed samples that remain in the heap may be reduced, so 
    // that their order may change. We choose a decreasing order to reduce 
    // the number of times that must be reduced during extrapolation 
    // from each fixed sample.
    while (!_heap.isEmpty()) {

      // Remove from the max-heap the fixed sample with largest time.
      TimeHeap2.Entry ef = _heap.remove();
      int k1 = ef.i1;
      int k2 = ef.i2;

      // The values to be extrapolated.
      float[] vk = _vk[k2][k1];

      // Prepare extrapolator to listen to time solver.
      ex.set(k1,k2,vk);

      // Zero the time at the fixed sample and compute times to neighbors.
      _tsol.zeroAt(k1,k2);
    }

    _tsol.removeListener(ex);
  }

  /**
   * Interpolates values from all fixed samples, using extrapolated samples.
   * After interpolation, all samples are either fixed or interpolated.
   */
  public void interpolate() {
  /*

    // Interpolation occurs in two stages. Both stages compute times by 
    // fast marching away from the sample to be interpolated. In stage 1, 
    // times and values for samples reached while marching are modified,
    // so that each sample interpolated in this first stage will affect
    // samples interpolated later. In stage 2, times (but not values) 
    // are again modified during marching, but are restored after marching. 
    // Therefore, values interpolated in stage 2 do not affect other 
    // interpolated values. 
    boolean stage1 = true;
    boolean stage2 = false;

    // In stage 2, times that must be restored are saved in a list while 
    // marching, and values interpolated are saved in a separate array so
    // that they will not affect other interpolated values. At the end of
    // stage 2, we will merge the two arrays of values. In stage 1, both
    // the time list and values array are null and not used.
    TimeList tl = null;
    float[][][] va = null;

    // Minimum number of samples to interpolate in stage 2, a fraction 
    // of the the total number of samples. This is an important parameter. 
    // Higher fractions close to one yield smoother interpolations, but
    // can be much more costly than lower fractions.
    int nstage2 = (int)(0.05*_n1*_n2);

    // Insert all extrapolated samples into the max-heap with their
    // current times. After the max-heap is built, the extrapolated
    // sample with largest time is at the top of the heap.
    _heap.clear();
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        if (_type[i2][i1]==EXTRA) {
          _heap.insert(i1,i2,_tk[i2][i1]);
        }
      }
    }

    // Interpolate all extrapolated (not-fixed) samples, one at a time, 
    // in order of decreasing time. The extrapolated sample with the 
    // largest time is at the top of the max-heap. In stage 1, as we 
    // interpolate this sample, times for other extrapolated samples 
    // that remain in the heap may be reduced, so that their order may 
    // change. We choose a decreasing order in stage 1 to reduce the 
    // number of extrapolated samples that are reached while marching.
    // In stage 2, the order will not matter.
    while (!_heap.isEmpty()) {

      // Switch to stage 2 when number of samples left is small enough. 
      if (stage1 && _heap.size()<nstage2) {
        stage1 = false;
        stage2 = true;
        tl = new TimeList();
        va = new float[_n2][_n1][];
      }

      // Remove from the max-heap the extrapolated sample with largest time.
      // This is the sample to be interpolated; the "interpolated sample".
      TimeHeap2.Entry te = _heap.remove();
      int k1 = te.i1;
      int k2 = te.i2;
      float tk = te.t;

      // The values to be interpolated. In stage 1, this array will be
      // assigned to all extrapolated samples during the march away from 
      // the interpolated sample. This assignment is one reason that an 
      // array of values is so useful, for we will not actually know the 
      // interpolated values until the march is complete. Later, when we 
      // have completed the computation of the interpolated values, those 
      // values will already be referenced by all extrapolated samples 
      // nearest to the interpolated sample.
      // In stage 2, interpolated values are not extrapolated, and must not
      // affect other interpolated values, so we store them in a separate
      // array of values to be merged later.
      float[] vk = Array.copy(_vk[k2][k1]);
      if (stage1) {
        _vk[k2][k1] = vk;
      } else {
        va[k2][k1] = vk;
      }

      // In stage 2, save the time for the interpolated sample.
      if (stage2) {
        tl.clear();
        tl.append(k1,k2,tk);
      }

      // Count of values accumulated for the interpolated sample.
      int nk = 1;

      // Mark all samples as far, set the type of the interpolated sample,
      // mark the interpolated sample as known with time zero, and update 
      // its neighbors.
      clearMarks();
      _hmin.clear();
      _type[k2][k1] = INTER;
      _mark[k2][k1] = _known;
      _k1[k2][k1] = k1;
      _k2[k2][k1] = k2;
      _tk[k2][k1] = 0.0f;
      updateNabors(k1,k2,tl);

      // March away from the interpolated sample to all extrapolated
      // samples that are nearer to the interpolated sample than to any 
      // fixed samples or samples previously interpolated in stage 1.
      // While marching, accumulate values needed for interpolation.
      while (!_hmin.isEmpty()) {

        // Get the extrapolated sample with minimum time.
        TimeHeap2.Entry e = _hmin.remove();
        int i1 = e.i1;
        int i2 = e.i2;
        float ti = e.t;

        // Accumulate values for the extrapolated sample.
        float[] vki = _vk[i2][i1];
        for (int iv=0; iv<_nv; ++iv)
          vk[iv] += vki[iv];
        ++nk;

        // Mark the extrapolated sample known. In stage 1, reduce it's
        // time in the max-heap. Also, in stage 1, it's values will be 
        // those of the interpolated sample. Continue marching by updating 
        // the nabor samples.
        _mark[i2][i1] = _known;
        if (stage1) {
          _heap.reduce(i1,i2,ti);
          _k1[i2][i1] = k1;
          _k2[i2][i1] = k2;
          _vk[i2][i1] = vk;
        }
        updateNabors(i1,i2,tl);
      }

      // The march is complete, and nearby values have been accumulated.
      // Now simply divide by the number of values accumulated.
      float vs = 1.0f/(float)nk;
      for (int iv=0; iv<_nv; ++iv)
        vk[iv] *= vs;

      // In stage 2, restore any times saved during marching.
      if (stage2) {
        int nl = tl.n;
        int[] k1l = tl.k1List;
        int[] k2l = tl.k2List;
        float[] tkl = tl.tkList;
        for (int il=0; il<nl; ++il)
          _tk[k2l[il]][k1l[il]] = tkl[il];
      }
    }

    // Finally, merge any interpolated values that were saved in stage 2.
    if (va!=null) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          float[] vi = va[i2][i1];
          if (vi!=null)
            _vk[i2][i1] = vi;
        }
      }
    }
  */
  }

  /**
   * Gets the array of times used in this painting.
   * These times are modified by extrapolation and interpolation.
   * @return array of times; by reference, not by copy.
   */
  public float[][] getTimes() {
    return _tk;
  }

  /**
   * Gets the painted value with index zero for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   */
  public float getValue(int i1, int i2) {
    return getValue(i1,i2,0);
  }

  /**
   * Gets the painted value for the specified sample.
   * Returns a default value if the specified sample is not painted.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param iv index of value to get.
   */
  public float getValue(int i1, int i2, int iv) {
    float[] vk = _vk[i2][i1];
    return (vk!=null)?vk[iv]:_dv[iv];
  }

  /**
   * Gets a copy of the painted values with index zero.
   * Returns default values for any samples not painted.
   * @return array of values; by copy, not by reference.
   */
  public float[][] getValues() {
    return getValues(0);
  }

  /**
   * Gets a copy of the painted values with specified index.
   * Returns default values for any samples not painted.
   * @param iv index of values to get.
   * @return array of values; by copy, not by reference.
   */
  public float[][] getValues(int iv) {
    float[][] v = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        v[i2][i1] = getValue(i1,i2,iv);
      }
    }
    return v;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // The value for times not yet computed. Also the value returned by
  // methods that compute times from nabor times when a valid time
  // cannot be computed. We use the maximum possible float so that
  // it will be larger than any valid times we compute.
  private static final float INFINITY = Float.MAX_VALUE;

  // Type of paint.
  private static final byte CLEAR = 0; // values null (not painted)
  private static final byte FIXED = 1; // values painted explicitly
  private static final byte EXTRA = 2; // values painted by extrapolation
  private static final byte INTER = 3; // values painted by interpolation

  private int _n1,_n2; // painting dimensions
  private int _nv; // number of values associated with each sample
  private float[] _dv; // default values for clear samples
  private float[][][] _vk; // painted values; null for clear samples
  private float[][] _tk; // time to nearest fixed or interpolated sample
  private int[][] _k1,_k2; // indices of nearest fixed or interp sample
  private byte[][] _type; // clear, fixed, extra, or inter
  private TimeHeap2 _heap; // max heap of times
  private TimeSolver2 _tsol; // computes times

  /**
   * Listens for times decreased by the time solver during extrapolation.
   */
  private class Extrapolator implements TimeSolver2.Listener {
    private int k1,k2;
    private float[] vk;
    void set(int k1, int k2, float[] vk) {
      this.k1 = k1;
      this.k2 = k2;
      this.vk = vk;
    }
    public void timeDecreased(int i1, int i2, float t) {
      if (_type[i2][i1]==FIXED) {
        if (i1!=k1 || i2!=k2)
          _heap.reduce(i1,i2,t);
      } else {
        _type[i2][i1] = EXTRA;
        _k1[i2][i1] = k1;
        _k2[i2][i1] = k2;
        _vk[i2][i1] = vk;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static void sleep(int ms) {
    try {
      Thread.currentThread().sleep(1000);
    } catch (InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static void plotImageTensors(float[][] x, EigenTensors2 et) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    //sp.setSize(650,600);
    //sp.setSize(950,900);
    sp.setSize(1130,820);
    PixelsView pv = sp.addPixels(x);
    //pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    int n1 = x[0].length;
    int n2 = x.length;
    float[][][] x12 = getTensorEllipses(n1,n2,10,et);
    float[][] x1 = x12[0];
    float[][] x2 = x12[1];
    PointsView ev = new PointsView(x1,x2);
    ev.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    ev.setLineColor(Color.YELLOW);
    sp.getPlotPanel().getTile(0,0).addTiledView(ev);
  }

  private static class Plot {
    Plot(float[][] f) {
      this(f,null);
    }
    Plot(float[][] f, IndexColorModel icm) {
      _sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
      //_sp.setSize(650,600);
      _sp.setSize(950,900);
      //_sp.setSize(1130,820);
      _pv = _sp.addPixels(f);
      if (icm==null) icm = ColorMap.JET;
      _pv.setColorModel(icm);
      //_pv.setInterpolation(PixelsView.Interpolation.NEAREST);
      _pv.setInterpolation(PixelsView.Interpolation.LINEAR);
    }
    void set(final float[][] f) {
      sleep(1000);
      SwingUtilities.invokeLater(new Runnable() {
        public void run() {
          _pv.set(f);
        }
      });
    }
    final private SimplePlot _sp;
    final private PixelsView _pv;
  }

  private static void plot(float[][] f) {
    plot(f,null);
  }

  private static void plot(float[][] f, IndexColorModel icm) {
    new Plot(f,icm);
  }

  private static float[][] readImage(int n1, int n2, String fileName) {
    try {
      //java.nio.ByteOrder bo = java.nio.ByteOrder.LITTLE_ENDIAN;
      java.nio.ByteOrder bo = java.nio.ByteOrder.BIG_ENDIAN;
      ArrayInputStream ais = new ArrayInputStream(fileName,bo);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      return null;
    }
  }

  private static class StructureTensors extends EigenTensors2 {
    StructureTensors(double sigma, float[][] x) {
      super(x[0].length,x.length);
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] u1 = new float[n2][n1];
      float[][] u2 = new float[n2][n1];
      float[][] su = new float[n2][n1];
      float[][] sv = new float[n2][n1];
      LocalOrientFilter lof = new LocalOrientFilter(sigma);
      lof.apply(x,null,u1,u2,null,null,su,sv,null);
      float[][] sc = Array.sub(1.0f,coherence(sigma,x));
      su = Array.mul(su,sc);
      sv = Array.mul(sv,sc);
      plot(su);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          setEigenvalues(i1,i2,su[i2][i1],sv[i2][i1]);
          setEigenvectorU(i1,i2,u1[i2][i1],u2[i2][i1]);
        }
      }
    }
  }

  private static float[][] coherence(double sigma, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    LocalOrientFilter lof1 = new LocalOrientFilter(sigma);
    LocalOrientFilter lof2 = new LocalOrientFilter(sigma*4);
    float[][] u11 = new float[n2][n1];
    float[][] u21 = new float[n2][n1];
    float[][] su1 = new float[n2][n1];
    float[][] sv1 = new float[n2][n1];
    float[][] u12 = new float[n2][n1];
    float[][] u22 = new float[n2][n1];
    float[][] su2 = new float[n2][n1];
    float[][] sv2 = new float[n2][n1];
    lof1.apply(x,null,u11,u21,null,null,su1,sv1,null);
    lof2.apply(x,null,u12,u22,null,null,su2,sv2,null);
    float[][] c = u11;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u11i = u11[i2][i1];
        float u21i = u21[i2][i1];
        float su1i = su1[i2][i1];
        float sv1i = sv1[i2][i1];
        float u12i = u12[i2][i1];
        float u22i = u22[i2][i1];
        float su2i = su2[i2][i1];
        float sv2i = sv2[i2][i1];
        float s111 = (su1i-sv1i)*u11i*u11i+sv1i;
        float s121 = (su1i-sv1i)*u11i*u21i     ;
        float s221 = (su1i-sv1i)*u21i*u21i+sv1i;
        float s112 = (su2i-sv2i)*u12i*u12i+sv2i;
        float s122 = (su2i-sv2i)*u12i*u22i     ;
        float s222 = (su2i-sv2i)*u22i*u22i+sv2i;
        float s113 = s111*s112+s121*s122;
        float s223 = s121*s122+s221*s222;
        float t1 = s111+s221;
        float t2 = s112+s222;
        float t3 = s113+s223;
        float t12 = t1*t2;
        c[i2][i1] = (t12>0.0f)?t3/t12:0.0f;
      }
    }
    return c;
  }

  private static float[][][] getTensorEllipses(
    int n1, int n2, int ns, EigenTensors2 et) 
  {
    int nt = 51;
    int m1 = (n1-1)/ns;
    int m2 = (n2-1)/ns;
    int j1 = (n1-1-(m1-1)*ns)/2;
    int j2 = (n2-1-(m2-1)*ns)/2;
    int nm = m1*m2;
    double r = 0.45*ns;
    float[][] x1 = new float[nm][nt];
    float[][] x2 = new float[nm][nt];
    double dt = 2.0*PI/(nt-1);
    double ft = 0.0f;
    for (int i2=j2,im=0; i2<n2; i2+=ns) {
      double y2 = i2+r;
      for (int i1=j1; i1<n1; i1+=ns,++im) {
        float[] u = et.getEigenvectorU(i1,i2);
        float[] s = et.getEigenvalues(i1,i2);
        double u1 = u[0];
        double u2 = u[1];
        double v1 = -u2;
        double v2 =  u1;
        double su = s[0];
        double sv = s[1];
        double a = r*sqrt(sv/su);
        double b = r;
        for (int it=0; it<nt; ++it) {
          double t = ft+it*dt;
          double cost = cos(t);
          double sint = sin(t);
          x1[im][it] = (float)(i1+a*cost*u1-b*sint*u2);
          x2[im][it] = (float)(i2+b*sint*u1+a*cost*u2);
        }
      }
    }
    return new float[][][]{x1,x2};
  }

  private static class SimpleTensors extends EigenTensors2 {
    SimpleTensors(int n1, int n2, double du, double dv, double v1) {
      super(n1,n2);
      float u2 = -(float)v1;
      float u1 = sqrt(1.0f-u2*u2);
      float au = (float)du;
      float av = (float)dv;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          setEigenvalues(i1,i2,au,av);
          setEigenvectorU(i1,i2,u1,u2);
        }
      }
    }
  }

  private static void testConstant() {
    int n1 = 301;
    int n2 = 301;
    int nv = 1;
    float du = 0.010f;
    float dv = 1.000f;
    float v1 = sin(-45.0f*FLT_PI/180.0f);
    SimpleTensors st = new SimpleTensors(n1,n2,du,dv,v1);
    Painting2 p = new Painting2(n1,n2,nv,st);
    p.paintAt(1*(n1-1)/4,1*(n2-1)/4,1.0f);
    p.paintAt(3*(n1-1)/4,3*(n2-1)/4,2.0f);
    p.extrapolate();
    //plotImageTensors(p.getTimes(),st);
    plot(p.getTimes(),ColorMap.PRISM);
    plot(p.getValues());
  }

  private static float[][] makeTargetImage(int n1, int n2) {
    float k = 0.3f;
    float[][] x = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float d2 = (float)(i2-n2/2);
      for (int i1=0; i1<n1; ++i1) {
        float d1 = (float)(i1-n1/2);
        x[i2][i1] = 10.0f*sin(k*sqrt(d1*d1+d2*d2));
      }
    }
    return x;
  }

  private static void testTarget() {
    int n1 = 315;
    int n2 = 315;
    int nv = 1;
    float[][] x = makeTargetImage(n1,n2);
    StructureTensors st = new StructureTensors(8,x);
    Painting2 p = new Painting2(n1,n2,nv,st);
    int m1 = 1;
    int m2 = 1;
    int nk = 1+(n2-1)/m2+1+(n1-1)/m1;
    int[] k1 = new int[nk];
    int[] k2 = new int[nk];
    float[] vk = new float[nk];
    int ik = 0;
    for (int i2=0; i2<n2; i2+=m2,++ik) {
      k1[ik] = n1/2;
      k2[ik] = i2;
      vk[ik] = x[k2[ik]][k1[ik]];
    }
    for (int i1=0; i1<n1; i1+=m1,++ik) {
      k1[ik] = i1;
      k2[ik] = n2/2;
      vk[ik] = x[k2[ik]][k1[ik]];
    }
    for (ik=0; ik<nk; ++ik)
      p.paintAt(k1[ik],k2[ik],vk[ik]);
    plotImageTensors(x,st);
    p.extrapolate();
    plot(p.getValues());
    p.interpolate();
    plot(p.getValues());
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testConstant();
        //testSeismic();
        //testChannels();
        //testTarget();
      }
    });
  }
}
