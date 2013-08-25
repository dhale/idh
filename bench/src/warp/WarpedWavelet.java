/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Warping without wavelet distortion.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.24
 */
public class WarpedWavelet {

  public interface Warping {
    public void apply(float[] x, float[] y);
  }

  public static class Nmo implements Warping {
    public Nmo(Sampling st, double offset, double vnmo) {
      this(st,offset,fillfloat((float)vnmo,st.getCount()));
    }
    public Nmo(Sampling st, double offset, float[] vnmo) {
      int nt = st.getCount();
      double dt = st.getDelta();
      double ft = st.getFirst();
      _st = st;
      _si = SincInterp.fromErrorAndFrequency(0.01,0.45);
      _ty = new float[nt];
      _ay = new float[nt];
      for (int it=0; it<nt; ++it) {
        double tx = ft+it*dt;
        _ty[it] = (float)sqrt(tx*tx+(offset*offset)/(vnmo[it]*vnmo[it]));
        _ay[it] = (float)((_ty[it]>0.0)?tx/_ty[it]:0.0);
      }
    }
    public void apply(float[] x, float[] y) {
      int nt = _st.getCount();
      double dt = _st.getDelta();
      double ft = _st.getFirst();
      _si.interpolate(nt,dt,ft,x,nt,_ty,y);
      for (int it=0; it<nt; ++it)
        y[it] *= _ay[it];
    }
    private Sampling _st;
    private SincInterp _si;
    private float[] _ty;
    private float[] _ay;
  }

  public WarpedWavelet(Warping w) {
    _w = w;
  }

  public float[] estimateInverse(int na, int ka, float[] x, float[] y) {

    // Temporary array holds shifted x[it-lag].
    int nt = x.length;
    float[] t = new float[nt];

    // Array of z = warped(x) - y for all lag = ka,ka+1,...,ka+na-1.
    float[][] z = new float[na][nt];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      int itlo = max(0,lag);   // 0 <= it-lag
      int ithi = min(nt,nt+lag); // it-lag < nt
      for (int it=0; it<itlo; ++it)
        t[it] = 0.0f;
      for (int it=itlo; it<ithi; ++it)
        t[it] = x[it-lag];
      for (int it=ithi; it<nt; ++it)
        t[it] = 0.0f;
      _w.apply(t,z[ia]);
      for (int it=itlo; it<ithi; ++it)
        z[ia][it] -= y[it-lag];
    }

    // The matrix R and right-hand-side vector b, for Ra = b. Because
    // a0 = a[-ka] = 1, only the na-1 coefficients of a for non-zero
    // lag are unknown.
    int ma = na-1;
    DMatrix r = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ir=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jr=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double rij = 0.0;
        for (int it=0; it<nt; ++it) {
          rij += z[ia][it]*z[ja][it];
        }
        r.set(ir,jr,rij);
        ++jr;
      }
      double bi = 0.0;
      for (int it=0; it<nt; ++it) {
        bi -= z[ia][it]*z[-ka][it];
      }
      b.set(ir,0,bi);
      ++ir;
    }
    //System.out.println("r=\n"+r);
    //System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of R.
    DMatrixChd chd = new DMatrixChd(r);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ir=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // a0 = 1
      } else {
        aa[ia] = (float)a.get(ir,0);
        ++ir;
      }
    }
    return aa;
  }

  private Warping _w;
}
