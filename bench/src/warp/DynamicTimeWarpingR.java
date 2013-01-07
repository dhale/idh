/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
//package edu.mines.jtk.dsp;
package warp;
import edu.mines.jtk.dsp.*;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic time warping.
 * <em>This fundamentally flawed version solves for r instead of u!</em>
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.10.31
 */
public class DynamicTimeWarpingR {

  public DynamicTimeWarpingR(Sampling st, double umin, double umax) {
    double fs = umin;
    double ds = 0.5*st.getDelta();
    int ns = 2+(int)((umax-umin)/ds);
    ds = (umax-umin)/(ns-1);
    _ss = new Sampling(ns,ds,fs);
    _st = st;
    _umin = (float)umin;
    _umax = (float)umax;
    _si = new SincInterp();
    _si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
  }

  public Sampling getShiftSampling() {
    return _ss;
  }

  public Sampling getTimeSampling() {
    return _st;
  }

  public float[][] computeErrors(
    Sampling sf, float[] f,
    Sampling sg, float[] g)
  {
    Sampling ss = _ss;
    Sampling st = _st;
    int nf = sf.getCount();
    int ng = sg.getCount();
    int nt = st.getCount();
    int ns = ss.getCount();
    float[][] e = new float[nt][ns];
    float[] fi = new float[nt];
    float[] gi = new float[nt];
    _si.interpolate(sf,f,st,fi);
    for (int is=0; is<ns; ++is) {
      double fts = st.getFirst()+ss.getValue(is);
      Sampling si = new Sampling(nt,st.getDelta(),fts);
      _si.interpolate(sg,g,si,gi);
      for (int it=0; it<nt; ++it)
        e[it][is] = error(fi[it],gi[it]);
    }
    return e;
  }

  public float[] findShifts(
    double u0, int nr, double rmin, double rmax, float[][] e) 
  {
    double dr = (nr>1)?(rmax-rmin)/(nr-1):(rmax-rmin);
    double fr = rmin;
    Sampling sr = new Sampling(nr,dr,fr);
    return findShifts(u0,sr,e);
  }

  public float[] findShifts(double u0, Sampling sr, float[][] e) {
    int nr = sr.getCount();
    int ns = _ss.getCount();
    int nt = _st.getCount();
    float dt = (float)_st.getDelta();
    double ds = _ss.getDelta();
    double fs = _ss.getFirst();
    float s0i = (float)u0;
    float e0i = _si.interpolate(ns,ds,fs,e[0],s0i);
    int[][] m = new int[nt][nr];
    float[][] s = new float[nt][nr];
    float[][] d = new float[nt][nr];
    for (int ir=0; ir<nr; ++ir) {
      s[0][ir] = s0i;
      d[0][ir] = e0i;
    }
    for (int it=1; it<nt; ++it) {
      for (int ir=0; ir<nr; ++ir) {
        float ridt = (float)sr.getValue(ir)*dt;
        int ir0 = ir;
        int irm = max(ir-1,0);
        int irp = min(ir+1,nr-1);
        float s0 = s[it-1][ir0]+ridt;
        float sm = s[it-1][irm]+ridt;
        float sp = s[it-1][irp]+ridt;
        s0 = max(_umin,min(_umax,s0));
        sm = max(_umin,min(_umax,sm));
        sp = max(_umin,min(_umax,sp));
        float e0 = _si.interpolate(ns,ds,fs,e[it],s0);
        float em = _si.interpolate(ns,ds,fs,e[it],sm);
        float ep = _si.interpolate(ns,ds,fs,e[it],sp);
        float d0 = d[it-1][ir0]+e0;
        float dm = d[it-1][irm]+em;
        float dp = d[it-1][irp]+ep;
        float dmin = min3(dm,d0,dp);
        float emin;
        float smin;
        int imin;
        if (dmin==d0) {
          smin = s0;
          imin = ir0;
        } else if (dmin==dm) {
          smin = sm;
          imin = irm;
        } else {
          smin = sp;
          imin = irp;
        }
        s[it][ir] = smin;
        m[it][ir] = imin;
        d[it][ir] = dmin;
      }
    }
    float[] u = new float[nt];
    int it = nt-1;
    int ir = (nr-1)/2;
    float di = d[it][ir];
    for (int jr=0; jr<nr; ++jr) {
      if (d[it][jr]<di) {
        ir = jr;
        di = d[it][jr];
      }
    }
    u[it] = s[it][ir];
    for (ir=m[it][ir],--it; it>=0; ir=m[it][ir],--it)
      u[it] = s[it][ir];
    return u;
  }

  public float[] applyShifts(
    Sampling sg, float[] g, 
    Sampling st, float[] u)
  {
    int ng = sg.getCount();
    int nt = st.getCount();
    float[] h = new float[nt];
    for (int it=0; it<nt; ++it) {
      double t = st.getValue(it)+u[it];
      h[it] = _si.interpolate(sg,g,t);
    }
    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private Sampling _ss,_st;
  private float _umin,_umax;
  private SincInterp _si;
  private float _epow = 2.0f;

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }
}
