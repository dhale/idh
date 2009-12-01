/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lss;

import java.awt.*;
import java.util.*;
import javax.swing.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Tests use of weighted semblance in velocity analysis.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.11.16
 */
public class Velan {

  public static float[][] velocitySpectrum(
    Sampling st, Sampling sx, float[][] p, 
    Sampling sv, double tsigma, boolean weighted)
  {
    int nv = sv.getCount();
    float[][] s = new float[nv][];
    RecursiveGaussianFilter tsmoother = new RecursiveGaussianFilter(tsigma);
    for (int iv=0; iv<nv; ++iv) {
      double v = sv.getValue(iv);
      float[][] q = nmo(v,st,sx,p);
      if (weighted)
        s[iv] = semblance(tsmoother,st,sx,v,q);
      else
        s[iv] = semblance(tsmoother,q);
    }
    return s;
  }

  public static float[] semblance(
    RecursiveGaussianFilter tsmoother, 
    Sampling st, Sampling sx, double vnmo, float[][] q)
  {
    int nx = q.length;
    int nt = q[0].length;
    float[] r = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        r[it] += q[ix][it];
      }
    }
    float[] arr = new float[nt];
    float[] arq = new float[nt];
    float[] aqq = new float[nt];
    float[] brr = new float[nt];
    float[] brq = new float[nt];
    float[] bqq = new float[nt];
    float gamma = (float)(1.0/vnmo*vnmo);
    for (int ix=0; ix<nx; ++ix) {
      float x = (float)sx.getValue(ix);
      float xx = x*x;
      float xxg = xx*gamma;
      for (int it=0; it<nt; ++it) {
        float t0 = (float)st.getValue(it);
        float ti = sqrt(t0*t0+xxg);
        float qi = q[ix][it];
        float ri = r[it];
        float ui = xx/ti;
        float rr = ri*ri;
        float rq = ri*qi;
        float qq = qi*qi;
        arr[it] += rr;
        arq[it] += rq;
        aqq[it] += qq;
        brr[it] += ui*rr;
        brq[it] += ui*rq;
        bqq[it] += ui*qq;
      }
    }
    tsmoother.apply0(arr,arr);
    tsmoother.apply0(arq,arq);
    tsmoother.apply0(aqq,aqq);
    tsmoother.apply0(brr,brr);
    tsmoother.apply0(brq,brq);
    tsmoother.apply0(bqq,bqq);
    float[] s = new float[nt];
    for (int it=0; it<nt; ++it) {
      double arri = arr[it];
      double arqi = arq[it];
      double aqqi = aqq[it];
      double brri = brr[it];
      double brqi = brq[it];
      double bqqi = bqq[it];
      double saab = arri*arqi*bqqi;
      double saba = arri*brqi*aqqi;
      double sbaa = brri*arqi*aqqi;
      double sabb = arri*brqi*bqqi;
      double sbab = brri*arqi*bqqi;
      double sbba = brri*brqi*aqqi;
      double b = 0.0;
      if (sabb<sbab && sbab<sbba || sbba<sbab && sbab<sabb) {
        double bnum = sbaa-2.0*saba+saab;
        double bden = sabb-2.0*sbab+sbba+bnum;
        b = (bden!=0.0)?bnum/bden:1.0;
      } else {
        double bnum = arqi;
        double bden = -brqi+bnum;
        b = (bden!=0.0)?bnum/bden:0.0;
      }
      if (b<0.0 || b>1.0) {
        double snuma = arqi*arqi;
        double sdena = arri*aqqi;
        double sa = (sdena>0.0)?(float)(snuma/sdena):0.0f;
        double snumb = brqi*brqi;
        double sdenb = brri*bqqi;
        double sb = (sdenb>0.0)?(float)(snumb/sdenb):0.0f;
        s[it] = (float)min(sa,sb);
      } else {
        double a = 1.0-b;
        double srri = a*arri+b*brri;
        double srqi = a*arqi+b*brqi;
        double sqqi = a*aqqi+b*bqqi;
        double snum = srqi*srqi;
        double sden = srri*sqqi;
        s[it] = (sden>0.0)?(float)(snum/sden):0.0f;
      }
    }
    return s;
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static float[] semblance(
    RecursiveGaussianFilter tsmoother, float[][] q)
  {
    int nx = q.length;
    int nt = q[0].length;
    float[] sn = new float[nt];
    float[] sd = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        float qi = q[ix][it];
        sn[it] += qi;
        sd[it] += qi*qi;
      }
    }
    mul(sn,sn,sn);
    tsmoother.apply0(sn,sn);
    tsmoother.apply0(sd,sd);
    float[] s = sn;
    for (int it=0; it<nt; ++it) {
      s[it] = sn[it]/(nx*sd[it]);
    }
    return s;
  }

  public static float[][] nmo(
    double vnmo, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(nt,dt,ft);
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
      double xxg = (x*x)/(vnmo*vnmo);
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        t[it] = (float)sqrt(t0*t0+xxg); 
      }
      si.setUniformSamples(p[ix]);
      si.interpolate(nt,t,q[ix]);
    }
    return q;
  }

  public static float[][] nmo(
    float[] vnmo, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(nt,dt,ft);
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        double v0 = vnmo[it];
        double xxg = (x*x)/(v0*v0);
        t[it] = (float)sqrt(t0*t0+xxg); 
      }
      si.setUniformSamples(p[ix]);
      si.interpolate(nt,t,q[ix]);
    }
    return q;
  }

  public static float[][] makeRickerGather(
    double fpeak, float[] vnmo, Sampling st, Sampling sx) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    //Random random = new Random(314159);
    Random random = new Random();
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    for (int jt=0; jt<nt; ++jt) {
      if (random.nextDouble()<0.1) // if reflection, ...
        continue;
      double t0 = st.getValue(jt);
      double a0 = 2.0*random.nextDouble()-1.0;
      double v0 = vnmo[st.indexOfNearest(t0)];
      double gamma = 1.0/(v0*v0);
      for (int ix=0; ix<nx; ++ix) {
        double x = (float)sx.getValue(ix);
        double t = sqrt(t0*t0+x*x*gamma);
        int itlo = max(0,(int)((t-thalf-ft)/dt));
        int ithi = min(nt-1,(int)((t+thalf-ft)/dt));
        for (int it=itlo; it<=ithi; ++it) {
          double twave = st.getValue(it)-t;
          p[ix][it] += (float)(a0*ricker(fpeak,twave));
        }
      }
    }
    return p;
  }

  private static float[] makeLinearVelocity(
    double vmin, double vmax, Sampling st) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    float[] v = new float[nt];
    for (int it=0; it<nt; ++it)
      v[it] = (float)((it*vmax+(nt-1-it)*vmin)/(nt-1));
    return v;
  }

  private static float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }

  private static void testGather() {
    Sampling st = new Sampling(1001,0.004,0.000); 
    Sampling sx = new Sampling(50,0.050,0.050);
    Sampling sv = new Sampling(101,0.020,1.5,2.5);
    double tsigma = 4.0;
    double fpeak = 25.0;
    float[] vp = makeLinearVelocity(2.00,3.00,st);
    float[][] p = makeRickerGather(fpeak,vp,st,sx);
    float[] vm = makeLinearVelocity(1.98,2.70,st);
    p = add(p,makeRickerGather(fpeak,vm,st,sx));
    float[][] q = nmo(vp,st,sx,p);
    SimplePlot spp = SimplePlot.asPixels(st,sx,p);
    spp.setHLabel("Offset (km)");
    spp.setVLabel("Time (s)");
    spp.setSize(400,900);
    SimplePlot spq = SimplePlot.asPixels(st,sx,q);
    spq.setHLabel("Offset (km)");
    spq.setVLabel("Time (s)");
    spq.setSize(400,900);
    for (boolean weighted:new boolean[]{true,false}) {
      float[][] s = velocitySpectrum(st,sx,p,sv,tsigma,weighted);
      System.out.println("s min="+min(s)+" max="+max(s));
      SimplePlot spv = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
      spv.setHLabel("Velocity (km/s)");
      spv.setVLabel("Time (km/s)");
      PixelsView pv = spv.addPixels(st,sv,s);
      pv.setColorModel(ColorMap.JET);
      pv.setInterpolation(PixelsView.Interpolation.NEAREST);
      pv.setClips(0.0f,1.0f);
      ContoursView cv = spv.addContours(st,sv,s);
      cv.setLineColor(Color.BLACK);
      cv.setContours(new Sampling(4,0.2,0.2));
      spv.setSize(400,900);
    }
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testGather();
      }
    });
  }
}
