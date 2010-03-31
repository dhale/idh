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
    return velocitySpectrum(st,sx,p,sv,tsigma,weighted,null);
  }
  public static float[][] velocitySpectrum(
    Sampling st, Sampling sx, float[][] p, 
    Sampling sv, double tsigma, boolean weighted,
    float[][] b)
  {
    int nv = sv.getCount();
    float[][] s = new float[nv][];
    for (int iv=0; iv<nv; ++iv) {
      double v = sv.getValue(iv);
      float[][] q = nmo(v,st,sx,p);
      if (weighted) {
        float[] biv = (b!=null)?b[iv]:null;
        s[iv] = semblance(st,sx,v,tsigma,q,biv);
      } else {
        s[iv] = semblance(tsigma,q);
      }
    }
    return s;
  }

  public static float[] semblance(
    Sampling st, Sampling sx, double vnmo, double tsigma, float[][] q)
  {
    return semblance(st,sx,vnmo,tsigma,q,null);
  }
  public static float[] semblance(
    Sampling st, Sampling sx, double vnmo, double tsigma, float[][] q,
    float[] bs)
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
    float gamma = (float)(1.0/(vnmo*vnmo));
    for (int ix=0; ix<nx; ++ix) {
      float x = (float)sx.getValue(ix);
      float xx = x*x;
      float xxg = xx*gamma;
      for (int it=0; it<nt; ++it) {
        float t0 = (float)st.getValue(it);
        float ti = sqrt(t0*t0+xxg);
        float qi = q[ix][it];
        if (qi==0.0f) continue;
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
    esmooth(tsigma,arr,arr);
    esmooth(tsigma,arq,arq);
    esmooth(tsigma,aqq,aqq);
    esmooth(tsigma,brr,brr);
    esmooth(tsigma,brq,brq);
    esmooth(tsigma,bqq,bqq);
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
        if (sa<=sb) {
          s[it] = (float)sa;
          if (bs!=null) bs[it] = 0.0f;
        } else {
          s[it] = (float)sb;
          if (bs!=null) bs[it] = 1.0f;
        }
      } else {
        double a = 1.0-b;
        double srri = a*arri+b*brri;
        double srqi = a*arqi+b*brqi;
        double sqqi = a*aqqi+b*bqqi;
        double snum = srqi*srqi;
        double sden = srri*sqqi;
        s[it] = (sden>0.0)?(float)(snum/sden):0.0f;
        if (bs!=null) bs[it] = (float)b;
      }
    }
    return s;
  }

  public static float[] semblance(double tsigma, float[][] q) {
    int nx = q.length;
    int nt = q[0].length;
    float[] sn = new float[nt];
    float[] sd = new float[nt];
    float[] sx = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        float qi = q[ix][it];
        if (qi==0.0f) continue;
        sn[it] += qi;
        sd[it] += qi*qi;
        sx[it] += 1.0f;
      }
    }
    mul(sn,sn,sn);
    esmooth(tsigma,sn,sn);
    esmooth(tsigma,sd,sd);
    float[] s = sn;
    for (int it=0; it<nt; ++it) {
      s[it] = (sd[it]>0.0)?sn[it]/(sx[it]*sd[it]):0.0f;
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
    Random random = new Random(314159);
    //Random random = new Random();
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

  public static float[][] addRandomNoise(float r, float[][] p) {
    int nt = p[0].length;
    int nx = p.length;
    //float pmax = max(abs(p)); // peak signal
    float prms = sqrt(sum(mul(p,p))/nt/nx); // rms of signal
    Random random = new Random(3);
    float[][] s = sub(randfloat(random,nt,nx),0.5f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply0X(s,s); // noise, bandlimited in time only
    float srms = sqrt(sum(mul(s,s))/nt/nx); // rms of noise
    return add(mul(prms/(srms*r),s),p); // r = rms-signal / rms-noise
  }

  public static float[][] makeSimpleGather(
    double fpeak, float[] vnmo, Sampling st, Sampling sx) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    int kt = nt/8; // spacing between reflections, in samples
    for (int jt=kt; jt<nt-1; jt+=kt) {
      double t0 = st.getValue(jt);
      double a0 = 1.0; // constant amplitude
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

  /**
   * Two-sided exponential smoothing with zero-slope boundary conditions.
   * Zero-slope means that values off the ends of the input array are 
   * assumed to equal the values at the ends. This boundary condition
   * ensures that this filter does nothing to an array of constant values.
   * <p>
   * Input and output arrays may the same array.
   * @param sigma filter half-width, in samples; approximating a Gaussian.
   * @param p array of input samples.
   * @param q array of output samples.
   */
  public static void esmooth(double sigma, float[] p, float[] q) {
    if (p==q) p = copy(p);
    float sigmas = (float)(sigma*sigma);
    float a = (sigmas>0.0)?(1.0f+sigmas-sqrt(1.0f+2.0f*sigmas))/sigmas:0.0f;
    float b = (1.0f-a)/(1.0f+a);
    int n = p.length;
    float qi = b/(1.0f-a)*p[0];
    q[0] = qi;
    for (int i=1; i<n; ++i) {
      qi = a*qi+b*p[i];
      q[i] = qi;
    }
    qi = a*b/(1.0f-a)*p[n-1];
    q[n-1] += qi;
    for (int i=n-2; i>=0; --i) {
      qi = a*(qi+b*p[i+1]);
      q[i] += qi;
    }
  }

  private static void testGather() {
    Sampling st = new Sampling(1001,0.004,0.000); 
    Sampling sx = new Sampling(60,0.050,0.050);
    Sampling sv = new Sampling(101,0.020,1.5,2.5);
    int nv = sv.getCount();
    int nt = st.getCount();
    double tsigma = 5.0;
    double fpeak = 25.0;
    double[] vps = {2.00,3.00};
    double[] vms = {1.98,2.65}; // lower velocities for multiples
    double[] ts =  {0.00,4.00};
    float[] vp = makeLinearVelocity(vps[0],vps[1],st);
    //float[][] p = makeSimpleGather(fpeak,vp,st,sx);
    float[][] p = makeRickerGather(fpeak,vp,st,sx);
    float[] vm = makeLinearVelocity(vms[0],vms[1],st);
    p = add(p,makeRickerGather(fpeak,vm,st,sx));
    float snratio = 1.0e6f; // = rms-signal / rms-noise
    p = addRandomNoise(snratio,p); // add noise
    float[][] q = nmo(vp,st,sx,p);
    SimplePlot spp = SimplePlot.asPixels(st,sx,p);
    spp.setVLimits(1.0,3.0);
    spp.setHLabel("Offset (km)");
    spp.setVLabel("Time (s)");
    spp.setSize(400,400);
    /*
    SimplePlot spq = SimplePlot.asPixels(st,sx,q);
    spq.setVLimits(1.0,3.0);
    spq.setHLabel("Offset (km)");
    spq.setVLabel("Time (s)");
    spq.setSize(400,400);
    */
    for (boolean weighted:new boolean[]{true,false}) {
      float[][] b = weighted?new float[nv][nt]:null;
      float[][] s = velocitySpectrum(st,sx,p,sv,tsigma,weighted,b);
      System.out.println("s min="+min(s)+" max="+max(s));
      SimplePlot spv = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
      spv.setHLimits(2.0,3.0);
      spv.setVLimits(1.0,3.0);
      spv.setHLabel("Velocity (km/s)");
      spv.setVLabel("Time (km/s)");
      PixelsView pv = spv.addPixels(st,sv,s);
      pv.setColorModel(ColorMap.JET);
      pv.setInterpolation(PixelsView.Interpolation.LINEAR);
      pv.setClips(0.0f,1.0f);
      ContoursView cv = spv.addContours(st,sv,s);
      cv.setLineColor(Color.BLACK);
      cv.setLineWidth(2.0f);
      cv.setContours(new Sampling(1,0.5,0.2));
      PointsView tvp = spv.addPoints(ts,vps);
      tvp.setLineColor(Color.BLACK);
      tvp.setLineWidth(2.0f);
      PointsView tvm = spv.addPoints(ts,vms);
      tvm.setLineColor(Color.BLACK);
      tvm.setLineWidth(2.0f);
      spv.setSize(400,400);
      if (b!=null) {
        SimplePlot spb = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
        spb.setHLimits(2.0,3.0);
        spb.setVLimits(1.0,3.0);
        spb.setHLabel("Velocity (km/s)");
        spb.setVLabel("Time (km/s)");
        PixelsView pvb = spb.addPixels(st,sv,b);
        pvb.setColorModel(ColorMap.JET);
        pvb.setInterpolation(PixelsView.Interpolation.LINEAR);
        pvb.setClips(0.0f,1.0f);
        PointsView tvpb = spb.addPoints(ts,vps);
        tvpb.setLineColor(Color.BLACK);
        tvpb.setLineWidth(2.0f);
        PointsView tvmb = spb.addPoints(ts,vms);
        tvmb.setLineColor(Color.BLACK);
        tvmb.setLineWidth(2.0f);
        spb.setSize(400,400);
      }
    }
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testGather();
      }
    });
  }
  private static void trace(String s) {
    System.out.println(s);
  }
}
