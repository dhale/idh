package tp;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import model.SeismicModel1D;

/**
 * Synthetic seismograms.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.03.28
 */
public class SynSeis {

  public static class Model {
    public double x1,x2,x3; // CSM coordinates of first sample
    public int nz; // number of sampled depths
    public Sampling sz; // sampling of depth, in km
    public float v0; // velocity for depths z < first sampled z
    public float[] v; // velocities, km/s
    public float[] d; // densities, gm/cc
    public float[] a; // acoustic impedance
    public float[] r; // reflectivity
    public float[] t; // two-way time
    public float tmin() {
      return t[0];
    }
    public float tmax() {
      return t[nz-1];
    }
  }

  public static class RickerWavelet {
    public RickerWavelet(double fpeak) {
      _a = PI*fpeak;
    }
    public float getValue(double t) {
      double b = _a*t;
      double c = b*b;
      return (float)((1.0-2.0*c)*exp(-c));
    }
    public float getWidth() {
      return (float)(6.0/_a);
    }
    private double _a;
  }

  public static float[] makeBetterSeismogram(
    Model model, double q, double fpeak, Sampling st)
  {
    Sampling sz = model.sz;
    int nz = sz.getCount();
    SeismicModel1D sm = new SeismicModel1D();
    sm.setSourceType(SeismicModel1D.SourceType.LAND_VIBROSEIS);
    sm.setSensorType(SeismicModel1D.SensorType.GEOPHONE);
    sm.setSurfaceReflectionCoefficient(1.0);
    sm.setRickerWavelet(fpeak);
    sm.addLayer(0.0,model.v0,2.0,q);
    sm.addSource(0.0,1.0);
    sm.addSensor(0.0);
    for (int iz=0; iz<nz; ++iz) {
      double zi = sz.getValue(iz);
      double vi = model.v[iz];
      double di = model.d[iz];
      sm.addLayer(zi,vi,di,q);
    }
    //sm.dumpLayers();
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int ntpad = nt+(int)(ft/dt);
    double fref = 10000.0;
    System.out.println("begin make ...");
    float[] f = sm.makeSeismograms(ntpad,dt,fref)[0];
    System.out.println("... done");
    float[] g = new float[nt];
    SincInterp si = new SincInterp();
    si.interpolate(ntpad,dt,0.0,f,nt,dt,ft,g);
    g = taperEnds(fpeak,st,g);
    return g;
  }

  public static float[] makeSimpleSeismogram(
    Model model, double fpeak, Sampling st)
  {
    Sampling sz = model.sz;
    int nz = sz.getCount();
    int nt = st.getCount();
    float dt = (float)(st.getDelta());
    float ft = (float)(st.getFirst());
    float[] f = new float[nt];
    float[] v = model.v;
    float[] r = model.r;
    float[] t = model.t;
    RickerWavelet rw = new RickerWavelet(fpeak);
    int kw = 1+(int)(0.5*rw.getWidth()/dt);
    for (int iz=0; iz<nz; ++iz) {
      float rz = r[iz];
      float tz = t[iz];
      int it = st.indexOfNearest(tz);
      int jtlo = max(0,it-kw);
      int jthi = min(nt,it+kw);
      for (int jt=jtlo; jt<jthi; ++jt) {
        float tj = ft+jt*dt;
        f[jt] += rz*rw.getValue(tj-tz);
      }
    }
    f = taperEnds(fpeak,st,f);
    return f;
  }

  private static float[] taperEnds(double fpeak, Sampling st, float[] f) {
    f = copy(f);
    int nt = st.getCount();
    double dt = st.getDelta();
    RickerWavelet rw = new RickerWavelet(fpeak);
    int kw = 1+(int)(0.5*rw.getWidth()/dt);
    for (int it=0; it<kw; ++it) {
      float w = (float)(0.5*(1.0-cos(PI*it/kw)));
      f[it] *= w;
      f[nt-1-it] *= w;
    }
    return f;
  }

  public static Model getModel(WellLog log) {

    // If depth sampling not uniform, then model is null.
    int nz = log.z.length;
    double[] z = new double[nz];
    for (int iz=0; iz<nz; ++iz)
      z[iz] = log.z[iz];
    Sampling sz = new Sampling(z);
    if (!sz.isUniform())
      return null;

    // Velocities and densities from log.
    float[] v = log.v;
    float[] d = log.d;

    // Index of first non-null velocity and density.
    int jz = 0;
    while (jz<nz && (v[jz]==WellLog.NULL_VALUE || d[jz]==WellLog.NULL_VALUE))
      ++jz;

    // Index of last non-null velocity and density.
    int lz = nz-1;
    while (lz>=0 && (v[lz]==WellLog.NULL_VALUE || d[lz]==WellLog.NULL_VALUE))
      --lz;

    // If any null values between first and last, then model is null.
    for (int iz=jz; iz<=lz; ++iz) {
      if (v[iz]==WellLog.NULL_VALUE || d[iz]==WellLog.NULL_VALUE)
        return null;
    }

    // If insufficient samples, then model is null.
    nz = 1+lz-jz;
    if (nz<=0)
      return null;

    // Uniformly sampled velocity, density, impedance and reflectivity.
    v = copy(nz,jz,v);
    d = copy(nz,jz,d);
    float[] a = new float[nz];
    float[] r = new float[nz];
    a[0] = v[0]*d[0];
    r[0] = 0.0f;
    for (int mz=0,iz=1; iz<nz; ++mz,++iz) {
      float ai = v[iz]*d[iz];
      float am = v[mz]*d[mz];
      a[iz] = ai;
      r[iz] = (am-ai)/(am+ai);
    }

    // Velocity for depths less than those sampled.
    float v0 = v[0]; // can we do better?

    // Sampling of depth in km for non-null values; 1 ft = 0.0003048 km
    float dz = (float)(sz.getDelta()*0.0003048);
    float fz = log.x1[jz];

    // Two-way time, a sampled function depth.
    float[] t = new float[nz];
    float tz = 2.0f*fz/v0;
    for (int iz=0; iz<nz; ++iz) {
      t[iz] = tz;
      tz += 2.0f*dz/v[iz];
    }

    // The model.
    Model m = new Model();
    m.x1 = log.x1[jz];
    m.x2 = log.x2[jz];
    m.x3 = log.x3[jz];
    m.nz = nz;
    m.sz = new Sampling(nz,dz,fz);
    m.v0 = v0;
    m.v = v;
    m.d = d;
    m.a = a;
    m.r = r;
    m.t = t;
    return m;
  }
}
