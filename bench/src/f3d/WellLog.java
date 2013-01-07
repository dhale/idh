package f3d;

import java.io.*;
import static java.lang.Math.abs;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

import edu.mines.jtk.mosaic.*;

/**
 * A well log from F3. 
 * Note that the coordinates x1 for these logs are in time, not depth.
 * For F3 we have only four wells with logs, and insufficient information
 * is available to accurately convert the seismic images from time to depth.
 * However, time-depth functions are provided for each of the four logs, so
 * we convert the logs to time.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.12.29
 */
public class WellLog {

  public String name; // unique well-name
  public double xe,yn,ze; // easting, northing, elevation
  public int n; // number of samples
  public float[] z; // measured depth (along wellbore)
  public float[] v,d,g,p; // velocity, density, gamma, porosity; may be null
  public float[] x1,x2,x3; // CSM resampled coordinates of well-bore (s,km)

  // This value represents missing values in all well curve data.
  public static final float NULL_VALUE = -999.2500f;

  /**
   * A collection of well log data.
   */
  public static class Data {

    /**
     * Constructs empty well log data.
     */
    public Data() {
    }

    /**
     * Constructs well log data from specified files.
     * @param wlDirName directory containing LAS files.
     */
    public Data(String wlDirName) {
      File wlDir = new File(wlDirName);
      load(wlDir,this);
      System.out.println("WellLog.Data: after initial loading, ...");
      System.out.println("  number of logs = "+size());
    }

    /**
     * Writes these well log data to a binary file.
     * @param fileName the file name.
     */
    public void writeBinary(String fileName) {
      try {
        ArrayOutputStream aos = new ArrayOutputStream(fileName);
        aos.writeInt(size());
        for (WellLog log:getAll()) {
          aos.writeUTF(log.name);
          aos.writeDouble(log.xe);
          aos.writeDouble(log.yn);
          aos.writeDouble(log.ze);
          aos.writeInt(log.n);
          int vdgp = 0;
          if (log.v!=null) vdgp |= 1;
          if (log.d!=null) vdgp |= 2;
          if (log.g!=null) vdgp |= 4;
          if (log.p!=null) vdgp |= 8;
          aos.writeInt(vdgp);
          aos.writeFloats(log.z);
          if (log.v!=null) aos.writeFloats(log.v);
          if (log.d!=null) aos.writeFloats(log.d);
          if (log.g!=null) aos.writeFloats(log.g);
          if (log.p!=null) aos.writeFloats(log.p);
          aos.writeFloats(log.x1);
          aos.writeFloats(log.x2);
          aos.writeFloats(log.x3);
        }
        aos.close();
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }
    
    /**
     * Reads well log data from a binary file.
     * @param fileName the file name.
     */
    public static Data readBinary(String fileName) {
      Data data = new Data();
      try {
        ArrayInputStream ais = new ArrayInputStream(fileName);
        int nlog = ais.readInt();
        for (int ilog=0; ilog<nlog; ++ilog) {
          WellLog log = new WellLog();
          log.name = ais.readUTF();
          log.xe = ais.readDouble();
          log.yn = ais.readDouble();
          log.ze = ais.readDouble();
          log.n = ais.readInt();
          int vdgp = ais.readInt();
          int n = log.n;
          log.z = new float[n];
          log.v = ((vdgp&1)!=0)?new float[n]:null;
          log.d = ((vdgp&2)!=0)?new float[n]:null;
          log.g = ((vdgp&4)!=0)?new float[n]:null;
          log.p = ((vdgp&8)!=0)?new float[n]:null;
          ais.readFloats(log.z);
          if (log.v!=null) ais.readFloats(log.v);
          if (log.d!=null) ais.readFloats(log.d);
          if (log.g!=null) ais.readFloats(log.g);
          if (log.p!=null) ais.readFloats(log.p);
          log.x1 = new float[n];
          log.x2 = new float[n];
          log.x3 = new float[n];
          ais.readFloats(log.x1);
          ais.readFloats(log.x2);
          ais.readFloats(log.x3);
          data.add(log);
        }
        ais.close();
        return data;
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

    /**
     * Adds the specified well log.
     * @param log the well log.
     */
    public void add(WellLog log) {
      _data.put(log.name,log);
    }

    /**
     * Returns the number of well logs.
     * @return the number of well logs.
     */
    public int size() {
      return _data.size();
    }

    /**
     * Clips these well log data by removing all logs outside bounds.
     * @param s2 sampling for bounds in 2nd dimension
     * @param s3 sampling for bounds in 3rd dimension
     */
    public void clip(Sampling s2, Sampling s3) {
      float f2 = (float)s2.getFirst();
      float f3 = (float)s3.getFirst();
      float l2 = (float)s2.getLast();
      float l3 = (float)s3.getLast();
      Map<String,WellLog> data = new HashMap<String,WellLog>();
      for (String name:_data.keySet()) {
        WellLog log = _data.get(name);
        float[] x2 = log.x2;
        float[] x3 = log.x3;
        int n = log.n;
        boolean ok = true;
        for (int i=0; i<n && ok; ++i)
          ok = f2<=x2[i] && x2[i]<=l2 && f3<=x3[i] && x3[i]<=l3;
        if (ok)
          data.put(name,log);
      }
      _data = data;
    }

    /**
     * Gets the well log for the specified well name.
     * @param name the well name.
     * @return the well log; null, if none.
     */
    public WellLog get(String name) {
      return _data.get(name);
    }

    /**
     * Gets all well logs.
     * @return list of logs.
     */
    public List<WellLog> getAll() {
      return new ArrayList<WellLog>(_data.values());
    }

    /**
     * Gets well logs with the specified curve.
     * @param curve the well log curve.
     * @return list of logs.
     */
    public List<WellLog> getLogsWith(String curve) {
      List<WellLog> list = new ArrayList<WellLog>();
      for (WellLog log:_data.values()) {
        if (log.getCurve(curve)!=null)
          list.add(log);
      }
      return list;
    }

    /**
     * Gets intersections of boreholes with constant-x1 plane.
     * @return array {x2,x3} of arrays of (x2,x3) coordinates.
     */
    public float[][] getIntersections(String curve, float x1) {
      FloatList x2l = new FloatList();
      FloatList x3l = new FloatList();
      for (WellLog log:getAll()) {
        float[] c = log.getCurve(curve);
        if (c!=null) {
          int n = log.n;
          float[] x1s = log.x1;
          float[] x2s = log.x2;
          float[] x3s = log.x3;
          for (int i=1; i<n; ++i) {
            if ((x1-x1s[i])*(x1-x1s[i-1])<=0.0) {
              float a = abs(x1-x1s[i-1])/abs(x1s[i]-x1s[i-1]);
              float x2 = a*x2s[i]+(1.0f-a)*x2s[i-1];
              float x3 = a*x3s[i]+(1.0f-a)*x3s[i-1];
              x2l.add(x2);
              x3l.add(x3);
            }
          }
        }
      }
      float x2[] = x2l.trim();
      float x3[] = x3l.trim();
      return new float[][]{x2,x3};
    }

    /**
     * Rasterizes all logs that have the specified curve.
     * This method assumes that all well log curves are sampled
     * more finely than the returned raster image. Sample values
     * in the image will be set to the average of all curve
     * values that are nearest to that image sample. Image samples
     * that are nearest to no curve values are set to the specified
     * null value.
     * in the raster image
     * @param curve the well log curve to rasterize.
     * @param fnull default raster image sample value.
     * @param s1 sampling in 1st dimension of raster image.
     * @param s2 sampling in 2nd dimension of raster image.
     * @param s3 sampling in 3rd dimension of raster image.
     * @return the sampled raster image.
     */
    public float[][][] rasterizeLogsWith(
      String curve, float fnull, 
      Sampling s1, Sampling s2, Sampling s3) 
    {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      double fx1 = s1.getFirst();
      double fx2 = s2.getFirst();
      double fx3 = s3.getFirst();
      double lx1 = s1.getLast();
      double lx2 = s2.getLast();
      double lx3 = s3.getLast();
      double sx1 = 1.0/s1.getDelta();
      double sx2 = 1.0/s2.getDelta();
      double sx3 = 1.0/s3.getDelta();
      float[][][] image = new float[n3][n2][n1];
      float[][][] count = new float[n3][n2][n1];
      for (WellLog log:getAll()) {
        float[] c = log.getCurve(curve);
        if (c!=null) {
          int n = log.n;
          float[] x1 = log.x1;
          float[] x2 = log.x2;
          float[] x3 = log.x3;
          for (int i=0; i<n; ++i) {
            if (c[i]!=NULL_VALUE) {
              double x1i = x1[i]; if (x1i<fx1 || lx1<x1i) continue;
              double x2i = x2[i]; if (x2i<fx2 || lx2<x2i) continue;
              double x3i = x3[i]; if (x3i<fx3 || lx3<x3i) continue;
              int i1 = (int)((x1i-fx1)*sx1+0.5);
              int i2 = (int)((x2i-fx2)*sx2+0.5);
              int i3 = (int)((x3i-fx3)*sx3+0.5);
              image[i3][i2][i1] += c[i];
              count[i3][i2][i1] += 1.0f;
            }
          }
        }
      }
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            if (count[i3][i2][i1]>0.0f) {
              image[i3][i2][i1] /= count[i3][i2][i1];
            } else {
              image[i3][i2][i1] = fnull;
            }
          }
        }
      }
      return image;
    }

    public void printCounts(float[][][] image, float fnull) {
      int n1 = image[0][0].length;
      int n2 = image[0].length;
      int n3 = image.length;
      int n3max = 0; 
      int i3max = -1;
      for (int i3=0; i3<n3; ++i3) {
        int ns = 0;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            if (image[i3][i2][i1]!=fnull)
              ++ns;
          }
        }
        System.out.println("i3="+i3+" ns="+ns);
        if (ns>n3max) { n3max = ns; i3max = i3; }
      }
      System.out.println("i3max="+i3max+" n3max="+n3max);
      int n2max = 0; 
      int i2max = -1;
      for (int i2=0; i2<n2; ++i2) {
        int ns = 0;
        for (int i3=0; i3<n3; ++i3) {
          for (int i1=0; i1<n1; ++i1) {
            if (image[i3][i2][i1]!=fnull)
              ++ns;
          }
        }
        System.out.println("i2="+i2+" ns="+ns);
        if (ns>n2max) { n2max = ns; i2max = i2; }
      }
      System.out.println("i2max="+i2max+" n2max="+n2max);
      int n1max = 0; 
      int i1max = -1;
      for (int i1=0; i1<n1; ++i1) {
        int ns = 0;
        for (int i3=0; i3<n3; ++i3) {
          for (int i2=0; i2<n2; ++i2) {
            if (image[i3][i2][i1]!=fnull)
              ++ns;
          }
        }
        System.out.println("i1="+i1+" ns="+ns);
        if (ns>n1max) { n1max = ns; i1max = i1; }
      }
      System.out.println("i1max="+i1max+" n1max="+n1max);
    }

    /**
     * Prints summary information for these well log data.
     */
    public void printInfo() {
      int nv = 0;
      int nd = 0;
      int ng = 0;
      int np = 0;
      for (WellLog log:getAll()) {
        if (log.v!=null) ++nv;
        if (log.d!=null) ++nd;
        if (log.g!=null) ++ng;
        if (log.p!=null) ++np;
      }
      System.out.println("number of logs = "+_data.size());
      System.out.println("number of velocity logs   = "+nv);
      System.out.println("number of density logs    = "+nd);
      System.out.println("number of gamma ray logs  = "+ng);
      System.out.println("number of porosity logs   = "+np);
    }

    private Map<String,WellLog> _data = new TreeMap<String,WellLog>();
    private static void load(File wlDir, WellLog.Data wldata) {
      Check.argument(wlDir.isDirectory(),wlDir+" is a directory");
      for (File file:wlDir.listFiles()) {
        if (file.isDirectory()) {
          load(file,wldata);
        } else if (file.isFile()) {
          WellLog log = WellLog.load(file);
          if (log!=null)
            wldata.add(log);
        }
      }
    }
  }

  /**
   * Loads one well log from the specified LAS file, if possible.
   * @param file the LAS file.
   * @return the well log; null, if the file format is not understood.
   */
  public static WellLog load(File file) {

    // Ignore files that are not LAS files.
    String fileName = file.getPath();
    System.out.println("WellLog.load: fileName="+fileName);
    if (!fileName.toLowerCase().endsWith(".las"))
      return null;

    // String data for each part of the LAS file that we need.
    Map<String,String> parts = getParts(file);
    String vdata = parts.get("VERSION");
    String wdata = parts.get("WELL");
    String cdata = parts.get("CURVE");
    String adata = parts.get("A");
    String value;
    Scanner scanner;

    // Ignore files with data wrapped into multiple lines.
    value = getValue(vdata,"WRAP.");
    if (value!=null && value.equals("YES"))
      return null;
    
    // Get well name.
    value = getValue(wdata,"WELL     .");
    if (value==null)
      return null;
    String name = value;

    // Get the well location.
    value = getValue(wdata,"LOC      .");
    if (value==null)
      return null;
    //X = 607903.0000 Y = 6077213.0000
    float xe = Float.parseFloat(value.substring(4,15));
    float yn = Float.parseFloat(value.substring(20,32));
    trace("xe="+xe+" yn="+yn);

    // Get the kelly bushing elevation in m.
    value = getValue(wdata,"EKB      .M");
    float ze = Float.parseFloat(value); // TODO: should this be zero?

    // Map well log curve names to column indices. Currently 
    // we process only depth in m, velocity, density, gamma, and porosity.
    // Measured depth is always first (index zero), as in this example:
    // 012345678901234567890123456789012345678901234567890123
    //  DEPTH    .M                  : 1     Name = Domain Curve ...
    //  CALI     .in                 : 2     Name = Caliper_1 ...
    //  RHOB     .g/cc               : 3     Name = Density_1 ...
    //  GR       .API                : 4     Name = Gamma Ray_math ...
    //  DT       .us/ft              : 5     Name = P-wave_1 ...
    //  DT       .us/ft              : 6     Name = P-wave_corr ...
    //  UNKNOWN  .fraction           : 7     Name = Porosity_1 ...
    scanner = new Scanner(cdata);
    int kz = -1, kv = -1, kd = -1, kg = -1, kp = -1;
    int index = 0;
    while (kz<0 && scanner.hasNextLine()) {
      String line = scanner.nextLine();
      if (line.startsWith(" DEPTH    .M"))
        kz = index;
    }
    while (scanner.hasNextLine()) {
      String line = scanner.nextLine();
      ++index;
      if (line.contains("RHO")) {
        kd = index;
      } else if (line.contains("GR")) {
        kg = index;
      } else if (line.contains("P-wave_1")) {
        kv = index;
      } else if (line.contains("UNKNOWN")) {
        kp = index;
      }
    }
    if (kz<0 || (kv<0 && kd<0 && kg<0 && kp<0))
      return null;

    // Read all well curve data as floats.
    FloatList zl = new FloatList();
    FloatList vl = new FloatList();
    FloatList dl = new FloatList();
    FloatList gl = new FloatList();
    FloatList pl = new FloatList();
    scanner = new Scanner(adata);
    int n = 0;
    while (scanner.hasNextLine()) {
      String line = scanner.nextLine().trim();
      if (line.startsWith("~"))
        break;
      String[] values = line.split("\\s+");
      if (values.length<1 || values[0].length()==0)
        break;
      zl.add(Float.parseFloat(values[kz]));
      if (kv>0) vl.add(Float.parseFloat(values[kv]));
      if (kd>0) dl.add(Float.parseFloat(values[kd]));
      if (kg>0) gl.add(Float.parseFloat(values[kg]));
      if (kp>0) pl.add(Float.parseFloat(values[kp]));
      ++n;
    }

    // Construct the log.
    WellLog log = new WellLog();
    log.name = name;
    log.xe = xe;
    log.yn = yn;
    log.ze = ze;
    log.n = n;
    log.z = zl.trim();
    log.v = vl.trim();
    log.d = dl.trim();
    log.g = gl.trim();
    log.p = pl.trim();

    // Convert velocity and density to SI units.
    siVelocity(log.v);
    //siDensity(log.d); // leave density in g/cc!

    // Clean up logs (or set to null, if too dirty).
    log.v = cleanVelocity(log.v);
    log.d = cleanDensity(log.d);
    log.g = cleanGamma(log.g);
    log.p = cleanPorosity(log.p);

    // If no valid curves, no log.
    if (log.v==null && log.d==null && log.g==null && log.p==null)
      return null;

    // Set resampled (x1,x2,x3) coordinates.
    log.setCoordinates();

    return log;
  }

  /**
   * Gets the well log curve with the specified name.
   * @param curve the curve name; e.g., "velocity".
   * @return array of curve values; null, if none.
   */
  public float[] getCurve(String curve) {
    if (curve.startsWith("v")) return v;
    if (curve.startsWith("d")) return d;
    if (curve.startsWith("g")) return g;
    if (curve.startsWith("p")) return p;
    return null;
  }

  /**
   * Applies a despking filter to all curves for this log.
   * Any null values are ignored and remain null during despiking.
   * @param nmed number of median-of-three filter passes.
   */
  public void despike(int nmed) {
    for (int imed=0; imed<nmed; ++imed) {
      despike(v);
      despike(d);
      despike(g);
      despike(p);
    }
  }
  private void despike(float[] f) {
    if (f==null) 
      return;
    int n = f.length;
    for (int i=1; i<n-1; ++i) {
      float fa = f[i-1];
      float fb = f[i  ];
      float fc = f[i+1];
      if (fa!=NULL_VALUE && fb!=NULL_VALUE && fc!=NULL_VALUE)
        f[i] = med3(fa,fb,fc);
    }
  }
  private float med3(float a, float b, float c) {
    return a<b ? 
           (b<c ? b : (a<c ? c : a)) : 
           (b>c ? b : (a>c ? c : a));
  }

  /**
   * Applies a Gaussian smoothing filter to all curves for this log.
   * Any null values are ignored and remain null during smoothing.
   * @param sigma the half-width of the Gaussian.
   */
  public void smooth(double sigma) {
    int lpad = 1+(int)(3.0*sigma);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    smooth(rgf,lpad,v);
    smooth(rgf,lpad,d);
    smooth(rgf,lpad,g);
    smooth(rgf,lpad,p);
  }
  private void smooth(
    RecursiveGaussianFilter rgf, int lpad, float[] f)
  {
    if (f==null) 
      return;

    // While more non-null samples may exist, ...
    int n = f.length;
    int j = 0;
    while (j<n) {
      
      // Find index j of next non-null sample.
      while (j<n && f[j]==NULL_VALUE)
        ++j;

      // If at least one more non-null sample, ...
      if (j<n) {

        // Find index k of next null sample.
        int k = j;
        while (k<n && f[k]!=NULL_VALUE)
          ++k;

        // Extract and pad sequence of non-null samples.
        float[] fpad = pad(lpad,j,k,f);

        // Apply the smoothing filter.
        rgf.apply0(fpad,fpad);

        // Replace input samples with smoothed samples.
        copy(k-j,lpad,fpad,j,f);

        // Index at which to begin next search.
        j = k;
      }
    }
  }
  private float[] pad(int lpad, int j, int k, float[] f) {
    int n = k-j;
    int npad = lpad+n+lpad;
    float[] fpad = new float[npad];
    for (int i=0; i<lpad; ++i)
      fpad[i] = f[j];
    copy(n,j,f,lpad,fpad);
    for (int i=lpad+n,l=k-1; i<npad; ++i)
      fpad[i] = f[l];
    return fpad;
  }

  /**
   * Gets non-null samples f(x1,x2,x3) for the specified well log curve.
   * @param curve the curve name; e.g., "velocity".
   * @return array of arrays {f,x1,x2,x3} of non-null samples; null if none.
   */
  public float[][] getSamples(String curve) {
    return getSamples(curve,null,null,null);
  }

  /**
   * Gets non-null samples f(x1,x2,x3) for the specified well log curve.
   * Includes only samples within specified sampling bounds.
   * @param curve the curve name; e.g., "velocity".
   * @param s1 sampling in 1st dimension; null, if unbounded.
   * @param s2 sampling in 2nd dimension; null, if unbounded.
   * @param s3 sampling in 3rd dimension; null, if unbounded.
   * @return array of arrays {f,x1,x2,x3} of non-null samples; null if none.
   */
  public float[][] getSamples(
    String curve, Sampling s1, Sampling s2, Sampling s3) 
  {
    double x1min,x1max,x2min,x2max,x3min,x3max;
    if (s1!=null) {
      x1min = s1.getFirst()-0.5*s1.getDelta();
      x1max = s1.getLast() +0.5*s1.getDelta();
    } else {
      x1min = -Float.MAX_VALUE; 
      x1max =  Float.MAX_VALUE;
    }
    if (s2!=null) {
      x2min = s2.getFirst()-0.5*s2.getDelta();
      x2max = s2.getLast() +0.5*s2.getDelta();
    } else {
      x2min = -Float.MAX_VALUE; 
      x2max =  Float.MAX_VALUE;
    }
    if (s3!=null) {
      x3min = s3.getFirst()-0.5*s3.getDelta();
      x3max = s3.getLast() +0.5*s3.getDelta();
    } else {
      x3min = -Float.MAX_VALUE; 
      x3max =  Float.MAX_VALUE;
    }
    float[] f = getCurve(curve);
    if (f==null)
      return null;
    FloatList fl = new FloatList();
    FloatList x1l = new FloatList();
    FloatList x2l = new FloatList();
    FloatList x3l = new FloatList();
    for (int i=0; i<n; ++i) {
      float fi = f[i];
      if (fi!=NULL_VALUE) {
        float x1i = x1[i];
        float x2i = x2[i];
        float x3i = x3[i];
        if (x1min<=x1i && x1i<=x1max && 
            x2min<=x2i && x2i<=x2max && 
            x3min<=x3i && x3i<=x3max) {
          fl.add(fi);
          x1l.add(x1i);
          x2l.add(x2i);
          x3l.add(x3i);
        }
      }
    }
    float[] fs = fl.trim();
    float[] x1s = x1l.trim();
    float[] x2s = x2l.trim();
    float[] x3s = x3l.trim();
    if (fs==null)
      return null;
    return new float[][]{fs,x1s,x2s,x3s};
  }

  /**
   * Returns an array of derivatives dx1/dz.
   * @return array of derivatives.
   */
  public float[] dx1dz() {
    float s = 1.0f/0.0003048f;
    float[] d = new float[n];
    for (int i=1; i<n; ++i)
      d[i] = s*(x1[i]-x1[i-1])/(z[i]-z[i-1]);
    d[0] = d[1];
    return d;
  }

  /**
   * Sets resampled (x1,x2,x3) coordinates for well bore.
   */
  private void setCoordinates() {
    x1 = new float[n];
    x2 = new float[n];
    x3 = new float[n];
    for (int i=0; i<n; ++i) {
      // Depths in depth-time pairs in other files are measured depths, 
      // which include the 30 m from the sea surface to the kelly bushing. 
      // In those files, 30 m <=> 0 s, so our first depth must be 30 m.
      //Coordinates.Odt odt = new Coordinates.Odt(xe,yn,ze-z[i]);
      Coordinates.Odt odt = new Coordinates.Odt(xe,yn,-z[i]);
      Coordinates.Csm csm = new Coordinates.Csm(odt);
      if (i==0) {
        trace("odt: ze="+odt.ze+" xe="+odt.xe+" yn="+odt.yn);
        trace("csm: x1="+csm.x1+" x2="+csm.x2+" x3="+csm.x3);
      }
      x1[i] = (float)csm.x1;
      x2[i] = (float)csm.x2;
      x3[i] = (float)csm.x3;
    }
    convertFromDepthToTime();
  }

  /**
   * Converts x1 coordinates from depth to time.
   */
  private void convertFromDepthToTime() {
    trace("converting depth to time for "+name);
    float[][] zt = getDepthTime(name);
    float[] z = zt[0];
    float[] t = zt[1];
    SimplePlot.asPoints(z,t);
    CubicInterpolator.Method method = CubicInterpolator.Method.LINEAR;
    CubicInterpolator ci = new CubicInterpolator(method,z,t);
    for (int i=0; i<n; ++i)
      x1[i] = ci.interpolate(x1[i]);
  }

  // Depth-time pairs from txt files in .../F3_Demo/WellInfo/DT_model/.
  // In those files depths are in meters and times are in seconds.
  private static String[] _wellNames = {"F02-01","F03-02","F03-04","F06-01"};
  private static double[][][] _wellDepthTimes = {
    { {30,0}, {553.6,0.544}, {612.9,0.607}, {683.31,0.675}, {716.65,0.712},
      {748.49,0.748}, {795.18,0.794}, {927.28,0.932}, {1025.42,1.031},
      {1048.84,1.051}, {1057.89,1.059}, {1068.74,1.068}, {1094.12,1.086},
      {1134.73,1.117}, {1174.62,1.146}, {1197.08,1.165}, {1252.26,1.213},
      {1285.09,1.242}, {1695,1.67}, {1872,1.861}, {2636,2.682},
      {2735,2.788}, {2997.5,3.07}, {3028,3.103}, {3150,3.234} },
    { {30,0}, {486.14,0.485}, {522.2,0.522}, {564.4,0.565}, {589.92,0.592}, 
      {610.9,0.613}, {676.43,0.675}, {715.02,0.71}, {835.93,0.82}, 
      {890.32,0.868}, {957.29,0.927}, {1036.94,0.995}, {1073.3,1.024}, 
      {1111.79,1.056}, {1249.72,1.178}, {1647,1.569}, {1885,1.698},
      {1932,1.732} },
    { {30,0}, {479.74,0.473}, {510.52,0.504}, {547.75,0.543}, {569.66,0.565},
      {611.39,0.61}, {664.95,0.662}, {697.11,0.692}, {856.08,0.84},
      {907.88,0.886}, {946.5,0.921}, {1024.39,0.991}, {1071.12,1.032},
      {1111.13,1.067}, {1200.93,1.15}, {1246,1.195}, {1359,1.304},
      {1785,1.74}, {1858,1.773}, {1970,1.873}, {2277,2.147}, {2295,2.163}, 
      {2317,2.182}, {2330,2.194}, {2362,2.222}, {2375,2.234}, {2736,2.556}, 
      {2806,2.618}, {2885,2.689}, {2923,2.722}, {2980,2.773}, {3092,2.873}, 
      {3111,2.89}, {3130,2.907} },
    { {30,0}, {589.14,0.593}, {670.54,0.677}, {697.12,0.705}, {725.25,0.735},
      {771.24,0.779}, {866.29,0.886}, {995.97,1.028}, {1015.24,1.042},
      {1030.08,1.053}, {1049.63,1.066}, {1077.92,1.087}, {1114.22,1.113},
      {1151.2,1.138}, {1170.07,1.152}, {1225.53,1.193}, {1261.34,1.221},
      // these last depth-time pairs were copied from nearby F02-01 above
      {1285.09,1.242}, {1695,1.67}, {1872,1.861}, {2636,2.682},
      {2735,2.788}, {2997.5,3.07}, {3028,3.103}, {3150,3.234} }
  };
  private static float[][] getDepthTime(String name) {
    for (int i=0; i<_wellNames.length; ++i) {
      if (_wellNames[i].equals(name)) {
        double[][] zt = _wellDepthTimes[i];
        int n = zt.length;
        float[] z = new float[n];
        float[] t = new float[n];
        for (int j=0; j<n; ++j) {
          z[j] = (float)(zt[j][0]*0.001);
          t[j] = (float)zt[j][1];
        }
        return new float[][]{z,t};
      }
    }
    return null;
  }

  /**
   * Each part of an LAS file begins with a line like "~VERSION INFORMATION".
   * In the returned map, the key is the part name = the first word in that
   * line (e.g., "VERSION"), and the value is all data for that part stored 
   * as a single string with embedded newlines. The single string of data is 
   * easy to scan.
   */
  private static Map<String,String> getParts(File file) {
    Map<String,String> parts = new HashMap<String,String>();
    try {
      BufferedReader br = new BufferedReader(new FileReader(file));
      String name = null;
      StringBuilder sb = null;
      String line = null;
      while ((line=br.readLine())!=null) {
        if (line.startsWith("~")) {
          if (name!=null)
            parts.put(name,sb.toString());
          int end = line.indexOf(" ");
          if (end<0) end = line.length();
          name = line.substring(1,end);
          name = name.toUpperCase();
          sb = new StringBuilder();
        } else {
          sb.append(line);
          sb.append('\n');
        }
      }
      if (name!=null)
        parts.put(name,sb.toString());
      br.close();
      return parts;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Searches a data string for a specified key and, if found, returns
   * the corresponding string value. The key-values pairs in LAS files
   * appear on separate lines like this one:
   * " STRT     .M        30.0000             : START DEPTH"
   * In this example, key is "STRT     .M" and value is "30.0000".
   * @return the string value; null if the key was not found.
   */
  private static String getValue(String data, String key) {
    Matcher m;
    if (key.equals("LOC      .")) {
      Pattern p = Pattern.compile(" "+key+"\\s+?(.+):");
      m = p.matcher(data);
    } else {
      Pattern p = Pattern.compile(" "+key+"\\s+?(\\S+)\\s*:");
      m = p.matcher(data);
    }
    return m.find()?m.group(1):null;
  }

  private static void siVelocity(float[] v) {
    if (v==null) return;
    int n = v.length;
    float sv = 304.8f; // from ft/us to km/s
    for (int i=0; i<n; ++i)
      if (v[i]!=NULL_VALUE) v[i] = sv/v[i]; // divide, since log is us/ft!
  }

  private static void siDensity(float[] d) {
    if (d==null) return;
    int n = d.length;
    float sd = 1000.0f; // from g/cc to kg/m3
    for (int i=0; i<n; ++i)
      if (d[i]!=NULL_VALUE) d[i] = sd*d[i];
  }

  private static float[] cleanVelocity(float[] v) {
    if (v==null) return null;
    int n = v.length;
    for (int i=0; i<n; ++i)
      if (v[i]!=NULL_VALUE && !validVelocity(v[i])) return null;
    return v;
  }
  private static boolean validVelocity(float v) {
    return 0.20f<v && v<20.0f;
  }

  private static float[] cleanDensity(float[] d) {
    if (d==null) return null;
    int n = d.length;
    for (int i=0; i<n; ++i)
      if (d[i]!=NULL_VALUE && !validDensity(d[i])) return null;
    return d;
  }
  private static boolean validDensity(float d) {
    return 0.5f<d && d<10.0f;
  }

  private static float[] cleanGamma(float[] g) {
    if (g==null) return null;
    int n = g.length;
    for (int i=0; i<n; ++i) {
      if (g[i]!=NULL_VALUE && !validGamma(g[i])) 
        return null;
    }
    return g;
  }
  private static boolean validGamma(float g) {
    if (g<0.0f) g = 0.0f;
    return 0.0f<=g && g<300.0f;
  }

  private static float[] cleanPorosity(float[] p) {
    if (p==null) return null;
    int n = p.length;
    for (int i=0; i<n; ++i)
      if (p[i]!=NULL_VALUE && !validPorosity(p[i])) return null;
    return p;
  }
  private static boolean validPorosity(float p) {
    return 0.0f<p && p<0.8f;
  }
  private static void trace(String s) {
    System.out.println(s);
  }
}
