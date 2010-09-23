package tp;

import java.io.*;
import static java.lang.Math.abs;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A well log from Teapot Dome.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.21
 */
public class WellLog {

  public long id; // unique 12-digit API well number
  public double xe,yn,ze; // easting, northing, elevation
  public int n; // number of samples
  public float[] z; // measured depth (along wellbore)
  public float[] v,d,g,p; // velocity, density, gamma, porosity; may be null
  public float[] x1,x2,x3; // CSM resampled coordinates of well-bore (km)

  // This value should be much lower than any reasonable elevation.
  // For Teapot Dome, elevations are approximately 5000 ft.
  // We ignore LAS files that have missing or invalid elevations.
  public static final float MIN_ELEVATION = 100.0f;

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
     * @param whFileName file containing well header data.
     * @param dsFileName file containing directional survey data.
     */
    public Data(String wlDirName, String whFileName, String dsFileName) {
      WellHeader.Data whdata = null;
      if (whFileName!=null)
        whdata = new WellHeader.Data(whFileName);
      DirectionalSurvey.Data dsdata = null;
      if (dsFileName!=null)
        dsdata = new DirectionalSurvey.Data(dsFileName);
      File wlDir = new File(wlDirName);
      ArrayList<Long> wldup = new ArrayList<Long>();
      load(wlDir,whdata,dsdata,this,wldup);
      System.out.println("WellLog.Data: after initial loading, ...");
      System.out.println("  number of logs = "+size());
      System.out.println("  number of duplicate ids = "+wldup.size());
      for (long id:wldup)
        _data.remove(id);
      System.out.println("WellLog.Data: after removing duplicates, ...");
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
          aos.writeLong(log.id);
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
          log.id = ais.readLong();
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
      _data.put(log.id,log);
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
      Map<Long,WellLog> data = new HashMap<Long,WellLog>();
      for (long id:_data.keySet()) {
        WellLog log = _data.get(id);
        float[] x2 = log.x2;
        float[] x3 = log.x3;
        int n = log.n;
        boolean ok = true;
        for (int i=0; i<n && ok; ++i)
          ok = f2<=x2[i] && x2[i]<=l2 && f3<=x3[i] && x3[i]<=l3;
        if (ok)
          data.put(id,log);
      }
      _data = data;
    }

    /**
     * Gets the well log for the specified well id.
     * @param id the well id.
     * @return the well log; null, if none.
     */
    public WellLog get(long id) {
      return _data.get(id);
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

    private Map<Long,WellLog> _data = new HashMap<Long,WellLog>();
    private static void load(
      File wlDir, WellHeader.Data whdata, DirectionalSurvey.Data dsdata,
      WellLog.Data wldata, ArrayList<Long> wldup)
    {
      Check.argument(wlDir.isDirectory(),wlDir+" is a directory");
      for (File file:wlDir.listFiles()) {
        if (file.isDirectory()) {
          load(file,whdata,dsdata,wldata,wldup);
        } else if (file.isFile()) {
          WellLog log = WellLog.load(file,whdata,dsdata);
          if (log!=null) {
            if (wldata.get(log.id)!=null) {
              wldup.add(log.id);
            } else {
              wldata.add(log);
            }
          }
        }
      }
    }
  }

  /**
   * Loads one well log from the specified LAS file, if possible.
   * @param file the LAS file.
   * @param whdata well header data; null, if none.
   * @param dsdata directional survey data; null, if none.
   * @return the well log; null, if the file format is not understood.
   */
  public static WellLog load(
    File file, WellHeader.Data whdata, DirectionalSurvey.Data dsdata) 
  {

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
    String pdata = parts.get("PARAMETER");
    String adata = parts.get("A");
    String value;
    Scanner scanner;

    // Ignore files with data wrapped into multiple lines.
    value = getValue(vdata,"WRAP.");
    if (value!=null && value.equals("YES"))
      return null;
    
    // Get first ten digits of the well API number; ignore file if not valid.
    value = getValue(wdata,"API .");
    if (value==null)
      return null;
    long id = idFromString(value);
    if (id<0)
      return null;

    // Get the log-measured-from elevation in ft.
    value = getValue(pdata,"LMF .");
    if (value==null)
      return null;
    if (value.equals("KB")) {
      value = getValue(pdata,"EKB .F");
    } else if (value.equals("PD")) {
      value = getValue(pdata,"EPD .F");
    } else if (value.equals("DF")) {
      value = getValue(pdata,"EDF .F");
    } else if (value.equals("GL")) {
      value = getValue(pdata,"EGL .F");
    } else {
      return null;
    }
    if (value==null)
      return null;
    float ze = Float.parseFloat(value);
    if (ze<MIN_ELEVATION) // check for wildly erroneous elevations
      return null;

    // Must have well header with map coordinates of well and the well
    // header elevation must match that in LAS file to within 0.5 ft.
    WellHeader wh = (whdata!=null)?whdata.get(id):null;
    if (wh==null)
      return null;
    if (abs(ze-wh.ze)>0.5)
      return null;

    // Map well log curve names (API codes) to column indices. Currently 
    // we process only depth in ft, velocity, density, gamma, and porosity.
    // Measured depth is always first (index zero), as in this example:
    // 012345678901234567890123456789012345678901234567890123
    //  DEPT.F                  :   
    //  RHOB.G/C3   45 350 01 00:   BULK DENSITY
    //  CNC .DEC    42 890 03 00:   NEUTRON POROSITY SANDSTONE
    scanner = new Scanner(cdata);
    int kz = -1, kv = -1, kd = -1, kg = -1, kp = -1;
    int index = 0;
    while (kz<0 && scanner.hasNextLine()) {
      String line = scanner.nextLine();
      if (line.startsWith(" DEPT.F"))
        kz = index;
    }
    while (scanner.hasNextLine()) {
      String line = scanner.nextLine();
      ++index;
      if (line.length()>22) {
        String api1 = line.substring(13,15);
        String api2 = line.substring(16,19);
        String api3 = line.substring(20,22);
        if (api2.equals("310")) {
          kg = index;
        } else if (api2.equals("250") ||
                   api2.equals("520")) {
          kv = index;
        } else if (api2.equals("890")) {
          kp = index;
        } else if (api2.equals("350")) {
          kd = index;
        }
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
    log.id = id;
    log.ze = ze;
    log.xe = wh.xe;
    log.yn = wh.yn;
    //log.ze = wh.ze;
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
    log.setCoordinates(id,dsdata);

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
   * Returns an array of derivatives (1/0.0003048)*dx1/dz.
   * The scale factor 0.0003048 (km/ft) compensates for the difference 
   * in units for dx1 (km) and dz (ft). This derivative should never 
   * exceed one, and for near-vertical wells it will be close to one.
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
   * Returns a 12-digit well log id corresponding to the specified string.
   * Any characters after the first 10 must be '0'. That is, this method
   * does not accept well log id strings for sidetrack wells.
   * @param s the string representing the well log id.
   * @return the well log id; -1, if string is not valid.
   */
  public static long idFromString(String s) {
    long id = -1;
    Matcher m = _idPattern.matcher(s);
    if (m.find()) {
      s = m.group(1);
      if (s.length()==10) {
        s = s+"00";
      } else if (s.length()==11) {
        s = s+"0";
      }
      id =  Long.parseLong(s);
      if (id%100!=0) // ignore sidetrack wells
        id = -1;
    }
    return id;
  }
  private static Pattern _idPattern = Pattern.compile("^(49025\\d{5,7}).*");

  /**
   * Sets resampled (x1,x2,x3) coordinates for well bore.
   * Assumes a vertical well unless a directional survey is available.
   * @param id the well id.
   * @param dsdata directional survey data.
   */
  private void setCoordinates(long id, DirectionalSurvey.Data dsdata) {
    DirectionalSurvey ds = (dsdata!=null)?dsdata.get(id):null;
    if (ds!=null) {
      ds.setCoordinates(this);
    } else {
      x1 = new float[n];
      x2 = new float[n];
      x3 = new float[n];
      for (int i=0; i<n; ++i) {
        Coordinates.Map m = new Coordinates.Map(xe,yn,ze-z[i]);
        Coordinates.Csm c = new Coordinates.Csm(m);
        x1[i] = (float)c.x1;
        x2[i] = (float)c.x2;
        x3[i] = (float)c.x3;
      }
    }
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
   * " STRT.F          200.0000:   START DEPTH"
   * In this example, the key is "STRT.F" and the value is "200.0000".
   * @return the string value; null if the key was not found.
   */
  private static String getValue(String data, String key) {
    Pattern p = Pattern.compile(" "+key+"\\s+?(\\S+):");
    Matcher m = p.matcher(data);
    return (m.find())?m.group(1):null;
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
    for (int i=0; i<n; ++i)
      if (g[i]!=NULL_VALUE && !validGamma(g[i])) return null;
    return g;
  }
  private static boolean validGamma(float g) {
    return 0.0f<g && g<300.0f;
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
}
