package tp;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import static java.lang.Math.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.util.*;

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
  public float[] z,v,d,g,p; // depth, velocity, density, gamma, porosity
  public float[] x1,x2,x3; // resampled coordinates of well-bore (km)

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
   * Returns an array of derivatives (1/0.0003048)*dx1/dz.
   * The scale factor 0.0003048 (km/ft) compensates for the difference 
   * in units for dz1 (km) and dz (ft). This derivative should never 
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
        Coordinates.Resampled r = new Coordinates.Resampled(m);
        x1[i] = (float)r.x1;
        x2[i] = (float)r.x2;
        x3[i] = (float)r.x3;
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
}
