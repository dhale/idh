package tp;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import static java.lang.Math.*;

import edu.mines.jtk.util.*;

/**
 * A well log from Teapot Dome.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.21
 */
public class WellLog {

  // This value should be much lower than the elevation of any Kelling bushing.
  // Some LAS files have missing or invalid elevations.
  public static final float MIN_ELEVATION = 100.0f;

  // This value represents missing values in all well curve data.
  public static final float NULL_VALUE = -999.2500f;

  public String id; // unique 10-character API well number
  public double xe,yn,ze; // easting, northing, elev of Kelling bushing (ft)
  public int n; // number of samples
  public float[] z,v,d,g,p; // depth, velocity, density, gamma, porosity
  public float[] x1,x2,x3; // resampled coordinates of well-bore (km)

  /**
   * Loads a well log from the specified LAS file.
   * @param file the LAS file
   * @return the well log; null, if the file format is not valid.
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
    if (value==null || value.length()<10 || !value.startsWith("49025"))
      return null;
    for (int i=10; i<value.length(); ++i)
      if (value.charAt(i)!='0')
        return null;
    String id = value.substring(0,10);

    // Get the elevation of the Kelling bushing in ft.
    // Assume small elevations are not valid.
    value = getValue(pdata,"EKB .F");
    if (value==null)
      return null;
    float ze = Float.parseFloat(value);
    if (ze<MIN_ELEVATION)
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

    // Read the well curve data as floats.
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
    log.n = n;
    log.z = zl.trim();
    log.v = vl.trim();
    log.d = dl.trim();
    log.g = gl.trim();
    log.p = pl.trim();
    log.setCoordinates(0.0,0.0);
    System.out.println("id="+id+" ze="+ze+" n="+n);
    return log;
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

  // Set (x1,x2,x3) coordinates using directional survey data.
  private static class Interpolator {
    public Interpolator(
      double xe, double yn, double ze, 
      double[] z, double[] t, double[] p) 
    {
      int n = z.length;

      // Sort directional data by increasing measured depth.
      int[] j = Array.rampint(0,1,n);
      Array.quickIndexSort(z,j);
      z = sort(z,j);
      t = sort(t,j);
      p = sort(p,j);

      // If smallest depth not zero, then insert data for depth zero.
      // Assume that both inclination and azimuth are zero at depth zero.
      if (z[0]!=0.0) {
        double[] ztemp = new double[n+1]; 
        double[] ttemp = new double[n+1]; 
        double[] ptemp = new double[n+1]; 
        System.arraycopy(z,0,ztemp,1,n);
        System.arraycopy(t,0,ttemp,1,n);
        System.arraycopy(p,0,ptemp,1,n);
        z = ztemp;
        t = ttemp;
        p = ptemp;
        ++n;
      }
      _n = n;

      // Parameters at measured depths required for interpolation.
      _a = new double[n];
      _zm = new double[n];
      _xp = new double[n];
      _yp = new double[n];
      _zp = new double[n];
      _xt = new double[n];
      _yt = new double[n];
      _zt = new double[n];
      _xp[0] = xe;
      _yp[0] = yn;
      _zp[0] = ze;
      for (int i=0,im1=0; i<n; ++i,im1=i-1) {
        double t1 = toRadians(t[im1]);
        double t2 = toRadians(t[i]);
        double p1 = toRadians(p[im1]);
        double p2 = toRadians(p[i]);
        double ct1 = cos(t1);
        double ct2 = cos(t2);
        double cp1 = cos(p1);
        double cp2 = cos(p2);
        double st1 = sin(t1);
        double st2 = sin(t2);
        double sp1 = sin(p1);
        double sp2 = sin(p2);
        double sdt = sin(0.5*(t2-t1));
        double sdp = sin(0.5*(p2-p1));
        double a = 2.0*asin(sqrt(sdt*sdt+st1*st2*sdp*sdp));
        double scale = 0.5*shapeFactor(a)*(z[i]-z[im1]);
        _a[i] = a;
        _zm[i] = z[i];
        _xt[i] = st2*sp2;
        _yt[i] = st2*cp2;
        _zt[i] = ct2;
        if (i==0) {
          _xp[i] = xe;
          _yp[i] = yn;
          _zp[i] = ze;
        } else {
          _xp[i] = _xp[im1]+scale*(_xt[i-1]+_xt[i]);
          _yp[i] = _yp[im1]+scale*(_yt[i-1]+_yt[i]);
          _zp[i] = _zp[im1]-scale*(_zt[i-1]+_zt[i]);
        }
      }
    }
    private int _n; // number of points with directional data
    private double[] _zm; // measured depths for all directional data
    private double[] _xp,_yp,_zp; // coordinates of points along wellbore
    private double[] _xt,_yt,_zt; // unit tangent vectors for all points
    private double[] _a; // angles of arc
                         // a[1] is for interval between zm[0] and zm[1]
                         // a[2] is for interval between zm[1] and zm[2]
                         // ...
    private static double shapeFactor(double a) {
      if (a<0.02) {
        double aa = a*a;
        return 1.0+aa/12.0*(1.0+aa/10.0*(1.0+aa/168.0*(1+31.0*aa/18.0)));
      } else {
        return tan(0.5*a)/(0.5*a);
      }
    }
    private static double[] sort(double[] x, int[] j) {
      int n = x.length;
      double[] y = new double[n];
      for (int i=0; i<n; ++i)
        y[i] = x[j[i]];
      return y;
    }
    public void interpolate(double zm, double[] x1x2x3) {
      /*
      int n = _n;

      // Find interval of measured depth.
      int i = Array.binarySearch(_zm,zm);
      if (i<0) i = -(i+1);

      // Coordinates of point corresponding to specified measured depth.
      double xp,yp,zp;

      // If extrapolating, ...
      if (i==0) {
        xp = _xp[0];
        yp = _yp[0];
        zp = _zp[0];
      } else if (i==n) {
        double dz = zm-_zm[n-1];
        xp = _xp[n-1]+dz*_xt[n-1];
        yp = _yp[n-1]+dz*_yt[n-1];
        zp = _zp[n-1]+dz*_zt[n-1];
      } 

      // Otherwise, if interpolating, ...
      else {
        int im1 = i-1;
        double zm1 = _zm[im1];
        double zm2 = _zm[i];
        double a12 = _a[im1];
        double a = a12*(zm-zm1)/(zm2-zm1); 
        if (a<0.02) {
        } else {
          double s1 = sin(a12-a);
          double s2 = sin(a);
          double s12 = sin(a12);
          xp1 = _xp[im1];
          yp1 = _yp[im1];
          zp1 = _zp[im1];
          xt1 = _xt[im1];
          yt1 = _yt[im1];
          zt1 = _zt[im1];
          xt2 = _xt[i];
          yt2 = _yt[i];
          zt2 = _zt[i];
          xt = (s1*xt1+s2*xt2)/s12;
          yt = (s1*yt1+s2*yt2)/s12;
          zt = (s1*zt1+s2*zt2)/s12;
          scale = 0.5*shapeFactor(a)*(zm-zm1);
          xp = xp1+scale*(xt1+xt);
          yp = yp1+scale*(yt1+yt);
          zp = zp1-scale*(zt1+zt);
        }
      }
      */
    }
  }
  public void setCoordinates(
    double xe, double yn, double[] z, double[] t, double[] p)
  {
    this.xe = xe;
    this.yn = yn;
    if (z[0]!=0.0) {
      int n = z.length;
      double[] zt = new double[n+1]; 
      double[] tt = new double[n+1]; 
      double[] pt = new double[n+1]; 
      System.arraycopy(z,0,zt,1,n);
      System.arraycopy(t,0,tt,1,n);
      System.arraycopy(p,0,pt,1,n);
    }
  }
  public void setCoordinates(double xe, double yn) {
    this.xe = xe;
    this.yn = yn;
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
  private static class FloatList {
    public int n;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
}
