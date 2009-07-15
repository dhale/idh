package tp;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A directional survey from Teapot Dome.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.26
 */
public class DirectionalSurvey {
  public long id; // unique 12-digit API well number
  public int n; // number of survey samples
  public double[] z,t,p; // arrays of measured depth, theta, and phi

  /**
   * Collection of directional survey data.
   */
  public static class Data { 

    /**
     * Constructs directional survey data from the specified file.
     * The file is in Excel CSV (comma-delimited) format.
     * @param fileName file containing directional survey data.
     */
    public Data(String fileName) {
      try {
        FileInputStream fis = new FileInputStream(fileName);
        Scanner s = new Scanner(fis);
        long id = -1;
        DoubleList zl = new DoubleList();
        DoubleList tl = new DoubleList();
        DoubleList pl = new DoubleList();
        while (s.hasNextLine()) {
          String line = s.nextLine();
          String[] fields = line.split(",");
          if (fields.length<4)
            continue;
          long idline = WellLog.idFromString(fields[0]);
          if (idline<0)
            continue;
          if (id<0) {
            id = idline;
          } else if (id!=idline || !s.hasNextLine()) {
            DirectionalSurvey ds = new DirectionalSurvey();
            ds.id = id;
            ds.z = zl.trim();
            ds.t = tl.trim();
            ds.p = pl.trim();
            if (id==490252304800L) { // Swap theta,phi for this survey!
              double[] temp = ds.t;
              ds.t = ds.p;
              ds.p = temp;
            }
            ds.n = ds.z.length;
            _data.put(id,ds);
            id = idline;
            zl = new DoubleList();
            tl = new DoubleList();
            pl = new DoubleList();
          }
          zl.add(Double.parseDouble(fields[1]));
          tl.add(Double.parseDouble(fields[2]));
          pl.add(Double.parseDouble(fields[3]));
        }
        s.close();
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

    /**
     * Adds the specified directional survey.
     * @param ds the directional survey.
     */
    public void add(DirectionalSurvey ds) {
      _data.put(ds.id,ds);
    }

    /**
     * Returns the number of directional surveys.
     * @return the number of directional surveys.
     */
    public int size() {
      return _data.size();
    }

    /**
     * Gets the directional survey for the specified well id.
     * @param id the well id.
     * @return the directional survey; null, if none.
     */
    public DirectionalSurvey get(long id) {
      return _data.get(id);
    }

    /**
     * Gets all directional surveys.
     * @return list of directional surveys.
     */
    public List<DirectionalSurvey> getAll() {
      return new ArrayList<DirectionalSurvey>(_data.values());
    }

    /**
     * Prints summary information for these directional survey data.
     */
    public void printInfo() {
      for (DirectionalSurvey ds:getAll())
        System.out.println("id="+ds.id+" n="+ds.n +
                           " first z="+ds.z[0]+" last z="+ds.z[ds.n-1]);
      System.out.println("number of surveys = "+_data.size());
    }
    private Map<Long,DirectionalSurvey> _data = 
      new HashMap<Long,DirectionalSurvey>();
  }

  /**
   * Sets resampled coordinates for the specified well log.
   * Uses this directional survey and the log's map coordinates (xe,yn,ze)
   * to compute the resampled coordinates (x1,x2,x3).
   * @param log the well log.
   */
  public void setCoordinates(WellLog log) {
    assert id==log.id:"id of well log and directional survey match";
    Locater locater = new Locater(log.xe,log.yn,log.ze,z,t,p);
    int n = log.z.length;
    log.x1 = new float[n];
    log.x2 = new float[n];
    log.x3 = new float[n];
    for (int i=0; i<n; ++i) {
      Coordinates.Csm csm = locater.locate(log.z[i]);
      log.x1[i] = (float)csm.x1;
      log.x2[i] = (float)csm.x2;
      log.x3[i] = (float)csm.x3;
    }
  }

  // Computes wellbore coordinates (x1,x2,x3) from directional survey data.
  // The directional survey data is (z,t,p), where z is measured depth,
  // t (theta) is the inclination angle, and p (phi) is the azimuthal angle.
  // Reference:
  // Sawaryn & Thorogood, 2005, A compendium of directional calculations 
  // based on the minimum curvature method: SPE Drilling & Completion.
  private static class Locater {
    public Locater(
      double xe, double yn, double ze, 
      double[] z, double[] t, double[] p) 
    {
      int n = z.length;

      // Sort directional data by increasing measured depth. This sort 
      // facilitates binary search for the measured depth interval that
      // contains a measured depth at which we compute well bore coordinates.
      int[] j = rampint(0,1,n);
      quickIndexSort(z,j);
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
    public Coordinates.Csm locate(double zmi) {
      int n = _n;

      // Find interval of measured depth.
      int i = binarySearch(_zm,zmi);
      if (i<0) i = -(i+1);

      // Coordinates of point corresponding to specified measured depth.
      double xp,yp,zp;

      // If beyond the ends of the directional survey, extrapolate.
      if (i==0) {
        xp = _xp[0];
        yp = _yp[0];
        zp = _zp[0];
      } else if (i==n) {
        double scale = zmi-_zm[n-1];
        xp = _xp[n-1]+scale*_xt[n-1];
        yp = _yp[n-1]+scale*_yt[n-1];
        zp = _zp[n-1]-scale*_zt[n-1];
      } 

      // Otherwise, interpolate using the minimum curvature method.
      else {
        int im1 = i-1;
        double zm1 = _zm[im1];
        double zm2 = _zm[i];
        double a12 = _a[i];
        double c = (zmi-zm1)/(zm2-zm1);
        double a1i = c*a12;
        double w1 = sinRatio(1.0-c,a12);
        double w2 = sinRatio(c,a12);
        double xp1 = _xp[im1];
        double yp1 = _yp[im1];
        double zp1 = _zp[im1];
        double xt1 = _xt[im1];
        double yt1 = _yt[im1];
        double zt1 = _zt[im1];
        double xt2 = _xt[i];
        double yt2 = _yt[i];
        double zt2 = _zt[i];
        double xti = (w1*xt1+w2*xt2);
        double yti = (w1*yt1+w2*yt2);
        double zti = (w1*zt1+w2*zt2);
        double scale = 0.5*shapeFactor(a1i)*(zmi-zm1);
        xp = xp1+scale*(xt1+xti);
        yp = yp1+scale*(yt1+yti);
        zp = zp1-scale*(zt1+zti);
      }

      // From map coordinates (xp,yp,zp) to resampled coordinates (x1,x2,x3).
      return new Coordinates.Csm(new Coordinates.Map(xp,yp,zp));
    }
    private int _n; // number of points with directional data
    private double[] _zm; // measured depths for all directional data
    private double[] _xp,_yp,_zp; // coordinates of points along wellbore
    private double[] _xt,_yt,_zt; // unit tangent vectors for all points
    private double[] _a; // angles of arc
                         // a[1] is for interval between zm[0] and zm[1]
                         // a[2] is for interval between zm[1] and zm[2]
                         // ...

    // Computes shape factor tan(a/2)/(a/2), for 0 <= a <= PI/2.
    private static double shapeFactor(double a) {
      if (a<0.02) {
        double aa = a*a;
        return 1.0+aa/12.0*(1.0+aa/10.0*(1.0+aa/168.0*(1+31.0*aa/18.0)));
      } else {
        return tan(0.5*a)/(0.5*a);
      }
    }

    // Returns sin(ca)/sin(a), for 0 <= c <= 1 and 0 <= a <= PI/2.
    private static double sinRatio(double c, double a) {
      if (a>=0.02) {
        return sin(c*a)/sin(a);
      } else {
        double cc = c*c;
        double aa = a*a;
        return c + 
          aa*(c*(1.0/6.0-cc/6.0) +
          aa*(c*(7.0/360.0-cc*(1.0/36.0-cc/120.0)) +
          aa*(c*(31.0/15120.0-cc*(7.0/2160.0-cc*(1.0/720-cc/5040.0))) +
          aa*c*(127.0/604800.0-cc*(31.0/90720.0 -
            cc*(7.0/43200.0-cc*(1.0/30240.0-cc/362880.0)))))));
      }
    }

    private static double[] sort(double[] x, int[] j) {
      int n = x.length;
      double[] y = new double[n];
      for (int i=0; i<n; ++i)
        y[i] = x[j[i]];
      return y;
    }
  }
}
