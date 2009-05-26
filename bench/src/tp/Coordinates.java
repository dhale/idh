package tp;

import static java.lang.Math.*;

/**
 * Coordinate transforms for Teapot Dome data.
 * Coordinates in many data files are map coordinates (xe,yn,ze), where 
 * xe is easting, yn is northing, and ze is elevation, all measured in ft.
 * In some files depths zd are specified instead of elevations ze.
 * 
 * This class defines two additional coordinate systems, one for seismic
 * survey coordinates (x1,x2,x3) = (depth,inline,crossline) and another
 * that is a resampled (translated and rotated) version of those seismic
 * coordinates. The resampled coordinate system enables many of the dead
 * traces in the seismic survey to be ignored while maintaining regular
 * sampling of the 3D seismic image.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.19
 */
public class Coordinates {

  // Conversion ft <=> km.
  static final double KM_PER_FT = 0.0003048;
  static final double FT_PER_KM = 1.0/KM_PER_FT;

  // Eastings and northings - (xe,yn) coordinates - of seismic rectangle.
  // Coordinates are measured in ft, obtained from teapot_3d_load.doc.
  static final double XLL_FT=788937.0, YLL_FT=938846.0; // (i2,i3) = (0,0)
  static final double XLR_FT=809502.0, YLR_FT=939334.0; // (i2,i3) = (n2-1,0)
  static final double XUR_FT=808604.0, YUR_FT=977163.0; // (i2,i3) = (n2-1,n3-1)
  static final double XUL_FT=788039.0, YUL_FT=976675.0; // (i2,i3) = (0,n3-1)

  // Seismic datum elevation measured in ft. Zero time or depth in seismic 
  // or resampled coordinates corresponds to this elevation in map coordinates.
  static final double DATUM_FT = 5500.0;

  // Rotation angle of seismic rectangle. This is the angle between
  // the seismic x2 (inline) axis and the map xe (easting) axis,
  // measured CCW in the map (xe,yn) coordinate system.
  static final double PHIS = atan2(YLR_FT-YLL_FT,XLR_FT-XLL_FT);
  static final double COSS = cos(PHIS);
  static final double SINS = sin(PHIS);

  // Resampled (x1,x2,x3) coordinates to trim off most dead traces.
  // This is a translation and rotation. The rotation center (in km) 
  // and angle (in radians) were determined visually. Seismic coordinates 
  // (x2s,x3s) and resampled coordinates (x2r,x3r) are measured in km and
  // related by:
  //   x2s = X2CR + x2r*cos(PHIR) + x3r*sin(PHIR)
  //   x3s = X3CR - x2r*sin(PHIR) + x3r*cos(PHIR)
  static final double X2CR = 4.192, X3CR = 0.935; // in km
  static final double PHIR = -0.485364; // in radians; equals -27.8093 degrees
  static final double COSR = cos(PHIR);
  static final double SINR = sin(PHIR);

  // Map coordinates xe (easting), yn (northing), and ze (elevation) in ft.
  public static class Map {
    public double xe,yn,ze;
    public Map(double xe, double yn) {
      this(xe,yn,0.0);
    }
    public Map(double xe, double yn, double ze) {
      this.xe = xe;
      this.yn = yn;
      this.ze = ze;
    }
    public Map(Seismic s) {
      double x1 = s.x1;
      double x2 = s.x2;
      double x3 = s.x3;
      x1 *= FT_PER_KM;
      x2 *= FT_PER_KM;
      x3 *= FT_PER_KM;
      xe = COSS*x2-SINS*x3;
      yn = SINS*x2+COSS*x3;
      ze = -x1;
      xe += XLL_FT;
      yn += YLL_FT;
      ze += DATUM_FT;
    }
    public Map(Resampled r) {
      this(new Seismic(r));
    }
  }

  // Seismic coordinates x1 (depth), x2 (inline) and x3 (crossline) in km.
  public static class Seismic {
    public double x1,x2,x3;
    public Seismic(double x2, double x3) {
      this(0.0,x2,x3);
    }
    public Seismic(double x1, double x2, double x3) {
      this.x1 = x1;
      this.x2 = x2;
      this.x3 = x3;
    }
    public Seismic(Map map) {
      double xe = map.xe;
      double yn = map.yn;
      double ze = map.ze;
      xe -= XLL_FT;
      yn -= YLL_FT;
      x1 = DATUM_FT-ze;
      x2 =  COSS*xe+SINS*yn;
      x3 = -SINS*xe+COSS*yn;
      x1 *= KM_PER_FT;
      x2 *= KM_PER_FT;
      x3 *= KM_PER_FT;
    }
    public Seismic(Resampled r) {
      double x1r = r.x1;
      double x2r = r.x2;
      double x3r = r.x3;
      x1 = x1r;
      x2 = X2CR+COSR*x2r+SINR*x3r;
      x3 = X3CR-SINR*x2r+COSR*x3r;
    }
  }

  // Resampled coordinates x1 (depth), x2 (inline) and x3 (crossline) in km.
  public static class Resampled {
    public double x1,x2,x3;
    public Resampled(double x2, double x3) {
      this(0.0,x2,x3);
    }
    public Resampled(double x1, double x2, double x3) {
      this.x1 = x1;
      this.x2 = x2;
      this.x3 = x3;
    }
    public Resampled(Seismic s) {
      double x1s = s.x1;
      double x2s = s.x2-X2CR;
      double x3s = s.x3-X3CR;
      x1 = x1s;
      x2 = COSR*x2s-SINR*x3s;
      x3 = SINR*x2s+COSR*x3s;
    }
    public Resampled(Map m) {
      this(new Seismic(m));
    }
  }

  public static void fromMapToSeismic(
    double[] xe, double yn[],
    double[] x2, double x3[])
  {
    int n = xe.length;
    for (int i=0; i<n; ++i) {
      Seismic s = new Seismic(new Map(xe[i],yn[i]));
      x2[i] = s.x2;
      x3[i] = s.x3;
    }
  }

  public static void fromMapToSeismic(
    double[] xe, double yn[], double[] ze,
    double[] x1, double[] x2, double x3[])
  {
    int n = xe.length;
    for (int i=0; i<n; ++i) {
      Seismic s = new Seismic(new Map(xe[i],yn[i],ze[i]));
      x1[i] = s.x1;
      x2[i] = s.x2;
      x3[i] = s.x3;
    }
  }
}
