package tp;

import static java.lang.Math.*;

/**
 * Coordinate transforms for Teapot Dome data.
 * Coordinates in many data files are map coordinates (xe,yn,ze), where 
 * xe is easting, yn is northing, and ze is elevation, all measured in ft.
 * In some files depths zd are specified instead of elevations ze.
 * 
 * This class defines two additional coordinate systems, one for DOE's
 * seismic coordinates (x1,x2,x3) = (time/depth,inline,crossline) and
 * another for CSM's resampled (translated and rotated) seismic coordinates. 
 * The CSM resampled coordinate system enables many of the dead traces in 
 * the DOE seismic survey to be excluded while maintaining regular sampling 
 * of the 3D seismic image.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.19
 */
public class Coordinates {

  // Seismic datum elevation measured in ft. Zero time or depth in seismic 
  // or resampled coordinates corresponds to this elevation in map coordinates.
  public static final double DATUM_FT = 5500.0;

  // Conversion ft <=> km.
  static final double KM_PER_FT = 0.0003048;
  static final double FT_PER_KM = 1.0/KM_PER_FT;

  // Eastings and northings - (xe,yn) coordinates - of seismic rectangle.
  // Coordinates are measured in ft, obtained from teapot_3d_load.doc.
  static final double XLL_FT=788937.0, YLL_FT=938846.0; // (i2,i3) = (0,0)
  static final double XLR_FT=809502.0, YLR_FT=939334.0; // (i2,i3) = (n2-1,0)
  static final double XUR_FT=808604.0, YUR_FT=977163.0; // (i2,i3) = (n2-1,n3-1)
  static final double XUL_FT=788039.0, YUL_FT=976675.0; // (i2,i3) = (0,n3-1)

  // Rotation angle of DOE seismic rectangle. This is the angle between
  // the DOE seismic x2 (inline) axis and the map xe (easting) axis,
  // measured CCW in the map (xe,yn) coordinate system.
  static final double PHID = atan2(YLR_FT-YLL_FT,XLR_FT-XLL_FT);
  static final double COSD = cos(PHID);
  static final double SIND = sin(PHID);

  // CSM resampled (x2,x3) coordinates to exclude most dead traces.
  // This is a translation and rotation. The rotation center (in km) 
  // and angle (in radians) were determined visually. DOE coordinates 
  // (x2d,x3d) and CSM coordinates (x2c,x3c) are measured in km and
  // related by:
  //   x2d = X2DC + x2c*cos(PHIC) - x3c*sin(PHIC)
  //   x3d = X3DC + x2c*sin(PHIC) + x3c*cos(PHIC)
  //static final double X2DC = 4.192; // km
  //static final double X3DC = 0.935; // km
  //static final double PHIC = 0.485364; // radians = 27.8093 degrees
  //static final double COSC = cos(PHIC);
  //static final double SINC = sin(PHIC);
  static final double X2DC =  0.039880; // km
  static final double X3DC =  8.807096; // km
  static final double PHIC = -1.085432; // radians = -62.1907 degrees
  static final double COSC = cos(PHIC);
  static final double SINC = sin(PHIC);

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
    public Map(Doe doe) {
      double x1 = doe.x1;
      double x2 = doe.x2;
      double x3 = doe.x3;
      x1 *= FT_PER_KM;
      x2 *= FT_PER_KM;
      x3 *= FT_PER_KM;
      xe = COSD*x2-SIND*x3;
      yn = SIND*x2+COSD*x3;
      ze = -x1;
      xe += XLL_FT;
      yn += YLL_FT;
      ze += DATUM_FT;
    }
    public Map(Csm csm) {
      this(new Doe(csm));
    }
  }

  /**
   * DOE coordinates x1 (time/depth), x2 (inline) and x3 (crossline) in km.
   */
  public static class Doe {
    public double x1,x2,x3;
    public Doe(double x2, double x3) {
      this(0.0,x2,x3);
    }
    public Doe(double x1, double x2, double x3) {
      this.x1 = x1;
      this.x2 = x2;
      this.x3 = x3;
    }
    public Doe(Map map) {
      double xe = map.xe;
      double yn = map.yn;
      double ze = map.ze;
      xe -= XLL_FT;
      yn -= YLL_FT;
      x1 = DATUM_FT-ze;
      x2 =  COSD*xe+SIND*yn;
      x3 = -SIND*xe+COSD*yn;
      x1 *= KM_PER_FT;
      x2 *= KM_PER_FT;
      x3 *= KM_PER_FT;
    }
    public Doe(Csm csm) {
      double x1c = csm.x1;
      double x2c = csm.x2;
      double x3c = csm.x3;
      x1 = x1c;
      x2 = X2DC+COSC*x2c-SINC*x3c;
      x3 = X3DC+SINC*x2c+COSC*x3c;
    }
  }

  /**
   * CSM coordinates x1 (time/depth), x2 (inline) and x3 (crossline) in km.
   */
  public static class Csm {
    public double x1,x2,x3;
    public Csm(double x2, double x3) {
      this(0.0,x2,x3);
    }
    public Csm(double x1, double x2, double x3) {
      this.x1 = x1;
      this.x2 = x2;
      this.x3 = x3;
    }
    public Csm(Doe doe) {
      double x1d = doe.x1;
      double x2d = doe.x2-X2DC;
      double x3d = doe.x3-X3DC;
      x1 = x1d;
      x2 =  COSC*x2d+SINC*x3d;
      x3 = -SINC*x2d+COSC*x3d;
    }
    public Csm(Map map) {
      this(new Doe(map));
    }
  }
}
