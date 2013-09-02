package nor;

import static java.lang.Math.*;

import edu.mines.jtk.dsp.Sampling;

/**
 * Coordinate transforms for Norne data.
 * Coordinates in many Norne files are map coordinates (xe,yn,ze), where 
 * xe is easting, yn is northing, and ze is elevation, all measured in m.
 * In some files depths zd are specified instead of elevations ze.
 * 
 * This class defines an additional CSM system of seismic-survey coordinates
 * (x1,x2,x3) = (time/depth,inline,crossline), with units of distance in km.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.31
 */
public class Coordinates {

  public static final Sampling SX = new Sampling( 321,0.0125,0.000);
  public static final Sampling SY = new Sampling(1001,0.0125,0.000);
  public static final Sampling SZ = new Sampling(1001,0.0040,0.000);
  public static final Sampling ST = new Sampling(1001,0.0040,0.000);

  // Seismic datum elevation measured in m. Zero time or depth in seismic 
  // or resampled coordinates corresponds to this elevation in map coordinates.
  public static final double DATUM_M = 0.0; // sea level

  // Conversion m <=> km.
  static final double KM_PER_M = 0.001;
  static final double M_PER_KM = 1.0/KM_PER_M;

  // Eastings and northings - (xe,yn) coordinates - of seismic rectangle.
  // Coordinates are measured in m, obtained from SEG-Y headers.
  // In the CSM files, n2 = 1001, n3 = 321.
  static final double XLL_M=456302.0, YLL_M=7317354.9; // (i2,i3) = (0,0)
  static final double XLR_M=464634.0, YLR_M=7326673.9; // (i2,i3) = (n2-1,0)
  static final double XUR_M=461652.0, YUR_M=7329340.2; // (i2,i3) = (n2-1,n3-1)
  static final double XUL_M=453320.0, YUL_M=7320021.1; // (i2,i3) = (0,n3-1)

  // Rotation angle of Norne seismic rectangle. This is the angle between
  // the CSM seismic x2 (inline) axis and the map xe (easting) axis,
  // measured CCW in the map (xe,yn) coordinate system.
  static final double PHIC = atan2(YLR_M-YLL_M,XLR_M-XLL_M);
  static final double COSC = cos(PHIC);
  static final double SINC = sin(PHIC);

  // Norne coordinates xe (easting), yn (northing), and ze (elevation) in m.
  public static class Norne {
    public double xe,yn,ze;
    public Norne(double xe, double yn) {
      this(xe,yn,0.0);
    }
    public Norne(double xe, double yn, double ze) {
      this.xe = xe;
      this.yn = yn;
      this.ze = ze;
    }
    public Norne(Csm csm) {
      double x1 = csm.x1;
      double x2 = csm.x2;
      double x3 = csm.x3;
      x1 *= M_PER_KM;
      x2 *= M_PER_KM;
      x3 *= M_PER_KM;
      xe = COSC*x2-SINC*x3;
      yn = SINC*x2+COSC*x3;
      ze = -x1;
      xe += XLL_M;
      yn += YLL_M;
      ze += DATUM_M;
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
    public Csm(Norne norne) {
      double xe = norne.xe;
      double yn = norne.yn;
      double ze = norne.ze;
      xe -= XLL_M;
      yn -= YLL_M;
      x1 = DATUM_M-ze;
      x2 =  COSC*xe+SINC*yn;
      x3 = -SINC*xe+COSC*yn;
      x1 *= KM_PER_M;
      x2 *= KM_PER_M;
      x3 *= KM_PER_M;
    }
  }
}
