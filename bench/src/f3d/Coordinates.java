package f3d;

import static java.lang.Math.*;

/**
 * Coordinate transforms for F3 data from opendtect.org (ODT).
 * Coordinates in many ODT files are map coordinates (xe,yn,ze), where 
 * xe is easting, yn is northing, and ze is elevation, all measured in m.
 * In some files depths zd are specified instead of elevations ze.
 * 
 * This class defines an additional CSM system of seismic-survey coordinates
 * (x1,x2,x3) = (time/depth,inline,crossline), with units of distance in km.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.12.29
 */
public class Coordinates {

  // Seismic datum elevation measured in m. Zero time or depth in seismic 
  // or resampled coordinates corresponds to this elevation in map coordinates.
  public static final double DATUM_M = 0.0; // sea level

  // Conversion m <=> km.
  static final double KM_PER_M = 0.001;
  static final double M_PER_KM = 1.0/KM_PER_M;

  // Eastings and northings - (xe,yn) coordinates - of seismic rectangle.
  // Coordinates are measured in m, obtained from SEG-Y headers.
  // In the CSM subset, n2 = 951, n3 = 591.
  static final double XLL_M=605835.5, YLL_M=6073556.4; // (i2,i3) = (0,0)
  static final double XLR_M=629576.3, YLR_M=6074219.9; // (i2,i3) = (n2-1,0)
  static final double XUR_M=629164.4, YUR_M=6088963.8; // (i2,i3) = (n2-1,n3-1)
  static final double XUL_M=605423.7, YUL_M=6088300.6; // (i2,i3) = (0,n3-1)

  // Rotation angle of ODT seismic rectangle. This is the angle between
  // the CSM seismic x2 (inline) axis and the map xe (easting) axis,
  // measured CCW in the map (xe,yn) coordinate system.
  static final double PHIC = atan2(YLR_M-YLL_M,XLR_M-XLL_M);
  static final double COSC = cos(PHIC);
  static final double SINC = sin(PHIC);

  // ODT coordinates xe (easting), yn (northing), and ze (elevation) in m.
  public static class Odt {
    public double xe,yn,ze;
    public Odt(double xe, double yn) {
      this(xe,yn,0.0);
    }
    public Odt(double xe, double yn, double ze) {
      this.xe = xe;
      this.yn = yn;
      this.ze = ze;
    }
    public Odt(Csm csm) {
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
    public Csm(Odt odt) {
      double xe = odt.xe;
      double yn = odt.yn;
      double ze = odt.ze;
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
