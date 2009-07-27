package gp210;

import java.awt.*;
import java.awt.image.IndexColorModel;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import javax.swing.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.SibsonInterpolator2;
import edu.mines.jtk.mesh.TriMesh;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.sgl.test.TestFrame;
import static edu.mines.jtk.util.ArrayMath.*;

public class Mag {

  private static float[][] readData() {

    // Scan data into one big list of doubles.
    // Find min latitude, longitude, elevation.
    double latmin = Double.MAX_VALUE;
    double lonmin = Double.MAX_VALUE;
    double elemin = Double.MAX_VALUE;
    ArrayList<Double> al = new ArrayList<Double>();
    try {
      Scanner s = new Scanner(new File("src/gp210/magneto0.txt"));
      while (s.hasNextLine()) {
        s.nextInt(); // skip sample index
        double lat = s.nextDouble(); // latitude
        double lon = s.nextDouble(); // longitude
        double ele = s.nextDouble(); // elevation
        latmin = min(latmin,lat);
        lonmin = min(lonmin,lon);
        elemin = min(elemin,ele);
        al.add(lat);
        al.add(lon);
        al.add(ele);
        s.nextDouble(); // skip
        s.nextDouble(); // skip
        s.nextDouble(); // skip
        al.add(s.nextDouble()); // mx
        al.add(s.nextDouble()); // my
        al.add(s.nextDouble()); // mz
        al.add(s.nextDouble()); // temp
        al.add(s.nextDouble()); // mx
        al.add(s.nextDouble()); // my
        al.add(s.nextDouble()); // mz
        al.add(s.nextDouble()); // temp
        s.nextLine();
      }
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }
    System.out.println("latmin="+latmin);
    System.out.println("lonmin="+lonmin);
    System.out.println("elemin="+elemin);

    // Meters per degree of latitude,longitude at minimum latitude.
    // (See Wikipedia entry for geographical coordinate system.)
    double a = 6378137.0;
    double b = 6356752.0;
    double c = cos(toRadians(latmin));
    double s = sin(toRadians(latmin));
    double aa = a*a, bb = b*b, cc = c*c, ss = s*s;
    double aaaa = aa*aa, bbbb = bb*bb;
    double lonToMeters = PI/180.0*c*sqrt((aaaa*cc+bbbb*ss)/(aa*cc+bb*ss));
    double latToMeters = 111300.0;
    System.out.println("lonToMeters="+lonToMeters);
    System.out.println("latToMeters="+latToMeters);

    // Split data into separate arrays of floats.
    // Convert global (lon,lat,ele) to local (x,y,z).
    // Compute gradient of magnitude of magnetic field.
    int n = al.size()/11; // kept 11 measurements above
    float[] x = new float[n];
    float[] y = new float[n];
    float[] z = new float[n];
    float[] g = new float[n];
    for (int i=0,j=0; i<n; ++i,j+=11) {
      double lat = al.get(j+0); // latitude
      double lon = al.get(j+1); // longitude
      double ele = al.get(j+2); // elevation
      double mx1 = al.get(j+3);
      double my1 = al.get(j+4);
      double mz1 = al.get(j+5);
      double mx2 = al.get(j+7);
      double my2 = al.get(j+8);
      double mz2 = al.get(j+9);
      double m1 = sqrt(mx1*mx1+my1*my1+mz1*mz1);
      double m2 = sqrt(mx2*mx2+my2*my2+mz2*mz2);
      x[i] = (float)((lon-lonmin)*lonToMeters);
      y[i] = (float)((lat-latmin)*latToMeters);
      z[i] = (float)(ele-elemin);
      g[i] = (float)(m1-m2);
    }
    float[][] data = new float[][]{x,y,z,g};
    return data;
  }

  private static Sampling[] makeSamplings(float[] x, float[] y) {
    float xmin = min(x);
    float xmax = max(x);
    float ymin = min(y);
    float ymax = max(y);
    double fx = xmin;
    double fy = ymin;
    //double dx = 0.5;
    double dx = 0.2;
    //double dx = 0.1;
    //double dx = 0.05; // makes discrete Sibson really slow!
    double dy = dx;
    int nx = 2+(int)((xmax-xmin)/dx);
    int ny = 2+(int)((ymax-ymin)/dy);
    System.out.println("nx="+nx+" dx="+dx+" fx="+fx);
    System.out.println("ny="+ny+" dy="+dy+" fy="+fy);
    Sampling sx = new Sampling(nx,dx,fx);
    Sampling sy = new Sampling(ny,dy,fy);
    return new Sampling[]{sx,sy};
  }

  private static float[][] interpolateSample(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    int n = x.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = fillfloat(-1.0f,nx,ny);
    for (int i=0; i<n; ++i) {
      int ix = sx.indexOfNearest(x[i]);
      int iy = sy.indexOfNearest(y[i]);
      zi[iy][ix] = z[i];
    }
    return zi;
  }

  private static float[][] interpolateSimple(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    int n = x.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    for (int iy=0; iy<ny; ++iy) {
      double yi = sy.getValue(iy);
      for (int ix=0; ix<nx; ++ix) {
        double xi = sx.getValue(ix);
        double num = 0.0;
        double den = 0.0;
        for (int j=0; j<n; ++j) {
          double dx = x[j]-xi;
          double dy = y[j]-yi;
          double ds = dx*dx+dy*dy;
          double wj = 1.0/(ds*ds);
          num += wj*z[j];
          den += wj;
        }
        zi[iy][ix] = (float)(num/den);
      }
    }
    return zi;
  }

  private static float[][] interpolateNearest(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    int n = x.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    TriMesh mesh = new TriMesh();
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    for (int i=0; i<n; ++i) {
      TriMesh.Node node = new TriMesh.Node(x[i],y[i]);
      if (mesh.addNode(node))
        zmap.put(node,new Float(z[i]));
    }
    for (int iy=0; iy<ny; ++iy) {
      float yi = (float)sy.getValue(iy);
      for (int ix=0; ix<nx; ++ix) {
        float xi = (float)sx.getValue(ix);
        TriMesh.Node node = mesh.findNodeNearest(xi,yi);
        zi[iy][ix] = (Float)zmap.get(node);
      }
    }
    return zi;
  }

  private static float[][] interpolateDiscreteSibson(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    DiscreteSibsonInterpolator dsi = new DiscreteSibsonInterpolator(sx,sy);
    dsi.setSmoothingIterations(0);
    return dsi.apply(x,y,z);
  }

  private static float[][] interpolateSibson(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    SibsonInterpolator2 si = new SibsonInterpolator2(z,x,y);
    si.setBounds(sx,sy);
    si.setNullValue(-1.0f);
    return si.interpolate(sx,sy);
  }

  private static float[][] interpolateSibson1(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    SibsonInterpolator2 si = new SibsonInterpolator2(z,x,y);
    si.setBounds(sx,sy);
    si.setNullValue(-1.0f);
    si.setGradientPower(1.0f);
    return si.interpolate(sx,sy);
  }

  private static void plot(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz,
    String title, String png) 
  {
    SimplePlot sp = new SimplePlot();
    //sp.setSize(1065,645);
    sp.setSize(1057,630);
    sp.setFontSize(36);
    sp.setTitle(title);
    sp.setHLabel("Easting (m)");
    sp.setVLabel("Northing (m)");
    sp.addColorBar();
    sp.getPlotPanel().setColorBarWidthMinimum(100);
    sp.setHLimits(sx.getFirst(),sx.getLast());
    sp.setVLimits(sy.getFirst(),sy.getLast());
    if (sz!=null) {
      PixelsView iv = sp.addPixels(sx,sy,sz);
      iv.setInterpolation(PixelsView.Interpolation.NEAREST);
      iv.setColorModel(makeColorModel());
      iv.setClips(-1.0f,125f);
    }
    if (x!=null && y!=null) {
      PointsView dv = sp.addPoints(x,y);
      dv.setLineStyle(PointsView.Line.NONE);
      dv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
      dv.setMarkSize(2.0f);
    }
    sp.paintToPng(300,6,png+".png");
  }

  private static IndexColorModel makeColorModel() {
    IndexColorModel icm = ColorMap.JET;
    int n = 256;
    byte[] r = new byte[n];
    byte[] g = new byte[n];
    byte[] b = new byte[n];
    byte[] a = new byte[n];
    icm.getReds(r);
    icm.getGreens(g);
    icm.getBlues(b);
    icm.getAlphas(a);
    r[0] = (byte)255;
    g[0] = (byte)255;
    b[0] = (byte)255;
    a[0] = (byte)0;
    icm = new IndexColorModel(8,n,r,g,b,a);
    return icm;
  }

  private static void plot3d(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz)
  {
    PointGroup pg = makePointGroup(x,y,z);
    TriangleGroup tg = makeTriangleGroup(sx,sy,sz);
    World world = new World();
    world.addChild(pg);
    world.addChild(tg);
    TestFrame frame = new TestFrame(world);
    OrbitView view = frame.getOrbitView();
    view.setAxesOrientation(AxesOrientation.XOUT_YRIGHT_ZUP);
    frame.setSize(new Dimension(1200,800));
    frame.setVisible(true);
  }
  private static PointGroup makePointGroup(float[] x, float[] y, float[] z) {
    int n = x.length;
    float[] xyz = new float[3*n];
    copy(n,0,1,x,0,3,xyz);
    copy(n,0,1,y,1,3,xyz);
    copy(n,0,1,z,2,3,xyz);
    float size = 0.2f;
    PointGroup pg = new PointGroup(size,xyz);
    StateSet states = new StateSet();
    ColorState cs = new ColorState();
    cs.setColor(Color.RED);
    states.add(cs);
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    states.add(lms);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setSpecular(Color.WHITE);
    ms.setShininess(100.0f);
    states.add(ms);
    pg.setStates(states);
    return pg;
  }
  private static TriangleGroup makeTriangleGroup(
    Sampling sx, Sampling sy, float[][] sz) 
  {
    sz = transpose(sz);
    TriangleGroup tg = new TriangleGroup(true,sx,sy,sz);
    StateSet states = new StateSet();
    ColorState cs = new ColorState();
    cs.setColor(Color.LIGHT_GRAY);
    states.add(cs);
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    states.add(lms);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setSpecular(Color.WHITE);
    ms.setShininess(100.0f);
    states.add(ms);
    tg.setStates(states);
    return tg;
  }

  private static float[][] mask(float maskValue, float[][] m, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) 
        y[i2][i1] = (m[i2][i1]==maskValue)?maskValue:x[i2][i1];
    return y;
  }

  private static void interpolate() {
    float[][] data = readData();
    float[] x = data[0], y = data[1], z = data[2]; //, g = data[3];
    mul(100.0f,z,z);
    Sampling[] s = makeSamplings(x,y);
    Sampling sx = s[0], sy = s[1];

    //float[][] z0 = interpolateSample(x,y,z,sx,sy);
    //float[][] z1 = interpolateSimple(x,y,z,sx,sy);
    float[][] z2 = interpolateSibson(x,y,z,sx,sy);
    //z1 = mask(-1.0f,z2,z1);
    //float[][] z3 = interpolateNearest(x,y,z,sx,sy);
    //float[][] z4 = interpolateDiscreteSibson(x,y,z,sx,sy);
    //plot(null,y,z,sx,sy,z0,"Elevation (cm)","elevm");
    //plot(null,y,z,sx,sy,z1,"Elevation (cm)","elevs");
    //plot(null,y,z,sx,sy,z2,"Natural neighbor elevation (cm)","elev2a");
    //plot(x,y,z,sx,sy,z2,"Elevation (cm)","elevi");
    //plot(null,y,z,sx,sy,z2,"Elevation (cm)","elevmi");
    //plot(x,y,z,sx,sy,z3,"Nearest neighbor elevation (cm)","elev3");
    //plot(null,y,z,sx,sy,z4,"Discrete Sibson elevation (cm)","elev4");

    //float[][] g1 = interpolateSimple(x,y,g,sx,sy);
    //float[][] g2 = interpolateSibson(x,y,g,sx,sy);
    //float[][] g3 = interpolateNearest(x,y,g,sx,sy);
    //float[][] g4 = interpolateDiscreteSibson(x,y,g,sx,sy);
    //plot(x,y,g,sx,sy,g1,"Simple magnetic gradient","grad1");
    //plot(x,y,g,sx,sy,g2,"Sibson magnetic gradient","grad2");
    //plot(x,y,g,sx,sy,g3,"Nearest magnetic gradient","grad3");
    //plot(x,y,g,sx,sy,g4,"Discrete Sibson magnetic gradient","grad4");

    //plot3d(x,y,z,sx,sy,z1);
    plot3d(x,y,z,sx,sy,z2);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        interpolate();
      }
    });
  }
}
