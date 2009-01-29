package gp210;

import java.awt.*;
import java.io.*;
import java.util.*;
import javax.swing.SwingUtilities;
import static java.lang.Math.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mesh.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.sgl.test.TestFrame;
import static edu.mines.jtk.ogl.Gl.*;

public class Data {

  private static float[][] readData() {

    // Scan data into one big list of floats.
    ArrayList<Float> al = new ArrayList<Float>();
    try {
      Scanner s = new Scanner(new File("src/gp210/"+DATA_NAME+".txt"));
      while (s.hasNextLine()) {
        al.add(s.nextFloat());
        al.add(s.nextFloat());
        al.add(s.nextFloat());
        s.nextLine();
      }
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);
    }

    // Split data into separate arrays of floats.
    int n = al.size()/3;
    float[] x = new float[n];
    float[] y = new float[n];
    float[] z = new float[n];
    for (int i=0,j=0; i<n; ++i,j+=3) {
      x[i] = al.get(j+0);
      y[i] = al.get(j+1);
      z[i] = al.get(j+2);
    }
    float[][] data = new float[][]{x,y,z};
    return data;
  }

  private static Sampling[] makeSamplings(float[] x, float[] y) {
    float xmin = Array.min(x);
    float xmax = Array.max(x);
    float ymin = Array.min(y);
    float ymax = Array.max(y);
    int nx = 401;
    int ny = 401;
    double fx = 0.000;
    double fy = 0.000;
    double dx = 1.0/(nx-1);
    double dy = 1.0/(ny-1);
    Sampling sx = new Sampling(nx,dx,fx);
    Sampling sy = new Sampling(ny,dy,fy);
    return new Sampling[]{sx,sy};
  }

  private static TriMesh makeMesh(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy, boolean extrap)
  {
    int n = x.length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    TriMesh mesh = new TriMesh();
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    for (int i=0; i<n; ++i) {
      TriMesh.Node node = new TriMesh.Node(x[i],y[i]);
      mesh.addNode(node);
      zmap.put(node,new Float(z[i]));
    }
    extrap = true;
    if (extrap) {
      for (int iy=0; iy<ny; iy+=ny-1) {
        float yi = (float)sy.getValue(iy);
        yi += (iy==0)?-0.005f:0.005f;
        for (int ix=0; ix<nx; ix+=nx-1) {
          float xi = (float)sx.getValue(ix);
          xi += (ix==0)?-0.005f:0.005f;
          TriMesh.Node near = mesh.findNodeNearest(xi,yi);
          TriMesh.Node node = new TriMesh.Node(xi,yi);
          mesh.addNode(node);
          zmap.put(node,(Float)zmap.get(near));
        }
      }
    }
    return mesh;
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

  private static float[][] interpolatePlanar(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    TriMesh mesh = makeMesh(x,y,z,sx,sy,true);
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    TriMesh.Tri triLast = null;
    float xa = 0.0f, ya = 0.0f, za = 0.0f, ax = 0.0f, ay = 0.0f;
    float xb,yb,zb,xc,yc,zc;
    float znull = 0.5f*(Array.min(z)+Array.max(z));
    for (int iy=0; iy<ny; ++iy) {
      float yi = (float)sy.getValue(iy);
      for (int ix=0; ix<nx; ++ix) {
        float xi = (float)sx.getValue(ix);
        TriMesh.PointLocation pl = mesh.locatePoint(xi,yi);
        TriMesh.Tri tri = pl.tri();
        if (pl.isInside()) {
          if (tri!=triLast && tri!=null) {
            triLast = tri;
            TriMesh.Node na = tri.nodeA();
            TriMesh.Node nb = tri.nodeB();
            TriMesh.Node nc = tri.nodeC();
            xa = na.x(); ya = na.y(); za = (Float)zmap.get(na);
            xb = nb.x(); yb = nb.y(); zb = (Float)zmap.get(nb);
            xc = nc.x(); yc = nc.y(); zc = (Float)zmap.get(nc);
            xb -= xa; yb -= ya; zb -= za;
            xc -= xa; yc -= ya; zc -= za;
            float odet = 1.0f/(xb*yc-xc*yb);
            ax = (yc*zb-yb*zc)*odet;
            ay = (xb*zc-xc*zb)*odet;
          }
          zi[iy][ix] = za+ax*(xi-xa)+ay*(yi-ya);
        } else {
          zi[iy][ix] = znull;
        }
      }
    }
    return zi;
  }

  private static float[][] interpolateNearest(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    TriMesh mesh = makeMesh(x,y,z,sx,sy,false);
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
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

  private static float[][] interpolateNatural(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] zi = new float[ny][nx];
    TriMesh mesh = makeMesh(x,y,z,sx,sy,true);
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    float znull = 0.5f*(Array.min(z)+Array.max(z));
    for (int iy=0; iy<ny; ++iy) {
      float yi = (float)sy.getValue(iy);
      for (int ix=0; ix<nx; ++ix) {
        float xi = (float)sx.getValue(ix);
        zi[iy][ix] = mesh.interpolateSibson(xi,yi,zmap,znull);
      }
    }
    return zi;
  }

  private static void plotMesh(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, 
    boolean labels, boolean tris, boolean polys,
    String title, String png) 
  {
    SimplePlot sp = new SimplePlot();
    sp.setSize(PLOT_MESH_WIDTH,PLOT_HEIGHT);
    sp.setTitle(title);
    sp.setHLabel("x (m)");
    sp.setVLabel("y (m)");
    PlotPanel pp = sp.getPlotPanel();
    pp.setLimits(-0.02,-0.02,1.02,1.02);
    TriMesh mesh = makeMesh(x,y,z,sx,sy,false);
    TriMeshView tv = new TriMeshView(mesh);
    tv.setTrisVisible(tris);
    tv.setPolysVisible(polys);
    tv.setMarkColor(Color.BLACK);
    tv.setTriColor(Color.RED);
    tv.setPolyColor(Color.BLUE);
    pp.getTile(0,0).addTiledView(tv);
    PointsView dv;
    if (labels) {
      dv = pp.addPoints(x,y,z);
    } else {
      dv = pp.addPoints(x,y);
    }
    dv.setLineStyle(PointsView.Line.NONE);
    dv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    dv.setMarkSize(4.0f);
    dv.setTextFormat("%2.0f");
    sp.paintToPng(300,6,png+"2.png");
  }

  private static void plot(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz,
    String title, String png) 
  {
    SimplePlot sp = new SimplePlot();
    sp.setSize(PLOT_WIDTH,PLOT_HEIGHT);
    sp.setTitle(title);
    sp.setHLabel("x (m)");
    sp.setVLabel("y (m)");
    sp.addColorBar();
    PlotPanel pp = sp.getPlotPanel();
    pp.setColorBarWidthMinimum(80);
    pp.setLimits(-0.02,-0.02,1.02,1.02);
    PixelsView iv = sp.addPixels(sx,sy,sz);
    iv.setInterpolation(PixelsView.Interpolation.NEAREST);
    iv.setColorModel(ColorMap.JET);
    PointsView dv = pp.addPoints(x,y,z);
    dv.setLineStyle(PointsView.Line.NONE);
    dv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    dv.setMarkSize(2.0f);
    dv.setTextFormat("%2.0f");
    sp.paintToPng(300,6,png+"_"+DATA_NAME+".png");
  }

  private static void plot3d(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz)
  {
    sz = Array.mul(0.01f,sz);
    z = Array.mul(0.01f,z);
    PointGroup pg = makePointGroup(x,y,z);
    TriangleGroup tg = makeTriangleGroup(sx,sy,sz);
    World world = new World();
    world.addChild(pg);
    world.addChild(tg);
    TestFrame frame = new TestFrame(world);
    OrbitView view = frame.getOrbitView();
    view.setScale(2.0f);
    view.setAxesOrientation(View.AxesOrientation.XOUT_YRIGHT_ZUP);
    frame.setSize(new Dimension(1200,800));
    frame.setVisible(true);
  }
  private static PointGroup makePointGroup(float[] x, float[] y, float[] z) {
    int n = x.length;
    float[] xyz = new float[3*n];
    Array.copy(n,0,1,x,0,3,xyz);
    Array.copy(n,0,1,y,1,3,xyz);
    Array.copy(n,0,1,z,2,3,xyz);
    float size = 0.01f;
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
    sz = Array.transpose(sz);
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

  private static void interpolate() {
    float[][] data = readData();
    float[] x = data[0], y = data[1], z = data[2];
    Sampling[] s = makeSamplings(x,y);
    Sampling sx = s[0], sy = s[1];
    //plotMesh(x,y,z,sx,sy,true,false,false,"Scattered data","sd");
    //plotMesh(x,y,z,sx,sy,false,true,false,"Delaunay triangles","dt");
    //plotMesh(x,y,z,sx,sy,false,false,true,"Voronoi polygons","vp");
    //plotMesh(x,y,z,sx,sy,false,true,true,"triangles & polygons","tp");
    //float[][] z1 = interpolateSimple(x,y,z,sx,sy);
    //float[][] z2 = interpolatePlanar(x,y,z,sx,sy);
    //float[][] z3 = interpolateNearest(x,y,z,sx,sy);
    float[][] z4 = interpolateNatural(x,y,z,sx,sy);
    //plot(x,y,z,sx,sy,z1,"Distance weighted","dw");
    //plot(x,y,z,sx,sy,z2,"Planar triangles","pt");
    //plot(x,y,z,sx,sy,z3,"Nearest neighbor","ne");
    plot(x,y,z,sx,sy,z4,"Natural neighbor","na");
    //plot3d(x,y,z,sx,sy,z1);
    //plot3d(x,y,z,sx,sy,z2);
    //plot3d(x,y,z,sx,sy,z3);
    //plot3d(x,y,z,sx,sy,z4);
  }

  //public static final String DATA_NAME = "data1";
  public static final String DATA_NAME = "data2";
  public static final int PLOT_HEIGHT = 785;
  public static final int PLOT_WIDTH = 875;
  public static final int PLOT_MESH_WIDTH = PLOT_WIDTH-132;

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        interpolate();
      }
    });
  }
}
