/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

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

/**
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.02.05
 */
public class SambridgeTest {

  private static float[][] makeData() {
    Random r = new Random();
    int seed = r.nextInt();
    System.out.println("seed="+seed);
    r.setSeed(seed);
    int n = 5000;
    float[] x = new float[n];
    float[] y = new float[n];
    float[] z = new float[n];
    for (int i=0,j=0; i<n; ++i,j+=3) {
      x[i] = r.nextFloat();
      y[i] = r.nextFloat();
      z[i] = 50.0f+25.0f*(float)(sin(PI*x[i])*sin(PI*y[i]));
    }
    float[][] data = new float[][]{x,y,z};
    return data;
  }

  private static Sampling[] makeSamplings(float[] x, float[] y) {
    int nx = 1001;
    int ny = 1001;
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

  private static float[][] interpolate(int method,
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
        if (method==1) {
          zi[iy][ix] = mesh.interpolateSibson(xi,yi,zmap,znull);
        } else if (method==2) {
          zi[iy][ix] = mesh.interpolateSambridge(xi,yi,zmap,znull);
        }
      }
    }
    return zi;
  }

  private static float[][] interpolateError(
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
        float z1 = mesh.interpolateSibson(xi,yi,zmap,znull);
        float z2 = mesh.interpolateSambridge(xi,yi,zmap,znull);
        zi[iy][ix] = z2-z1;
      }
    }
    return zi;
  }

  private static void interpolateErrorUlp(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    //System.out.print("makeMesh ... ");
    TriMesh mesh = makeMesh(x,y,z,sx,sy,true);
    //System.out.println("done");
    TriMesh.NodePropertyMap zmap = mesh.getNodePropertyMap("z");
    float znull = 0.5f*(Array.min(z)+Array.max(z));
    TriMesh.NodeIterator ni = mesh.getNodes();
    TriMesh.NodeList nl = new TriMesh.NodeList();
    long nxy = 0;
    while (ni.hasNext()) {
      TriMesh.Node na = ni.next();
      float xa = na.x();
      float ya = na.y();
      if (xa<0.1f || xa>0.9f) continue;
      if (ya<0.1f || ya>0.9f) continue;
      nl.clear();
      mesh.getNodeNabors(na,nl);
      int nnabor = nl.nnode();
      TriMesh.Node[] nabors = nl.nodes();
      for (int inabor=0; inabor<nnabor; ++inabor) {
        TriMesh.Node nb = nabors[inabor];
        float xb = nb.x();
        float yb = nb.y();
        if (xb<0.1f || xb>0.9f) continue;
        if (yb<0.1f || yb>0.9f) continue;
        int ns = 1000;
        float dx = (xb-xa)/(float)(ns-1);
        float dy = (yb-ya)/(float)(ns-1);
        for (int is=0; is<ns; ++is) {
          float xi = xa+(float)is*dx;
          float yi = ya+(float)is*dy;
          ++nxy;
          //if (nxy%100000==0)
          //  System.out.println("nxy="+nxy+" xi="+xi+" yi="+yi);
          float z1 = mesh.interpolateSibson(xi,yi,zmap,znull);
          float z2 = mesh.interpolateSambridge(xi,yi,zmap,znull);
          float ze = 100.0f*abs(z2-z1)/z1;
          if (ze>0.1f) {
            System.out.println("-----------------------------------------");
            System.out.println("xi="+xi);
            System.out.println("yi="+yi);
            System.out.println("z1="+z1);
            System.out.println("z2="+z2);
            System.out.println("ze="+ze);
            System.out.println("-----------------------------------------");
          }
          if (ze>1.0f)
            System.exit(-1);
        }
      }
    }
    System.out.println("nxy="+nxy);
  }

  private static long nx(float xmin, float xmax) {
    long n = 0;
    for (float x=xmin; x<=xmax; x+=ulp(x),++n)
      ;
    return n;
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
    //sp.paintToPng(300,6,png+".png");
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
    //sp.paintToPng(300,6,png+".png");
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
    float[][] data = makeData();
    float[] x = data[0], y = data[1], z = data[2];
    Sampling[] s = makeSamplings(x,y);
    Sampling sx = s[0], sy = s[1];
    //plotMesh(x,y,z,sx,sy,true,false,false,"Scattered data","sd");
    //plotMesh(x,y,z,sx,sy,false,true,false,"Delaunay triangles","dt");
    //plotMesh(x,y,z,sx,sy,false,false,true,"Voronoi polygons","vp");
    //plotMesh(x,y,z,sx,sy,false,true,true,"triangles & polygons","tp");
    //float[][] z1 = interpolate(1,x,y,z,sx,sy);
    //float[][] z2 = interpolate(2,x,y,z,sx,sy);
    //float[][] ze = interpolateError(x,y,z,sx,sy);
    //System.out.println("error = "+Array.max(Array.abs(ze)));
    //plot(x,y,z,sx,sy,ze,"error","se");
    //plot(x,y,z,sx,sy,z1,"Sibson interpolation","si");
    //plot(x,y,z,sx,sy,z2,"Sambridge interpolation","sa");
    //plot3d(x,y,z,sx,sy,z1);
    //plot3d(x,y,z,sx,sy,z2);
    //System.out.println("nx [0:1] = "+nx(0.0f,1.0f));
    //System.out.println("nx [0.4:0.6] = "+nx(0.4f,0.6f));
    interpolateErrorUlp(x,y,z,sx,sy);
  }

  public static final int PLOT_HEIGHT = 785;
  public static final int PLOT_WIDTH = 875;
  public static final int PLOT_MESH_WIDTH = PLOT_WIDTH-132;

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        while(true) {
          System.out.println("#########################################");
          interpolate();
        }
      }
    });
  }
}
