/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import javax.swing.*;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.SibsonInterpolator2;
import edu.mines.jtk.mesh.TriMesh;
import edu.mines.jtk.mesh.TriMeshView;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.sgl.SimpleFrame;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Illustrates significant error in Sambridge's implementation of
 * natural-neighbor interpolation.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.03.28
 */
public class SambridgeError {

  // Parameters that will produce significant error.
  /*
  private static final int NNODE = 5000;
  private static final int SEED = 268695554;
  private static final float XBAD = 0.49457437f;
  private static final float YBAD = 0.7369756f;

  // 6% (near the middle)
  private static final int NNODE = 50;
  private static final int SEED = 1754078249;
  private static final float XBAD = 0.6789419f;
  private static final float YBAD = 0.64597934f;

  // 15.9% (but negative and too near the edge)
  private static final int NNODE = 50;
  private static final int SEED = 282375838;
  private static final float XBAD = 0.76774025f;
  private static final float YBAD = 0.13308403f;

  // 13.5% (near center but negative) z1=73.94699 z2=64.0 ze=13.45
  private static final int NNODE = 50;
  private static final int SEED = 1679885684;
  private static final float XBAD = 0.5697848f;
  private static final float YBAD = 0.49317682f;

  // 7.0% (near center and positive) z1=73.5013 z2=78.676605 ze=7.0411115
  private static final int NNODE = 50;
  private static final int SEED = 674706574;
  private static final float XBAD = 0.57742673f;
  private static final float YBAD = 0.44822684f;

  // 4.2% (near center and positive) z1=71.683334 z2=74.666664 ze=4.161818
  private static final int NNODE = 50;
  private static final int SEED = 853812808;
  private static final float XBAD = 0.4459424f;
  private static final float YBAD = 0.36336002f;
  */

  private static float[] XHULL = {
    0.42878067f,0.59729666f,0.41604030f,0.58106000f,0.59000f};
  private static float[] YHULL = {
    0.46198040f,0.44638836f,0.41495920f,0.49765354f,0.35000f};
  private static float[] ZHULL = {
    74.2031700f,73.5038100f,73.2791300f,74.1930850f,73.3907f};
  private static final float XBAD = 0.57742673f;
  private static final float YBAD = 0.44822684f;
  private static final float XGOOD = 0.56f;
  private static final float YGOOD = 0.41f;
  private static final float XNODE = XBAD;
  private static final float YNODE = YBAD;
  private static final float XMIN = 0.41f;
  private static final float XMAX = 0.61f;
  private static final float YMIN = 0.32f;
  private static final float YMAX = 0.52f;

  /*
  private static float[][] makeDataRandom() {
    Random r = new Random(SEED);
    int n = NNODE;
    float[] x = new float[n];
    float[] y = new float[n];
    float[] z = new float[n];
    for (int i=0; i<n; ++i) {
      x[i] = r.nextFloat();
      y[i] = r.nextFloat();
      z[i] = z(x[i],y[i]);
    }
    return new float[][]{x,y,z};
  }
  */
  private static float[][] makeData() {
    float[] x = XHULL;
    float[] y = YHULL;
    float[] z = ZHULL;
    for (int i=0; i<x.length; ++i)
      System.out.println("x="+x[i]+" y="+y[i]+" z="+z[i]);
    return new float[][]{x,y,z};
  }
  private static float z(float x, float y) {
      return 50.0f+25.0f*(float)(sin(PI*x)*sin(PI*y));
  }

  private static float[][] windowData(
    float xmin, float xmax, float ymin, float ymax, float[][] d)
  {
    float[] x = d[0], y = d[1], z = d[2];
    int n = x.length;
    int m = 0;
    for (int i=0; i<n; ++i) {
      if (xmin<=x[i] && x[i]<=xmax && ymin<=y[i] && y[i]<=ymax) {
        System.out.println("x="+x[i]+" y="+y[i]);
        ++m;
      }
    }
    //m += 4;
    float[] xw = new float[m];
    float[] yw = new float[m];
    float[] zw = new float[m];
    float zmin =  Float.MAX_VALUE;
    float zmax = -Float.MAX_VALUE;
    for (int i=0,j=0; i<n; ++i) {
      if (xmin<=x[i] && x[i]<=xmax && ymin<=y[i] && y[i]<=ymax) {
        xw[j] = x[i];
        yw[j] = y[i];
        zw[j] = z[i];
        if (zw[j]<zmin) zmin = zw[j];
        if (zw[j]>zmax) zmax = zw[j];
        ++j;
      }
    }
    /*
    zmin -= 3.0f*(zmax-zmin);
    xw[m-4] = xmin; yw[m-4] = ymin; zw[m-4] = zmin;
    xw[m-3] = xmin; yw[m-3] = ymax; zw[m-3] = zmin;
    xw[m-2] = xmax; yw[m-2] = ymin; zw[m-2] = zmin;
    xw[m-1] = xmax; yw[m-1] = ymax; zw[m-1] = zmin;
    */
    return new float[][]{xw,yw,zw};
  }

  private static Sampling[] makeSamplings(float[] x, float[] y) {
    int nx = 201;
    int ny = 201;
    double dx = (XMAX-XMIN)/(nx-1);
    double dy = (YMAX-YMIN)/(ny-1);
    double fx,fy;
    for (fx=XNODE; fx>=XMIN+dx; fx-=dx);
    for (fy=YNODE; fy>=YMIN+dy; fy-=dy);
    Sampling sx = new Sampling(nx,dx,fx);
    Sampling sy = new Sampling(ny,dy,fy);
    return new Sampling[]{sx,sy};
  }

  private static TriMesh makeMesh(
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy,
    boolean extrap, boolean withNode)
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
        yi += (iy==0)?-0.001f:0.001f;
        for (int ix=0; ix<nx; ix+=nx-1) {
          float xi = (float)sx.getValue(ix);
          xi += (ix==0)?-0.001f:0.001f;
          TriMesh.Node near = mesh.findNodeNearest(xi,yi);
          TriMesh.Node node = new TriMesh.Node(xi,yi);
          mesh.addNode(node);
          zmap.put(node,(Float)zmap.get(near));
        }
      }
    }
    if (withNode) {
      float xnode = XNODE, ynode = YNODE;
      System.out.println("node at x="+xnode+" y="+ynode);
      TriMesh.Node node = new TriMesh.Node(xnode,ynode);
      mesh.addNode(node);
      zmap.put(node,new Float(z(xnode,ynode)));
      TriMesh.Tri[] tris = mesh.getTriNabors(node);
      for (int itri=0; itri<tris.length; ++itri) {
        double[] c = tris[itri].centerCircle();
        System.out.println("center x="+c[0]+" y="+c[1]);
      }
    }
    return mesh;
  }

  private static float[][] interpolate(int method,
    float[] x, float[] y, float[] z, Sampling sx, Sampling sy)
  {
    SibsonInterpolator2.Method m = (method==1) ?
      SibsonInterpolator2.Method.HALE_LIANG :
      SibsonInterpolator2.Method.WATSON_SAMBRIDGE;
    SibsonInterpolator2 si = new SibsonInterpolator2(m,z,x,y);
    si.setBounds(sx,sy);
    return si.interpolate(sx,sy);
  }

  private static void plotMesh(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, 
    boolean nodes, boolean tris, boolean polys, boolean withNode,
    String title, String png) 
  {
    SimplePlot sp = new SimplePlot();
    sp.setSize(PLOT_WIDTH,PLOT_HEIGHT);
    sp.setTitle(title);
    sp.setHLabel("x");
    sp.setVLabel("y");
    PlotPanel pp = sp.getPlotPanel();
    pp.setLimits(XMIN,YMIN,XMAX,YMAX);
    if (withNode) {
      TriMesh mesh = makeMesh(x,y,z,sx,sy,false,true);
      TriMeshView tv = new TriMeshView(mesh);
      tv.setNodesVisible(nodes);
      tv.setTrisVisible(tris);
      tv.setPolysVisible(polys);
      tv.setMarkColor(Color.RED);
      tv.setTriColor(Color.RED);
      tv.setPolyColor(Color.RED);
      pp.getTile(0,0).addTiledView(tv);
    }
    TriMesh mesh = makeMesh(x,y,z,sx,sy,false,false);
    TriMeshView tv = new TriMeshView(mesh);
    tv.setNodesVisible(nodes);
    tv.setTrisVisible(tris);
    tv.setPolysVisible(polys);
    tv.setMarkColor(Color.BLUE);
    tv.setTriColor(Color.BLUE);
    tv.setPolyColor(Color.BLUE);
    pp.getTile(0,0).addTiledView(tv);
    //sp.paintToPng(300,6,png+".png");
  }

  private static void plot3d(
    float[] x, float[] y, float[] z,
    Sampling sx, Sampling sy, float[][] sz)
  {
    z = mul(0.01f,z);
    sz = mul(0.01f,sz);
    TriMesh mesh = makeMesh(x,y,z,sx,sy,false,false);
    PointGroup pg = makePointGroup(x,y,z);
    TriangleGroup tg = makeTriangleGroup(mesh,sx,sy,sz);
    World world = new World();
    world.addChild(pg);
    world.addChild(tg);
    //System.out.println("bs = "+world.getBoundingSphere(true));
    SimpleFrame frame = new SimpleFrame(world);
    OrbitView view = frame.getOrbitView();
    double xs = 0.5*(XMIN+XMAX);
    double ys = 0.5*(XMIN+XMAX);
    double zs = 0.7;
    double rs = 0.25;
    view.setWorldSphere(new BoundingSphere(xs,ys,zs,rs));
    view.setScale(2.3f);
    view.setAzimuth(-90.0f);
    view.setElevation(18.0f);
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
    float size = 0.003f;
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
    TriMesh mesh, Sampling sx, Sampling sy, float[][] sz) 
  {
    sz = transpose(sz);
    TriangleGroup tg = new TriangleGroup(true,makeVertices(mesh,sx,sy,sz));
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
  private static float[] makeVertices(
    TriMesh mesh, Sampling sx, Sampling sy, float[][] z) 
  {
    int nx = sx.getCount()-1;
    int ny = sy.getCount()-1;
    float[] xyz = new float[3*6*nx*ny];
    int ntri = 0;
    for (int ix=0,i=0; ix<nx; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        if (inBounds(x0,XMIN,XMAX) && inBounds(x1,XMIN,XMAX) &&
            inBounds(y0,YMIN,YMAX) && inBounds(y1,YMIN,YMAX)) {
          xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
          xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
          xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
          ++ntri;
          xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
          xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
          xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
          ++ntri;
        }
        /*
        if (inHull(x0,y0,mesh) && inHull(x0,y1,mesh) && inHull(x1,y0,mesh)) {
          xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
          xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
          xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
          ++ntri;
        }
        if (inHull(x1,y0,mesh) && inHull(x0,y1,mesh) && inHull(x1,y1,mesh)) {
          xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
          xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
          xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
          ++ntri;
        }
        */
      }
    }
    return copy(3*3*ntri,xyz);
  }
  private static boolean inHull(float x, float y, TriMesh mesh) {
    TriMesh.PointLocation pl = mesh.locatePoint(x,y);
    return pl.isInside();
  }
  private static boolean inBounds(float x, float xmin, float xmax) {
    return xmin<=x && x<=xmax;
  }

  private static void interpolate() {
    //float[][] data = windowData(XMIN,XMAX,YMIN,YMAX,makeData());
    float[][] data = makeData();
    float[] x = data[0], y = data[1], z = data[2];
    Sampling[] s = makeSamplings(x,y);
    Sampling sx = s[0], sy = s[1];
    plotMesh(x,y,z,sx,sy,true,false,false,false,"Scattered data","sd");
    plotMesh(x,y,z,sx,sy,true,true,false,false,"Delaunay triangles","dt");
    plotMesh(x,y,z,sx,sy,true,true,false,true,"Delaunay triangles","dt");
    plotMesh(x,y,z,sx,sy,true,false,true,false,"Voronoi polygons","vp");
    plotMesh(x,y,z,sx,sy,true,false,true,true,"Voronoi polygons","vp");
    float[][] z1 = interpolate(1,x,y,z,sx,sy);
    float[][] z2 = interpolate(2,x,y,z,sx,sy);
    plot3d(x,y,z,sx,sy,z1);
    plot3d(x,y,z,sx,sy,z2);
  }

  private static void showRoundingError() {
    double a00,f00,af00,a10,f10,af10,a20,f20,af20;
    double a01,f01,af01,a11,f11,af11,a21,f21,af21;
    double a02,f02,af02,a12,f12,af12,a22,f22,af22;
    a00=0.0013129787198019422;  f00=73.27913;  af00=0.09621393701980535;
    a10=4.566615736033324E12;   f10=73.50381;  af10=3.356636420144976E14;
    a20=6.104319161389225E11;   f20=74.20317;  af20=4.5295983720601516E13;
    a01=0.0012717433624584266;  f01=74.193085; af01=0.09435456302890216;
    a11=-6.104319161389172E11;  f11=74.20317;  af11=-4.529598372060113E13;
    a21=-4.566615736033315E12;  f21=73.50381;  af21=-3.3566364201449694E14;
    a02=5.650321386785432E-4;   f02=73.3907;   af02=0.041468104911236044;
    a12=-0.009599257234373463;  f12=73.50381;  af12=-0.7055819517502528;
    a22=-0.0011378416603291027; f22=73.27913;  af22=-0.08338004584105453;

    // Add/subtract in Sambridge's order.
    double asBad =  ( (a00+ a10+ a20)+ (a01+ a11+ a21))+ (a02+ a12+ a22);
    double afsBad = ((af00+af10+af20)+(af01+af11+af21))+(af02+af12+af22);
    System.out.println("asBad  = "+asBad);
    System.out.println("afsBad = "+afsBad);
    System.out.println("fBad   = "+afsBad/asBad);

    // Add/subtract in order of increasing magnitude.
    double asLow =   a02+ a22+ a01+ a00+ a12+ a11+ a20+ a21+ a10;
    double afsLow = af02+af22+af01+af00+af12+af11+af20+af21+af10;
    System.out.println("asLow  = "+asLow);
    System.out.println("afsLow = "+afsLow);
    System.out.println("fLow   = "+afsLow/asLow);
  }

  public static final int PLOT_HEIGHT = 797;
  public static final int PLOT_WIDTH = 765;

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        interpolate();
        //showRoundingError();
      }
    });
  }
}
