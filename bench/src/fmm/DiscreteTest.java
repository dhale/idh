/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import java.util.Random;
import javax.swing.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.SibsonInterpolator2;
import edu.mines.jtk.mesh.TriMesh;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.sgl.test.TestFrame;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.11.20
 */
public class DiscreteTest {

  ///////////////////////////////////////////////////////////////////////////
  // 1D

  private static float[][] getDistanceNearest(float[] x, float[] f, int n1) {
    int nx = x.length;
    float[] d = new float[n1];
    float[] p = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      float xi = (float)i1;
      int ki = binarySearch(x,xi);
      if (ki<0) ki = -1-ki;
      if (ki>=nx) ki = nx-1;
      float di = abs(x[ki]-xi);
      if (ki>0) {
        float dm = abs(x[ki-1]-xi);
        if (dm<di) {
          di = dm;
          ki -= 1;
        }
      }
      d[i1] = di;
      p[i1] = f[ki];
    }
    return new float[][]{d,p};
  }

  private static void smoothGaussian(
    float c, float[] s, float[] x, float[] y) 
  {
    int n1 = x.length;
    float csmax = c* max(s);
    int niter = 1+4*(int)csmax;
    trace("niter="+niter);
    float[] z = copy(x);
    float[] w = copy(x);
    c *= -0.5f/(float)niter;
    for (int jiter=0; jiter<niter; ++jiter) {
      float gi;
      for (int i1=1; i1<n1; ++i1) {
        gi  = z[i1  ];
        gi -= z[i1-1];
        gi *= c*(s[i1]+s[i1-1]);
        w[i1-1] -= gi;
        w[i1  ] += gi;
      }
      float[] t = w; w = z; z = t;
      copy(z,w);
    }
    copy(z,y);
  }

  private static float[] interpolateGaussian(float[] d, float[] p) {
    int n1 = d.length;
    float[] q = new float[n1];
    float[] s = pow(d,2.0f);
    smoothGaussian(0.5f,s,p,q);
    return q;
  }

  private static float[] interpolateSmooth(float[] d, float[] p) {
    int n1 = d.length;
    float[] q = new float[n1];
    //float[] s = add(0.0f,pow(d,2.0f));
    float[] s = add(64.0f,pow(d,2.0f));
    SimplePlot.asPoints(s);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(0.5f,s,p,q);
    return q;
  }

  private static void testInterpolation1() {
    int n1 = 315;
    float[] x = { 60.0f, 100.0f, 170.0f, 200.0f, 250.0f};
    float[] f = {  1.0f,   2.0f,   2.7f,   3.0f,   2.0f};
    //float[] x = { 10.0f, 100.0f, 170.0f, 200.0f, 250.0f};
    //float[] f = {  1.0f,   2.0f,   1.2f,   2.3f,   1.5f};
    float[][] dp = getDistanceNearest(x,f,n1);
    float[] d = dp[0];
    float[] p = dp[1];
    //d = fillfloat(10.0f,n1);
    //p = zerofloat(n1);
    for (int i=0; i<x.length; ++i)
      p[(int)x[i]] = f[i];
    //float[] q = interpolateGaussian(d,p);
    float[] q = interpolateSmooth(d,p);
    SimplePlot sp = new SimplePlot();
    PointsView pvf = sp.addPoints(x,f);
    pvf.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE);
    pvf.setLineStyle(PointsView.Line.NONE);
    PointsView pvp = sp.addPoints(p);
    pvp.setLineStyle(PointsView.Line.DASH);
    PointsView pvq = sp.addPoints(q);
  }

  ///////////////////////////////////////////////////////////////////////////
  // 2D

  private static float[][] makeSamples(int mode, int n, int n1, int n2) {
    float[] f = new float[n];
    float[] g1 = new float[n];
    float[] g2 = new float[n];
    float[] x1 = new float[n];
    float[] x2 = new float[n];
    float scale = 1.0f/(float)(n1+n2);
    double radius = 0.25*min(n1,n2);
    double dtheta = 4.0*PI/(n-1);
    Random r = new Random(314159);
    for (int i=0; i<n; ++i) {
      int i1,i2;
      if (mode==4) {
        if (i==0) {
          i1 = n1/2;
          i2 = n2/2;
        } else {
          double y1,y2,theta;
          if (i<=n/2) {
            theta = i*dtheta;
            y1 = radius*cos(theta);
            y2 = radius*sin(theta);
          } else {
            theta = (i-n/2)*dtheta;
            y1 = 1.5*radius*cos(theta);
            y2 = 1.5*radius*sin(theta);
          }
          i1 = (int)(n1/2+y1+0.5);
          i2 = (int)(n2/2+y2+0.5);
        }
      } else if (mode==3 && i==0) {
        i1 = n1/2;
        i2 = n2/2;
      } else {
        i1 = r.nextInt(n1);
        i2 = r.nextInt(n2);
      }
      x1[i] = i1;
      x2[i] = i2;
      if (mode==0) {
        f[i] = (float)(0.5+0.25*sin(PI*i1/n1)*sin(PI*i2/n2));
        g1[i] = (float)(0.25*PI/n1*cos(PI*i1/n1)*sin(PI*i2/n2));
        g2[i] = (float)(0.25*PI/n2*sin(PI*i1/n1)*cos(PI*i2/n2));
      } else if (mode==1) {
        f[i] = (float)(i1+i2)*scale;
        g1[i] = scale;
        g2[i] = scale;
      } else if (mode==2) {
        f[i] = 0.25f+0.50f*r.nextFloat();
      } else if (mode==3 || mode==4) {
        f[i] = (i==0)?0.9f:0.5f;
        g1[i] = 0.0f;
        g2[i] = 0.0f;
      }
    }
    return new float[][]{f,g1,g2,x1,x2};
  }

  private static float[][][] getDistanceNearest(
    float[] f, float[] g1, float[] g2, float[] x1, float[] x2, int n1, int n2) 
  {
    int nx = x1.length;
    float[][] d = new float[n2][n1];
    float[][] p = new float[n2][n1];
    TriMesh tm = new TriMesh();
    for (int ix=0; ix<nx; ++ix) {
      TriMesh.Node node = new TriMesh.Node(x1[ix],x2[ix]);
      node.index = ix;
      tm.addNode(node);
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float y1i = (float)i1;
        float y2i = (float)i2;
        TriMesh.Node node = tm.findNodeNearest(y1i,y2i);
        int index = node.index;
        float fi = f[index];
        float g1i = (g1!=null)?g1[index]:0.0f;
        float g2i = (g2!=null)?g2[index]:0.0f;
        float x1i = x1[index];
        float x2i = x2[index];
        float d1i = y1i-x1i;
        float d2i = y2i-x2i;
        d[i2][i1] = sqrt(d1i*d1i+d2i*d2i);
        p[i2][i1] = fi+g1i*d1i+g2i*d2i;
      }
    }
    return new float[][][]{d,p};
  }

  private static float[][] interpolateLaplace(float[][] d, float[][] p) {
    int n1 = d[0].length;
    int n2 = d.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[][] q = copy(p);
    float dmax = max(d);
    int niter = 1+(int)(dmax*dmax);
    trace("interpolateLaplace: niter="+niter);
    for (int jiter=0; jiter<niter; ++jiter) {
      for (int i2=0; i2<n2; ++i2) {
        int i2m = (i2==0  )?i2:i2-1;
        int i2p = (i2==n2m)?i2:i2+1;
        for (int i1=0; i1<n1; ++i1) {
          int i1m = (i1==0  )?i1:i1-1;
          int i1p = (i1==n1m)?i1:i1+1;
          if (d[i2][i1]!=0.0f) {
            float q0m = q[i2 ][i1m];
            float q0p = q[i2 ][i1p];
            float qm0 = q[i2m][i1 ];
            float qp0 = q[i2p][i1 ];
            float qav = 0.25f*(q0m+q0p+qm0+qp0);
            q[i2][i1] = qav;
          }
        }
      }
    }
    return q;
  }

  private static float[][] interpolateBiLaplace(float[][] d, float[][] p) {
    int n1 = d[0].length;
    int n2 = d.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float a =  8.0f/20.0f;
    float b = -2.0f/20.0f;
    float c = -1.0f/20.0f;
    float[][] q = copy(p);
    float dmax = max(d);
    float pmax = max(p);
    float pmin = min(p);
    float dpeps = 0.000001f*(pmax-pmin);
    float dseps = dpeps*dpeps;
    float dsmax = Float.MAX_VALUE;
    int maxiter = 1+(int)(dmax*dmax*dmax);
    trace("interpolateBiLaplace: maxiter="+maxiter+" dseps="+dseps);
    for (int niter=0; niter<maxiter && dsmax>dseps; ++niter) {
      if (niter%100==0)
        trace("interpolateBiLaplace: niter="+niter+" dsmax="+dsmax);
      dsmax = 0.0f;
      for (int i2=0; i2<n2; ++i2) {
        int i2m = (i2==0  )?i2:i2-1;
        int i2p = (i2==n2m)?i2:i2+1;
        int i2mm = (i2m==0  )?i2m:i2m-1;
        int i2pp = (i2p==n2m)?i2p:i2p+1;
        for (int i1=0; i1<n1; ++i1) {
          int i1m = (i1==0  )?i1:i1-1;
          int i1p = (i1==n1m)?i1:i1+1;
          int i1mm = (i1m==0  )?i1m:i1m-1;
          int i1pp = (i1p==n1m)?i1p:i1p+1;
          if (d[i2][i1]!=0.0f) {
            float aq = a*(q[i2 ][i1m]+q[i2 ][i1p]+q[i2m][i1 ]+q[i2p][i1 ]);
            float bq = b*(q[i2m][i1m]+q[i2m][i1p]+q[i2p][i1m]+q[i2p][i1p]);
            float cq = c*(q[i2][i1mm]+q[i2][i1pp]+q[i2mm][i1]+q[i2pp][i1]);
            float qn = aq+bq+cq;
            float qc = q[i2][i1];
            float dq = qn-qc;
            float ds = dq*dq;
            if (ds>dsmax) 
              dsmax = ds;
            q[i2][i1] = qn;
          }
        }
      }
    }
    return q;
  }

  private static float[][] interpolateApproxSibson(float[][] d, float[][] p) {
    int n1 = d[0].length;
    int n2 = d.length;
    PaintingX px = new PaintingX(n1,n2,1);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (d[i2][i1]==0.0f)
          px.paintAt(i1,i2,p[i2][i1]);
      }
    }
    px.extrapolate();
    px.interpolate();
    return px.getValues();
  }

  private static float[][] interpolateDiscrSibson(float[][] d, float[][] p) {
    int n1 = d[0].length;
    int n2 = d.length;
    float[][] q = new float[n2][n1];
    float[][] c = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float pi = p[i2][i1];
        float di = d[i2][i1];
        if (di==0.0f) {
          q[i2][i1] = pi;
          c[i2][i1] = 1.0f;
        } else {
          float ds = di*di;
          int j1lo = max( 0,(int)(i1-di)-1);
          int j1hi = min(n1,(int)(i1+di)+2);
          int j2lo = max( 0,(int)(i2-di)-1);
          int j2hi = min(n2,(int)(i2+di)+2);
          for (int j2=j2lo; j2<j2hi; ++j2) {
            for (int j1=j1lo; j1<j1hi; ++j1) {
              float d1 = (float)(j1-i1);
              float d2 = (float)(j2-i2);
              if (d1*d1+d2*d2<ds) {
                q[j2][j1] += pi;
                c[j2][j1] += 1.0f;
              }
            }
          }
        }
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        q[i2][i1] /= c[i2][i1];
      }
    }
    return q;
  }

  private static float[][] interpolateExactSibson(float[][] d, float[][] p) {
    int n1 = d[0].length;
    int n2 = d.length;
    int n = 0;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        if (d[i2][i1]==0.0f) ++n;
    float[] f = new float[n];
    float[] x1 = new float[n];
    float[] x2 = new float[n];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (d[i2][i1]==0.0f) {
          f[i] = p[i2][i1];
          x1[i] = i1;
          x2[i] = i2;
        }
      }
    }
    SibsonInterpolator2 si = new SibsonInterpolator2(f,x1,x2);
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    si.setBounds(s1,s2);
    return si.interpolate(s1,s2);
  }

  private static float[][] sinsin(int n1, int n2) {
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[i2][i1] = (float)(0.5+0.25*sin(PI*i1/n1)*sin(PI*i2/n2));
      }
    }
    return f;
  }

  private static float[][] smooth(double sigma, float[][] x) {
    float[][] y = copy(x);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply00(x,y);
    return y;
  }

  private static void shift(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    for (int i2=n2-1; i2>0; --i2) {
      for (int i1=n1-1; i1>0; --i1) {
        s[i2][i1] = 0.25f*(s[i2][i1]+s[i2-1][i1]+s[i2][i1-1]+s[i2-1][i1-1]);
      }
    }
  }

  private static float[][] interpolateSmooth(float[][] d, float[][] p) {
    float[][] q = copy(p);
    float[][] s = mul(d,d);
    shift(s);
    float c = 0.500f;
    Tensors2 t2 = new IdentityTensors2();
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.0001,10000);
    lsf.apply(t2,c,s,p,q);
    return q;
  }
  private static class IdentityTensors2 implements Tensors2 {
    public void getTensor(int i1, int i2, float[] d) {
      d[0] = 1.0f;
      d[1] = 0.0f;
      d[2] = 1.0f;
    }
  }

  private static void plot(String title, float[][] f) {
    SimplePlot sp = new SimplePlot();
    sp.setSize(880,775);
    sp.setTitle(title);
    sp.addColorBar();
    PixelsView pv = sp.addPixels(f);
    pv.setColorModel(ColorMap.JET);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  }

  private static void printSamples(float[][] s) {
    float[] x1 = s[3];
    float[] x2 = s[4];
    float[] x3 = s[0];
    int n = x1.length;
    for (int i=0; i<n; ++i) {
      System.out.printf("%10.2f %10.2f %10.2f%n",x1[i],x2[i],x3[i]);
    }
  }

  private static void adjustSamples(int n1, int n2, float[][] s) {
    float[] x1 = s[3];
    float[] x2 = s[4];
    float[] x3 = s[0];
    mul(1.0f/(float)(n1-1),x1,x1);
    mul(1.0f/(float)(n2-1),x2,x2);
    mul(100.0f,x3,x3);
  }

  private static void plotSamples(float[][] s) {
    float[] x1 = s[3];
    float[] x2 = s[4];
    float[] x3 = s[0];
    PlotPanel panel = new PlotPanel(1,1);
    double pad = 0.05;
    panel.setLimits(-pad,-pad,1.0+pad,1.0+pad);
    PointsView pv = panel.addPoints(x1,x2,x3);
    pv.setLineStyle(PointsView.Line.NONE);
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv.setTextFormat("%2.0f");
    PlotFrame frame = new PlotFrame(panel);
    frame.setSize(800,800);
    frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
    frame.setVisible(true);
    frame.paintToPng(300,6,"junk.png");
  }

  private static void plot(float[][] s, float[][] f) {
    float[] x1 = s[3];
    float[] x2 = s[4];
    int n1 = f[0].length;
    int n2 = f.length;
    PlotPanel panel = new PlotPanel(1,1);
    double pad = 0.05;
    panel.setLimits(-pad,-pad,1.0+pad,1.0+pad);
    panel.setColorBarWidthMinimum(60);
    panel.addColorBar();
    Sampling s1 = new Sampling(n1,1.0/(n1-1),0.0);
    Sampling s2 = new Sampling(n2,1.0/(n2-1),0.0);
    PixelsView fv = panel.addPixels(s1,s2,f);
    fv.setColorModel(ColorMap.GRAY);
    PointsView sv = panel.addPoints(x1,x2);
    sv.setLineStyle(PointsView.Line.NONE);
    sv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    sv.setMarkColor(Color.RED);
    PlotFrame frame = new PlotFrame(panel);
    frame.setSize(910,800);
    frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
    frame.setVisible(true);
    //frame.paintToPng(300,6,"junk.png");
  }

  private static PointGroup makePointGroup(
    float[][] s, Sampling s1, Sampling s2) 
  {
    float f1 = (float)s1.getFirst();
    float f2 = (float)s2.getFirst();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    int np = s[0].length;
    float[] xyz = new float[3*np];
    for (int ip=0,i=0; ip<np; ++ip) {
      xyz[i++] =  f2+d2*s[4][ip];
      xyz[i++] =  f1+d1*s[3][ip];
      xyz[i++] = -s[0][ip];
    }
    float size = 2.0f*min(d1,d2);
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
    Sampling s1, Sampling s2, float[][] f) 
  {
    f = neg(f);
    TriangleGroup tg = new TriangleGroup(true,s2,s1,f);
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

  private static void plot3d(
    float[][] s, Sampling s1, Sampling s2, float[][] f) 
  {
    PointGroup pg = makePointGroup(s,s1,s2);
    TriangleGroup tg = makeTriangleGroup(s2,s1,f);
    World world = new World();
    world.addChild(pg);
    world.addChild(tg);
    TestFrame frame = new TestFrame(world);
    OrbitView view = frame.getOrbitView();
    view.setScale(2.0f);
    view.setElevation(30.0f); // good for sinsin points
    view.setAzimuth(-70.0f);
    //view.setElevation(40.714287f); // good for random points
    //view.setAzimuth(-130.72289f);
    //view.setElevation(30.0f); // good for impulse and circle points
    //view.setAzimuth(18.0f);
    view.setWorldSphere(new BoundingSphere(0.5,0.5,-0.5,1.0));
    //BoundingSphere bs = world.getBoundingSphere(true);
    //trace("bs center: "+bs.getCenter());
    //trace("bs radius: "+bs.getRadius());
    frame.setSize(new Dimension(1200,800));
    frame.setVisible(true);
  }

  private static void testInterpolation(int mode) {
    int n1 = 315;
    int n2 = 315;
    int n = (mode==4)?51:50;
    Sampling s1 = new Sampling(n1,1.0/(n1-1),0.0);
    Sampling s2 = new Sampling(n2,1.0/(n2-1),0.0);
    float[][] s = makeSamples(mode,n,n1,n2);
    float[] f = s[0], g1 = s[1], g2 = s[2], x1 = s[3], x2 = s[4];
    g1 = g2 = null;
    float[][][] dp = getDistanceNearest(f,g1,g2,x1,x2,n1,n2);
    float[][] d = dp[0];
    float[][] p = dp[1];
    float[][] qf = (mode==0)?sinsin(n1,n2):zerofloat(n1,n2);
    //float[][] q0 = interpolateApproxSibson(d,p);
    //float[][] q1 = interpolateDiscrSibson(d,p);
    //float[][] q2 = interpolateExactSibson(d,p);
    float[][] q3 = interpolateSmooth(d,p);
    //float[][] q4 = interpolateLaplace(d,p);
    //float[][] q5 = interpolateBiLaplace(d,p);
    //plot3d(s,s1,s2,p);
    //plot3d(s,s1,s2,qf);
    //plot3d(s,s1,s2,q0);
    //plot3d(s,s1,s2,q1);
    //plot3d(s,s1,s2,q2);
    plot3d(s,s1,s2,q3);
    //plot3d(s,s1,s2,q4);
    //plot3d(s,s1,s2,q5);
    /*
    adjustSamples(n1,n2,s);
    plotSamples(s);
    printSamples(s);
    plot(s,p);
    */
  }
  private static void testSinSin() {
    testInterpolation(0);
  }
  private static void testLinear() {
    testInterpolation(1);
  }
  private static void testRandom() {
    testInterpolation(2);
  }
  private static void testImpulse() {
    testInterpolation(3);
  }
  private static void testCircle() {
    testInterpolation(4);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //testInterpolation1();
        //testRandom();
        //testSinSin();
        testLinear();
        //testImpulse();
        //testCircle();
      }
    });
  }
}
