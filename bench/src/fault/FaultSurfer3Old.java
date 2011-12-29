/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import java.util.ArrayList;
import java.util.HashSet;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Finds fault surfaces using fault likelihoods and orientations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.11.23
 */
public class FaultSurfer3Old {

  public static class Node {
    public float fl,fp,ft; // fault likelihood, strike and dip
    public float x1,x2,x3; // point on fault
    public float u1,u2,u3; // normal vector
    public float v1,v2,v3; // strike vector
    public float w1,w2,w3; // dip vector
    Node left,right,above,below; // nabors of this node; may be null
    Node(float x1, float x2, float x3, float fl, float fp, float ft) {
      float p = toRadians(fp);
      float t = toRadians(ft);
      float cp = cos(p); 
      float sp = sin(p);
      float ct = cos(t);
      float st = sin(t);
      this.fl = fl;
      this.fp = fp;
      this.ft = ft;
      this.x1 = x1;
      this.x2 = x2;
      this.x3 = x3;
      this.u1 = -st;
      this.u2 = -sp*ct;
      this.u3 =  cp*ct;
      this.v1 = 0.0f;
      this.v2 = cp;
      this.v3 = sp;
      this.w1 = ct;
      this.w2 = -sp*st;
      this.w3 =  cp*st;
    }
    public String toString() {
      return "("+x1+","+x2+","+x3+"):("+fl+","+fp+","+ft+")";
    }
  }

  public static class Quad {
    Node na,nb,nc,nd; // four nodes in this quad
    Quad left,right,above,below; // four nabors
    Quad(Node na, Node nb, Node nc, Node nd) {
      this.na = na;
      this.nb = nb;
      this.nc = nc;
      this.nd = nd;
    }
    public boolean equals(Quad q) {
      return q.na==na && q.nb==nb || 
             q.nb==na && q.na==nb;
    }
    public int hashCode() {
      return na.hashCode()^nb.hashCode();
    }
  }

  /**
   * Constructs a fault surfer for specified likelihoods and orientations.
   * @param fpt array {f,p,t} of fault likelihoods, strikes and dips.
   */
  public FaultSurfer3Old(float[][][][] fpt) {
    _f = fpt[0];
    _p = fpt[1];
    _t = fpt[2];
    _n1 = _f[0][0].length;
    _n2 = _f[0].length;
    _n3 = _f.length;
  }

  /**
   * Sets a lower bound on fault likelihoods.
   * The default lower bounds is 0.1.
   * @param fmin the lower bound
   */
  public void setThreshold(double fmin) {
    _fmin = (float)fmin;
  }

  public static float[][] slice1(int i1, float[][][] f) {
    int n2 = f[0].length;
    int n3 = f.length;
    float[][] fs = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        fs[i3][i2] = f[i3][i2][i1];
    return fs;
  }

  public static Node[] slice1(int i1, Node[] nodes) {
    int n = nodes.length;
    ArrayList<Node> nodeList = new ArrayList<Node>(n);
    for (int i=0; i<n; ++i) {
      Node node = nodes[i];
      if (node.x1==i1)
        nodeList.add(node);
    }
    return nodeList.toArray(new Node[0]);
  }

  public Node[] findNodes() {
    final int edge = 7;
    int n1 = _n1; int n2 = _n2, n3 = _n3;
    float[][] g2 = new float[n3][n2];
    float[][] g3 = new float[n3][n2];
    float[][] h2 = g2;
    float[][] h3 = g3;
    float[][] h22 = new float[n3][n2];
    float[][] h23 = new float[n3][n2];
    float[][] h33 = new float[n3][n2];
    ArrayList<Node> nodeList = new ArrayList<Node>();
    for (int i1=0; i1<n1; ++i1) {
      float[][] f = slice1(i1,_f);
      float[][] p = slice1(i1,_p);
      float[][] t = slice1(i1,_t);
      _rgf.apply10(f,g2);
      _rgf.apply01(f,g3);
      _rgf.apply20(f,h22);
      _rgf.apply11(f,h23);
      _rgf.apply02(f,h33);
      float[][] a = new float[2][2];
      float[][] z = new float[2][2];
      float[] e = new float[2];
      float etiny = 0.01f;
      for (int i3=0+edge; i3<n3-edge; ++i3) {
        for (int i2=0+edge; i2<n2-edge; ++i2) {
          float h2i = 0.0f;
          float h3i = 0.0f;
          a[0][0] = h22[i3][i2];
          a[0][1] = h23[i3][i2];
          a[1][0] = h23[i3][i2];
          a[1][1] = h33[i3][i2];
          Eigen.solveSymmetric22(a,z,e);
          float eu = e[0];
          float ev = e[1];
          if (ev<0.0f) {
            float v2 = z[1][0]; // normal vector from
            float v3 = z[1][1]; // fault likelihoods
            float fp = toRadians(p[i3][i2]);
            float u2 = -sin(fp); // normal vector from
            float u3 =  cos(fp); // fault strikes
            float uv = u2*v2+u3*v3;
            if (uv*uv>0.5f) { // if normal vectors consistent, ...
              float uscale = 1.0f;
              if (eu-ev<=etiny) {
                uscale = 1.0f-(eu-ev)/etiny;
                uscale = 1.0f-uscale*uscale;
              }
              float ug = uscale*(u2*g2[i3][i2]+u3*g3[i3][i2]);
              h2i = ug*u2;
              h3i = ug*u3;
            }
          }
          h2[i3][i2] = h2i;
          h3[i3][i2] = h3i;
        }
      }
      for (int i3=0+edge; i3<n3-edge; ++i3) {
        int i3p = min(n3-1,i3+1);
        for (int i2=0+edge; i2<n2-edge; ++i2) {
          int i2p = min(n2-1,i2+1);
          if (f[i3][i2]<_fmin && f[i3][i2p]<_fmin && f[i3p][i2]<_fmin)
            continue;
          float h2i,h3i,h2p,h3p,hsi,hsp;
          h2i = h2[i3][i2];
          h3i = h3[i3][i2];
          h2p = h2[i3][i2p];
          h3p = h3[i3][i2p];
          if (h2i*h2p+h3i*h3p<0.0f) {
            hsi = sqrt(h2i*h2i+h3i*h3i);
            hsp = sqrt(h2p*h2p+h3p*h3p);
            float wi = hsp/(hsi+hsp);
            float wp = 1.0f-wi;
            float x2 = wi*i2+wp*i2p;
            float x3 = i3;
            float fl = wi*f[i3][i2]+wp*f[i3][i2p];
            //float fl = (wi>=wp)?f[i3][i2]:f[i3][i2p];
            float fp = (wi>=wp)?p[i3][i2]:p[i3][i2p];
            float ft = (wi>=wp)?t[i3][i2]:t[i3][i2p];
            Node node = new Node(i1,x2,x3,fl,fp,ft);
            nodeList.add(node);
          }
          h2p = h2[i3p][i2];
          h3p = h3[i3p][i2];
          if (h2i*h2p+h3i*h3p<0.0f) {
            hsi = sqrt(h2i*h2i+h3i*h3i);
            hsp = sqrt(h2p*h2p+h3p*h3p);
            float wi = hsp/(hsi+hsp);
            float wp = 1.0f-wi;
            float x2 = i2;
            float x3 = wi*i3+wp*i3p;
            float fl = wi*f[i3][i2]+wp*f[i3p][i2];
            //float fl = (wi>=wp)?f[i3][i2]:f[i3p][i2];
            float fp = (wi>=wp)?p[i3][i2]:p[i3p][i2];
            float ft = (wi>=wp)?t[i3][i2]:t[i3p][i2];
            Node node = new Node(i1,x2,x3,fl,fp,ft);
            nodeList.add(node);
          }
        }
      }
    }
    return nodeList.toArray(new Node[0]);
  }

  public float[][] getXyzUvwRgb(Quad[] quads) {
    int nq = quads.length;
    float[] xyz = new float[3*4*nq];
    float[] uvw = new float[3*4*nq];
    float[] fc = new float[4*nq];
    for (int iq=0,iv=0,in=0,ic=0; iq<nq; ++iq) {
      Quad qi = quads[iq];
      xyz[iv++] = qi.na.x3; 
      xyz[iv++] = qi.na.x2; 
      xyz[iv++] = qi.na.x1;
      xyz[iv++] = qi.nb.x3; 
      xyz[iv++] = qi.nb.x2; 
      xyz[iv++] = qi.nb.x1;
      xyz[iv++] = qi.nc.x3; 
      xyz[iv++] = qi.nc.x2; 
      xyz[iv++] = qi.nc.x1;
      xyz[iv++] = qi.nd.x3; 
      xyz[iv++] = qi.nd.x2; 
      xyz[iv++] = qi.nd.x1;
      uvw[in++] = qi.na.u3;
      uvw[in++] = qi.na.u2;
      uvw[in++] = qi.na.u1;
      uvw[in++] = qi.nb.u3;
      uvw[in++] = qi.nb.u2;
      uvw[in++] = qi.nb.u1;
      uvw[in++] = qi.nc.u3;
      uvw[in++] = qi.nc.u2;
      uvw[in++] = qi.nc.u1;
      uvw[in++] = qi.nd.u3;
      uvw[in++] = qi.nd.u2;
      uvw[in++] = qi.nd.u1;
      fc[ic++] = qi.na.fl;
      fc[ic++] = qi.nb.fl;
      fc[ic++] = qi.nc.fl;
      fc[ic++] = qi.nd.fl;
    }
    float fcmin = min(fc);
    float fcmax = max(fc);
    ColorMap cmap = new ColorMap(fcmin,fcmax,ColorMap.JET);
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz,uvw,rgb};
  }
  
  private static KdTree makeKdTree(Node[] nodes) {
    int nnode = nodes.length;
    float[] x1 = new float[nnode];
    float[] x2 = new float[nnode];
    float[] x3 = new float[nnode];
    for (int inode=0; inode<nnode; ++inode) {
      Node nodei = nodes[inode];
      x1[inode] = nodei.x1;
      x2[inode] = nodei.x2;
      x3[inode] = nodei.x3;
    }
    float[][] x = new float[][]{x1,x2,x3};
    return new KdTree(x);
  }

  public Node[] mergeNodes(Node[] nodes) {
    int nnode = nodes.length;
    KdTree kdt = makeKdTree(nodes);
    ArrayList<Node> nodeList = new ArrayList<Node>(nnode);
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    float[][] a = new float[2][2];
    float[][] v = new float[2][2];
    float[] d = new float[2];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          xmin[0] = i1-0.5f; xmax[0] = i1+0.5f;
          xmin[1] = i2-0.5f; xmax[1] = i2+0.5f;
          xmin[2] = i3-0.5f; xmax[2] = i3+0.5f;
          int[] i = kdt.findInRange(xmin,xmax);
          int nn = i.length;
          Node ni = null;
          if (nn==1) {
            ni = nodes[i[0]];
          } else if (nn>1) {
            float fl = 0.0f, fp = 0.0f, ft = 0.0f;
            float x1 = 0.0f, x2 = 0.0f, x3 = 0.0f;
            float sn = 1.0f/nn;
            float ccp = 0.0f;
            float ssp = 0.0f;
            float csp = 0.0f;
            for (int j=0; j<nn; ++j) {
              Node nj = nodes[i[j]];
              x1 += sn*nj.x1;
              x2 += sn*nj.x2;
              x3 += sn*nj.x3;
              fl += sn*nj.fl;
              ft += sn*nj.ft;
              float cpj = nj.v2;
              float spj = nj.v3;
              ccp += cpj*cpj;
              ssp += spj*spj;
              csp += cpj*spj;
            }
            a[0][0] =  ccp; a[1][0] = -csp;
            a[0][1] = -csp; a[1][1] =  ssp;
            Eigen.solveSymmetric22(a,v,d);
            fp = toDegrees(atan2(v[1][0],v[1][1]));
            if (fp<-90.0f) fp += 180.0f;
            if (fp> 90.0f) fp -= 180.0f;
            ni = new Node(x1,x2,x3,fl,fp,ft);
          }
          if (ni!=null)
            nodeList.add(ni);
        }
      }
    }
    return nodeList.toArray(new Node[0]);
  }

  public Node[] linkNodes(Node[] nodes) {
    int nn = nodes.length;
    KdTree kdt = makeKdTree(nodes);
    float[] xmin = new float[3];
    float[] xmax = new float[3];

    // Linking threshold for dot products of dip and strike vectors.
    final float dlink = 0.5f;

    // Link with nabors left and right.
    for (int in=0; in<nn; ++in) {
      Node ni = nodes[in];
      float x1i = ni.x1;
      float x2i = ni.x2;
      float x3i = ni.x3;
      xmin[0] = x1i-0.5f; xmax[0] = x1i+0.5f;
      xmin[1] = x2i-1.5f; xmax[1] = x2i+1.5f;
      xmin[2] = x3i-1.5f; xmax[2] = x3i+1.5f;
      int[] i = kdt.findInRange(xmin,xmax);
      int jmin = -1;
      int jmax = -1;
      float dmin = 0.0f;
      float dmax = 0.0f;
      float vi2 = ni.v2;
      float vi3 = ni.v3;
      for (int j=0; j<i.length; ++j) {
        Node nj = nodes[i[j]];
        if (nj!=ni && nj.x1==ni.x1) {
          float vj2 = nj.x2-ni.x2;
          float vj3 = nj.x3-ni.x3;
          float d = (vi2*vj2+vi3*vj3)/sqrt(vj2*vj2+vj3*vj3);
          if (d<dmin) { jmin = j; dmin = d; }
          if (d>dmax) { jmax = j; dmax = d; }
        }
      }
      ni.left  = (dmin<-dlink)?nodes[i[jmin]]:null;
      ni.right = (dmax> dlink)?nodes[i[jmax]]:null;
    }

    // Link with best nabors above and below.
    for (int in=0; in<nn; ++in) {
      Node ni = nodes[in];
      float x1i = ni.x1;
      float x2i = ni.x2;
      float x3i = ni.x3;
      xmin[0] = x1i-1.5f; xmax[0] = x1i+1.5f;
      xmin[1] = x2i-1.5f; xmax[1] = x2i+1.5f;
      xmin[2] = x3i-1.5f; xmax[2] = x3i+1.5f;
      int[] i = kdt.findInRange(xmin,xmax);
      int jmin = -1;
      int jmax = -1;
      float dmin = 0.0f;
      float dmax = 0.0f;
      float wi1 = ni.w1;
      float wi2 = ni.w2;
      float wi3 = ni.w3;
      for (int j=0; j<i.length; ++j) {
        Node nj = nodes[i[j]];
        if (nj!=ni && nj.x1!=ni.x1) {
          float wj1 = nj.x1-ni.x1;
          float wj2 = nj.x2-ni.x2;
          float wj3 = nj.x3-ni.x3;
          float d = (wi1*wj1+wi2*wj2+wi3*wj3)/sqrt(wj1*wj1+wj2*wj2+wj3*wj3);
          if (d<dmin) { jmin = j; dmin = d; }
          if (d>dmax) { jmax = j; dmax = d; }
        }
      }
      ni.above = (dmin<-dlink)?nodes[i[jmin]]:null;
      ni.below = (dmax> dlink)?nodes[i[jmax]]:null;
    }

    // Remove links to node nabors that are not mutual.
    // Return only those nodes with at least one nabor.
    ArrayList<Node> nodeList = new ArrayList<Node>(nn);
    for (int in=0; in<nn; ++in) {
      Node ni = nodes[in];
      if (ni.left!=null && ni!=ni.left.right && ni!=ni.left.left)
        ni.left = null;
      if (ni.right!=null && ni!=ni.right.left && ni!=ni.right.right)
        ni.right = null;
      if (ni.above!=null && ni!=ni.above.below)
        ni.above = null;
      if (ni.below!=null && ni!=ni.below.above)
        ni.below = null;
      if (ni.left!=null || ni.right!=null || ni.above!=null || ni.below!=null)
        nodeList.add(ni);
    }
    trace("linkNodes: nn in="+nn+" out="+nodeList.size());
    return nodeList.toArray(new Node[0]);
  }

  public static Quad[] findQuads(Node[] nodes) {
    int nn = nodes.length;
    HashSet<Quad> quadSet = new HashSet<Quad>();
    for (int in=0; in<nn; ++in) {
      Node ni = nodes[in];
      Node na = ni.above;
      Node nb = ni.below;
      Node nl = ni.left;
      Node nr = ni.right;
      if (nr!=null && na!=null) {
        Node nra = nr.above;
        Node nar = na.right;
        Node nal = na.left;
        if (nra!=null && (nra==nar || nra==nal))
          quadSet.add(new Quad(ni,nr,nra,na));
      }
      if (nl!=null && na!=null) {
        Node nla = nl.above;
        Node nar = na.right;
        Node nal = na.left;
        if (nla!=null && (nla==nar || nla==nal))
          quadSet.add(new Quad(ni,nl,nla,na));
      }
    }
    return quadSet.toArray(new Quad[0]);
  }

  public static float[][][] slice1(int i1, Quad[] quads) {
    float x1 = i1;
    int nq = quads.length;
    ArrayList<float[]> x2List = new ArrayList<float[]>();
    ArrayList<float[]> x3List = new ArrayList<float[]>();
    for (int iq=0; iq<nq; ++iq) {
      Quad qi = quads[iq];
      if (qi.na.x1==x1) {
        float x2a = qi.na.x2;
        float x3a = qi.na.x3;
        float x2b = qi.nb.x2;
        float x3b = qi.nb.x3;
        x2List.add(new float[]{x2a,x2b});
        x3List.add(new float[]{x3a,x3b});
      } else if (qi.nc.x1==x1) {
        float x2c = qi.nc.x2;
        float x3c = qi.nc.x3;
        float x2d = qi.nd.x2;
        float x3d = qi.nd.x3;
        x2List.add(new float[]{x2c,x2d});
        x3List.add(new float[]{x3c,x3d});
      }
    }
    float[][] x2 = x2List.toArray(new float[0][]);
    float[][] x3 = x3List.toArray(new float[0][]);
    return new float[][][]{x2,x3};
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static class FloatList {
    public int n;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

  public float[][][] findFaults() {
    final int n1 = _n1, n2 = _n2, n3 = _n3;
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    final float[][][] h2 = g2;
    final float[][][] h3 = g3;
    final float[][][] h22 = new float[n3][n2][n1];
    final float[][][] h23 = new float[n3][n2][n1];
    final float[][][] h33 = new float[n3][n2][n1];
    _rgf.apply010(_f,g2);
    _rgf.apply001(_f,g3);
    _rgf.apply020(_f,h22);
    _rgf.apply011(_f,h23);
    _rgf.apply002(_f,h33);
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[][] a = new float[2][2];
        float[][] z = new float[2][2];
        float[] e = new float[2];
        float etiny = 0.01f;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float h2i = 0.0f;
            float h3i = 0.0f;
            a[0][0] = h22[i3][i2][i1];
            a[0][1] = h23[i3][i2][i1];
            a[1][0] = h23[i3][i2][i1];
            a[1][1] = h33[i3][i2][i1];
            Eigen.solveSymmetric22(a,z,e);
            float eu = e[0];
            float ev = e[1];
            if (ev<0.0f) {
              float v2 = z[1][0];
              float v3 = z[1][1];
              float fp = toRadians(_p[i3][i2][i1]);
              float u2 = -sin(fp);
              float u3 =  cos(fp);
              float uv = u2*v2+u3*v3;
              if (uv*uv>0.5f) {
                float uscale = 1.0f;
                if (eu-ev<=etiny) {
                  uscale = 1.0f-(eu-ev)/etiny;
                  uscale = 1.0f-uscale*uscale;
                }
                float ug = uscale*(u2*g2[i3][i2][i1]+u3*g3[i3][i2][i1]);
                h2i = ug*u2;
                h3i = ug*u3;
              }
            }
            h2[i3][i2][i1] = h2i;
            h3[i3][i2][i1] = h3i;
          }
        }
      }
    });
    final float[][][] c = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        int i3p = i3+1;
        for (int i2=0; i2<n2; ++i2) {
          int i2p = i2+1;
          for (int i1=0; i1<n1; ++i1) {
            int i1p = i1+1;
            if (_f[i3][i2][i1]<_fmin)
              continue;
            float h2i,h3i,h2p,h3p,hsi,hsp;
            h2i = h2[i3][i2][i1];
            h3i = h3[i3][i2][i1];
            if (h2i==0.0f && h3i==0.0f)
              continue;
            if (i1p<n1) {
              h2p = h2[i3][i2][i1p];
              h3p = h3[i3][i2][i1p];
              if (h2i*h2p+h3i*h3p<0.0f) {
                hsi = h2i*h2i+h3i*h3i;
                hsp = h2p*h2p+h3p*h3p;
                if (hsi<=hsp) {
                  c[i3][i2][i1 ] = _f[i3][i2][i1 ];
                } else {
                  c[i3][i2][i1p] = _f[i3][i2][i1p];
                }
              }
            }
            if (i2p<n2) {
              h2p = h2[i3][i2p][i1];
              h3p = h3[i3][i2p][i1];
              if (h2i*h2p+h3i*h3p<0.0f) {
                hsi = h2i*h2i+h3i*h3i;
                hsp = h2p*h2p+h3p*h3p;
                if (hsi<=hsp) {
                  c[i3][i2 ][i1] = _f[i3][i2 ][i1];
                } else {
                  c[i3][i2p][i1] = _f[i3][i2p][i1];
                }
              }
            }
            if (i3p<n3) {
              h2p = h2[i3p][i2][i1];
              h3p = h3[i3p][i2][i1];
              if (h2i*h2p+h3i*h3p<0.0f) {
                hsi = h2i*h2i+h3i*h3i;
                hsp = h2p*h2p+h3p*h3p;
                if (hsi<=hsp) {
                  c[i3 ][i2][i1] = _f[i3 ][i2][i1];
                } else {
                  c[i3p][i2][i1] = _f[i3p][i2][i1];
                }
              }
            }
          }
        }
      }
    });
    return c;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _f,_p,_t;
  private int _n1,_n2,_n3;
  private float _fmin = 0.1f;
  private RecursiveGaussianFilter _rgf = new RecursiveGaussianFilter(1.0);

  private float[] vectorU(int i1, int i2, int i3) {
    float phi = toRadians(_p[i3][i2][i1]);
    float theta = toRadians(_t[i3][i2][i1]);
    float cp = cos(phi); 
    float sp = sin(phi);
    float ct = cos(theta);
    float st = sin(theta);
    float u1i = -st;
    float u2i = -sp*ct;
    float u3i =  cp*ct;
    return new float[]{u1i,u2i,u3i};
  }

  private float[] getNormalVector(Node node, float[] u) {
    float fp = toRadians(node.fp);
    float ft = toRadians(node.ft);
    float cp = cos(fp); 
    float sp = sin(fp);
    float ct = cos(ft);
    float st = sin(ft);
    float u1 = -st;
    float u2 = -sp*ct;
    float u3 =  cp*ct;
    if (u==null) u = new float[3];
    u[0] = u1;
    u[1] = u2;
    u[2] = u3;
    return u;
  }
  private void uvw(int i1, int i2, int i3, float[][] uvw) {
    float phi = toRadians(_p[i3][i2][i1]);
    float theta = toRadians(_t[i3][i2][i1]);
    float cp = cos(phi); 
    float sp = sin(phi);
    float ct = cos(theta);
    float st = sin(theta);
    float u1i = -st;
    float u2i = -sp*ct;
    float u3i =  cp*ct;
    float v1i = 0.0f;
    float v2i = cp;
    float v3i = sp;
    float w1i = ct;
    float w2i = -sp*st;
    float w3i =  cp*st;
    uvw[0][0] = u1i; uvw[0][1] = u2i; uvw[0][2] = u3i;
    uvw[1][0] = v1i; uvw[1][1] = v2i; uvw[1][2] = v3i;
    uvw[2][0] = w1i; uvw[2][1] = w2i; uvw[2][2] = w3i;
  }
}
