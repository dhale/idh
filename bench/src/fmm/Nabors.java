/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.util.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

/**
 * Some ideas in natural-neighbor interpolation.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.10.11
 */
public class Nabors {

  public static void main(String[] args) {
    int n1 = 101;
    int n2 = 101;
    ArrayList<Node> nodeList = new ArrayList<Node>();
    nodeList.add(new Node(   0,   0,0.0f));
    nodeList.add(new Node(n1-1,   0,0.0f));
    nodeList.add(new Node(   0,n2-1,0.0f));
    nodeList.add(new Node(n1-1,n2-1,0.0f));
    nodeList.add(new Node(n1/2,n2/2,1.0f));
    float[][] f = interpolate(n1,n2,nodeList);
    plot(f);
    DistanceMap dmap = new DistanceMap(n1,n2,nodeList);
    //plot(dmap.dists);
    //plot(interpolateNearest(dmap));
    Random random = new Random();
    //int nran = (int)Math.sqrt(n1*n2);
    //int nran = 128;
    int nran = n1*n2/64;
    for (int iran=0; iran<nran; ++iran) {
      Node node = interpolateAndUpdate(dmap);
    }
    float[][] g = interpolate(dmap);
    plot(g);
    //plot(dmap.dists);
  }

  private static float[][] interpolateNearest(DistanceMap dmap) {
    int n1 = dmap.n1;
    int n2 = dmap.n2;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[i2][i1] = dmap.nodes[i2][i1].f;
      }
    }
    return f;
  }

  private static void plot(float[][] f) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setSize(650,600);
    PixelsView pv = sp.addPixels(f);
    pv.setColorModel(ColorMap.JET);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    int n = f.length;
    float[] g = new float[n];
    for (int i=0; i<n; ++i)
      g[i] = f[i][i];
    SimplePlot.asPoints(g);
  }

  private static float[][] interpolate(int n1, int n2, List<Node> nodeList) {
    DistanceMap dmap = new DistanceMap(n1,n2,nodeList);
    return interpolate(dmap);
  }

  private static float[][] interpolate(DistanceMap dmap) {
    int n1 = dmap.n1;
    int n2 = dmap.n2;
    float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[i2][i1] = interpolate(dmap,i1,i2);
      }
    }
    return f;
  }

  private static Node interpolateAndUpdate(DistanceMap dmap) {
    int n1 = dmap.n1;
    int n2 = dmap.n2;
    Node[][] nodes = dmap.nodes;
    float[][] dists = dmap.dists;

    // find indices of point farthest from any existing node
    int j1 = -1;
    int j2 = -1;
    float dj = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (dists[i2][i1]>dj) {
          j1 = i1;
          j2 = i2;
          dj = dists[i2][i1];
        }
      }
    }

    // make a new node at that point
    Node node = new Node(j1,j2,0.0f);

    // update distance map while accumulating node value
    float fsum = nodes[j2][j1].f;
    int nsum = 1;
    dists[j2][j1] = 0.0f;
    nodes[j2][j1] = node;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float dist = distanceSquared(node,i1,i2);
        if (dist<dists[i2][i1]) {
          fsum += nodes[i2][i1].f;
          nsum += 1;
          dists[i2][i1] = dist;
          nodes[i2][i1] = node;
        }
      }
    }
    node.f = fsum/(float)nsum;
    return node;
  }

  private static float interpolate(DistanceMap dmap, int i1, int i2) {
    int n1 = dmap.n1;
    int n2 = dmap.n2;
    Node[][] nodes = dmap.nodes;
    float[][] dists = dmap.dists;
    float fsum = nodes[i2][i1].f;
    int nsum = 1;
    for (int j2=0; j2<n2; ++j2) {
      for (int j1=0; j1<n1; ++j1) {
        float dij = distanceSquared(i1,i2,j1,j2);
        if (dij<dists[j2][j1]) {
          fsum += nodes[j2][j1].f;
          nsum += 1;
        }
      }
    }
    return fsum/(float)nsum;
  }

  private static Node getNearestNode(List<Node> nodeList, int i1, int i2) {
    float dmin = Float.MAX_VALUE;
    Node nmin = null;
    for (Node node:nodeList) {
      float d = distanceSquared(node,i1,i2);
      if (d<dmin) {
        dmin = d;
        nmin = node;
      }
    }
    return nmin;
  }

  private static float distanceSquared(Node node, int i1, int i2) {
    return distanceSquared(i1,i2,node.i1,node.i2);
  }

  private static float distanceSquared(int i1, int i2, int j1, int j2) {
    float d1 = (float)(i1-j1);
    float d2 = (float)(i2-j2);
    return d1*d1+d2*d2;
  }

  private static class DistanceMap {
    int n1,n2;
    Node[][] nodes;
    float[][] dists;
    DistanceMap(int n1, int n2, List<Node> nodeList) {
      this.n1 = n1;
      this.n2 = n2;
      nodes = new Node[n2][n1];
      dists = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          nodes[i2][i1] = getNearestNode(nodeList,i1,i2);
          dists[i2][i1] = distanceSquared(nodes[i2][i1],i1,i2);
        }
      }
    }
  }

  private static class Node {
    int i1,i2;
    float f;
    Node(int i1, int i2, float f) {
      this.i1 = i1; 
      this.i2 = i2; 
      this.f = f;
    }
  }
}
