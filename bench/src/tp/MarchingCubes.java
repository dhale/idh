/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package tp;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.MathPlus.sqrt;

/**
 * The marching cubes algorithm for contour surfaces of f(x1,x2,x3).
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.08.15
 */
public class MarchingCubes {

  /**
   * A contour surface computed via the marching cubes algorithm.
   * The surface is represented by lists of vertices and triangles.
   */
  public static class Contour {

    /**
     * Vertex array with packed coordinates (x1,x2,x3).
     * The number of vertices equals x.length/3.
     */
    public float[] x;

    /**
     * Optional normal vector array with packed components (u1,u2,u3).
     * If not null, the number of normal vectors equals u.length/3.
     * This number also equals the number of vertices.
     */
    public float[] u;

    /**
     * Triangle array of packed vertex indices (i1,i2,i3). When multiplied by 
     * 3, each index references the first coordinate x1 of a vertex (x1,x2,x3) 
     * stored in the packed vertex array. The number of triangles equals 
     * i.length/3. A vertex may be referenced by more than one triangle.
     */
    public int[] i;
  }

  /**
   * Constructs marching cubes for an image with specified sampling.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param s3 sampling of 3rd dimension.
   * @param f 3D array of image values; referenced (not copied).
   */
  public MarchingCubes(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _f = f;
  }

  /**
   * Sets the sampling for the image.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param s3 sampling of 3rd dimension.
   */
  public void setSampling(Sampling s1, Sampling s2, Sampling s3) {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
  }

  /**
   * Determines whether contours have normal vectors. The default is
   * true, so that normal vectors are computed for every vertex.
   * @param normals true, for normal vectors; false, otherwise.
   */
  public void setNormals(boolean normals) {
    _normals = normals;
  }

  /**
   * Determines whether contours have vertex coordinates x1 and x3 swapped.
   * If true, any corresponding components u1 and u3 of normal vectors are 
   * also swapped. The default is true, so that coordinates are swapped.
   * <p>
   * This swapping is useful when vertex coordinates (x1,x2,x3) correspond 
   * to (z,y,x) instead of (x,y,z) in 3D graphical views of contours.
   * @param swap13 true, to swap x1 and x3; false, otherwise.
   */
  public void setSwap13(boolean swap13) {
    _swap13 = swap13;
  }

  /**
   * Gets a contour surface for a specified image value.
   * @param c the image value to contour.
   * @param withNormals true, for normal vectors; false, otherwise.
   */
  public Contour getContour(float c) {
    IntList tlist = new IntList();
    FloatList xlist = new FloatList();
    FloatList ulist = _normals?new FloatList():null;
    march(_s1.getCount(),_s2.getCount(),_s3.getCount(),
          _s1.getDelta(),_s2.getDelta(),_s3.getDelta(),
          _s1.getFirst(),_s2.getFirst(),_s3.getFirst(),
          _f,c,tlist,xlist,ulist);
    Contour contour = new Contour();
    contour.i = tlist.trim();
    contour.x = xlist.trim();
    contour.u = _normals?ulist.trim():null;
    if (_swap13) {
      float[] x = contour.x;
      float[] u = contour.u;
      for (int i=x.length-3; i>=0; i-=3) {
        float x1 = x[i  ];
        float x3 = x[i+2];
        x[i  ] = x3;
        x[i+2] = x1;
        if (u!=null) {
          float u1 = u[i  ];
          float u3 = u[i+2];
          u[i  ] = u3;
          u[i+2] = u1;
        }
      }
    }
    return contour;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _s1,_s2,_s3;
  private float[][][] _f;
  private boolean _normals = true;
  private boolean _swap13 = true;

  // In the code below, corners and edges for each cube have indices:
  //         .7----6---.6
  //      10  |     11  |
  //    .     7   .     5
  //  3----2----2       |
  //  |      .4-|--4---.5
  //  3    8    1    9
  //  | .       | . 
  //  0----0----1

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
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

  private static class IntList {
    public int n;
    public int[] a = new int[1024];
    public void add(int f) {
      if (n==a.length) {
        int[] t = new int[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public int[] trim() {
      int[] t = new int[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

  private static void computeNormalVector(
    int j1, int j2, int j3, int k1, int k2, int k3,
    int n1, int n2, int n3, double d1, double d2, double d3,
    double dx, float[][][] f, float[] u)
  {
    double u1,u2,u3;
    double v1,v2,v3;
    if (j1==0) {
      u1 = (f[j3][j2][j1+1]-f[j3][j2][j1  ]);
    } else if (j1==n1-1) {
      u1 = (f[j3][j2][j1  ]-f[j3][j2][j1-1]);
    } else {
      u1 = (f[j3][j2][j1+1]-f[j3][j2][j1-1])*0.5;
    }
    if (k1==0) {
      v1 = (f[k3][k2][k1+1]-f[k3][k2][k1  ]);
    } else if (k1==n1-1) {
      v1 = (f[k3][k2][k1  ]-f[k3][k2][k1-1]);
    } else {
      v1 = (f[k3][k2][k1+1]-f[k3][k2][k1-1])*0.5;
    }
    if (j2==0) {
      u2 = (f[j3][j2+1][j1]-f[j3][j2  ][j1]);
    } else if (j2==n2-1) {
      u2 = (f[j3][j2  ][j1]-f[j3][j2-1][j1]);
    } else {
      u2 = (f[j3][j2+1][j1]-f[j3][j2-1][j1])*0.5;
    }
    if (k2==0) {
      v2 = (f[k3][k2+1][k1 ]-f[k3][k2  ][k1]);
    } else if (k2==n2-1) {
      v2 = (f[k3][k2  ][k1 ]-f[k3][k2-1][k1]);
    } else {
      v2 = (f[k3][k2+1][k1]-f[k3][k2-1][k1])*0.5;
    }
    if (j3==0) {
      u3 = (f[j3+1][j2][j1]-f[j3  ][j2][j1]);
    } else if (j3==n3-1) {
      u3 = (f[j3  ][j2][j1]-f[j3-1][j2][j1]);
    } else {
      u3 = (f[j3+1][j2][j1]-f[j3-1][j2][j1])*0.5;
    }
    if (k3==0) {
      v3 = (f[k3+1][k2][k1]-f[k3  ][k2][k1]);
    } else if (k3==n3-1) {
      v3 = (f[k3  ][k2][k1]-f[k3-1][k2][k1]);
    } else {
      v3 = (f[k3+1][k2][k1]-f[k3-1][k2][k1])*0.5;
    }
    u1 = (u1+(v1-u1)*dx)/d1;
    u2 = (u2+(v2-u2)*dx)/d2;
    u3 = (u3+(v3-u3)*dx)/d3;
    double us = 1.0/sqrt(u1*u1+u2*u2+u3*u3);
    u[0] = (float)(u1*us);
    u[1] = (float)(u2*us);
    u[2] = (float)(u3*us);
  }

  private static void march(
    int n1, int n2, int n3,
    double d1, double d2, double d3,
    double f1, double f2, double f3,
    float[][][] f, float c,
    IntList tlist, FloatList xlist, FloatList ulist)
  {
    float[] u = new float[3];

    // Arrays of indices of vertices computed only once and stored in lists.
    // The two arrays ixa and ixb contain indices for one 2D slab of cubes.
    // A non-negative index in either array points to a computed vertex, and
    // a triangle is represented by three such indices. Here we initialize
    // all of the indices to -1, since no vertices have yet been computed.
    int[][][] ixa = new int[n2][n1][3];
    int[][][] ixb = new int[n2][n1][3];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int kk=0; kk<3; ++kk) {
          ixa[i2][i1][kk] = -1;
          ixb[i2][i1][kk] = -1;
        }
      }
    }

    // Numbers of vertices (and normal vectors) and triangle.
    int nx = 0;
    int nt = 0;

    // For all slabs of cubes, ...
    for (int i3=0; i3<n3-1; ++i3) {

      // For all cubes in this slab, ...
      for (int i2=0; i2<n2-1; ++i2) {
        for (int i1=0; i1<n1-1; ++i1) {

          // Eight corner values for this cube.
          float c0 = f[i3  ][i2  ][i1  ];
          float c1 = f[i3  ][i2  ][i1+1];
          float c2 = f[i3  ][i2+1][i1+1];
          float c3 = f[i3  ][i2+1][i1  ];
          float c4 = f[i3+1][i2  ][i1  ];
          float c5 = f[i3+1][i2  ][i1+1];
          float c6 = f[i3+1][i2+1][i1+1];
          float c7 = f[i3+1][i2+1][i1  ];

          // Case index for this cube.
          int caseIndex = 0;
          if (c0>c) caseIndex +=   1;
          if (c1>c) caseIndex +=   2;
          if (c2>c) caseIndex +=   4;
          if (c3>c) caseIndex +=   8;
          if (c4>c) caseIndex +=  16;
          if (c5>c) caseIndex +=  32;
          if (c6>c) caseIndex +=  64;
          if (c7>c) caseIndex += 128;

          // If at least one triangle for this case, ...
          if (caseIndex>0 && caseIndex<255) {

            // Edges intersected by contour.
            int[] edges = _edges[caseIndex];
            int ne = edges.length;

            // For all triangles (triplets of edge intersections), ...
            for (int ie=0; ie<ne; ie+=3) {

              // For each of three triangle vertices, ...
              for (int je=0; je<3; ++je) {

                // Decode edge j->k into sample indices of j and axis to k.
                int edge = edges[ie+je];
                float cj,ck;
                int j1,j2,j3,kk;
                switch(edge) {
                case 0: // 0->1
                  cj = c0;
                  ck = c1;
                  j1 = i1;
                  j2 = i2;
                  j3 = i3;
                  kk = 0;
                  break;
                case 1: // 1->2
                  cj = c1;
                  ck = c2;
                  j1 = i1+1;
                  j2 = i2;
                  j3 = i3;
                  kk = 1;
                  break;
                case 2: // 3->2
                  cj = c3;
                  ck = c2;
                  j1 = i1;
                  j2 = i2+1;
                  j3 = i3;
                  kk = 0;
                  break;
                case 3: // 0->3
                  cj = c0;
                  ck = c3;
                  j1 = i1;
                  j2 = i2;
                  j3 = i3;
                  kk = 1;
                  break;
                case 4: // 4->5
                  cj = c4;
                  ck = c5;
                  j1 = i1;
                  j2 = i2;
                  j3 = i3+1;
                  kk = 0;
                  break;
                case 5: // 5->6
                  cj = c5;
                  ck = c6;
                  j1 = i1+1;
                  j2 = i2;
                  j3 = i3+1;
                  kk = 1;
                  break;
                case 6: // 7->6
                  cj = c7;
                  ck = c6;
                  j1 = i1;
                  j2 = i2+1;
                  j3 = i3+1;
                  kk = 0;
                  break;
                case 7: // 4->7
                  cj = c4;
                  ck = c7;
                  j1 = i1;
                  j2 = i2;
                  j3 = i3+1;
                  kk = 1;
                  break;
                case 8: // 0->4
                  cj = c0;
                  ck = c4;
                  j1 = i1;
                  j2 = i2;
                  j3 = i3;
                  kk = 2;
                  break;
                case 9: // 1->5
                  cj = c1;
                  ck = c5;
                  j1 = i1+1;
                  j2 = i2;
                  j3 = i3;
                  kk = 2;
                  break;
                case 10: // 3->7
                  cj = c3;
                  ck = c7;
                  j1 = i1;
                  j2 = i2+1;
                  j3 = i3;
                  kk = 2;
                  break;
                default: // 2->6
                  cj = c2;
                  ck = c6;
                  j1 = i1+1;
                  j2 = i2+1;
                  j3 = i3;
                  kk = 2;
                }

                // Index of vertex, if already computed; or -1, if not yet.
                int[] ixjj = (j3==i3)?ixa[j2][j1]:ixb[j2][j1];
                int ix = ixjj[kk];

                // If vertex not yet computed, compute and store coordinates,
                // and optionally compute and store normal vector components.
                if (ix<0) {
                  int k1,k2,k3;
                  double x1,x2,x3;
                  float dx = (c-cj)/(ck-cj);
                  switch(kk) {
                  case 0: // edge aligned with axis 1
                    k1 = j1+1;
                    k2 = j2;
                    k3 = j3;
                    x1 = f1+d1*(j1+dx);
                    x2 = f2+d2*(j2   );
                    x3 = f3+d3*(j3   );
                    break;
                  case 1: // edge aligned with axis 2
                    k1 = j1;
                    k2 = j2+1;
                    k3 = j3;
                    x1 = f1+d1*(j1   );
                    x2 = f2+d2*(j2+dx);
                    x3 = f3+d3*(j3   );
                    break;
                  default: // edge aligned with axis 3
                    k1 = j1;
                    k2 = j2;
                    k3 = j3+1;
                    x1 = f1+d1*(j1   );
                    x2 = f2+d2*(j2   );
                    x3 = f3+d3*(j3+dx);
                  }
                  ix = ixjj[kk] = nx;
                  xlist.add((float)x1);
                  xlist.add((float)x2);
                  xlist.add((float)x3);
                  ++nx;
                  if (ulist!=null) {
                    computeNormalVector(j1,j2,j3,k1,k2,k3,
                                        n1,n2,n3,d1,d2,d3,
                                        dx,f,u);
                    ulist.add(u[0]);
                    ulist.add(u[1]);
                    ulist.add(u[2]);
                  }
                }

                // Append index of vertex to triangle list.
                tlist.add(ix);
              }

              // Number of triangles.
              ++nt;
            }
          }
        }
      }

      // Swap the index arrays ixa and ixb, and re-initialize all ixb to -1.
      int[][][] ixt = ixa;
      ixa = ixb;
      ixb = ixt;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int kk=0; kk<3; ++kk) {
            ixb[i2][i1][kk] = -1;
          }
        }
      }
    }
  }

  // Edges intersected. Each group of three indices represents a triangle. 
  // For the eight sample values in each cube, there are 256 cases. However, 
  // most of those 256 cases are complements or rotations of 16 base cases.
  // Comments at end of each line are case number and base-case number.
  // This table was adapted from one in VTK, the Visualization ToolKit.
  private static final int[][] _edges = { 
    {}, // 0 0
    { 0, 3, 8}, // 1 1
    { 0, 9, 1}, // 2 1
    { 1, 3, 8, 9, 1, 8}, // 3 2
    { 1, 11, 2}, // 4 1
    { 0, 3, 8, 1, 11, 2}, // 5 3
    { 9, 11, 2, 0, 9, 2}, // 6 2
    { 2, 3, 8, 2, 8, 11, 11, 8, 9}, // 7 5
    { 3, 2, 10}, // 8 1
    { 0, 2, 10, 8, 0, 10}, // 9 2
    { 1, 0, 9, 2, 10, 3}, // 10 3
    { 1, 2, 10, 1, 10, 9, 9, 10, 8}, // 11 5
    { 3, 1, 11, 10, 3, 11}, // 12 2
    { 0, 1, 11, 0, 11, 8, 8, 11, 10}, // 13 5
    { 3, 0, 9, 3, 9, 10, 10, 9, 11}, // 14 5
    { 9, 11, 8, 11, 10, 8}, // 15 8
    { 4, 8, 7}, // 16 1
    { 4, 0, 3, 7, 4, 3}, // 17 2
    { 0, 9, 1, 8, 7, 4}, // 18 3
    { 4, 9, 1, 4, 1, 7, 7, 1, 3}, // 19 5
    { 1, 11, 2, 8, 7, 4}, // 20 4
    { 3, 7, 4, 3, 4, 0, 1, 11, 2}, // 21 7
    { 9, 11, 2, 9, 2, 0, 8, 7, 4}, // 22 7
    { 2, 9, 11, 2, 7, 9, 2, 3, 7, 7, 4, 9}, // 23 14
    { 8, 7, 4, 3, 2, 10}, // 24 3
    {10, 7, 4, 10, 4, 2, 2, 4, 0}, // 25 5
    { 9, 1, 0, 8, 7, 4, 2, 10, 3}, // 26 6
    { 4, 10, 7, 9, 10, 4, 9, 2, 10, 9, 1, 2}, // 27 9
    { 3, 1, 11, 3, 11, 10, 7, 4, 8}, // 28 7
    { 1, 11, 10, 1, 10, 4, 1, 4, 0, 7, 4, 10}, // 29 11
    { 4, 8, 7, 9, 10, 0, 9, 11, 10, 10, 3, 0}, // 30 12
    { 4, 10, 7, 4, 9, 10, 9, 11, 10}, // 31 5
    { 9, 4, 5}, // 32 1
    { 9, 4, 5, 0, 3, 8}, // 33 3
    { 0, 4, 5, 1, 0, 5}, // 34 2
    { 8, 4, 5, 8, 5, 3, 3, 5, 1}, // 35 5
    { 1, 11, 2, 9, 4, 5}, // 36 3
    { 3, 8, 0, 1, 11, 2, 4, 5, 9}, // 37 6
    { 5, 11, 2, 5, 2, 4, 4, 2, 0}, // 38 5
    { 2, 5, 11, 3, 5, 2, 3, 4, 5, 3, 8, 4}, // 39 9
    { 9, 4, 5, 2, 10, 3}, // 40 4
    { 0, 2, 10, 0, 10, 8, 4, 5, 9}, // 41 7
    { 0, 4, 5, 0, 5, 1, 2, 10, 3}, // 42 7
    { 2, 5, 1, 2, 8, 5, 2, 10, 8, 4, 5, 8}, // 43 11
    {11, 10, 3, 11, 3, 1, 9, 4, 5}, // 44 7
    { 4, 5, 9, 0, 1, 8, 8, 1, 11, 8, 11, 10}, // 45 12
    { 5, 0, 4, 5, 10, 0, 5, 11, 10, 10, 3, 0}, // 46 14
    { 5, 8, 4, 5, 11, 8, 11, 10, 8}, // 47 5
    { 9, 8, 7, 5, 9, 7}, // 48 2
    { 9, 0, 3, 9, 3, 5, 5, 3, 7}, // 49 5
    { 0, 8, 7, 0, 7, 1, 1, 7, 5}, // 50 5
    { 1, 3, 5, 3, 7, 5}, // 51 8
    { 9, 8, 7, 9, 7, 5, 11, 2, 1}, // 52 7
    {11, 2, 1, 9, 0, 5, 5, 0, 3, 5, 3, 7}, // 53 12
    { 8, 2, 0, 8, 5, 2, 8, 7, 5, 11, 2, 5}, // 54 11
    { 2, 5, 11, 2, 3, 5, 3, 7, 5}, // 55 5
    { 7, 5, 9, 7, 9, 8, 3, 2, 10}, // 56 7
    { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 10, 7}, // 57 14
    { 2, 10, 3, 0, 8, 1, 1, 8, 7, 1, 7, 5}, // 58 12
    {10, 1, 2, 10, 7, 1, 7, 5, 1}, // 59 5
    { 9, 8, 5, 8, 7, 5, 11, 3, 1, 11, 10, 3}, // 60 10
    { 5, 0, 7, 5, 9, 0, 7, 0, 10, 1, 11, 0, 10, 0, 11}, // 61 7
    {10, 0, 11, 10, 3, 0, 11, 0, 5, 8, 7, 0, 5, 0, 7}, // 62 7
    {10, 5, 11, 7, 5, 10}, // 63 2
    {11, 5, 6}, // 64 1
    { 0, 3, 8, 5, 6, 11}, // 65 4
    { 9, 1, 0, 5, 6, 11}, // 66 3
    { 1, 3, 8, 1, 8, 9, 5, 6, 11}, // 67 7
    { 1, 5, 6, 2, 1, 6}, // 68 2
    { 1, 5, 6, 1, 6, 2, 3, 8, 0}, // 69 7
    { 9, 5, 6, 9, 6, 0, 0, 6, 2}, // 70 5
    { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2}, // 71 11
    { 2, 10, 3, 11, 5, 6}, // 72 3
    {10, 8, 0, 10, 0, 2, 11, 5, 6}, // 73 7
    { 0, 9, 1, 2, 10, 3, 5, 6, 11}, // 74 6
    { 5, 6, 11, 1, 2, 9, 9, 2, 10, 9, 10, 8}, // 75 12
    { 6, 10, 3, 6, 3, 5, 5, 3, 1}, // 76 5
    { 0, 10, 8, 0, 5, 10, 0, 1, 5, 5, 6, 10}, // 77 14
    { 3, 6, 10, 0, 6, 3, 0, 5, 6, 0, 9, 5}, // 78 9
    { 6, 9, 5, 6, 10, 9, 10, 8, 9}, // 79 5
    { 5, 6, 11, 4, 8, 7}, // 80 3
    { 4, 0, 3, 4, 3, 7, 6, 11, 5}, // 81 7
    { 1, 0, 9, 5, 6, 11, 8, 7, 4}, // 82 6
    {11, 5, 6, 1, 7, 9, 1, 3, 7, 7, 4, 9}, // 83 12
    { 6, 2, 1, 6, 1, 5, 4, 8, 7}, // 84 7
    { 1, 5, 2, 5, 6, 2, 3, 4, 0, 3, 7, 4}, // 85 10
    { 8, 7, 4, 9, 5, 0, 0, 5, 6, 0, 6, 2}, // 86 12
    { 7, 9, 3, 7, 4, 9, 3, 9, 2, 5, 6, 9, 2, 9, 6}, // 87 7
    { 3, 2, 10, 7, 4, 8, 11, 5, 6}, // 88 6
    { 5, 6, 11, 4, 2, 7, 4, 0, 2, 2, 10, 7}, // 89 12
    { 0, 9, 1, 4, 8, 7, 2, 10, 3, 5, 6, 11}, // 90 13
    { 9, 1, 2, 9, 2, 10, 9, 10, 4, 7, 4, 10, 5, 6, 11}, // 91 6
    { 8, 7, 4, 3, 5, 10, 3, 1, 5, 5, 6, 10}, // 92 12
    { 5, 10, 1, 5, 6, 10, 1, 10, 0, 7, 4, 10, 0, 10, 4}, // 93 7
    { 0, 9, 5, 0, 5, 6, 0, 6, 3, 10, 3, 6, 8, 7, 4}, // 94 6
    { 6, 9, 5, 6, 10, 9, 4, 9, 7, 7, 9, 10}, // 95 3
    {11, 9, 4, 6, 11, 4}, // 96 2
    { 4, 6, 11, 4, 11, 9, 0, 3, 8}, ///97 7
    {11, 1, 0, 11, 0, 6, 6, 0, 4}, // 98 5
    { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 11, 1}, // 99 14
    { 1, 9, 4, 1, 4, 2, 2, 4, 6}, // 100 5
    { 3, 8, 0, 1, 9, 2, 2, 9, 4, 2, 4, 6}, // 101 12
    { 0, 4, 2, 4, 6, 2}, // 102 8
    { 8, 2, 3, 8, 4, 2, 4, 6, 2}, // 103 5
    {11, 9, 4, 11, 4, 6, 10, 3, 2}, // 104 7
    { 0, 2, 8, 2, 10, 8, 4, 11, 9, 4, 6, 11}, // 105 10
    { 3, 2, 10, 0, 6, 1, 0, 4, 6, 6, 11, 1}, // 106 12
    { 6, 1, 4, 6, 11, 1, 4, 1, 8, 2, 10, 1, 8, 1, 10}, // 107 7
    { 9, 4, 6, 9, 6, 3, 9, 3, 1, 10, 3, 6}, // 108 11
    { 8, 1, 10, 8, 0, 1, 10, 1, 6, 9, 4, 1, 6, 1, 4}, // 109 7
    { 3, 6, 10, 3, 0, 6, 0, 4, 6}, // 110 5
    { 6, 8, 4, 10, 8, 6}, // 111 2
    { 7, 6, 11, 7, 11, 8, 8, 11, 9}, // 112 5
    { 0, 3, 7, 0, 7, 11, 0, 11, 9, 6, 11, 7}, // 113 11
    {11, 7, 6, 1, 7, 11, 1, 8, 7, 1, 0, 8}, // 114 9
    {11, 7, 6, 11, 1, 7, 1, 3, 7}, // 115 5
    { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6}, // 116 14
    { 2, 9, 6, 2, 1, 9, 6, 9, 7, 0, 3, 9, 7, 9, 3}, // 117 7
    { 7, 0, 8, 7, 6, 0, 6, 2, 0}, // 118 5
    { 7, 2, 3, 6, 2, 7}, // 119 2
    { 2, 10, 3, 11, 8, 6, 11, 9, 8, 8, 7, 6}, // 120 12
    { 2, 7, 0, 2, 10, 7, 0, 7, 9, 6, 11, 7, 9, 7, 11}, // 121 7
    { 1, 0, 8, 1, 8, 7, 1, 7, 11, 6, 11, 7, 2, 10, 3}, // 122 6
    {10, 1, 2, 10, 7, 1, 11, 1, 6, 6, 1, 7}, // 123 3
    { 8, 6, 9, 8, 7, 6, 9, 6, 1, 10, 3, 6, 1, 6, 3}, // 124 7
    { 0, 1, 9, 10, 7, 6}, // 125 4
    { 7, 0, 8, 7, 6, 0, 3, 0, 10, 10, 0, 6}, // 126 3
    { 7, 6, 10}, // 127 1
    { 7, 10, 6}, // 128 1
    { 3, 8, 0, 10, 6, 7}, // 129 3
    { 0, 9, 1, 10, 6, 7}, // 130 4
    { 8, 9, 1, 8, 1, 3, 10, 6, 7}, // 131 7
    {11, 2, 1, 6, 7, 10}, // 132 3
    { 1, 11, 2, 3, 8, 0, 6, 7, 10}, // 133 6
    { 2, 0, 9, 2, 9, 11, 6, 7, 10}, // 134 7
    { 6, 7, 10, 2, 3, 11, 11, 3, 8, 11, 8, 9}, // 135 12
    { 7, 3, 2, 6, 7, 2}, // 136 2
    { 7, 8, 0, 7, 0, 6, 6, 0, 2}, // 137 5
    { 2, 6, 7, 2, 7, 3, 0, 9, 1}, // 138 7
    { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7}, // 139 14
    {11, 6, 7, 11, 7, 1, 1, 7, 3}, // 140 5
    {11, 6, 7, 1, 11, 7, 1, 7, 8, 1, 8, 0}, // 141 9
    { 0, 7, 3, 0, 11, 7, 0, 9, 11, 6, 7, 11}, // 142 11
    { 7, 11, 6, 7, 8, 11, 8, 9, 11}, // 143 5
    { 6, 4, 8, 10, 6, 8}, // 144 2
    { 3, 10, 6, 3, 6, 0, 0, 6, 4}, // 145 5
    { 8, 10, 6, 8, 6, 4, 9, 1, 0}, // 146 7
    { 9, 6, 4, 9, 3, 6, 9, 1, 3, 10, 6, 3}, // 147 11
    { 6, 4, 8, 6, 8, 10, 2, 1, 11}, // 148 7
    { 1, 11, 2, 3, 10, 0, 0, 10, 6, 0, 6, 4}, // 149 12
    { 4, 8, 10, 4, 10, 6, 0, 9, 2, 2, 9, 11}, // 150 10
    {11, 3, 9, 11, 2, 3, 9, 3, 4, 10, 6, 3, 4, 3, 6}, // 151 7
    { 8, 3, 2, 8, 2, 4, 4, 2, 6}, // 152 5
    { 0, 2, 4, 4, 2, 6}, // 153 8
    { 1, 0, 9, 2, 4, 3, 2, 6, 4, 4, 8, 3}, // 154 12
    { 1, 4, 9, 1, 2, 4, 2, 6, 4}, // 155 5
    { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 11}, // 156 14
    {11, 0, 1, 11, 6, 0, 6, 4, 0}, // 157 5
    { 4, 3, 6, 4, 8, 3, 6, 3, 11, 0, 9, 3, 11, 3, 9}, // 158 7
    {11, 4, 9, 6, 4, 11}, // 159 2
    { 4, 5, 9, 7, 10, 6}, // 160 3
    { 0, 3, 8, 4, 5, 9, 10, 6, 7}, // 161 6
    { 5, 1, 0, 5, 0, 4, 7, 10, 6}, // 162 7
    {10, 6, 7, 8, 4, 3, 3, 4, 5, 3, 5, 1}, // 163 12
    { 9, 4, 5, 11, 2, 1, 7, 10, 6}, // 164 6
    { 6, 7, 10, 1, 11, 2, 0, 3, 8, 4, 5, 9}, // 165 13
    { 7, 10, 6, 5, 11, 4, 4, 11, 2, 4, 2, 0}, // 166 12
    { 3, 8, 4, 3, 4, 5, 3, 5, 2, 11, 2, 5, 10, 6, 7}, // 167 6
    { 7, 3, 2, 7, 2, 6, 5, 9, 4}, // 168 7
    { 9, 4, 5, 0, 6, 8, 0, 2, 6, 6, 7, 8}, // 169 12
    { 3, 2, 6, 3, 6, 7, 1, 0, 5, 5, 0, 4}, // 170 10
    { 6, 8, 2, 6, 7, 8, 2, 8, 1, 4, 5, 8, 1, 8, 5}, // 171 7
    { 9, 4, 5, 11, 6, 1, 1, 6, 7, 1, 7, 3}, // 172 12
    { 1, 11, 6, 1, 6, 7, 1, 7, 0, 8, 0, 7, 9, 4, 5}, // 173 6
    { 4, 11, 0, 4, 5, 11, 0, 11, 3, 6, 7, 11, 3, 11, 7}, // 174 7
    { 7, 11, 6, 7, 8, 11, 5, 11, 4, 4, 11, 8}, // 175 3
    { 6, 5, 9, 6, 9, 10, 10, 9, 8}, // 176 5
    { 3, 10, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9}, // 177 9
    { 0, 8, 10, 0, 10, 5, 0, 5, 1, 5, 10, 6}, // 178 14
    { 6, 3, 10, 6, 5, 3, 5, 1, 3}, // 179 5
    { 1, 11, 2, 9, 10, 5, 9, 8, 10, 10, 6, 5}, // 180 12
    { 0, 3, 10, 0, 10, 6, 0, 6, 9, 5, 9, 6, 1, 11, 2}, // 181 6
    {10, 5, 8, 10, 6, 5, 8, 5, 0, 11, 2, 5, 0, 5, 2}, // 182 7
    { 6, 3, 10, 6, 5, 3, 2, 3, 11, 11, 3, 5}, // 183 3
    { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8}, // 184 11
    { 9, 6, 5, 9, 0, 6, 0, 2, 6}, // 185 5
    { 1, 8, 5, 1, 0, 8, 5, 8, 6, 3, 2, 8, 6, 8, 2}, // 186 7
    { 1, 6, 5, 2, 6, 1}, // 187 2
    { 1, 6, 3, 1, 11, 6, 3, 6, 8, 5, 9, 6, 8, 6, 9}, // 188 7
    {11, 0, 1, 11, 6, 0, 9, 0, 5, 5, 0, 6}, // 189 3
    { 0, 8, 3, 5, 11, 6}, // 190 4
    {11, 6, 5}, // 191 1
    {10, 11, 5, 7, 10, 5}, // 192 2
    {10, 11, 5, 10, 5, 7, 8, 0, 3}, // 193 7
    { 5, 7, 10, 5, 10, 11, 1, 0, 9}, // 194 7
    {11, 5, 7, 11, 7, 10, 9, 1, 8, 8, 1, 3}, // 195 10
    {10, 2, 1, 10, 1, 7, 7, 1, 5}, // 196 5
    { 0, 3, 8, 1, 7, 2, 1, 5, 7, 7, 10, 2}, // 197 12
    { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 10}, // 198 14
    { 7, 2, 5, 7, 10, 2, 5, 2, 9, 3, 8, 2, 9, 2, 8}, // 199 7
    { 2, 11, 5, 2, 5, 3, 3, 5, 7}, // 200 5
    { 8, 0, 2, 8, 2, 5, 8, 5, 7, 11, 5, 2}, // 201 11
    { 9, 1, 0, 5, 3, 11, 5, 7, 3, 3, 2, 11}, // 202 12
    { 9, 2, 8, 9, 1, 2, 8, 2, 7, 11, 5, 2, 7, 2, 5}, // 203 7
    { 1, 5, 3, 3, 5, 7}, // 204 8
    { 0, 7, 8, 0, 1, 7, 1, 5, 7}, // 205 5
    { 9, 3, 0, 9, 5, 3, 5, 7, 3}, // 206 5
    { 9, 7, 8, 5, 7, 9}, // 207 2
    { 5, 4, 8, 5, 8, 11, 11, 8, 10}, // 208 5
    { 5, 4, 0, 5, 0, 10, 5, 10, 11, 10, 0, 3}, // 209 14
    { 0, 9, 1, 8, 11, 4, 8, 10, 11, 11, 5, 4}, // 210 12
    {11, 4, 10, 11, 5, 4, 10, 4, 3, 9, 1, 4, 3, 4, 1}, // 211 7
    { 2, 1, 5, 2, 5, 8, 2, 8, 10, 4, 8, 5}, // 212 11
    { 0, 10, 4, 0, 3, 10, 4, 10, 5, 2, 1, 10, 5, 10, 1}, // 213 7
    { 0, 5, 2, 0, 9, 5, 2, 5, 10, 4, 8, 5, 10, 5, 8}, // 214 7
    { 9, 5, 4, 2, 3, 10}, // 215 4
    { 2, 11, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8}, // 216 9
    { 5, 2, 11, 5, 4, 2, 4, 0, 2}, // 217 5
    { 3, 2, 11, 3, 11, 5, 3, 5, 8, 4, 8, 5, 0, 9, 1}, // 218 6
    { 5, 2, 11, 5, 4, 2, 1, 2, 9, 9, 2, 4}, // 219 3
    { 8, 5, 4, 8, 3, 5, 3, 1, 5}, // 220 5
    { 0, 5, 4, 1, 5, 0}, // 221 2
    { 8, 5, 4, 8, 3, 5, 9, 5, 0, 0, 5, 3}, // 222 3
    { 9, 5, 4}, // 223 1
    { 4, 7, 10, 4, 10, 9, 9, 10, 11}, // 224 5
    { 0, 3, 8, 4, 7, 9, 9, 7, 10, 9, 10, 11}, // 225 12
    { 1, 10, 11, 1, 4, 10, 1, 0, 4, 7, 10, 4}, // 226 11
    { 3, 4, 1, 3, 8, 4, 1, 4, 11, 7, 10, 4, 11, 4, 10}, // 227 7
    { 4, 7, 10, 9, 4, 10, 9, 10, 2, 9, 2, 1}, // 228 9
    { 9, 4, 7, 9, 7, 10, 9, 10, 1, 2, 1, 10, 0, 3, 8}, // 229 6
    {10, 4, 7, 10, 2, 4, 2, 0, 4}, // 230 5
    {10, 4, 7, 10, 2, 4, 8, 4, 3, 3, 4, 2}, // 231 3
    { 2, 11, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4}, // 232 14
    { 9, 7, 11, 9, 4, 7, 11, 7, 2, 8, 0, 7, 2, 7, 0}, // 233 7
    { 3, 11, 7, 3, 2, 11, 7, 11, 4, 1, 0, 11, 4, 11, 0}, // 234 7
    { 1, 2, 11, 8, 4, 7}, // 235 4
    { 4, 1, 9, 4, 7, 1, 7, 3, 1}, // 236 5
    { 4, 1, 9, 4, 7, 1, 0, 1, 8, 8, 1, 7}, // 237 3
    { 4, 3, 0, 7, 3, 4}, // 238 2
    { 4, 7, 8}, // 239 1
    { 9, 8, 11, 11, 8, 10}, // 240 8
    { 3, 9, 0, 3, 10, 9, 10, 11, 9}, // 241 5
    { 0, 11, 1, 0, 8, 11, 8, 10, 11}, // 242 5
    { 3, 11, 1, 10, 11, 3}, // 243 2
    { 1, 10, 2, 1, 9, 10, 9, 8, 10}, // 244 5
    { 3, 9, 0, 3, 10, 9, 1, 9, 2, 2, 9, 10}, // 245 3
    { 0, 10, 2, 8, 10, 0}, // 246 2
    { 3, 10, 2}, // 247 1
    { 2, 8, 3, 2, 11, 8, 11, 9, 8}, // 248 5
    { 9, 2, 11, 0, 2, 9}, // 249 2
    { 2, 8, 3, 2, 11, 8, 0, 8, 1, 1, 8, 11}, // 250 3
    { 1, 2, 11}, // 251 1
    { 1, 8, 3, 9, 8, 1}, // 252 2
    { 0, 1, 9}, // 253 1
    { 0, 8, 3}, // 254 1
    {} // 255 0
  };

  private static void trace(String s) {
    System.out.println(s);
  }
}
