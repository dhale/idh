/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fault;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Tracks faults using thinned fault likelihoods and orientations.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.10.25
 */
public class FaultTracker3 {

  /**
   * Constructs a fault tracker for specified likelihoods and orientations.
   * @param fpt array {f,p,t} of thinned fault likelihoods, strikes and dips.
   */
  public FaultTracker3(float[][][][] fpt) {
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

  /**
   * Returns arrays {i1,i2,i3} of coordinates of points on a fault.
   * @param i1 index of seed in 1st dimension.
   * @param i2 index of seed in 2nd dimension.
   * @param i3 index of seed in 3rd dimension.
   * @return arrays {i1,i2,i3} of coordinate of points
   */
  public int[][] track(int i1, int i2, int i3) {
    int n1 = _n1, n2 = _n2, n3 = _n3;

    // Lists for samples in the fault.
    IntList i1List = new IntList();
    IntList i2List = new IntList();
    IntList i3List = new IntList();

    // Stacks for candidate samples that have been found.
    // Candidates will eventually be included in the fault, 
    // but their neighbors may have not yet been processed.
    // The found flags ensure that no sample is processed
    // (added to the stacks) more than once.
    IntList i1Stack = new IntList();
    IntList i2Stack = new IntList();
    IntList i3Stack = new IntList();
    boolean[][][] found = new boolean[n3][n2][n1];

    // First candidate is the specified seed.
    i1Stack.push(i1);
    i2Stack.push(i2);
    i3Stack.push(i3);
    found[i3][i2][i1] = true;

    // While candidate samples exist, ...
    while (!i1Stack.isEmpty()) {

      // Get candidate sample i from the stack.
      i1 = i1Stack.pop();
      i2 = i2Stack.pop();
      i3 = i3Stack.pop();

      // If candidate sample i is in the fault, ...
      if (sampleInFault(i1,i2,i3)) {

        // Add the candiate sample i to the fault list.
        i1List.add(i1);
        i2List.add(i2);
        i3List.add(i3);

        // Normal and strike vectors for sample i.
        float phi = toRadians(_p[i3][i2][i1]);
        float theta = toRadians(_t[i3][i2][i1]);
        float cp = cos(phi); 
        float sp = sin(phi);
        float ct = cos(theta);
        float st = sin(theta);
        float u1i = -st;
        float u2i = -sp*ct;
        float u3i =  cp*ct;

        // For all neighbor samples, ...
        for (int k3=0,j3=i3-1; k3<3; ++k3,++j3) {
          for (int k2=0,j2=i2-1; k2<3; ++k2,++j2) {
            for (int k1=0,j1=i1-1; k1<3; ++k1,++j1) {

              // Skip neighbors our of bounds or already in stack.
              if (j1<0 || j1>=n1 ||
                  j2<0 || j2>=n2 ||
                  j3<0 || j3>=n3 ||
                  found[j3][j2][j1])
                continue;

              // If neighbor could be in the fault, push it onto the stack.
              if (naborMayBeInFault(u1i,u2i,u3i,i1,i2,i3,k1,k2,k3)) {
                i1Stack.push(j1);
                i2Stack.push(j2);
                i3Stack.push(j3);
                found[j3][j2][j1] = true;
              }
            }
          }
        }
      }
    }
    int[] i1s = i1List.trim();
    int[] i2s = i2List.trim();
    int[] i3s = i3List.trim();
    return new int[][]{i1s,i2s,i3s};
  }
  public static float[] xyz(int[][] i123) {
    int[] i1 = i123[0];
    int[] i2 = i123[1];
    int[] i3 = i123[2];
    int n = i1.length;
    float[] xyz = new float[3*n];
    for (int i=0,j=0; i<n; ++i) {
      xyz[j++] = i3[i];
      xyz[j++] = i2[i];
      xyz[j++] = i1[i];
    }
    return xyz;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _f,_p,_t;
  private int _n1,_n2,_n3;
  private float _fmin = 0.1f;

  private static final float HUGE = Float.MAX_VALUE;

  private static float[][][] VK1 = new float[3][3][3];
  private static float[][][] VK2 = new float[3][3][3];
  private static float[][][] VK3 = new float[3][3][3];
  static {
    for (int k3=0; k3<3; ++k3) {
      float v3 = k3-1.0f;
      for (int k2=0; k2<3; ++k2) {
        float v2 = k2-1.0f;
        for (int k1=0; k1<3; ++k1) {
          float v1 = k1-1.0f;
          VK1[k3][k2][k1] = v1;
          VK2[k3][k2][k1] = v2;
          VK3[k3][k2][k1] = v3;
        }
      }
    }
  }

  private class IntList {
    public int n;
    public int[] a = new int[1024];
    public void add(int i) {
      if (n==a.length) {
        int[] t = new int[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = i;
    }
    public int[] trim() {
      if (n==0)
        return null;
      int[] t = new int[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
    public void push(int i) {
      add(i);
    }
    public int pop() {
      return a[--n];
    }
    public boolean isEmpty() {
      return n==0;
    }
  }

  private float distanceSquaredToPlane(
    float u1, float u2, float u3,
    int k1, int k2, int k3)
  {
    float d = u1*VK1[k3][k2][k1] +
              u2*VK2[k3][k2][k1] +
              u3*VK3[k3][k2][k1];
    return d*d;
  }

  private boolean naborMayBeInFault(
    float u1, float u2, float u3,
    int i1, int i2, int i3,
    int k1, int k2, int k3)
  {
    int j1 = i1+k1-1; 
    int j2 = i2+k2-1; 
    int j3 = i3+k3-1; 
    if (j1==i1 && j2==i2 && j3==i3)
      return false;
    if (j1<0 || j1==_n1 || 
        j2<0 || j2==_n2 || 
        j3<0 || j3==_n3)
      return false;
    float fj = _f[j3][j2][j1];
    if (fj<=_fmin)
      return false;
    float ds = distanceSquaredToPlane(u1,u2,u3,k1,k2,k3);
    if (ds>0.7f)
      return false;
    float[] uj = vectorU(j1,j2,j3);
    float u1j = uj[0], u2j = uj[1], u3j = uj[2];
    float uu = u1*u1j+u2*u2j+u3*u3j;
    float uus = uu*uu;
    if (uus<0.97)
      return false;
    return true;
  }

  private boolean sampleInFault(int i1, int i2, int i3) {
    float[] ui = vectorU(i1,i2,i3);
    float u1i = ui[0], u2i = ui[1], u3i = ui[2];
    float dsi = 0.0f;
    float wsi = 0.0f;
    for (int k3=0,j3=i3-1; k3<3; ++k3,++j3) {
      if (j3<0 || j3>=_n3) continue;
      for (int k2=0,j2=i2-1; k2<3; ++k2,++j2) {
        if (j2<0 || j2>=_n2) continue;
        for (int k1=0,j1=i1-1; k1<3; ++k1,++j1) {
          if (j1<0 || j1>=_n1) continue;
          float[] uj = vectorU(j1,j2,j3);
          float u1j = uj[0], u2j = uj[1], u3j = uj[2];
          float uu = u1i*u1j+u2i*u2j+u3i*u3j;
          float uus = uu*uu;
          float wj = uus*_f[j3][j2][j1];
          dsi += wj*distanceSquaredToPlane(u1i,u2i,u3i,k1,k2,k3);
          wsi += wj;
        }
      }
    }
    float e = sqrt(dsi/wsi);
    return e<0.7f;
  }

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
