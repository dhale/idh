/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package sw;

import java.io.*;
import java.util.InputMismatchException;
import java.util.Scanner;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Reads horizons for sw data.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.09.19
 */
public class HorizonReader {

  private static final int X_FIRST=1325;
  private static final int Y_FIRST=1316;
  private static final int X_DELTA=1;
  private static final int Y_DELTA=2;

  private static final int NX = 367;
  private static final float DX = 0.025f;
  private static final float FX = 0.0f;

  private static final int NY = 623;
  private static final float DY = 0.025f;
  private static final float FY = 0.0f;

  private static final float BOGUS = 4.99999f;

  private static final String FILENAME = "horizons.txt";

  public static void main(String[] args) {
    float[][][] tb = readTopBase();
    float[][] t = tb[0];
    float[][] b = tb[1];
    float tmin = min(t);
    float tmax = max(t);
    float bmin = min(b);
    float bmax = max(b);
    System.out.println("tmin="+tmin+" tmax="+tmax+" bmin="+bmin+" bmax="+bmax);
  }

  public static float[] xyz(float[][] t) {
    int nx = NX;
    int ny = NY;
    float dx = DX;
    float dy = DY;
    float fx = FX;
    float fy = FY;
    int nmax = 3*3*2*(nx-1)*(ny-1);
    float[] xyz = new float[nmax];
    int n = 0;
    for (int ix=1; ix<nx; ++ix) {
      float xm = fx+(float)(ix-1)*dx;
      float xp = fx+(float)(ix  )*dx;
      for (int iy=1; iy<ny; ++iy) {
        float ym = fy+(float)(iy-1)*dy;
        float yp = fy+(float)(iy  )*dy;
        float a = t[ix-1][iy-1];
        float b = t[ix-1][iy  ];
        float c = t[ix  ][iy  ];
        float d = t[ix  ][iy-1];
        if (a!=BOGUS && b!=BOGUS && c!=BOGUS) {
          xyz[n++] = xm;
          xyz[n++] = ym;
          xyz[n++] = a;
          xyz[n++] = xm;
          xyz[n++] = yp;
          xyz[n++] = b;
          xyz[n++] = xp;
          xyz[n++] = yp;
          xyz[n++] = c;
        } else if (a!=BOGUS && b!=BOGUS && d!=BOGUS) {
          xyz[n++] = xm;
          xyz[n++] = ym;
          xyz[n++] = a;
          xyz[n++] = xm;
          xyz[n++] = yp;
          xyz[n++] = b;
          xyz[n++] = xp;
          xyz[n++] = ym;
          xyz[n++] = d;
        }
        if (a!=BOGUS && c!=BOGUS && d!=BOGUS) {
          xyz[n++] = xm;
          xyz[n++] = ym;
          xyz[n++] = a;
          xyz[n++] = xp;
          xyz[n++] = yp;
          xyz[n++] = c;
          xyz[n++] = xp;
          xyz[n++] = ym;
          xyz[n++] = d;
        } else if (b!=BOGUS && c!=BOGUS && d!=BOGUS) {
          xyz[n++] = xm;
          xyz[n++] = yp;
          xyz[n++] = b;
          xyz[n++] = xp;
          xyz[n++] = yp;
          xyz[n++] = c;
          xyz[n++] = xp;
          xyz[n++] = ym;
          xyz[n++] = d;
        }
      }
    }
    return copy(n,xyz);
  }

  public static float[][][] readTopBase() {
    float[][][] tb = new float[2][NX][NY];
    fill(BOGUS,tb);
    Scanner s;
    try {
      s = new Scanner(new BufferedReader(new FileReader(FILENAME)));
      s.nextLine();
    } catch (IOException ioe) {
      throw new RuntimeException("Cannot open file: "+FILENAME);
    }
    try {
      while (s.hasNextLine()) {
        int kx = s.nextInt();
        int ky = s.nextInt();
        float t = s.nextFloat();
        float b = s.nextFloat();
        //System.out.println("kx="+kx+" ky="+ky+" t="+t+" b="+b);
        int ix = (kx-X_FIRST)/X_DELTA;
        int iy = (ky-Y_FIRST)/Y_DELTA;
        //System.out.println("ix="+ix+" iy="+iy);
        tb[0][ix][iy] = (t==99999.0f)?BOGUS:0.001f*t;
        tb[1][ix][iy] = (b==99999.0f)?BOGUS:0.001f*b;
        s.nextLine();
      }
      s.close();
      return tb;
    } catch (InputMismatchException ime) {
      throw new RuntimeException("Unknown format of file: "+FILENAME);
    }
  }
}
