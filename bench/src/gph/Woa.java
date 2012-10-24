/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gph;

import java.io.*;
import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Utilities for working with World Ocean Atlas (WOA) data.
 * @author Dave Hale, Colorado School of Mines.
 * @version 2012.09.20
 */
public class Woa {

  public static Sampling getLonSampling() {
    return new Sampling(NLON,DLON,FLON);
  }

  public static Sampling getLatSampling() {
    return new Sampling(NLAT,DLAT,FLAT);
  }

  public static float[][] readData(String fileName, int level) {
    float[][] f = null;
    try {
      Scanner s = new Scanner(new FileInputStream(fileName));
      int nskip = (level-1)*NLON*NLAT/10; // 10 samples per line of file
      for (int i=0; i<nskip; ++i)
        s.nextLine();
      f = new float[NLAT][NLON];
      int n = 0;
      String line = null;
      for (int ilat=0; ilat<NLAT; ++ilat) {
        for (int ilon=0; ilon<NLON; ++ilon,++n) {
          int i = n%10;
          if (i==0)
            line = s.nextLine();
          String field = line.substring(i*8,(i+1)*8).trim();
          f[ilat][(ilon+180)%360] = Float.parseFloat(field);
        }
      }
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    return f;
  }

  public static boolean[][] mask(float[][] f, float fnull) {
    int n1 = f[0].length;
    int n2 = f.length;
    boolean[][] m = new boolean[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        m[i2][i1] = (f[i2][i1]!=fnull)?true:false;
      }
    }
    return m;
  }
  public static boolean[][] not(boolean[][] m) {
    int n1 = m[0].length;
    int n2 = m.length;
    m = copy(m);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        m[i2][i1] = !m[i2][i1];
      }
    }
    return m;
  }
  public static int count(boolean[][] m) {
    int n1 = m[0].length;
    int n2 = m.length;
    int n = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (m[i2][i1]) ++n;
      }
    }
    return n;
  }
  public static boolean[][] copy(boolean[][] m) {
    int n1 = m[0].length;
    int n2 = m.length;
    boolean[][] c = new boolean[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        c[i2][i1] = m[i2][i1];
      }
    }
    return c;
  }
  public static void copy(boolean[][] m, float[][] f, float[][] g) {
    int n1 = m[0].length;
    int n2 = m.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (m[i2][i1])
          g[i2][i1] = f[i2][i1];
      }
    }
  }
  public static float min(boolean[][] m, float[][] f) {
    int n1 = m[0].length;
    int n2 = m.length;
    float fmin = Float.MAX_VALUE;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fi = f[i2][i1];
        if (fi<fmin && m[i2][i1])
          fmin = fi;
      }
    }
    return fmin;
  }
  public static float max(boolean[][] m, float[][] f) {
    int n1 = m[0].length;
    int n2 = m.length;
    float fmax = -Float.MAX_VALUE;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fi = f[i2][i1];
        if (fi>fmax && m[i2][i1])
          fmax = fi;
      }
    }
    return fmax;
  }

  public static float[][] sdfix(float[][] n, float[][] sd) {
    int n1 = n[0].length;
    int n2 = n.length;
    float[][] td = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (n[i2][i1]>=2) {
          float ni = n[i2][i1];
          float si = sd[i2][i1];
          td[i2][i1] = sqrt(si*si*(ni-1.0f)/(ni-1.5f));
        }
      }
    }
    float td1 = nmed(2,10,n,td);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (n[i2][i1]==1)
          td[i2][i1] = td1;
      }
    }
    return td;
  }

  public static float[][] sefix(float[][] n, float[][] se) {
    float se1 = nmed(2,5,n,se);
    int n1 = n[0].length;
    int n2 = n.length;
    float[][] te = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (n[i2][i1]==1) {
          te[i2][i1] = se1;
        } else {
          te[i2][i1] = se[i2][i1];
        }
      }
    }
    return te;
  }

  /**
   * Returns median of values f for which n is in [nmin,nmax].
   */
  public static float nmed(
    float nmin, float nmax, float[][] n, float[][] f) 
  {
    // Count samples for which nmin <= n <= nmax.
    int n1 = n[0].length;
    int n2 = n.length;
    int m = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (nmin<=n[i2][i1] && n[i2][i1]<=nmax)
          ++m;
      }
    }

    // Get samples for which nmin <= n <= nmax.
    float[] g = new float[m];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (nmin<=n[i2][i1] && n[i2][i1]<=nmax)
          g[i++] = f[i2][i1];
      }
    }

    // Median.
    quickPartialSort(m/2,g);
    return g[m/2];
  }

  public static float[] clips(
    float pmin, float pmax, 
    boolean[][] m, float[][] f) 
  {
    if (pmin==0.0f && pmax==100.0f) {
      return new float[]{min(m,f),max(m,f)};
    } else {
      int n1 = m[0].length;
      int n2 = m.length;
      int n = count(m);
      float[] g = new float[n];
      for (int i2=0,j=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (m[i2][i1])
            g[j++] = f[i2][i1];
        }
      }
      Clips clips = new Clips(g);
      clips.setPercentiles(pmin,pmax);
      return new float[]{clips.getClipMin(),clips.getClipMax()};
    }
  }
  public static boolean[][] subset(float p, boolean[][] m) {
    return subset(new Random(),p,m);
  }
  public static boolean[][] subset(Random r, float p, boolean[][] m) {
    int n1 = m[0].length;
    int n2 = m.length;
    boolean[][] ms = new boolean[n2][n1];

    // Count number of live samples.
    int nm = count(m);
    if (nm==0)
      return ms;

    // Number of live samples in subset.
    int ns = 1+(int)(nm*p/100.0f);

    // Get arrays of indices of live samples.
    int[] m1 = new int[nm];
    int[] m2 = new int[nm];
    for (int i2=0,i=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (m[i2][i1]) {
          m1[i] = i1;
          m2[i] = i2;
          ++i;
        }
      }
    }

    // Make the first ns indices a random subset.
    for (int is=0; is<ns; ++is) {
      int js = is+r.nextInt(nm-is);
      int m1t = m1[is]; m1[is] = m1[js]; m1[js] = m1t;
      int m2t = m2[is]; m2[is] = m2[js]; m2[js] = m2t;
    }

    // Mark live samples in subset.
    for (int is=0; is<ns; ++is)
      ms[m2[is]][m1[is]] = true;
    return ms;
  }

  public static float[][] scattered(boolean[][] m, float[][] f) {
    Sampling slon = getLonSampling();
    Sampling slat = getLatSampling();
    int nlon = slon.getCount();
    int nlat = slat.getCount();
    int n = count(m);
    float[] fs = new float[n];
    float[] lons = new float[n];
    float[] lats = new float[n];
    for (int ilat=0,i=0; ilat<nlat; ++ilat) {
      float lati = (float)slat.getValue(ilat);
      for (int ilon=0; ilon<nlon; ++ilon) {
        float loni = (float)slon.getValue(ilon);
        if (m[ilat][ilon]) {
          fs[i] = f[ilat][ilon];
          lons[i] = loni;
          lats[i] = lati;
          ++i;
        }
      }
    }
    return new float[][]{fs,lons,lats};
  }

  public static EigenTensors2 makeTensors(boolean[][] m) {
    Sampling slon = getLonSampling();
    Sampling slat = getLatSampling();
    int nlon = slon.getCount();
    int nlat = slat.getCount();
    float[][] au = new float[nlat][nlon];
    float[][] av = new float[nlat][nlon];
    float[][] u1 = new float[nlat][nlon];
    float[][] u2 = new float[nlat][nlon];
    for (int ilat=0; ilat<nlat; ++ilat) {
      float lati = (float)slat.getValue(ilat);
      float cosl = cos(toRadians(lati));
      if (cosl<0.1f) cosl = 0.1f;
      float auc = 4.0f/(cosl*cosl);
      //System.out.println("cos="+cosl+" auc="+auc);
      for (int ilon=0; ilon<nlon; ++ilon) {
        float loni = (float)slon.getValue(ilon);
        float scale = m[ilat][ilon]?1.0f:0.001f;
        //au[ilat][ilon] = 1.0f;
        //av[ilat][ilon] = 1.0f;
        au[ilat][ilon] = scale*auc;
        av[ilat][ilon] = scale;
        u1[ilat][ilon] = 1.0f;
        u2[ilat][ilon] = 0.0f;
      }
    }
    return new EigenTensors2(u1,u2,au,av);
  }

  private static final int NLON = 360;
  private static final int NLAT = 180;
  private static final double DLON = 1.0;
  private static final double DLAT = 1.0;
  private static final double FLON = -179.5;
  private static final double FLAT = -89.5;
}
