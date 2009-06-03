package tp;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import static java.lang.Math.*;

import edu.mines.jtk.io.*;
import edu.mines.jtk.util.*;

/**
 * A gridded horizon from Teapot Dome.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.28
 */
public class Horizon {

  public int ns; // number of samples
  public float[] x1,x2,x3; // sample coordinates
  public int nt; // number of triangles
  public int[] ia,ib,ic; // triangle indices

  /**
   * Reads a horizon from a binary file with the specified name.
   * @param fileName the file name.
   * @return the horizon.
   */
  public static Horizon readBinary(String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      int ns = ais.readInt();
      float[] x1 = new float[ns]; ais.readFloats(x1);
      float[] x2 = new float[ns]; ais.readFloats(x2);
      float[] x3 = new float[ns]; ais.readFloats(x3);
      int nt = ais.readInt();
      int[] ia = new int[nt]; ais.readInts(ia);
      int[] ib = new int[nt]; ais.readInts(ib);
      int[] ic = new int[nt]; ais.readInts(ic);
      ais.close();
      return new Horizon(ns,x1,x2,x3,nt,ia,ib,ic);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a horizon from a text file with the specified name.
   * The file format is that used by Transform Software.
   * @param fileName the file name.
   * @return the horizon.
   */
  public static Horizon readText(String fileName) {
    double zdatum = Coordinates.DATUM_FT;
    double nullValue = -999.99;
    int m = 500, n = 500; // (m,n) = (inline,xline) dimensions
    int[][] grid = Array.fillint(-1,n,m);
    int ns = 0;
    FloatList x1List = new FloatList();
    FloatList x2List = new FloatList();
    FloatList x3List = new FloatList();
    try {
      BufferedReader br = new BufferedReader(new FileReader(fileName));
      String name = null;
      StringBuilder sb = null;
      String line = null;
      while ((line=br.readLine())!=null) {
        if (line.startsWith("#"))
          continue;
        String[] fields = line.split("\\s+");
        int i = Integer.parseInt(fields[0]);
        int j = Integer.parseInt(fields[1]);
        double xe = Double.parseDouble(fields[2]);
        double yn = Double.parseDouble(fields[3]);
        double zd = Double.parseDouble(fields[4]);
        assert grid[i][j]<0:"no duplicate samples";
        if (xe==nullValue || yn==nullValue || zd==nullValue)
          continue;
        Coordinates.Map xyz = new Coordinates.Map(xe,yn,zdatum-zd);
        Coordinates.Resampled r = new Coordinates.Resampled(xyz);
        x1List.add((float)r.x1);
        x2List.add((float)r.x2);
        x3List.add((float)r.x3);
        grid[i][j] = ns++;
      }
      br.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    float[] x1 = x1List.trim();
    float[] x2 = x2List.trim();
    float[] x3 = x3List.trim();
    IntList iaList = new IntList();
    IntList ibList = new IntList();
    IntList icList = new IntList();
    for (int i=1; i<m; ++i) {
      for (int j=1; j<n; ++j) {
        int imm = grid[i-1][j-1];
        int imp = grid[i-1][j  ];
        int ipm = grid[i  ][j-1];
        int ipp = grid[i  ][j  ];
        if (imm>=0 && imp>=0 && ipm>=0 && ipp>=0) {
          iaList.add(imp); ibList.add(imm); icList.add(ipm);
          iaList.add(imp); ibList.add(ipm); icList.add(ipp);
        } else if (  imp>=0 &&        imm>=0 &&        ipm>=0) {
          iaList.add(imp); ibList.add(imm); icList.add(ipm);
        } else if (  ipp>=0 &&        imp>=0 &&        imm>=0) {
          iaList.add(ipp); ibList.add(imp); icList.add(imm);
        } else if (  ipm>=0 &&        ipp>=0 &&        imp>=0) {
          iaList.add(ipm); ibList.add(ipp); icList.add(imp);
        } else if (  imm>=0 &&        ipm>=0 &&        ipp>=0) {
          iaList.add(imm); ibList.add(ipm); icList.add(ipp);
        }
      }
    }
    int[] ia = iaList.trim();
    int[] ib = ibList.trim();
    int[] ic = icList.trim();
    int nt = ia.length;
    return new Horizon(ns,x1,x2,x3,nt,ia,ib,ic);
  }

  /**
   * Writes this horizon to a binary file with the specified name.
   * @param fileName the file name.
   */
  public void writeBinary(String fileName) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      aos.writeInt(ns);
      aos.writeFloats(x1);
      aos.writeFloats(x2);
      aos.writeFloats(x3);
      aos.writeInt(nt);
      aos.writeInts(ia);
      aos.writeInts(ib);
      aos.writeInts(ic);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Gets packed triplets (ia,ib,ic) of triangle indices.
   * @return array[3*nt] of triangle indices.
   */
  public int[] getIABC() {
    int[] iabc = new int[3*nt];
    for (int it=0,m=0; it<nt; ++it) {
      iabc[m++] = ia[it];
      iabc[m++] = ib[it];
      iabc[m++] = ic[it];
    }
    return iabc;
  }

  /**
   * Gets packed triplets (x3,x2,x1) of sample coordinates.
   * @return array[3*ns] of sample coordinates.
   */
  public float[] getX321() {
    float[] x321 = new float[3*ns];
    for (int is=0,m=0; is<ns; ++is) {
      x321[m++] = x3[is];
      x321[m++] = x2[is];
      x321[m++] = x1[is];
    }
    return x321;
  }

  private Horizon(
    int ns, float[] x1, float[] x2, float[] x3,
    int nt, int[] ia, int[] ib, int[] ic) 
  {
    this.ns = ns; this.x1 = x1; this.x2 = x2; this.x3 = x3;
    this.nt = nt; this.ia = ia; this.ib = ib; this.ic = ic;
  }
}
