package util;

import java.nio.ByteOrder;

import edu.mines.jtk.io.*;
import edu.mines.jtk.util.ArrayMath;

public class Floats {

  public static float[] readSubset(String fileName, int n1) {
    return null;
  }

  public static void subset(
    String infile, String outfile,
    int n1i, int n2i, int n3i,
    int n1o, int n2o, int n3o,
    int j1o, int j2o, int j3o,
    int k1o, int k2o, int k3o)
  {
    float[] fi = new float[n1i-j1o];
    float[] fo = new float[n1o];
    long offset1 = j1o*4;
    try {
      ArrayFile afi = new ArrayFile(infile,"r");
      ArrayFile afo = new ArrayFile(outfile,"rw");
      for (int i3o=0; i3o<n3o; ++i3o) {
        int i3i = j3o+i3o*k3o;
        long offset3 = i3i*n2i*n1i*4;
        for (int i2o=0; i2o<n2o; ++i2o) {
          int i2i = j2o+i2o*k2o;
          long offset2 = i2i*n1i*4;
          long offset = offset1+offset2+offset3;
          afi.seek(offset);
          if (k1o>1) {
            afi.readFloats(fi);
            ArrayMath.copy(n1o,0,k1o,fi,0,1,fo);
          } else {
            afi.readFloats(fo);
          }
          afo.writeFloats(fo);
        }
      }
      afi.close();
      afo.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static float[] read(String fileName, int n1) {
    return read(ByteOrder.BIG_ENDIAN,fileName,n1);
  }

  public static float[][] read(String fileName, int n1, int n2) {
    return read(ByteOrder.BIG_ENDIAN,fileName,n1,n2);
  }

  public static float[][][] read(String fileName, int n1, int n2, int n3) {
    return read(ByteOrder.BIG_ENDIAN,fileName,n1,n2,n3);
  }

  public static float[] readLittleEndian(String fileName, int n1) {
    return read(ByteOrder.LITTLE_ENDIAN,fileName,n1);
  }

  public static float[][] readLittleEndian(
    String fileName, int n1, int n2) 
  {
    return read(ByteOrder.LITTLE_ENDIAN,fileName,n1,n2);
  }

  public static float[][][] readLittleEndian(
    String fileName, int n1, int n2, int n3) 
  {
    return read(ByteOrder.LITTLE_ENDIAN,fileName,n1,n2,n3);
  }

  public static void write(String fileName, float[] f) {
    write(ByteOrder.BIG_ENDIAN,fileName,f);
  }

  public static void write(String fileName, float[][] f) {
    write(ByteOrder.BIG_ENDIAN,fileName,f);
  }

  public static void write(String fileName, float[][][] f) {
    write(ByteOrder.BIG_ENDIAN,fileName,f);
  }

  public static void writeLittleEndian(String fileName, float[] f) {
    write(ByteOrder.LITTLE_ENDIAN,fileName,f);
  }

  public static void writeLittleEndian(String fileName, float[][] f) {
    write(ByteOrder.LITTLE_ENDIAN,fileName,f);
  }

  public static void writeLittleEndian(String fileName, float[][][] f) {
    write(ByteOrder.LITTLE_ENDIAN,fileName,f);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Floats() {
  }

  private static float[] read(ByteOrder bo, String fileName, int n1) {
    float[] f = new float[n1];
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,bo);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return f;
  }

  private static float[][] read(
    ByteOrder bo, String fileName, int n1, int n2) 
  {
    float[][] f = new float[n2][n1];
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,bo);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return f;
  }

  private static float[][][] read(
    ByteOrder bo, String fileName, int n1, int n2, int n3) 
  {
    float[][][] f = new float[n3][n2][n1];
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,bo);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return f;
  }

  private static void write(ByteOrder bo, String fileName, float[] f) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,bo);
      aos.writeFloats(f);
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  private static void write(ByteOrder bo, String fileName, float[][] f) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,bo);
      aos.writeFloats(f);
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  private static void write(ByteOrder bo, String fileName, float[][][] f) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,bo);
      aos.writeFloats(f);
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }
}
