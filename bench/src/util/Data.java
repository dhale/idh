package util;

import edu.mines.jtk.io.*;

public class Data {

  public static float[][] readFloats(String fileName, int n1, int n2) {
    return readFloats(DataFile.ByteOrder.BIG_ENDIAN,fileName,n1,n2);
  }

  public static float[][] readFloatsLittleEndian(
    String fileName, int n1, int n2) 
  {
    return readFloats(DataFile.ByteOrder.LITTLE_ENDIAN,fileName,n1,n2);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Data() {
  }

  private static float[][] readFloats(
    DataFile.ByteOrder bo, String fileName, int n1, int n2) 
  {
    float[][] f = new float[n2][n1];
    try {
      DataFile df = new DataFile(fileName,"r",bo);
      df.readFloats(f);
      df.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return f;
  }
}
