package ysun;

import java.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

public class First {
  public static void main(String[] args) {
    int n1 = 101;
    int n2 = 103;
    int n3 = 105;
    float[][][] f = makeRandom(n1,n2,n3);
    System.out.println("f min="+ ArrayMath.min(f)+" max="+ ArrayMath.max(f));
    SimplePlot.asPixels(f[n3/2]);
    float[][][] g = stretch(f);
    SimplePlot.asPixels(g[n3/2]);
    float sigma = 6.0f;
    float[][][] u = ArrayMath.zerofloat(n1,n2,n3); // new float[n3][n2][n1];
    LocalShiftFinder lsf = new LocalShiftFinder(sigma);
    int min1 = -2;
    int max1 =  2;
    lsf.find1(min1,max1,g,f,u);
    SimplePlot.asPixels(u[n3/2]);
    System.out.println("u min="+ ArrayMath.min(u)+" max="+ ArrayMath.max(u));
    /*
    writeFile("junkf.dat",f);
    writeFile("junkg.dat",g);
    float[][][] f2 = readFile("junkf.dat",n1,n2,n3);
    float[][][] g2 = readFile("junkg.dat",n1,n2,n3);
    System.out.println("max |f - f2| ="+ArrayMath.max(ArrayMath.abs(ArrayMath.sub(f,f2))));
    */
  }
  private static float[][][] makeRandom(int n1, int n2, int n3) {
    float[][][] f = ArrayMath.sub(0.5f, ArrayMath.randfloat(n1,n2,n3));
    float sigma = 1.0f;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply000(f,f);
    return f;
  }
  private static float[][][] makeT(int n1, int n2, int n3) {
    float a0 = 0.000f;
    float a1 = 0.100f;
    float a2 = 0.100f;
    float a3 = 0.100f;
    float[][][] t = 
      ArrayMath.mul(1.0f,
        ArrayMath.sin(
          ArrayMath.rampfloat(a0,a1,a2,a3,n1,n2,n3)));
    SimplePlot.asPixels(t[n3/2]);
    t = ArrayMath.add(t, ArrayMath.rampfloat(0.0f,1.0f,0.0f,0.0f,n1,n2,n3));
    return t;
  }
  private static float[][][] stretch(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] t = makeT(n1,n2,n3);
    float[][][] g = ArrayMath.zerofloat(n1,n2,n3);
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
            si.setUniform(n1,1.0f,0.0f,f[i3][i2]);
            si.interpolate(n1,t[i3][i2],g[i3][i2]);
        }
    }
    return g;
  }
  private static void writeFile(String filename, float[][][] f) {
    try {
        ArrayFile af = new ArrayFile(filename,"rw");
        af.writeFloats(f);
        af.close();
    } catch (IOException ioe) {
      // ...
    }
  }
  private static float[][][] readFile(String filename, int n1, int n2, int n3) {
    float[][][] f = new float[n3][n2][n1];
    try {
        ArrayFile af = new ArrayFile(filename,"r");
        af.readFloats(f);
        af.close();
    } catch (IOException ioe) {
      // ...
    }
    return f;
  }
}
