package neil;

import java.io.IOException;
import javax.swing.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Demonstrates use of LocalShiftFinder.
 */
public class LsfDemo {

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        demo(); // demo has graphics, so run it on the Swing thread
      }
    });
  }

  public static void demo() {

    // Small array dimensions for testing.
    int n1 = 101, n2 = 103, n3 = 105;

    // Array g[n3][n2][n1] of random values.
    float[][][] g = makeRandom(n1,n2,n3);
    System.out.println("g min="+ min(g)+" max="+ max(g));

    // Array s1[n3][n2][n1] of displacements in 1st dimension.
    float[][][] s1 = makeShifts(n1,n2,n3);
    System.out.println("s1 min="+ min(s1)+" max="+ max(s1));

    // Warp g[n3][n2][n1] to get f[n3][n2][n1].
    float[][][] f = warp1(s1,g);
    System.out.println("f min="+ min(f)+" max="+ max(f));

    // Construct a local shift finder.
    float sigma = 4.0f;
    LocalShiftFinder lsf = new LocalShiftFinder(sigma);

    // Arrays for estimated displacements, initially zero.
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];

    // Iterative search for shifts in all dimensions.
    float[][][] du = new float[n3][n2][n1];
    float[][][] h = copy(g);
    int min1 = -2, max1 =  2;
    int min2 = -2, max2 =  2;
    int min3 = -2, max3 =  2;
    for (int iter=0; iter<2; ++iter) {
      lsf.find1(min1,max1,f,h,du);
      System.out.println("1: du min="+ min(du)+" max="+ max(du));
      lsf.shift1(du,u1,u2,u3,h);
      lsf.find2(min2,max2,f,h,du);
      System.out.println("2: du min="+ min(du)+" max="+ max(du));
      lsf.shift2(du,u1,u2,u3,h);
      lsf.find3(min3,max3,f,h,du);
      System.out.println("3: du min="+ min(du)+" max="+ max(du));
      lsf.shift3(du,u1,u2,u3,h);
    }
    System.out.println("u1 min="+ min(u1)+" max="+ max(u1));
    System.out.println("u2 min="+ min(u2)+" max="+ max(u2));
    System.out.println("u3 min="+ min(u3)+" max="+ max(u3));

    // Plots.
    plot("g = random image",g);
    plot("f = warped g",f);
    plot("actual shifts s1",s1);
    plot("estimated u1",u1);
    plot("estimated u2",u2);
    plot("estimated u3",u3);

    /*
    // Write arrays f and g to binary files. Then read them into new arrays 
    // f2 and g2. The max abs difference between f and f2 should be zero.
    writeFile("junkf.dat",f);
    writeFile("junkg.dat",g);
    float[][][] f2 = readFile("junkf.dat",n1,n2,n3);
    float[][][] g2 = readFile("junkg.dat",n1,n2,n3);
    System.out.println("max |f - f2| ="+max(abs(sub(f,f2))));
    */
  }

  // Plots one slice for i3 = n3/2 of a 3D array f[n3][n2][n1].
  private static void plot(String title, float[][][] f) {
    int n3 = f.length;
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setVLabel("i1");
    sp.setHLabel("i2");
    sp.setTitle(title);
    sp.setSize(600,600);
    PixelsView pv = sp.addPixels(f[n3/2]);
    pv.setClips(-1.0f,1.0f);
  }

  // Returns a 3D array r[n3][n2][n1] of Gaussian-filtered random values.
  private static float[][][] makeRandom(int n1, int n2, int n3) {
    float[][][] r = mul(5.0f, sub(0.5f, randfloat(n1,n2,n3)));
    float sigma = 1.0f;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply000(r,r);
    return r;
  }

  // Returns a displacement field u(x1,x2,x3) (an array u[n3][n2][n1]).
  private static float[][][] makeShifts(int n1, int n2, int n3) {
    float a0 = 0.000f;
    float a1 = 0.100f;
    float a2 = 0.100f;
    float a3 = 0.100f;
    float[][][] s = sin(rampfloat(a0,a1,a2,a3,n1,n2,n3));
    return s;
  }

  // Returns f[i3][i2][i1] = g[i3][i2](i1+s1[i3][i2][i1]), where
  // parentheses implies sinc interpolation in the 1st dimension.
  private static float[][][] warp1(float[][][] s1, float[][][] g) {
    int n1 = g[0][0].length, n2 = g[0].length, n3 = g.length;
    float[][][] x1 = rampfloat(0.0f,1.0f,0.0f,0.0f,n1,n2,n3);
    x1 = add(x1,s1);
    float[][][] f = zerofloat(n1,n2,n3);
    SincInterp si = new SincInterp();
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        si.interpolate(n1,1.0f,0.0f,g[i3][i2],n1,x1[i3][i2],f[i3][i2]);
      }
    }
    return f;
  }

  // Write a 3D array[n3][n2][n1] to a file.
  private static void writeFile(String filename, float[][][] f) {
    try {
        ArrayFile af = new ArrayFile(filename,"rw");
        af.writeFloats(f);
        af.close();
    } catch (IOException ioe) {
      // ...
    }
  }

  // Read a 3D array[n3][n2][n1] from a file.
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
