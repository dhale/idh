package fault;

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING ONLY!
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.Stopwatch;
import java.util.*;
import javax.swing.*;

/** 
 * A Euclidean distance (and closest point) transform.
 * <p>
 * Based on the distance transform developed by 
 * Felzenszwalb, P.F. and D.P. Huttenlocher, 2004, 
 * Distance transforms of sampled functions: 
 * Cornell Computing and Information Science 
 * Technical Report TR2004-1963.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.01.17
 */
public class DistanceTransform {

  /**
   * Constructs a distance transform.
   */
  public DistanceTransform() {
    this(1.0);;
  }

  /**
   * Constructs a distance transform with specified scale factor.
   * @param s scale factor for all dimensions.
   */
  public DistanceTransform(double s) {
    this(s,s);
  }

  /**
   * Constructs a distance transform with specified scale factors.
   * @param s1 scale factor for 1st dimension.
   * @param s2 scale factor for 2nd and higher dimensions.
   */
  public DistanceTransform(double s1, double s2) {
    this(s1,s2,s2);
  }

  /**
   * Constructs a distance transform with specified scale factors.
   * @param s1 scale factor for 1st dimension.
   * @param s2 scale factor for 2nd dimension.
   * @param s3 scale factor for 3rd dimension.
   */
  public DistanceTransform(double s1, double s2, double s3) {
    _s1 = (float)s1;
    _s2 = (float)s2;
    _s3 = (float)s3;
  }

  /**
   * Applies this transform for a specified image.
   * @param fnull value used to mark samples in f that are unknown.
   * @param f input array of values with at least one value known.
   * @param d output of distances to the nearest known value.
   * @param i1 output array of 1st indices of nearest known value.
   * @param i2 output array of 2nd indices of nearest known value.
   */
  public void apply(
    final float fnull, final float[][] f, final float[][] d, 
    final short[][] i1, final short[][] i2) 
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    Parallel.loop(n2,new Parallel.LoopInt() { // axis 1
    public void compute(int j2) {
      float[] f1 = new float[n1];
      float[] d1 = new float[n1];
      short[] k1 = new short[n1];
      for (int j1=0; j1<n1; ++j1)
        f1[j1] = (f[j2][j1]==fnull)?HUGE:0.0f;
      dt(_s1,f1,d1,k1);
      for (int j1=0; j1<n1; ++j1) {
         d[j2][j1] = d1[j1];
        i1[j2][j1] = k1[j1];
      }
    }});
    Parallel.loop(n1,new Parallel.LoopInt() { // axis 2
    public void compute(int j1) {
      int[] ij = new int[n2];
      float[] f2 = new float[n2];
      float[] d2 = new float[n2];
      short[] k1 = new short[n2];
      short[] k2 = new short[n2];
      for (int j2=0; j2<n2; ++j2)
        f2[j2] = d[j2][j1];
      dt(_s2,f2,d2,k2);
      for (int j2=0; j2<n2; ++j2) {
        d[j2][j1] = sqrt(d2[j2]);
        k1[j2] = i1[k2[j2]][j1];
      }
      for (int j2=0; j2<n2; ++j2) {
        i1[j2][j1] = k1[j2];
        i2[j2][j1] = k2[j2];
      }
    }});
  }

  /**
   * Applies this transform for a specified image.
   * @param fnull value used to mark samples in f that are unknown.
   * @param f input array of values with at least one value known.
   * @param d output of distances to the nearest known value.
   * @param i1 output array of 1st indices of nearest known value.
   * @param i2 output array of 2nd indices of nearest known value.
   * @param i3 output array of 3rd indices of nearest known value.
   */
  public void apply(
    final float fnull, final float[][][] f, final float[][][] d, 
    final short[][][] i1, final short[][][] i2, final short[][][] i3) 
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int j3) {
      { // axis 1
        float[] f1 = new float[n1];
        float[] d1 = new float[n1];
        short[] k1 = new short[n1];
        for (int j2=0; j2<n2; ++j2) {
          for (int j1=0; j1<n1; ++j1)
              f1[j1] = (f[j3][j2][j1]==fnull)?HUGE:0.0f;
          dt(_s1,f1,d1,k1);
          for (int j1=0; j1<n1; ++j1) {
             d[j3][j2][j1] = d1[j1];
            i1[j3][j2][j1] = k1[j1];
          }
        }
      }
      { // axis 2
        float[] f2 = new float[n2];
        float[] d2 = new float[n2];
        short[] k1 = new short[n2];
        short[] k2 = new short[n2];
        for (int j1=0; j1<n1; ++j1) {
          for (int j2=0; j2<n2; ++j2)
            f2[j2] = d[j3][j2][j1];
          dt(_s2,f2,d2,k2);
          for (int j2=0; j2<n2; ++j2)
            k1[j2] = i1[j3][k2[j2]][j1];
          for (int j2=0; j2<n2; ++j2) {
             d[j3][j2][j1] = d2[j2];
            i1[j3][j2][j1] = k1[j2];
            i2[j3][j2][j1] = k2[j2];
          }
        }
      }
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int j2) {
      // axis 3
      float[] f3 = new float[n3];
      float[] d3 = new float[n3];
      short[] k1 = new short[n3];
      short[] k2 = new short[n3];
      short[] k3 = new short[n3];
      for (int j1=0; j1<n1; ++j1) {
        for (int j3=0; j3<n3; ++j3)
          f3[j3] = d[j3][j2][j1];
        dt(_s3,f3,d3,k3);
        for (int j3=0; j3<n3; ++j3) {
          k1[j3] = i1[k3[j3]][j2][j1];
          k2[j3] = i2[k3[j3]][j2][j1];
        }
        for (int j3=0; j3<n3; ++j3) {
           d[j3][j2][j1] = sqrt(d3[j3]);
          i1[j3][j2][j1] = k1[j3];
          i2[j3][j2][j1] = k2[j3];
          i3[j3][j2][j1] = k3[j3];
        }
      }
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _s1,_s2,_s3;

  private static float HUGE = 1.0e20f;

  /**
   * Algorithm DT(f) of Felzenszwalb and Huttenlocher (2004),
   * augmented to compute the index i of the nearest sample.
   * @param s the distance between two adjacent samples.
   * @param f input array of accumulated distances squared.
   * @param d output array of distances squared to nearest samples.
   * @param i output array of indices of nearest samples.
   */
  private static void dt(float s, float[] f, float[] d, short[] i) {
    int n = f.length;
    int[] v = new int[n];
    float[] z = new float[n+1];
    z[0] = -HUGE;
    z[1] =  HUGE;
    for (int q=1,k=0; q<n; ++q) {
      float r = ((f[q]+q*q)-(f[v[k]]+v[k]*v[k]))/(2*q-2*v[k]);
      while (r<=z[k]) {
        --k; 
        r = ((f[q]+q*q)-(f[v[k]]+v[k]*v[k]))/(2*q-2*v[k]);
      }
      ++k;
      v[k] = q;
      z[k] = r;
      z[k+1] = HUGE; 
    }
    for (int q=0,k=0; q<n; ++q) {
      while (z[k+1]<q)
        ++k;
      float qv = s*(q-v[k]);
      d[q] = qv*qv+f[v[k]];
      i[q] = (short)v[k];
    }
  }

  /*
  private static void dt1(float s, float[][] f, float[][] d, short[][] i) {
    int n = f.length;
    for (int p=0; p<n; ++p)
      dt1(s,f[p],d[p],i[p]);
  }

  private static void dt2(float s, float[][] f, float[][] d, short[][] i) {
    int m = f[0].length;
    int n = f.length;
    int[][] v = new int[n][m];
    float[][] z = new float[n+1][m];
    for (int p=0; p<m; ++p) {
      z[0][p] = -HUGE;
      z[1][p] =  HUGE;
    }
    for (int q=1,k=0; q<n; ++q) {
      for (int p=0; p<m; ++p) {
        int vk = v[k][p];
        float fq = f[q][p];
        float fv = f[vk][p];
        float r = ((fq+q*q)-(fv+vk*vk))/(2*q-2*vk);
        while (r<=z[k][p]) {
          --k; 
          vk = v[k][p];
          fv = f[vk][p];
          r = ((fq+q*q)-(fv+vk*vk))/(2*q-2*vk);
        }
        ++k;
        v[k][p] = q;
        z[k][p] = r;
        z[k+1][p] = HUGE; 
      }
    }
    for (int q=0,k=0; q<n; ++q) {
      for (int p=0; p<m; ++p) {
        while (z[k+1][p]<q)
          ++k;
        int vk = v[k][p];
        float qv = s*(q-vk);
        d[q][p] = qv*qv+f[vk][p];
        i[q][p] = (short)vk;
      }
    }
  }
  */

  ///////////////////////////////////////////////////////////////////////////
  // test

  private void applySlow(
    float fnull, float[][] f, float[][] d, 
    short[][] i1, short[][] i2) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    int n = 0;
    for (int j2=0; j2<n2; ++j2) {
      for (int j1=0; j1<n1; ++j1) {
        if (f[j2][j1]!=fnull)
          ++n;
      }
    }
    int[] k1 = new int[n];
    int[] k2 = new int[n];
    float[] fk = new float[n];
    for (int j2=0,j=0; j2<n2; ++j2) {
      for (int j1=0; j1<n1; ++j1) {
        if (f[j2][j1]!=fnull) {
          fk[j] = f[j2][j1];
          k1[j] = j1;
          k2[j] = j2;
          ++j;
        }
      }
    }
    for (int j2=0; j2<n2; ++j2) {
      for (int j1=0; j1<n1; ++j1) {
        float dsmin = Float.MAX_VALUE;
        int i1min = -1;
        int i2min = -1;
        for (int j=0; j<n; ++j) {
          float d1 = j1-k1[j];
          float d2 = j2-k2[j];
          float ds = d1*d1+d2*d2;
          if (ds<dsmin) {
            dsmin = ds;
            i1min = k1[j];
            i2min = k2[j];
          }
        }
        d[j2][j1] = sqrt(dsmin);
        i1[j2][j1] = (short)i1min;
        i2[j2][j1] = (short)i2min;
      }
    }
  }
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //test2Random();
        test3Random();
      }
    });
  }
  private static void test2Random() {
    int n1 = 1001;
    int n2 = 1002;
    int n = 100;
    Random r = new Random();
    int seed = r.nextInt();
    System.out.println("seed="+seed);
    r = new Random(seed);
    float[][] f = zerofloat(n1,n2);
    for (int i=0; i<n; ++i) {
      int i1 = r.nextInt(n1);
      int i2 = r.nextInt(n2);
      f[i2][i1] = 1.0f+i;
    }
    float[][] d = zerofloat(n1,n2);
    float[][] ds = zerofloat(n1,n2);
    short[][] i1 = zeroshort(n1,n2);
    short[][] i2 = zeroshort(n1,n2);
    DistanceTransform dt = new DistanceTransform();
    dt.applySlow(0.0f,f,ds,i1,i2);
    dt.apply(0.0f,f,d,i1,i2);
    System.out.println("error max="+max(abs(sub(ds,d))));
    float[][] p = zerofloat(n1,n2);
    for (int j2=0; j2<n2; ++j2) {
      for (int j1=0; j1<n1; ++j1) {
        p[j2][j1] = f[i2[j2][j1]][i1[j2][j1]];
        //p[j2][j1] = (float)i2[j2][j1]*n1+i1[j2][j1];
        //p[j2][j1] = (float)i1[j2][j1];
      }
    }
    SimplePlot spd = new SimplePlot();
    PixelsView pvd = spd.addPixels(d);
    pvd.setInterpolation(PixelsView.Interpolation.NEAREST);
    pvd.setColorModel(ColorMap.JET);
    spd.addColorBar();
    SimplePlot spp = new SimplePlot();
    PixelsView pvp = spp.addPixels(p);
    pvp.setInterpolation(PixelsView.Interpolation.NEAREST);
    pvp.setColorModel(ColorMap.JET);
    spp.addColorBar();
  }
  private static void test3Random() {
    int n1 = 401;
    int n2 = 402;
    int n3 = 403;
    int n = 10000;
    Random r = new Random();
    int seed = r.nextInt();
    System.out.println("seed="+seed);
    r = new Random(seed);
    float[][][] f = zerofloat(n1,n2,n3);
    for (int i=0; i<n; ++i) {
      int i1 = r.nextInt(n1);
      int i2 = r.nextInt(n2);
      int i3 = r.nextInt(n3);
      f[i3][i2][i1] = 1.0f+i;
    }
    float[][][] d = zerofloat(n1,n2,n3);
    float[][][] ds = zerofloat(n1,n2,n3);
    short[][][] i1 = zeroshort(n1,n2,n3);
    short[][][] i2 = zeroshort(n1,n2,n3);
    short[][][] i3 = zeroshort(n1,n2,n3);
    DistanceTransform dt = new DistanceTransform();
    //dt.applySlow(0.0f,f,ds,i1,i2,i3);
    Stopwatch sw = new Stopwatch();
    sw.restart();
    dt.apply(0.0f,f,d,i1,i2,i3);
    sw.stop();
    System.out.println("dt time = "+sw.time());
    //System.out.println("error max="+max(abs(sub(ds,d))));
    float[][][] p = zerofloat(n1,n2,n3);
    for (int j3=0; j3<n3; ++j3) {
      for (int j2=0; j2<n2; ++j2) {
        for (int j1=0; j1<n1; ++j1) {
          int i3j = i3[j3][j2][j1];
          int i2j = i2[j3][j2][j1];
          int i1j = i1[j3][j2][j1];
          p[j3][j2][j1] = f[i3j][i2j][i1j];
        }
      }
    }
    SimpleFrame sfd = new SimpleFrame();
    ImagePanelGroup ipd = sfd.addImagePanels(d);
    ipd.setColorModel(ColorMap.JET);
    SimpleFrame sfp = new SimpleFrame();
    ImagePanelGroup ipp = sfp.addImagePanels(p);
    ipp.setColorModel(ColorMap.JET);
  }
}
