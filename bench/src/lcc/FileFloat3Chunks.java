/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import java.io.File;
import java.io.IOException;

import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.util.ArrayMath;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.MathPlus.max;
import static edu.mines.jtk.util.MathPlus.min;

/**
 * A filter for 3-D arrays of floats in files too large to fit in memory.
 * This filter works by processing 3-D arrays in chunks that are not larger
 * than a specified maximum size. Input chunks may be overlapped and output 
 * chunks may be merged so that the output arrays in the files are seamless.
 * <p>
 * The files must support random access, so that floats may be read and
 * written in any order. In other words, the order may not be sequential.
 * <p>
 * The filter transforms zero or more 3-D input chunks into zero or more
 * 3D output chunks. The input and output chunks passed to the filter reside 
 * in arrays in memory, and may be smaller (but never larger) than the arrays 
 * stored in the corresponding files. All chunks passed to the filter have 
 * the same dimensions.
 * <p>
 * The filter is assumed to have an impulse response with finite-duration.
 * This property enables chunks to be filtered independently.
 * <p>
 * Most non-trivial filters compute each output sample from a 3-D window of 
 * input samples. To avoid seams in filtered output files, we must know the 
 * extent of this window. This extent is determined for each of the three 
 * array dimensions by two parameters: left and right overlaps. 
 * <p>
 * For the 1st dimension, the left overlap l1 is the number of input samples 
 * with indices less than i1 that are required to compute the output at index 
 * i1. Likewise, the right overlap r1 is the number of input samples with 
 * indices greater than index i1 required to compute the output sample at 
 * index i1. Corresponding left and right overlaps must be specified for
 * the 2nd and 3rd dimensions as well.
 * <p>
 * Chunks will overlap by these specified amounts, which implies that some 
 * input samples will be filtered more than once. The amount of wasted 
 * computation is proportional to the sum of left and right overlaps, which
 * may vary for different array dimensions. Chunk dimensions are computed 
 * to minimize this wasted computation.
 * <p>
 * Currently chunks are constructed by slicing only one array dimension.
 * Therefore, the maximum chunk size must be sufficiently large to contain 
 * complete slabs of the 3-D array, where the minimum slab thickness 
 * equals one plus the left and right overlaps. For efficiency (to 
 * reduce wasted computation) slabs should be much thicker than this.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.11.27
 */
public class FileFloat3Chunks {

  /**
   * Interface for a filter that processes zero or more 3-D input chunks 
   * to obtain zero or more 3-D output chunks.
   */
  public interface Filter {

    /**
     * Applies this filter to specified chunks.
     * @param i1 start index of chunks in 1st dimension.
     * @param i2 start index of chunks in 2nd dimension.
     * @param i3 start index of chunks in 3rd dimension.
     * @param x array of 3-D input chunks.
     * @param y array of 3-D output chunks.
     */
    public void apply(
      int i1, int i2, int i3, 
      float[][][][] x, float[][][][] y);
  }

  /**
   * Constructs a filter.
   * @param maxChunkSize the maximum size in floats of an input or output 
   *  chunk to be passed to the filter. The actual chunks (array) passed
   *  to the filter may be smaller but never larger than this size.
   * @param n1 number of samples in 1st dimension of arrays in files.
   * @param l1 number of samples of left overlap required in 1st dimension.
   * @param r1 number of samples of right overlap required in 1st dimension.
   * @param n2 number of samples in 2nd dimension of arrays in files.
   * @param l2 number of samples of left overlap required in 2nd dimension.
   * @param r2 number of samples of right overlap required in 2nd dimension.
   * @param n3 number of samples in 3rd dimension of arrays in files.
   * @param l3 number of samples of left overlap required in 3rd dimension.
   * @param r3 number of samples of right overlap required in 3rd dimension.
   * @exception IllegalArgumentException if the maxChunkSize is too small.
   */
  public FileFloat3Chunks(
    long maxChunkSize,
    int n1, int l1, int r1, 
    int n2, int l2, int r2, 
    int n3, int l3, int r3)
  {
    _mc = maxChunkSize;
    _n1 = n1;  _l1 = l1;  _r1 = r1;
    _n2 = n2;  _l2 = l2;  _r2 = r2;
    _n3 = n3;  _l3 = l3;  _r3 = r3;

    // Choose chunked dimension to minimize waste, while attempting to not
    // exceed the max chunk size limit. In case of ties, favor chunking of 
    // the 3rd dimension for which file I/O will be sequential. 
    int w1 = (l1+r1)*n2*n3; // waste in 1st dimension
    int w2 = (l2+r2)*n1*n3; // waste in 2nd dimension
    int w3 = (l3+r3)*n1*n2; // waste in 3rd dimension
    int m1 = (int)((_mc-w1)/(n2*n3)); // (l1+m1+r1)*n2*n3 <= mc
    int m2 = (int)((_mc-w2)/(n1*n3)); // (l2+m2+r2)*n1*n3 <= mc
    int m3 = (int)((_mc-w3)/(n1*n2)); // (l3+m3+r3)*n1*n2 <= mc
    if (m3>0 && (w3<=w2 || m2<=0) && (w3<=w1 || m1<=0)) {
      int mc = (l3+m3+r3)*n1*n2;
      Check.state(mc<=_mc,"chunk size <= max chunk size");
      _mc = mc;
      _m3 = m3;
    } else if (m2>0 && (w2<=w1 || m1<=0) && (w2<=w3 || m3<=0)) {
      int mc = (l2+m2+r2)*n1*n3;
      Check.state(mc<=_mc,"chunk size <= max chunk size");
      _mc = mc;
      _m2 = m2;
    } else if (m1>0 && (w1<=w2 || m2<=0) && (w1<=w3 || m3<=0)) {
      int mc = (l1+m1+r1)*n2*n3;
      Check.state(mc<=_mc,"chunk size <= max chunk size");
      _mc = mc;
      _m1 = m1;
    }
    Check.argument(_m1>0 || _m2>0 || _m3>0,"max chunk size is large enough");
  }

  /**
   * Transforms chunks in the specified files using the specified filter.
   * @param filter the filter.
   * @param xf array of files for input chunks.
   * @param yf array of files for output chunks.
   * @param ip indices for which in-place chunks y[i] == x[ip[i]];
   *  ip[i]&lt;0 indicates that chunk y[i] should be a new array of zeros.
   */
  public void apply(Filter filter, ArrayFile[] xf, ArrayFile[] yf, int[] ip) 
    throws IOException 
  {
    ip = (ip!=null)?ip: ArrayMath.fillint(-1,yf.length);
    if (_m1>0) {
      apply1(filter,xf,yf,ip);
    } else if (_m2>0) {
      apply2(filter,xf,yf,ip);
    } else if (_m3>0) {
      apply3(filter,xf,yf,ip);
    }
  }

  /**
   * Transforms chunks in the specified files using the specified filter.
   * @param filter the filter.
   * @param xf array of files for input chunks.
   * @param yf array of files for output chunks.
   */
  public void apply(Filter filter, ArrayFile[] xf, ArrayFile[] yf) 
    throws IOException 
  {
    apply(filter,xf,yf,null);
  }

  /**
   * Gets the array size (the chunk size).
   * @return the array size.
   */
  public long getChunkSize() {
    return _mc;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  
  private long _mc;
  private int _n1,_m1,_l1,_r1;
  private int _n2,_m2,_l2,_r2;
  private int _n3,_m3,_l3,_r3;

  private long byteOffset(int i1, int i2, int i3) {
    return 4*((long)i1+(long)i2*_n1+(long)i3*_n1*_n2);
  }

  private void apply1(Filter filter, ArrayFile[] xf, ArrayFile[] yf, int[] ip)
    throws IOException 
  {
    int nx = xf.length;
    int ny = yf.length;
    float[][][][] x = new float[nx][][][];
    float[][][][] y = new float[ny][][][];
    for (int i1=0; i1<_n1; i1+=_m1) {
      int l1 = min(_l1,i1); // 0 <= i1-l1
      int m1 = min(_m1,_n1-i1); // i1+m1 <= n1
      int r1 = min(_r1,_n1-i1-m1); // i1+m1+r1 <= n1
      int s1 = l1+m1+r1; // chunk size
      int j1 = i1-l1; // chunk index
      for (int ix=0; ix<nx; ++ix) {
        float[][][] xi = x[ix] = new float[_n3][_n2][s1];
        for (int i3=0; i3<_n3; ++i3) {
          for (int i2=0; i2<_n2; ++i2) {
            xf[ix].seek(byteOffset(j1,i2,i3));
            xf[ix].readFloats(xi[i3][i2]);
          }
        }
      }
      for (int iy=0; iy<ny; ++iy)
        y[iy] = ip[iy]<0?new float[_n3][_n2][s1]:x[ip[iy]];
      filter.apply(i1,0,0,x,y);
      for (int iy=0; iy<ny; ++iy) {
        float[][][] yi = y[iy];
        for (int i3=0; i3<_n3; ++i3) {
          for (int i2=0; i2<_n2; ++i2) {
            yf[iy].seek(byteOffset(i1,i2,i3));
            yf[iy].writeFloats(yi[i3][i2],l1,m1);
          }
        }
      }
    }
  }

  private void apply2(Filter filter, ArrayFile[] xf, ArrayFile[] yf, int[] ip)
    throws IOException 
  {
    int nx = xf.length;
    int ny = yf.length;
    float[][][][] x = new float[nx][][][];
    float[][][][] y = new float[ny][][][];
    for (int i2=0; i2<_n2; i2+=_m2) {
      int l2 = min(_l2,i2); // 0 <= i2-l2
      int m2 = min(_m2,_n2-i2); // i2+m2 <= n2
      int r2 = min(_r2,_n2-i2-m2); // i2+m2+r2 <= n2
      int s2 = l2+m2+r2; // chunk size
      int j2 = i2-l2; // chunk index
      for (int ix=0; ix<nx; ++ix) {
        float[][][] xi = x[ix] = new float[_n3][s2][_n1];
        for (int i3=0; i3<_n3; ++i3) {
          xf[ix].seek(byteOffset(0,j2,i3));
          xf[ix].readFloats(xi[i3]);
        }
      }
      for (int iy=0; iy<ny; ++iy)
        y[iy] = ip[iy]<0?new float[_n3][s2][_n1]:x[ip[iy]];
      filter.apply(0,i2,0,x,y);
      for (int iy=0; iy<ny; ++iy) {
        float[][][] yi = y[iy];
        for (int i3=0; i3<_n3; ++i3) {
          yf[iy].seek(byteOffset(0,i2,i3));
          for (int k2=0; k2<m2; ++k2)
            yf[iy].writeFloats(yi[i3][l2+k2]);
        }
      }
    }
  }

  private void apply3(Filter filter, ArrayFile[] xf, ArrayFile[] yf, int[] ip)
    throws IOException 
  {
    int nx = xf.length;
    int ny = yf.length;
    float[][][][] x = new float[nx][][][];
    float[][][][] y = new float[ny][][][];
    for (int i3=0; i3<_n3; i3+=_m3) {
      int l3 = min(_l3,i3); // 0 <= i3-l3
      int m3 = min(_m3,_n3-i3); // i3+m3 <= n3
      int r3 = min(_r3,_n3-i3-m3); // i3+m3+r3 <= n3
      int s3 = l3+m3+r3; // chunk size
      int j3 = i3-l3; // chunk index
      for (int ix=0; ix<nx; ++ix) {
        float[][][] xi = x[ix] = new float[s3][_n2][_n1];
        xf[ix].seek(byteOffset(0,0,j3));
        xf[ix].readFloats(xi);
      }
      for (int iy=0; iy<ny; ++iy)
        y[iy] = ip[iy]<0?new float[s3][_n2][_n1]:x[ip[iy]];
      filter.apply(0,0,i3,x,y);
      for (int iy=0; iy<ny; ++iy) {
        float[][][] yi = y[iy];
        yf[iy].seek(byteOffset(0,0,i3));
        for (int k3=0; k3<m3; ++k3)
          yf[iy].writeFloats(yi[l3+k3]);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // unused

  private void apply1(Filter filter, ArrayFile[] xf, ArrayFile[] yf)
    throws IOException 
  {
    int nx = xf.length;
    int ny = yf.length;
    float[][][][] x = new float[nx][][][];
    float[][][][] y = new float[ny][][][];
    for (int i1=0; i1<_n1; i1+=_m1) {
      int l1 = min(_l1,i1); // 0 <= i1-l1
      int m1 = min(_m1,_n1-i1); // i1+m1 <= n1
      int r1 = min(_r1,_n1-i1-m1); // i1+m1+r1 <= n1
      int s1 = l1+m1+r1; // chunk size
      int j1 = i1-l1; // chunk index
      for (int ix=0; ix<nx; ++ix) {
        float[][][] xi = x[ix] = new float[_n3][_n2][s1];
        for (int i3=0; i3<_n3; ++i3) {
          for (int i2=0; i2<_n2; ++i2) {
            xf[ix].seek(byteOffset(j1,i2,i3));
            xf[ix].readFloats(xi[i3][i2]);
          }
        }
      }
      for (int iy=0; iy<ny; ++iy)
        y[iy] = yarray(x,xf,yf[iy]);
      filter.apply(i1,0,0,x,y);
      for (int iy=0; iy<ny; ++iy) {
        float[][][] yi = y[iy];
        for (int i3=0; i3<_n3; ++i3) {
          for (int i2=0; i2<_n2; ++i2) {
            yf[iy].seek(byteOffset(i1,i2,i3));
            yf[iy].writeFloats(yi[i3][i2],l1,m1);
          }
        }
      }
    }
  }

  private void apply2(Filter filter, ArrayFile[] xf, ArrayFile[] yf)
    throws IOException 
  {
    int nx = xf.length;
    int ny = yf.length;
    float[][][][] x = new float[nx][][][];
    float[][][][] y = new float[ny][][][];
    for (int i2=0; i2<_n2; i2+=_m2) {
      int l2 = min(_l2,i2); // 0 <= i2-l2
      int m2 = min(_m2,_n2-i2); // i2+m2 <= n2
      int r2 = min(_r2,_n2-i2-m2); // i2+m2+r2 <= n2
      int s2 = l2+m2+r2; // chunk size
      int j2 = i2-l2; // chunk index
      for (int ix=0; ix<nx; ++ix) {
        float[][][] xi = x[ix] = new float[_n3][s2][_n1];
        for (int i3=0; i3<_n3; ++i3) {
          xf[ix].seek(byteOffset(0,j2,i3));
          xf[ix].readFloats(xi[i3]);
        }
      }
      for (int iy=0; iy<ny; ++iy)
        y[iy] = yarray(x,xf,yf[iy]);
      filter.apply(0,i2,0,x,y);
      for (int iy=0; iy<ny; ++iy) {
        float[][][] yi = y[iy];
        for (int i3=0; i3<_n3; ++i3) {
          yf[iy].seek(byteOffset(0,i2,i3));
          for (int k2=0; k2<m2; ++k2)
            yf[iy].writeFloats(yi[i3][l2+k2]);
        }
      }
    }
  }

  private void apply3(Filter filter, ArrayFile[] xf, ArrayFile[] yf)
    throws IOException 
  {
    int nx = xf.length;
    int ny = yf.length;
    float[][][][] x = new float[nx][][][];
    float[][][][] y = new float[ny][][][];
    for (int i3=0; i3<_n3; i3+=_m3) {
      int l3 = min(_l3,i3); // 0 <= i3-l3
      int m3 = min(_m3,_n3-i3); // i3+m3 <= n3
      int r3 = min(_r3,_n3-i3-m3); // i3+m3+r3 <= n3
      int s3 = l3+m3+r3; // chunk size
      int j3 = i3-l3; // chunk index
      for (int ix=0; ix<nx; ++ix) {
        float[][][] xi = x[ix] = new float[s3][_n2][_n1];
        xf[ix].seek(byteOffset(0,0,j3));
        xf[ix].readFloats(xi);
      }
      for (int iy=0; iy<ny; ++iy)
        y[iy] = yarray(x,xf,yf[iy]);
      filter.apply(0,0,i3,x,y);
      for (int iy=0; iy<ny; ++iy) {
        float[][][] yi = y[iy];
        yf[iy].seek(byteOffset(0,0,i3));
        for (int k3=0; k3<m3; ++k3)
          yf[iy].writeFloats(yi[l3+k3]);
      }
    }
  }

  private float[][][] yarray(float[][][][] x, ArrayFile[] xf, ArrayFile yf) {
    int n = x.length;
    for (int i=0; i<n; ++i) {
      if (yf==xf[i])
        return x[i];
    }
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    return new float[n3][n2][n1];
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  // Implements a simple moving average filter.
  private static class TestFilter implements Filter {
    TestFilter(int l1, int r1, int l2, int r2, int l3, int r3) {
      _l1 = l1;  _r1 = r1;
      _l2 = l2;  _r2 = r2;
      _l3 = l3;  _r3 = r3;
    }
    int _l1,_l2,_l3;
    int _r1,_r2,_r3;
    public void apply(
      int i1, int i2, int i3, 
      float[][][][] x, float[][][][] y) 
    {
      float[][][] xi = x[0];
      float[][][] yi = y[0];
      int n1 = xi[0][0].length;
      int n2 = xi[0].length;
      int n3 = xi.length;
      System.out.println("TestFilter.apply: i1="+i1+" i2="+i2+" i3="+i3);
      System.out.println("                  n1="+n1+" n2="+n2+" n3="+n3);
      for (int j3=0; j3<n3; ++j3) {
        int k3lo = max(0,j3-_l3);
        int k3hi = min(n3-1,j3+_r3);
        for (int j2=0; j2<n2; ++j2) {
          int k2lo = max(0,j2-_l2);
          int k2hi = min(n2-1,j2+_r2);
          for (int j1=0; j1<n1; ++j1) {
            int k1lo = max(0,j1-_l1);
            int k1hi = min(n1-1,j1+_r1);
            float sum = 0.0f;
            for (int k3=k3lo; k3<=k3hi; ++k3) {
              for (int k2=k2lo; k2<=k2hi; ++k2) {
                for (int k1=k1lo; k1<=k1hi; ++k1) {
                  sum += xi[k3][k2][k1];
                }
              }
            }
            yi[j3][j2][j1] = sum;
          }
        }
      }
    }
  }

  private static void testFilter(
    int maxChunkSize,
    int n1, int l1, int r1,
    int n2, int l2, int r2,
    int n3, int l3, int r3) 
    throws IOException 
  {
    float[][][] x = ArrayMath.randfloat(n1,n2,n3);
    float[][][] y = ArrayMath.zerofloat(n1,n2,n3);
    float[][][] z = ArrayMath.zerofloat(n1,n2,n3);
    TestFilter tf = new TestFilter(l1,r1,l2,r2,l3,r3);
    File xfile = null;
    File yfile = null;
    ArrayFile xaf = null;
    ArrayFile yaf = null;
    try {
      tf.apply(0,0,0,new float[][][][]{x},new float[][][][]{z});
      xfile = File.createTempFile("junkx","dat");
      yfile = File.createTempFile("junky","dat");
      xaf = new ArrayFile(xfile,"rw");
      yaf = new ArrayFile(yfile,"rw");
      xaf.writeFloats(x);
      FileFloat3Chunks ff3c =
        new FileFloat3Chunks(maxChunkSize,n1,l1,r1,n2,l2,r2,n3,l3,r3);
      ff3c.apply(tf,new ArrayFile[]{xaf},new ArrayFile[]{yaf});
      yaf.seek(0);
      yaf.readFloats(y);
      float xsum = ArrayMath.sum(x);
      float ysum = ArrayMath.sum(y);
      float zsum = ArrayMath.sum(z);
      float emax = ArrayMath.max(ArrayMath.abs(ArrayMath.sub(y,z)));
      System.out.println("xsum = "+xsum);
      System.out.println("ysum = "+ysum);
      System.out.println("zsum = "+zsum);
      System.out.println("emax = "+emax);
    } finally {
      if (xaf!=null) xaf.close();
      if (yaf!=null) yaf.close();
      if (xfile!=null) xfile.delete();
      if (yfile!=null) yfile.delete();
    }
  }

  private static void test1() throws IOException {
    int n1 = 12,  l1 = 1,  r1 = 2;
    int n2 = 11,  l2 = 2,  r2 = 3;
    int n3 = 10,  l3 = 3,  r3 = 4;
    int maxChunkSize = (l1+n1/4+r1)*n2*n3;
    System.out.println("test1:");
    testFilter(maxChunkSize,n1,l1,r1,n2,l2,r2,n3,l3,r3);
  }

  private static void test2() throws IOException {
    int n1 = 11,  l1 = 2,  r1 = 3;
    int n2 = 12,  l2 = 1,  r2 = 2;
    int n3 = 10,  l3 = 3,  r3 = 4;
    int maxChunkSize = (l2+n2/4+r2)*n1*n3;
    System.out.println("test2:");
    testFilter(maxChunkSize,n1,l1,r1,n2,l2,r2,n3,l3,r3);
  }

  private static void test3() throws IOException {
    int n1 = 11,  l1 = 2,  r1 = 3;
    int n2 = 10,  l2 = 3,  r2 = 4;
    int n3 = 12,  l3 = 1,  r3 = 2;
    int maxChunkSize = (l3+n3/4+r3)*n1*n2;
    System.out.println("test3:");
    testFilter(maxChunkSize,n1,l1,r1,n2,l2,r2,n3,l3,r3);
  }

  public static void main(String[] args) throws IOException {
    test1();
    test2();
    test3();
  }
}
