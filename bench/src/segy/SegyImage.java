/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package segy;

import java.io.*;
import java.nio.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A uniformly sampled seismic image in SEG-Y format.
 * The image may be two- or three-dimensional. The 1st dimension is
 * assumed to be time (or depth), and the 2nd dimension corresponds 
 * to the inline direction. For 3D images the 3rd dimension corresponds 
 * to the crossline direction.
 * <p>
 * All time and space coordinates are converted (if necessary) and 
 * stored as seconds and kilometers.
 * <p>
 * The image sampling grid has three indices i1, i2, and i3 in ranges
 * [0,n1-1], [i2min,i2max], and [i3min,i3max], corresponding to time, 
 * inline, and crossline dimensions, respectively. The minimum inline 
 * and crossline indices i2min and i3min in SEG-Y trace headers are 
 * typically non-zero. 
 * <p>
 * Image traces may be obtained from the image by specifying either a 
 * sequential trace index itrace or a pair of grid indices (i2,i3). 
 * For grid indices at which a trace does not exist, sample values 
 * are set to zero.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.06.18
 */
public class SegyImage {

  /**
   * Constructs an image with specified SEG-Y file name.
   * Assumes byte order is BIG_ENDIAN, which is the SEG-Y standard.
   * <p>
   * This method may 
   * @param fileName name of the SEG-Y file that contains the image.
   */
  public SegyImage(String fileName) {
    this(fileName,ByteOrder.BIG_ENDIAN);
  }

  /**
   * Constructs an image with specified SEG-Y file name and byte order.
   * @param fileName name of the SEG-Y file that contains the image.
   * @param byteOrder byte order, either BIG_ENDIAN or LITTLE_ENDIAN.
   */
  public SegyImage(String fileName, ByteOrder byteOrder) {
    _fileName = fileName;
    _byteOrder = byteOrder;
    openArrayFile();
    loadBinaryHeaderInfo();
    loadTraceHeaderInfo();
  }

  /**
   * Closes the SEG-Y file corresponding to this image.
   */
  public void close() {
    closeArrayFile();
  }

  /**
   * Prints information derived from binary file and trace headers.
   * To ensure that this information makes sense, this method is often 
   * the first one called after constructing an SEG-Y image. Nonsensical 
   * information may be caused by non-standard byte order or non-standard 
   * SEG-Y headers.
   */
  public void printInfo() {
    System.out.println("number of bytes = "+_nbyte);
    System.out.println("number of traces = "+_ntrace);
    if (_format==1) {
      System.out.println("format = 1 (4-byte IBM floating point)");
    } else if (_format==2) {
      System.out.println("format = 2 (4-byte two's complement integer)");
    } else if (_format==3) {
      System.out.println("format = 3 (2-byte two's complement integer)");
    } else if (_format==4) {
      System.out.println("format = 4 (4-byte fixed-point with gain)");
    } else if (_format==5) {
      System.out.println("format = 5 (4-byte IEEE floating point)");
    } else if (_format==8) {
      System.out.println("format = 8 (1-byte two's complement integer)");
    } else {
      System.out.println("format is unknown!");
    }
    System.out.println(
      "units for spatial coordinates: "+(_feet?"ft":"m")+
      " (will be converted to km)");
    System.out.printf("n1 = %5d (number of samples per trace)%n",_n1);
    System.out.printf("n2 = %5d (number of traces in inline direction)%n",_n2);
    System.out.printf("n3 = %5d (number of traces in crossline direction)%n",
      _n3);
    System.out.printf("d1 = %8.6f (time sampling interval, in s)%n",_d1);
    System.out.printf("d2 = %8.6f (inline sampling interval, in km)%n",_d2);
    System.out.printf("d3 = %8.6f (crossline sampling interval, in km)%n",_d3);
    System.out.printf("i2min = %5d, i2max = %5d (inline index bounds)%n",
      _i2min,_i2max);
    System.out.printf("i3min = %5d, i3max = %5d (crossline index bounds)%n",
      _i3min,_i3max);
    System.out.printf(
      "xmin = %10.6f, xmax = %10.6f (x coordinate bounds, in km)%n",
      _xmin,_xmax);
    System.out.printf(
      "ymin = %10.6f, ymax = %10.6f (y coordinate bounds, in km)%n",
      _ymin,_ymax);
    System.out.printf("grid azimuth = %6.2f degrees%n",toDegrees(_azim));
    System.out.println("grid reference point:");
    System.out.printf("  i2ref = %5d, i3ref = %5d, x = %10.6f, y = %10.6f%n",
      _i2ref,_i3ref,_xref,_yref);
    System.out.println("grid corner points:");
    System.out.printf("  i2min = %5d, i3min = %5d, x = %10.6f, y = %10.6f%n",
      _i2min,_i3min,getX(_i2min,_i3min),getY(_i2min,_i3min));
    System.out.printf("  i2max = %5d, i3min = %5d, x = %10.6f, y = %10.6f%n",
      _i2max,_i3min,getX(_i2max,_i3min),getY(_i2max,_i3min));
    System.out.printf("  i2min = %5d, i3max = %5d, x = %10.6f, y = %10.6f%n",
      _i2min,_i3max,getX(_i2min,_i3max),getY(_i2min,_i3max));
    System.out.printf("  i2max = %5d, i3max = %5d, x = %10.6f, y = %10.6f%n",
      _i2max,_i3max,getX(_i2max,_i3max),getY(_i2max,_i3max));
  }

  /**
   * Returns the number of bytes in this image.
   * @return the number of bytes.
   */
  public long countBytes() {
    return _nbyte;
  }

  /**
   * Returns the number of traces in this image.
   * @return the number of traces.
   */
  public int countTraces() {
    return _ntrace;
  }

  /**
   * Gets the format code for this image.
   * The code is 1 for IBM floats, 5 for IEEE floats.
   * @return the format code.
   */
  public int getFormat() {
    return _format;
  }

  /**
   * Sets the format code for this image.
   * Overrides the format found in the binary file header.
   * @param format the format code.
   */
  public void setFormat(int format) {
    _format = format;
  }

  /**
   * Gets the sampling of the 1st (time) dimension, in seconds.
   * @return the sampling.
   */
  public Sampling getSampling1() {
    return new Sampling(_n1,_d1,_f1);
  }

  /**
   * Gets the sampling of the 2nd (inline) dimension, in kilometers.
   * @return the sampling.
   */
  public Sampling getSampling2() {
    return new Sampling(_n2,_d2,_f2);
  }

  /**
   * Gets the sampling of the 3rd (crossline) dimension, in kilometers.
   * @return the sampling.
   */
  public Sampling getSampling3() {
    return new Sampling(_n3,_d3,_f3);
  }

  /**
   * Gets the number of samples in 1st dimension of image sampling grid.
   * @return the number of samples.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of samples in 2nd dimension of image sampling grid.
   * @return the number of samples.
   */
  public int getN2() {
    return _n2;
  }

  /**
   * Gets the number of samples in 3rd dimension of image sampling grid.
   * @return the number of samples.
   */
  public int getN3() {
    return _n3;
  }

  /**
   * Gets the sampling interval in 1st dimension of image sampling grid.
   * @return the sampling interval.
   */
  public double getD1() {
    return _d1;
  }

  /**
   * Gets the sampling interval in 2nd dimension of image sampling grid.
   * @return the sampling interval.
   */
  public double getD2() {
    return _d2;
  }

  /**
   * Gets the sampling interval in 3rd dimension of image sampling grid.
   * @return the sampling interval.
   */
  public double getD3() {
    return _d3;
  }

  /**
   * Sets the sampling interval for the 1st dimension.
   * Overrides the interval found in the binary file header.
   * @param d1 sampling interval for the 1st dimension.
   */
  public void setD1(double d1) {
    _d1 = d1;
  }

  /**
   * Sets the sampling interval for the 2nd dimension.
   * Overrides the interval computed from trace headers.
   * @param d2 sampling interval for the 2nd dimension.
   */
  public void setD2(double d2) {
    _d2 = d2;
  }

  /**
   * Sets the sampling interval for the 3rd dimension.
   * Overrides the interval computed from trace headers.
   * @param d3 sampling interval for the 3rd dimension.
   */
  public void setD3(double d3) {
    _d3 = d3;
  }

  /**
   * Gets the trace with specified index.
   * @param i the trace index.
   * @return array of trace samples.
   */
  public float[] getTrace(int i) {
    float[] f = new float[_n1];
    getTrace(i,f);
    return f;
  }

  /**
   * Gets the trace with specified grid indices.
   * If a trace does not exist for those indices, all samples are zero.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @return array of trace samples.
   */
  public float[] getTrace(int i2, int i3) {
    float[] f = new float[_n1];
    getTrace(i2,i3,f);
    return f;
  }

  /**
   * Gets the trace with specified index.
   * @param i the trace index.
   * @param f output array to fill with trace samples.
   */
  public void getTrace(int i, float[] f) {
    checkTraceIndex(i);
    try {
      _af.seek(offset(i));
      if (_format==1) { // 4-byte IBM floats
        _af.readInts(_ibuf);
        ibmToFloat(_ibuf,f);
      } else if (_format==2) { // 4-byte integers
        _af.readInts(_ibuf);
        intToFloat(_ibuf,f);
      } else if (_format==3) { // 2-byte integers
        _af.readShorts(_sbuf);
        shortToFloat(_sbuf,f);
      } else if (_format==4) { // 4-byte fixed-point with gain (obsolete)
        Check.state(_format!=4,"data sample format != 4 (obsolete)");
      } else if (_format==5) { // 4-byte IEEE floats
        _af.readFloats(f);
      } else if (_format==8) { // 1-byte integers
        _af.readBytes(_bbuf);
        byteToFloat(_bbuf,f);
      }
    } catch (IOException e) {
      throw new RuntimeException("cannot read trace at index "+i+" ("+e+")");
    }
  }

  /**
   * Gets the trace with specified grid indices.
   * If a trace does not exist for those indices, all samples are zero.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @param f output array to fill with trace samples.
   */
  public void getTrace(int i2, int i3, float[] f) {
    int i = index(i2,i3);
    if (i>=0) {
      getTrace(i,f);
    } else {
      for (int i1=0; i1<_n1; ++i1)
        f[i1] = 0.0f;
    }
  }

  /**
   * Determines if a trace exists with specified grid indices.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @return true, if the trace exists; false, otherwise.
   */
  public boolean hasTraceAt(int i2, int i3) {
    return e(i2,i3);
  }

  /**
   * Gets the x coordinate for the specified grid indices.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @return the x coordinate.
   */
  public double getX(int i2, int i3) {
    return _xref+(i2-_i2ref)*_d2*_cosa-(i3-_i3ref)*_d3*_sina;
  }

  /**
   * Gets the y coordinate for the specified grid indices.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @return the y coordinate.
   */
  public double getY(int i2, int i3) {
    return _yref+(i2-_i2ref)*_d2*_sina+(i3-_i3ref)*_d3*_cosa;
  }

  /**
   * Gets array of i2 indices, one for each trace.
   * @return array of i2 coordinates.
   */
  public int[] getI2s() {
    return copy(_i2s);
  }

  /**
   * Gets array of i3 indices, one for each trace.
   * @return array of i3 coordinates.
   */
  public int[] getI3s() {
    return copy(_i3s);
  }

  /**
   * Gets array of i2 indices, one for each trace, as floats.
   * @return array of i2 indices, as floats.
   */
  public float[] getI2sAsFloats() {
    return asFloats(_i2s);
  }

  /**
   * Gets array of i3 indices, one for each trace, as floats.
   * @return array of i3 indices, as floats.
   */
  public float[] getI3sAsFloats() {
    return asFloats(_i3s);
  }

  /**
   * Gets array of x coordinates, one for each trace.
   * @return array of x coordinates.
   */
  public double[] getXs() {
    return copy(_xs);
  }

  /**
   * Gets array of y coordinates, one for each trace.
   * @return array of y coordinates.
   */
  public double[] getYs() {
    return copy(_ys);
  }

  /**
   * Gets the minimum grid index for the 2nd dimension.
   * @return the minimum grid index.
   */
  public int getI2Min() {
    return _i2min;
  }

  /**
   * Gets the maximum grid index for the 2nd dimension.
   * @return the maximum grid index.
   */
  public int getI2Max() {
    return _i2max;
  }

  /**
   * Gets the minimum grid index for the 3rd dimension.
   * @return the minimum grid index.
   */
  public int getI3Min() {
    return _i3min;
  }

  /**
   * Gets the maximum grid index for the 3rd dimension.
   * @return the maximum grid index.
   */
  public int getI3Max() {
    return _i3max;
  }

  /**
   * Writes this image to a simple file of floats.
   * Writes zeros for any missing traces.
   */
  public void writeFloats(String fileName) {
    writeFloats(fileName,0,_n1-1,_i2min,_i2max,_i3min,_i3max);
  }

  /**
   * Writes a subset of this image to a simple file of floats.
   * Writes zeros for any missing traces; however, all minimum and 
   * maximum sample indices must be in the bounds for the grid.
   * @param fileName the file name.
   * @param i1min minimum sample index in 1st dimension.
   * @param i1max maximum sample index in 1st dimension.
   * @param i2min minimum sample index in 2nd dimension.
   * @param i2max maximum sample index in 2nd dimension.
   * @param i3min minimum sample index in 3rd dimension.
   * @param i3max maximum sample index in 3rd dimension.
   */
  public void writeFloats(
    String fileName,
    int i1min, int i1max,
    int i2min, int i2max,
    int i3min, int i3max)
  {
    checkSampleIndex(i1min);
    checkSampleIndex(i1max);
    checkGridIndices(i2min,i3min);
    checkGridIndices(i2max,i3max);
    Check.argument(i1min<=i1max,"i1min<=i1max");
    Check.argument(i2min<=i2max,"i2min<=i2max");
    Check.argument(i3min<=i3max,"i3min<=i3max");
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      int m1 = 1+i1max-i1min;
      float[] f = new float[_n1];
      float[] g = new float[m1];
      int mb = (int)(0.5+4.0*m1*(1+i2max-i2min)*(1+i3max-i3min)/1.0e6);
      System.out.print("writing floats ("+mb+" MB) ");
      for (int i3=i3min; i3<=i3max; ++i3) {
        if ((i3-i3min)%((i3max-i3min)/10)==0) {
          double perc = (int)((i3-i3min)*100.0/(i3max-i3min));
          System.out.print(".");
        }
        for (int i2=i2min; i2<=i2max; ++i2) {
          getTrace(i2,i3,f);
          copy(m1,i1min,f,0,g);
          aos.writeFloats(g);
        }
      }
      aos.close();
      System.out.println(" done");
    } catch (IOException e) {
      throw new RuntimeException("cannot write trace ("+e+")");
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private ArrayFile _af; // array file with random access
  private String _fileName; // SEG-Y file name
  private ByteOrder _byteOrder; // BIG_ENDIAN or LITTLE_ENDIAN
  private int _format; // sample format code
  private boolean _feet; // true, if feet; false if meters.
  private int _bytesPerSample; // number of bytes per sample
  private int _ntrace; // number of traces in file
  private long _nbyte; // number of bytes in file
  private Sampling _s1,_s2,_s3; // samplings
  private int _n1,_n2,_n3; // numbers of samples
  private double _d1,_d2,_d3; // sampling intervals
  private double _f1,_f2,_f3; // first sample values
  private int[] _i2s; // array[ntrace] of indices i2
  private int[] _i3s; // array[ntrace] of indices i3
  private double[] _xs; // array[ntrace] of x coordinates
  private double[] _ys; // array[ntrace] of y coordinates
  private int[][] _itrace; // array[n3][n2] of trace indices; -1 if no trace
  private int _i2min,_i2max; // bounds on index i2
  private int _i3min,_i3max; // bounds on index i3
  private double _xmin,_xmax; // bounds on coordinate x
  private double _ymin,_ymax; // bounds on coordinate y
  private int _i2ref,_i3ref; // grid indices for reference trace
  private double _xref,_yref; // coordinates for reference trace
  private double _azim = 0.0; // azimuth of sampling grid
  private double _cosa = 1.0; // cosine of grid azimuth
  private double _sina = 0.0; // sine of grid azimuth
  private int[] _ibuf; // buffer for trace samples as ints
  private short[] _sbuf; // buffer for trace samples as shorts
  private byte[] _bbuf; // buffer for trace samples as bytes

  private void println(String s) {
    System.out.println(s);
  }

  private float[] asFloats(int[] i) {
    int n = i.length;
    float[] f = new float[n];
    for (int j=0; j<n; ++j)
      f[j] = i[j];
    return f;
  }

  private boolean formatInt1() {
    return _format==8;
  }
  private boolean formatInt2() {
    return _format==3;
  }
  private boolean formatIbmFloat4() {
    return _format==1;
  }
  private boolean formatIeeeFloat4() {
    return _format==5;
  }

  private void openArrayFile() {
    try {
      if (_af==null)
        _af = new ArrayFile(_fileName,"r",_byteOrder,_byteOrder);
    } catch (IOException e) {
      throw new RuntimeException("cannot open SEG-Y file "+_fileName);
    }
  }
  private void closeArrayFile() {
    try {
      if (_af!=null)
        _af.close();
    } catch (IOException e) {
      throw new RuntimeException("cannot close SEG-Y file "+_fileName);
    }
    _af = null;
  }

  private void loadBinaryHeaderInfo() {
    short[] hs = new short[400/2]; // 400 bytes as 2-byte shorts
    try {
      _af.seek(3200L); // skip text header (3200 bytes)
      _af.readShorts(hs); // read binary header (400 bytes)
    } catch (IOException e) {
      throw new RuntimeException("cannot read binary header");
    }
    _format = hs[12]; // data sample format code
    if (_format==8) {
      _bytesPerSample = 1;
    } else if (_format==3) {
      _bytesPerSample = 2;
    } else if (_format==1 || _format==2 || _format==4 || _format==5) {
      _bytesPerSample = 4;
    } else {
      throw new RuntimeException("unknown data format: "+_format);
    }
    double d1 = hs[8]*1.0e-6; // sampling interval, in seconds
    if (d1<0.0 || d1>1.0) 
      throw new RuntimeException("invalid sampling interval: "+d1);
    int n1 = hs[10]; // number of samples per trace
    if (n1<0 || n1>100000) 
      throw new RuntimeException("invalid number of samples: "+n1);
    _feet = hs[27]==2; // feet or meters
    _n1 = n1;
    _d1 = d1;
    _nbyte = new File(_fileName).length();
    _ntrace = (int)((_nbyte-3200-400)/(240+_bytesPerSample*n1));
    _ibuf = new int[_n1];
    _sbuf = new short[_n1];
    _bbuf = new byte[_n1];
  }

  private void checkTraceIndex(int i) {
    Check.argument(i>=0,"trace index i is non-negative");
    Check.argument(i<_ntrace,"trace index i is less than number of traces");
  }

  private void checkSampleIndex(int i1) {
    Check.argument(i1>=0,"sample index i1 is non-negative");
    Check.argument(i1<_n1,"sample index i1 is less than number of samples");
  }

  private void checkGridIndices(int i2, int i3) {
    Check.argument(_i2min<=i2,"grid index i2 >= i2min");
    Check.argument(_i3min<=i3,"grid index i3 >= i3min");
    Check.argument(i2<=_i2max,"grid index i2 <= i2max");
    Check.argument(i3<=_i3max,"grid index i3 <= i3max");
  }

  private int index(int i2, int i3) {
    checkGridIndices(i2,i3);
    return _itrace[i3-_i3min][i2-_i2min];
  }
  private int indexForTrace(int i2, int i3) {
    int i = index(i2,i3);
    Check.state(i>=0,"have trace for indices ("+i2+","+i3+")");
    return i;
  }
  private boolean e(int i2, int i3) {
    return index(i2,i3)>=0;
  }
  private double x(int i2, int i3) {
    return _xs[indexForTrace(i2,i3)];
  }
  private double y(int i2, int i3) {
    return _ys[indexForTrace(i2,i3)];
  }
  private long offset(int i) {
    return 3600L+i*(240L+_bytesPerSample*_n1)+240L;
  }
  private long offset(int i2, int i3) {
    int i = index(i2,i3);
    return (i>=0)?offset(i):-1;
  }
  private void loadTraceHeaderInfo() {
    _i2min =  Integer.MAX_VALUE;
    _i2max = -Integer.MAX_VALUE;
    _i3min =  Integer.MAX_VALUE;
    _i3max = -Integer.MAX_VALUE;
    _xmin =  Double.MAX_VALUE;
    _xmax = -Double.MAX_VALUE;
    _ymin =  Double.MAX_VALUE;
    _ymax = -Double.MAX_VALUE;
    double uxy = 0.001*(_feet?0.3048:1.0);

    // Read indices (i2,i3) and coordinates (x,y) from trace headers.
    _i2s = new int[_ntrace];
    _i3s = new int[_ntrace];
    _xs = new double[_ntrace];
    _ys = new double[_ntrace];
    try {
      _af.seek(3600L); // skip text and binary file headers
      int[] hi = new int[240/4]; // 240-byte trace header as 4-byte ints
      System.out.print("reading "+_ntrace+" headers ");
      for (int itrace=0; itrace<_ntrace; ++itrace) {
        if ((itrace%(_ntrace/10))==0) {
          int perc = (int)(itrace*100.0/_ntrace);
          System.out.print(".");
        }
        _af.readInts(hi); // read the trace header
        _af.skipBytes(_n1*_bytesPerSample); // skip the trace samples
        int pxy = hi[17]&0xffff; // power of 10 is a 2-byte short
        double sxy = uxy*pow(10.0,pxy); // scale factor for x and y
        double x = hi[45]*sxy; // x coordinate
        double y = hi[46]*sxy; // y coordinate
        int i2 = hi[48]; // xline number
        int i3 = hi[47]; // iline number
        if (x<_xmin) _xmin = x;
        if (x>_xmax) _xmax = x;
        if (y<_ymin) _ymin = y;
        if (y>_ymax) _ymax = y;
        if (i2<_i2min) _i2min = i2;
        if (i2>_i2max) _i2max = i2;
        if (i3<_i3min) _i3min = i3;
        if (i3>_i3max) _i3max = i3;
        _xs[itrace] = x;
        _ys[itrace] = y;
        _i2s[itrace] = i2;
        _i3s[itrace] = i3;
      }
      System.out.println(" done");
    } catch (IOException e) {
      throw new RuntimeException("cannot read trace headers");
    }

    // Build mapping from grid indices to trace indices. With these
    // indices we can determine whether or not a trace exists and,
    // if it does, find x and y coordinates for any grid indices i2
    // and i3. A trace exists if the trace index is non-negative. A 
    // non-negative trace index can also be used to compute the byte 
    // offset for the corresponding trace in the SEG-Y file.
    _n2 = 1+_i2max-_i2min;
    _n3 = 1+_i3max-_i3min;
    _itrace = new int[_n3][_n2];
    for (int i3=0; i3<_n3; ++i3)
      for (int i2=0; i2<_n2; ++i2)
        _itrace[i3][i2] = -1;
    for (int itrace=0; itrace<_ntrace; ++itrace) {
      int i2 = _i2s[itrace]-_i2min;
      int i3 = _i3s[itrace]-_i3min;
      _itrace[i3][i2] = itrace;
    }

    // Compute sampling intervals d2 and d3 using largest spans 
    // of grid indices i2 and i3 for which traces exist. Using
    // the largest spans minimizes errors caused by rounding of
    // x and y coordinates. Compute also the grid azimuth and
    // grid reference points that enable computation of (x,y)
    // coordinates for any pair of grid indices (i2,i3).
    _d2 = 1.0;
    int k2 = 0;
    for (int j3=_i3min; j3<=_i3max; ++j3) {
      int j2lo = _i2min;
      while (j2lo<=_i2max && !e(j2lo,j3))
        ++j2lo;
      int j2hi = _i2max;
      while (j2hi>=_i2min && !e(j2hi,j3))
        --j2hi;
      if (j2hi-j2lo>k2) {
        k2 = j2hi-j2lo;
        double dx = x(j2hi,j3)-x(j2lo,j3);
        double dy = y(j2hi,j3)-y(j2lo,j3);
        _d2 = sqrt(dx*dx+dy*dy)/k2;
        _i2ref = j2lo;
        _i3ref = j3;
        _xref = x(_i2ref,_i3ref);
        _yref = y(_i2ref,_i3ref);
        _azim = atan2(dy,dx);
        _cosa = cos(_azim);
        _sina = sin(_azim);
      }
    }
    _d3 = 1.0;
    int k3 = 0;
    for (int j2=_i2min; j2<=_i2max; ++j2) {
      int j3lo = _i3min;
      while (j3lo<=_i3max && !e(j2,j3lo))
        ++j3lo;
      int j3hi = _i3max;
      while (j3hi>=_i3min && !e(j2,j3hi))
        --j3hi;
      if (j3hi-j3lo>k3) {
        k3 = j3hi-j3lo;
        double dx = x(j2,j3hi)-x(j2,j3lo);
        double dy = y(j2,j3hi)-y(j2,j3lo);
        _d3 = sqrt(dx*dx+dy*dy)/k3;
      }
    }
  }

  private static void byteToFloat(byte[] b, float[] f) {
    int n = b.length;
    for (int i=0; i<n; ++i)
      f[i] = b[i];
  }
  private static void shortToFloat(short[] s, float[] f) {
    int n = s.length;
    for (int i=0; i<n; ++i)
      f[i] = s[i];
  }
  private static void intToFloat(int[] i, float[] f) {
    int n = i.length;
    for (int j=0; j<n; ++j)
      f[j] = i[j];
  }
  private static void ieeeToFloat(int[] ieee, float[] f) {
    int n = ieee.length;
    for (int i=0; i<n; ++i)
      f[i] = ieeeToFloat(ieee[i]);
  }
  private static void ibmToFloat(int[] ibm, float[] f) {
    int n = ibm.length;
    for (int i=0; i<n; ++i)
      f[i] = ibmToFloat(ibm[i]);
  }
  private static float ieeeToFloat(int ieee) {
    return Float.intBitsToFloat(ieee);
  }
  private static float ibmToFloat(int ibm) {
    // 1) Extract sign bit, exponent, and mantissa.
    // 2) Convert exponent: subtract 64, multiply by 4, subtract 1, add 127
    // 3) Convert mantissa:
    //    while high mantissa bit is zero {
    //      shift mantissa left (low order bits are zeroed)
    //      decrement exponent
    // 4) Put sign and exponent bits back in number
    // 5) Reverse order of bytes?
    int s = 0x80000000&ibm; // sign bit
    int e = 0x7f000000&ibm; // exponent
    int m = 0x00ffffff&ibm; // mantissa
    int ieee = 0;
    if (m!=0) {
      e = (e>>22)-130; // = ((e>>24)-64)*4 - 1 + 127
      while ((m&0x00800000)==0) {
        m <<= 1;
        --e;
      }
      if (e<=0) {
        ieee = 0;
      } else if (e>=255) {
        ieee = s|0x7f7fffff;
      } else {
        ieee = s|(e<<23)|(m&0x007fffff);
      }
    }
    return Float.intBitsToFloat(ieee);
  }
  private static float ibmToFloatFromSuCurrentlyUnused(int ibm) {
    int fconv = ((ibm    )     )<<24 |
                ((ibm>> 8)&0xff)<<16 |
                ((ibm>>16)&0xff)<< 8 |
                ((ibm>>24)&0xff);
    int fmant = 0x00ffffff&fconv;
    if (fconv!=0 && fmant!=0) {
      int t = ((0x7f000000&fconv)>>>22)-130;
      while ((fmant&0x00800000)==0) {
        --t;
        fmant <<= 1;
      }
      if (t>254) {
        fconv = (0x80000000&fconv)|0x7f7fffff;
      } else if (t<=0) {
        fconv = 0;
      } else {
        fconv = (0x80000000&fconv)|(t<<23)|(0x007fffff&fmant);
      }
    }
    return Float.intBitsToFloat(fconv);
  }
}
