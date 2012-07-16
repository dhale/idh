/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package segy;

import java.io.*;
import java.nio.*;
import java.util.Random;

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
  }

  /**
   * Closes the SEG-Y file corresponding to this image.
   */
  public void close() {
    closeArrayFile();
  }

  /**
   * Prints summary information derived from the binary file header.
   * To ensure that this information makes sense, this method is often 
   * the first one called after constructing an SEG-Y image. Nonsensical 
   * information may be caused by non-standard byte order or non-standard 
   * SEG-Y headers.
   * <p>
   * For SEG-Y files with either IBM or IEEE floats, this method reads a 
   * small number of traces to guess the format code (because this code 
   * is often set incorrectly in the binary file header), and prints a 
   * warning if the guess does not match the format in the header.
   */
  public void printSummaryInfo() {
    loadBinaryHeaderInfo();
    System.out.println("****** beginning of SEG-Y file summary info ******");
    printBasicInfo();
    System.out.printf("n1 = %5d (number of samples per trace)%n",_n1);
    System.out.printf("d1 = %8.6f (time sampling interval, in s)%n",_d1);
    System.out.println("****** end of SEG-Y file summary info ******");
  }
  private void printBasicInfo() {
    System.out.println("file name = "+_fileName);
    System.out.println("byte order = "+_byteOrder);
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
    if (_formatGuess==1 && _format==5) {
      System.out.println("WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      System.out.println("WARNING: format may actually be 1 (IBM float)");
      System.out.println("WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    } else if (_formatGuess==5 && _format==1) {
      System.out.println("WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      System.out.println("WARNING: format may actually be 5 (IEEE float)");
      System.out.println("WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    }
    System.out.println(
      "units for spatial coordinates: "+(_feet?"ft":"m")+
      " (will be converted to km)");
  }

  /**
   * Prints all information derived from binary file and trace headers.
   * This method reads the entire SEG-Y file, which may take a long time.
   */
  public void printAllInfo() {
    loadTraceHeaderInfo();
    System.out.println("****** beginning of SEG-Y file info ******");
    printBasicInfo();
    System.out.println("indices and coordinates from trace headers:");
    System.out.printf(
      "  i2min = %5d, i2max = %5d (inline indices)%n",
      _i2min,_i2max);
    System.out.printf(
      "  i3min = %5d, i3max = %5d (crossline indices)%n",
      _i3min,_i3max);
    System.out.printf(
      "  xmin = %11.6f, xmax = %11.6f (x coordinates, in km)%n",
      _xmin,_xmax);
    System.out.printf(
      "  ymin = %11.6f, ymax = %11.6f (y coordinates, in km)%n",
      _ymin,_ymax);
    System.out.println("grid sampling:");
    System.out.printf(
      "  n1 = %5d (number of samples per trace)%n",_n1);
    System.out.printf(
      "  n2 = %5d (number of traces in inline direction)%n",_n2);
    System.out.printf(
      "  n3 = %5d (number of traces in crossline direction)%n",
      _n3);
    System.out.printf(
      "  d1 = %8.6f (time sampling interval, in s)%n",_d1);
    System.out.printf(
      "  d2 = %8.6f (inline sampling interval, in km)%n",_d2);
    System.out.printf(
      "  d3 = %8.6f (crossline sampling interval, in km)%n",_d3);
    System.out.println("grid corner points:");
    System.out.printf(
      "  i2min = %5d, i3min = %5d, x = %11.6f, y = %11.6f%n",
      _i2min,_i3min,getX(_i2min,_i3min),getY(_i2min,_i3min));
    System.out.printf(
      "  i2max = %5d, i3min = %5d, x = %11.6f, y = %11.6f%n",
      _i2max,_i3min,getX(_i2max,_i3min),getY(_i2max,_i3min));
    System.out.printf(
      "  i2min = %5d, i3max = %5d, x = %11.6f, y = %11.6f%n",
      _i2min,_i3max,getX(_i2min,_i3max),getY(_i2min,_i3max));
    System.out.printf(
      "  i2max = %5d, i3max = %5d, x = %11.6f, y = %11.6f%n",
      _i2max,_i3max,getX(_i2max,_i3max),getY(_i2max,_i3max));
    System.out.printf("grid azimuth: %5.2f degrees%n",_azim);
    System.out.println("****** end of SEG-Y file info ******");
  }

  /**
   * Prints the binary file header. Useful for diagnosing problems.
   * See the SEG-Y standard for complete explanations of all fields.
   */
  public void printBinaryHeader() {
    loadBinaryHeaderInfo();
    int[] h = new int[400/4];
    try {
      _af.seek(3200);
      _af.readInts(h);
    } catch (IOException e) {
      throw new RuntimeException("cannot read binary header");
    }
    System.out.println("****** beginning of binary file header ******");
    for (String field:_binaryHeaderFields) {
      int ibyte = Integer.parseInt(field.substring(0,4));
      int nbyte = Integer.parseInt(field.substring(5,9))-ibyte+1;
      int index = (ibyte-1-3200)/4;
      int ihalf = (ibyte-1-3200)%4;
      int value = h[index]; 
      if (nbyte==2)
        value = (ihalf==0)?h2(value):l2(value);
      String bytes = field.substring(0,11);
      String descr = field.substring(11);
      if (descr.equals("unassigned")) {
        System.out.printf(bytes+"           (unassigned)%n");
      } else {
        System.out.printf(bytes+"%10d ("+descr+")%n",value);
      }
    }
    System.out.println("****** end of binary file header ******");
  }
  
  /**
   * Prints a specified trace header. Useful for diagnosing problems.
   * See the SEG-Y standard for complete explanations of all fields.
   * @param i the (zero-based) trace index.
   */
  public void printTraceHeader(int i) {
    loadBinaryHeaderInfo();
    int[] h = new int[240/4];
    try {
      _af.seek(headerOffset(i));
      _af.readInts(h);
    } catch (IOException e) {
      throw new RuntimeException("cannot read header at index "+i+" ("+e+")");
    }
    System.out.println("****** beginning of header for trace "+i+" ******");
    for (String field:_traceHeaderFields) {
      int ibyte = Integer.parseInt(field.substring(0,4).trim());
      int nbyte = Integer.parseInt(field.substring(5,9).trim())-ibyte+1;
      int index = (ibyte-1)/4;
      int ihalf = (ibyte-1)%4;
      int value = h[index]; 
      if (nbyte==2)
        value = (ihalf==0)?h2(value):l2(value);
      String bytes = field.substring(0,11);
      String descr = field.substring(11);
      if (descr.equals("skipping")) {
        System.out.printf(bytes+"           (skipping)%n");
      } else {
        System.out.printf(bytes+"%10d ("+descr+")%n",value);
      }
    }
    System.out.println("****** end of header for trace "+i+" ******");
  }

  /**
   * Returns the number of bytes in this image.
   * @return the number of bytes.
   */
  public long countBytes() {
    loadBinaryHeaderInfo();
    return _nbyte;
  }

  /**
   * Returns the number of traces in this image.
   * @return the number of traces.
   */
  public int countTraces() {
    loadBinaryHeaderInfo();
    return _ntrace;
  }

  /**
   * Gets the format code for this image.
   * The code is 1 for IBM floats, 5 for IEEE floats.
   * @return the format code.
   */
  public int getFormat() {
    loadBinaryHeaderInfo();
    return _format;
  }

  /**
   * Sets the format code for this image.
   * Overrides the format found in the binary file header.
   * @param format the format code.
   */
  public void setFormat(int format) {
    _format = format;
    _formatSet = true;
    if (_format==8) {
      _bytesPerSample = 1;
    } else if (_format==3) {
      _bytesPerSample = 2;
    } else if (_format==1 || _format==2 || _format==4 || _format==5) {
      _bytesPerSample = 4;
    } else {
      throw new RuntimeException("unknown data format: "+_format);
    }
  }

  /**
   * Sets bytes used to find inline and xline numbers in trace headers.
   * Both inline and xline numbers are represented by four-byte integers, 
   * but these integers are not always placed in the standard bytes
   * 189-192 (for inline) and 193-196 (for xline) of trace headers. 
   * This method allows these standard byte locations to be overridden.
   * <p>
   * Note that the first byte is 1, not 0, as in the documentation 
   * for SEG-Y format.
   * @param inlineByte first byte of inline number.
   * @param xlineByte first byte of xline number.
   */
  public void setInlineXlineBytes(int inlineByte, int xlineByte) {
    _i2hi = (xlineByte-1)/4;
    _i3hi = (inlineByte-1)/4;
  }

  /**
   * Returns a guess for the format code, if a guess is possible.
   * Currently attempts to guess only if either IBM or IEEE floats.
   * <p>
   * The guess is obtained by analyzing amplitude spectra of a 
   * random sample of a specified number of traces, based on the
   * observation that using the wrong format when reading traces 
   * tends to introduce noise near the Nyquist frequency.
   * <p>
   * This method <emph>does not</emph> alter the format used to get
   * traces, which is obtained from the binary file header, and may 
   * be overridden only by explicitly setting it.
   * @param na number of traces to analyze.
   * @return the guess; 1 if IBM, 5 if IEEE, 0 if no guess.
   */ 
  public int guessFormat(int na) {
    loadBinaryHeaderInfo();

    // If neither IBM nor IEEE format, cannot guess.
    if (_format!=1 && _format!=5)
      return 0;

    // Time (t) and frequency (w) sampling.
    int nt = _n1;
    Fft fft = new Fft(nt);
    int nw = fft.getFrequencySampling1().getCount();

    // Guess will be positive if IBM, negative if IEEE, zero if unknown.
    int guess = 0;

    // Randomly select traces to analyze.
    Random r = new Random(3);

    // Maximum number of traces to read, including any all-zero traces.
    int nr = 10*na;

    // While more traces to analyze, and not too many already read, ...
    while (na>0 && nr>0) {

      // Read traces in both IBM and IEEE formats.
      int format = _format;
      int itrace = r.nextInt(_ntrace);
      _format = 1;
      float[] f1 = getTrace(itrace);
      _format = 5;
      float[] f5 = getTrace(itrace);
      _format = format;
      --nr;

      // Compute weighted sums of amplitude spectra of traces. 
      // One sum emphasizes low frequencies, while the other sum
      // emphasizes high frequencies. The correct format should 
      // minimize the ratio of the high- and low-frequency sums.
      // The weights for the high-frequency sum are nearly zero 
      // except near the Nyquist frequency, where they increase
      // rapidly.
      float[] a1 = cabs(fft.applyForward(f1));
      float[] a5 = cabs(fft.applyForward(f5));
      float swl = 0.0f;
      float swh = 0.0f;
      float s1l = 0.0f;
      float s5l = 0.0f;
      float s1h = 0.0f;
      float s5h = 0.0f;
      for (int iw=0; iw<nw; ++iw) {
        float wl = pow(1.0f-(float)iw/(nw-1),0.125f);
        float wh = 1.0f-wl; // ~ 0 except near Nyquist
        swl += wl;
        swh += wh;
        s1l += wl*a1[iw];
        s1h += wh*a1[iw];
        s5l += wl*a5[iw];
        s5h += wh*a5[iw];
      }
      s1l /= swl;
      s1h /= swh;
      s5l /= swl;
      s5h /= swh;

      // If trace not all zeros, compute ratios and update the guess.
      if (s1l>0.0f && s5l>0.0f) {
        float r1 = s1h/s1l;
        float r5 = s5h/s5l;
        if (r1<r5) {
          ++guess;
        } else if (r1>r5) {
          --guess;
        }
        --na;
      }
    }

    // Guess the format.
    if (guess>0) {
      return 1;
    } else if (guess<0) {
      return 5;
    } else {
      return 0;
    }
  }

  /**
   * Gets the number of samples in 1st dimension of image sampling grid.
   * @return the number of samples.
   */
  public int getN1() {
    loadBinaryHeaderInfo();
    return _n1;
  }

  /**
   * Gets the number of samples in 2nd dimension of image sampling grid.
   * @return the number of samples.
   */
  public int getN2() {
    loadTraceHeaderInfo();
    return _n2;
  }

  /**
   * Gets the number of samples in 3rd dimension of image sampling grid.
   * @return the number of samples.
   */
  public int getN3() {
    loadTraceHeaderInfo();
    return _n3;
  }

  /**
   * Gets the sampling interval in 1st dimension of image sampling grid.
   * @return the sampling interval.
   */
  public double getD1() {
    loadBinaryHeaderInfo();
    return _d1;
  }

  /**
   * Gets the sampling interval in 2nd dimension of image sampling grid.
   * @return the sampling interval.
   */
  public double getD2() {
    loadTraceHeaderInfo();
    return _d2;
  }

  /**
   * Gets the sampling interval in 3rd dimension of image sampling grid.
   * @return the sampling interval.
   */
  public double getD3() {
    loadTraceHeaderInfo();
    return _d3;
  }

  /**
   * Sets the sampling interval for the 1st dimension.
   * Overrides the interval found in the binary file header.
   * @param d1 sampling interval for the 1st dimension.
   */
  public void setD1(double d1) {
    _d1 = d1;
    _d1Set = true;
  }

  /**
   * Sets the sampling interval for the 2nd dimension.
   * Overrides the interval computed from trace headers.
   * @param d2 sampling interval for the 2nd dimension.
   */
  public void setD2(double d2) {
    _d2 = d2;
    _d2Set = true;
  }

  /**
   * Sets the sampling interval for the 3rd dimension.
   * Overrides the interval computed from trace headers.
   * @param d3 sampling interval for the 3rd dimension.
   */
  public void setD3(double d3) {
    _d3 = d3;
    _d3Set = true;
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
    loadBinaryHeaderInfo();
    checkTraceIndex(i);
    try {
      _af.seek(traceOffset(i));
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
    loadTraceHeaderInfo();
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
    loadTraceHeaderInfo();
    return e(i2,i3);
  }

  /**
   * Gets the x coordinate for the specified grid indices.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @return the x coordinate.
   */
  public double getX(int i2, int i3) {
    loadTraceHeaderInfo();
    return _xref+(i2-_i2ref)*_d2*_sin2+(i3-_i3ref)*_d3*_sin3;
  }

  /**
   * Gets the y coordinate for the specified grid indices.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3rd dimension.
   * @return the y coordinate.
   */
  public double getY(int i2, int i3) {
    loadTraceHeaderInfo();
    return _yref+(i2-_i2ref)*_d2*_cos2+(i3-_i3ref)*_d3*_cos3;
  }

  /**
   * Gets array of i2 indices, one for each trace.
   * @return array of i2 coordinates.
   */
  public int[] getI2s() {
    loadTraceHeaderInfo();
    return copy(_i2s);
  }

  /**
   * Gets array of i3 indices, one for each trace.
   * @return array of i3 coordinates.
   */
  public int[] getI3s() {
    loadTraceHeaderInfo();
    return copy(_i3s);
  }

  /**
   * Gets array of i2 indices, one for each trace, as floats.
   * @return array of i2 indices, as floats.
   */
  public float[] getI2sAsFloats() {
    loadTraceHeaderInfo();
    return asFloats(_i2s);
  }

  /**
   * Gets array of i3 indices, one for each trace, as floats.
   * @return array of i3 indices, as floats.
   */
  public float[] getI3sAsFloats() {
    loadTraceHeaderInfo();
    return asFloats(_i3s);
  }

  /**
   * Gets array of x coordinates, one for each trace.
   * @return array of x coordinates.
   */
  public double[] getXs() {
    loadTraceHeaderInfo();
    return copy(_xs);
  }

  /**
   * Gets array of y coordinates, one for each trace.
   * @return array of y coordinates.
   */
  public double[] getYs() {
    loadTraceHeaderInfo();
    return copy(_ys);
  }

  /**
   * Gets the minimum grid index for the 2nd dimension.
   * @return the minimum grid index.
   */
  public int getI2Min() {
    loadTraceHeaderInfo();
    return _i2min;
  }

  /**
   * Gets the maximum grid index for the 2nd dimension.
   * @return the maximum grid index.
   */
  public int getI2Max() {
    loadTraceHeaderInfo();
    return _i2max;
  }

  /**
   * Gets the minimum grid index for the 3rd dimension.
   * @return the minimum grid index.
   */
  public int getI3Min() {
    loadTraceHeaderInfo();
    return _i3min;
  }

  /**
   * Gets the maximum grid index for the 3rd dimension.
   * @return the maximum grid index.
   */
  public int getI3Max() {
    loadTraceHeaderInfo();
    return _i3max;
  }

  /**
   * Gets the minimum x coordinate.
   * @return the minimum x coordinate.
   */
  public double getXMin() {
    loadTraceHeaderInfo();
    return _xmin;
  }

  /**
   * Gets the maximum x coordinate.
   * @return the maximum x coordinate.
   */
  public double getXMax() {
    loadTraceHeaderInfo();
    return _xmax;
  }

  /**
   * Gets the minimum y coordinate.
   * @return the minimum y coordinate.
   */
  public double getYMin() {
    loadTraceHeaderInfo();
    return _ymin;
  }

  /**
   * Gets the maximum y coordinate.
   * @return the maximum y coordinate.
   */
  public double getYMax() {
    loadTraceHeaderInfo();
    return _ymax;
  }

  /**
   * Writes this image to a simple file of floats.
   * Writes zeros for any missing traces.
   */
  public void writeFloats(String fileName) {
    loadTraceHeaderInfo();
    writeFloats(fileName,1.0,0,_n1-1,_i2min,_i2max,_i3min,_i3max);
  }

  /**
   * Writes a subset of this image to a simple file of floats.
   * Writes zeros for any missing traces; however, all minimum and 
   * maximum sample indices must be in the bounds for the grid.
   * @param fileName the file name.
   * @param scaleFactor scaling applied to each sample.
   * @param i1min minimum sample index in 1st dimension.
   * @param i1max maximum sample index in 1st dimension.
   * @param i2min minimum sample index in 2nd dimension.
   * @param i2max maximum sample index in 2nd dimension.
   * @param i3min minimum sample index in 3rd dimension.
   * @param i3max maximum sample index in 3rd dimension.
   */
  public void writeFloats(
    String fileName,
    double scaleFactor,
    int i1min, int i1max,
    int i2min, int i2max,
    int i3min, int i3max)
  {
    loadTraceHeaderInfo();
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
      int m2 = 1+i2max-i2min;
      int m3 = 1+i3max-i3min;
      int mb = (int)(0.5+4.0*m1*m2*m3/1.0e6);
      System.out.print(
        "writing "+m1+"*"+m2+"*"+m3+" floats ("+mb+" MB) ");
      float[] f = new float[_n1];
      float[] g = new float[m1];
      float s = (float)scaleFactor;
      for (int i3=i3min; i3<=i3max; ++i3) {
        if (i3min<i3max && (i3-i3min)%((i3max-i3min)/10)==0) {
          double perc = (int)((i3-i3min)*100.0/(i3max-i3min));
          System.out.print(".");
        }
        for (int i2=i2min; i2<=i2max; ++i2) {
          getTrace(i2,i3,f);
          copy(m1,i1min,f,0,g);
          if (s!=1.0f) {
            for (int i1=0; i1<m1; ++i1)
              g[i1] *= s;
          }
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
  private int _formatGuess; // guess for format; zero if none
  private boolean _formatSet; // true, if format set explicitly
  private boolean _feet; // true, if feet; false if meters.
  private int _bytesPerSample; // number of bytes per sample
  private int _ntrace; // number of traces in file
  private long _nbyte; // number of bytes in file
  private Sampling _s1,_s2,_s3; // samplings
  private int _n1,_n2,_n3; // numbers of samples
  private double _d1,_d2,_d3; // sampling intervals
  private boolean _d1Set,_d2Set,_d3Set; // true, if set explicitly
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
  private double _azim = 90.0; // azimuth of sampling grid, in degrees
  private double _sin2 = 1.0; // sin(azimuth) for 2nd dimension
  private double _sin3 = 1.0; // sin(azimuth) for 3rd dimension
  private double _cos2 = 0.0; // cos(azimuth) for 2nd dimension
  private double _cos3 = 0.0; // cos(azimuth) for 3rd dimension
  private int[] _ibuf; // buffer for trace samples as ints
  private short[] _sbuf; // buffer for trace samples as shorts
  private byte[] _bbuf; // buffer for trace samples as bytes
  private boolean _infoBH; // true, if binary header info has been loaded
  private boolean _infoTH; // true, if trace header info has been loaded
  private int _i2hi = 48; // index in trace header of integer xline number
  private int _i3hi = 47; // index in trace header of integer iline number

  private static String[] _binaryHeaderFields = {
    "3201-3204: job identification number",
    "3205-3208: line number",
    "3209-3212: reel number",
    "3213-3214: number of data traces per ensemble",
    "3215-3216: number of auxiliary traces per ensemble",
    "3217-3218: sample interval in microseconds",
    "3219-3220: original sample interval in microseconds",
    "3221-3222: number of samples per data trace",
    "3223-3224: original number of samples per data trace",
    "3225-3226: data sample format code",
    "3227-3228: ensemble fold",
    "3229-3230: trace sorting code, 2=CDP, ...",
    "3231-3232: vertical sum code",
    "3233-3234: sweep frequency at start, in Hz",
    "3235-3236: sweep frequency at end, in Hz",
    "3237-3238: sweep length, in ms",
    "3239-3240: sweep type code",
    "3241-3242: trace number of sweep channel",
    "3243-3244: sweep taper at start, in ms",
    "3245-3246: sweep taper at end, in ms",
    "3247-3248: taper type",
    "3249-3250: correlated data traces",
    "3251-3252: binary gain recovered",
    "3253-3254: amplitude recovery method",
    "3255-3256: measurement system, 1=m, 2=ft",
    "3257-3258: impulse signal polarity, 1=neg, 2=pos",
    "3259-3260: vibratory polarity code",
    "3261-3500: unassigned",
    "3501-3502: SEG-Y format revision number",
    "3503-3504: fixed length trace flag",
    "3505-3506: number of 3200-byte header extensions",
    "3507-3600: unassigned",
  };

  private static String[] _traceHeaderFields = {
    "  1 -   4: trace sequence number within line",
    "  5 -   8: trace sequence number within file",
    "  9 -  12: original field record number",
    " 13 -  16: trace number in original field record",
    " 17 -  20: energy source point number",
    " 21 -  24: ensemble number",
    " 25 -  28: trace number within ensemble",
    " 29 -  30: trace identification code",
    " 31 -  32: number of vertically summed traces",
    " 33 -  34: number of horizontally stacked traces",
    " 35 -  36: data use: 1 = production, 2 = test",
    " 37 -  40: distance from source to receiver",
    " 41 -  44: receiver group elevation",
    " 45 -  48: surface elevation at source",
    " 49 -  52: source depth below surface",
    " 53 -  56: datum elevation at receiver group",
    " 57 -  60: datum elevation at source",
    " 61 -  64: water depth at source",
    " 65 -  68: water depth at group",
    " 69 -  70: scalar applied to elevations and depths",
    " 71 -  72: scalar applied to other coordinates",
    " 73 -  76: source coordinate x",
    " 77 -  80: source coordinate y",
    " 81 -  84: group coordinate x",
    " 85 -  88: group coordinate y",
    " 89 -  90: coordinate units, 1=length, ...",
    " 91 -  92: weathering velocity, in ft/s or m/s",
    " 93 -  94: subweathering velocity, in ft/s or m/s",
    " 95 -  96: uphole time at source, in ms",
    " 97 -  98: uphole time at group, in ms",
    " 99 - 100: source static correction, in ms",
    "101 - 102: group static correction, in ms",
    "103 - 104: total static applied, in ms",
    "105 - 106: lag time A, in ms",
    "107 - 108: lag time B, in ms",
    "109 - 110: delay recording time, in ms",
    "111 - 112: mute time start, in ms",
    "113 - 114: mute time end, in ms",
    "115 - 116: number of samples in this trace",
    "117 - 118: sample interval, in microseconds",
    "119 - 120: gain type of field instruments",
    "121 - 122: instrument gain constant, in dB",
    "123 - 124: instrument early or initial gain, in dB",
    "125 - 126: correlated, 1=no, 2=yes",
    "127 - 128: sweep frequency at start, in Hz",
    "129 - 130: sweep frequency at end, in Hz",
    "131 - 132: sweep length, in ms",
    "133 - 134: sweep type, 1=linear, ...",
    "135 - 136: sweep trace taper at start, in ms",
    "137 - 138: sweep trace taper at end, in ms",
    "139 - 140: taper type, 1=linear, ...",
    "141 - 142: alias filter frequency, in Hz",
    "143 - 144: alias filter slope, in dB/octave",
    "145 - 146: notch filter frequency, in Hz",
    "147 - 148: notch filter slope, in dB/octave",
    "149 - 150: low-cut frequency, in Hz",
    "151 - 152: high-cut frequency, in Hz",
    "153 - 154: low-cut slope, in dB/octave",
    "155 - 156: high-cut slope, in dB/octave",
    "157 - 158: year data recorded",
    "159 - 160: day of year",
    "161 - 162: hour of day",
    "163 - 164: minute of hour",
    "165 - 166: second of minute",
    "167 - 168: time basis code",
    "169 - 170: trace weighting factor",
    "171 - 172: geophone group number, roll position one",
    "173 - 174: geophone group number, first trace",
    "175 - 176: geophone group number, last trace",
    "177 - 178: gap size = total groups dropped",
    "179 - 180: over travel, 1=down/behind, 2=up/ahead",
    "181 - 184: x coordinate of ensemble position",
    "185 - 188: y coordinate of ensemble position",
    "189 - 192: 3D inline number, the line number",
    "193 - 196: 3D crossline number, the trace within line",
    "197 - 200: shotpoint location nearest to ensemble",
    "201 - 202: scalar for shotpoint location",
    "203 - 204: trace value measurement unit",
    "205 - 240: skipping",
  };

  private short l2(int i) {
    if (_byteOrder==ByteOrder.BIG_ENDIAN) {
      return (short)(i&0xffff);
    } else {
      return (short)((i>>16)&0xffff);
    }
  }
  private short h2(int i) {
    if (_byteOrder==ByteOrder.BIG_ENDIAN) {
      return (short)((i>>16)&0xffff);
    } else {
      return (short)(i&0xffff);
    }
  }

  private static float[] asFloats(int[] i) {
    int n = i.length;
    float[] f = new float[n];
    for (int j=0; j<n; ++j)
      f[j] = i[j];
    return f;
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
    if (_infoBH) 
      return;
    _infoBH = true;
    short[] hs = new short[400/2]; // 400 bytes as 2-byte shorts
    try {
      _af.seek(3200); // skip text header (3200 bytes)
      _af.readShorts(hs); // read binary header (400 bytes)
    } catch (IOException e) {
      throw new RuntimeException("cannot read binary header");
    }
    if (!_formatSet)
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
    if (!_d1Set) 
      _d1 = d1;
    _nbyte = new File(_fileName).length();
    _ntrace = (int)((_nbyte-3600)/(240+_bytesPerSample*n1));
    if (_ntrace*(240L+_bytesPerSample*n1)!=(_nbyte-3600L))
      throw new RuntimeException("invalid file length or binary header");
    _ibuf = new int[_n1];
    _sbuf = new short[_n1];
    _bbuf = new byte[_n1];
    _formatGuess = guessFormat(10);
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
  private long traceOffset(int i2, int i3) {
    int i = index(i2,i3);
    return (i>=0)?traceOffset(i):-1;
  }
  private long traceOffset(int i) {
    return headerOffset(i)+240L;
  }
  private long headerOffset(int i) {
    return 3600L+i*(240L+_bytesPerSample*_n1);
  }
  private void loadTraceHeaderInfo() {
    if (_infoTH) 
      return;
    _infoTH = true;
    loadBinaryHeaderInfo();
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
      _af.seek(3600); // skip text and binary file headers
      int[] hi = new int[240/4]; // 240-byte trace header as 4-byte ints
      System.out.print("reading "+_ntrace+" trace headers ");
      for (int itrace=0; itrace<_ntrace; ++itrace) {
        if ((itrace%(_ntrace/10))==0) {
          int perc = (int)(itrace*100.0/_ntrace);
          System.out.print(".");
        }
        _af.readInts(hi); // read the trace header
        _af.skipBytes(_n1*_bytesPerSample); // skip the trace samples
        double sxy = uxy; // scale factor for x and y
        int pxy = l2(hi[17]); // scale factor is a short
        if (pxy>0) // if positive, multiply
          sxy *= pxy;
        else if (pxy<0) // if negative, divide
          sxy /= -pxy;
        double x = hi[45]*sxy; // x coordinate
        double y = hi[46]*sxy; // y coordinate
        int i2 = hi[_i2hi]; // xline number
        if (i2==0) // if no xline number, ... 
          i2 = hi[5]; // try the CDP number
        int i3 = hi[_i3hi]; // iline number
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
        _i2ref = j2lo;
        _i3ref = j3;
        _xref = x(_i2ref,_i3ref);
        _yref = y(_i2ref,_i3ref);
        double dx = x(j2hi,j3)-x(j2lo,j3);
        double dy = y(j2hi,j3)-y(j2lo,j3);
        double dd = dx*dx+dy*dy;
        if (dd>0.0) {
          if (!_d2Set)
            _d2 = sqrt(dx*dx+dy*dy)/k2;
          double az = atan2(dx,dy);
          _sin2 = sin(az);
          _cos2 = cos(az);
          _azim = toDegrees(az);
        }
      }
    }
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
        double dd = dx*dx+dy*dy;
        if (dd>0.0) {
          if (!_d3Set)
            _d3 = sqrt(dx*dx+dy*dy)/k3;
          double az = atan2(dx,dy);
          _sin3 = sin(az);
          _cos3 = cos(az);
        }
      }
    }
  }
  
  private static float computeRatio(float[] a) {
    int len = a.length;
    float hSum  = 0;
    float lSum  = 0;
    float hwSum = 0;
    float lwSum = 0;
    for (int i = 0; i < len; i++) {
      float radians = (float)Math.toRadians((float)i / (float)(len-1) * 90f);
      float hw = (float)Math.pow(Math.sin(radians), 2);
      hwSum += hw;
      hSum  += hw*a[i];
			
      float lw = (float)Math.pow(Math.cos(radians), 2);
      lwSum += lw;
      lSum  += lw*a[i];
    }
    
    float h = hSum/hwSum;
    float l = lSum/lwSum;
    return h/l;
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
