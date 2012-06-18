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

/**
 * A file containing seismic data in the SEG-Y format.
 * The file is currently assumed to contain a post-stack 
 * seismic image.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.06.18
 */
public class SegyFile {

  public SegyFile(String fileName, ByteOrder byteOrder) {
    _fileName = fileName;
    _byteOrder = byteOrder;
    testSegyFileAccess();
    loadBinaryHeaderInfo();
  }

  /**
   * Gets the sampling of the 1st (time) dimension, in seconds.
   * @return the sampling.
   */
  public Sampling getSampling1() {
    return _s1;
  }

  /**
   * Gets the sampling of the 2nd (inline) dimension, in kilometers.
   * @return the sampling.
   */
  public Sampling getSampling2() {
    return _s2;
  }

  /**
   * Gets the sampling of the 3rd (crossline) dimension, in kilometers.
   * @return the sampling.
   */
  public Sampling getSampling3() {
    return _s3;
  }

  public void getTrace(float[] f) {
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private String _fileName; // SEG-Y file name
  private ByteOrder _byteOrder; // BIG_ENDIAN or LITTLE_ENDIAN
  private int _format; // sample format code
  private boolean _feet; // true, if feet; false if meters.
  private int _bytesPerSample; // number of bytes per sample
  private int _ntrace; // number of traces in file
  private long _nbyte; // number of bytes in file
  private Sampling _s1,_s2,_s3; // samplings for time, inline and crossline
  private boolean[][] _exists; // true if trace exists; false, otherwise

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

  private void testSegyFileAccess() {
    try {
      ArrayInputStream ais = new ArrayInputStream(_fileName,_byteOrder);
      ais.skipBytes(3200);
      ais.close();
    } catch (IOException e) {
      throw new RuntimeException("Cannot read SEG-Y file: "+_fileName);
    }
  }

  private void loadBinaryHeaderInfo() {
    short[] hs = new short[400/2]; // 400 bytes as 2-byte shorts
    try {
      ArrayInputStream ais = new ArrayInputStream(_fileName,_byteOrder);
      ais.skipBytes(3200); // skip text header (3200 bytes)
      ais.readShorts(hs); // read binary header (400 bytes)
      ais.close();
    } catch (IOException e) {
      throw new RuntimeException("Cannot read binary header.");
    }
    _format = hs[12]; // data sample format code
    if (_format==8) {
      _bytesPerSample = 1;
    } else if (_format==3) {
      _bytesPerSample = 2;
    } else if (_format==2 || _format==4 || _format==5) {
      _bytesPerSample = 4;
    } else {
      throw new RuntimeException("Invalid data format: "+_format);
    }
    double d1 = hs[8]*1.0e-6; // sampling interval, in seconds
    if (d1<0.0 || d1>1.0) 
      throw new RuntimeException("Invalid sampling interval: "+d1);
    int n1 = hs[10]; // number of samples per trace
    if (n1<0 || n1>100000) 
      throw new RuntimeException("Invalid number of samples: "+n1);
    _feet = hs[27]==2; // feet or meters
    _s1 = new Sampling(n1,d1,0.0); // time sampling
    _nbyte = new File(_fileName).length();
    _ntrace = (int)((_nbyte-3200-400)/(240+_bytesPerSample*n1));
  }

  private static void shortToFloat(short[] s, float[] f) {
    int n = s.length;
    for (int i=0; i<n; ++i)
      f[i] = s[i];
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
