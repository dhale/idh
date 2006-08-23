/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gp404;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.MathPlus.*;

/**
 * Reads files containing sequences - sampled functions of one variable.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.08.22
 */
public class SequenceReader {

  public static void main(String[] args) {
    String fileName = "hunter1.txt";
    int column = 1;
    if (args.length>0)
      fileName = args[0];
    if (args.length>1)
      column = Integer.parseInt(args[1]);
    Sequence x = read(fileName,column);
    Sampling s = x.getSampling();
    int nt = s.getCount();
    double dt = s.getDelta();
    double ft = s.getFirst();
    System.out.println("In "+fileName);
    System.out.println("  nt="+nt+" dt="+dt+" ft="+ft);
  }

  /**
   * Reads a sequence from the specified file.
   * @param fileName name of file containing the sequence.
   * @return the sequence.
   */
  public static Sequence read(String fileName) {
    return read(fileName,1);
  }

  /**
   * Reads a sequence from the specified column of the specified file.
   * @param fileName name of file containing the sequence.
   * @return the sequence.
   */
  public static Sequence read(String fileName, int column) {
    BufferedReader br = null;
    try {
      java.net.URL url = SequenceReader.class.getResource("data/"+fileName);
      Reader reader;
      if (url!=null) {
        reader = new InputStreamReader(url.openStream());
      } else {
        reader = new FileReader(fileName);
      }
      br = new BufferedReader(reader);
    } catch (IOException ioe) {
      throw new RuntimeException("Cannot open file: "+fileName);
    }
    try {
      Scanner s = new Scanner(br);
      s.findWithinHorizon("nt=",0);
      int nt = s.nextInt();
      s.findWithinHorizon("dt=",0);
      double dt = s.nextDouble();
      s.findWithinHorizon("ft=",0);
      double ft = s.nextDouble();
      s.findWithinHorizon("x=",0);
      s.nextLine();
      float[] x = new float[nt];
      for (int it=0; it<nt; ++it) {
        for (int skip=1; skip<column; ++skip)
          s.next();
        x[it] = s.nextFloat();
        s.nextLine();
      }
      s.close();
      return new Sequence(nt,dt,ft,x);
    } catch (InputMismatchException ime) {
      throw new RuntimeException("Unknown format of file: "+fileName);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Not used, yet.

  /*
  scanf()               Regular Expression 
  %c 	                .
  %5c 	                .{5}
  %d 	                [-+]?\d+
  %e, %E, %f, %g 	[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?
  %i 	                [-+]?(0[xX][\dA-Fa-f]+|0[0-7]*|\d+)
  %o 	                0[0-7]*
  %s 	                \S+
  %u 	                \d+
  %x, %X 	        0[xX][\dA-Fa-f]+
  */
  // quoted string: "[^"\\\r\n]*(\\.[^"\\\r\n]*)*"

  private static class Named {
    public String name;
    public boolean found;
    public Named(String name) {
      this(name,"\\s");
    }
    protected Named(String name, String pattern) {
      this.name = name;
      _pattern = Pattern.compile(name+"=("+pattern+")");
    }
    protected String findString(String s) {
      Matcher m = _pattern.matcher(s);
      if (!m.find())
        return null;
      found = true;
      return m.group(1);
    }
    private Pattern _pattern;

    private static class Integer extends Named {
      public int value;
      public Integer(String name, int value) {
        super(name,"([-+]?\\d+)");
      }
      public boolean find(String s) {
        String v = super.findString(s);
        if (v==null)
          return false;
        value = java.lang.Integer.parseInt(v);
        found = true;
        return true;
      }
    }

    private static class Double extends Named {
      public double value;
      public Double(String name, double value) {
        super(name,"[-+]?(\\d+(\\.\\d*)?|\\d*\\.\\d+)([eE][-+]?\\d+)?");
      }
      public boolean find(String s) {
        String v = super.findString(s);
        if (v==null)
          return false;
        value = java.lang.Double.parseDouble(v);
        found = true;
        return true;
      }
    }
  }
}
