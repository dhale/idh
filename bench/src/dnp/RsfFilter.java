/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dnp;

import java.io.*;
import java.util.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.*;

// for testing only
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Filters arrays of floats using Madagascar (RSF) programs.
 * Applies this filter through the following sequence of steps:
 * <ol><li> make temporary RSF input and output files
 * </li><li> write from the input array into the RSF input file
 * </li><li> associate the RSF input file with stdin
 * </li><li> associate the RSF output file with stdout
 * </li><li> run the RSF program; reads from stdin and writes to stdout
 * </li><li> read from the RSF output file into the output array
 * </li><li> delete the temporary files
 * </li></ol>
 * The RSF program runs in the environment in which the JVM is running.
 * That environment includes the current directory, path, etc.
 * @author Dave Hale and Elias Arias, Colorado School of Mines
 * @version 2014.01.25
 */
public class RsfFilter {

  /**
   * Constructs a filter for specified RSF program and arguments.
   * The program name and arguments are specified as they might be in a shell
   * command. The arguments can be specified in an array of strings or with a
   * comma-separated sequence of strings.
   * @param name the RSF program name.
   * @param args the program arguments.
   */
  public RsfFilter(String name, String... args) {
     _argList = new ArrayList<String>();
     _argList.add(name);
     for (String arg:args) 
       _argList.add(arg);
  }

  /**
   * Sets the directory in which temporary RSF files are created.
   * The default directory is that set in the property java.io.tmpdir.
   * @param tempDirName name of the directory used for temporary files.
   */
  public void setTempDir(String tempDirName) {
    if (tempDirName==null) {
      _tempDir = null;
    } else {
      _tempDir = new File(tempDirName);
    }
  }

  /**
   * Specifies whether to keep (or delete) temporary RSF files.
   * This method is useful when trying to determine why an RSF program may
   * have failed. The default is false, so that all temporary files are
   * deleted. 
   * @param keep true, to keep temporary files; false, to delete them.
   */
  public void setTempKeep(boolean keep) {
    _tempKeep = keep;
  }

  /**
   * Applies this filter to an array with default samplings.
   * @param x the input array.
   * @return the output array.
   */
  public float[][] apply(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    apply(x,y);
    return y;
  }

  /**
   * Applies this filter for arrays with default samplings.
   * @param x the input array.
   * @param y the output array.
   */
  public void apply(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    apply(new Sampling(n1),new Sampling(n2),x,y);
  }

  /**
   * Applies this filter to an array with specified samplings.
   * @param s1 the sampling for the 1st dimension.
   * @param s2 the sampling for the 2nd dimension.
   * @param x the input array.
   * @return the output array.
   */
  public float[][] apply(Sampling s1, Sampling s2, float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    apply(s1,s2,x,y);
    return y;
  }

  /**
   * Applies this filter for arrays with specified samplings.
   * @param s1 the sampling for the 1st dimension.
   * @param s2 the sampling for the 2nd dimension.
   * @param x the input array.
   * @param y the output array.
   */
  public void apply(Sampling s1, Sampling s2, float[][] x, float[][] y) {
    apply(new Sampling[]{s1,s2},x,y);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private ArrayList<String> _argList;
  private File _tempDir;
  private boolean _tempKeep;

  private void apply(Sampling[] ss, Object x, Object y) {
    File rsfInFile = null;
    File rsfOutFile = null;
    File rsfErrFile = null;
    File binInFile = null;
    File binOutFile = null;

    // If any exceptions, rethrow them as RuntimeExceptions.
    try {

      // Make temporary RSF files.
      rsfInFile = File.createTempFile("rsfi",".rsf",_tempDir);
      rsfOutFile = File.createTempFile("rsfo",".rsf",_tempDir);
      rsfErrFile = File.createTempFile("rsfe",".txt",_tempDir);
      binInFile = new File(rsfInFile+"@");
      binOutFile = new File(rsfOutFile+"@");

      // Write floats to RSF input file.
      ArrayOutputStream aos = new ArrayOutputStream(binInFile);
      if (x instanceof float[]) {
        aos.writeFloats((float[])x);
      } else if (x instanceof float[][]) {
        aos.writeFloats((float[][])x);
      } else if (x instanceof float[][][]) {
        aos.writeFloats((float[][][])x);
      }
      aos.close();
      PrintWriter pw = new PrintWriter(rsfInFile);
      pw.println("in=\""+binInFile+"\"");
      for (int is=0; is<ss.length; ++is) {
        Sampling si = ss[is];
        int js = is+1;
        pw.print(" n"+js+"="+si.getCount());
        pw.print(" d"+js+"="+si.getDelta());
        pw.print(" o"+js+"="+si.getFirst());
        pw.println();
      }
      pw.println("esize=4 type=float data_format=\"xdr_float\"");
      pw.close();

      // Add the RSF output file to the command. Use of "out=" will
      // override the RSF datapath, if any, that may have been set.
      ArrayList<String> command = new ArrayList<String>(_argList);
      command.add("out="+binOutFile);

      // Build the process with redirected standard in,out,err.
      ProcessBuilder pb = new ProcessBuilder(command);
      pb.redirectInput(rsfInFile);
      pb.redirectOutput(rsfOutFile);
      pb.redirectError(rsfErrFile);

      // Start the process and wait for it to finish. If the RSF program fails
      // for any reason, print standard error, and throw a RuntimeException.
      try {
        int result = pb.start().waitFor();
        if (result!=0) {
          Scanner s = new Scanner(rsfErrFile);
          while (s.hasNextLine()) {
            String line = s.nextLine();
            System.err.println(line);
          }
          s.close();
          throw new RuntimeException("RSF command failed");
        }
      } catch (InterruptedException e) {
        throw new RuntimeException(e);
      }

      // Read floats from the RSF output file.
      ArrayInputStream ais = new ArrayInputStream(binOutFile);
      if (y instanceof float[]) {
        ais.readFloats((float[])y);
      } else if (y instanceof float[][]) {
        ais.readFloats((float[][])y);
      } else if (y instanceof float[][][]) {
        ais.readFloats((float[][][])y);
      }
      ais.close();

    } catch (Exception e) {
      throw new RuntimeException(e);
    } finally {
      if (!_tempKeep) {
        if (rsfInFile!=null) rsfInFile.delete();
        if (rsfOutFile!=null) rsfOutFile.delete();
        if (rsfErrFile!=null) rsfErrFile.delete();
        if (binInFile!=null) binInFile.delete();
        if (binOutFile!=null) binOutFile.delete();
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    int n1 = 5;
    int n2 = 6;
    float[][] x = cos(rampfloat(0.0f,0.1f,0.02f,n1,n2));
    RsfFilter rf = new RsfFilter("sfdip","pmin=-0.5","pmax=0.5");
    //rf.setTempDir(".");
    //rf.setTempKeep(true);
    float[][]y = rf.apply(x);
    dump(x);
    dump(y);
  }
}
