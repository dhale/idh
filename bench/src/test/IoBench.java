package test;

import java.io.File;
import java.io.IOException;

import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.util.Stopwatch;
import static edu.mines.jtk.util.ArrayMath.*;

public class IoBench {
  public static void main(String[] args) {
    //benchCopy();
    benchReadWrite(23);
    //benchReadWrite(32);
  }

  private static void benchReadWrite(int order) {
    String file;
    for(int iter=0; iter<3; ++iter) {
      file = dira+"junk.dat"; 
      remove(file);
      makeFile(file);
      benchReadWrite(file,order);
      file = dirb+"junk.dat"; 
      remove(file);
      makeFile(file);
      benchReadWrite(file,order);
    }
  }
  private static void benchReadWrite(String file, int order) {
    trace("benchReadWrite: file="+file+" order="+order);
    Stopwatch sw = new Stopwatch();
    sw.start();
    readWrite(file,order);
    sw.stop();
    int rate = (int)(nbytes*1.0e-6/sw.time());
    trace("time = "+sw.time()+" s rate = "+rate+" MB/s");
  }

  private static void benchCopy() {
    String afile,bfile;
    for(int iter=0; iter<3; ++iter) {
      /*
      afile = dira+"junka.dat"; 
      bfile = dira+"junkb.dat"; 
      remove(afile);
      remove(bfile);
      makeFile(afile);
      benchCopy(afile,bfile);
      */
      afile = dirb+"junk.dat"; 
      bfile = dirc+"junk.dat"; 
      //remove(afile);
      //remove(bfile);
      //makeFile(afile);
      benchCopy(afile,bfile);
    }
  }
  private static void benchCopy(String afile, String bfile) {
    trace("benchCopy: afile="+afile+" bfile="+bfile);
    Stopwatch sw = new Stopwatch();
    sw.start();
    copy(afile,bfile);
    sw.stop();
    int rate = (int)(nbytes*1.0e-6/sw.time());
    trace("time = "+sw.time()+" s rate = "+rate+" MB/s");
  }

  private static void benchRandom() {
    // TODO
  }

  private static final String dira = "/data/tmp/";
  private static final String dirb = "/data/tmp/";
  private static final String dirc = "/data/tmp/";
  private static final int n1 = 1000;
  private static final int n2 = 100;
  private static final int n3 = 100;
  private static final long nbytes = n1*n2*n3*4L;

  private static void makeFile(String name) {
    trace("makeFile: begin");
    try {
      ArrayFile af = new ArrayFile(name,"rw");
      float[] a = randfloat(n1);
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          af.writeFloats(a);
        }
      }
      af.close();
    } catch(IOException e) {
      throw new RuntimeException(e);
    }
    trace("done");
  }

  private static void copy(String afile, String bfile) {
    try {
      ArrayFile afa = new ArrayFile(afile,"r");
      ArrayFile afb = new ArrayFile(bfile,"rw");
      float[] a = new float[n1];
      for (int i3=0; i3<n3; ++i3) {
        if (i3%100==0)
          trace("copy: i3="+i3);
        for (int i2=0; i2<n2; ++i2) {
          afa.readFloats(a);
          afb.writeFloats(a);
        }
      }
      afa.close();
      afb.close();
    } catch(IOException e) {
      throw new RuntimeException(e);
    }
  }

  private static void readWrite(String file, int order) {
    try {
      ArrayFile af = new ArrayFile(file,"rw");
      float[] a = new float[n1];
      if (order==23) {
        for (int i3=0; i3<n3; ++i3) {
          for (int i2=0; i2<n2; ++i2) {
            af.seek(4*n1*i2*i3);
            af.readFloats(a);
            af.seek(4*n1*i2*i3);
            af.writeFloats(a);
          }
        }
      } else if (order==32) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i3=0; i3<n3; ++i3) {
            af.seek(4*n1*i2*i3);
            af.readFloats(a);
            af.seek(4*n1*i2*i3);
            af.writeFloats(a);
          }
        }
      }
      af.close();
    } catch(IOException e) {
      throw new RuntimeException(e);
    }
  }

  private static void remove(String fileName) {
    File file = new File(fileName);
    if (file.exists())
      file.delete();
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
