package test;

import edu.mines.jtk.util.Stopwatch;

public class Run {
  public static void main(String[] args) {
    Stopwatch s = new Stopwatch();
    s.start();
    while (s.time()<1)
      ;
    s.stop();
    System.out.println("time = "+s.time());
  }
}
