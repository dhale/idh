package tp;

public class DoubleList {
  public int n;
  public double[] a = new double[1024];
  public void add(double d) {
    if (n==a.length) {
      double[] t = new double[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = d;
  }
  public double[] trim() {
    if (n==0)
      return null;
    double[] t = new double[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
}
