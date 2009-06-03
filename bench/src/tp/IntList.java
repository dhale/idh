package tp;

public class IntList {
  public int n;
  public int[] a = new int[1024];
  public void add(int i) {
    if (n==a.length) {
      int[] t = new int[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = i;
  }
  public int[] trim() {
    if (n==0)
      return null;
    int[] t = new int[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
}
