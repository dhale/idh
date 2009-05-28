package tp;

public class FloatList {
  public int n;
  public float[] a = new float[1024];
  public void add(float f) {
    if (n==a.length) {
      float[] t = new float[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = f;
  }
  public float[] trim() {
    if (n==0)
      return null;
    float[] t = new float[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
}
