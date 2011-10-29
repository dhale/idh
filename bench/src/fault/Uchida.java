package fault;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Uchida's approximation to 2D warping.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.09.28
 */
public class Uchida {

  public static float[][][] getCosts(
    int lmin, int lmax, float[][] f, float[][] g) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    int nl = 1+lmax-lmin;
    float[][][] c = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      int lag = lmin+il;
      float ulag = lag;
      for (int i2=0; i2<n2; ++i2) {
        float[] f2 = f[i2];
        float[] g2 = g[i2];
        float[] cl2 = c[il][i2];
        int i1lo = min(n1-1,max(0,-lag));
        int i1hi = max(0,min(n1,n1-lag));
        float clo = cost(f2[i1lo],g2[0]);
        float chi = cost(f2[i1hi-1],g2[n1-1]);
        for (int i1=0; i1<i1lo; ++i1)
          cl2[i1] = clo;
        for (int i1=i1lo; i1<i1hi; ++i1)
          cl2[i1] = cost(f2[i1],g2[i1+lag]);
        for (int i1=i1hi; i1<n1; ++i1)
          cl2[i1] = chi;
      }
    }
    return c;
  }

  public static float[][] getShifts(
    int m2, int lmin, int lmax, float[][][] c) 
  {
    int n1 = c[0][0].length;
    int n2 = c[0].length;
    int nl = c.length;
    Transitions t = new Transitions(n1,m2,nl);
    float[][][] cj = new float[nl][m2][];
    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      int k2lo = max(0,min(n2-m2,i2-m2/2));
      int k2hi = k2lo+m2;
      for (int il=0; il<nl; ++il) {
        for (int k2=k2lo,j2=0; k2<k2hi; ++k2,++j2) {
          cj[il][j2] = c[il][k2];
        }
      }
      float[][][][] dj = getDistances(t,cj);
      float[][] sj = getShifts(lmin,lmax,t,dj);
      for (int k2=k2lo,j2=0; k2<k2hi; ++k2,++j2) {
        copy(n1,sj[j2],s[k2]);
      }
    }
    return s;
  }

  private static float[][][][] getDistances(Transitions t, float[][][] c) {
    int n1 = t.n1();
    int n2 = t.n2();
    int nk = t.nk();
    int nl = t.nl();
    float[][][][] d = new float[n1][n2][nk][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        int j1 = t.j1(i1,i2);
        int j2 = t.j2(i1,i2);
        for (int ik=0; ik<nk; ++ik) {
          for (int il=0; il<nl; ++il) {
            float dmin = 0.0f;
            int nj = t.nj(i2,ik,il);
            if ((i1>0 || i2>0) && nj>0) {
              dmin = Float.MAX_VALUE;
              for (int j=0; j<nj; ++j) {
                int jk = t.jk(i2,ik,il,j);
                int jl = t.jl(i2,ik,il,j);
                float dj = d[j1][j2][jk][jl];
                if (dj<dmin)
                  dmin = dj;
              }
            }
            d[i1][i2][ik][il] = dmin+c[il][i2][i1];
          }
        }
      }
    }
    return d;
  }

  private static float[][] getShifts(
    int lmin, int lmax, Transitions t, float[][][][] d) 
  {
    int n1 = t.n1();
    int n2 = t.n2();
    int nk = t.nk();
    int nl = t.nl();

    // Shifts correspond to lags to be found while backtracking.
    float[][] s = new float[n2][n1];

    // For sample (n1-1,n2-1), find state (ik,il) with smallest distance.
    float dmin = Float.MAX_VALUE;
    int ik = -1;
    int il = -1;
    for (int jk=0; jk<nk; ++jk) {
      for (int jl=0; jl<nl; ++jl) {
        float dj = d[n1-1][n2-1][jk][jl];
        if (dj<dmin && t.nj(n2-1,jk,jl)>0) {
          dmin = dj;
          ik = jk;
          il = jl;
        }
      }
    }
    assert ik>=0:"ik>=0";
    assert il>=0:"il>=0";

    // For all sample indices (i1,i2) in backtrack reverse order, ...
    for (int i1=n1-1; i1>=0; --i1) {
      for (int i2=n2-1; i2>=0; --i2) {

        // Shift at (i1,i2) corresponds to lag index il in state (ik,il).
        s[i2][i1] = lmin+il;

        // Indices (j1,j2) of next sample in backtrack.
        int j1 = t.j1(i1,i2);
        int j2 = t.j2(i1,i2);
        if (j1<0) 
          break;

        // For sample (j1,j2), find state (ik,il) with smallest distance.
        dmin = Float.MAX_VALUE;
        int jkmin = -1;
        int jlmin = -1;
        int nj = t.nj(i2,ik,il);
        for (int j=0; j<nj; ++j) {
          int jk = t.jk(i2,ik,il,j);
          int jl = t.jl(i2,ik,il,j);
          float dj = d[j1][j2][jk][jl];
          if (dj<dmin) {
            dmin = dj;
            jkmin = jk;
            jlmin = jl;
          }
        }
        ik = jkmin;
        il = jlmin;
      }
    }
    return s;
  }

  private static float cost(float x, float y) {
    float d = x-y;
    return d*d; // square of difference
    //return d>=0.0f?d:-d; // absolute value of difference
  }

  private static class Transitions {

    private Transitions(int n1, int n2, int nl) {

      // Total number of states is nk*nl, where nk = 3^(n2-1).
      int nk = 1;
      for (int i2=1; i2<n2; ++i2)
        nk *= 3;
      //System.out.println("table memory = "+(2*n2*nk*nl*3*4));

      // Compute lag[n2][nk][nl] for nk*nl states, each with n2 lags.
      int[][][] lag = new int[n2][nk][nl];
      for (int ik=0; ik<nk; ++ik) {
        for (int il=0; il<nl; ++il) {
          lag[0][ik][il] = il;
        }
      }
      for (int i2=1,mk=1; i2<n2; ++i2,mk*=3) {
        for (int ik=0; ik<nk; ++ik) {
          for (int il=0; il<nl; ++il) {
            int lagi = lag[i2-1][ik][il]+(ik/mk)%3-1;
            lag[i2][ik][il] = (0<=lagi && lagi<nl)?lagi:-1;
          }
        }
      }
      /*
      // Print lags.
      for (int ik=0; ik<nk; ++ik) {
        for (int il=0; il<nl; ++il) {
          System.out.print("ik="+ik+" il="+il+" lag=");
          for (int i2=0; i2<n2; ++i2) {
            System.out.print(lag[i2][ik][il]);
            if (i2==n2-1) {
              System.out.println();
            } else {
              System.out.print(",");
            }
          }
        }
      }
      */

      // For all states (ik,il), compute indices (jk,jl) of ancestors.
      // The indices (jk,jl) depend on ik, il and i2. The maximum 
      // number of ancestors for any state (ik,il) is three.
      int[][][][] jk3 = new int[n2][nk][nl][3];
      int[][][][] jl3 = new int[n2][nk][nl][3];
      for (int i2=0; i2<n2; ++i2) {
        for (int ik=0; ik<nk; ++ik) {
          for (int il=0; il<nl; ++il) {
            for (int j=0; j<3; ++j) {
              jk3[i2][ik][il][j] = -1; // initialize to -1, so we know
              jl3[i2][ik][il][j] = -1; // which (jk,jl) are set below
            }
            int j = 0;
            int lagii = lag[i2][ik][il];
            if (lagii<0) continue;
            for (int jk=0; jk<nk; ++jk) {
              for (int jl=0; jl<nl; ++jl) {
                int lagij = lag[i2][jk][jl];
                if (lagij<0) continue;
                boolean transition = abs(lagii-lagij)<=1; // transition ok?
                for (int j2=0; transition && j2<n2; ++j2) {
                  if (j2==i2) continue;
                  int lagji = lag[j2][ik][il];
                  int lagjj = lag[j2][jk][jl];
                  transition = lagji>=0 && lagjj>=0 && lagji==lagjj;
                }
                if (transition) { // if transition is ok, ...
                  //System.out.println(
                  //  "i2="+i2+" ik="+ik+" il="+il+" jk="+jk+" jl="+jl+" j="+j);
                  jk3[i2][ik][il][j] = jk;
                  jl3[i2][ik][il][j] = jl;
                  ++j;
                }
              }
            }
          }
        }
      }
      _n1 = n1;
      _n2 = n2;
      _nk = nk;
      _nl = nl;
      _jk3 = jk3;
      _jl3 = jl3;
    }

    // Returns index j1 of sample that precedes sample with indices (i1,i2).
    int j1(int i1, int i2) {
      return (i2>0)?i1:i1-1;
    }

    // Returns index j2 of sample that precedes sample with indices (i1,i2).
    int j2(int i1, int i2) {
      return (i2>0)?i2-1:_n2-1;
    }

    // Returns number of samples in 1st dimension.
    int n1() {
      return _n1;
    }

    // Returns number of samples in 2nd dimension.
    int n2() {
      return _n2;
    }

    // Returns factor 3^(n2-1) in number of states nk*nl.
    int nk() {
      return _nk;
    }

    // Returns number of lags in number of states nk*nl.
    int nl() {
      return _nl;
    }

    // Returns number of states j preceding the specified state i.
    int nj(int i2, int ik, int il) {
      int nj = 0;
      int[] jk = _jk3[i2][ik][il];
      while (nj<3 && jk[nj]>=0)
        ++nj;
      return nj;
    }

    // Returns index jk for j'th state that precedes the specified state i.
    int jk(int i2, int ik, int il, int j) {
      return _jk3[i2][ik][il][j];
    }

    // Returns index jl for j'th state that precedes the specified state i.
    int jl(int i2, int ik, int il, int j) {
      return _jl3[i2][ik][il][j];
    }

    private int _n1,_n2,_nk,_nl;
    private int[][][][] _jk3;
    private int[][][][] _jl3;
  }
}
