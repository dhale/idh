package fault;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Utilities.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.09.01
 */
public class Util {

  private static float cost(float x, float y) {
    float d = x-y;
    return d*d; // square of difference
    //return d>=0.0f?d:-d; // absolute value of difference
  }

  public static float[][][] dtwCost1(
    int lmin, int lmax, double sigma, float[][] f, float[][] g) 
  {
    int n1 = f[0].length;
    int n2 = f.length;
    int nl = 1+lmax-lmin;
    float[][][] c = new float[nl][n2][n1];
    RecursiveExponentialFilter ref = null;
    if (sigma>0.0)
      ref = new RecursiveExponentialFilter(sigma);
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
      if (ref!=null)
        ref.apply2(c[il],c[il]);
    }
    return c;
  }

  public static float[][][] dtwCost12(float[][][] c, int limit) {
    float[][][] d2f = dtwDistance2F(c,limit);
    float[][][] d2r = dtwDistance2R(c,limit);
    return sub(add(d2f,d2r),c);
  }

  public static float[][][] dtwDistance2F(float[][][] c, int limit) {
    int n1 = c[0][0].length;
    int n2 = c[0].length;
    int nl = c.length;
    int nlm1 = nl-1;
    float[][][] d = new float[nl][n2][n1];
    int i2 = 0;
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        d[il][i2][i1] = c[il][i2][i1];
      }
    }
    if (limit>1) {
      i2 = 1;
      for (int il=0; il<nl; ++il) {
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        for (int i1=0; i1<n1; ++i1) {
          d[il][i2][i1] = c[il][i2][i1]+min(d[im][i2-1][i1],
                                            d[il][i2-1][i1],
                                            d[ip][i2-1][i1]);
        }
      }
    }
    if (limit>2) {
      i2 = 2;
      for (int il=0; il<nl; ++il) {
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        for (int i1=0; i1<n1; ++i1) {
          d[il][i2][i1] = c[il][i2][i1]+min(d[im][i2-2][i1]+c[im][i2-1][i1],
                                            d[il][i2-1][i1],
                                            d[ip][i2-2][i1]+c[ip][i2-1][i1]);
        }
      }
    }
    if (limit==1) {
      for (i2=1; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          float[] cii = c[il][i2];
          float[] dii = d[il][i2];
          float[] dmm = d[im][i2-1];
          float[] dim = d[il][i2-1];
          float[] dpm = d[ip][i2-1];
          for (int i1=0; i1<n1; ++i1) {
            dii[i1] = cii[i1]+min(dmm[i1],
                                  dim[i1],
                                  dpm[i1]);
          }
        }
      }
    } else if (limit==2) {
      for (i2=2; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          float[] cii = c[il][i2  ];
          float[] cmm = c[im][i2-1];
          float[] cpm = c[ip][i2-1];
          float[] dii = d[il][i2  ];
          float[] dmm = d[im][i2-2];
          float[] dim = d[il][i2-1];
          float[] dpm = d[ip][i2-2];
          for (int i1=0; i1<n1; ++i1) {
            dii[i1] = cii[i1]+min(dmm[i1]+cmm[i1],
                                  dim[i1],
                                  dpm[i1]+cpm[i1]);
          }
        }
      }
    } else if (limit==3) {
      for (i2=3; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          float[] cii =  c[il][i2  ];
          float[] cmm1 = c[im][i2-1];
          float[] cpm1 = c[ip][i2-1];
          float[] cmm2 = c[im][i2-2];
          float[] cpm2 = c[ip][i2-2];
          float[] dii  = d[il][i2  ];
          float[] dmm3 = d[im][i2-3];
          float[] dim1 = d[il][i2-1];
          float[] dpm3 = d[ip][i2-3];
          for (int i1=0; i1<n1; ++i1) {
            dii[i1] = cii[i1]+min(dmm3[i1]+cmm1[i1]+cmm2[i1],
                                  dim1[i1],
                                  dpm3[i1]+cpm1[i1]+cpm2[i1]);
          }
        }
      }
    }
    return d;
  }

  public static int[][][][][] getTransitions(int nl, int n2) {
    int nk = 1;
    for (int i2=1; i2<n2; ++i2)
      nk *= 3;

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

    // Print lags.
    /*
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
    System.out.println("table memory = "+(2*n2*nk*nl*3*4));
    return new int[][][][][]{jk3,jl3};
  }

  public static float[][][] dtwDistance2R(float[][][] c, int limit) {
    int n1 = c[0][0].length;
    int n2 = c[0].length;
    int nl = c.length;
    int nlm1 = nl-1;
    float[][][] d = new float[nl][n2][n1];
    int i2 = n2-1;
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        d[il][i2][i1] = c[il][i2][i1];
      }
    }
    if (limit>1) {
      i2 = n2-2;
      for (int il=0; il<nl; ++il) {
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        for (int i1=0; i1<n1; ++i1) {
          d[il][i2][i1] = c[il][i2][i1]+min(d[im][i2+1][i1],
                                            d[il][i2+1][i1],
                                            d[ip][i2+1][i1]);
        }
      }
    }
    if (limit>2) {
      i2 = n2-3;
      for (int il=0; il<nl; ++il) {
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        for (int i1=0; i1<n1; ++i1) {
          d[il][i2][i1] = c[il][i2][i1]+min(d[im][i2+2][i1]+c[im][i2+1][i1],
                                            d[il][i2+1][i1],
                                            d[ip][i2+2][i1]+c[ip][i2+1][i1]);
        }
      }
    }
    if (limit==1) {
      for (i2=n2-2; i2>=0; --i2) {
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          float[] cii = c[il][i2];
          float[] dii = d[il][i2];
          float[] dmp = d[im][i2+1];
          float[] dip = d[il][i2+1];
          float[] dpp = d[ip][i2+1];
          for (int i1=0; i1<n1; ++i1) {
            dii[i1] = cii[i1]+min(dmp[i1],
                                  dip[i1],
                                  dpp[i1]);
          }
        }
      }
    } else if (limit==2) {
      for (i2=n2-3; i2>=0; --i2) {
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          float[] cii = c[il][i2  ];
          float[] cmp = c[im][i2+1];
          float[] cpp = c[ip][i2+1];
          float[] dii = d[il][i2  ];
          float[] dmp = d[im][i2+2];
          float[] dip = d[il][i2+1];
          float[] dpp = d[ip][i2+2];
          for (int i1=0; i1<n1; ++i1) {
            dii[i1] = cii[i1]+min(dmp[i1]+cmp[i1],
                                  dip[i1],
                                  dpp[i1]+cpp[i1]);
          }
        }
      }
    } else if (limit==3) {
      for (i2=n2-4; i2>=0; --i2) {
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          float[] cii =  c[il][i2  ];
          float[] cmp1 = c[im][i2+1];
          float[] cpp1 = c[ip][i2+1];
          float[] cmp2 = c[im][i2+2];
          float[] cpp2 = c[ip][i2+2];
          float[] dii  = d[il][i2  ];
          float[] dmp3 = d[im][i2+3];
          float[] dip1 = d[il][i2+1];
          float[] dpp3 = d[ip][i2+3];
          for (int i1=0; i1<n1; ++i1) {
            dii[i1] = cii[i1]+min(dmp3[i1]+cmp1[i1]+cmp2[i1],
                                  dip1[i1],
                                  dpp3[i1]+cpp1[i1]+cpp2[i1]);
          }
        }
      }
    }
    return d;
  }

  public static float[][][] dtwCostU1(
    int lmin, int lmax, float[][][] c,
    float uweight, float[][] u) 
  {
    int n1 = c[0][0].length;
    int n2 = c[0].length;
    int nl = c.length;

    // Update a copy of the specified costs.
    c = copy(c);

    // Scale the specified weight by the average cost.
    uweight *= sum(c)/nl/n2/n1;
    System.out.println("uweight="+uweight);

    // Add costs for derivatives of u along 2nd dimension.
    for (int il=0; il<nl; ++il) {
      int lag = lmin+il;
      float ulag = lag;
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(0,i2-1);
        int i2p = min(n2-1,i2+1);
        float[] u2m = u[i2m];
        float[] u2p = u[i2p];
        float[] cl2  = c[il][i2];
        for (int i1=0; i1<n1; ++i1) {
          float dum = ulag-u2m[i1];
          float dup = ulag-u2p[i1];
          cl2[i1] += uweight*(abs(dum)+abs(dup));
        }
      }
    }
    return c;
  }

  public static float[][] dtwShifts1(
    int lmin, int lmax, float[][][] c, int limit) 
  {
    int n1 = c[0][0].length;
    int n2 = c[0].length;
    int nl = c.length;
    float[][] u = new float[n2][];
    float[][] c2 = new float[nl][];
    for (int i2=0; i2<n2; ++i2) {
      for (int il=0; il<nl; ++il)
        c2[il] = c[il][i2];
      float[][] d2 = dtwDistance(c2,limit);
      u[i2] = dtwShifts(lmin,lmax,d2);
    }
    return u;
  }

  public static float[][] dtwCost(int lmin, int lmax, float[] f, float[] g) {
    int n1 = f.length;
    int nl = 1+lmax-lmin;
    float[][] c = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      int lag = lmin+il;
      float[] cl = c[il];
      int i1lo = min(n1-1,max(0,-lag));
      int i1hi = max(0,min(n1,n1-lag));
      float clo = cost(f[i1lo],g[0]);
      float chi = cost(f[i1hi-1],g[n1-1]);
      for (int i1=0; i1<i1lo; ++i1)
        cl[i1] = clo;
      for (int i1=i1lo; i1<i1hi; ++i1)
        cl[i1] = cost(f[i1],g[i1+lag]);
      for (int i1=i1hi; i1<n1; ++i1)
        cl[i1] = chi;
    }
    return c;
  }

  public static float[][] dtwDistance(float[][] c, int limit) {
    int n1 = c[0].length;
    int nl = c.length;
    int nlm1 = nl-1;
    float[][] d = new float[nl][n1];
    for (int il=0; il<nl; ++il) // i1 = 0
      d[il][0] = c[il][0];
    if (limit==1) {
      for (int i1=1; i1<n1; ++i1) { // i1 = 1 to n1-1
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          d[il][i1] = c[il][i1]+min(d[im][i1-1],
                                    d[il][i1-1],
                                    d[ip][i1-1]);
        }
      }
    } else if (limit==2) {
      for (int il=0; il<nl; ++il) { // i1 = 1
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        d[il][1] = c[il][1]+min(d[im][0],
                                d[il][0],
                                d[ip][0]);
      }
      for (int i1=2; i1<n1; ++i1) { // i1 = 2 to n1-1
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          d[il][i1] = c[il][i1]+min(d[im][i1-2]+c[im][i1-1],
                                    d[il][i1-1],
                                    d[ip][i1-2]+c[ip][i1-1]);
        }
      }
    } else if (limit==3) {
      for (int il=0; il<nl; ++il) { // i1 = 1
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        d[il][1] = c[il][1]+min(d[im][0],
                                d[il][0],
                                d[ip][0]);
      }
      for (int il=0; il<nl; ++il) { // i1 = 2
        int im = max(0,il-1);
        int ip = min(nlm1,il+1);
        d[il][2] = c[il][2]+min(d[im][0]+c[im][1],
                                d[il][1],
                                d[ip][0]+c[ip][1]);
      }
      for (int i1=3; i1<n1; ++i1) { // i1 = 3 to n1-1
        for (int il=0; il<nl; ++il) {
          int im = max(0,il-1);
          int ip = min(nlm1,il+1);
          d[il][i1] = c[il][i1]+min(d[im][i1-3]+c[im][i1-1]+c[im][i1-2],
                                    d[il][i1-1],
                                    d[ip][i1-3]+c[ip][i1-1]+c[ip][i1-2]);
        }
      }
    }
    return d;
  }

  public static float[] dtwShifts(int lmin, int lmax, float[][] d) {
    int n1 = d[0].length;
    int nl = 1+lmax-lmin;
    float[] u = new float[n1];
    int lm = 0;
    float dm = d[lm][n1-1];
    for (int il=1; il<nl; ++il) {
      if (d[il][n1-1]<dm) {
        dm = d[il][n1-1];
        lm = il;
      }
    }
    u[n1-1] = lm+lmin;
    for (int i1=n1-2; i1>=0; --i1) {
      int la = lm>0?lm-1:lm;
      int lb = lm;
      int lc = lm<nl-1?lm+1:lm;
      float da = d[la][i1];
      float db = d[lb][i1];
      float dc = d[lc][i1];
      dm = min(da,db,dc);
      if (db==dm) {
        lm = lb;
      } else if (da==dm) {
        lm = la;
      } else {
        lm = lc;
      }
      u[i1] = lm+lmin;
    }
    return u;
  }

  public static float[][] dtwWarpMatrix(float[][] d) {
    int n = d.length;
    int i = n-1;
    int j = n-1;
    /*
    float dij = d[i][j];
    for (int k=0; k<n; ++k) {
      if (d[k][n-1]<dij) {
        i = k; j = n-1;
        dij = d[i][j];
      }
    }
    for (int k=0; k<n; ++k) {
      if (d[n-1][k]<dij) {
        i = n-1; j = k;
        dij = d[i][j];
      }
    }
    */
    ArrayList<Integer> il = new ArrayList<Integer>(n);
    ArrayList<Integer> jl = new ArrayList<Integer>(n);
    il.add(i);
    jl.add(j);
    while (i>0 || j>0) {
      if (i==0) {
        --j;
      } else if (j==0) {
        --i;
      } else {
        float dmm = d[i-1][j-1];
        float d0m = d[i  ][j-1];
        float dm0 = d[i-1][j  ];
        float dmin = min(dmm,d0m,dm0);
        if (dmin==dmm) {
          --i;
          --j;
        } else if (dmin==dm0) {
          --i;
        } else {
          --j;
        }
      }
      il.add(i);
      jl.add(j);
    }
    n = il.size();
    float[] x = new float[n];
    float[] y = new float[n];
    for (int k=0; k<n; ++k) {
      x[k] = il.get(k);
      y[k] = jl.get(k);
    }
    return new float[][]{x,y};
  }

  public static float[][] dtwDistanceMatrix(float[][] c) {
    int n = c.length;
    float[][] d = fillfloat(Float.MAX_VALUE,n,n);
    d[0][0] = 0.0f;
    for (int i=1; i<n; ++i)
      d[i][0] = d[i-1][0]+c[i][0];
    for (int j=1; j<n; ++j)
      d[0][j] = d[0][j-1]+c[0][j];
    for (int i=1; i<n; ++i) {
      for (int j=1; j<n; ++j) {
        d[i][j] = c[i][j]+min(d[i-1][j-1],d[i][j-1],d[i-1][j]);
      }
    }
    return d;
  }

  public static float[] findShiftsSmooth(
    double s, double sigma, int lmin, int lmax, float[][] c) 
  {
    int n1 = c[0].length;
    float[] g = new float[n1];
    float[] h = new float[n1];
    float[] u = new float[n1];
    float[] v = new float[n1];
    float[] d = new float[n1];
    //u = mul(4.0f,sin(rampfloat(0.0f,2.0f*FLT_PI/(n1-1),n1)));
    //v = copy(u);
    VecArrayFloat1 vg = new VecArrayFloat1(g);
    VecArrayFloat1 vd = new VecArrayFloat1(d);
    Sampling sl = new Sampling(1+lmax-lmin,1.0,lmin);
    CgSolver cs = new CgSolver(0.01,4);
    for (int iter=0; iter<64; ++iter) {
      // Solve (sI-(1-s)S'C''S)dv = (1-s)S'c'-sv = g
      gH(sl,c,u,g,h); // g = c', h = c''
      RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
      ref.apply(g,g); // g = S'c'
      float sg = (float)s;
      float tg = 1.0f-sg;
      for (int i1=0; i1<n1; ++i1)
        g[i1] = tg*g[i1]-sg*v[i1]; // g = (1-s)S'c'-sv
      A a = new A(s,h,ref);
      zero(d);
      cs.solve(a,vg,vd);
      add(d,v,v);
      ref.apply(v,u);
      //ref.apply(d,d);
      //add(d,u,u);
    }
    return u;
  }
  private static void gH(
    Sampling sl, float[][] c, float[] u, float[] g, float[] h) 
  {
    int n1 = c[0].length;
    int nl = c.length;
    for (int i1=0; i1<n1; ++i1) {
      float ui = u[i1];
      int il = sl.indexOfNearest(ui);
      if (il==0   ) il = 1;
      if (il==nl-1) il = nl-2;
      double du = ui-sl.getValue(il);
      double ai = c[il-1][i1];
      double bi = c[il  ][i1];
      double ci = c[il+1][i1];
      //double q0 = bi;
      double q1 = 0.5*(ci-ai);
      double q2 = 0.5*(ci+ai)-bi;
      //double ci = q0+du*(q1+du*q2);
      g[i1] = (float)(q1+2.0*du*q2);
      h[i1] = (float)(2.0*q2);
    }
    fill(min(h),h);
    //System.out.println("max g="+max(g)+" min h="+min(h));
    //System.out.println("g="+g[3*n1/4]+" h="+h[3*n1/4]);
    // c = 1-(x-4)^2, c' = -2*(x-4), c'' = -2
  }
  private static class A implements CgSolver.A {
    A(double s, float[] h, RecursiveExponentialFilter ref) {
      _s = s;
      _h = h;
      _ref = ref;
    }
    public void apply(Vec vx, Vec vy) {
      float[] x = ((VecArrayFloat1)vx).getArray();
      float[] y = ((VecArrayFloat1)vy).getArray();
      int n = x.length;
      //(sI-(1-s)S'C''S)dv = (1-s)S'c'-sv
      //(sI-(1-s)S'C''S)x = y
      _ref.apply(x,y);
      mul(_h,y,y);
      _ref.apply(y,y);
      float s = (float)_s;
      float t = 1.0f-s;
      for (int i=0; i<n; ++i)
        y[i] = s*x[i]-t*y[i];
    }
    private double _s;
    private float[] _h;
    private RecursiveExponentialFilter _ref;
  }
}
/*
E = s(u'R'Ru)/2 - (1-s)c(u)
g = sR'Ru - (1-s)c' = gradient
H = sR'R - (1-s)C'' = Hessian
(sR'R-(1-s)C'')du = (1-s)c'-sR'Ru => s = 1 -> du = -u -> u = u+du = 0
(sI-(1-s)S'C''S)dv = (1-s)S'c'-sv => s = 1 -> dv = -v -> v = v+dv = 0
*/
