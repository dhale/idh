package fault;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/*
TODO: 
Implement constrained accumulate2 methods for specified i2c
Begin as for unconstrained methods:
  for all i2
    accumulate1
  for all i1
    accumulate2
choose i2 = i2c variance of u1 in x2 direction is small
for i2 = i2c, accumulateForward1 and findShiftsReverse1 to get u[i2c][i1]
for i2 = i2c+1, i2c+2, ...
  set e[i2][i1][il] = emax for all i1 outside of u[i2-1][i1] +- 1
  accumulateForward1 and findShiftsReverse1 to get u[i2][i1]
for i2 = i2c-1, i2c-2, ...
  set e[i2][i1][il] = emax for all i1 outside of u[i2-1][i1] +- 1
  accumulateForward1 and findShiftsReverse1 to get u[i2][i1]
Problem: these last loops over i2 do not enforce bstretch2 on du/d2!
*/

/**
 * Dynamic warping of sequences or images.
 * <p>
 * For sequences f and g, dynamic warping finds a 1D array of 
 * integer shifts u[i1] such that f[i1] ~ g[i1+u[i1]], subject to 
 * constraints |u[i1]-u[i1-1]|&lt;b1, where b1 is a bound on 
 * the local average slope of u[i1].
 * <p>
 * A positive slope u[i1]-u[i1-1] = 1 implies that g[i1] between 
 * indices i1-1 and i1 is a stretched version of f[i1] ~ g[i1+u[i1]].
 * (E.g., f[i1-1] ~ g[i1-1+u[i1-1]] and f[i1] ~ g[i1+u[i1-1]+1].)
 * Values in f for indices i1 and i1-1 are one sample apart. Those
 * same values in g are two samples apart, which implies stretching 
 * by 100%. Likewise, a negative slope u[i1]-u[i1-1] = -1 implies 
 * squeezing by 100%, such that both values f[i1] and f[i1-1]
 * correspond to g[i1+u1[i1-1]-1].
 * <p>
 * In practice, stretching or squeezing by as much as 100% is extreme.
 * Therefore, the constraint b1 on the slope of u[i1] may be smaller
 * than one. For example, if b1 = 0.5, then the absolute value of the 
 * local average slope of u[i1] is bounded by 0.5. This bound is 
 * complicated by the fact that the shifts u[i1] are integers. Again, 
 * for b1 = 0.5, a local average slope of 0.5 would correspond to 
 * u[i1-1]-u[i1-2] = 0 and u[i1]-u[i1-1] = 1.
 * <p>
 * For images f and g, dynamic warping finds a 2D array of integer 
 * shifts u[i2][i1] such that f[i2][i1] ~= g[i2][i1+u[i2][i1]], subject 
 * to |u[i2][i1]-u[i2][i1-1]|&le;b1 and |u[i2][i1]-u[i2-1][i1]|&le;b2,
 * where b1 and b2 are again bounds on local average slopes.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.09.30
 */
public class DynamicWarping {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u in f[i1] ~ g[i1+u].
   * @param shiftMax upper bound on shift u in f[i1] ~ g[i1+u].
   */
  public DynamicWarping(int shiftMin, int shiftMax) {
    Check.argument(shiftMax-shiftMin>1,"shiftMax-shiftMin>1");
    _lmin = shiftMin;
    _lmax = shiftMax;
    _nl = 1+_lmax-_lmin;
    _bstretch1 = 1;
    _bstretch2 = 1;
  }

  /**
   * Sets bound on stretching (and squeezing) for all dimensions.
   * @param stretchMax the bound, a value less than or equal to one.
   */
  public void setStretchMax(double stretchMax) {
    Check.argument(stretchMax<=1.0,"stretchMax<=1.0");
    Check.argument(stretchMax>0.0,"stretchMax>0.0");
    setStretchMax(stretchMax,stretchMax);
  }

  /**
   * Sets bound on stretching (and squeezing) in 1st and 2nd dimensions.
   * @param stretchMax1 the bound for the 1st dimension.
   * @param stretchMax2 the bound for the 2nd dimension.
   */
  public void setStretchMax(double stretchMax1, double stretchMax2) {
    Check.argument(stretchMax1<=1.0,"stretchMax1<=1.0");
    Check.argument(stretchMax2<=1.0,"stretchMax2<=1.0");
    Check.argument(stretchMax1>0.0,"stretchMax1>0.0");
    Check.argument(stretchMax2>0.0,"stretchMax2>0.0");
    _bstretch1 = (int)(1.0/stretchMax1+0.5);
    _bstretch2 = (int)(1.0/stretchMax2+0.5);
  }

  /**
   * Returns squared errors for all samples and lags.
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax]. Squared errors are then computed as
   * (f[i1]-g[i1+il+shiftMin])^2.
   * @param f array[n1] for the sequence f[i1].
   * @param g array[n1] for the sequence g[i1].
   * @return array[n1][nl] of squared errors.
   */
  public float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    int nl = _nl;
    int n1m = n1-1;
    float[][] e = new float[n1][nl];
    // 0 <= il < nl, where il is index for lag
    // 0 <= i1 < n1, where i1 is index for sequence f
    // 0 <= j1 < n1, where j1 index for sequence g
    // j1 = i1+il+lmin, where il+lmin = lag
    // 0 <= i1+il+lmin < n1, so that j1 is in bounds
    // max(0,-lmin-i1) <= il < min(nl,n1-lmin-i1)
    // max(0,-lmin-il) <= i1 < min(n1,n1-lmin-il)
    // j1 = 0    => i1 = -lmin-il
    // j1 = n1-1 => i1 = n1-1-lmin-il
    for (int i1=0; i1<n1; ++i1) {
      int illo = min(nl-1,max(0,-_lmin-i1)); // see notes
      int ilhi = max(0,min(nl,n1-_lmin-i1)); // above
      for (int il=0,j1=i1+il+_lmin; il<illo; ++il,++j1)
        e[i1][il] = (j1>=0) ?
          error(f[i1],g[j1]) : 
          error(f[-_lmin-il],g[0]);
      for (int il=illo,j1=i1+il+_lmin; il<ilhi; ++il,++j1)
        e[i1][il] = error(f[i1],g[j1]);
      for (int il=ilhi,j1=i1+il+_lmin; il<nl; ++il,++j1)
        e[i1][il] = (j1<n1) ?
          error(f[i1],g[j1]) : 
          error(f[n1-1-_lmin-il],g[n1-1]);
    }
    return e;
  }

  /**
   * Returns squared errors for all samples and lags.
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax]. Squared errors are then computed as
   * (f[i2][i1]-g[i2][i1+il+shiftMin])^2.
   * @param f array[n2][n1] for the image f[i2][i1].
   * @param g array[n2][n1] for the sequence g[i2][i1].
   * @return array[n2][n1][nl] of squared errors.
   */
  public float[][][] computeErrors(float[][] f, float[][] g) {
    int n2 = f.length;
    float[][][] e = new float[n2][][];
    for (int i2=0; i2<n2; ++i2)
      e[i2] = computeErrors(f[i2],g[i2]);
    return e;
  }

  /**
   * Returns the sum of errors for specified shifts.
   * @param e array[n1][nl] of errors.
   * @param u array[n1] of shifts.
   * @return the sum of errors.
   */
  public float sumErrors(float[][] e, float[] u) {
    int n1 = e.length;
    int nl = e[0].length;
    float ul = 0.5f-_lmin;
    double sum = 0.0;
    for (int i1=0; i1<n1; ++i1) {
      int il = (int)(u[i1]+ul);
      il = max(0,min(nl-1,il));
      sum += e[i1][il];
    }
    return (float)sum;
  }

  /**
   * Returns the sum of errors for specified shifts.
   * @param e array[n2][n1][nl] of errors.
   * @param u array[n2][n1] of shifts.
   * @return the sum of errors.
   */
  public float sumErrors(float[][][] e, float[][] u) {
    int n2 = e.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2)
      sum += sumErrors(e[i2],u[i2]);
    return (float)sum;
  }

  public float[][] accumulate(float[][] e) {
    float[][] ea = new float[e.length][e[0].length];
    accumulate(e,ea);
    return ea;
  }
  public float[][] accumulateForward(float[][] e) {
    float[][] ea = new float[e.length][e[0].length];
    accumulateForward(e,ea);
    return ea;
  }
  public float[][] accumulateReverse(float[][] e) {
    float[][] ea = new float[e.length][e[0].length];
    accumulateReverse(e,ea);
    return ea;
  }
  public float[][][] accumulate1(float[][][] e) {
    float[][][] ea = new float[e.length][e[0].length][e[0][0].length];
    accumulate1(e,ea);
    return ea;
  }
  public float[][][] accumulateForward1(float[][][] e) {
    float[][][] ea = new float[e.length][e[0].length][e[0][0].length];
    accumulateForward1(e,ea);
    return ea;
  }
  public float[][][] accumulateReverse1(float[][][] e) {
    float[][][] ea = new float[e.length][e[0].length][e[0][0].length];
    accumulateReverse1(e,ea);
    return ea;
  }
  public float[][][] accumulate2(float[][][] e) {
    float[][][] ea = new float[e.length][e[0].length][e[0][0].length];
    accumulate2(e,ea);
    return ea;
  }
  public float[][][] accumulateForward2(float[][][] e) {
    float[][][] ea = new float[e.length][e[0].length][e[0][0].length];
    accumulateForward2(e,ea);
    return ea;
  }
  public float[][][] accumulateReverse2(float[][][] e) {
    float[][][] ea = new float[e.length][e[0].length][e[0][0].length];
    accumulateReverse2(e,ea);
    return ea;
  }

  public void accumulate(float[][] e, float[][] ea) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulateForward(e,ef);
    accumulateReverse(e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        ea[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }
  public void accumulateForward(float[][] e, float[][] ea) {
    accumulate( 1,_bstretch1,e,ea);
  }
  public void accumulateReverse(float[][] e, float[][] ea) {
    accumulate(-1,_bstretch1,e,ea);
  }
  public void accumulate1(float[][][] e, float[][][] ea) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulate(e[i2],ea[i2]);
  }
  public void accumulateForward1(float[][][] e, float[][][] ea) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateForward(e[i2],ea[i2]);
  }
  public void accumulateReverse1(float[][][] e, float[][][] ea) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateReverse(e[i2],ea[i2]);
  }
  public void accumulate2(float[][][] e, float[][][] ea) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] eai1 = new float[n2][];
    float[][] ef = new float[n2][nl];
    float[][] er = new float[n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
         ei1[i2] =  e[i2][i1];
        eai1[i2] = ea[i2][i1];
        for (int il=0; il<nl; ++il) {
          ef[i2][il] = 0.0f;
          er[i2][il] = 0.0f;
        }
      }
      accumulate( 1,_bstretch2,ei1,ef);
      accumulate(-1,_bstretch2,ei1,er);
      for (int i2=0; i2<n2; ++i2)
        for (int il=0; il<nl; ++il)
          eai1[i2][il] = ef[i2][il]+er[i2][il]-ei1[i2][il];
    }
  }
  public void accumulateForward2(float[][][] e, float[][][] ea) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] eai1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
         ei1[i2] =  e[i2][i1];
        eai1[i2] = ea[i2][i1];
      }
      accumulate( 1,_bstretch2,ei1,eai1);
    }
  }
  public void accumulateReverse2(float[][][] e, float[][][] ea) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] eai1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
         ei1[i2] =  e[i2][i1];
        eai1[i2] = ea[i2][i1];
      }
      accumulate(-1,_bstretch2,ei1,eai1);
    }
  }

  public float[] findShiftsReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    findShiftsReverse(d,e,u);
    return u;
  }
  public float[][] findShiftsReverse1(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    findShiftsReverse1(d,e,u);
    return u;
  }
  public float[][] findShiftsReverse2(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    findShiftsReverse2(d,e,u);
    return u;
  }

  public void findShiftsReverse(float[][] d, float[][] e, float[] u) {
    findShifts(-1,_bstretch1,_lmin,d,e,u);
  }
  public void findShiftsReverse1(float[][][] d, float[][][] e, float[][] u) {
    int n2 = d.length;
    for (int i2=0; i2<n2; ++i2)
      findShiftsReverse(d[i2],e[i2],u[i2]);
  }
  public void findShiftsReverse2(float[][][] d, float[][][] e, float[][] u) {
    int n1 = d[0].length;
    int n2 = d.length;
    float[][] di1 = new float[n2][];
    float[][] ei1 = new float[n2][];
    float[] ui1 = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        di1[i2] = d[i2][i1];
        ei1[i2] = e[i2][i1];
      }
      findShifts(-1,_bstretch2,_lmin,di1,ei1,ui1);
      for (int i2=0; i2<n2; ++i2)
        u[i2][i1] = ui1[i2];
    }
  }

  public static float[][] transposeLag(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] t = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        t[il][i1] = e[i1][il];
      }
    }
    return t;
  }
  public static float[][][] transposeLag(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[il][i2][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _lmin,_lmax,_nl;
  private int _bstretch1,_bstretch2;

  private static float error(float f, float g) {
    float e = f-g;
    return e*e;
  }

  private static void xaccumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int n1 = e.length;
    int nlm1 = nl-1;
    int n1m1 = n1-1;
    int i1b = (dir>0)?0:n1m1;
    int i1e = (dir>0)?n1:-1;
    int i1s = (dir>0)?1:-1;
    for (int i1=i1b; i1!=i1e; i1+=i1s) {
      int j1i = max(0,min(n1m1,i1-i1s));
      int j1b = max(0,min(n1m1,i1-i1s*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1<=-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[j1b][ilm1];
        float di = d[j1i][il  ];
        float dp = d[j1b][ilp1];
        for (int k1b=j1i; k1b!=j1b; k1b-=i1s) {
          dm += e[k1b][ilm1];
          dp += e[k1b][ilp1];
        }
        d[i1][il] = min3(dm,di,dp)+e[i1][il];
      }
    }
  }

  private static void xfindShifts(
    int dir, int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    int nl = d[0].length;
    int n1 = d.length;
    int nlm1 = nl-1;
    int n1m1 = n1-1;
    int i1b = (dir>0)?0:n1m1;
    int i1e = (dir>0)?n1m1:0;
    int i1s = (dir>0)?1:-1;
    int i1 = i1b;
    int il = 0;
    float dl = d[i1][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[i1][jl]<dl) {
        dl = d[i1][jl];
        il = jl;
      }
    }
    u[i1] = il+lmin;
    while (i1!=i1e) {
      int j1i = max(0,min(n1m1,i1+i1s));
      int j1b = max(0,min(n1m1,i1+i1s*b));
      int ilm1 = il-1; if (ilm1<=-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[j1b][ilm1];
      float di = d[j1i][il  ];
      float dp = d[j1b][ilp1];
      for (int k1b=j1i; k1b!=j1b; k1b+=i1s) { // i1-1
        dm += e[k1b][ilm1];
        dp += e[k1b][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      i1 += i1s;
      u[i1] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        for (int k1b=j1i; k1b!=j1b; k1b+=i1s) {
          i1 += i1s;
          u[i1] = il+lmin;
        }
      }
    }
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1<=-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param lmin minimum lag corresponding to lag index zero.
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private static void findShifts(
    int dir, int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?nim1:0;
    int is = (dir>0)?1:-1;
    int ii = ib;
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il+lmin;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1<=-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=is) { // ii-1
        dm += e[kb][ilm1];
        dp += e[kb][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      ii += is;
      u[ii] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        for (int kb=ji; kb!=jb; kb+=is) {
          ii += is;
          u[ii] = il+lmin;
        }
      }
    }
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  /////////////////////////////////////////////////////////////////////////////
  // Listings in paper
  void xxaccumulate(int dir, int b, 
                  float[][] e, float[][] d) 
  {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int i=ib; i!=ie; i+=is) {
      int ji = max(0,min(nim1,i-is));
      int jb = max(0,min(nim1,i-is*b));
      for (int l=0; l<nl; ++l) {
        int lm1 = l-1; if (lm1<=-1) lm1 = 0;
        int lp1 = l+1; if (lp1==nl) lp1 = nlm1;
        float dm = d[jb][lm1];
        float di = d[ji][l  ];
        float dp = d[jb][lp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][lm1];
          dp += e[kb][lp1];
        }
        d[i][l] = min3(dm,di,dp)+e[i][l];
      }
    }
  }
  void xxbacktrack(int dir, int b, int lmin, 
    float[][] d, float[][] e, float[] u) 
  {
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?nim1:0;
    int is = (dir>0)?1:-1;
    int ii = ib;
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il+lmin;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1<=-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=is) {
        dm += e[kb][ilm1];
        dp += e[kb][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      ii += is;
      u[ii] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        for (int kb=ji; kb!=jb; kb+=is) {
          ii += is;
          u[ii] = il+lmin;
        }
      }
    }
  }

}
