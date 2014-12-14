/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warp;

import java.util.Random;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.util.MedianFinder;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING ONLY!
import java.awt.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.mosaic.*;

import dnp.CgSolver;
import dnp.Vec;
import dnp.VecArrayFloat2;

/**
 * Dynamic warping for alignment of well logs. 
 * <p>
 * This application of dynamic warping differs from others in that it must
 * account for missing log values, and the fact that values in different logs
 * are measured at different depths.
 * <p>
 * The alignment of two log sequences f[i] and g[j] is represented by a
 * sequence of index pairs (i,j). This representation is the same as that
 * proposed by Sakoe and Chiba (1978) in their description of what is now
 * known as dynamic time warping (DTW). However, unlike Sakoe and Chiba, we
 * need not assume that the first and last samples of the sequence f[i] are
 * aligned with the first and last samples of the sequence g[j]. Indeed,
 * that assumption is rarely valid for well log sequences.
 * <p>
 * As for DTW, the first step is to compute alignment errors for all (i,j),
 * subject to the constraint that |j-i|&le;lmax. The difference l = j-i is
 * called lag (or shift), and one can use specified constraints on geologic
 * dip and distances between wells to compute the maximum lag lmax.
 * <p>
 * As noted above, conventional DTW assumes that lag l is zero for the first
 * and last samples of the optimal path. To permit the optimal path to begin
 * and end with non-zero lag, alignment errors are computed in a rotated
 * coordinate system: e[k,l] = pow(abs(f[i]-g[j]),epow), where k = j+i and l =
 * j-i. Here, k = imin+jmin,...,imax+jmax, and l = -lmax,...,lmax. (The
 * zero-based array index for any lag l is simply l+lmax.) When accumulating
 * alignment errors, we begin at kmin = imin+jmin and end at kmax = imax+jmax.
 * <p>
 * Half of the samples in the array of alignment errors e[k,l] are unused,
 * because i = (k-l)/2 and j = (k+l)/2 are integers for only even values of
 * k+l. For display purposes only, errors e[k,l] for which k+l is an odd
 * number may be computed by linear interpolation of the other errors.
 * <p>
 * Alignment errors e[k,l] for two log sequences f[i] and g[j] are computed
 * only for a range of k = i+j for which at least one of the two logs has a
 * non-null value. Within this range, where either f[i] or g[j] is null, the
 * null values are replaced by a non-null value randomly selected from the
 * sequence. This replacement does not actually alter the log sequences, as it
 * occurs only during the computation of alignment errors. Outside of this
 * range, alignment errors are set to a null error. When accumulating
 * alignment errors, any null errors are ignored. Accumulation of alignment
 * errors begins and ends with the first and last indices k for which
 * alignment errors are not null.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.12.13
 */
public class WellLogWarping {

  /**
   * Sets the maximum shift (lag).
   * @param lmax the maximum lag.
   */
  public void setMaxShift(int lmax) {
    _lmax = lmax;
  }

  /**
   * Sets the exponent (the power) used to compute alignment errors.
   * @param epow the exponent.
   */
  public void setPowError(double epow) {
    _epow = (float)epow;
  }

  /**
   * Sets the alignment error that represents no computed error. 
   * @param enull the null error.
   */
  public void setNullError(float enull) {
    _enull = enull;
  }

  /**
   * Sets the log value that represents no measured value.
   * The default null value is -999.2500.
   * @param vnull the null value.
   */
  public void setNullValue(float vnull) {
    _vnull = vnull;
  }

  /**
   * Returns an array of alignment errors e[k,l] for two sequences.
   * @param f array of values f[i] in 1st log sequence.
   * @param g array of values g[j] in 2nd log sequence.
   * @return array of alignment errors e[k,l].
   */
  public float[][] computeErrors(float[] f, float[] g) {
    int ni = f.length;
    int nj = g.length;
    int lmax = min(max(ni,nj)-1,_lmax);
    int lmin = -lmax;
    int nl = 1+lmax-lmin;
    int nk = ni+nj-1;
    float[][] e = fillfloat(_enull,nl,nk);
    int[] igood = findGood(f);
    int[] jgood = findGood(g);
    int imin = min(igood);
    int imax = max(igood);
    int jmin = min(jgood);
    int jmax = max(jgood);
    int kmin = imin+jmin;
    int kmax = imax+jmax;
    //trace("computeErrors: kmin="+kmin+" kmax="+kmax);
    Random random = new Random(314159);
    for (int k=kmin; k<=kmax; ++k) {
      for (int l=lmin,ll=l-lmin; l<=lmax; ++l,++ll) {
        if ((k+l)%2==0) {
          int i = (k-l)/2;
          int j = (k+l)/2;
          float fi = value(random,igood,f,i);
          float gj = value(random,jgood,g,j);
          e[k][ll] = error(fi,gj);
        }
      }
    }
    return e;
  }

  /**
   * Returns accumulated errors d[k,l].
   * Any null errors in e[k,l] will be null in d[k,l].
   * @param e array of alignment errors e[k,l].
   * @return array of accumulated errors d[k,l].
   */
  public float[][] accumulateErrors(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    int lmax = (nl-1)/2;
    int lmin = -lmax;
    float[][] d = fillfloat(_enull,nl,nk);

    // Range of k for which to accumulate errors.
    int kmin = kminNotNull(e);
    int kmax = kmaxNotNull(e);
    //trace("accumulateErrors: kmin="+kmin+" kmax="+kmax);

    // Special case: k = kmin.
    int k=kmin,km1,km2;
    for (int l=lmin,ll=l-lmin; l<=lmax; ++l,++ll) {
      if ((k+l)%2==0)
        d[k][ll] = e[k][ll];
    }

    // Special case: k = kmin+1.
    k = kmin+1; km1 = k-1;
    for (int l=lmin,ll=l-lmin,lm=ll-1,lp=ll+1; l<=lmax; ++l,++lm,++ll,++lp) {
      if ((k+l)%2==0) {
        float da = lm>=0?d[km1][lm]:FLT_MAX;
        float dc = lp<nl?d[km1][lp]:FLT_MAX;
        d[k][ll] = min(da,dc)+e[k][ll];
      }
    }

    // General case: k = kmin+2 to kmax.
    for (k=kmin+2,km1=k-1,km2=k-2; k<=kmax; ++k,++km1,++km2) {
      for (int l=lmin,ll=l-lmin,lm=ll-1,lp=ll+1; l<=lmax; ++l,++lm,++ll,++lp) {
        if ((k+l)%2==0) {
          float da = lm>=0?d[km1][lm]:FLT_MAX;
          float db =       d[km2][ll];
          float dc = lp<nl?d[km1][lp]:FLT_MAX;
          float dm = dmin(da,db,dc);
          d[k][ll] = dm+e[k][ll];
        }
      }
    }
    return d;
  }

  /**
   * Returns the optimal warping path as pairs of sample indices (k,l).
   * The pairs are returned as an array of two arrays, one array for the
   * indices k and the other array for the corresponding indices l.
   * @param d array of accumulated errors.
   * @return array[2][] containing indices (k,l). The array[0] will
   *  contain the indices k, and the array[1] will contain the indices l. The
   *  lengths of these two arrays will equal the number of pairs on the
   *  warping path; this number is unknown before the optimal path has been
   *  found.
   */
  public int[][] findWarping(float[][] d) {
    int nk = d.length;
    int nl = d[0].length;
    int lmax = (nl-1)/2;
    int lmin = -lmax;

    // Range of k for which to find the optimal warping path.
    int kmin = kminNotNull(d);
    int kmax = kmaxNotNull(d);
    //trace("findWarping: kmin="+kmin+" kmax="+kmax);

    // Initially empty arrays for (k,l) pairs.
    int nw = 0;
    int[] kw = new int[nk];
    int[] lw = new int[nk];

    // Find lag l with minimum accumulated error at k = kmax.
    int kp = kmax;
    int lp = -1;
    float dp = FLT_MAX;
    for (int l=lmin,ll=0; l<=lmax; ++l,++ll) {
      if ((kp+l)%2==0) {
        if (d[kp][ll]<dp) {
          lp = l;
          dp = d[kp][ll];
        }
      }
    }

    // Add the corresponding pair (k,l) to the path.
    kw[0] = kp;
    lw[0] = lp;
    nw += 1;

    // While the path is not yet complete, backtrack.
    while (kp>kmin) {
      int ll = lp-lmin;
      float da = lp>lmin  ?d[kp-1][ll-1]:FLT_MAX;
      float db = kp>kmin+1?d[kp-2][ll  ]:FLT_MAX;
      float dc = lp<lmax  ?d[kp-1][ll+1]:FLT_MAX;
      float dm = dmin(da,db,dc);
      if (dm==db) {
        kp -= 2;
      } else if (dm==da) {
        kp -= 1;
        lp -= 1;
      } else {
        kp -= 1;
        lp += 1;
      }
      kw[nw] = kp;
      lw[nw] = lp;
      nw += 1;
    }

    // Remove any wasted space from the arrays of indices, while reordering
    // the indices (k,l) so that k are increasing, not decreasing.
    int[] kt = new int[nw];
    int[] lt = new int[nw];
    for (int mw=0; mw<nw; ++mw) {
      kt[mw] = kw[nw-1-mw];
      lt[mw] = lw[nw-1-mw];
    }
    return new int[][]{kt,lt};
  }

  /**
   * Returns an optimal warping path converted from (k,l) to (i,j).
   * Omits any indices (i,j) for which either f[i] is null or g[j] is null.
   * @param kl array[2][] containing indices (k,l). The array[0] contains
   *  the indices k, and the array[1] contains the indices l.
   * @param f array of values f[i] in 1st log sequence.
   * @param g array of values g[j] in 2nd log sequence.
   * @return array[2][] containing indices (i,j). The array[0] contains
   *  the indices i, and the array[1] contains the indices j.
   */
  public int[][] convertWarping(int[][] kl, float[] f, float[] g) {
    int ni = f.length;
    int nj = g.length;
    int[] ks = kl[0];
    int[] ls = kl[1];
    int nkl = ks.length;

    // Initially empty arrays for index pairs (i,j).
    int nij = 0;
    int[] is = new int[nkl];
    int[] js = new int[nkl];

    // Collect index pairs (i,j) for all non-null values.
    for (int ikl=0; ikl<nkl; ++ikl) {
      int k = ks[ikl];
      int l = ls[ikl];
      int i = (k-l)/2;
      int j = (k+l)/2;
      if (0<=i && i<ni && 0<=j && j<nj && f[i]!=_vnull && g[j]!=_vnull) {
        is[nij] = i;
        js[nij] = j;
        ++nij;
      }
    }

    // Return trimmed arrays of index pairs (i,j).
    is = copy(nij,is);
    js = copy(nij,js);
    return new int[][]{is,js};
  }

  /**
   * Apply warping to specified sequences f[i] and g[j].
   * The returned sequences f and g will both be indexed by k = i+j.
   * and will include null values for any missing values.
   * @param kl array[2][] containing indices (k,l). The array[0] contains
   *  the indices k, and the array[1] contains the indices l.
   * @param f array of values f[i] in 1st log sequence.
   * @param g array of values g[j] in 2nd log sequence.
   * @return array[2][] containing warped sequences f[k] and g[k].
   */
  public float[][] applyWarping(int[][] kl, float[] f, float[] g) {
    int[] ks = kl[0];
    int[] ls = kl[1];
    int nkl = ks.length;
    int ni = f.length;
    int nj = g.length;
    int nk = ni+nj-1;
    float[] fk = fillfloat(_vnull,nk);
    float[] gk = fillfloat(_vnull,nk);
    for (int ikl=0; ikl<nkl; ++ikl) {
      int k = ks[ikl];
      int l = ls[ikl];
      int i = (k-l)/2;
      int j = (k+l)/2;
      if (0<=i && i<ni && f[i]!=_vnull && 
          0<=j && j<nj && g[j]!=_vnull) {
        fk[k] = f[i];
        gk[k] = g[j];
      }
      if (ikl<nkl-1 && ks[ikl+1]==k+2) {
        fk[k+1] = fk[k];
        gk[k+1] = gk[k];
      }
    }
    return new float[][]{fk,gk};
  }

  /**
   * Interpolates alignment (or accumulated) errors for odd k+l.
   * Does not modify errors e[k,l] for which k+l is even.
   * <p>
   * Errors for odd k+l are never used. This interpolation is useful only for
   * displays of errors.
   * @param e array of errors e[k,l].
   */
  public void interpolateOddErrors(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    int lmax = (nl-1)/2;
    int lmin = -lmax;
    for (int k=0; k<nk; ++k) {
      for (int l=lmin; l<=lmax; ++l) {
        if ((k+l)%2!=0) {
          int km = k-1; if (km<   0) km += 2;
          int kp = k+1; if (kp>= nk) kp -= 2;
          int lm = l-1; if (lm<lmin) lm += 2;
          int lp = l+1; if (lp>lmax) lp -= 2;
          float ekm = e[km][l-lmin];
          float ekp = e[kp][l-lmin];
          float elm = e[k][lm-lmin];
          float elp = e[k][lp-lmin];
          if (ekm==_enull || ekp==_enull || elm==_enull || elp==_enull) {
            e[k][l-lmin] = _enull;
          } else {
            e[k][l-lmin] = 0.25f*(ekm+ekp+elm+elp);
          }
        }
      }
    }
  }

  /**
   * Returns shifts for each specified log. Logs are assumed to have been
   * resampled so that every log has the same depth sampling. The returned
   * shifts are in units of samples, but may have non-zero fractional parts.
   * @param fs array[nl][nz] of log values, nz values for each of nl logs.
   * @return array[nl][nt] of shifts, where nt = nz is the number of RGTs.
   */
  public float[][] findShifts(float[][] fs) {
    int nl = fs.length;
    int nz = fs[0].length;
    int nt = nz;

    // For each pair of logs, find all pairs (i,j) = [(iz,jz);(il,jl)] of
    // corresponding log samples. Each (i,j) pair corresponds to one equation
    // t[il][iz] = iz + s[il][iz]  ~=  t[jl][jz] = jz + s[jl][jz]
    // where t denotes relative geologic age (time). Corresponding log samples
    // should have equal t, but we have many such equations, so we seek a
    // least-squares solution of
    // s[il][iz] - s[jl][jz]  ~=  jz - iz
    // To prevent unnecessary shifting, squeezing, or stretching, we instead
    // solve for shifts r defined by z[il][it] = it + r[il][it], where z
    // denotes depth. Shifts r and s are related by 
    // s[il][z[il][it]] = r[il][t[il][iz]]. 
    // In least-squares solutions for shifts r we can impose the constraint
    // that the sum over il of r[il][it] = 0. which leads to a least-squares
    // solution of
    // r[il][t[il][iz]] - r[jl][t[jl][jz]]  ~=  jz - iz
    // Note that we need the mapping t[il][iz] from depth z to time t. 
    // So we first initialize r[il][iz] = 0, t[il][iz] = iz. Then,
    // while shifts r not yet converged {
    //   use t[il][iz] to solve for shifts r
    //   compute z[il][it] = it - r[il][it]
    //   interpolate t[il][it] from z[il][it]
    // }

    // Use dynamic warping to find pairs of corresponding log samples.
    int np = nl*(nl-1)/2; // number of pairs of logs
    Pairs[] ps = new Pairs[np];
    for (int il=0,ip=0; il<nl; ++il) {
      for (int jl=il+1; jl<nl; ++jl,++ip) {
        float[] fi = fs[il];
        float[] gj = fs[jl];
        float[][] e = computeErrors(fi,gj);
        float[][] d = accumulateErrors(e);
        int[][] kl = findWarping(d);
        int[][] ij = convertWarping(kl,fi,gj);
        ps[ip] = new Pairs(il,jl,ij[0],ij[1]);
      }
    }

    // Initial mapping t[il][iz] from depth z to time t for shifts r = 0.
    float[][] r = new float[nl][nt];
    int[][] t = new int[nl][nz];
    computeTzFromShifts(r,t);

    // Outer iterations over CG solutions.
    int maxouter = 10;
    boolean converged = false;
    for (int nouter=0; nouter<maxouter && !converged; ++nouter) {
      //zero(r); // zero shifts before each set of inner CG iterations?

      // Inner preconditioned conjugate-gradient iterations.
      int ninner = (nouter<maxouter-1)?10:50;
      CgSolver cs = new CgSolver(0.001f,ninner);
      float[][] b = makeRhs(ps,t);
      A a = new A(ps,t,10.0);
      M m = new M(100.0);
      VecArrayFloat2 vb = new VecArrayFloat2(b);
      VecArrayFloat2 vr = new VecArrayFloat2(r);
      cs.solve(a,m,vb,vr);
      //plotPoints(r);

      // Ensure monotonically increasing z[il][it] = it-r[il][it].
      cleanShifts(r);

      // Update t[il][iz] by inverse interpolation of z[il][it].
      computeTzFromShifts(r,t);
    }

    return r;
  }
  /**
   * Applies specified shifts to specified logs.
   * @param f array[nl][nz] of log values.
   * @param r array[nl][nt] of shifts.
   * @return array[nl][nt] of shifted log values.
   */
  public float[][] applyShifts(float[][] f, float[][] r) {
    int nl = f.length;
    int nz = f[0].length;
    int nt = r[0].length;
    float[][] g = fillfloat(_vnull,nt,nl);
    for (int il=0; il<nl; ++il) { // for all logs, ...
      for (int it=0; it<nt; ++it) { // for all times, ...
        float zi = it-r[il][it]; // depth at which to interpolate
        int iz = (int)(zi+0.5f); // index of nearest log sample
        if (0<=iz && iz<nz) {
          g[il][it] = f[il][iz]; // nearest-neighbor interpolation
        }
      }
    }
    return g;
  }

  /**
   * Sorts well indices to approximately minimize distances between wells.
   * Specifically, this method approximately minimizes the total distance
   * traveled from well to well in a sequential iteration over all well
   * locations, in which each well is visited only once. Because exactly 
   * minimizing this total distance would require a costly solution to the
   * traveling-salesman problem, this method instead uses a simple greedy
   * solution.
   * <p>
   * This method is useful primarily in 2D displays of logs. Logs displayed as
   * adjacent pixels or curves are likely to appear more correlated than they
   * would be for an arbitrary ordering of wells.
   * @param x array of x coordinates of well locations.
   * @param y array of y coordinates of well locations.
   * @return the sorted array of well indices.
   */
  public static int[] sortWells(double[] x, double[] y) {
    int nw = x.length;
    Random r = new Random(314159);
    double dsmin = DBL_MAX; // the minimized distance
    int[] ksmin = null;
    for (int mw=0; mw<nw; ++mw) { // for all possible first-well indices
      boolean[] bs = new boolean[nw]; // flags for visited wells
      int[] ks = new int[nw]; // indices of visited wells
      int kw = mw; // index of the first well
      double xk = x[kw]; // x-coordinate
      double yk = y[kw]; // y-coordinate
      double ds = 0.0f; // distance sum
      int iw = 0;
      ks[iw] = kw; // first index in the list
      bs[kw] = true; // have now visited well with index kw
      for (iw=1; iw<nw; ++iw) { // for all other wells, ...
        int jmin = -1;
        double dmin = DBL_MAX;
        for (int jw=0; jw<nw; ++jw) { // find nearest not yet visited
          if (!bs[jw]) { // if not yet visited, ...
            double xj = x[jw];
            double yj = y[jw];
            double dx = xk-xj;
            double dy = yk-yj;
            double dj = dx*dx+dy*dy; // distance-squared
            if (dj<dmin) { // if nearest so far, ...
              dmin = dj;
              jmin = jw;
            }
          }
        }
        kw = jmin; // visit the nearest well
        xk = x[kw];
        yk = y[kw];
        ds += dmin;
        ks[iw] = kw;
        bs[kw] = true; // mark this well as visited
      }
      //trace("sortWells: ds="+ds);
      if (ds<dsmin) { // if this path has less distance, remember it
        dsmin = ds;
        ksmin = ks;
      }
    }
    //trace("sortWells: dsmin="+dsmin);
    return ksmin;
  }

  /**
   * Gets a uniform sampling of depth for specified depths and logs.
   * Considers only depths for which log values are non-null. The returned
   * sampling will include the shallowest and the deepest depths logged.
   * @param z array of arrays of depths; one array for each log.
   * @param f array of arrays of log values; one array for each log.
   * @return the uniform sampling.
   */
  public Sampling getDepthSampling(float[][] z, float[][] f) {
    int nl = z.length;

    // Total number of depths specified.
    int nlz = 0;
    for (int il=0; il<nl; ++il)
      nlz += z[il].length;

    // Array for depth increments, and counter for number of increments.
    float[] dz = new float[nlz];
    int ndz = 0;

    // Find min and max depths, while storing depth increments.
    // Consider only log samples with non-null values.
    float zmin =  FLT_MAX;
    float zmax = -FLT_MAX;
    for (int il=0; il<nl; ++il) {
      int nz = z[il].length;
      float zi = z[il][0];
      float fi = f[il][0];
      for (int iz=0; iz<nz; ++iz) {
        float zim1 = zi;
        float fim1 = fi;
        zi = z[il][iz];
        fi = f[il][iz];
        if (fi!=_vnull && zi<zmin)
          zmin = zi;
        if (fi!=_vnull && zi>zmax)
          zmax = zi;
        if (iz>0 && fi!=_vnull && fim1!=_vnull)
          dz[ndz++] = zi-zim1;
      }
    }

    // Depth interval is median of all depth increments.
    dz = copy(ndz,dz);
    MedianFinder mf = new MedianFinder(ndz);
    float zdel = mf.findMedian(dz);

    // Uniform sampling.
    int nz = 1+(int)ceil((zmax-zmin)/zdel);
    return new Sampling(nz,zdel,zmin);
  }

  /**
   * Resamples a log with specified depths and values.
   * The desired depth sampling for the output array of values can be
   * different from (e.g., coarser than) that implied by the input array of
   * depths.
   * @param s the desired uniform sampling of depths.
   * @param z array of depths for which values are provided.
   * @param f array of values; some values may be null.
   * @return array of values for uniformly sampled depths.
   */
  public float[] resampleLog(Sampling s, float[] z, float[] f) {
    int nz = s.getCount();
    float zmin = (float)s.getFirst();
    float zmax = (float)s.getLast();
    int n = z.length;
    float[] g = new float[nz];
    float[] c = new float[nz];
    for (int i=0; i<n; ++i) {
      float zi = z[i];
      float fi = f[i];
      if (zmin<=zi && zi<=zmax && fi!=_vnull) {
        int iz = s.indexOfNearest(zi);
        g[iz] += fi;
        c[iz] += 1.0f;
      }
    }
    for (int iz=0; iz<nz; ++iz) {
      if (c[iz]>0.0f) {
        g[iz] /= c[iz];
      } else {
        g[iz] = _vnull;
      }
    }
    return g;
  }

  /**
   * Resamples multiple logs with specified depths and values.
   * This method simply calls the method 
   * {@link #resampleLog(Sampling,float[],float[])}
   * for all specified arrays of depths and values.
   * @param s the desired uniform sampling of depths.
   * @param z array of depths for which values are provided.
   * @param f array of values; some values may be null.
   * @return array of values for uniformly sampled depths.
   */
  public float[][] resampleLogs(Sampling s, float[][] z, float[][] f) {
    int n = z.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      g[i] = resampleLog(s,z[i],f[i]);
    return g;
  }

  /**
   * Replaces any null values found in a log with a specified value.
   * Typically used only for displays.
   * @param f array of log values.
   * @param freplace value used to replace any null values.
   * @return array of logs values with null values replaced.
   */
  public float[] replaceNulls(float[] f, float freplace) {
    int n = f.length;
    float[] g = new float[n];
    for (int i=0; i<n; ++i) {
      if (f[i]!=_vnull) {
        g[i] = f[i];
      } else {
        g[i] = freplace;
      }
    }
    return g;
  }
  /**
   * Replaces any null values found in multiple logs with a specified value.
   * Typically used only for displays.
   * @param f array of log values.
   * @param freplace value used to replace any null values.
   * @return array of logs values with null values replaced.
   */
  public float[][] replaceNulls(float[][] f, float freplace) {
    int n = f.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      g[i] = replaceNulls(f[i],freplace);
    return g;
  }

  public float[] toFloat(int[] i) {
    int n = i.length;
    float[] f = new float[n];
    for (int j=0; j<n; ++j)
      f[j] = (float)i[j];
    return f;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
 
  private int _lmax = Integer.MAX_VALUE;
  private float _epow = 1.0f;
  private float _enull = -FLT_MIN;
  private float _vnull = -999.2500f;

  /**
   * Arrays of pairs of depth sample indices (iz,jz) with weights w. The
   * arrays of depth sample indices were computed by warping a pair of logs
   * with indices (il,jl). All weights are initialized to one.
   */
  private static class Pairs {
    Pairs(int il, int jl, int[] izs, int[] jzs) {
      this.il = il;
      this.jl = jl;
      this.nzp = izs.length;
      this.izs = izs;
      this.jzs = jzs;
      this.ws = fillfloat(1.0f,nzp);
    }
    int il,jl; // a pair of log indices
    int nzp; // number of depth index pairs
    int[] izs,jzs; // arrays of pairs of depth indices
    float[] ws; // array of weights
  }

  /**
   * Conjugate-gradient operator A and preconditioner M.
   * The preconditioner smooths shifts r(t,l) along time t, and 
   * subracts the mean over logs (indexed by l) of r(t,l).
   */
  private static class A implements CgSolver.A {
    A(Pairs[] ps, int[][] t, double eps) {
      _ps = ps;
      _t = t;
      _eps = (float)eps;
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      applyLhs(_ps,_t,_eps,x,y);
    }
    private Pairs[] _ps;
    private int[][] _t;
    private float _eps;
  }
  private static class M implements CgSolver.A {
    public M(double sigma) {
      _ref = new RecursiveExponentialFilter(sigma);
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      int nl = x.length;
      int nt = x[0].length;
      copy(x,y);
      subtractMeanOverLogs(y);
      smoothOverTime(y);
      subtractMeanOverLogs(y);
    }
    private RecursiveExponentialFilter _ref;
    private void smoothOverTime(float[][] r) {
      _ref.apply1(r,r);
    }
    private void subtractMeanOverLogs(float[][] r) {
      int nl = r.length;
      int nt = r[0].length;
      float[] rsum = new float[nt];
      for (int il=0; il<nl; ++il)
        add(r[il],rsum,rsum); // sum over logs
      mul(1.0f/nl,rsum,rsum); // mean over logs
      for (int il=0; il<nl; ++il)
        sub(r[il],rsum,r[il]); // subtract mean
    }
  }

  private static float[][] makeRhs(Pairs[] ps, int[][] t) {
    int np = ps.length; // number of log pairs
    int nl = t.length; // number of logs
    int nz = t[0].length; // number of depths (and times)
    int nt = nz;
    float[][] y = new float[nl][nt]; // output rhs
    for (int ip=0; ip<np; ++ip) { // for all log pairs, ...
      Pairs p = ps[ip];
      int il = p.il;
      int jl = p.jl;
      int nzp = p.nzp;
      int[] izs = p.izs;
      int[] jzs = p.jzs;
      float[] ws = p.ws;
      float[] yi = y[il];
      float[] yj = y[jl];
      for (int kzp=0; kzp<nzp; ++kzp) { // for all index pairs (i,j), ...
        int iz = izs[kzp];
        int jz = jzs[kzp];
        int it = t[il][iz];
        int jt = t[jl][jz];
        if (0<=it && it<nt && 0<=jt && jt<nt) {
          float wk = ws[kzp];
          float scl = wk*wk;
          float dif = jz-iz;
          dif *= scl;
          yi[it] += dif;
          yj[jt] -= dif;
        }
      }
    }
    return y;
  }

  private static void applyLhs(
      Pairs[] ps, int[][] t, float eps, float[][] x, float[][] y) {
    int np = ps.length; // number of log pairs
    int nl = x.length; // number of logs
    int nt = x[0].length; // number of times (and depths)
    zero(y); // zero y before accumulating below
    for (int ip=0; ip<np; ++ip) { // for all log pairs, ...
      Pairs p = ps[ip];
      int il = p.il;
      int jl = p.jl;
      float[] xi = x[il];
      float[] xj = x[jl];
      float[] yi = y[il];
      float[] yj = y[jl];
      int nzp = p.nzp;
      int[] izs = p.izs;
      int[] jzs = p.jzs;
      float[] ws = p.ws;
      for (int kzp=0; kzp<nzp; ++kzp) {
        int iz = izs[kzp];
        int jz = jzs[kzp];
        int it = t[il][iz];
        int jt = t[jl][jz];
        if (0<=it && it<nt && 0<=jt && jt<nt) {
          float wk = ws[kzp];
          float scl = wk*wk;
          float dif = 0.0f;
          dif += xi[it];
          dif -= xj[jt];
          dif *= scl;
          yi[it] += dif;
          yj[jt] -= dif;
        }
      }
      float scl = eps*eps;
      for (int kl=0; kl<nl; ++kl) {
        for (int kt=1; kt<nt; ++kt) {
          float dif = 0.0f;
          dif += x[kl][kt  ];
          dif -= x[kl][kt-1];
          dif *= scl;
          y[kl][kt  ] += dif;
          y[kl][kt-1] -= dif;
        }
      }
    }
  }

  private static float dmin(float a, float b, float c) {
    float d = b;
    if (a<d) d = a;
    if (c<d) d = c;
    return d;
  }

  private float error(float f, float g) {
    float d = f-g;
    if (_epow==2.0f) {
      return d*d;
    } else if (_epow==1.0f) {
      return abs(d);
    } else {
      if (d<0.0f) d = -d;
      return pow(d,_epow);
    }
  }

  private int[] findGood(float[] f) {
    int n = f.length;
    int[] igood = new int[n];
    int ngood = 0;
    for (int i=0; i<n; ++i) {
      if (f[i]!=_vnull) {
        igood[ngood] = i;
        ++ngood;
      }
    }
    return copy(ngood,igood);
  }

  private float value(Random random, int[] igood, float[] f, int i) {
    int n = f.length;
    if (0<=i && i<n && f[i]!=_vnull) {
      return f[i];
    } else {
      int j = random.nextInt(igood.length);
      return f[igood[j]];
    }
  }

  private int kminNotNull(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    for (int k=0; k<nk; ++k) {
      if (e[k][0]!=_enull || e[k][1]!=_enull)
        return k;
    }
    return nk;
  }

  private int kmaxNotNull(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    for (int k=nk-1; k>=0; --k) {
      if (e[k][0]!=_enull || e[k][1]!=_enull)
        return k;
    }
    return -1;
  }

  /**
   * Linear interpolation and constant extrapolation of known shifts s. Known
   * shifts are marked by non-zero counts c. That is, c[k]==0.0f if and only
   * if s[k] is unknown. This method assumes that s contains at least one
   * known shift.
   */
  private static float[] interpolateShifts(float[] c, float[] s) {
    int n = s.length;
    float[] t = new float[n];
    int klo = -1; // lower index of known shift
    int kup = -1; // upper index of known shift
    for (int k=0; k<n; ++k) {
      if (c[k]!=0.0f) { // if shift is known, ...
        t[k] = s[k]; // copy it
        klo = k; // update index of last known shift
      } else { // else, if shift is unknown, ...
        if (kup<k) { // if necessary, find next known shift
          for (kup=k+1; kup<n && c[kup]==0.0f; ++kup)
            ;
        }
        if (klo<0) { // if no known shift with index klo < k
          t[k] = s[kup]; // copy shift at upper index kup
        } else if (kup>=n) { // else if no known shift with index kup > k
          t[k] = s[klo]; // copy shift with lower index klo
        } else { // else, linearly interpolate two known shifts
          float whi = (float)(k-klo)/(float)(kup-klo);
          float wlo = 1.0f-whi;
          t[k] = wlo*s[klo]+whi*s[kup];
        }
      }
    }
    return t;
  }
  private static float[][] interpolateShifts(float[][] c, float[][] s) {
    int n = s.length;
    float[][] t = new float[n][];
    for (int l=0; l<n; ++l)
      t[l] = interpolateShifts(c[l],s[l]);
    return t;
  }

  private static float[][] interpolateShiftsSlow(float[][] c, float[][] s) {
    int nl = s.length;
    int nk = s[0].length;
    float[][] si = new float[nl][nk];
    for (int il=0; il<nl; ++il) {
      for (int ik=0; ik<nk; ++ik) {
        if (c[il][ik]!=0.0f) {
          si[il][ik] = s[il][ik];
        } else {
          int iklo = ik-1;
          while (iklo>=0 && c[il][iklo]==0.0f)
            --iklo;
          int ikhi = ik+1;
          while (ikhi<nk && c[il][ikhi]==0.0f)
            ++ikhi;
          if (iklo<0) {
            si[il][ik] = s[il][ikhi];
          } else if (ikhi>=nk) {
            si[il][ik] = s[il][iklo];
          } else {
            float whi = (float)(ik-iklo)/(ikhi-iklo);
            float wlo = 1.0f-whi;
            si[il][ik] = wlo*s[il][iklo]+whi*s[il][ikhi];
          }
        }
      }
    }
    return si;
  }

  private static void cleanShifts(float[] r) {
    int nt = r.length;
    for (int it=1; it<nt; ++it) {
      if (r[it]>r[it-1]+0.99f)
        r[it] = r[it-1]+0.99f;
    }
    r[nt-1] = r[nt-2]; // ensures constant extrapolation of shifts
  }
  private static void cleanShifts(float[][] r) {
    for (float[] ri:r)
      cleanShifts(ri);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  /**
   * Linear interpolation of a function y(x).
   * @param xi input array of monotonically increasing x.
   * @param yi input array of y(x), for all x in array xi.
   * @param xo input array of x for which to interpolate y(x).
   * @param yo output array of y(x), for all x in array xo.
   */
  private static void interpolate(
      float[] xi, float[] yi, float[] xo, float[] yo) {
    CubicInterpolator ci = 
      new CubicInterpolator(CubicInterpolator.Method.LINEAR,xi,yi);
    ci.interpolate(xo,yo);
  }

  /**
   * Computes t(z,l) from shifts r(t,l).
   */
  private static void computeTzFromShifts(float[][] r, int[][] t) {
    int nl = r.length;
    int nt = r[0].length;
    int nz = t[0].length;
    float[] uz = rampfloat(0.0f,1.0f,nz); // z = 0,1,2,...,nz-1
    float[] ut = rampfloat(0.0f,1.0f,nt); // t = 0,1,2,...,nt-1
    float[] zt = new float[nt]; // for z[il][it]
    float[] tz = new float[nz]; // for t[il][iz]
    for (int il=0; il<nl; ++il) {
      for (int it=0; it<nt; ++it) {
        zt[it] = it-r[il][it];
      }
      interpolate(zt,ut,uz,tz);
      for (int iz=0; iz<nz; ++iz) {
        float tzi = tz[iz];
        t[il][iz] = (int)(tzi+0.5f);
      }
    }
  }

  private static float[][] shiftErrors(Pairs[] ps, float[][] s) {
    int np = ps.length;
    int nl = s.length;
    int nz = s[0].length;
    float[][] e = new float[nl][nz];
    float[][] c = new float[nl][nz];
    for (Pairs p:ps) {
      int il = p.il;
      int jl = p.jl;
      int nzp = p.nzp;
      int[] izs = p.izs;
      int[] jzs = p.jzs;
      for (int kzp=0; kzp<nzp; ++kzp) {
        int iz = izs[kzp];
        int jz = jzs[kzp];
        float ez = (iz+s[il][iz])-(jz+s[jl][jz]);
        float es = ez*ez;
        e[il][iz] += es; c[il][iz] += 1.0f;
        e[jl][jz] += es; c[jl][jz] += 1.0f;
      }
    }
    for (int il=0; il<nl; ++il) {
      for (int iz=0; iz<nz; ++iz) {
        e[il][iz] = c[il][iz]>0.0f ? sqrt(e[il][iz]/c[il][iz]) : 0.0f;
      }
    }
    return e;
  }

  private static void plotPoints(float[][] s) {
    Color[] colors = {
      Color.BLACK,
      Color.RED,Color.GREEN,Color.BLUE,
      Color.CYAN,Color.MAGENTA,Color.YELLOW,
    };
    int ncolor = colors.length;
    SimplePlot sp = new SimplePlot();
    int ns = s.length;
    for (int is=0; is<ns; ++is) {
      PointsView pv = sp.addPoints(s[is]);
      pv.setLineColor(colors[is%ncolor]);
    }
  }

  private static void plotPixels(float[][] s) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    PixelsView pv = sp.addPixels(s);
    pv.setColorModel(ColorMap.JET);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    sp.addColorBar();
  }
    

}
