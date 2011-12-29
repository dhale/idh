package fault;

/** 
 * A k-d tree of samples in a k-dimensional space. Provides efficient 
 * searches for samples nearest to specified query points, and for
 * samples with coordinates within specified bounds.
 * <p>
 * Computational complexity to construct the tree is O(k*n*log(n)), 
 * where k is the number of coordinate values in each sample, and 
 * n is the number of samples. After the tree has been constructed, 
 * the average computational complexity for each search is O(log(n)). 
 * <p>
 * Storage required is O(n), not counting the storage of the k*n 
 * floats required for the sample coordinates. Those coordinates 
 * are specified in an array x[k][n] that is referenced (not copied) 
 * in the k-d tree. Searches return indices of samples stored in this
 * referenced array.
 * <p>
 * This implementation of the k-d tree is based on that described by 
 * Friedman, Bentley and Finkel (1977), An algorithm for finding best 
 * matches in logarithmic expected time: ACM Transactions on Mathematical
 * Software, v. 3, n. 3, p. 209-226.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.11.26
 */
public class KdTree {

  /**
   * A measure of distance between two points with k coordinates.
   * <p>
   * Distance is comprised of two functions. The first function f1 
   * computes a 1-dimensional coordinate distance for a specified 
   * coordinate. These non-negative coordinate distances f1 must 
   * satisfy properties of symmetry and monotonicity.
   * <p>
   * Given a sum of coordinate distances f1, the second function f2 
   * completes the computation of distance. The function f2 must be 
   * a monotonic (strictly increasing) function of its argument, the 
   * sum of coordinate distances f1.
   * <p>
   * In summary, distance = f2(f1(1,a1,b1)+f1(2,a2,b2)+...+f1(k,ak,bk)).
   * For example, simple Euclidean distance corresponds to functions
   * f1(i,ai,bi) = (ai-bi)*(ai-bi) and f2(f1sum) = sqrt(f1sum).
   * <p>
   * In the work of Bentley et al. (1977, see reference above), the 
   * first function f1 is called fi, and the second function f2 is 
   * called F.
   */
  public interface Distance {

    /**
     * Returns the i'th coordinate distance between two points a and b.
     * @param i the coordinate index in [1,k].
     * @param ai the i'th coordinate of point a.
     * @param bi the i'th coordinate of point b.
     * @return the i'th coordinate distance between points a and b.
     */
    public float f1(int i, float ai, float bi);

    /**
     * Returns the distance for a specified sum of coordinate distances.
     * This method completes the computation of distance between two
     * points a and b, each with k coordinates.
     * @param f1sum sum of coordinate distances f1 between a and b.
     * @return the distance between points a and b.
     */
    public float f2(float f1sum);
  }

  /**
   * Simple Euclidean distance between two points.
   */
  public static class EuclideanDistance implements Distance {
    public float f1(int i, float ai, float bi) {
      float ab = ai-bi;
      return ab*ab;
    }
    public float f2(float f1sum) {
      return (float)Math.sqrt(f1sum);
    }
  }

  /**
   * Constructs a k-d tree for specified sample coordinates.
   * Uses simple Euclidean distance.
   * @param x array[k][n] of n samples, each with k coordinates.
   */
  public KdTree(float[][] x) {
    this(x,new EuclideanDistance());
  }

  /**
   * Constructs a k-d tree for specified sample coordinates and distance.
   * @param x array[k][n] of n samples, each with k coordinates.
   * @param d the measure of distance.
   */
  public KdTree(float[][] x, Distance d) {
    _n = x[0].length;
    _k = x.length;
    _x = x;
    _d = d;
    _i = new int[_n];
    for (int i=0; i<_n; ++i)
      _i[i] = i;
    _root = new Node(_d,_x,0,_n-1,_i);
  }

  /**
   * Returns the index of the sample nearest to the specified point.
   * @param x array {x1,x2,...,xk} of point coordinates.
   * @return the index of the nearest sample.
   */
  public int findNearest(float[] x) {
    Search search = new Search(x);
    findNearest(_root,search);
    return search.i;
  }

  /**
   * Returns an array of indices of samples in the specified range.
   * @param xmin array of lower bounds of point coordinates.
   * @param xmax array of upper bounds of point coordinates.
   * @return array of sample indices.
   */
  public int[] findInRange(float[] xmin, float[] xmax) {
    BoxSearch search = new BoxSearch(xmin,xmax);
    findInRange(_root,search);
    return search.ilist.trim();
  }

  /**
   * Returns the distance from the i'th sample to the specified point.
   * @param i the index of the i'th sample in this tree.
   * @param x array {x1,x2,...,xk} of point coordinates.
   * @return the distance.
   */
  public float distance(int i, float[] x) {
    return _d.f2(f1sum(i,x));
  }

  /**
   * FOR TESTING ONLY.
   */
  private int findNearestSlow(float[] x) {
    float dmin = Float.MAX_VALUE;
    int imin = -1;
    for (int i=0; i<_n; ++i) {
      float d = distance(i,x);
      if (d<dmin) {
        dmin = d;
        imin = i;
      }
    }
    return imin;
  }

  /**
   * FOR TESTING ONLY.
   */
  private int[] findInRangeSlow(float[] xmin, float[] xmax) {
    IntList ilist = new IntList();
    for (int i=0; i<_n; ++i) {
      boolean outside = false;
      for (int j=0; j<_k && !outside; ++j) {
        float xji = _x[j][i];
        outside = xji<xmin[j] || xji>xmax[j];
      }
      if (!outside)
        ilist.add(i);
    }
    return ilist.trim();
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _k; // number of values per sample
  private int _n; // number of samples in tree
  private float[][] _x; // array[k][n] of sample values
  private Distance _d; // measure of distance
  private int[] _i; // array of sample indices
  private Node _root; // root node in tree

  private static final int NLEAF = 12; // min # of samples in leaf nodes

  private float f1sum(int i, float[] x) {
    float ds = 0.0f;
    for (int j=0; j<_k; ++j)
      ds += _d.f1(j,_x[j][i],x[j]);
    return ds;
  }

  /**
   * Returns true iff the search ball lies within the bounds box.
   * The bounds must contain the query point.
   */
  private boolean ballWithinBounds(Search search) {
    float[] x = search.x;
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;
    float sds = search.ds;
    for (int j=0; j<_k; ++j) {
      if (_d.f1(j,x[j],xmin[j])<=sds ||
          _d.f1(j,x[j],xmax[j])<=sds) 
        return false;
    }
    return true;
  }

  /**
   * Returns true iff the bounds box overlaps with the search ball.
   * The bounds box is assumed to not contain the search query point.
   * This method requires some careful explanation, partly because of
   * an error in the paper by Friedman et al (FBF, 1997, see above).
   * This method works by summing coordinate distances from the query 
   * point to the bounds box. If/when that sum exceeds the same sum
   * for the search ball, then we know that the bounds box is too
   * far away and cannot overlap the search ball, so we can quickly
   * return false without accumulating more coordinate distances.
   */
  private boolean boundsOverlapBall(Search search) {
    float[] x = search.x;
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;
    float sds = search.ds;
    float sum = 0.0f;
    for (int j=0; j<_k; ++j) {
      if (x[j]<xmin[j]) {
        float ds = _d.f1(j,x[j],xmin[j]); 
        sum += ds; 
        if (sum>sds) 
          return false; // not true <- FBF error
      } else if (x[j]>xmax[j]) {
        float ds = _d.f1(j,x[j],xmax[j]); 
        sum += ds; 
        if (sum>sds) 
          return false; // not true <- FBF error
      }
    }
    return true; // not false <- FBF error
  }

  /**
   * Returns the index of the sample nearest to a query point in a search.
   * After searching, the index in the search equals that returned by this
   * method.
   */
  private int findNearest(Search search) {
    findNearest(_root,search);
    return search.i;
  }

  private void findInRange(Node node, BoxSearch search) {
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;

    // If leaf node, append any samples in bounds box to search result.
    if (node.isLeaf()) {
      int p = node._p;
      int q = node._q;
      for (int i=p; i<=q; ++i) {
        boolean outside = false;
        for (int j=0; j<_k && !outside; ++j) {
          float xji = _x[j][_i[i]];
          outside = xji<xmin[j] || xji>xmax[j];
        }
        if (!outside)
          search.ilist.add(_i[i]);
      }
    }
    
    // Else if non-leaf node, search left and/or right children.
    else {
      int jv = node._j;
      float xv = node._x;
      if (xmin[jv]<=xv) 
        findInRange(node._left,search);
      if (xmax[jv]>=xv) 
        findInRange(node._right,search);
    }
  }

  /**
   * Searches recursively the specified node.
   * Returns true if the search is complete; false, otherwise.
   */
  private boolean findNearest(Node node, Search search) {

    // If leaf node, ...
    if (node.isLeaf()) {
      int p = node._p;
      int q = node._q;
      float[] x = search.x;

      // Examine its samples, updating nearest sample found.
      for (int i=p; i<=q; ++i) {
        float ds = f1sum(_i[i],x);
        if (ds<search.ds) {
          search.ds = ds;
          search.i = _i[i];
        }
      }

      // If ball is within bounds, then done.
      return ballWithinBounds(search);
    }

    // Dimension partitioned and value of non-leaf node.
    int jv = node._j;
    float xv = node._x;
    float[] x = search.x;
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;

    // Recursive call on the child that contains the query point.
    if (x[jv]<=xv) {
      float xmaxjv = xmax[jv]; 
      xmax[jv] = xv;
      if (findNearest(node._left,search))
        return true;
      xmax[jv] = xmaxjv;
    } else {
      float xminjv = xmin[jv]; 
      xmin[jv] = xv;
      if (findNearest(node._right,search))
        return true;
      xmin[jv] = xminjv;
    }

    // Recursive call on the other child, only if necessary.
    if (x[jv]<=xv) {
      float xminjv = xmin[jv];
      xmin[jv] = xv;
      if (boundsOverlapBall(search) && findNearest(node._right,search))
        return true;
      xmin[jv] = xminjv;
    } else {
      float xmaxjv = xmax[jv];
      xmax[jv] = xv;
      if (boundsOverlapBall(search) && findNearest(node._left,search))
        return true;
      xmax[jv] = xmaxjv;
    }
    return ballWithinBounds(search);
  }

  /**
   * A node in the tree.
   */
  private static class Node {

    Node(Distance d, float[][] x, int p, int q, int[] i) {

      // If number of samples in [p:q] is small, ...
      if (q-p<NLEAF) {

        // This is a leaf node with no children.
        _p = p;
        _q = q;
        _left = null;
        _right = null;
      }

      // Otherwise, if the number of samples in [p:q] is not small, ...
      else {

        // Determine the dimension j with maximum spread for samples in [p:q].
        // Here, spread is simply the coordinate distance f1 between the
        // minimum and maximum values for each of the k coordinates.
        int k = x.length;
        int n = x[0].length;
        int jsmax = 0;
        float xsmax = 0.0f;
        for (int j=0; j<k; ++j) {
          float xmin = x[j][i[p]];
          float xmax = xmin;
          for (int m=p+1; m<=q; ++m) {
            float xm = x[j][i[m]];
            if (xm<xmin) xmin = xm;
            if (xm>xmax) xmax = xm;
          }
          float xs = d.f1(j,xmin,xmax);
          if (xs>xsmax) {
            xsmax = xs;
            jsmax = j;
          }
        }
        _j = jsmax;

        // Split samples in [p:q] by median value in j'th dimension.
        int m = medianSplit(p,q,x[_j],i);
        _x = x[_j][i[m]];

        // Recursively make left and right children.
        _left = new Node(d,x,p,m,i);
        _right = new Node(d,x,m+1,q,i);
      }
    }

    boolean isLeaf() {
      return _left==null;
    }

    int _j; // dimension split by this node, if not a leaf node
    float _x; // the split value, if not a leaf node
    int _p,_q; // range [p:q] of indices spanned by this node
    Node _left,_right; // left and right children; null, if leaf

    /**
     * Partially sorts indices i[p:q] so that the median value in x[p:q] 
     * is x[i[m]], where the index m of the median is = (p+q)/2.
     */
    int medianSplit(int p, int q, float[] x, int[] i) {

      // Index of median is halfway between p and q.
      int m = (p+q)/2;

      // Partially sort the sample indices i such that
      // x[i[l]] <= x[i[m]], for p <= l <= m
      // x[i[l]] >= x[i[m]], for m <= l <= q
      while (p<q) {

        // Choose a pivot element between p and q.
        int pivot = (p+q)/2;
        float xpivot = x[i[pivot]];

        // Partition indices in the subarray i[p:q] so that there
        // exist integers r and s with the following properties:
        // p <= r < s <= q
        // x[i[l]] <= xpivot, for p <= l <= r
        // x[i[l]] == xpivot, for r <  l <  s
        // x[i[l]] >= xpivot, for s <= l <= q
        int s = p;
        int r = q;
        for (;;) {
          while (x[i[s]]<=xpivot && s<q) ++s;
          while (x[i[r]]>=xpivot && r>p) --r;
          if (s<r) {
            int is = i[s];
            i[s++] = i[r];
            i[r--] = is;
          } else {
            break;
          }
        }
        if (s<pivot) {
          int is = i[s];
          i[s++] = i[pivot];
          i[pivot] = is;
        } else if (pivot<r) {
          int ir = i[r];
          i[r--] = i[pivot];
          i[pivot] = ir;
        }

        // If median is in lower/upper subarray, partition that
        // subarray again; else x[i[m]] is the median value.
        if (m<=r) {
          q = r;
        } else if (m>=s) {
          p = s;
        } else {
          break;
        }
      }

      // Return index of median value.
      return m;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // May someday want to make this inner class public, to enable searches for 
  // m>1 nearest samples and/or incremental searching. For now we hide it.
  ///////////////////////////////////////////////////////////////////////////
  /**
   * A search within this tree. A search consists of the coordinates of
   * the most recent query point and search results.
   * <p> 
   * This class facilitates parallel (multithreaded) search within a k-d 
   * tree. A search cannot be shared by multiple threads. However, a k-d 
   * tree can be accessed by multiple threads, each with its own search.
   */
  private class Search {
    private float[] x = null; // coordinates of query point
    private int i = -1; // index of nearest point
    private float ds = Float.MAX_VALUE; // f1sum for nearest point
    private float[] xmin; // current lower bounds on coordinates
    private float[] xmax; // current upper bounds on coordinates
    private Search(float[] x) {
      this.x = x;
      this.i = -1;
      ds = Float.MAX_VALUE;
      if (xmin==null) {
        xmin = new float[x.length];
        xmax = new float[x.length];
      }
      for (int j=0; j<x.length; ++j) {
        xmin[j] = -Float.MAX_VALUE;
        xmax[j] =  Float.MAX_VALUE;
      }
    }
  }

  private class BoxSearch {
    private float[] xmin; // current lower bounds on coordinates
    private float[] xmax; // current upper bounds on coordinates
    private IntList ilist; // list of sample indices
    private BoxSearch(float[] xmin, float[] xmax) {
      this.xmin = xmin;
      this.xmax = xmax;
      this.ilist = new IntList();
    }
  }
  private static class IntList {
    public int n;
    public int[] a = new int[16];
    public void add(int i) {
      if (n==a.length) {
        int[] t = new int[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = i;
    }
    public int[] trim() {
      int[] t = new int[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
}
