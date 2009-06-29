package dnp;

import edu.mines.jtk.sgl.Point3;
import edu.mines.jtk.sgl.Vector3;

/*
*		  		
*  Triangle-Triangle Overlap Test Routines				
*  July, 2002                                                          
*  Updated December 2003                                                
*                                                                       
*  This file contains C implementation of algorithms for                
*  performing two and three-dimensional triangle-triangle intersection test 
*  The algorithms and underlying theory are described in                    
*                                                                           
* "Fast and Robust Triangle-Triangle Overlap Test 
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*                                                 
*  Journal of Graphics Tools, 8(1), 2003                                    
*                                                                           
*  Several geometric predicates are defined.  Their parameters are all      
*  points.  Each point is an array of two or three double precision         
*  floating point numbers. The geometric predicates implemented in          
*  this file are:                                                            
*                                                                           
*    int tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)                         
*    int tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)                         
*                                                                           
*    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
*                                     coplanar,source,target)               
*                                                                           
*       is a version that computes the segment of intersection when            
*       the triangles overlap (and are not coplanar)                        
*                                                                           
*    each function returns 1 if the triangles (including their              
*    boundary) intersect, otherwise 0                                       
*                                                                           
*                                                                           
*  Other information are available from the Web page                        
*  http://www.acm.org/jgt/papers/GuigueDevillers03/                         
*                                                                           
*/
public class Intersect {

  /**
   * Computes the line-segment intersection of two triangles in 3D. 
   * If the triangles are co-planar, they are assumed to not intersect.
   * @param p1 vertex p of 1st triangle.
   * @param q1 vertex q of 1st triangle.
   * @param r1 vertex r of 1st triangle.
   * @param p2 vertex p of 2nd triangle.
   * @param q2 vertex q of 2nd triangle.
   * @param r2 vertex r of 2nd triangle.
   * @param a endpoint of line segment.
   * @param b endpoint of line segment.
   * @return true, if the triangles intersect; false, otherwise.
   */
  public static boolean triTri(
    Point3 p1, Point3 q1, Point3 r1, 
    Point3 p2, Point3 q2, Point3 r2,
    Point3 a, Point3 b)
  {
    double dp1,dq1,dr1,dp2,dq2,dr2;

    // Temporary vectors.
    Vector3 v1 = new Vector3();
    Vector3 v2 = new Vector3();
    Vector3 v = new Vector3();
    Vector3 n1 = new Vector3();
    Vector3 n2 = new Vector3();
    Vector3 n = new Vector3();

    // Distance signs of p1, q1 and r1 to the plane of triangle(p2,q2,r2).
    sub(p2,r2,v1);
    sub(q2,r2,v2);
    cross(v1,v2,n2);
    sub(p1,r2,v1);
    dp1 = dot(v1,n2);
    sub(q1,r2,v1);
    dq1 = dot(v1,n2);
    sub(r1,r2,v1);
    dr1 = dot(v1,n2);
    if (((dp1*dq1)>0.0) && ((dp1*dr1)>0.0)) 
      return false;

    // Distance signs of p2, q2 and r2 to the plane of triangle(p1,q1,r1).
    sub(q1,p1,v1);
    sub(r1,p1,v2);
    cross(v1,v2,n1);
    sub(p2,r1,v1);
    dp2 = dot(v1,n1);
    sub(q2,r1,v1);
    dq2 = dot(v1,n1);
    sub(r2,r1,v1);
    dr2 = dot(v1,n1);
    if (((dp2*dq2)>0.0) && ((dp2*dr2)>0.0)) 
      return false;

    // Permutation in a canonical form of T1's vertices
    if (dp1>0.0) {
      if (dq1>0.0) 
        return triTri(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2,v,v1,v2,n,n1,n2,a,b);
      else if (dr1>0.0)
        return triTri(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2,v,v1,v2,n,n1,n2,a,b);
      else 
        return triTri(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2,v,v1,v2,n,n1,n2,a,b);
    } else if (dp1<0.0) {
      if (dq1<0.0) 
        return triTri(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2,v,v1,v2,n,n1,n2,a,b);
      else if (dr1<0.0)
        return triTri(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2,v,v1,v2,n,n1,n2,a,b);
      else 
        return triTri(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2,v,v1,v2,n,n1,n2,a,b);
    } else {
      if (dq1<0.0) {
        if (dr1>=0.0)
          return triTri(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2,v,v1,v2,n,n1,n2,a,b);
        else 
          return triTri(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2,v,v1,v2,n,n1,n2,a,b);
      }
      else if (dq1>0.0) {
        if (dr1>0.0)
          return triTri(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2,v,v1,v2,n,n1,n2,a,b);
        else 
          return triTri(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2,v,v1,v2,n,n1,n2,a,b);
      }
      else  {
        if (dr1>0.0)
          return triTri(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2,v,v1,v2,n,n1,n2,a,b);
        else if (dr1<0.0)
          return triTri(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2,v,v1,v2,n,n1,n2,a,b);
        else {
          // triangles are co-planar
          return false;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static void cross(Vector3 v1, Vector3 v2, Vector3 v3) {
    v3.x = v1.y*v2.z-v1.z*v2.y;
    v3.y = v1.z*v2.x-v1.x*v2.z;
    v3.z = v1.x*v2.y-v1.y*v2.x;
  }

  private static double dot(Vector3 v1, Vector3 v2) {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  }

  private static void scale(double s, Vector3 v1, Vector3 v2) {
    v2.x = s*v1.x;
    v2.y = s*v1.y;
    v2.z = s*v1.z;
  }

  private static void sub(Point3 p1, Point3 p2, Vector3 v3) {
    v3.x = p1.x-p2.x;
    v3.y = p1.y-p2.y;
    v3.z = p1.z-p2.z;
  }

  private static void sub(Point3 p1, Vector3 v1, Point3 p2) {
    p2.x = p1.x-v1.x;
    p2.y = p1.y-v1.y;
    p2.z = p1.z-v1.z;
  }

  /**
   * Computes the segment [a,b] of intersection of two triangles.
   */
  private static boolean constructIntersection(
    Point3 p1, Point3 q1, Point3 r1, 
    Point3 p2, Point3 q2, Point3 r2,
    Vector3 v, Vector3 v1, Vector3 v2,
    Vector3 n, Vector3 n1, Vector3 n2,
    Point3 a, Point3 b)
  {
    double alpha;
    sub(q1,p1,v1);
    sub(r2,p1,v2);
    cross(v1,v2,n);
    sub(p2,p1,v);
    if (dot(v,n)>0.0) {
      sub(r1,p1,v1);
      cross(v1,v2,n);
      if (dot(v,n)<=0.0) {
        sub(q2,p1,v2);
        cross(v1,v2,n);
        if (dot(v,n)>0.0) {
          sub(p1,p2,v1);
          sub(p1,r1,v2);
          alpha = dot(v1,n2)/dot(v2,n2);
          scale(alpha,v2,v1);
          sub(p1,v1,a);
          sub(p2,p1,v1);
          sub(p2,r2,v2);
          alpha = dot(v1,n1)/dot(v2,n1);
          scale(alpha,v2,v1);
          sub(p2,v1,b);
          return true;
        } else {
          sub(p2,p1,v1);
          sub(p2,q2,v2);
          alpha = dot(v1,n1)/dot(v2,n1);
          scale(alpha,v2,v1);
          sub(p2,v1,a);
          sub(p2,p1,v1);
          sub(p2,r2,v2);
          alpha = dot(v1,n1)/dot(v2,n1);
          scale(alpha,v2,v1);
          sub(p2,v1,b);
          return true;
        }
      } else {
        return false;
      }
    } else {
      sub(q2,p1,v2);
      cross(v1,v2,n);
      if (dot(v,n)<0.0) {
        return false;
      } else {
        sub(r1,p1,v1);
        cross(v1,v2,n);
        if (dot(v,n)>=0.0) {
          sub(p1,p2,v1);
          sub(p1,r1,v2);
          alpha = dot(v1,n2)/dot(v2,n2);
          scale(alpha,v2,v1);
          sub(p1,v1,a);
          sub(p1,p2,v1);
          sub(p1,q1,v2);
          alpha = dot(v1,n2)/dot(v2,n2);
          scale(alpha,v2,v1);
          sub(p1,v1,b);
          return true;
        } else {
          sub(p2,p1,v1);
          sub(p2,q2,v2);
          alpha = dot(v1,n1)/dot(v2,n1);
          scale(alpha,v2,v1);
          sub(p2,v1,a);
          sub(p1,p2,v1);
          sub(p1,q1,v2);
          alpha = dot(v1,n2)/dot(v2,n2);
          scale(alpha,v2,v1);
          sub(p1,v1,b);
          return true;
        }
      }
    }
  }

  private static boolean triTri(
    Point3 p1, Point3 q1, Point3 r1, 
    Point3 p2, Point3 q2, Point3 r2,
    double dp2, double dq2, double dr2,
    Vector3 v, Vector3 v1, Vector3 v2,
    Vector3 n, Vector3 n1, Vector3 n2,
    Point3 a, Point3 b)
  {
    if (dp2>0.0) {
       if (dq2>0.0) 
         return constructIntersection(p1,r1,q1,r2,p2,q2,v,v1,v2,n,n1,n2,a,b);
       else if (dr2>0.0) 
         return constructIntersection(p1,r1,q1,q2,r2,p2,v,v1,v2,n,n1,n2,a,b);
       else 
         return constructIntersection(p1,q1,r1,p2,q2,r2,v,v1,v2,n,n1,n2,a,b);
    } else if (dp2<0.0) {
      if (dq2<0.0) 
        return constructIntersection(p1,q1,r1,r2,p2,q2,v,v1,v2,n,n1,n2,a,b);
      else if (dr2<0.0) 
        return constructIntersection(p1,q1,r1,q2,r2,p2,v,v1,v2,n,n1,n2,a,b);
      else 
        return constructIntersection(p1,r1,q1,p2,q2,r2,v,v1,v2,n,n1,n2,a,b);
    } else {
      if (dq2<0.0) {
        if (dr2>=0.0) 
          return constructIntersection(p1,r1,q1,q2,r2,p2,v,v1,v2,n,n1,n2,a,b);
        else 
          return constructIntersection(p1,q1,r1,p2,q2,r2,v,v1,v2,n,n1,n2,a,b);
      }
      else if (dq2>0.0) {
        if (dr2>0.0) 
          return constructIntersection(p1,r1,q1,p2,q2,r2,v,v1,v2,n,n1,n2,a,b);
        else 
          return constructIntersection(p1,q1,r1,q2,r2,p2,v,v1,v2,n,n1,n2,a,b);
      }
      else  {
        if (dr2>0.0) 
          return constructIntersection(p1,q1,r1,r2,p2,q2,v,v1,v2,n,n1,n2,a,b);
        else if (dr2<0.0) 
          return constructIntersection(p1,r1,q1,r2,p2,q2,v,v1,v2,n,n1,n2,a,b);
        else
          return false;
      }
    }
  }
}
