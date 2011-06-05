/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
*****************************************************************************
Benchmark nested parallel loops with Scala parallel collections.
@author Dave Hale, Colorado School of Mines
@version 2011.06.04
****************************************************************************/
import java.lang.System
import scala.math._
import scala.util.Random

object ParallelNest {

  type Float1 = Array[Float]  // 1-D array of floats
  type Float2 = Array[Float1] // 2-D array of floats
  type Float3 = Array[Float2] // 3-D array of floats

  val maxtest = 3
  val maxtime = 3.0

  def main(args:Array[String]) = {
    println("Parallel benchmark in Scala:")
    val ns = List(
      (1000, 50000,     5),
      (1000,  5000,    50),
      (1000,   500,   500),
      (1000,    50,  5000),
      (1000,     5, 50000))
    for ((n1,n2,n3) <- ns)
      benchSolr(n1,n2,n3)
  }

  def benchSolr(n1:Int, n2:Int, n3:Int) = {
    println("Array solr: n1="+n1+" n2="+n2+" n3="+n3)
    val a1 = -1.8f
    val a2 = 0.81f
    val b0 = 2.0f
    val b1 = -3.2f
    val b2 = 1.28f
    val x = Array.fill(n3,n2,n1){1.0f}
    val ys = Array.fill(n3,n2,n1){1.0f}
    val yp = Array.fill(n3,n2,n1){1.0f}
    var s = new Stopwatch
    val mflop2 = 9.0e-6*n1*n2
    val mflop3 = 9.0e-6*n1*n2*n3
    var niter = 0
    for (ntest <- 0 until maxtest) {
      s.restart
      niter = 0; while (s.time<maxtime) {
        solrS(a1,a2,b0,b1,b2,x(niter%n3),ys(niter%n3))
        niter += 1
      }
      s.stop
      println("2D S: rate = "+(niter*mflop2/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        solrP(a1,a2,b0,b1,b2,x(niter%n3),yp(niter%n3))
        niter += 1
      }
      s.stop
      println("2D P: rate = "+(niter*mflop2/s.time).toInt)
      assertEqual(ys(0),yp(0))
      s.restart
      niter = 0; while (s.time<maxtime) {
        solrS(a1,a2,b0,b1,b2,x,ys)
        niter += 1
      }
      s.stop
      println("3D S: rate = "+(niter*mflop3/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        solrP(a1,a2,b0,b1,b2,x,yp)
        niter += 1
      }
      s.stop
      println("3D P: rate = "+(niter*mflop3/s.time).toInt)
      assertEqual(ys,yp)
    }
  }
  def solrS(
    a1:Float, a2:Float, b0:Float, b1:Float, b2:Float,
    x:Float1, y:Float1):Unit = 
  {
    val n = x.length
    var yim2 = 0.0f
    var yim1 = 0.0f
    var xim2 = 0.0f
    var xim1 = 0.0f
    var i = 0
    while (i<n) {
      val xi = x(i)
      val yi = b0*xi+b1*xim1+b2*xim2-a1*yim1-a2*yim2
      y(i) = yi
      yim2 = yim1
      yim1 = yi
      xim2 = xim1
      xim1 = xi
      i += 1
    }
  }
  private def solrS(
    a1:Float, a2:Float, b0:Float, b1:Float, b2:Float,
    x:Float2, y:Float2):Unit = 
  {
    x.indices.foreach(i=>solrS(a1,a2,b0,b1,b2,x(i),y(i)))
  }
  private def solrS(
    a1:Float, a2:Float, b0:Float, b1:Float, b2:Float,
    x:Float3, y:Float3):Unit = 
  {
    x.indices.foreach(i=>solrS(a1,a2,b0,b1,b2,x(i),y(i)))
  }
  private def solrP(
    a1:Float, a2:Float, b0:Float, b1:Float, b2:Float,
    x:Float2, y:Float2):Unit = 
  {
    if (x.length<10)
      x.indices.foreach(j=>solrS(a1,a2,b0,b1,b2,x(j),y(j)))
    else
      x.indices.par.foreach(j=>solrS(a1,a2,b0,b1,b2,x(j),y(j)))
  }
  private def solrP(
    a1:Float, a2:Float, b0:Float, b1:Float, b2:Float,
    x:Float3, y:Float3):Unit = 
  {
    if (x.length<10)
      x.indices.foreach(j=>solrP(a1,a2,b0,b1,b2,x(j),y(j)))
    else
      x.indices.par.foreach(j=>solrP(a1,a2,b0,b1,b2,x(j),y(j)))
  }

  def benchArraySqr(n1:Int, n2:Int, n3:Int) = {
    println("Array sqr: n1="+n1+" n2="+n2+" n3="+n3)
    val a = Array.fill(n3,n2,n1){1.0f}
    val bs = Array.fill(n3,n2,n1){0.0f}
    val bp = Array.fill(n3,n2,n1){0.0f}
    var s = new Stopwatch
    val mflop2 = 1.0e-6*n1*n2
    val mflop3 = 1.0e-6*n1*n2*n3
    var niter = 0
    for (ntest <- 0 until maxtest) {
      s.restart
      niter = 0; while (s.time<maxtime) {
        sqrS(a(niter%n3),bs(niter%n3))
        niter += 1
      }
      s.stop
      println("2D S: rate = "+(niter*mflop2/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        sqrP(a(niter%n3),bp(niter%n3))
        niter += 1
      }
      s.stop
      println("2D P: rate = "+(niter*mflop2/s.time).toInt)
      assertEqual(bs(0),bp(0))
      s.restart
      niter = 0; while (s.time<maxtime) {
        sqrS(a,bs)
        niter += 1
      }
      s.stop
      println("3D S: rate = "+(niter*mflop3/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        sqrP(a,bp)
        niter += 1
      }
      s.stop
      println("3D P: rate = "+(niter*mflop3/s.time).toInt)
      assertEqual(bs,bp)
    }
  }
  private def sqrS(a:Float1,b:Float1):Unit = {
    val n = a.length
    var i = 0; while (i<n) {
      b(i) = a(i)*a(i)
      i += 1
    }
  }
  private def sqrS(a:Float2, b:Float2):Unit = {
    for (i <- a.indices)
      sqrS(a(i),b(i))
  }
  private def sqrS(a:Float3, b:Float3):Unit = {
    for (i <- a.indices)
      sqrS(a(i),b(i))
  }
  private def sqrP(a:Float2, b:Float2):Unit = {
    for (i <- a.indices.par)
      sqrS(a(i),b(i))
  }
  private def sqrP(a:Float3, b:Float3):Unit = {
    for (i <- a.indices.par)
      sqrP(a(i),b(i))
  }

  def benchArraySum(n1:Int, n2:Int, n3:Int) = {
    println("Array sum: n1="+n1+" n2="+n2+" n3="+n3)
    val a = Array.fill(n3,n2,n1){1.0f}
    var s = new Stopwatch
    val mflop2 = 1.0e-6*n1*n2
    val mflop3 = 1.0e-6*n1*n2*n3
    var niter = 0
    for (ntest <- 0 until maxtest) {
      var ss = 0.0f
      var sp = 0.0f
      s.restart
      niter = 0; while (s.time<maxtime) {
        ss = sumS(a(niter%n3))
        niter += 1
      }
      s.stop
      println("2D S: rate = "+(niter*mflop2/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        sp = sumP(a(niter%n3))
        niter += 1
      }
      s.stop
      println("2D P: rate = "+(niter*mflop2/s.time).toInt)
      assert(abs(ss-sp)<0.0001*ss);
      ss = 0.0f
      sp = 0.0f
      s.restart
      niter = 0; while (s.time<maxtime) {
        ss = sumS(a)
        niter += 1
      }
      s.stop
      println("3D S: rate = "+(niter*mflop3/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        sp = sumP(a)
        niter += 1
      }
      s.stop
      println("3D P: rate = "+(niter*mflop3/s.time).toInt)
      assert(abs(ss-sp)<0.0001*ss);
    }
  }
  private def sumS(a:Float1):Float = {
    //a.sum // too slow
    val n = a.length
    var s = 0.0f
    var i = 0; while (i<n) {
      s += a(i)
      i += 1
    }
    s
  }
  private def sumS(a:Float2):Float = {
    a.aggregate(0.0f)(_+sumS(_),_+_)
  }
  private def sumS(a:Float3):Float = {
    a.aggregate(0.0f)(_+sumS(_),_+_)
  }
  private def sumP(a:Float2):Float = {
    //if (a.length<400) sumS(a) else a.par.aggregate(0.0f)(_+sumS(_),_+_)
    a.par.aggregate(0.0f)(_+sumS(_),_+_)
  }
  private def sumP(a:Float3):Float = {
    a.par.aggregate(0.0f)(_+sumP(_),_+_)
  }

  def benchArraySaxpy(n1:Int, n2:Int, n3:Int) = {
    println("Array saxpy: n1="+n1+" n2="+n2+" n3="+n3)
    val a = 3.14f;
    val x = Array.fill(n3,n2,n1){1.0f}
    val ys = Array.fill(n3,n2,n1){1.0f}
    val yp = Array.fill(n3,n2,n1){1.0f}
    var s = new Stopwatch
    val mflop2 = 2.0e-6*n1*n2
    val mflop3 = 2.0e-6*n1*n2*n3
    var niter = 0
    for (ntest <- 0 until maxtest) {
      s.restart
      niter = 0; while (s.time<maxtime) {
        saxpyS(a,x(niter%n3),ys(niter%n3))
        niter += 1
      }
      s.stop
      println("2D S: rate = "+(niter*mflop2/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        saxpyP(a,x(niter%n3),yp(niter%n3))
        niter += 1
      }
      s.stop
      println("2D P: rate = "+(niter*mflop2/s.time).toInt)
      //assertEqual(ys(0),yp(0))
      s.restart
      niter = 0; while (s.time<maxtime) {
        saxpyS(a,x,ys)
        niter += 1
      }
      s.stop
      println("3D S: rate = "+(niter*mflop3/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        saxpyP(a,x,yp)
        niter += 1
      }
      s.stop
      println("3D P: rate = "+(niter*mflop3/s.time).toInt)
      //assertEqual(ys,yp)
    }
  }
  private def saxpyS(a:Float, x:Float1, y:Float1):Unit = {
    val n = x.length
    var i = 0; while (i<n) {
      y(i) += a*x(i)
      i += 1
    }
  }
  private def saxpyS(a:Float, x:Float2, y:Float2):Unit = {
    for (i <- x.indices)
      saxpyS(a,x(i),y(i))
  }
  private def saxpyS(a:Float, x:Float3, y:Float3):Unit = {
    for (i <- x.indices)
      saxpyS(a,x(i),y(i))
  }
  private def saxpyP(a:Float, x:Float2, y:Float2):Unit = {
    for (i <- x.indices.par)
      saxpyS(a,x(i),y(i))
  }
  private def saxpyP(a:Float, x:Float3, y:Float3):Unit = {
    for (i <- x.indices.par)
      saxpyP(a,x(i),y(i))
  }

  def benchArraySdot(n1:Int, n2:Int, n3:Int) = {
    println("Array sdot: n1="+n1+" n2="+n2+" n3="+n3)
    val x = Array.fill(n3,n2,n1){1.0f}
    val y = Array.fill(n3,n2,n1){2.0f}
    var s = new Stopwatch
    val mflop2 = 2.0e-6*n1*n2
    val mflop3 = 2.0e-6*n1*n2*n3
    var niter = 0
    for (ntest <- 0 until maxtest) {
      var ss = 0.0f
      var sp = 0.0f
      s.restart
      niter = 0; while (s.time<maxtime) {
        ss = sdotS(x(niter%n3),y(niter%n3))
        niter += 1
      }
      s.stop
      println("2D S: rate = "+(niter*mflop2/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        sp = sdotP(x(niter%n3),y(niter%n3))
        niter += 1
      }
      s.stop
      println("2D P: rate = "+(niter*mflop2/s.time).toInt)
      assert(abs(ss-sp)<0.0001*ss);
      ss = 0.0f
      sp = 0.0f
      s.restart
      niter = 0; while (s.time<maxtime) {
        ss = sdotS(x,y)
        niter += 1
      }
      s.stop
      println("3D S: rate = "+(niter*mflop3/s.time).toInt)
      s.restart
      niter = 0; while (s.time<maxtime) {
        sp = sdotP(x,y)
        niter += 1
      }
      s.stop
      println("3D P: rate = "+(niter*mflop3/s.time).toInt)
      assert(abs(ss-sp)<0.0001*ss);
    }
  }
  private def sdotS(x:Float1, y:Float1):Float = {
    val n = x.length
    var s = 0.0f
    var i = 0; while (i<n) {
      s += x(i)*y(i)
      i += 1
    }
    s
  }
  private def sdotS(x:Float2, y:Float2):Float = {
    x.indices.aggregate(0.0f)((s,i)=>s+sdotS(x(i),y(i)),_+_)
  }
  private def sdotS(x:Float3, y:Float3):Float = {
    x.indices.aggregate(0.0f)((s,i)=>s+sdotS(x(i),y(i)),_+_)
  }
  private def sdotP(x:Float2, y:Float2):Float = {
    x.indices.par.aggregate(0.0f)((s,i)=>s+sdotS(x(i),y(i)),_+_)
  }
  private def sdotP(x:Float3, y:Float3):Float = {
    x.indices.par.aggregate(0.0f)((s,i)=>s+sdotP(x(i),y(i)),_+_)
  }

  private def assertEqual(a:Float1, b:Float1):Unit = {
    assert(a.length==b.length,"same dimensions")
    val n = a.length
    var i = 0; while (i<n) {
      assert(a(i)==b(i),"same elements")
      i += 1
    }
  }
  private def assertEqual(a:Float2, b:Float2):Unit = {
    assert(a.length==b.length,"same dimensions")
    val n = a.length
    var i = 0; while (i<n) {
      assertEqual(a(i),b(i))
      i += 1
    }
  }
  private def assertEqual(a:Float3, b:Float3):Unit = {
    assert(a.length==b.length,"same dimensions")
    val n = a.length
    var i = 0; while (i<n) {
      assertEqual(a(i),b(i))
      i += 1
    }
  }

  private class Stopwatch {
    def start:Unit = {
      if (!_running) {
        _running = true
        _start = System.nanoTime
      }
    }
    def stop:Unit = {
      if (_running) {
        _time += System.nanoTime-_start
        _running = false
      }
    }
    def reset:Unit = {
      stop
      _time = 0
    }
    def restart:Unit = {
      reset
      start
    }
    def time:Double = {
      if (_running)
        1.0e-9*(_time+System.nanoTime-_start)
      else
        1.0e-9*_time
    }
    private var _running:Boolean = false
    private var _start:Long = 0L
    private var _time:Long = 0L
  }
}
