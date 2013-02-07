#############################################################################
# Dynamic warping for 2D images

from imports import *
from warp import DynamicWarpingS

#############################################################################

def main(args):
  f,g = getSinoImages()
  nt,nx = len(f[0]),len(f)
  ni,nl = 501,351
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  SimplePlot.asPixels(f)
  SimplePlot.asPixels(g)
  ht = 50
  for mx in [1,5,10,20,40,80]:
    ixl = nx/2-mx/2
    ixu = ixl+mx
    nj = 1+(ni-1)/ht
    k = rampint(0,ht,nj)
    e = zerofloat(nl,ni)
    for ix in range(ixl,ixu):
      ei = DynamicWarpingS.computeErrors(nl,f[ix],g[ix])
      add(ei,e,e)
    u1 = DynamicWarpingS.findShifts(0.0,2.0,e)
    uj = DynamicWarpingS.findShiftsK(0.0,2.0,e,k)
    tj = rampfloat(0.0,ht,nj)
    ti = rampfloat(0.0,1.0,ni)
    ci = CubicInterpolator(CubicInterpolator.Method.MONOTONIC,tj,uj)
    u2 = ci.interpolate(ti)
    plote(e)
    plote(e,u1,u2,tj,uj)

def plote(e,ua=None,ub=None,tj=None,uj=None):
  nl,ni = len(e[0]),len(e)
  sp = SimplePlot()
  pv = sp.addPixels(pow(transpose(e),1.0))
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  pv.setPercentiles(2,98)
  if ua:
    pv = sp.addPoints(ua)
    pv.setLineColor(Color.WHITE)
    pv.setLineStyle(PointsView.Line.DASH)
    pv.setLineWidth(3)
  if ub:
    pv = sp.addPoints(ub)
    pv.setLineColor(Color.WHITE)
    pv.setLineWidth(3)
  if tj and uj:
    pv = sp.addPoints(tj,uj)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.WHITE)
  sp.setLimits(0,0,ni-1,nl-1)
  sp.setSize(800,585)

def getSinoImages():
  dataDir = "/data/seis/sino/"
  n1,d1,f1 = 2001,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"z260.dat",n1,n2)
  g = readImage(dataDir+"x260.dat",n1,n2)
  gain(100,f)
  gain(100,g)
  return f,g

#############################################################################
# utilities

def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
