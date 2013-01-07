"""
Displays F3 data.
"""
from f3utils import *

#############################################################################
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
#ltype="v"; lmin=1.9; lmax=2.4
#ltype="d"; lmin=2.0; lmax=2.3
ltype="g"; lmin=20.0; lmax=90.0
#ltype="p"; lmin=0.25; lmax=0.4
clip = 1.0

#############################################################################
def main(args):
  display2DThruWells()
  #displayWithWells(ltype,lmin,lmax)
  #display3d()
  #displaySlice3()

def wellLocations():
  wldata = WellLog.Data.readBinary(wellLogsDir+"f3dwell.dat")
  wldata.printInfo()
  x2,x3 = [],[]
  for log in wldata.getAll():
    print log.name
    x2.append(log.x2[0])
    x3.append(log.x3[0])
  nw = len(x2) # assume nw = 4; make last well first
  #x2t = x2[3]; x2[3] = x2[2]; x2[2] = x2[1]; x2[1] = x2[0]; x2[0] = x2t
  #x3t = x3[3]; x3[3] = x3[2]; x3[2] = x3[1]; x3[1] = x3[0]; x3[0] = x3t
  x2t = x2[1]; x2[1] = x2[3]; x2[3] = x2t;
  x3t = x3[1]; x3[1] = x3[3]; x3[3] = x3t;
  ns = nw
  x2s = zerofloat(nw); copy(x2,x2s)
  x3s = zerofloat(nw); copy(x3,x3s)
  return x2s,x3s

def curveThruPoints(x2s,x3s):
  ns = len(x2s); print "ns =",ns
  for i in range(2):
    ds = zerofloat(ns)
    ds[0] = 0.0
    for js in range(1,ns):
      ds[js] = ds[js-1]+hypot(x2s[js]-x2s[js-1],x3s[js]-x3s[js-1])
    ci2 = CubicInterpolator(ds,x2s)
    ci3 = CubicInterpolator(ds,x3s)
    smin,smax = ds[0],ds[-1]
    smin -= 0.250
    smax += 0.250
    ns = 1+int((smax-smin)/s2.delta); print "ns =",ns
    ds = (smax-smin)/(ns-1)
    sj = rampfloat(smin,ds,ns)
    x2s = zerofloat(ns)
    x3s = zerofloat(ns)
    ci2.interpolate(sj,x2s)
    ci3.interpolate(sj,x3s)
  return x2s,x3s

def getImageAlongCurve(x2s,x3s):
  setupForSubset("all")
  s1,s2,s3 = getSamplings()
  f = readImage("f3dall")
  ns = len(x2s)
  g = zerofloat(n1,ns)
  si = SincInterp()
  for js in range(ns):
    x2j = x2s[js]
    x3j = x3s[js]
    for j1 in range(n1):
      x1j = s1.getValue(j1)
      g[js][j1] = si.interpolate(s1,s2,s3,f,x1j,x2j,x3j)
  return f,g
  
def display2DThruWells():
  setupForSubset("all")
  s1,s2,s3 = getSamplings()
  x2w,x3w = wellLocations()
  x2s,x3s = curveThruPoints(x2w,x3w)
  f,fs = getImageAlongCurve(x2s,x3s)
  n2,n3 = s2.count,s3.count
  fj = zerofloat(n2,n3)
  j1 = int((1.244-s1.first)/s1.delta+0.5)
  for i3 in range(n3):
    for i2 in range(n2):
      fj[i3][i2] = f[i3][i2][j1]
  sp = SimplePlot()
  pv = sp.addPixels(s2,s3,fj)
  pv.setClips(-clip,clip)
  pv = sp.addPoints(x2s,x3s)
  pv.setLineColor(Color.YELLOW)
  pv = sp.addPoints(x2w,x3w)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkColor(Color.YELLOW)
  #sp.setHLimits(s2.first,s2.last)
  #sp.setVLimits(s3.first,s3.last)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  sp.paintToPng(300,3.3,"png/curvy23.png")
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  ss = Sampling(len(fs),s2.delta,0.0)
  pv = sp.addPixels(s1,ss,fs)
  pv.setClips(-clip,clip)
  sp.setVLimits(0.3,s1.last)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  sp.paintToPng(300,3.3,"png/curvy12.png")

def displayWithWells(ltype,lmin,lmax):
  world = World()
  setupForSubset("all")
  x = readImage("f3dall")
  #x = readImage("f3d")
  addImageToWorld(world,x,cmin=-clip,cmax=clip)
  addLogsToWorld(world,ltype,cmin=lmin,cmax=lmax)
  makeFrame(world)

def display3d():
  x = readImage("f3d")
  print "x min =",min(x)," max =",max(x)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,x)
  ipg.setClips(-clip,clip)
  view = frame.getOrbitView()
  view.setScale(3.0)
  view.setAxesScale(1.0,1.0,8.0)
  frame.setSize(1200,900)
  frame.setVisible(True)

def displaySlice3():
  k3 = 75
  x = readImage("f3d")
  y = x[k3]
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(y)
  pv.setClips(-clip,clip)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
