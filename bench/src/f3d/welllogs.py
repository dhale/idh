"""
Reads, reformats, displays well logs for F3.
Author: Dave Hale, Colorado School of Mines
Version: 2012.12.29
"""
from f3utils import *

#############################################################################
setupForSubset("all")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
ss = Sampling(1176,d2,0.0)
ns,ds,fs = ss.count,ss.delta,ss.first
f3dDataDir = getF3dDataDir()
odtWellLogsDir = f3dDataDir+"odt/"
csmWellLogsDir = f3dDataDir
odtWellLogs = odtWellLogsDir
csmWellLogs = csmWellLogsDir+getF3dBaseName()+"well.dat"
clip = 1.0

ltype,lmin,lmax,llabel = None,None,None,None
curvePars = {
  "v":(1.9,2.4,"Velocity (km/s)"),
  "d":(2.0,2.3,"Density (gm/cc)"),
  "g":(20.0,90.0,"Gamma ray (API)"),
  "p":(0.25,0.40,"Porosity")
}
def setupForLogCurve(curve):
  global ltype,lmin,lmax,llabel
  ltype = curve
  lmin,lmax,llabel = curvePars[curve]
setupForLogCurve("g")

#############################################################################
def main(args):
  #makeBinaryWellLogs()
  #makeSliceThruWells()
  displaySlicesThruWells()
  #displayWithWells(ltype,lmin,lmax)
  #viewWellCurves("velocity")
  #viewWellCurves("density")
  #viewWellCurves("gamma")
  #viewWellCurves("porosity")

def wellLogImagesOnCurve(ltype,lmin,lmax,x2s,x3s):
  def clipNonNull(fnull,fmin,fmax,f):
    n = len(f)
    f = copy(f)
    for i in range(n):
      if f[i]!=fnull:
        if f[i]<fmin: f[i] = fmin
        if f[i]>fmax: f[i] = fmax
    return f
  def replaceNull(fnull,fval,f):
    n = len(f)
    f = copy(f)
    for i in range(n):
      if f[i]==fnull:
        f[i] = fval
    return f
  wldata = readWellLogData()
  wlis = []
  fnull = -Float.MAX_VALUE
  fmin = lmin+(lmax-lmin)/256
  fmax = lmax
  for log in wldata.getAll():
    f = log.resample(ltype,fnull,s1)
    f = clipNonNull(fnull,fmin,fmax,f)
    f = replaceNull(fnull,lmin,f)
    g = zerofloat(n1,2)
    copy(f,g[0])
    copy(f,g[1])
    i = indexOfNearestPoint(log.x2[0],log.x3[0],x2s,x3s)
    sw = Sampling(2,4*ds,ss.getValue(i)-2*ds)
    wlis.append((s1,sw,g))
  return wlis

def displaySlicesThruWells():
  fsw = readImage2("f3dsw",n1,ns)
  fk1 = readImage2(getF3dSlice1Name(309),n2,n3)
  x2w,x3w = getWellLocations()
  x2s,x3s = getCurveThruPoints(x2w,x3w)
  sp = SimplePlot()
  pv = sp.addPixels(s2,s3,fk1)
  pv.setClips(-clip,clip)
  pv = sp.addPoints(x2s,x3s)
  pv.setLineColor(Color.YELLOW)
  pv.setLineWidth(3.0)
  pv = sp.addPoints(x2w,x3w)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(16.0)
  pv.setMarkColor(Color.YELLOW)
  #sp.setHLimits(s2.first,s2.last)
  #sp.setVLimits(s3.first,s3.last)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  sp.setSize(1100,800)
  #sp.paintToPng(300,10.0,"png/sw23.png")
  alpha = fillfloat(1.0,256); alpha[0] = 0.0
  cmap = ColorMap.setAlpha(ColorMap.JET,alpha)
  for curve in ["v","d","g"]:
    setupForLogCurve(curve)
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPixels(s1,ss,fsw)
    pv.setClips(-clip,clip)
    wlis = wellLogImagesOnCurve(ltype,lmin,lmax,x2s,x3s)
    for wli in wlis:
      s1w,s2w,fw = wli
      pv = sp.addPixels(s1w,s2w,fw)
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv.setClips(lmin,lmax)
      pv.setColorModel(cmap)
    sp.addColorBar(llabel)
    sp.plotPanel.setColorBarWidthMinimum(100)
    sp.setVLimits(0.0,s1.last)
    sp.setHLimits(ss.first,ss.last)
    sp.setHLabel("Distance (km)")
    sp.setVLabel("Time (s)")
    sp.setSize(1400,800)
    #sp.paintToPng(300,10.0,"png/sw1"+curve+".png")

def displayWithWells(ltype,lmin,lmax):
  world = World()
  x = readF3dImage()
  addImageToWorld(world,x,cmin=-clip,cmax=clip)
  addLogsToWorld(world,ltype,cmin=lmin,cmax=lmax)
  makeFrame(world)

# Gets well locations in CSM (x2,x3) coordinates
def getWellLocations():
  wldata = readWellLogData()
  #wldata.printInfo()
  x2,x3 = [],[]
  for log in wldata.getAll():
    print log.name
    x2.append(log.x2[0])
    x3.append(log.x3[0])
  nw = len(x2)
  if nw==4: # reorder for curve through wells
    x2t = x2[1]; x2[1] = x2[3]; x2[3] = x2t;
    x3t = x3[1]; x3[1] = x3[3]; x3[3] = x3t;
  ns = nw
  x2s = zerofloat(nw); copy(x2,x2s)
  x3s = zerofloat(nw); copy(x3,x3s)
  return x2s,x3s

# Gets finely sampled curve through points (well locations).
def getCurveThruPoints(x2s,x3s):
  ns = len(x2s)
  for i in range(2): # 2 iterations should be sufficient
    ds = zerofloat(ns)
    ds[0] = 0.0
    for js in range(1,ns):
      ds[js] = ds[js-1]+hypot(x2s[js]-x2s[js-1],x3s[js]-x3s[js-1])
    ci2 = CubicInterpolator(ds,x2s)
    ci3 = CubicInterpolator(ds,x3s)
    smin,smax = ds[0],ds[-1]
    smin -= 0.250
    smax += 0.250
    ns = 1+int((smax-smin)/s2.delta)
    ds = (smax-smin)/(ns-1)
    sj = rampfloat(smin,ds,ns)
    x2s = zerofloat(ns)
    x3s = zerofloat(ns)
    ci2.interpolate(sj,x2s)
    ci3.interpolate(sj,x3s)
  return x2s,x3s

# Gets 2D seismic image along specified curve.
def getImageAlongCurve(x2s,x3s):
  f = readF3dImage()
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

def indexOfNearestPoint(x2,x3,x2s,x3s):
  def distanceSquared(xa,ya,xb,yb):
    dx = xb-xa
    dy = yb-ya
    return dx*dx+dy*dy
  ns = len(x2s)
  jsmin = ns
  dsmin = Double.MAX_VALUE
  for js in range(ns):
    ds = distanceSquared(x2,x3,x2s[js],x3s[js])
    if ds<dsmin:
      dsmin = ds
      jsmin = js
  return jsmin

def makeBinaryWellLogs():
  wldata = WellLog.Data(odtWellLogs)
  print "well log data"
  wldata.printInfo()
  #wldata.clip(s2,s3)
  #print "after clipping"
  #wldata.printInfo()
  wldata.writeBinary(csmWellLogs)
  
# Writes a 2D image slice through wells
def makeSliceThruWells():
  x2w,x3w = getWellLocations()
  x2s,x3s = getCurveThruPoints(x2w,x3w)
  f,fs = getImageAlongCurve(x2s,x3s)
  print "slice through well has",len(fs),"traces"
  writeImage2("f3dsw",fs)

def viewWellCurves(curve):
  wldata = WellLog.Data.readBinary(csmWellLogs)
  flist,zlist = [],[]
  colors = [Color.RED,Color.GREEN,Color.BLUE,Color.MAGENTA]
  icolor = 0
  sp = SimplePlot()
  sp.setHLabel("Time (s)")
  if curve=="velocity":
    sp.setVLabel("Velocity (km/s)")
  elif curve=="density":
    sp.setVLabel("Density (gm/cc)")
  elif curve=="gamma":
    sp.setVLabel("Gamma Ray")
  elif curve=="porosity":
    sp.setVLabel("Porosity")
  for log in wldata.getLogsWith(curve):
    #log.despike(3)
    print log.name,"",colors[icolor]
    log.smooth(25)
    f,z,y,x = log.getSamples(curve)
    tv = PointsView(z,f)
    tv.setLineColor(colors[icolor])
    icolor = (icolor+1)%len(colors)
    sp.add(tv)

def viewPointsPixels(x,y,z,sx,sy,zs,zlabel):
  sp = SimplePlot()
  sp.setHLabel("X (km)")
  sp.setVLabel("Y (km)")
  pp = sp.getPlotPanel()
  pp.setColorBarWidthMinimum(100)
  pv = sp.addPixels(sx,sy,zs)
  pv.setColorModel(ColorMap.JET)
  pv = sp.addPoints(x,y)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(3)
  sp.addColorBar(zlabel)
  sp.setSize(970,440)

#############################################################################
run(main)
