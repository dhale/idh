"""
Displays images through conical faults in F3.
Author: Dave Hale, Colorado School of Mines
Version: 2013.03.02
"""
from f3utils import *

#############################################################################
setupForDataSet("seta")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
f3dDataDir = getF3dDataDir()
clip = 1.0

pngDir = "png/"
#pngDir = None

# P-wave velocities for sediments with cones ~ 2.2 km/s => dz ~ 1.1*dt 
# So, for 1:1 plots, make pixel_ratio = 1.1*value_ratio
# for vertical exaggeration by s, scale pixel_ratio by 1/s

#############################################################################
def main(args):
  goFigures()

def goFigures():
  for name in ["g","gs8"]:
    #makeSliceThruPoints(name,"c1")
    #displaySlicesThruCones(name,"c1")
    displayCones(name)
  
def displayCones(name):
  g = readImage(name)
  nt = 126 # sampling of times with cones present
  dt = 0.004
  ft = 1.040
  st = Sampling(nt,dt,ft)
  nd = 73  # sampling of distance d from center of cone
  dd = s2.delta
  fd = -dd*(nd-1)/2
  sd = Sampling(nd,dd,fd)
  na = 6 # sampling of azimuth a; north = 0 degrees
  da = 30.0
  fa = 0.0
  sa = Sampling(na,da,fa)
  c2s,c3s = getConeLocations()
  for ic in [0,1,2,3]:
  #for ic in [0]:
    c2,c3 = c2s[ic],c3s[ic]
    h = makeConeSlices(c2,c3,st,sd,sa,g)
    def plota(ia,aspect11):
      sia = str(int(sa.getValue(ia)))
      sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
      #sp.setTitle("Azimuth = "+sia+" degrees")
      sp.setHLabel("Radial distance (km)")
      sp.setVLabel("Time (s)")
      #sp.setHLimits(sd.first,sd.last)
      #sp.setVLimits(1.05,1.55)
      if aspect11:
        sp.setSize(1400,580) # 1:1 scale
      else:
        #sp.setSize(590,765) # 1:5 scale (vertical exaggeration with title)
        sp.setSize(590,712) # 1:5 scale (vertical exaggeration)
      wpt = 200.0
      sp.setFontSizeForPrint(8,wpt)
      pv = sp.addPixels(st,sd,h[ia])
      pv.setClips(-1.0,1.0)
      pv = sp.addPoints([st.first,st.last],[0,0])
      pv.setLineColor(Color.WHITE)
      pngFile="cone"+str(ic)+name
      if len(sia)==1:
        sia = "00"+sia
      elif len(sia)==2:
        sia = "0"+sia
      if aspect11:
        pngFile += "asa"+sia+".png"
      else:
        pngFile += "asb"+sia+".png"
      if pngDir:
        sp.paintToPng(360,wpt/72,pngDir+pngFile)
    for ia in range(na):
      plota(ia,False)

def displaySlicesThruCones(name,points):
  ss = getSamplingS(name,points)
  ns,ds,fs = ss.count,ss.delta,ss.first
  fcp = readImage2("slices/"+name+points,n1,ns)
  fk1 = readImage2(getF3dSlice1Name(name,309),n2,n3)
  x2c,x3c = getConeLocations()
  x2p,x3p = getPointSet(points)
  x2s,x3s = getCurveThruPoints(x2p,x3p)
  sp = SimplePlot()
  pv = sp.addPixels(s2,s3,fk1)
  pv.setClips(-clip,clip)
  pv = sp.addPoints(x2s,x3s)
  pv.setLineColor(Color.WHITE)
  pv.setLineWidth(1.0)
  pv = sp.addPoints(x2c,x3c)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  pv.setMarkSize(64.0)
  pv.setMarkColor(Color.WHITE)
  #sp.setHLimits(s2.first,s2.last)
  #sp.setVLimits(s3.first,s3.last)
  sp.setHLimits(0.0,9.0)
  sp.setVLimits(0.0,15.0)
  sp.setHInterval(2.0)
  sp.setVInterval(2.0)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  wpt = 240.0
  sp.setFontSizeForPrint(8,wpt)
  sp.setSize(460,710)
  if pngDir:
    sp.paintToPng(360,wpt/72,pngDir+name+points+"23.png")
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(s1,ss,fcp)
  pv.setClips(-clip,clip)
  x1,xs = linesThruCones(x2c,x3c,x2s,x3s,ss)
  pv = sp.plotPanel.addPoints(x1,xs)
  pv.setLineColor(Color.WHITE)
  sp.setVLimits(1.0,s1.last)
  sp.setHLimits(ss.first,ss.last)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  wpt = 504.0
  sp.setFontSizeForPrint(8,wpt)
  #sp.setSize(775,440)
  sp.setSize(775,255)
  if pngDir:
    sp.paintToPng(360,wpt/72,pngDir+name+points+"1s.png")

# Gets cone locations in CSM (x2,x3) coordinates
def getConeLocations():
  #x2c = [ 0.910, 2.387, 4.533, 6.225, 7.002]
  #x3c = [11.896, 6.229, 6.140, 1.678, 1.013]
  x2c = [ 0.910, 2.387, 4.533, 7.002]
  x3c = [11.896, 6.229, 6.140, 1.013]
  return x2c,x3c

def getPointSet(points):
  x2c,x3c = getConeLocations()
  if (points=="c1"):
    x2p = [x2c[0],x2c[1],x2c[3]]
    x3p = [x3c[0],x3c[1],x3c[3]]
  return x2p,x3p

# Gets finely sampled curve through points.
def getCurveThruPoints(x2s,x3s):
  ns = len(x2s)
  for i in range(2): # 2 iterations should be sufficient
    ds = zerofloat(ns)
    ds[0] = 0.0
    for js in range(1,ns):
      ds[js] = ds[js-1]+hypot(x2s[js]-x2s[js-1],x3s[js]-x3s[js-1])
    method = CubicInterpolator.Method.SPLINE
    ci2 = CubicInterpolator(method,ds,x2s)
    ci3 = CubicInterpolator(method,ds,x3s)
    smin,smax = ds[0],ds[-1]
    smin -= 0.500
    smax += 0.500
    ns = 1+int((smax-smin)/s2.delta)
    ds = (smax-smin)/(ns-1)
    sj = rampfloat(smin,ds,ns)
    x2s = zerofloat(ns)
    x3s = zerofloat(ns)
    ci2.interpolate(sj,x2s)
    ci3.interpolate(sj,x3s)
  return x2s,x3s

# Gets 2D seismic image along specified curve.
def getImageAlongCurve(name,x2s,x3s):
  f = readImage(name)
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

# Writes a 2D image slice that passes through point set
def makeSliceThruPoints(name,points):
  x2p,x3p = getPointSet(points)
  x2s,x3s = getCurveThruPoints(x2p,x3p)
  f,fs = getImageAlongCurve(name,x2s,x3s)
  print "slice through points has",len(fs),"traces"
  writeImage2("slices/"+name+points,fs)

# Sampling of 2nd dimension for slices
def getSamplingS(name,points):
  f = File(getF3dDataSetDir()+"slices/"+name+points+".dat")
  n2 = f.length()/n1/4
  return Sampling(n2,d2,0.0)

def linesThruCones(x2c,x3c,x2s,x3s,ss):
  nc = len(x2c)
  x1,xs = [],[]
  for ic in range(nc):
    i = indexOfNearestPoint(x2c[ic],x3c[ic],x2s,x3s)
    if i>=0:
      xi = ss.getValue(i)
      x1.append([1.04,1.54])
      xs.append([xi,xi])
  return x1,xs

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
  if dsmin>d2:
    jsmin = -1
  return jsmin

def makeConeSlices(c2,c3,st,sd,sa,g):
  nt,nd,na = st.count,sd.count,sa.count
  dt,dd,da = st.delta,sd.delta,sa.delta
  ft,fd,fa = st.first,sd.first,sa.first
  # output array is a 3D image[na][nd][n1]
  h = zerofloat(nt,nd,na)
  # survey azimuth is 88.4 degrees
  ps = toRadians(88.4)
  # interpolate slices
  si = SincInterp()
  for ja in range(na):
    aj = sa.getValue(ja)
    pj = toRadians(aj)
    cj = cos(pj-ps)
    sj = sin(pj-ps)
    for jd in range(nd):
      dj = sd.getValue(jd)
      x2 = c2+cj*dj;
      x3 = c3-sj*dj;
      for jt in range(nt):
        x1 = st.getValue(jt)
        h[ja][jd][jt] = si.interpolate(s1,s2,s3,g,x1,x2,x3)
  return h

#############################################################################
run(main)
