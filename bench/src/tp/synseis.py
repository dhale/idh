"""
Synthetic seismograms for Teapot Dome.
Author: Dave Hale, Colorado School of Mines
Version: 2013.03.28
"""
from imports import *
from model import *
from warp import *
from dnp import *

#############################################################################
def main(args):
  wellId = 490252305400
  #goModel(wellId)
  #goSimple(wellId)
  goTie(wellId)

def goSeis(wellId):
  setGlobals("all")
  wldata = WellLog.Data.readBinary(csmWellLogs)
  log = wldata.get(wellId)
  model = SynSeis.getModel(log)
  ds = 0.002
  fs = model.tmin()
  ls = model.tmax()
  fs = int(fs/ds)*ds
  ls = int(ls/ds)*ds
  ns = 1+int((ls-fs)/ds)
  fpeak = 35.0
  q = 100.0
  ss = Sampling(ns,ds,fs)
  fs = SynSeis.makeSimpleSeismogram(model,fpeak,ss)
  fb = SynSeis.makeBetterSeismogram(model,q,fpeak,ss)
  def plot(s,f):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPoints(s,f)
    sp.setSize(300,800)
  plot(ss,fs)
  plot(ss,fb)

def goTie(wellId):
  setGlobals("all")
  wldata = WellLog.Data.readBinary(csmWellLogs)
  log = wldata.get(wellId)
  model = SynSeis.getModel(log)
  ds = 0.002
  fs = model.tmin()
  ls = model.tmax()
  fs = int(fs/ds)*ds
  ls = int(ls/ds)*ds
  ns = 1+int((ls-fs)/ds)
  fpeak = 35.0
  q = 100.0
  ss = Sampling(ns,ds,fs)
  #f = SynSeis.makeSimpleSeismogram(model,fpeak,ss)
  f = SynSeis.makeBetterSeismogram(model,q,fpeak,ss)
  f = normalize(f)
  n1,n2,n3 = s1.count,s2.count,s3.count
  ga = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(csmSeismict)
  ais.readFloats(ga)
  ais.close()
  i2 = s2.indexOfNearest(model.x2)
  i3 = s3.indexOfNearest(model.x3)
  x2 = s2.getValue(i2)
  x3 = s3.getValue(i3)
  print "i2 =",i2," i3 =",i3
  print "x2 =",x2," x3 =",x3
  #phaseScan(ss,f,i2,i3,ga); return
  #f = phaseShift(f,150.0) # simple
  f = phaseShift(f,250.0) # better, q = 80 to 150
  #f = phaseShift(f,160.0) # better, q = infinity
  def plot(s,f):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPoints(s,f)
    sp.setSize(300,800)
    sp.setVLimits(0.1,1.1)
  g = ga[i3][i2]
  g = normalize(g)
  esum,e,u = findShifts(ss,f,i2,i3,ga)
  st,r = invertShifts(ss,u)
  nt,dt,ft = st.count,st.delta,st.first
  s = sub(rampfloat(ft,dt,nt),r)
  si = SincInterpolator()
  h = zerofloat(nt)
  si.interpolate(ns,ds,fs,f,nt,s,h)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(ss,f)
  pv.setLineColor(Color.RED)
  sp.setVLimits(0.1,1.1)
  sp.setSize(300,800)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(s1,g)
  pv = sp.addPoints(st,h)
  pv.setLineColor(Color.RED)
  sp.setVLimits(0.1,1.1)
  sp.setSize(300,800)
  plotSynOnSeis(x2,x3,ss,f)
  plotSynOnSeis(x2,x3,st,h)

def plotSynOnSeis(x2,x3,ss,f):
  ns = ss.count
  n1,n2,n3 = s1.count,s2.count,s3.count
  i2 = s2.indexOfNearest(x2)
  i3 = s3.indexOfNearest(x3)
  ga = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(csmSeismict)
  ais.readFloats(ga)
  ais.close()
  g = ga[i3]
  for j2 in range(i2-1,i2+2):
    for js in range(ns):
      tj = ss.getValue(js)
      j1 = s1.indexOfNearest(tj)
      g[j2][j1] = f[js]
  for i2 in range(n2):
    g[i2] = normalize(g[i2])
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  sp.setHLimits(3.0,4.4)
  sp.setVLimits(0.1,1.1)
  sp.setSize(720,800)

def findShifts(sf,f,i2,i3,ga):
  umin,umax = -0.10,0.20
  rmin,rmax = -0.02,0.02
  dmin = 0.1
  dw = DynamicWarpingR(umin,umax,sf)
  dw.setStrainLimits(rmin,rmax)
  dw.setSmoothness(dmin)
  ss = dw.samplingS
  e = zerofloat(ss.count,sf.count)
  n2,n3 = s2.count,s3.count
  hw23 = 2
  for j3 in range(i3-2,i3+3):
    if j3<0 or j3>=n3: continue
    for j2 in range(i2-hw23,i2+hw23+1):
      if j3<0 or j3>=n3: continue
      g = normalize(ga[j3][j2])
      e = add(e,dw.computeErrors(sf,f,s1,g))
  u = dw.findShifts(e)
  esum = 0.0
  si = SincInterpolator()
  for jf in range(sf.count):
    esum += si.interpolate(ss,e[jf],u[jf])
  plot = True
  if plot:
    sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
    pv = sp.addPixels(sf,ss,pow(transpose(e),0.5))
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv = sp.addPoints(sf,u)
    pv.setLineColor(Color.WHITE)
  return esum,e,u

def phaseScan(sf,f,i2,i3,ga):
  sp = Sampling(36,10.0,0.0)
  ep = []
  for ip in range(sp.count):
    p = sp.getValue(ip)
    fp = phaseShift(f,p)
    esum,e,u = findShifts(sf,fp,i2,i3,ga)
    ep.append(esum)
    print "phase =",p," error =",esum
  SimplePlot.asPoints(sp,ep)

def phaseShift(f,phase):
  """
  Applies a constant-phase shift, for phase specified in degrees.
  """
  n = len(f)
  g = zerofloat(n)
  htf = HilbertTransformFilter()
  htf.apply(n,f,g)
  cp = cos(toRadians(phase))
  sp = sin(toRadians(phase))
  return add(mul(cp,f),mul(sp,g))

def invertShifts(ss,u):
  """ 
  Given a uniformly sampled u(s) such that s+u(s) increases monotonically, 
  computes a uniformly sampled r(t) = u(t-r(t)) such that t-r(t) increases 
  monotonically. Uses the following 3-step algorithm:
  (1) computes t(s) = s+u(s)
  (2) computes s(t) = by inverse linear interpolation of t(s)
  (3) computes r(t) = t-s(t)
  Returns the sampling of time st and the sequence of sampled r(t).
  """
  ns,ds,fs,ls = ss.count,ss.delta,ss.first,ss.last
  dt = ds # make sampling intervals equal
  ft = fs+u[0] # ft = time t of first sample
  lt = ls+u[ns-1] # lt = time of last sample
  ft = int(ft/dt)*dt # force ft to be a multiple of interval dt
  nt = 1+int((lt-ft)/dt) # number of t samples
  st = Sampling(nt,dt,ft) # sampling of t
  t = add(rampfloat(fs,ds,ns),u) # t(s) = s+u(s)
  s = zerofloat(nt)
  ii = InverseInterpolator(ss,st)
  ii.invert(t,s) # both t(s) and s(t) increase monotonically
  r = sub(rampfloat(ft,dt,nt),s) # r(t) = t-s(t)
  return st,r

def normalize(f):
  sigma = 100
  n = len(f)
  g = zerofloat(n)
  ref = RecursiveExponentialFilter(sigma)
  g = mul(f,f)
  ref.apply(g,g)
  g = div(f,sqrt(g))
  return g

def goModel(wellId):
  setGlobals("all")
  wldata = WellLog.Data.readBinary(csmWellLogs)
  log = wldata.get(wellId)
  model = SynSeis.getModel(log)
  sz = model.sz;
  nz,dz,fz = sz.count,sz.delta,sz.first
  print "nz =",nz," dz =",dz," fz =",fz
  x1,x2,x3 = model.x1,model.x2,model.x3
  print "r: x1 =",x1," x2 =",x2," x3 =",x3
  f,x1,x2,x3 = log.getSamples("v")
  print "v: x1 =",x1[0]," x2 =",x2[0]," x3 =",x3[0]
  f,x1,x2,x3 = log.getSamples("d")
  print "d: x1 =",x1[0]," x2 =",x2[0]," x3 =",x3[0]
  print "v0 =",model.v0
  v = model.v
  d = model.d
  a = model.a
  r = model.r
  t = model.t
  print "r: min =",min(r)," max =",max(r)
  def plot(f):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPoints(sz,f)
    sp.setSize(300,800)
  #plot(model.v)
  #plot(model.d)
  #plot(model.a)
  plot(model.r)
  plot(model.t)

#############################################################################
# Files for well logs and seismic images.
tpDir = "/data/seis/tpd/"
csmWellLogsDir = tpDir+"csm/welllogs/"
csmSeismiczDir = tpDir+"csm/seismicz/"
csmSeismictDir = tpDir+"csm/seismict/"
csmWellLogs = ""
csmSeismicz = csmSeismiczDir+"tpsz.dat"
csmSeismict = csmSeismictDir+"tpst.dat"

# Sampling for 3D seismic images.
sz = Sampling(2762,0.002,0.000)
st = Sampling(1501,0.002,0.000)
s1 = st
s2 = Sampling(357,0.025,0.000)
s3 = Sampling(161,0.025,0.000)

def setGlobals(what):
  global csmWellLogs
  if what=="all":
    csmWellLogs = csmWellLogsDir+"tpwa.dat"
  elif what=="deep":
    csmWellLogs = csmWellLogsDir+"tpwd.dat"
  elif what=="shallow":
    csmWellLogs = csmWellLogsDir+"tpws.dat"
  elif what=="test":
    csmWellLogs = csmWellLogsDir+"tpwt.dat"

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
