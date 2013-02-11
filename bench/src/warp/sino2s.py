#############################################################################
# Dynamic warping for 2D images
"""
two noisy images
pp ps
ppz1 psz1

smoother shifts
eh001 
eh025
eh050
eh100
shows subsampling of shifts
effect of too little smoothing
effect of too much smoothing
smoothing rough shifts will not work

more accurate shifts
e01z1 e01uz1
e05z1 e05uz1
e11z1 e11uz1
e21z1 e21uz1 (show zoom box for e21uz2 below)
ps image is not shifted version of pp image
  due to noise, differences in reflection coefficients
horizontal averaging of alignment errors helps
  unless shifts vary laterally
  within averaging window

subsampling shifts for better resolution of changes in shift
e21uz2
e21z2 with stencil for h=1
e21z2 with stencil for h=10 (for illustration only)

1D image warping with 1D shifts (mx=121)
pwh01z pwh50z
pmh01z pmh50z

shifts and Vp/Vs
estimates of shifts u
estimates of vp/vs
"""

from imports import *
from warp import DynamicWarpingR,DynamicWarpingS

#############################################################################
nt,nx = 2001,721
ni,nl = 501,351
st,sx = Sampling(nt),Sampling(nx)
si,sl = Sampling(ni),Sampling(nl)
#pngDir = "./png/sinos/"
pngDir = None

def main(args):
  #goImages()
  #goErrors()
  #goShifts1()
  goShifts2()

def vpvs(u):
  n1,n2 = len(u[0]),len(u)
  h = [1,-1]
  v = copy(u)
  for i2 in range(n2):
    Conv.conv(2,0,h,n1,0,u[i2],n1,0,v[i2])
    add(1.0,mul(2.0,v[i2]),v[i2])
    v[i2][0] = v[i2][1]
  print "v: min =",min(v)," max =",max(v)
  return v

def goShifts2():
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  h = copy(f)
  sf = Sampling(ni)
  sg = Sampling(ni+nl)
  dw = DynamicWarpingR(sl.first,sl.last,si,sx)
  dw.setStrainLimits(0.0,2.0,-0.1,0.1)
  #dw.setStrainLimits(0.0,2.0,-0.0,0.0)
  dw.setSmoothness(50,20)
  u = dw.findShifts(sf,f,sg,g)
  v = vpvs(u)
  for ix in range(nx):
    h[ix] = dw.applyShifts(sg,g[ix],u[ix])
  d = sub(h,f)
  print "rms(h-f) =",sqrt(sum(mul(d,d))/ni/nx)
  zoom = True
  plotImage(f,fmax=5,pp=True,zoom=zoom)
  plotImage(h,fmax=5,pp=True,zoom=zoom)
  plotImage(d,fmax=5,pp=True,zoom=zoom)
  if zoom:
    plotImage(v,fmin=1.9,fmax=2.3,pp=True,zoom=zoom)
    plotImage(u,fmin=150,fmax=250,pp=True,zoom=zoom,cv=True)
  else:
    plotImage(v,fmin=2.0,fmax=3.0,pp=True,zoom=zoom)
    plotImage(u,fmin=50,fmax=350,pp=True,zoom=zoom,cv=True)

def goShifts1():
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  h = copy(f)
  e = computeErrors(121,f,g)
  sf = Sampling(ni)
  sg = Sampling(ni+nl)
  ix = nx/2
  dw = DynamicWarpingR(sl.first,sl.last,si)
  dw.setStrainLimits(0.0,2.0)
  plotImage(f,fmax=5,pp=True,zoom=True)
  for hs in [1,50]:
    dw.setSmoothness(hs)
    u = dw.findShifts(e)
    for ix in range(nx):
      h[ix] = dw.applyShifts(sg,g[ix],u)
    d = sub(h,f)
    print "hs =",hs," rms(h-f) =",sqrt(sum(mul(d,d))/ni/nx)
    plotImage(h,fmax=5,pp=True,zoom=True)
    plotImage(d,fmax=5,pp=True,zoom=True)

def goErrors():
  zoom = 0
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  e = computeErrors(5,f,g)
  for ht in [1,25,50,100]:
    nj = 1+(ni-1)/ht
    i = rampint(0,ht,nj)
    u1 = DynamicWarpingS.findShifts(0.0,2.0,e)
    uj = DynamicWarpingS.findShiftsI(0.0,2.0,e,i)
    tj = rampfloat(0.0,ht,nj)
    ti = rampfloat(0.0,1.0,ni)
    ci = CubicInterpolator(CubicInterpolator.Method.MONOTONIC,tj,uj)
    u2 = ci.interpolate(ti)
    zoom = 0
    plotErrors(e,ua=None,ub=u2,tj=tj,uj=uj,zoom=zoom)
  zoom = 1
  ht = 50
  for mx in [1,5,11,21]:
    ixl = nx/2-mx/2
    ixu = ixl+mx
    nj = 1+(ni-1)/ht
    i = rampint(0,ht,nj)
    e = computeErrors(mx,f,g)
    u1 = DynamicWarpingS.findShifts(0.0,2.0,e)
    uj = DynamicWarpingS.findShiftsI(0.0,2.0,e,i)
    tj = rampfloat(0.0,ht,nj)
    ti = rampfloat(0.0,1.0,ni)
    ci = CubicInterpolator(CubicInterpolator.Method.MONOTONIC,tj,uj)
    u2 = ci.interpolate(ti)
    plotErrors(e,zoom=1)
    plotErrors(e,ua=u1,ub=u2,tj=tj,uj=uj,zoom=1)
    if mx==21:
      plotErrors(e,zoom=2)
      plotErrors(e,ua=u1,ub=u2,tj=tj,uj=uj,zoom=2)

def goImages():
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  plotImage(f,pp=True,zoom=False,png="pp")
  plotImage(g,pp=False,zoom=False,png="ps")
  plotImage(f,pp=True,zoom=True,png="ppz")
  plotImage(g,pp=False,zoom=True,png="psz")

#############################################################################

def computeErrors(mx,f,g):
  ixl = nx/2-mx/2
  ixu = ixl+mx
  e = zerofloat(nl,ni)
  for ix in range(ixl,ixu):
    ei = DynamicWarpingS.computeErrors(nl,f[ix],g[ix])
    add(ei,e,e)
  return e

def plotImage(f,fmin=None,fmax=None,pp=True,zoom=False,cv=False,png=None):
  wpt = 252 # width (in points) of plot
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  if fmax:
    if fmin:
      pv.setClips(fmin,fmax)
    else:
      pv.setClips(-fmax,fmax)
  if cv:
    cv = sp.addContours(f)
    cv.setContours(Sampling(36,10,0))
    cv.setLineColor(Color.YELLOW)
  sp.setFontSizeForPrint(8,wpt)
  sp.setHLabel("trace index")
  sp.setVLabel("sample index")
  if zoom:
    if pp:
      sp.setLimits(290,145,410,305)
    else:
      sp.setLimits(290,290,410,560)
  sp.setSize(430,490)
  if pngDir and png:
    sp.paintToPng(720,wpt/72.0,pngDir+png+".png")

def plotErrors(e,ua=None,ub=None,tj=None,uj=None,zoom=0,png=None):
  wpt = 252 # width in points
  nl,ni = len(e[0]),len(e)
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,wpt)
  sp.setHLabel("sample index")
  sp.setVLabel("lag index")
  pv = sp.addPixels(transpose(e))
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
  if tj and uj and len(tj)<ni:
    pv = sp.addPoints(tj,uj)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.WHITE)
  if zoom==0:
    sp.setLimits(0,0,ni-1,nl-1)
  elif zoom==1:
    sp.setLimits(145,145,305,255)
  elif zoom==2:
    #sp.setLimits(235,205,251,217)
    sp.setLimits(227,199,243,211)
  sp.setSize(800,585)
  if pngDir and png:
    sp.paintToPng(720,wpt/72.0,pngDir+png+".png")

def goTest():
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  SimplePlot.asPixels(f)
  SimplePlot.asPixels(g)
  mx = 4
  ht = 50
  #for mx in [1,5,10,20,40,80]:
  #for mx in [1,4,16,64]:
  for ht in [1,25,50,100]:
    ixl = nx/2-mx/2
    ixu = ixl+mx
    nj = 1+(ni-1)/ht
    i = rampint(0,ht,nj)
    e = zerofloat(nl,ni)
    for ix in range(ixl,ixu):
      ei = DynamicWarpingS.computeErrors(nl,f[ix],g[ix])
      add(ei,e,e)
    u1 = DynamicWarpingS.findShifts(0.0,2.0,e)
    uj = DynamicWarpingS.findShiftsI(0.0,2.0,e,i)
    tj = rampfloat(0.0,ht,nj)
    ti = rampfloat(0.0,1.0,ni)
    ci = CubicInterpolator(CubicInterpolator.Method.MONOTONIC,tj,uj)
    u2 = ci.interpolate(ti)
    plotErrors(e)
    if ht>1:
      plotErrors(e,None,u2,tj,uj)
    else:
      plotErrors(e,None,u2)

def goNoisy():
  pass

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
