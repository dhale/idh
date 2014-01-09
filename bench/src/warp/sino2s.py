#############################################################################
# Dynamic warping for 2D images
"""
two noisy images
pp ps (no zoom)
ppz1 psz1 (zoomed)

smoother shifts
e005h001
e005h025
e005h050
e005h100
shows subsampling of shifts
effect of too little smoothing
effect of too much smoothing
smoothing rough shifts will not work

more accurate shifts (h = 50)
e001z1 e001uz1
e005z1 e005uz1
e011z1 e011uz1
e721z1 e721uz1 (show zoom box for e721uz2 below)
ps image is not shifted version of pp image
  due to noise, differences in reflection coefficients
horizontal averaging of alignment errors can help
  unless shifts vary laterally within averaging window
  important to see

subsampling shifts for better resolution of changes in shift
e721uz2
e721z2 with stencil for h=1
e721z2 with stencil for h=10 (for illustration only)
shows computational details
accumulator acts as a (non-linear) anti-alias filter

1D and 2D image warping
pw1z pw2z (warped ps)
pd1z pd2z (differences between pp and warped ps)
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
  goShifts1()
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
  dw.setSmoothness(50,20)
  u = dw.findShifts(sf,f,sg,g)
  v = vpvs(u)
  h = dw.applyShifts(sg,g,u)
  d = sub(h,f)
  print "nrms(h,f) =",nrms2(h,f)
  zoom = True
  #plotImage(f,fmax=5,pp=True,zoom=zoom)
  plotImage(h,fmax=5,pp=True,zoom=zoom,png="pw2z")
  plotImage(d,fmax=5,pp=True,zoom=zoom,png="pd2z")
  """
  if zoom:
    plotImage(v,fmin=1.9,fmax=2.3,pp=True,zoom=zoom)
    plotImage(u,fmin=150,fmax=250,pp=True,zoom=zoom,cv=True)
  else:
    plotImage(v,fmin=2.0,fmax=3.0,pp=True,zoom=zoom)
    plotImage(u,fmin=50,fmax=350,pp=True,zoom=zoom,cv=True)
    """
  
def goShifts1():
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  h = copy(f)
  e = computeErrors(nx,f,g)
  sf = Sampling(ni)
  sg = Sampling(ni+nl)
  ix = nx/2
  dw = DynamicWarpingR(sl.first,sl.last,si)
  dw.setStrainLimits(0.0,2.0)
  plotImage(f,fmax=5,pp=True,zoom=True)
  for hs in [50]:
    dw.setSmoothness(hs)
    u = dw.findShifts(e)
    for ix in range(nx):
      h[ix] = dw.applyShifts(sg,g[ix],u)
    d = sub(h,f)
    print "nrms(h,f) =",nrms2(h,f)
    plotImage(h,fmax=5,pp=True,zoom=True,png="pw1z")
    plotImage(d,fmax=5,pp=True,zoom=True,png="pd1z")

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
    hts = str(ht)
    while len(hts)<3: hts = "0"+hts
    plotErrors(e,ua=None,ub=u2,tj=tj,uj=uj,zoom=zoom,png="eh"+hts)
  zoom = 1
  ht = 50
  for mx in [1,5,11,nx]:
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
    mxs = str(mx)
    while len(mxs)<3: mxs = "0"+mxs
    plotErrors(e,zoom=1,png="e"+mxs+"z1")
    plotErrors(e,ua=u1,ub=u2,tj=tj,uj=uj,zoom=1,png="e"+mxs+"uz1")
    if mx==nx:
      plotErrors(e,zoom=2,png="e"+mxs+"z2")
      plotErrors(e,ua=u1,ub=u2,tj=tj,uj=uj,zoom=2,png="e"+mxs+"uz2")

def goImages():
  f,g = getSinoImages()
  f = copy(ni,nx,f)
  g = copy(ni+nl,nx,g)
  plotImage(f,pp=True,zoom=False,png="pp")
  plotImage(g,pp=False,zoom=False,png="ps")
  plotImage(f,pp=True,zoom=True,png="ppz")
  plotImage(g,pp=False,zoom=True,png="psz")

#############################################################################

def rms2(x):
  n1,n2 = len(x[0]),len(x)
  return sqrt(sum(mul(x,x))/n1/n2)

def nrms2(x,y):
  return 2*rms2(sub(x,y))/(rms2(x)+rms2(y))

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
