#############################################################################
# Dynamic time warping for 2D images

from imports import *
from util import FakeData

#############################################################################

#pngDir = "./png"
pngDir = None

seed = abs(Random().nextInt()/1000000)
seed = 580
seed = 588
print "seed =",seed

nrms = 2.00
npass = 3
stretchMax = 0.25
showLsf = False
smoothShifts = True

def main(args):
  #goTestImages()
  #goFaultImages()
  #goTestShifts() #smax = 0.20, nrms = 2.0
  goFaultShifts() #smax = 0.25, npass = 3

def goTestShifts():
  shift = 16
  sigma = shift
  ml = 2*shift
  uclips = (-shift,shift)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax)
  f,g,s = makeTestImages()
  fclips = (min(f),max(f))
  plot(f,"f",fclips)
  plot(g,"g",fclips)
  e = dw.computeErrors(f,g)
  if npass>=1:
    u = shifts1(dw,e);
    print "errors      u1: sum =",dw.sumErrors(e,u)
    plot(u,"u1",uclips)
  if npass>=2:
    u = shifts12(dw,e)
    print "errors     u12: sum =",dw.sumErrors(e,u)
    plot(u,"u12",uclips)
  if npass>=3:
    u = shifts121(dw,e)
    print "errors    u121: sum =",dw.sumErrors(e,u)
    plot(u,"u121",uclips)
  if npass>=4:
    u = shifts1212(dw,e)
    print "errors   u1212: sum =",dw.sumErrors(e,u)
    plot(u,"u1212",uclips)
  if npass>=5:
    u = shifts12121(dw,e)
    print "errors  u12121: sum =",dw.sumErrors(e,u)
    plot(u,"u12121",uclips)
  if showLsf:
    v = copy(u)
    LocalShiftFinder(ml,sigma).find1(-ml,ml,f,g,v)
    print "errors   u lsf: sum =",dw.sumErrors(e,v)
    plot(v,"u lsf",uclips)
  if s:
    print "errors       s: sum =",dw.sumErrors(e,s)
    plot(s,"s")
  h = align(u,f,g)
  plot(h,"h",fclips)

def goFaultShifts():
  shift = 10
  ml = 2*shift
  uclips = (0,8)
  dw = DynamicWarping(0,ml)
  dw.setStretchMax(stretchMax)
  f,g = makeFaultImages()
  fclips = (min(f),max(f))
  plot(f,"f",fclips)
  plot(g,"g",fclips)
  e = dw.computeErrors(f,g)
  if npass>=1:
    u = shifts1(dw,e);
    print "errors      u1: sum =",dw.sumErrors(e,u)
    plot(u,"u1",uclips)
  if npass>=2:
    u = shifts12(dw,e)
    print "errors     u12: sum =",dw.sumErrors(e,u)
    plot(u,"u12",uclips)
  if npass>=3:
    u = shifts121(dw,e)
    print "errors    u121: sum =",dw.sumErrors(e,u)
    plot(u,"u121",uclips)
  if npass>=4:
    u = shifts1212(dw,e)
    print "errors   u1212: sum =",dw.sumErrors(e,u)
    plot(u,"u1212",uclips)
  if npass>=5:
    u = shifts12121(dw,e)
    print "errors  u12121: sum =",dw.sumErrors(e,u)
    plot(u,"u12121",uclips)
  if showLsf:
    v = copy(u)
    LocalShiftFinder(ml,shift).find1(0,ml,f,g,v)
    print "errors   u lsf: sum =",dw.sumErrors(e,v)
    plot(v,"u lsf",uclips)
  h = align(u,f,g)
  plot(h,"h",fclips)

def goTestImages():
  f,g,s = makeTestImages()
  plot(f,"f")
  plot(g,"g")
  plot(s,"s")

def goFaultImages():
  f,g = makeFaultImages()
  plot(f,"f")
  plot(g,"g")

def makeTestImages():
  dip = 30.0
  shift = 16
  n1,n2 = 501,501; f = FakeData.seismic2d2011A(n1,n2,dip)
  #n1,n2 = 462,951; f = readImage("/data/seis/f3d/f3d75.dat",n1,n2)
  f = sub(f,sum(f)/n1/n2)
  #w = Warp2.constant(shift,0.0,n1,n2)
  w = Warp2.sinusoid(shift,0.0,n1,n2)
  g = w.warp(f)
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  s = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      s[i2][i1] = w.u1x(i1,i2)
  return f,g,s

def makeFaultImages():
  n1,n2 = 222,220; f = readImage("/data/seis/f3d/faults/s1gfm.dat",n1,n2)
  n1,n2 = 222,220; g = readImage("/data/seis/f3d/faults/s1gfp.dat",n1,n2)
  return f,g

#############################################################################
# shifts

def smooth(u):
  v = copy(u)
  if smoothShifts:
    rgf = RecursiveGaussianFilter(8); rgf.apply00(u,v)
  return v

def normalize(x):
  xmin = min(x)
  xmax = max(x)
  return mul(sub(x,xmin),1.0/(xmax-xmin))

def align(u,f,g):
  n1,n2 = len(u[0]),len(u)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  h = copy(g)
  r = rampfloat(0.0,1.0,n1)
  for i2 in range(n2):
    t = add(r,u[i2])
    si.setUniformSamples(g[i2])
    si.interpolate(n1,t,h[i2])
  return h

def shifts1(dw,e):
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

def shifts12(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  d = dw.accumulateForward2(e)
  u = dw.findShiftsReverse2(d,e)
  return smooth(u)

def shifts121(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

def shifts1212(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  e = dw.accumulate1(e)
  e = normalize(e)
  d = dw.accumulateForward2(e)
  u = dw.findShiftsReverse2(d,e)
  return smooth(u)

def shifts12121(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

#############################################################################
# utilities

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def addRickerWavelet(fpeak,f):
  n1,n2 = len(f[0]),len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    Conv.conv(nh,-ih,h,n1,0,f[i2],n1,0,g[i2])
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,f,seed=0):
  n1,n2 = len(f[0]),len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  g = mul(2.0,sub(randfloat(r,n1,n2),0.5))
  g = addRickerWavelet(0.125,g) # same wavelet used in signal
  rgf = RecursiveGaussianFilter(1.0)
  rgf.applyX0(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2)
  #frms = max(abs(f))
  grms = sqrt(sum(mul(g,g))/n1/n2)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

def plot(f,title=None,clips=None):
  #sp = SimplePlot.asPixels(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(f)
  if clips:
    pv.setClips(clips[0],clips[1])
  if title:
    sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(900,900)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
