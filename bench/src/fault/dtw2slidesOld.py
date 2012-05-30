#############################################################################
# Dynamic time warping for 2D images

import sys
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fault import *
from fault.Util import *
from util import FakeData

#############################################################################

pngDir = "./png"
#pngDir = None

seed = abs(Random().nextInt()/1000000)
seed = 580
seed = 588
print "seed =",seed

nrms = 2.00
npass = 5
stretchMax = 0.20
showLsf = True
smoothShifts = True
smoothSigma = 8

def main(args):
  #goTestImages()
  goTestShifts() #smax = 0.20, nrms = 2.0
  #goFaultImages()
  #goFaultShifts()

def goTestShifts():
  shift = 16
  sigma = shift
  ml = 2*shift
  uclips = (-shift,shift)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax)
  f,g,s = makeTestImages()
  fclips = (min(f),max(f))
  plot(f,fclips,label="Amplitude",png="asynf")
  plot(g,fclips,label="Amplitude",png="asyng")
  e = dw.computeErrors(f,g)
  if npass>=1:
    u = shifts1(dw,e);
    print "errors      u1: sum =",dw.sumErrors(e,u)
    plot(u,uclips,label="Lag (samples)",png="asynu1")
  if npass>=2:
    u = shifts12(dw,e)
    print "errors     u12: sum =",dw.sumErrors(e,u)
    plot(u,uclips,label="Lag (samples)",png="asynu2")
  if npass>=3:
    u = shifts121(dw,e)
    print "errors    u121: sum =",dw.sumErrors(e,u)
    plot(u,uclips,label="Lag (samples)",png="asynu3")
  if npass>=4:
    u = shifts1212(dw,e)
    print "errors   u1212: sum =",dw.sumErrors(e,u)
    plot(u,uclips,label="Lag (samples)",png="asynu4")
  if npass>=5:
    u = shifts12121(dw,e)
    print "errors  u12121: sum =",dw.sumErrors(e,u)
    plot(u,uclips,label="Lag (samples)",png="asynu5")
  if showLsf:
    v = copy(u)
    LocalShiftFinder(ml,sigma).find1(-ml,ml,f,g,v)
    print "errors   u lsf: sum =",dw.sumErrors(e,v)
    plot(v,uclips,label="Lag (samples)",png="asynv")
  if s:
    print "errors       s: sum =",dw.sumErrors(e,s)
    plot(s,uclips,label="Lag (samples)",png="asyns")
  h = align(u,f,g)
  plot(h,fclips,label="Amplitude",png="asynh")

def goFaultShifts():
  global smoothShifts; smoothShifts = True
  global smoothSigma; smoothSigma = 2.0
  shift = 10
  ml = 2*shift
  uclips = (0,8)
  dw = DynamicWarping(0,ml)
  dw.setStretchMax(0.25)
  f,g = makeFaultImages()
  fclips = (-3.0,3.0)
  plot(f,fclips,label="Amplitude",png="faultf")
  plot(g,fclips,label="Amplitude",png="faultg")
  e = dw.computeErrors(f,g)
  u = shifts121(dw,e)
  uclips = (0.0,8*4)
  plot(mul(4,u),uclips,label="Fault throw (ms)",png="faultu")
  #v = copy(u)
  #LocalShiftFinder(ml,shift).find1(0,ml,f,g,v)
  #plot(mul(4,v),uclips,label="Fault throw (ms)",png="faultv")
  h = align(u,f,g)
  plot(h,fclips,label="Amplitude",png="faulth")

def goTestImages():
  f,g,s = makeTestImages()
  plot(f,title="f")
  plot(g,title="g")
  plot(s,title="s")

def goFaultImages():
  f,g = makeFaultImages()
  clips = (-3,3)
  plot(f,clips,label="Amplitude",png="faultf")
  plot(g,clips,label="Amplitude",png="faultg")

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def makeTestImages():
  dip = 30.0
  shift = 16
  n1,n2 = 501,501; f = FakeData.seismic2d2011A(n1,n2,dip)
  #n1,n2 = 462,951; f = readImage("/data/seis/f3d/f3d75.dat",n1,n2)
  f = sub(f,sum(f)/n1/n2)
  #f = slog(f)
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
    rgf = RecursiveGaussianFilter(smoothSigma); rgf.apply00(u,v)
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

""" sampling of 3D image subset used to extract fault images
  j1,j2,j3 = 240, 50,100
  n1,n2,n3 = 222,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.964,0.000,0.000
"""

def plot(f,clips=None,title=None,label=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  sp.plotPanel.setColorBarWidthMinimum(160)
  s12 = len(f)==220
  if s12:
    s1 = Sampling(222,0.004,0.964)
    s2 = Sampling(220,0.025,0.000)
    pv = sp.addPixels(s1,s2,f)
  else:
    pv = sp.addPixels(f)
  if clips:
    pv.setClips(clips[0],clips[1])
  if title:
    sp.setTitle(title)
  if label:
    sp.addColorBar(label)
  if s12:
    sp.setHLabel("Distance along fault strike (km)")
    sp.setVLabel("Time (s)")
    if clips[0]==0.0:
      pv.setColorModel(ColorMap.JET)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(1150,815)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(400,3.2,pngDir+"/"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
