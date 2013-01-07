#############################################################################
# Dynamic time warping for 2D images
"""
Figures
teaser,nrms=1,npass=5,s=40: diw for noisy large shifts
  f,g,s,u zoomed

nrms=0,npass=1,s=40: dtw for clean large shifts (1c, 2x2)
  f,g,s,u

nrms=0,npass=1: dtw for clean range of shifts (1c, 2x2)
  u01,v01
  u10,v10

nrms=1,npass=1,s=40: dtw for noisy large shifts (1c, 2x2)
  f,g,s,u

nrms=1,s=40: diw for noisy large shift and different passes (2c, 1x4)
  u2,up3,up4,up5
nrms=1,npass=5: diw for noisy range of shifts (2c, 2x4)
  u01,u10,u20,u40
  v01,v10,v20,v40
"""

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

#pngDir = "./png"
pngDir = None

seed = abs(Random().nextInt()/1000000)
seed = 580
seed = 588
print "seed =",seed

nrms = 0.00
npass = 1
shiftMax = 40
stretchMax1 = 1.00
stretchMax2 = 1.00
smoothSigma1 = 2.0/stretchMax1
smoothSigma2 = 2.0/stretchMax2
smoothShifts = True
showLsf = True
teaser = False

s1,s2 = None,None
label1,label2 = None,None
fclips = (-3.0,3.0)

def main(args):
  #goTestImages()
  #goTestShifts() #smax = 0.20, nrms = 2.0
  #goFaultImages()
  #goFaultShifts()
  #goOzImages()
  #goOzTeaser()
  #goOzFiguresFgus()
  #goOzFiguresSmooth()
  #goOzFiguresNpass()
  goOzFiguresDtwUv()

def goOzTeaser():
  global stretchMax1,stretchMax2,smoothShifts,nrms,teaser
  stretchMax1 = 0.25
  stretchMax2 = 1.00
  smoothShifts = True
  smoothSigma1 = 2.0/stretchMax1
  smoothSigma2 = 2.0/stretchMax2
  nrms = 1.00
  teaser = True
  shift = 40
  ml = 3*shift
  uclips = (0,shift*8)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax1,stretchMax2)
  f,g,s = makeOzImages(shift)
  n1,n2 = s1.count,s2.count
  fclips = (-3.0,3.0)
  plot(f,fclips,label="Amplitude",png="toz02f")
  plot(g,fclips,label="Amplitude",png="toz02g")
  e = dw.computeErrors(f,g)
  u = shifts12121(dw,e)
  plot(mul(4,u),uclips,label="Shift (ms)",png="toz02u")
  h = align(u,g)
  v = zerofloat(n1,n2)
  LocalShiftFinder(80,8).find1(-ml,ml,f,g,v)
  plot(mul(4,v),uclips,label="Shift (ms)",png="toz02v")
  plot(mul(4,s),uclips,label="Shift (ms)",png="toz02s")

def goOzFiguresFgus():
  global smoothShifts,smoothSigma1,smoothSigma2,nrms,npass,shiftMax
  stretchMax1 = 0.25
  stretchMax2 = 1.00
  smoothShifts = True
  smoothSigma1 = 2.0/stretchMax1
  smoothSigma2 = 0.0/stretchMax2
  nrms = 0.00
  shift = shiftMax
  npass = 1
  ml = 3*shift
  #uclips = (-shift*4,shift*4)
  uclips = (0,shift*8)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax1,stretchMax2)
  for nrms in [0.0,1.0]:
    f,g,s = makeOzImages(shift)
    n1,n2 = s1.count,s2.count
    fclips = (-3.0,3.0)
    plot2(f,mul(4,s),uclips,label="Shift (ms)",png="oz02fs")
    plot2(g,mul(4,s),uclips,label="Shift (ms)",png="oz02gs")
    e = dw.computeErrors(f,g)
    if npass==1: u = shifts1(dw,e)
    elif npass==2: u = shifts12(dw,e)
    elif npass==3: u = shifts121(dw,e)
    elif npass==4: u = shifts1212(dw,e)
    elif npass==5: u = shifts12121(dw,e)
    h = align(u,g)
    plot2(g,mul(4,u),uclips,label="Shift (ms)",png="oz02gu")
    plot2(h,mul(4,u),uclips,label="Shift (ms)",png="oz02hu")
    plot2(f,mul(4,u),uclips,label="Shift (ms)",png="oz02fu")

def goOzFiguresNpass():
  global smoothShifts,smoothSigma1,smoothSigma2,nrms,npass,shiftMax
  stretchMax1 = 0.25
  stretchMax2 = 1.00
  smoothShifts = True
  smoothSigma1 = 2.0/stretchMax1
  smoothSigma2 = 2.0/stretchMax2
  nrms = 1.00
  shift = 40
  ml = 3*shift
  uclips = (0,shift*8)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax1,stretchMax2)
  f,g,s = makeOzImages(shift)
  e = dw.computeErrors(f,g)
  u0 = shifts1(dw,e)
  u1 = shifts121(dw,e)
  u2 = shifts12121(dw,e)
  h0 = align(u0,g)
  h1 = align(u1,g)
  h2 = align(u2,g)
  plot2(f,mul(4,u0),uclips,label="Shift (ms)",png="oz02fu0")
  plot2(f,mul(4,u1),uclips,label="Shift (ms)",png="oz02fu1")
  plot2(f,mul(4,u2),uclips,label="Shift (ms)",png="oz02fu2")
  plot2(g,mul(4,u1),uclips,label="Shift (ms)",png="oz02gu1")
  plot2(g,mul(4,u0),uclips,label="Shift (ms)",png="oz02gu0")
  plot2(g,mul(4,u2),uclips,label="Shift (ms)",png="oz02gu2")
  plot2(h0,mul(4,u0),uclips,label="Shift (ms)",png="oz02hu0")
  plot2(h1,mul(4,u1),uclips,label="Shift (ms)",png="oz02hu1")
  plot2(h2,mul(4,u2),uclips,label="Shift (ms)",png="oz02hu2")

def goOzFiguresSmooth():
  global smoothShifts,smoothSigma1,smoothSigma2,nrms,npass,shiftMax
  stretchMax1 = 0.25
  stretchMax2 = 1.00
  smoothShifts = True
  smoothSigma1 = 2.0/stretchMax1
  smoothSigma2 = 2.0/stretchMax2
  nrms = 1.00
  shift = 40
  ml = 3*shift
  uclips = (0,shift*8)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax1,stretchMax2)
  f,g,s = makeOzImages(shift)
  d1 = s1.getDelta()
  sl = Sampling(1+2*ml,d1,-ml*d1)
  e0 = dw.computeErrors(f,g)
  e0 = normalize(e0)
  plot3(dw.transposeLag(e0),None,s1,s2,sl,perc=95,png="oze0")
  plot3(dw.transposeLag(e0),s,s1,s2,sl,perc=95,png="oze0s")
  e1 = dw.accumulate1(e0)
  e1 = normalize(e1)
  plot3(dw.transposeLag(e1),None,s1,s2,sl,perc=95,png="oze1")
  plot3(dw.transposeLag(e1),s,s1,s2,sl,perc=95,png="oze1s")
  e2 = dw.accumulate2(e1)
  e2 = normalize(e2)
  plot3(dw.transposeLag(e2),None,s1,s2,sl,perc=95,png="oze2")
  plot3(dw.transposeLag(e2),s,s1,s2,sl,perc=95,png="oze2s")

def goOzFiguresDtwUv():
  global smoothShifts,smoothSigma1,smoothSigma2,nrms,npass,shiftMax
  smoothShifts = True
  for nrms in [0.0]:
    for shiftMax in [40]:
      stretchMax1 = min(1.00,0.25*shiftMax/40.0)
      stretchMax2 = min(1.00,1.00*shiftMax/40.0)
      smoothSigma1 = 1.0/stretchMax1
      if nrms==0.0:
        npass = 1
        smoothSigma2 = 0.0
      else:
        npass = 5
        smoothSigma2 = 1.0/stretchMax2
      shift = shiftMax
      ml = 3*shift
      dw = DynamicWarping(-ml,ml)
      dw.setStretchMax(stretchMax1,stretchMax2)
      uclips = (0,shift*8)
      f,g,s = makeOzImages(shift)
      n1,n2 = s1.count,s2.count
      e = dw.computeErrors(f,g)
      if npass==1: u = shifts1(dw,e)
      elif npass==2: u = shifts12(dw,e)
      elif npass==3: u = shifts121(dw,e)
      elif npass==4: u = shifts1212(dw,e)
      elif npass==5: u = shifts12121(dw,e)
      v = zerofloat(n1,n2)
      lsf = LocalShiftFinder(80,8)
      if nrms==0.0:
        for i2 in range(n2):
          lsf.find1(-ml,ml,f[i2],g[i2],v[i2])
      else:
        lsf.find1(-ml,ml,f,g,v)
      plot2(align(u,g),mul(4,u),uclips,label="Shift (ms)",png="oz02hu")
      plot2(align(v,g),mul(4,v),uclips,label="Shift (ms)",png="oz02hv")

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
  h = align(u,g)
  plot(h,fclips,label="Amplitude",png="asynh")

def goFaultShifts():
  global smoothShifts,smoothSigma1,smoothSigma2; 
  smoothShifts = True
  smoothSigma1 = 2.0
  smoothSigma2 = 2.0
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
  h = align(u,g)
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

def goOzImages():
  f,g,s = makeOzImages(20)
  clips = (-3,3)
  plot(f,clips,label="Amplitude",png="oz02f")
  plot(g,clips,label="Amplitude",png="oz02g")

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def tpow(power,st,f):
  """Applies t^power gain."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  tp = pow(rampfloat(ft,dt,nt),power) # sampled times raised to power
  g = zerofloat(nt,nx)
  for ix in range(nx):
    mul(tp,f[ix],g[ix])
  return g

def gpow(power,f):
  """Applies sgn(f)*|f|^power gain."""
  return mul(sgn(f),pow(abs(f),power))

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

def makeOzImages(shift):
  n1,d1,f1 = 1025,0.004,0.004
  n2,d2,f2 = 127,0.030480,-0.999
  #n2,d2,f2 = 127,0.030480,-1.005840 # correct
  #n2,d2,f2 = 127,0.030480,-2.834640
  fileName = "/data/seis/oz/oz02.F" # suffix F implies floats
  f = readImage(fileName,n1,n2)
  global s1,s2,label1,label2
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  label1 = "Time (s)"
  label2 = "Offset (km)"
  for i2 in range(n2/2):
    fi2 = f[i2]
    f[i2] = f[n2-i2-1]
    f[n2-i2-1] = fi2;
  f = mul(1.0e-12,f)
  f = tpow(2,s1,f)
  f = slog(f)
  w = Warp2.sinusoid(shift,0.0,shift,0.0,n1,n2)
  g = w.warp(f)
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  s = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      s[i2][i1] = w.u1x(i1,i2)
  return f,g,s

""" sampling of 3D image subset used to extract fault images
  j1,j2,j3 = 240, 50,100
  n1,n2,n3 = 222,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.964,0.000,0.000
"""
def makeFaultImages():
  n1,n2 = 222,220
  global s1,s2,label1,label2
  s1 = Sampling(n1,0.004,0.964)
  s2 = Sampling(n2,0.025,0.000)
  label1 = "Time (s)"
  label2 = "Distance along fault strike (km)"
  f = readImage("/data/seis/f3d/faults/s1gfm.dat",n1,n2)
  g = readImage("/data/seis/f3d/faults/s1gfp.dat",n1,n2)
  return f,g

#############################################################################
# shifts

def smooth(u):
  v = copy(u)
  if smoothShifts:
    #RecursiveGaussianFilter(smoothSigma1,smoothSigma2).apply00(u,v)
    RecursiveExponentialFilter(smoothSigma1,smoothSigma2).apply(u,v)
  return v

def normalize(x):
  xmin = min(x)
  xmax = max(x)
  return mul(sub(x,xmin),1.0/(xmax-xmin))

def align(u,g):
  n1,n2 = len(u[0]),len(u)
  si = SincInterp()
  h = copy(g)
  r = rampfloat(0.0,1.0,n1)
  for i2 in range(n2):
    t = add(r,u[i2])
    si.interpolate(n1,1.0,0.0,g[i2],n1,t,h[i2])
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

def plot(f,clips=None,title=None,label=None,png=None):
  n1,n2 = len(f[0]),len(f)
  if n1>510:
    if teaser:
      width,height = 500,415
      cbwm = 70
    else:
      width,height = 610,815
      cbwm = 145
  else:
    width,height = 1150,815
    cbwm = 160
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  global s1,s2
  if s1==None: 
    s1 = Sampling(n1,1.0,0.0)
  if s2==None: 
    s2 = Sampling(n2,1.0,0.0)
  pv = sp.addPixels(s1,s2,f)
  if teaser:
    sp.setVLimits(0.1,2.1)
    if clips[0]==0:
      clips = (160,320)
  #gv = sp.addGrid("H-")
  #gv.setColor(Color.YELLOW)
  if clips:
    pv.setClips(clips[0],clips[1])
  #title = png
  if title:
    sp.setTitle(title)
  if label:
    sp.addColorBar(label)
  sp.setVLabel(label1)
  sp.setHLabel(label2)
  #if clips[0]==0.0:
  #  pv.setColorModel(ColorMap.JET)
  sp.setFontSizeForSlide(0.5,0.8)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    if nrms>0.0:
      png += "n"
    png += str(int(shiftMax))
    sp.paintToPng(720,1.25,pngDir+"/"+png+".png")

def plot2(f,g,gclips=None,label=None,png=None):
  width,height,cbwm = 1095,815,200
  n1,n2 = len(f[0]),len(f)
  global s1,s2
  if s1==None: s1 = Sampling(n1,1.0,0.0)
  if s2==None: s2 = Sampling(n2,1.0,0.0)
  panel = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.mosaic.setWidthTileSpacing(10);
  pv0 = panel.addPixels(0,0,s1,s2,f)
  pv1 = panel.addPixels(0,1,s1,s2,g)
  pv0.setClips(fclips[0],fclips[1])
  if gclips:
    pv1.setClips(gclips[0],gclips[1])
  pv1.setColorModel(ColorMap.JET)
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setVLabel(0,label1)
  panel.setHLabel(0,label2)
  panel.setHLabel(1,label2)
  panel.setColorBarWidthMinimum(cbwm)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    if nrms>0.0:
      png += "n"
    if shiftMax<10:
      png += "0"+str(int(shiftMax))
    else:
      png += str(int(shiftMax))
    frame.paintToPng(720,3.3,pngDir+"/"+png+".png")

def plot3(e,s,s1,s2,sl,perc=None,png=None):
  width,height,cbwm = 860,825,200
  n1,n2,nl = s1.count,s2.count,sl.count
  k1,k2,kl = 300,n2/2,nl/2
  orient = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT;
  axespl = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(orient,axespl,s1,s2,sl,e)
  if perc:
    panel.setPercentiles(0.0,perc)
  panel.mosaic.setWidthElastic(0,100)
  panel.mosaic.setWidthElastic(1,75)
  panel.mosaic.setHeightElastic(0,75)
  panel.mosaic.setHeightElastic(1,100)
  panel.setSlice23(k1)
  panel.setSlice13(k2)
  panel.setSlice12(kl)
  if s:
    st1,sl1 = zerofloat(n1),zerofloat(n1)
    sx2,sl2 = zerofloat(n2),zerofloat(n2)
    for i1 in range(n1): 
      st1[i1] = s1.getValue(i1)
      sl1[i1] = s[k2][i1]*sl.delta
    for i2 in range(n2): 
      sx2[i2] = s2.getValue(i2)
      sl2[i2] = s[i2][k1]*sl.delta
  #panel.setSlice13(70)
  panel.setLabel1("Time (s)")
  panel.setLabel2("Offset (km)")
  panel.setLabel3("Lag (s)")
  panel.setInterval2(2.0)
  panel.setInterval3(0.4)
  panel.setColorModel(ColorMap.JET)
  panel.setLineColor(Color.WHITE)
  if s:
    pv1 = PointsView(st1,sl1)
    pv1.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
    pv1.setLineColor(Color.WHITE)
    pv1.setLineWidth(5)
    pv1.setLineStyle(PointsView.Line.DASH)
    pv2 = PointsView(sx2,sl2)
    pv2.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
    pv2.setLineColor(Color.WHITE)
    pv2.setLineWidth(5)
    pv2.setLineStyle(PointsView.Line.DASH)
    panel.addTiledView(1,1,pv1)
    panel.addTiledView(0,0,pv2)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(1,s1.first,s1.last)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    if nrms>0.0:
      png += "n"
    if shiftMax<10:
      png += "0"+str(int(shiftMax))
    else:
      png += str(int(shiftMax))
    frame.paintToPng(720,3.3,pngDir+"/"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
