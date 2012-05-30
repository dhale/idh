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

from imports import *
from util import FakeData

#############################################################################

pngDir = "./png"
#pngDir = None

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

def main(args):
  #goTestImages()
  #goTestShifts() #smax = 0.20, nrms = 2.0
  #goFaultImages()
  #goFaultShifts()
  #goOzImages()
  #goOzTeaser()
  #goOzFiguresFgus()
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
    plot2(f,g,fclips,label="Amplitude",png="oz02fg")
    e = dw.computeErrors(f,g)
    if npass==1: u = shifts1(dw,e)
    elif npass==2: u = shifts12(dw,e)
    elif npass==3: u = shifts121(dw,e)
    elif npass==4: u = shifts1212(dw,e)
    elif npass==5: u = shifts12121(dw,e)
    plot2(mul(4,u),mul(4,s),uclips,label="Shift (ms)",png="oz02us")

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
  u1 = shifts1(dw,e)
  u2 = shifts12(dw,e)
  u3 = shifts121(dw,e)
  u4 = shifts1212(dw,e)
  u5 = shifts12121(dw,e)
  plot2(mul(4,u2),mul(4,u3),uclips,label="Shift (ms)",png="oz02p23")
  plot2(mul(4,u4),mul(4,u5),uclips,label="Shift (ms)",png="oz02p45")
  print "u1-u5: min/max =",min(sub(u1,u5)),max(sub(u1,u5))
  print "u2-u5: min/max =",min(sub(u2,u5)),max(sub(u2,u5))
  print "u3-u5: min/max =",min(sub(u3,u5)),max(sub(u3,u5))
  print "u4-u5: min/max =",min(sub(u4,u5)),max(sub(u4,u5))

def goOzFiguresDtwUv():
  global smoothShifts,smoothSigma1,smoothSigma2,nrms,npass,shiftMax
  smoothShifts = True
  for nrms in [1.0]:
    for shiftMax in [2,20,40]:
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
      plot2(mul(4,u),mul(4,v),uclips,label="Shift (ms)",png="oz02uv")

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
  n2,d2,f2 = 127,0.030480,-1.005840
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
  sp.setFontSizeForPrint(8,120)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    if nrms>0.0:
      png += "n"
    png += str(int(shiftMax))
    sp.paintToPng(720,1.25,pngDir+"/"+png+".png")

def plot2(f,g,clips=None,label=None,png=None):
  width,height,cbwm = 1095,815,145
  n1,n2 = len(f[0]),len(f)
  global s1,s2
  if s1==None: s1 = Sampling(n1,1.0,0.0)
  if s2==None: s2 = Sampling(n2,1.0,0.0)
  panel = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.mosaic.setWidthTileSpacing(10);
  pv0 = panel.addPixels(0,0,s1,s2,f)
  pv1 = panel.addPixels(0,1,s1,s2,g)
  if clips:
    pv0.setClips(clips[0],clips[1])
    pv1.setClips(clips[0],clips[1])
  panel.addColorBar()
  if label:
    panel.addColorBar(label)
  panel.setVLabel(0,label1)
  panel.setHLabel(0,label2)
  panel.setHLabel(1,label2)
  panel.setColorBarWidthMinimum(cbwm)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSizeForPrint(8,240)
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

def plot4(u,clips=None,label=None,png=None):
  width,height,cbwm = 1095,815,145
  n1,n2 = len(u[0][0]),len(u[0])
  global s1,s2
  if s1==None: s1 = Sampling(n1,1.0,0.0)
  if s2==None: s2 = Sampling(n2,1.0,0.0)
  panel = PlotPanel(1,4,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.mosaic.setWidthTileSpacing(10);
  panel.setVLabel(0,label1)
  for i in range(len(u)):
    pv = panel.addPixels(0,i,s1,s2,u[i])
    if clips:
      pv.setClips(clips[0],clips[1])
  panel.setHLabel(i,label2)
  panel.addColorBar()
  if label:
    panel.addColorBar(label)
  panel.setColorBarWidthMinimum(cbwm)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSizeForPrint(8,504)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    if nrms>0.0:
      png += "n"
    frame.paintToPng(720,3.3,pngDir+"/"+png+".png")


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
