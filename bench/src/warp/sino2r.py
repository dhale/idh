#############################################################################
# Dynamic warping for 2D images

from imports import *
from warp import DynamicWarpingR

#############################################################################

sf,sg,s1,s2 = None,None,None,None

# Different time windows for plotting
ilims = ["0","1","2"]
flims = [(0.0,5.333),(0.8,2.8),(2.8,4.8)] # Vp/Vs = 2 => Tpp/Tps = 3/2
glims = [(0.0,8.000),(1.2,4.2),(4.2,7.2)] # Vp/Vs = 2 => Tpp/Tps = 3/2
#flims = [(0.0,6.000),(0.8,2.8),(2.8,4.8)]
#glims = [(0.0,8.000),(1.2,4.2),(4.2,7.2)]

def main(args):
  goSinoImages()
  #goSinoWarp()

def goSinoImages():
  f,g = getSinoImages()
  clips = (-2.0,2.0)
  for i in [0]: #range(len(ilims)):
    fpng = "si"+ilims[i]+"fpp"
    gpng = "si"+ilims[i]+"gps"
    plot(f,sf,clips,flims[i],title="PP image",cbar="Amplitude",png=fpng)
    plot(g,sg,clips,flims[i],title="PS image",cbar="Amplitude",png=gpng)
  warp1(f,g)

def goSinoWarp():
  fclips = (-2.0,2.0)
  fcbar = "Amplitude"
  ucbar = "Shift (ms)"
  psbar = "Vp/Vs"
  f,g = getSinoImages()
  u1,h1 = warp1(f,g)
  u2,h2 = warp2(f,h1)
  u = addShifts(u1,u2)
  c = sg.delta/sf.delta
  psa = vpvs(u,c,True)
  psi = vpvs(u,c,False)
  u  = mul(1000.0*sf.delta,u)
  u1 = mul(1000.0*sf.delta,u1)
  u2 = mul(1000.0*sf.delta,u2)
  uclips = (225,245)
  #uclips = (215,245)
  for i in [1]:
    flim = flims[i]
    glim = glims[i]
    pre = "si"+ilims[i]
    plot(g ,sf,fclips,flim,title="PS image",cbar=fcbar,png=pre+"g")
    #plot(h1,sf,fclips,flim,title="PS 1st warp",cbar=fcbar,png=pre+"h1")
    plot(h2,sf,fclips,flim,title="PS 2nd warp",cbar=fcbar,png=pre+"h2")
    plot(f ,sf,fclips,flim,title="PP image",cbar=fcbar,png=pre+"f")
    #plot(u1,sf,uclips,flim,title="1st shifts",cmap=jet,cbar=ucbar,png=pre+"u1")
    #plot(u2,sf,uclips,flim,title="2nd shifts",cmap=jet,cbar=ucbar,png=pre+"u2")
    plot(u ,sf,uclips,flim,title="Shifts",cmap=jet,cbar=ucbar,png=pre+"u")
    #plot(psa,sf,(2.0,3.0),flim,title="Vp/Vs (average)",
    #     cmap=jet,cbar=psbar,png=pre+"psa")
    #plot(psi,sf,(1.5,2.5),flim,title="Vp/Vs (interval)",
    #     cmap=jet,cbar=psbar,png=pre+"psi")

def addShifts(u1,u2):
  n1,n2 = len(u1[0]),len(u1)
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniformSampling(n1,1.0,0.0)
  t1 = rampfloat(0.0,1.0,n1)
  s1 = zerofloat(n1)
  y1 = zerofloat(n1)
  us = zerofloat(n1,n2)
  for i2 in range(n2):
    add(u2[i2],t1,s1)
    li.setUniformSamples(u1[i2])
    li.interpolate(n1,s1,y1)
    add(y1,u2[i2],us[i2])
  return us

def vpvs(u,c,avg=False):
  n1,n2 = len(u[0]),len(u)
  if avg:
    ut = div(u,rampfloat(1.0,1.0,0.0,n1,n2))
  else:
    ut = zerofloat(n1,n2)
    rgf = RecursiveGaussianFilter(1.0)
    rgf.apply1X(u,ut)
  ut = add(2.0*c-1.0,mul(2.0*c,ut))
  smoothX(3.0,ut)
  return ut

def smoothX(sigma,x):
  n = 8
  ref = RecursiveExponentialFilter(float(sigma)/sqrt(n))
  for i in range(n):
    ref.apply(x,x)

def warp2(f,g):
  #esmooth,usmooth = 0,0.0
  esmooth,usmooth = 2,8.0
  rsmooth = 101
  strainMax1 = 0.125
  strainMax2 = 0.125
  shiftMax = 10
  shiftMin = -shiftMax
  dw = DynamicWarping(shiftMin,shiftMax)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.REFLECT)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setShiftSmoothing(usmooth)
  e = dw.computeErrors(f,g)
  for ismooth in range(esmooth):
    dw.smoothErrors(e,e)
  d = dw.accumulateForward1(e)
  u = dw.backtrackReverse1(d,e)
  u = dw.smoothShifts(u)
  h = dw.applyShifts(u,g)
  print "warp2: u min =",min(u)," max =",max(u)
  return u,h

def warp1(f,g):
  smin,smax = 0.0,3.0
  dw = DynamicWarpingR(smin,smax,s1,s2)
  ss = dw.getSamplingS()
  ns,n1,n2 = ss.count,s1.count,s2.count
  e = zerofloat(ns,n1)
  for i2 in range(n2/2-11,n2/2+12):
    ei = dw.computeErrors(sf,f[i2],sg,g[i2])
    add(ei,e,e)
  dw.normalizeErrors(e)
  def plotShifts():
    sp = SimplePlot()
    sp.setSize(1400,500)
    pv = sp.addPixels(s1,ss,pow(transpose(e),1.0))
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.GRAY)
    pv.setPercentiles(2,98)
    #pv = sp.addPoints(sf,mul(sf.delta,u1))
  plotShifts()

def samplingPsToPp(vpvs,s1):
  c = 0.5*(1.0+vpvs)
  return Sampling(s1.count,s1.delta/c,s1.first/c)

def getSinoImages():
  global sf,sg,s1,s2
  dataDir = "/data/seis/sino/"
  n1,d1,f1 = 2001,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"z260.dat",n1,n2)
  g = readImage(dataDir+"x260.dat",n1,n2)
  sf = Sampling(n1,d1,f1)
  sg = samplingPsToPp(1.0,sf)
  c = sf.delta/sg.delta
  n1 = int((sf.last-sf.first)/(c*sf.delta))
  s1 = Sampling(n1,sf.delta,sf.first)
  s2 = Sampling(n2,d2,f2)
  gain(100,f)
  gain(100*c,g)
  return f,g

#############################################################################
# utilities

def lowpass(f3db,f):
  """ low-pass filter with specified 3dB frequency, in cycles per sample """
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

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
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET

def plotp(f,s1,clips=None,limits=None,title=None,
         cmap=None,cbar=None,png=None):
  n1,n2 = len(f[0]),len(f)
  width,height,cbwm = 610,815,60
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #sp.addGrid("H-").setColor(Color.BLACK)
  if clips:
    pv.setClips(clips[0],clips[1])
  if limits:
    sp.setVLimits(limits[0],limits[1])
  if cmap:
    pv.setColorModel(cmap)
  if cbar:
    cone = cbar=="Amplitude"
    if cone:
      cbar = "Amplitude (normalized)"
    cbar = sp.addColorBar(cbar)
    if cone:
      cbar.setInterval(1)
  sp.setVLabel("PP time (s)")
  sp.setHLabel("Distance (km)")
  sp.setFontSizeForPrint(8,240)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,3.33333,pngDir+png+".png")

def plots(f,s1,clips=None,limits=None,title=None,
         cmap=None,cbar=None,png=None):
  n1,n2 = len(f[0]),len(f)
  #width,height,cbwm = 610,815,145
  width,height,cbwm = 900,900,180
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #sp.addGrid("H-").setColor(Color.YELLOW)
  if clips:
    pv.setClips(clips[0],clips[1])
  if limits:
    sp.setVLimits(limits[0],limits[1])
  if title:
    sp.setTitle(title)
  if cmap:
    pv.setColorModel(cmap)
  if cbar:
    cone = cbar=="Amplitude"
    cbar = sp.addColorBar(cbar)
    if cone:
      cbar.setInterval(1)
  sp.setVInterval(1.0)
  if s1==sf:
    sp.setVLabel("Z time (s)")
  else:
    sp.setVLabel("X time (s)")
  sp.setHLabel("Distance (km)")
  sp.setFontSizeForPrint(8,150)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,2.0,pngDir+png+".png")

#pngDir = "./png/sino/"
pngDir = None
plot = plotp

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
