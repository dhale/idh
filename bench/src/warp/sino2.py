#############################################################################
# Dynamic warping for 2D images

from imports import *

from warp import DynamicWarpingX

#############################################################################

#pngDir = "./png/sino/"
pngDir = None

s1f,s1g,s2 = None,None,None

# Different time windows for plotting
ilims = ["0","1","2"]
#flims = [(0.0,5.0),(1.0,3.0),(3.0,5.0)]
#glims = [(0.0,8.0),(1.6,4.8),(4.8,8.0)]
flims = [(0.0,5.0),(0.8,2.8),(2.8,4.8)]
glims = [(0.0,8.0),(1.6,4.8),(4.8,8.0)]

def main(args):
  #goSinoImages()
  goSinoWarp()

def goSinoImages():
  f,g = getSinoImages()
  clips = (-2.0,2.0)
  for i in [0]: #range(len(ilims)):
    fpng = "si"+ilims[i]+"fz"
    gpng = "si"+ilims[i]+"gx"
    plot(f,s1f,clips,flims[i],title="Z component",cbar="Amplitude",png=fpng)
    plot(g,s1g,clips,glims[i],title="X component",cbar="Amplitude",png=gpng)

def goSinoWarp():
  fclips = (-2.0,2.0)
  fcbar = "Amplitude"
  ucbar = "Shift (ms)"
  psbar = "Vp/Vs"
  f,g = getSinoImages()
  u1,h1 = warp1(f,g)
  u2,h2 = warp2(f,h1)
  u = addShifts(u1,u2)
  c = s1g.delta/s1f.delta
  psa = vpvs(u,c,True)
  psi = vpvs(u,c,False)
  """
  SimplePlot.asPoints(u1[350])
  SimplePlot.asPoints(u2[350])
  SimplePlot.asPoints(u[350])
  SimplePlot.asPoints(clip(1.5,2.9,psi[350]))
  """
  u  = mul(1000.0*s1f.delta,u)
  u1 = mul(1000.0*s1f.delta,u1)
  u2 = mul(1000.0*s1f.delta,u2)
  #for i in range(len(ilims)):
  for i in [0]:
    flim = flims[i]
    pre = "si"+ilims[i]
    plot(g ,s1f,fclips,flim,title="X component",cbar=fcbar,png=pre+"g")
    plot(h1,s1f,fclips,flim,title="X 1st warp",cbar=fcbar,png=pre+"h1")
    plot(h2,s1f,fclips,flim,title="X 2nd warp",cbar=fcbar,png=pre+"h2")
    plot(f ,s1f,fclips,flim,title="Z component",cbar=fcbar,png="f")
    plot(u1,s1f,None,flim,title="1st shifts",cmap=jet,cbar=ucbar,png=pre+"u1")
    plot(u2,s1f,None,flim,title="2nd shifts",cmap=jet,cbar=ucbar,png=pre+"u2")
    plot(u ,s1f,None,flim,title="Shifts",cmap=jet,cbar=ucbar,png=pre+"u")
    plot(psa,s1f,(2.0,3.5),flim,title="Vp/Vs (average)",
         cmap=jet,cbar=psbar,png=pre+"psa")
    plot(psi,s1f,(1.5,3.0),flim,title="Vp/Vs (interval)",
         cmap=jet,cbar=psbar,png=pre+"psi")

def addShifts(u1,u2):
  dw = DynamicWarpingX(-1,1)
  return add(u2,dw.applyShifts(u2,u1))

def vpvs(u,c,avg=False):
  n1,n2 = len(u[0]),len(u)
  if avg:
    ut = div(u,rampfloat(1.0,1.0,0.0,n1,n2))
  else:
    ut = zerofloat(n1,n2)
    rgf = RecursiveGaussianFilter(1.0)
    rgf.apply1X(u,ut)
  ut = add(2.0*c-1.0,mul(2.0*c,ut))
  RecursiveGaussianFilter(2.0).apply00(ut,ut)
  return ut

def warp2(f,g):
  #esmooth,usmooth = 0,0.0
  #esmooth,usmooth = 0,1.0
  esmooth,usmooth = 2,1.0
  strainMax1 = 0.125
  strainMax2 = 0.125
  shiftMax = 5
  shiftMin = -shiftMax
  dw = DynamicWarpingX(-shiftMax,shiftMax)
  dw.setErrorExtrapolation(DynamicWarpingX.ErrorExtrapolation.AVERAGE)
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
  usmooth = 4.0
  strainMax1 = 0.125
  shiftMax = 80
  shiftMin = -shiftMax
  dw = DynamicWarpingX(-shiftMax,shiftMax)
  dw.setErrorExtrapolation(DynamicWarpingX.ErrorExtrapolation.AVERAGE)
  dw.setStrainMax(strainMax1)
  dw.setShiftSmoothing(usmooth)
  e = dw.computeErrors(f,g)
  nl,n1,n2 = len(e[0][0]),len(e[0]),len(e)
  e1 = zerofloat(nl,n1)
  for i2 in range(n2):
    add(e[i2],e1,e1)
  dw.normalizeErrors(e1)
  d1 = dw.accumulateForward(e1)
  u1 = dw.backtrackReverse(d1,e1)
  u1 = dw.smoothShifts(u1)
  if False:
    sp = SimplePlot()
    sp.setSize(1800,500)
    sl = Sampling(nl,1.0,shiftMin)
    s1 = Sampling(n1,1.0,0.0)
    pv = sp.addPixels(s1,sl,pow(transpose(e1),0.25))
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.JET)
    pv.setPercentiles(2,98)
    pv = sp.addPoints(u1)
  h = zerofloat(n1,n2)
  u = zerofloat(n1,n2)
  for i2 in range(n2):
    copy(u1,u[i2])
    dw.applyShifts(u1,g[i2],h[i2])
  print "warp1: u min =",min(u)," max =",max(u)
  return u,h

def getSinoImages():
  dataDir = "/data/seis/sino/"
  n1f,d1f,f1f = 2001,0.0025,0.0 # z component, 0 to 5 s
  n1g,d1g,f1g = 2001,0.0040,0.0 # x component, 0 to 8 s
  n2,d2,f2 =  721,0.0150,0.000
  global s1f,s1g,s2
  s1f = Sampling(n1f,d1f,f1f)
  s1g = Sampling(n1g,d1g,f1g)
  s2 = Sampling(n2,d2,f2)
  f = readImage(dataDir+"z260.dat",n1f,n2)
  g = readImage(dataDir+"x260.dat",n1g,n2)
  #n1f = 1201; f = copy(n1f,n2,f)
  #n1g = 1201; g = copy(n1g,n2,g)
  stretch(d1g/d1f,f)
  gain(100,f)
  gain(100,g)
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

def stretch(c,f):
  """ stretch (supersample) by specified factor c time sampling of image f """
  n1,n2 = len(f[0]),len(f)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  g = zerofloat(n1)
  for i2 in range(n2):
    si.setUniformSamples(f[i2])
    si.interpolate(n1,1.0/c,0.0,g)
    copy(g,f[i2])

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

def plot(f,s1,clips=None,limits=None,title=None,
         cmap=None,cbar=None,png=None):
  n1,n2 = len(f[0]),len(f)
  #width,height,cbwm = 610,815,145
  width,height,cbwm = 900,900,150
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
  if s1==s1f:
    sp.setVLabel("Z time (s)")
  else:
    sp.setVLabel("X time (s)")
  sp.setHLabel("Distance (km)")
  sp.setFontSizeForPrint(8,150)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,2.0,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
