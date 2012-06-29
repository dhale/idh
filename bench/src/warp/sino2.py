#############################################################################
# Dynamic time warping for 2D images

from imports import *

from warp import DynamicWarpingX

#############################################################################

pngDir = "./png"
#pngDir = None

s1f,s1g,s2 = None,None,None

# Different time windows for plotting
ilims = ["0","1","2"]
flims = [(0.0,5.0),(1.0,3.0),(3.0,5.0)]
glims = [(0.0,8.0),(1.6,4.8),(4.8,8.0)]

def main(args):
  goSinoImages()
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
  f,g = getSinoImages()
  u1,h1 = warp1(f,g)
  u2,h2 = warp2(f,h1)
  u = addShifts(u1,u2)
  u = mul(1000.0*s1f.getDelta(),u)
  u1 = mul(1000.0*s1f.getDelta(),u1)
  u2 = mul(1000.0*s1f.getDelta(),u2)
  for i in range(len(ilims)):
    flim = flims[i]
    pre = "si"+ilims[i]
    fpng = pre+"f"
    gpng = pre+"g"
    upng = pre+"u"
    u1png = pre+"u1"
    u2png = pre+"u2"
    h1png = pre+"h1"
    h2png = pre+"h2"
    plot(g ,s1f,fclips,flim,title="X component",cbar=fcbar,png=gpng)
    plot(h1,s1f,fclips,flim,title="X, 1st warping",cbar=fcbar,png=h1png)
    plot(h2,s1f,fclips,flim,title="X, 2nd warping",cbar=fcbar,png=h2png)
    plot(f ,s1f,fclips,flim,title="Z component",cbar=fcbar,png=fpng)
    plot(u1,s1f,None,flim,title="1st warping",cmap=jet,cbar=ucbar,png=u1png)
    plot(u2,s1f,None,flim,title="2nd warping",cmap=jet,cbar=ucbar,png=u1png)
    plot(u ,s1f,None,flim,title="Total warping",cmap=jet,cbar=ucbar,png=upng)

def addShifts(u1,u2):
  dw = DynamicWarpingX(-1,1)
  return add(u2,dw.applyShifts(u2,u1))

def warp2(f,g):
  esmooth = 0
  usmooth = 0.0
  strainMax1 = 0.125
  strainMax2 = 0.125
  shiftMax = 5
  shiftMin = -shiftMax
  dw = DynamicWarpingX(-shiftMax,shiftMax)
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
  usmooth = 1.0
  strainMax1 = 0.125
  shiftMax = 80
  shiftMin = -shiftMax
  dw = DynamicWarpingX(-shiftMax,shiftMax)
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
    pv.setColorModel(ColorMap.PRISM)
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
    cbar = sp.addColorBar(cbar)
    if clips and clips[1]<10:
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
    sp.paintToPng(720,2.0,pngDir+"/"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
