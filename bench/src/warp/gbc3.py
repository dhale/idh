#############################################################################
# Dynamic warping for 3D images

from imports import *

from warp import DynamicWarpingX as DynamicWarping

#############################################################################

pngDir = "./png/gbc/"
#pngDir = None

datDir = "/data/seis/gbc/dat/"
n1f,d1f,f1f = 2000,0.00150,0.00150 # p wave, 3 s
n1g,d1g,f1g = 2000,0.00200,0.00200 # s1/s2 (fast/slow) shear wave, 4 s
n2,d2,f2 =  150,0.033531,0.0
n3,d3,f3 =  145,0.033531,0.0
global s1f,s1g,s2,s3
s1f = Sampling(n1f,d1f,f1f)
s1g = Sampling(n1g,d1g,f1g)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)

def main(args):
  #goGbcImages()
  goGbcWarp(True)

def goGbcImages():
  f,g1,g2 = getGbcImages()
  clips = (-1.0,1.0)
  plot3(f,s1f,clips,title="PP", cbar="Amplitude",png="gbcpp")
  plot3(g1,s1g,clips,title="PS1",cbar="Amplitude",png="gbcps1")
  plot3(g2,s1g,clips,title="PS2",cbar="Amplitude",png="gbcps2")

def goGbcWarp(doWarp2):
  fclips = (-1.0,1.0)
  uclips = (-40.0,50.0)
  fcbar = "Amplitude"
  ucbar = "Shift (ms)"
  f,g1,g2 = getGbcImages()
  s12 = {g1:"ps1",g2:"ps2"}
  for g in [g1,g2]:
    ss = s12[g]
    ts = ss.upper()
    u1,h1 = warp1(f,g)
    if doWarp2:
      u2,h2 = warp2(f,h1)
      u = addShifts(u1,u2)
      u = mul(1000.0*s1f.getDelta(),u)
      u2 = mul(1000.0*s1f.getDelta(),u2)
    u1 = mul(1000.0*s1f.getDelta(),u1)
    pre = "gbc"
    fpng = pre+"pp"
    gpng = pre+ss
    upng = pre+ss+"u"
    u1png = pre+ss+"u1"
    u2png = pre+ss+"u2"
    h1png = pre+ss+"w1"
    h2png = pre+ss+"w"
    plot3(g ,s1f,fclips,title=ts,cbar=fcbar,png=gpng)
    plot3(h1,s1f,fclips,title=ts+": 1st warped",cbar=fcbar,png=h1png)
    if doWarp2:
      plot3(h2,s1f,fclips,title=ts+": warped",cbar=fcbar,png=h2png)
    plot3(f ,s1f,fclips,title="PP",cbar=fcbar,png=fpng)
    plot3(u1,s1f,uclips,title=ts+": 1st shifts",cmap=jet,cbar=ucbar,png=u1png)
    if doWarp2:
      plot3(u ,s1f,uclips,title=ts+": shifts",cmap=jet,cbar=ucbar,png=upng)

def addShifts(u1,u2):
  dw = DynamicWarping(-1,1)
  return add(u2,dw.applyShifts(u2,u1))

def warp1(f,g):
  usmooth = 1.0
  #strainMax1 = 0.25 # Vp/Vs = 5/3 + 8/3*0.25 = 1.67 +- 0.67
  strainMax1 = 0.125 # Vp/Vs = 5/3 + 8/3*0.125 = 1.67 +- 0.33
  shiftMax = 60
  shiftMin = -shiftMax
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(strainMax1)
  dw.setShiftSmoothing(usmooth)
  u1 = dw.findShifts1(f,g)
  nl,n1,n2,n3 = 1+shiftMax-shiftMin,len(f[0][0]),len(f[0]),len(f)
  h = zerofloat(n1,n2,n3)
  u = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      copy(u1,u[i3][i2])
      dw.applyShifts(u1,g[i3][i2],h[i3][i2])
  print "warp1: u min =",min(u)," max =",max(u)
  return u,h

def warp2(f,g):
  esmooth = 2
  usmooth = 1.0
  strainMax1 = 0.200
  strainMax2 = 0.200
  strainMax3 = 0.200
  shiftMax = 10
  shiftMin = -shiftMax
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(strainMax1,strainMax2,strainMax3)
  dw.setErrorSmoothing(esmooth)
  dw.setShiftSmoothing(usmooth)
  u = dw.findShifts(f,g)
  h = dw.applyShifts(u,g)
  print "warp2: u min =",min(u)," max =",max(u)
  return u,h

def getGbcImages():
  f = readImage( "pp",n1f,n2,n3)
  g1 = readImage("ps1",n1g,n2,n3)
  g2 = readImage("ps2",n1g,n2,n3)
  #n1f = 1201; f = copy(n1f,n2,f)
  #n1g = 1201; g1 = copy(n1g,n2,g1)
  #n1g = 1201; g2 = copy(n1g,n2,g2)
  stretch(d1g/d1f,f)
  #gain(100,f)
  #gain(100,g1)
  #gain(100,g2)
  return f,g1,g2

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
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  g = zerofloat(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      si.setUniformSamples(f[i3][i2])
      si.interpolate(n1,1.0/c,0.0,g)
      copy(g,f[i3][i2])

def readImage(fileName,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datDir+fileName+".dat")
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(datDir+fileName+".dat")
  aos.writeFloats(x)
  aos.close()
 
#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET

def plot3(f,s1,clips=None,limits=None,title=None,
          cmap=None,cbar=None,png=None):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  width,height,cbwm = 745,900,100
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,f)
  pp.setInterpolation(PixelsView.Interpolation.NEAREST)
  pp.setSlices(633,n2/2,n3/2)
  pp.setColorBarWidthMinimum(cbwm)
  if clips:
    pp.setClips(clips[0],clips[1])
  #if limits:
  #  pp.setVLimits(limits[0],limits[1])
  if title:
    pp.setTitle(title)
  if cmap:
    pp.setColorModel(cmap)
  if cbar:
    cbar = pp.addColorBar(cbar)
  #pp.setVInterval(1.0)
  if s1==s1f:
    pp.setVLabel(1,"PP time (s)")
  else:
    pp.setVLabel(1,"PS time (s)")
  pp.setVLabel(0,"Crossline (km)")
  pp.setHLabel(0,"Inline (km)")
  pp.setHLabel(1,"Crossline (km)")
  pp.setVLimits(1,0.0,2.0)
  mosaic = pp.getMosaic()
  mosaic.setWidthElastic(0,50)
  mosaic.setWidthElastic(1,50)
  mosaic.setHeightElastic(0,50)
  mosaic.setHeightElastic(1,100)
  pf = PlotFrame(pp)
  #pf.setFontSizeForPrint(8,150)
  pf.setSize(width,height)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(720,2.0,pngDir+png+".png")

def plot2(f,s1,clips=None,limits=None,title=None,
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
    cbar = sp.addColorBar(cbar)
    if clips and clips[1]<10:
      cbar.setInterval(1)
  sp.setVInterval(1.0)
  if s1==s1f:
    sp.setVLabel("PP time (s)")
  else:
    sp.setVLabel("PS time (s)")
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
