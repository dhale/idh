#############################################################################
# Dynamic time warping for 2D images

from imports import *

#############################################################################

#pngDir = "./png"
pngDir = None

s1f,s1g,s2 = None,None,None

def main(args):
  #goSinoImages()
  goSinoWarp()

def goSinoWarp():
  esmooth = 2
  usmooth = 0.0
  shiftMax = 50
  strainMax1 = 0.25
  strainMax2 = 0.05
  fclips = (-2.0,2.0)
  uclips = (-shiftMax,shiftMax)
  f,g = getSinoImages()
  plot(f,"f",fclips,label="Amplitude",png="sinof")
  plot(g,"f",fclips,label="Amplitude",png="sinog")
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setErrorSmoothing(esmooth)
  dw.setShiftSmoothing(usmooth,0.0)
  u = dw.findShifts(f,g)
  h = dw.applyShifts(u,g)
  plot(u,"f",uclips,label="Shift (samples)",png="sinou")
  plot(h,"f",fclips,label="Amplitude",png="sinoh")

def goSinoImages():
  f,g = getSinoImages()
  clips = (-2,2)
  plot(f,"f",clips,label="Amplitude",png="sinof")
  plot(g,"g",clips,label="Amplitude",png="sinog")

def getSinoImages():
  n1f,d1f,f1f = 2001,0.0025,0.0
  n1g,d1g,f1g = 2001,0.0040,0.0
  n2,d2,f2 =  721,0.0150,0.000
  dataDir = "/data/seis/sino/"
  f = readImage(dataDir+"z260.dat",n1f,n2)
  g = readImage(dataDir+"x260.dat",n1g,n2)
  global s1f,s1g,s2
  s1f = Sampling(n1f,d1f,f1f)
  s1g = Sampling(n1g,d1g,f1g)
  s2 = Sampling(n2,d2,f2)
  stretch(1.6,f) # 1.6 = d1g/d1f
  gain(500,f)
  gain(500,g)
  return f,g

#############################################################################
# utilities

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

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def stretch(c,f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0.0,1.0/c,n1)
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

def plot(f,what="f",clips=None,title=None,label=None,png=None):
  n1,n2 = len(f[0]),len(f)
  #width,height,cbwm = 610,815,145
  width,height,cbwm = 700,900,160
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  global s1f,s1g,s2
  if what=="f":
    s1 = s1f
  else:
    s1 = s1g
  pv = sp.addPixels(s1,s2,f)
  #gv = sp.addGrid("H-")
  #gv.setColor(Color.YELLOW)
  if clips:
    pv.setClips(clips[0],clips[1])
  #title = png
  if title:
    sp.setTitle(title)
  if label:
    sp.addColorBar(label)
  sp.setVLabel("Time (s)")
  sp.setHLabel("Distance (km)")
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
  pv0 = panel.addPixels(0,0,s1f,s2,f)
  pv1 = panel.addPixels(0,1,s1g,s2,g)
  pv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  if clips:
    pv0.setClips(clips[0],clips[1])
    pv1.setClips(clips[0],clips[1])
  panel.addColorBar()
  if label:
    panel.addColorBar(label)
  panel.setVLabel(0,"Time (s)")
  panel.setHLabel(0,"Distance (km)")
  panel.setHLabel(1,"Distance (km)")
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

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
