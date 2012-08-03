#############################################################################
# Dynamic warping for 2D images

from imports import *
from util import FakeData

#############################################################################

#pngDir = "./png"
pngDir = None

seed = abs(Random().nextInt()/1000000)
seed = 580
seed = 588
print "seed =",seed

nrms = 0.00
usmooth = 0.0
esmooth = 0
shiftMax = 40
strainMax1 = 1.00
strainMax2 = 1.00

s1,s2 = None,None
label1,label2 = None,None

def main(args):
  goOzImages()
  #goOzShifts()

def goOzShifts():
  global nrms,esmooth,usmooth,shiftMax
  nrms = 0.0
  esmooth = 2
  usmooth = 1.0
  strainMax1 = 0.25
  strainMax2 = 1.00
  shift = shiftMax
  mlag = 3*shift
  uclips = (0,shift*8)
  dw = DynamicWarping(-mlag,mlag)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setErrorSmoothing(esmooth)
  for nrms in [0.0,1.0]:
    f,g,s = makeOzImages(shift)
    n1,n2 = s1.count,s2.count
    fclips = (-3.0,3.0)
    plot2(f,g,fclips,label="Amplitude",png="oz02fg")
    if esmooth==0:
      dw.setShiftSmoothing(usmooth,0.0)
    else:
      dw.setShiftSmoothing(usmooth)
    u = dw.findShifts(f,g)
    h = dw.applyShifts(u,g)
    plot2(h,g,fclips,label="Amplitude",png="oz02hg")
    plot2(mul(4,u),mul(4,s),uclips,label="Shift (ms)",png="oz02us")

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
  for i2 in range(n2/2): # flip horizontally
    fi2 = f[i2]
    f[i2] = f[n2-i2-1]
    f[n2-i2-1] = fi2;
  f = mul(1.0e-12,f)
  f = tpow(2,s1,f)
  f = slog(f)
  w = WarpFunction2.constantPlusSinusoid(shift,0.0,shift,0.0,n1,n2)
  g = w.warp(f)
  s = w.u1x()
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  return f,g,s

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
  width,height,cbwm = 610,815,145
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  global s1,s2
  if s1==None: 
    s1 = Sampling(n1,1.0,0.0)
  if s2==None: 
    s2 = Sampling(n2,1.0,0.0)
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
  pv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv1.setInterpolation(PixelsView.Interpolation.NEAREST)
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
