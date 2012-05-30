#############################################################################
# Dynamic time warping for 3D images

from imports import *
from util import FakeData

#############################################################################

#pngDir = "./png"
pngDir = None

seed = abs(Random().nextInt()/1000000)
seed = 588
print "seed =",seed

n1,n2,n3 = 101,101,101
d1,d2,d3 = 0.004,0.025,0.025
f1,f2,f3 = 0.000,0.000,0.000
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
label1,label2,label3 = "Time (s)","Inline (km)","Crossline (km)"
dataDir = "/data/seis/fake/"

def main(args):
  #goFakeImages()
  goFakeShifts()

def goFakeImages():
  f,g,s = makeFakeImages(6,0.0) # max strain is 6*2*pi/200 ~ 0.188
  print "s: min =",min(s),"max =",max(s)
  writeImage(dataDir+"fakef.dat",f)
  writeImage(dataDir+"fakeg.dat",g)
  writeImage(dataDir+"fakes.dat",s)
  #show(f)
  #show(g)
  #show(s)

def goFakeShifts():
  #f = readImage(dataDir+"fakef.dat",n1,n2,n3)
  #g = readImage(dataDir+"fakeg.dat",n1,n2,n3)
  #s = readImage(dataDir+"fakes.dat",n1,n2,n3)
  f,g,s = makeFakeImages(6,0.5)
  smin,smax = min(s),max(s)
  print "s: min =",smin,"max =",smax
  esmooth = 2
  usmooth = 0.5
  strainMax1 = 0.50
  strainMax2 = 0.25
  strainMax3 = 0.25
  mlag = 4+int(max(-smin,smax))
  dw = DynamicWarping(-mlag,mlag)
  dw.setTempFileDirectory(dataDir);
  dw.setStrainMax(strainMax1,strainMax2,strainMax3)
  dw.setErrorSmoothing(esmooth)
  dw.setShiftSmoothing(usmooth)
  u = dw.findShifts(f,g)
  h = dw.applyShifts(u,g)
  print "s: min =",min(s),"max =",max(s)
  print "u: min =",min(u),"max =",max(u)
  show(s)
  show(u)
  show(g)
  show(f)
  show(h)

def makeFakeImages(shift,nrms):
  f = FakeData.seismic3d2010A(n1,n2,n3,20.0,10.0,30.0,0.5,0.0);
  w = Warp3.sinusoid(shift,0.0,0.0,shift,0.0,0.0,n1,n2,n3)
  #w = Warp3.constant(shift,0.0,0.0,n1,n2,n3)
  g = w.warp1(f)
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  s = w.u1x()
  return f,g,s

#############################################################################
# utilities

def readImage(fileName,n1,n2,n3):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(x)
  aos.close()

def addRickerWavelet(fpeak,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      Conv.conv(nh,-ih,h,n1,0,f[i3][i2],n1,0,g[i3][i2])
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,f,seed=0):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  g = mul(2.0,sub(randfloat(r,n1,n2,n3),0.5))
  g = addRickerWavelet(0.125,g) # bandlimited in time
  rgf = RecursiveGaussianFilter(1.0)
  rgf.applyX0X(g,g)
  rgf.applyXX0(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2/n3)
  grms = sqrt(sum(mul(g,g))/n1/n2/n3)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

def show(f):
  frame = SimpleFrame()
  ip = frame.addImagePanels(f)
  ip.setSlices(n1-1,n2/2,n3/2)
  frame.getOrbitView().setScale(2.0)
  frame.setSize(900,900)

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
