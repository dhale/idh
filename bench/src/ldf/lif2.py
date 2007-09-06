import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

from ldf import *

True = 1
False = 0

#############################################################################
# parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

paintBar = None

n1 = 315
n2 = 315
aniso = 10
small = 0.001
niter = 1000
lof = LocalOrientFilter(8)
lof.setGradientSmoothing(1)

#############################################################################
# functions

def main(args):
  global paintBar
  paintBar = makePaintBar()
  #doImage();
  goInterp()
  return

def doImage():
  x = readImage()
  #x = Array.transpose(x)
  #x = makeTargetImage()
  #x = flip2(x)
  plot(x,10.0,gray,"x")
  return x

def goInterp():
  x1 = readImage()
  #x2 = Array.transpose(x1)
  #x3 = makeTargetImage()
  #for x,s in [(x1,"_1"),(x2,"_2"),(x3,"_3")]:
  for x,s in [(x1,"_1")]:
    #x = bigger(bigger(x))
    plot(x,10.0,gray,"x"+s)
    doInterp(x,"d"+s)

def doInterp(x,png):
  n1,n2 = len(x[0]),len(x)
  v1,v2,ds = getV(x)
  #plot2(x,ds)
  #v1,v2 = makeVectorsRadial(n1,n2)
  #v1,v2 = makeVectors45(n1,n2)
  lif = LocalInterpolationFilter(aniso,small,niter)
  ds = None
  es = None
  f = Array.zerobyte(n1,n2)
  y = Array.zerofloat(n1,n2)
  z = Array.zerofloat(n1,n2)
  for i2 in [1,n2/2,n2-2]:
  #for i2 in [1,1*n2/6,2*n2/6,3*n2/6,4*n2/6,5*n2/6,n2-2]:
  #for i2 in [2*n2/4]:
    for i1 in range(n1):
      f[i2][i1] = 1
      y[i2][i1] = i1
      #y[i2][i1] = x[i2][i1]
  #plot(y,0.0,gray,"y"+png)
  #plot(y,10.0,jet,"y"+png)
  z = Array.copy(y)
  lif.applyLinear(ds,es,v1,f,z)
  #plot(z,0.0,gray,"z"+png)
  #plot(z,10.0,jet,"z"+png)
  plot2(x,z,0,315,paintBar)

def bigger(x):
  m1 = len(x[0])
  m2 = len(x)
  n1 = 2*m1
  n2 = 2*m2
  y = Array.zerofloat(n1,n2)
  Array.copy(m1,m2,0,0,x, 0, 0,y)
  Array.copy(m1,m2,0,0,x,m1, 0,y)
  Array.copy(m1,m2,0,0,x, 0,m2,y)
  Array.copy(m1,m2,0,0,x,m1,m2,y)
  return y

def makeBlock(n1,n2):
  ds = Array.fillfloat(0.0,n1,n2);
  for i2 in range(n2/5,4*n2/5):
    for i1 in range(n1/5,4*n1/5):
  #for i2 in range(n2/2,n2):
  #  for i1 in range(n1):
      ds[i2][i1] = 1.0
  return ds

def smooth(x):
  n1 = len(x[0])
  n2 = len(x)
  t = Array.zerofloat(n1,n2)
  y = Array.zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply0X(x,t)
  rgf.applyX0(t,y)
  return y

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def makeTargetImage():
  k = 0.3
  c1 = n1/2
  c2 = n2/2
  f = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    d2 = i2-c2
    for i1 in range(n1):
      d1 = i1-c1
      f[i2][i1] = 10.0*sin(k*sqrt(d1*d1+d2*d2))
  return f

def flip2(f):
  n1 = len(f[0])
  n2 = len(f)
  g = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    Array.copy(f[n2-1-i2],g[i2])
  return g

def getV(x):
  n1 = len(x[0])
  n2 = len(x)
  v1 = Array.zerofloat(n1,n2)
  v2 = Array.zerofloat(n1,n2)
  el = Array.zerofloat(n1,n2)
  lof.apply(x,None,None,None,v1,v2,None,None,el)
  #clips = Clips(0,75,el)
  #elmin = clips.getClipMin()
  #elmax = clips.getClipMax()
  #ds = Array.zerofloat(n1,n2)
  #for i2 in range(n2):
  #  for i1 in range(n1):
  #    ds[i2][i1] = min(elmax,max(elmin,el[i2][i1]))
  ds = el
  return v1,v2,ds
 
#############################################################################
# plot

def plot(f,clip=0.0,cmap=ColorMap.GRAY,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  if clip!=0.0:
    pv.setClips(-clip,clip)
    #pv.setClips(0,clip)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  frame(p,png)

def plot2(f,g,cmin=0,cmax=0,pb=None,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.getGray())
  pv.setClips(-10,10)
  pv = p.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin!=cmax:
    pv.setClips(cmin,cmax)
  pv.setColorModel(ColorMap.getJet(0.3))
  if pb:
    pv.addColorMapListener(pb)
  pb.setRange(cmin+0.25*(cmax-cmin),cmin+0.75*(cmax-cmin))
  frame(p,png)

fontSize = 24
#width = 500
#height = 520
width = 700
height = 600
widthColorBar = 80

def panel():
  #p = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  p.addColorBar()
  p.setColorBarWidthMinimum(widthColorBar)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(100,6,pngDir+"/"+png+".png")
  return frame

def makePaintBar():
  pb = PaintBar()
  frame = JFrame()
  frame.add(pb,BorderLayout.CENTER)
  frame.setSize(80,500)
  frame.setVisible(True)
  return pb

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
