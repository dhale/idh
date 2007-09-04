import sys
from math import *
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

fontSize = 24
width = 500
height = 520
widthColorBar = 80
dataDir = "/data"
#pngDir = "./png"
pngDir = None

n1 = 315
n2 = 315
small = 0.001
niter = 1000
lof = LocalOrientFilter(16)

#############################################################################
# functions

def main(args):
  #doImage();
  goInterp()
  return

def doImage():
  x = readImage()
  #x = Array.transpose(x)
  #x = makePlaneImage(63.435)
  #x = makePlaneImage(30)
  #x = makeTargetImage()
  #x = flip2(x)
  plot(x,10.0,"x")
  return x

def goInterp():
  x1 = readImage()
  #x2 = Array.transpose(x1)
  #x3 = makeTargetImage()
  #for x,s in [(x1,"_1"),(x2,"_2"),(x3,"_3")]:
  for x,s in [(x1,"_1")]:
    #x = bigger(bigger(x))
    plot(x,10.0,"x"+s)
    doInterp(x,"d"+s)

def doInterp(x,png):
  n1,n2 = len(x[0]),len(x)
  v1,v2 = getV(x)
  #v1,v2 = makeVectorsRadial(n1,n2)
  #v1,v2 = makeVectors45(n1,n2)
  lif = LocalInterpolationFilter(small,niter)
  ds = None
  #ds = makeBlock(n1,n2)
  xf = Array.zerobyte(n1,n2)
  y = Array.zerofloat(n1,n2)
  #for i2 in [1,1*n2/4,2*n2/4,3*n2/4,n2-2]:
  for i2 in [2*n2/4]:
    for i1 in range(n1):
      xf[i2][i1] = 1
      #y[i2][i1] = x[i2][i1]
      y[i2][i1] = i1
  plot(y,n1,"y"+png)
  lif.applyLinear(ds,v1,xf,y)
  plot(y,n1,"y"+png)

def makeVectorsRadial(n1,n2):
  k1,k2 = n1/2,n2/2
  v1 = Array.zerofloat(n1,n2)
  v2 = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    v2i = float(i2-k2)
    for i1 in range(n1):
      v1i = 0.1+float(i1-k1)
      sv = 1.0/sqrt(v1i*v1i+v2i*v2i)
      if v2i<0.0:
        sv = -sv
      v1[i2][i1] = sv*v1i
      v2[i2][i1] = sv*v2i
  return v1,v2

def makeVectors45(n1,n2):
  k1,k2 = n1/2,n2/2
  v1 = Array.zerofloat(n1,n2)
  v2 = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    v2i = 0.1+float(i2-k2)
    v2i = v2i/abs(v2i)
    for i1 in range(n1):
      v1i = 0.1+float(i1-k1)
      v1i = v1i/abs(v1i)
      sv = 1.0/sqrt(v1i*v1i+v2i*v2i)
      if v2i<0.0:
        sv = -sv
      v1[i2][i1] = sv*v1i
      v2[i2][i1] = sv*v2i
  return v1,v2

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

def makePlaneImage(angle):
  a = angle*pi/180.0
  k = 0.3
  c = k*cos(a)
  s = k*sin(a)
  return Array.mul(10.0,Array.sin(Array.rampfloat(0.0,c,s,n1,n2)))

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
  lof.applyForNormal(x,v2,v1)
  Array.neg(v1,v1);
  return v1,v2

#############################################################################
# plot

def plot(f,clip=0.0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  if n1<50 and n2<50:
    s1 = Sampling(n1,1,-(n1-1)/2)
    s2 = Sampling(n2,1,-(n2-1)/2)
  pv = p.addPixels(s1,s2,f)
  if clip!=0.0:
    #pv.setClips(-clip,clip)
    pv.setClips(0,clip)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  frame(p,png)

def plotf(f,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  pv = p.addPixels(f)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(p,png)

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  #p.addColorBar()
  #p.setColorBarWidthMinimum(widthColorBar)
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

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
