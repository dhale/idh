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
sigma = 32
small = 0.1
niter = 1000

#############################################################################
# functions

def main(args):
  goSmooth()
  return

def goSmooth():
  x1 = readImage()
  for x,s in [(x1,"_1")]:
    doSmooth(x,"d"+s)

def doSmooth(x,png):
  n1,n2 = len(x[0]),len(x)
  r = makeRandom(n1,n2)
  s = Array.zerofloat(n1,n2)
  y = Array.zerofloat(n1,n2)
  lof = LocalOrientFilter(8)
  t = lof.applyForTensors(x)
  eu = Array.fillfloat(0.0,n1,n2)
  ev = Array.fillfloat(1.0,n1,n2)
  #ev = makeBlock(n1,n2)
  #ev = makeNotch(n1,n2)
  #ev = makeCircle(n1,n2)
  t.setEigenvalues(eu,ev);
  scale = 0.5*sigma*sigma;
  lsf = LocalSmoothingFilterX(scale,small,niter)
  Array.zero(x)
  x[(n2-1)/2][(n1-1)/2] = 1.0
  lsf.apply(t,x,y);
  print "y min =",Array.min(y)," max =",Array.max(y)
  lsf.apply(t,r,s);
  plot(s, 0.0,"s"+png)
  plot(x, 0.0,"x"+png)
  plot(y, 0.0,"y"+png)

def makeIdeal(angle):
  n1 = 105
  n2 = 105
  j1 = (n1-1)/2
  j2 = (n2-1)/2
  s1 = 1.0/(n1-1)
  s2 = 1.0/(n2-1)
  a = angle*pi/180.0
  v1 = sin(a)
  v2 = cos(a)
  ai = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    k2 = (i2-j2)*s2
    for i1 in range(n1):
      k1 = (i1-j1)*s1
      if k1==0.0 and k2==0.0:
        ai[i2][i1] = 0.0
      else:
        s = v1*k1+v2*k2
        t = 0.001+k1*k1+k2*k2
        ai[i2][i1] = 1.0-0.05/(0.05+s*s/t)
  return ai

def makeRandom(n1,n2):
  r = Random(314159)
  return smooth(Array.sub(Array.randfloat(r,n1,n2),0.5))

def makeBlock(n1,n2):
  ds = Array.fillfloat(0.0,n1,n2);
  for i2 in range(n2):
    for i1 in range(n1):
      #if i1>i2: ds[i2][i1] = 1.0
      if n2/5<=i2<=4*n2/5 and n1/5<=i1<=4*n1/5: ds[i2][i1] = 1.0
  return ds

def makeCircle(n1,n2):
  s = Array.fillfloat(0.0,n1,n2);
  k1 = (n1-1)/2
  k2 = (n2-1)/2
  ds = (0.4*min(n1,n2))**2
  for i2 in range(n2):
    d2 = (i2-k2)**2
    for i1 in range(n1):
      d1 = (i1-k1)**2
      if d1+d2>ds: s[i2][i1] = 1.0
  return s

def makeNotch(n1,n2):
  ds = Array.fillfloat(0.0,n1,n2);
  m2 = (n2-1)/2
  for i2 in range(n2):
    for i1 in range(n1):
      if i2<m2 or i2>m2: ds[i2][i1] = 1.0
  return ds

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

def smooth(x):
  n1,n2 = len(x[0]),len(x)
  y = Array.zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply00(x,y)
  return y

def flip2(f):
  n1,n2 = len(x[0]),len(x)
  g = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    Array.copy(f[n2-1-i2],g[i2])
  return g

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
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
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
