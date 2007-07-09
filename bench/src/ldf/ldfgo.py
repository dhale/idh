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
sigma = 14
small = 0.00001
niter = 100
nbefore = 2
ncycle = 2
nafter = 2
lof = LocalOrientFilter(8)

#############################################################################
# functions

def main(args):
  #doImage();
  goDiff()
  #goAmpDiff()
  return

def goIdeal():
  for dip in [20,40,60,80]:
    suffix = str(dip)
    ai = makeIdeal(dip)
    plotf(ai,"ai"+suffix)

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

def makeImpulse(angle):
  n1 = 105
  n2 = 105
  a = angle*pi/180.0
  c = cos(a)
  s = sin(a)
  x = Array.zerofloat(n1,n2)
  x[(n2-1)/2][(n1-1)/2] = 1.0
  u1 = Array.fillfloat(c,n1,n2)
  u2 = Array.fillfloat(-s,n1,n2)
  return x,u1,u2

def doImage():
  x = readImage()
  #x = Array.transpose(x)
  #x = makePlaneImage(63.435)
  #x = makePlaneImage(30)
  #x = makeTargetImage()
  #x = flip2(x)
  plot(x,10.0,"x")
  return x

def goAmpDiff():
  #for dip in [20,40,60,80]:
  for dip in [20]:
    suffix = str(dip)
    doAmpDiff(dip,"ahs"+suffix)

def goDiff():
  x1 = readImage()
  x2 = Array.transpose(x1)
  x3 = makeTargetImage()
  for x,s in [(x1,"_1"),(x2,"_2"),(x3,"_3")]:
  #for x,s in [(x1,"_1")]:
    plot(x,10.0,"x"+s)
    doDiff(x,"d"+s)

def makeRandom():
  r = Random(314159)
  return Array.sub(Array.randfloat(r,n1,n2),0.5)

def doDiff(x,png):
  y = Array.zerofloat(n1,n2)
  z = Array.zerofloat(n1,n2)
  s = Array.zerofloat(n1,n2)
  r = makeRandom()
  r = smooth(r)
  v1,v2 = getV(x)
  #ldf = LocalDiffusionFilter(sigma)
  ldf = LocalDiffusionFilterCg(sigma,small,niter)
  #ldf = LocalDiffusionFilterMg(sigma,small,niter,nbefore,ncycle,nafter)
  ds = makeBlock()
  #ds = None
  ldf.applyInlineKill(ds,v1,x,y)
  ldf.applyInlinePass(ds,v1,x,z)
  ldf.applyInlinePass(ds,v1,r,s)
  plot(y, 2.0,"y"+png)
  plot(z,10.0,"z"+png)
  plot(s, 0.0,"s"+png)

def makeBlock():
  ds = Array.fillfloat(0.0,n1,n2);
  for i2 in range(n2/5,4*n2/5):
    for i1 in range(n1/5,4*n1/5):
      ds[i2][i1] = 1.0
  return ds

def doAmpDiff(dip,png=None):
  x,u1,u2 = makeImpulse(dip)
  h = Array.copy(x)
  su = Array.copy(x)
  sv = Array.copy(x)
  Array.fill(0.0,su)
  Array.fill(1.0,sv)
  sigma = 10
  ldf = LocalDiffusionFilter(sigma)
  #ldf.applyLineSmoothing(u2,x,h)
  ldf.apply(su,sv,u2,x,h)
  h = Array.sub(x,h)
  ah = frequencyResponse(h)
  plotf(ah,png)

def frequencyResponse(x):
  n1 = len(x[0])
  n2 = len(x)
  n1 = FftComplex.nfftSmall(n1)
  n2 = FftComplex.nfftSmall(n2)
  xr = Array.copy(n1,n2,x)
  xi = Array.zerofloat(n1,n2)
  cx = Array.cmplx(xr,xi)
  fft1 = FftComplex(n1)
  fft2 = FftComplex(n2)
  fft1.complexToComplex1(1,n2,cx,cx)
  fft2.complexToComplex2(1,n1,cx,cx)
  ax = Array.cabs(cx)
  a = Array.zerofloat(n1,n2)
  j1 = n1/2
  j2 = n2/2
  Array.copy(n1-j1,n2-j2,0,0,ax,j1,j2,a)
  Array.copy(j1,j2,n1-j1,n2-j2,ax,0,0,a)
  Array.copy(n1-j1,j2,0,n2-j2,ax,j1,0,a)
  Array.copy(j1,n2-j2,n1-j1,0,ax,0,j2,a)
  return a

def smooth(x):
  n1 = len(x[0])
  n2 = len(x)
  t = Array.zerofloat(n1,n2)
  y = Array.zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply0X(x,t)
  rgf.applyX0(t,y)
  """
  bf = ButterworthFilter(0.25,2,ButterworthFilter.Type.LOW_PASS)
  bf.apply1Forward(x,t)
  bf.apply1Reverse(t,y)
  bf.apply2Forward(y,t)
  bf.apply2Reverse(t,y)
  """
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
