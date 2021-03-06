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
from edu.mines.jtk.util.ArrayMath import *

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
small = 0.01
niter = 80
nbefore = 2
ncycle = 2
nafter = 2
lof = LocalOrientFilter(8)

#############################################################################
# functions

def main(args):
  #doImage();
  #goDiff()
  #goAmpDiff()
  goAmpLaplacianDiffusionFilter()
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
  ai = zerofloat(n1,n2)
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
  x = zerofloat(n1,n2)
  x[(n2-1)/2][(n1-1)/2] = 1.0
  v1 = fillfloat(s,n1,n2)
  v2 = fillfloat(c,n1,n2)
  return x,v1,v2

def makeVectorsRadial(n1,n2):
  k1,k2 = n1/2,n2/2
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
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
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
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

def doImage():
  x = readImage()
  #x = transpose(x)
  #x = makePlaneImage(63.435)
  #x = makePlaneImage(30)
  #x = makeTargetImage()
  #x = flip2(x)
  plot(x,10.0,"x")
  return x

def goAmpDiff():
  #for dip in [20,40,60,80]:
  for dip in [-20,20]:
    suffix = str(dip)
    doAmpDiff(dip,"ahs"+suffix)

def bigger(x):
  m1 = len(x[0])
  m2 = len(x)
  n1 = 2*m1
  n2 = 2*m2
  y = zerofloat(n1,n2)
  copy(m1,m2,0,0,x, 0, 0,y)
  copy(m1,m2,0,0,x,m1, 0,y)
  copy(m1,m2,0,0,x, 0,m2,y)
  copy(m1,m2,0,0,x,m1,m2,y)
  return y

def goDiff():
  x1 = readImage()
  x2 = transpose(x1)
  x3 = makeTargetImage()
  #for x,s in [(x1,"_1"),(x2,"_2"),(x3,"_3")]:
  for x,s in [(x1,"_1")]:
    #x = bigger(bigger(x))
    plot(x,10.0,"x"+s)
    doDiff(x,"d"+s)

def makeRandom(n1,n2):
  r = Random(314159)
  return sub(randfloat(r,n1,n2),0.5)

def doDiff(x,png):
  n1,n2 = len(x[0]),len(x)
  y = zerofloat(n1,n2)
  z = zerofloat(n1,n2)
  s = zerofloat(n1,n2)
  r = makeRandom(n1,n2)
  r = smooth(r)
  #v1,v2 = getV(x)
  #v1,v2 = makeVectorsRadial(n1,n2)
  v1,v2 = makeVectors45(n1,n2)
  #ldf = LocalDiffusionFilter(sigma)
  ldf = LocalDiffusionFilterMp(sigma)
  #ldf = LocalDiffusionFilterCg(sigma,small,niter)
  #ldf = LocalDiffusionFilterMg(sigma,small,niter,nbefore,ncycle,nafter)
  ds = None
  #ds = makeBlock(n1,n2)
  ldf.applyLinearKill(ds,v1,x,y)
  ldf.applyLinearPass(ds,v1,x,z)
  ldf.applyLinearPass(ds,v1,r,s)
  plot(y,10.0,"y"+png)
  plot(z,10.0,"z"+png)
  plot(s, 0.0,"s"+png)

def makeBlock(n1,n2):
  ds = fillfloat(0.0,n1,n2);
  for i2 in range(n2/5,4*n2/5):
    for i1 in range(n1/5,4*n1/5):
  #for i2 in range(n2/2,n2):
  #  for i1 in range(n1):
      ds[i2][i1] = 1.0
  return ds

def goAmpLaplacianDiffusionFilter():
  for dip in [0,20,40,60]:
    for rs in [0,8,25]:
      suffix = str(dip)+str(rs)
      rs /= 100.0
      doAmpLaplacianDiffusionFilter(rs,dip,"ahs"+suffix)

def doAmpLaplacianDiffusionFilter(rs,dip,png=None):
  ldf = LaplacianDiffusionFilter(rs)
  x,v1,v2 = makeImpulse(dip)
  ldt = LocalDiffusionTensors2(1.0,10.0,v1)
  h = copy(x)
  ldf.apply(ldt,x,h)
  ah = frequencyResponse(h)
  plotf(ah,png)

def doAmpDiff(dip,png=None):
  ldf = LocalDiffusionFilterCg(sigma,small,niter)
  #ldf = LocalDiffusionFilterMp(sigma)
  x,v1,v2 = makeImpulse(dip)
  h = copy(x)
  ldf.applyLinearPass(None,v1,x,h)
  ah = frequencyResponse(h)
  plotf(ah,png)

def frequencyResponse(x):
  n1 = len(x[0])
  n2 = len(x)
  n1 = FftComplex.nfftSmall(n1)
  n2 = FftComplex.nfftSmall(n2)
  xr = copy(n1,n2,x)
  xi = zerofloat(n1,n2)
  cx = cmplx(xr,xi)
  fft1 = FftComplex(n1)
  fft2 = FftComplex(n2)
  fft1.complexToComplex1(1,n2,cx,cx)
  fft2.complexToComplex2(1,n1,cx,cx)
  ax = cabs(cx)
  a = zerofloat(n1,n2)
  j1 = n1/2
  j2 = n2/2
  copy(n1-j1,n2-j2,0,0,ax,j1,j2,a)
  copy(j1,j2,n1-j1,n2-j2,ax,0,0,a)
  copy(n1-j1,j2,0,n2-j2,ax,j1,0,a)
  copy(j1,n2-j2,n1-j1,0,ax,0,j2,a)
  return a

def smooth(x):
  n1 = len(x[0])
  n2 = len(x)
  t = zerofloat(n1,n2)
  y = zerofloat(n1,n2)
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
  f = zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def makeTargetImage():
  k = 0.3
  c1 = n1/2
  c2 = n2/2
  f = zerofloat(n1,n2)
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
  return mul(10.0,sin(rampfloat(0.0,c,s,n1,n2)))

def flip2(f):
  n1 = len(f[0])
  n2 = len(f)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    copy(f[n2-1-i2],g[i2])
  return g

def getV(x):
  n1 = len(x[0])
  n2 = len(x)
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof.applyForNormal(x,v2,v1)
  neg(v1,v1);
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
