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

from lcc import *

#############################################################################
# parameters

fontSize = 24
width = 640
height = 505
widthColorBar = 80
dataDir = "/data"
pngDir = None

n1 = 315
n2 = 315
sigma = 8
type = LocalPlaneFilter.Type.HALE1

#############################################################################
# functions

def main(args):
  goPlane()
  #goPef()
  return

def goPlane():
  x = doImage()
  doPlane(x,sigma,LocalPlaneFilter.Type.HALE1)
  #doPlane(x,sigma,LocalPlaneFilter.Type.FOMEL1)
  #doPlane(x,sigma,LocalPlaneFilter.Type.HALE2)
  doPlaneX(x,sigma,LocalPlaneFilter.Type.HALE3)
  #doPlane(x,sigma,LocalPlaneFilter.Type.QUAD)
  doPlane(x,sigma,LocalPlaneFilter.Type.FOMEL1)
  #doPlane(x,sigma,LocalPlaneFilter.Type.FOMEL2)

def goPef():
  x = doImage()
  doPef(x,sigma,type)

def doImage():
  #x = readImage()
  #x = Array.transpose(x)
  #x = makePlaneImage(63.435)
  #x = makePlaneImage(30)
  x = makeTargetImage()
  #x = flip2(x)
  plot(x,10.0,"x")
  return x

def doPlane(x,sigma,type):
  lpf = LocalPlaneFilter(sigma,type)
  p = lpf.find(x);
  #plot(p[0],0.0,None)
  #plot(p[1],0.0,None)
  #plot(p[2],0.0,None)
  y = Array.zerofloat(n1,n2)
  lpf.applyForward(p,x,y);
  plot(y,2.0,"y")
  z = Array.zerofloat(n1,n2)
  lpf.applyInverse(p,y,z);
  plot(z,10.0,"z")
  print "max |z-x| =",Array.max(Array.abs(Array.sub(z,x)))
  #plot(Array.sub(z,x))
  r = Array.sub(Array.randfloat(n1,n2),0.5)
  r = smooth(r)
  #plot(r)
  s = Array.zerofloat(n1,n2)
  lpf.applyInverse(p,r,s);
  plot(s)

def doPlaneX(x,sigma,type):
  lpf = LocalPlaneFilter(sigma,type)
  p = lpf.find(x);
  #plot(p[0],0.0,None)
  #plot(p[1],0.0,None)
  #plot(p[2],0.0,None)
  y = Array.zerofloat(n1,n2)
  #lpf.applyForwardF(p,x,y);
  lpf.applyForwardX(p,x,y);
  plot(y,2.0)
  z = Array.zerofloat(n1,n2)
  lpf.applyInverseX(p,y,z);
  plot(z,10,0)
  print "max |z-x| =",Array.max(Array.abs(Array.sub(z,x)))
  #plot(Array.sub(z,x))
  r = Array.sub(Array.randfloat(n1,n2),0.5)
  r = smooth(r)
  #plot(r)
  s = Array.zerofloat(n1,n2)
  lpf.applyInverseX(p,r,s);
  plot(s)

def doPef(x,sigma,type):
  lpf = LocalPlaneFilter(sigma,type)
  p = lpf.find(x);
  #plot(p[0],0.0,None)
  #plot(p[1],0.0,None)
  #plot(p[2],0.0,None)
  y = Array.zerofloat(n1,n2)
  lpf.applyForward(p,x,y);
  plot(y,1.0,"y")
  z = Array.sub(x,y)
  plot(z,10.0,"z")

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

def doTransposeTest():
  lpf = LocalPlaneFilter(sigma)
  xa = readImage()
  pa = lpf.find(xa);
  xb = Array.transpose(xa)
  pb = lpf.find(xb);
  pb[0] = Array.transpose(pb[0]);
  pb[1] = Array.transpose(pb[1]);
  pb[2] = Array.transpose(pb[2]);
  print "p0 error =",Array.max(Array.sub(pa[0],pb[0]))
  n1 = len(xa[0])
  n2 = len(xa)
  for i2 in range(n2):
    for i1 in range(n1):
      p1 = pb[2][i2][i1]
      p2 = pb[1][i2][i1]
      if p1<0.0:
        p1 = -p1
        p2 = -p2
      pb[1][i2][i1] = p1
      pb[2][i2][i1] = p2
  print "p1 error =",Array.max(Array.sub(pa[1],pb[1]))
  print "p2 error =",Array.max(Array.sub(pa[2],pb[2]))

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def makeTargetImage():
  k = 0.2
  c1 = n1/2
  c2 = n2/2
  f = Array.zerofloat(n1,n2);
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

def panel():
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
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
    frame.paintToPng(200,6,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
