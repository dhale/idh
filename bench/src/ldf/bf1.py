# Test 1D bilateral filter.
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
from edu.mines.jtk.util.ArrayMath import *

from ldf import *

#############################################################################
# functions

def main(args):
  goRandomBlocks()
  #goShowRangeFunctions()
  return

def rms(x):
  n = len(x)
  m = sum(x)/n
  x = sub(x,m)
  return sqrt(sum(mul(x,x))/n)

def mad(x):
  m = Quantiler.estimate(0.5,x)
  y = abs(sub(x,m))
  return Quantiler.estimate(0.5,y)

def qqd(x):
  return 0.5*(Quantiler.estimate(0.75,x)-Quantiler.estimate(0.25,x))

def goRandomBlocks():
  n1 = 801
  x = makeBlocks(n1)
  y = add(x,makeRandom(3141,8.0,3.0,n1))
  yrms,ymad,yqqd = rms(y),mad(y),qqd(y)
  print "y rms =",yrms," mad =",ymad," qqd =",yqqd
  plot2(x,y,-9,9)
  sigmaS = 20.0
  z = zerofloat(n1)
  for sigmaX in [1.0*yqqd,10.0*yqqd]:
    bf = BilateralFilter(sigmaS,sigmaX)
    bf.apply(y,z)
    plot2(x,z,-9,9)

def plot2(x,y,ymin=0.0,ymax=0.0):
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,240)
  sp.setHLabel("Sample index")
  sp.setVLabel("Amplitude")
  sp.setSize(720,300)
  pv = sp.addPoints(x) 
  pv.setLineWidth(3)
  pv.setLineStyle(PointsView.Line.DOT)
  pv = sp.addPoints(y) 
  pv.setLineWidth(3)
  #pv.setLineColor(Color.RED)
  if ymin<ymax:
    sp.setVLimits(ymin,ymax)

gauss = BilateralFilter.Type.GAUSS
huber = BilateralFilter.Type.HUBER
tukey = BilateralFilter.Type.TUKEY

def goShowRangeFunctions():
  xmin,xmax,sigma = -6.0,6.0,2.0
  nx = 101
  dx = (xmax-xmin)/(nx-1)
  fx = xmin
  sx = Sampling(nx,dx,fx)
  yg = BilateralFilter.sampleRangeFunction(gauss,sigma,sx)
  yh = BilateralFilter.sampleRangeFunction(huber,sigma,sx)
  yt = BilateralFilter.sampleRangeFunction(tukey,sigma,sx)
  #x = rampfloat(fx,dx,nx)
  #yg = mul(x,yg)
  #yh = mul(x,yh)
  #yt = mul(x,yt)
  sp = SimplePlot()
  pv = sp.addPoints(sx,yg); pv.setLineColor(Color.BLACK); pv.setLineWidth(3)
  pv = sp.addPoints(sx,yh); pv.setLineColor(Color.RED); pv.setLineWidth(3)
  pv = sp.addPoints(sx,yt); pv.setLineColor(Color.GREEN); pv.setLineWidth(3)

def plot(x,ymin=0.0,ymax=0.0):
  sp = SimplePlot.asPoints(x)
  if ymin<ymax:
    sp.setVLimits(ymin,ymax)

def makeBlocks(n1):
  nb = 17
  db = 1.0
  pb = 8.0
  sb = 1.0
  xb = pb
  m1 = 1+n1/nb
  x = zerofloat(n1)
  for i1 in range(n1):
    if (i1+1)%m1==0:
      pb = pb-db
      sb = -sb
    x[i1] = sb*pb
  return x

def makeRandom(seed,scale,sigma,n1):
  r = Random(seed)
  x = mul(2.0*scale,sub(randfloat(r,n1),0.5))
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply0(x,x)
  return x

def makeImpulse(n1):
  x = zerofloat(n1)
  x[n1/2] = 1.0
  return x

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
