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
  #goShowEdgeFilters()
  #goShowRangeFunctions()
  goFilter()
  return

def goShowEdgeFilters():
  n = 201
  fft = Fft(n)
  x = zerofloat(n)
  x[n/2] = 1.0
  lsf = LocalSmoothingFilter()
  y1 = zerofloat(n)
  y2 = zerofloat(n)
  sigma1 = 1.0
  sigma2 = 2.0
  lsf.apply(0.5*sigma1*sigma1,x,y1)
  lsf.apply(0.5*sigma2*sigma2,x,y2)
  y = sub(y1,y2)
  ay = cabs(fft.applyForward(y))
  #sp = SimplePlot.asPoints(x)
  #sp = SimplePlot.asPoints(y)
  SimplePlot.asPoints(ay).setVLimits(0.0,0.8)
  rgf1 = RecursiveGaussianFilter(sigma1)
  rgf2 = RecursiveGaussianFilter(sigma2)
  rgf1.apply0(x,y1)
  rgf2.apply0(x,y2)
  y = sub(y1,y2)
  ay = cabs(fft.applyForward(y))
  #sp = SimplePlot.asPoints(y)
  SimplePlot.asPoints(ay).setVLimits(0.0,0.8)

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

def goFilter():
  n1 = 401
  x = makeStepUp(n1/2,n1)
  x = add(x,mul(1.0,mul(0.5,randfloat(n1))))
  plot(x)
  y = zerofloat(n1)
  sigmaS = 20.0
  sigmaX = 0.5
  bf = BilateralFilter(sigmaS,sigmaX)
  for type in [gauss,huber,tukey]:
    bf.setType(type)
    bf.apply(x,y)
    plot(y)

def plot(x):
  sp = SimplePlot.asPoints(x)
  sp.setVLimits(0.0,2.0)

def makeBlock(n1):
  ds = fillfloat(0.3,n1)
  for i1 in range(1*n1/3,2*n1/3):
    ds[i1] = 1.0
  return ds

def makeStepUp(k1,n1):
  ds = zerofloat(n1)
  for i1 in range(k1,n1):
    ds[i1] = 1.0
  return ds

def makeStepDown(k1,n1):
  ds = fillfloat(1.0,n1)
  for i1 in range(k1,n1):
    ds[i1] = 0.0
  return ds

def makeImpulse(n1):
  x = zerofloat(n1)
  x[n1/2] = 1.0
  return x

def makeRandom(n1):
  #r = Random(314159)
  #return sub(randfloat(r,n1),0.5)
  return sub(randfloat(n1),0.5)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
