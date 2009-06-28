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

n1 = 301
sigma = 0.5
small = 0.001
niter = 100

#############################################################################
# functions

def main(args):
  goInterp()
  return

def goInterp():
  f = zerobyte(n1)
  x = zerofloat(n1)
  f[1*n1/4] = 1; x[1*n1/4] = 0.5
  f[2*n1/4] = 1; x[2*n1/4] = 1.5
  f[3*n1/4] = 1; x[3*n1/4] = 0.5
  y = copy(x)
  #ds = fillfloat(1.0,n1)
  #ds = randfloat(n1)
  ds = makeRaisedCosine(n1)
  #lif = LocalInterpolationFilter1(sigma,0.001,niter)
  #lif.apply(ds,f,y)
  lif = LocalInterpolationFilter1(sigma,0.01,niter)
  lif.applyExact(ds,f,y)
  SimplePlot.asSequence(y);

def reset(f,x,y):
  n = len(f)
  for i in range(n):
    if f[i]:
      y[i] = x[i]

def makeRaisedCosine(n1):
  ds = zerofloat(n1)
  for i1 in range(n1):
    ds[i1] = 0.51+0.49*cos(pi*(2*i1-n1)/n1)
  return ds

def makeBlocky(n1):
  ds = fillfloat(0.2,n1)
  for i1 in range(1*n1/3,2*n1/3):
    ds[i1] = 1.0
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
