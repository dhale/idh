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

n1 = 201
sigma = 1.0
small = 0.001
niter = 100

#############################################################################
# functions

def main(args):
  goInterp()
  return

def goInterp():
  f = Array.zerobyte(n1)
  x = Array.zerofloat(n1)
  f[1*n1/4] = 1; x[1*n1/4] = 0.5
  f[2*n1/4] = 1; x[2*n1/4] = 1.5
  f[3*n1/4] = 1; x[3*n1/4] = 0.5
  y = Array.copy(x)
  ds = None
  #ds = makeBlock(n1)
  #ds = makeStepUp(n1/2+1,n1)
  #ds = makeStepDown(n1/2+1,n1)
  lif = LocalInterpolationFilter1(sigma,small,niter)
  lif.apply(ds,f,y)
  SimplePlot.asSequence(x);
  SimplePlot.asSequence(y);

def makeBlock(n1):
  ds = Array.fillfloat(0.3,n1)
  for i1 in range(1*n1/3,2*n1/3):
    ds[i1] = 1.0
  return ds

def makeStepUp(k1,n1):
  ds = Array.zerofloat(n1)
  for i1 in range(k1,n1):
    ds[i1] = 1.0
  return ds

def makeStepDown(k1,n1):
  ds = Array.fillfloat(1.0,n1)
  for i1 in range(k1,n1):
    ds[i1] = 0.0
  return ds

def makeImpulse(n1):
  x = Array.zerofloat(n1)
  x[n1/2] = 1.0
  return x

def makeRandom(n1):
  #r = Random(314159)
  #return Array.sub(Array.randfloat(r,n1),0.5)
  return Array.sub(Array.randfloat(n1),0.5)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
