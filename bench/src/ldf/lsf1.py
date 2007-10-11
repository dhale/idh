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
sigma = 16
small = 0.01
niter = 80

#############################################################################
# functions

def main(args):
  #goTestSymmetric()
  #goTestTranspose()
  #goSmooth()
  goSmooth1()
  return

def goTestTranspose():
  # A is symmetric, then x'(Ay) = (Ay)'x = y'(A'x)
  n = 11
  x = makeRandom(n)
  y = makeRandom(n)
  ds = Array.randfloat(n)
  ax = Array.zerofloat(n)
  ay = Array.zerofloat(n)
  lsf = LocalSmoothingFilter1(sigma,n,ds)
  lsf.apply(y,ay)
  lsf.applyTranspose(x,ax)
  xay = 0.0
  yax = 0.0
  for i in range(n):
    xay += x[i]*ay[i]
    yax += y[i]*ax[i]
  print "xay =",xay," yax =",yax

def goTestSymmetric():
  # if A is symmetric, then x'Ay = (Ay)'x = y'A'x = y'Ax
  n = 11
  x = makeRandom(n)
  y = makeRandom(n)
  ax = Array.zerofloat(n)
  ay = Array.zerofloat(n)
  lsf = LocalSmoothingFilter(sigma)
  lsf.applyPass(None,x,ax)
  lsf.applyPass(None,y,ay)
  xay = 0.0
  yax = 0.0
  for i in range(n):
    xay += x[i]*ay[i]
    yax += y[i]*ax[i]
  print "xay =",xay," yax =",yax

def goSmooth():
  #x = makeRandom(n1)
  x = makeImpulse(n1)
  y = Array.zerofloat(n1)
  #ds = None
  ds = makeBlock(n1)
  #ds = makeStepUp(n1/2+1,n1)
  #ds = makeStepDown(n1/2+1,n1)
  lsf = LocalSmoothingFilter(sigma)
  lsf.applyPass(ds,x,y)
  SimplePlot.asSequence(x);
  SimplePlot.asSequence(y);

def goSmooth1():
  #x = makeRandom(n1)
  #x = makeImpulse(n1)
  x = Array.zerofloat(n1)
  x[5] = x[n1-6] = x[n1/2] = 1.0
  y = Array.zerofloat(n1)
  #ds = None
  ds = makeBlock(n1)
  #ds = makeStepUp(n1/2+1,n1)
  #ds = makeStepDown(n1/2+1,n1)
  lsf = LocalSmoothingFilter1(sigma,n1,ds)
  lsf.applyTranspose(x,y)
  lsf.apply(y,y)
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
