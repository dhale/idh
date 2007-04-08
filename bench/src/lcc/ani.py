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

#############################################################################
# functions


def main(args):
  testInverseForward()

def testInverseForward():
  n = 315
  x = Array.zerofloat(n)
  y = Array.zerofloat(n)
  z = Array.zerofloat(n)
  x[n/2] = 1.0
  sigma1 = 4
  sigma2 = sqrt(sigma1*sigma1-2)
  f1 = makeFilter(sigma1)
  f1.applyInverse(x,y)
  f1.applyInverseTranspose(y,y)
  f2 = makeFilter(sigma2)
  f2.applyTranspose(y,z)
  f2.apply(z,z)
  sp = SimplePlot()
  #sp.addPoints(x)
  sp.addPoints(y)
  sp.addPoints(z)

def makeFilter(sigma):
  dr1 = 1.12075
  di1 = 1.27788
  dr2 = 1.76952
  di2 = 0.46611
  r1 = 1.0/sqrt(dr1*dr1+di1*di1)
  r2 = 1.0/sqrt(dr2*dr2+di2*di2)
  s1 = atan2(di1,dr1)
  s2 = atan2(di2,dr2)
  a1 = -2.0*pow(r1,2.0/sigma)*cos(2.0*s1/sigma)
  a2 = pow(r1,4.0/sigma)
  b1 = -2.0*pow(r2,2.0/sigma)*cos(2.0*s2/sigma)
  b2 = pow(r2,4.0/sigma)
  lag1 = (0,1,2,3,4)
  a = Array.zerofloat(5)
  Array.copy((1.0,a1+b1,a2+a1*b1+b2,a1*b2+b1*a2,a2*b2),a)
  a = Array.mul(1.0/Array.sum(a),a)
  return CausalFilter(lag1,a)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
