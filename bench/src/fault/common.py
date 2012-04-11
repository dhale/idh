import sys
from org.python.util import PythonObjectInputStream
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fault import *

#############################################################################
# globals
n1,n2,n3 = 0,0,0
s1,s2,s3 = None,None,None
dataDir,dataSub = None,None
 
#############################################################################
# utilities

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))
 
#############################################################################
# data read/write

def readImage(fileName):
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+dataSub+fileName+".dat")
  ais.readFloats(f)
  ais.close()
  print "readImage: min =",min(f)," max =",max(f)
  return f
def writeImage(f,fileName):
  aos = ArrayOutputStream(dataDir+dataSub+fileName+".dat")
  aos.writeFloats(f)
  aos.close()
def samplingS1():
  global n1,n2,n3
  global s1,s2,s3
  global dataDir,dataSub
  dataDir = "/data/seis/f3d/faults/"
  dataSub = "s1/"
  n1,n2,n3 = 222,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.964,1.250,2.500
  s1,s2,s3 = samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3)
  return s1,s2,s3
def samplingS1A():
  global n1,n2,n3
  global s1,s2,s3
  global dataDir,dataSub
  dataDir = "/data/seis/f3d/faults/"
  dataSub = "s1a/"
  n1,n2,n3 = 90,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 1.484,1.250,2.500
  s1,s2,s3 = samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3)
  return s1,s2,s3
def samplingS1B():
  global n1,n2,n3
  global s1,s2,s3
  global dataDir,dataSub
  dataDir = "/data/seis/f3d/faults/"
  dataSub = "s1b/"
  n1,n2,n3 = 90,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 1.024,1.250,2.500
  s1,s2,s3 = samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3)
  return s1,s2,s3
def samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3):
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  print "s1:",s1.count,s1.delta,s1.first
  print "s2:",s2.count,s2.delta,s2.first
  print "s3:",s3.count,s3.delta,s3.first
  return s1,s2,s3
def sampling():
  return dataSub[:-1]

#############################################################################
# plotting

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def jetNotch(alpha):
  return ColorMap.setAlpha(ColorMap.JET,
    clip(0.0,1.0,rampfloat(0.0,25.0*alpha/256,256)))

  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
    """
    if i<96:
      a[i] = 1.0
    elif i<128:
      a[i] = alpha*(128.0-i)/32.0
    elif i<160:
      a[i] = alpha*(i-127.0)/32.0
    else:
      a[i] = 1.0
    """
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
    """
    if i<96:
      a[i] = 1.0
    elif i<128:
      a[i] = alpha*(128.0-i)/32.0
    elif i<160:
      a[i] = alpha*(i-127.0)/32.0
    else:
      a[i] = 1.0
    """
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
