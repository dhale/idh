import sys
from math import *

from java.io import *
from java.lang import *
from java.nio import *
from java.awt.image import *
from javax.imageio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *


def main(args):
  doTp()

def doTp():
  n1,n2=251,357
  dataDir = "/data/seis/tp/"
  f = readFloats(n1,n2,dataDir+"tp73.dat")
  b = makeBytes(f)
  p = makeImage(n2,n1,b)
  writeImage(p,dataDir+"tp73.png")

def doAtw():
  n1,n2=500,500
  dataDir = "/data/seis/atw/"
  f = readFloats(n1,n2,dataDir+"atwj1s.dat")
  b = makeBytes(f)
  p = makeImage(n2,n1,b)
  writeImage(p,dataDir+"atwj1s.png")

def readFloats(n1,n2,fileName):
  ais = ArrayInputStream(fileName)
  x = zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

def makeBytes(f,cmin=0.0,cmax=0.0):
  f = transpose(f)
  n1,n2 = len(f[0]),len(f)
  nx,ny = n1,n2
  n = n1*n2
  b = zerobyte(n)
  if cmin==cmax:
    cmin = min(f)
    cmax = max(f)
  s = 255.0/(cmax-cmin)
  i = 0
  for i2 in range(n2):
    for i1 in range(n1):
      ci = s*(f[i2][i1]-cmin);
      if ci<0.0:
        bi = 0
      elif ci>255.0:
        bi = 255
      else:
        bi = int(ci)
      if bi>127:
        bi -= 256
      b[i] = bi
      i += 1
  return b

def makeImage(nx,ny,b):
  bi = BufferedImage(nx,ny,BufferedImage.TYPE_BYTE_GRAY)
  wr = bi.getRaster()
  wr.setDataElements(0,0,nx,ny,b)
  return bi

def writeImage(bi,fileName):
  file = File(fileName)
  ImageIO.write(bi,"png",file)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
