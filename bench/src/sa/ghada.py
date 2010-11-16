# Process Ghada's CT image

import sys
from java.awt import *
from java.lang import *
from java.io import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

dataDir = "/data/sa/ghada/"
n1,d1,f1 = 401,1.0,0.0
n2,d2,f2 = 401,1.0,0.0
n3,d3,f3 = 401,1.0,0.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)

def main(args):
  #convert()
  #threshold()
  slice2()

def slice2():
  x = readImage("ct")
  x = x[n3/2]
  writeImage("ct"+str(n3/2),x)
  display2(x)
  stats(x)
  y = sub(x,0.50)
  y = sgn(y)
  y = clip(0.0,1.0,y)
  display2(y)

def threshold():
  x = readImage("ct")
  stats(x)
  y = sub(x,0.50)
  y = sgn(y)
  y = clip(0.0,1.0,y)
  display3(x)
  display3(y)

def stats(x):
  print "min =",min(x)," max =",max(x)
  h = Histogram(flatten(x))
  SimplePlot.asPoints(h.getBinSampling(),h.getDensities())

def convert():
  x = readGhadaImage()
  x = copy(n1,n2,n3,x)
  x = mul(1.0/32768.0,x)
  writeImage("ct",x)

def readGhadaImage():
  name = "S_52_B_size_700X700X700.bin"
  x = zerofloat(700,700,700)
  ais = ArrayInputStream(dataDir+name,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x

def readImage(name):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+name+".dat")
  ais.readFloats(x)
  ais.close()
  x = copy(n1,n2,n3,x)
  return x

def writeImage(name,x):
  aos = ArrayOutputStream(dataDir+name+".dat")
  aos.writeFloats(x)
  aos.close()

def display2(x):
  sp = SimplePlot()
  sp.setSize(929,823)
  sp.addColorBar()
  pv = sp.addPixels(x)
  pv.setClips(0.0,1.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)

def display3(x):
  sf = SimpleFrame.asImagePanels(x)
  sf.orbitView.setScale(2.0)
  sf.setSize(900,900)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

