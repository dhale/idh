import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from dnp import *

seismicDir = "/data/seis/tp/csm/oldslices/"
ffile = "tp73"
s1 = Sampling(251,0.004,0.500)
s2 = Sampling(357,0.025,0.000)
n1,n2 = s1.count,s2.count


def main(args):
  flatten()

def flatten():
  f = readImage(ffile)
  plot(f)
  #for fl in [FlattenerF(8.0,0.1),Flattener(8.0)]:
  for fl in [Flattener(8.0)]:
    s = fl.findShifts(f)
    g = fl.applyShifts(f,s)
    plot(g)
    plot(s,jet)

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(x,cmap=ColorMap.GRAY):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar();
  sp.setSize(1012,676)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)

#############################################################################
# read/write files

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
