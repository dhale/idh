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

seismicDir = "/data/seis/tpd/csm/oldslices/"
ffile = "tp73"
s1 = Sampling(251,0.004,0.500)
s2 = Sampling(357,0.025,0.000)
n1,n2 = s1.count,s2.count
fmin,fmax = -1.0,1.0
emin,emax = -1.0,1.0
pmin,pmax = -0.5,0.5


def main(args):
  slopesLsf()
  slopesPwd()

def slopesPwd():
  f = readImage(ffile)
  pwd = PlaneWaveDestructor(-4.0,4.0)
  pwd.setSmoothness(8,2)
  pwd.setLateralBias(0.0)
  pwd.setOuterIterations(5)
  p = pwd.findSlopes(f)
  g = pwd.applyFilter(p,f)
  plot(f,cmin=fmin,cmax=fmax,title="input image")
  #plot(g,cmin=emin,cmax=emax,title="output image")
  plot(p,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: pwd")

def slopesLsf():
  f = readImage(ffile)
  sigma1 = 4.0
  sigma2 = 2.0
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(sigma1,sigma2,4.0)
  lsf.findSlopes(f,p2,el)
  #plot(f,cmin=fmin,cmax=fmax,title="input image")
  plot(p2,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: lsf")
  plot(el,cmap=jet,title="linearities: lsf")

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(x,cmap=ColorMap.GRAY,cmin=0,cmax=0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
      sp.setTitle(title);
  sp.addColorBar();
  sp.setSize(1000,800)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if cmin<cmax:
    pv.setClips(cmin,cmax)

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
