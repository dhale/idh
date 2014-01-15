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


def main(args):
  #slopesLsf()
  slopesPsf()

def slopesPsf():
  f = readImage(ffile)
  psf = PefSlopeFinder(-2,2)
  psf.setSmoothness(1,1)
  p0 = zerofloat(n1,n2)
  e0 = psf.applyPef(p0,f)
  p1 = copy(p0)
  psf.updateSlopes(f,p1)
  e1 = psf.applyPef(p1,f)
  p2 = copy(p1)
  psf.updateSlopes(f,p2)
  e2 = psf.applyPef(p2,f)
  p3 = copy(p2)
  psf.updateSlopes(f,p3)
  e3 = psf.applyPef(p3,f)
  fmin,fmax = -1.0,1.0
  emin,emax = -0.2,0.2
  pmin,pmax = -0.5,0.5
  g = sub(f,e3)
  plot(f,cmin=fmin,cmax=fmax,title="input image")
  plot(g,cmin=fmin,cmax=fmax,title="predicted image")
  plot(e0,cmin=emin,cmax=emax,title="prediction errors: 0")
  plot(e1,cmin=emin,cmax=emax,title="prediction errors: 1")
  plot(e2,cmin=emin,cmax=emax,title="prediction errors: 2")
  plot(e3,cmin=emin,cmax=emax,title="prediction errors: 3")
  plot(p1,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: 1")
  plot(p2,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: 2")
  plot(p3,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: 3")

def slopesLsf():
  f = readImage(ffile)
  pmax = 10.0
  sigma1 = 20.0
  sigma2 = 1.0
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,el)
  plot(f,title="f")
  plot(p2,cmap=jet,title="p2")
  plot(el,cmap=jet,title="el")

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
