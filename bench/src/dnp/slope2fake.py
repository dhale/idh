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
from util import FakeData

s1 = Sampling(501,1.000,0.000)
s2 = Sampling(501,1.000,0.000)
n1,n2 = s1.count,s2.count
fmin,fmax = -3.0,3.0
emin,emax = -3.0,3.0
pmin,pmax = -3.0,3.0
seismicDir = "/Users/dhale/Desktop/"

def main(args):
  #slopesLsf()
  #slopesPwd()
  #compare()
  sergey()

def sergey():
  for nrms in [0.0,0.5,1.0]:
    f,p = FakeData.seismicAndSlopes2d2014A(nrms)
    plot(f,cmin=fmin,cmax=fmax,title="input image")
    plot(p,cmin=pmin,cmax=pmax,title="known slopes")
    fname = "f"+str(int(nrms*10))
    pname = "p"
    writeImage(fname,f)
    writeImage(pname,p)

def compare():
  sigma1,sigma2 = 12,2
  eps1,eps2 = 2.0,0.5
  nrms = 0.50
  f,p = FakeData.seismicAndSlopes2d2014A(nrms)
  pwd = PlaneWaveDestructor(-8,8)
  pwd.setOuterIterations(10) # 5 is the default
  pwd.setSmoothness(eps1,eps2)
  pa = pwd.findSlopes(f)
  ea = sub(pa,p)
  ga = pwd.applyFilter(pa,f)
  lsf = LocalSlopeFinder(sigma1,sigma2,8)
  pb = zerofloat(n1,n2)
  lsf.findSlopes(f,pb,None)
  eb = sub(pb,p)
  gb = pwd.applyFilter(pb,f)
  plot(p ,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: syn")
  plot(pa,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: pwd")
  plot(pb,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: lsf")
  plot(ea,cmin=0.5*pmin,cmax=0.5*pmax,cmap=jet,title="errors: pwd")
  plot(eb,cmin=0.5*pmin,cmax=0.5*pmax,cmap=jet,title="errors: lsf")
  plot(ga,cmin=emin,cmax=emax,title="output image: pwd")
  plot(gb,cmin=emin,cmax=emax,title="output image: lsf")
  plot(f,cmin=fmin,cmax=fmax,title="input image")

def slopesPwd():
  f,p = FakeData.seismicAndSlopes2d2014A(nrms)
  pmax = max(abs(min(p)),abs(max(p)))
  pwd = PlaneWaveDestructor(pmin,pmax)
  pwd.setSmoothness(1.5,0.5)
  pwd.setLateralBias(0.0) # 0.0 is the default
  pwd.setOuterIterations(10) # 5 is the default
  q = pwd.findSlopes(f)
  g = pwd.applyFilter(q,f)
  plot(f,cmin=fmin,cmax=fmax,title="input image")
  plot(g,cmin=emin,cmax=emax,title="output image")
  plot(p,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: pwd")

def slopesLsf():
  f = readImage(ffile)
  sigma1 = 3.0
  sigma2 = 1.0
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,el)
  #plot(f,cmin=fmin,cmax=fmax,title="input image")
  plot(p2,cmin=pmin,cmax=pmax,cmap=jet,title="slopes: lsf")
  #plot(el,cmap=jet,title="linearities: lsf")

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
