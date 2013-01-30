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

seismicDir = "/data/seis/tpd/csm/oldslices/"
ffile = "tp73"
s1 = Sampling(251,0.004,0.500)
s2 = Sampling(357,0.025,0.000)
n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta

def main(args):
  flatten()

def flatten():
  #f = FakeData.seismic2d2011A(n1,n2,30)
  f = readImage(ffile)
  #sigma = 1.0 # good for fake data
  sigma = 8.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(sigma,pmax)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,el)
  p2 = mul(d1/d2,p2)
  el = pow(el,6)
  plot(s1,s2,el,cmap=jet)
  plot(s1,s2,p2,cmap=jet,cmin=-0.1,cmax=0.1)
  fl = Flattener2()
  fm = fl.getMappingsFromSlopes(s1,s2,p2,el)
  g = fm.flatten(f)
  h = fm.unflatten(g)
  s = fm.getShiftsS()
  plot(s1,s2,f,u=fm.u1)
  plot(s1,s2,g)
  plot(s1,s2,h)
  plot(s1,s2,s,cmap=jet)
  print "average shift =",sum(s)/(n1*n2),"samples"

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(s1,s2,x,u=None,cmap=ColorMap.GRAY,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar();
  sp.setSize(600,900)
  sp.plotPanel.setColorBarWidthMinimum(80)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    cv = sp.addContours(s1,s2,u)
    cv.setLineColor(Color.YELLOW)

#############################################################################
# utilities

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
