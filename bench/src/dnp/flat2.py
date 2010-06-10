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

#seismicDir = "/data/seis/tp/csm/oldslices/"
seismicDir = "/Users/dhale/Home/box/jtk/trunk/data/"
ffile = "tp73"
s1 = Sampling(251,0.004,0.500)
s2 = Sampling(357,0.025,0.000)
n1,n2 = s1.count,s2.count


def main(args):
  #slopes()
  flatten()
  #flattenTest()

def slopes():
  f = readImage(ffile)
  pmax = 10.0
  sigma1 = 20.0
  sigma2 = 10.0
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf = LocalSlopeFinderS(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,el)
  plot(f)
  plot(p2,jet)
  plot(el,jet)

def flatten():
  f = readImage(ffile)
  plot(f)
  sigma = 8.0
  pmax = 10.0
  lsf = LocalSlopeFinder(sigma,pmax)
  sigma1 = 6.0
  sigma2 = 12.0
  for fl in [FlattenerCg(sigma1,sigma2)]:
    p2 = zerofloat(n1,n2)
    el = zerofloat(n1,n2)
    lsf.findSlopes(f,p2,el)
    el = pow(el,6)
    #plot(el,gray)
    #plot(p2,gray,-1,1)
    s = fl.findShifts(p2,el)
    g = fl.applyShifts(f,s)
    plot(g)
    plot(s,jet)
    print "average shift =",sum(s)/(n1*n2),"samples"

def flattenTest():
  """Test for t(tau,x) = tau*(1+a*sin(b*x))"""
  x = rampfloat(0,0,1,n1,n2)
  t = rampfloat(0,1,0,n1,n2)
  smax = 5.0
  a = smax/(n1-1)
  b = 2*PI/(n2-1)
  bx = mul(b,x)
  bt = mul(b,t)
  cosbx = cos(bx)
  sinbx = sin(bx)
  acosbx = mul(a,cosbx)
  asinbx = mul(a,sinbx)
  p2 = div(mul(bt,acosbx),add(1,asinbx))
  el = fillfloat(1,n1,n2)
  fl = FlattenerS(8.0,0.01)
  sf = fl.findShifts(p2,el) # found shifts
  se = neg(mul(t,asinbx)) # exact shifts
  plot(sf,jet,-smax,smax)
  plot(se,jet,-smax,smax)

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(x,cmap=ColorMap.GRAY,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar();
  sp.setSize(600,900)
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
