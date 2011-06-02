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
from dnp.FlattenerVS import *
from util import FakeData

#seismicDir = "/data/seis/tp/csm/oldslices/"
seismicDir = "/Users/dhale/Home/box/jtk/trunk/data/"
ffile = "tp73"
#s1 = Sampling(251,0.004,0.500)
#s2 = Sampling(357,0.025,0.000)
#n1,n2 = s1.count,s2.count

s1 = Sampling(251,1.000,0.500)
s2 = Sampling(501,1.000,0.000)
n1,n2 = s1.count,s2.count

def main(args):
  showFake()
  #flatten()
  #flattenTest()
  #slopes()

def showFake():
  clip = 0.85
  f,g,s1,s2,r1,r2 = FakeData.seismicAndShifts2d2011A(n1,n2,45)
  plot(f,cmin=-clip,cmax=clip)
  plot(g,cmin=-clip,cmax=clip)
  plot(r1,jet)
  plot(r2,jet)
  d = getDeterminantsFromShifts([r1,r2])
  plot(d,jet)
  a = getAFromShifts([r1,r2])
  plot(a,jet,cmin=-1,cmax=2)
  a = applyInverseShiftsL([s1,s2],a)
  plot(a,jet,cmin=-1,cmax=2)


def flatten():
  clip = 0.85
  #f = readImage(ffile)
  #sigma = 8.0
  f,g,s1,s2,r1,r2 = FakeData.seismicAndShifts2d2011A(n1,n2,45)
  f = g
  d = getDeterminantsFromShifts([r1,r2])
  plot(f,cmin=-clip,cmax=clip)
  #plot(r1,jet)
  #plot(r2,jet)
  #plot(d,jet,cmin=0.65,cmax=1.35)
  sigma,pmax = 1.0,10.0
  lsf = LocalSlopeFinder(sigma,pmax)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,el)
  el = pow(el,6)
  sigma1,sigma2 = 12.0,12.0
  fl = FlattenerVS(sigma1,sigma2)
  #plot(el,gray)
  #plot(p2,gray,-1,1)
  slist = []
  for rotate in [0.0,1.0]:
    s = fl.findShifts(rotate,p2,el)
    #s = fl.findShiftsA(p2,el)
    g = fl.applyShifts(s,f)
    plot(g,cmin=-clip,cmax=clip)
    s1,s2 = s[0],s[1]
    #plot(s1,jet)
    #plot(s2,jet)
    #d = getDeterminantsFromShifts(s)
    #plot(d,jet,cmin=0.65,cmax=1.35)
    print "average s1 =",sum(s1)/n1/n2,"samples"
    print "average s2 =",sum(s2)/n1/n2,"samples"
    slist.append(s)
  sv,sr = slist[0],slist[1]
  a = getAFromShifts(sr,sv)
  print "a min =",min(a),"max =",max(a)
  plot(a,jet)

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
  fl = FlattenerCg(8.0,0.01)
  sf = fl.findShifts(p2,el) # found shifts
  se = neg(mul(t,asinbx)) # exact shifts
  plot(sf,jet,-smax,smax)
  plot(se,jet,-smax,smax)

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(x,cmap=ColorMap.GRAY,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(600,900)
  sp.setSize(1530,770)
  sp.plotPanel.setColorBarWidthMinimum(80)
  sp.addColorBar();
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
