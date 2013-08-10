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
from fault import *
from util import FakeData

seismicDir = "/data/seis/tpd/csm/oldslices/"
ffile = "tp73"
#s1 = Sampling(251,0.004,0.500)
#s2 = Sampling(357,0.025,0.000)
s1 = Sampling(251,1.0,0.0)
s2 = Sampling(357,1.0,0.0)
n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta

def main(args):
  flatten()

def flatten():
  k1s = [[ 44, 40],[190,181],[160,157]]
  k2s = [[210,260],[ 90,190],[120,180]]
  #f = FakeData.seismic2d2011A(n1,n2,30)
  f = readImage(ffile)
  plot(s1,s2,f,title="Input",png="f")
  plot(s1,s2,f,c=(k1s,k2s),title="Constraints",png="c")
  #sigma = 1.0 # good for fake data
  sigma = 8.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(sigma,pmax)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,wp)
  p2 = mul(d1/d2,p2)
  fl = faults(f)
  wp = pow(sub(1.0,fl),8)
  plot(s1,s2,wp,cmap=jet,title="Weights",png="w")
  plot(s1,s2,p2,cmap=jet,title="Slopes",cmin=-0.5*d1/d2,cmax=0.5*d1/d2,png="p")
  fl = Flattener2C()
  fl.setWeight1(0.02)
  fl.setIterations(0.01,1000)
  fl.setSmoothings(4.0,8.0)
  for cs in [(k1s,k2s),None]:
    if cs:
      nc = len(cs[0])
      psuffix = str(nc)
      tsuffix = " ("+str(nc)+" constraints)"
    else:
      psuffix = "0"
      tsuffix = " (no constraints)"
    fm = fl.getMappingsFromSlopes(s1,s2,p2,wp,cs)
    g = fm.flatten(f)
    h = fm.unflatten(g)
    s = fm.getShiftsS()
    plot(s1,s2,f,u=fm.u1,title="Horizons"+tsuffix,png="fu"+psuffix)
    plot(s1,s2,g,title="Flattened"+tsuffix,png="g"+psuffix)
    plot(s1,s2,h,title="Unflattened"+tsuffix,png="h"+psuffix)
    #plot(s1,s2,s,cmap=jet,title="Shifts"+tsuffix,png="s"+psuffix)
    print "average shift =",sum(s)/(n1*n2),"samples"

def flip2(f):
  n2 = len(f)
  for i2 in range(n2/2):
    fi2 = f[i2]
    f[i2] = f[n2-1-i2]
    f[n2-1-i2] = fi2

def faults(f):
  fs = FaultSemblance()
  p2 = fs.slopes(f)
  snd = fs.semblanceNumDen(p2,f)
  fs = FaultScanner2(20,snd,FaultScanner2.Smoother.SHEAR)
  (fl,ft) = fs.scan(-20,20)
  #(fl,ft) = fs.thin((fl,ft))
  return fl

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
pngDir = None
#pngDir = "./png/"
def plot(s1,s2,x,u=None,c=None,cmap=ColorMap.GRAY,clab=None,cmin=0,cmax=0,
         title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(600,900)
  sp.plotPanel.setColorBarWidthMinimum(80)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    cv = sp.addContours(s1,s2,u)
    cv.setContours(50)
    cv.setLineColor(Color.YELLOW)
  if c:
    k1s,k2s = c
    pv = PointsView(k1s,k2s)
    pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
    pv.setLineColor(Color.WHITE)
    pv.setLineWidth(2)
    pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
    pv.setMarkColor(Color.WHITE)
    pv.setMarkSize(10)
    sp.add(pv)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")

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
