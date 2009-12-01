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
  filter()
  #slopes()

def filter():
  f = readImage(ffile)
  pmax = 8.0
  sigma1 = pmax
  sigma2 = 20.0
  gp = zerofloat(n1,n2)
  gm = zerofloat(n1,n2)
  pp = zerofloat(n1,n2)
  pm = zerofloat(n1,n2)
  cp = zerofloat(n1,n2)
  cm = zerofloat(n1,n2)
  lsf = LocalSlopeFinderS(sigma1,sigma2,pmax)
  lsf.filter( 1,f,gp,pp,cp)
  lsf.filter(-1,f,gm,pm,cm)
  a = (sigma2-0.5*sqrt(2))/(sigma2+0.*sqrt(2))
  g = mul(1/(1+a),sub(add(gp,gm),mul(1-a,f)))
  plot(f,gray,-4,4,"input")
  plot(g,gray,-4,4,"filtered symmetric")
  plot(gp,gray,-4,4,"filtered left to right")
  plot(gm,gray,-4,4,"filtered right to left")
  plot(pp,jet,-1,1,"slopes left to right")
  plot(pm,jet,-1,1,"slopes right to left")
  plot(cp,jet,0,1,"correlations left to right")
  plot(cm,jet,0,1,"correlations right to left")

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

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(x,cmap=ColorMap.GRAY,cmin=0,cmax=0,title=""):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
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
