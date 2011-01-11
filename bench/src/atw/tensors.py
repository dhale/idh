import sys

from java.awt import *
from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

# John Mathewson's subsets, resampled to 20 x 20 m trace spacing
# atwj1.dat = John's area1 file
# atwj1s.dat = horizontal slice of John's atwj1 for i1=40
# atwj3.dat = John's area3 file
#n1= 129; d1=0.0040; f1=0.0000
#n2= 500; d2=0.0200; f2=0.0000
#n3= 500; d3=0.0200; f3=0.0000

def main(args):
  goTensors()

def goTensors():
  s1 = Sampling(500,0.02,0.0)
  s2 = Sampling(500,0.02,0.0)
  g = readImage("atwj1s",s1,s2)
  d1 = makeTensors1(g)
  d2 = makeTensors2(g)
  d3 = makeTensors3(g)
  plot(g,s1,s2,d1,dscale=2)
  plot(g,s1,s2,d2,dscale=1)
  plot(g,s1,s2,d3,dscale=2)

def makeTensors1(g):
  lof = LocalOrientFilter(4.0)
  d = lof.applyForTensors(g)
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av) # eigenvalues are gradients squared; au >= av
  au = pow(au,1.0) # amplify gradients to emphasize edges in image
  av = pow(av,1.0) # amplify gradients to emphasize edges in image
  amax = max(au) # largest eigenvalue
  aeps = 0.000001*amax # a relatively small eigenvalue
  au = add(aeps,au) # avoid divide by zero
  av = add(aeps,av) # avoid divide by zero
  amin = min(av) # smallest eigenvalue
  au = div(amin,au) # au <= av <= 1
  av = div(amin,av) # au <= av <= 1
  d.setEigenvalues(au,av)
  return d

def makeTensors2(g):
  lof = LocalOrientFilter(4.0)
  d = lof.applyForTensors(g)
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av) # eigenvalues are gradients squared; au >= av
  au = pow(au,1.0) # amplify gradients to emphasize edges in image
  av = pow(av,1.0) # amplify gradients to emphasize edges in image
  amax = max(au) # largest eigenvalue
  aeps = 0.000001*amax # a relatively small eigenvalue
  au = add(aeps,au) # avoid divide by zero
  av = add(aeps,av) # avoid divide by zero
  au = div(av,au) # au <= 1
  av = div(av,av) # av = 1
  d.setEigenvalues(au,av)
  return d

def makeTensors3(g):
  lof = LocalOrientFilter(4.0)
  d = lof.applyForTensors(g)
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av) # eigenvalues are gradients squared; au >= av
  amax = max(au) # largest eigenvalue
  aeps = 0.000001*amax # a relatively small eigenvalue
  au = add(aeps,au) # avoid divide by zero
  av = add(aeps,av) # avoid divide by zero
  amin = min(av) # smallest eigenvalue
  aiso = div(av,au) # isotropy
  au = div(amin,au) # au <= av <= 1
  av = div(amin,av) # au <= av <= 1
  au = mul(aiso,au) # scale by isotropy
  d.setEigenvalues(au,av)
  return d

#############################################################################
# plotting

#pngDir = "png/atw/" # where to put PNG images of plots
pngDir = None # for no PNG images

def plot(g,s1,s2,d=None,dscale=1,cmin=0,cmax=0,png=None):
  sp = SimplePlot()
  sp.setHLabel("Crossline distance (km)")
  sp.setVLabel("Inline distance (km)")
  sp.setHInterval(2.0)
  sp.setVInterval(2.0)
  sp.setFontSizeForPrint(8,240)
  sp.setSize(910,945)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(1)
    tv.setEllipsesDisplayed(30)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir and png:
    sp.paintToPng(600,3,pngDir+png+".png")

#############################################################################
# data input/output

dataDir = "/data/seis/atw/"

def readImage(fileName,s1,s2):
  n1,n2 = s1.count,s2.count
  ais = ArrayInputStream(dataDir+fileName+".dat")
  x = zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(dataDir+fileName+".dat")
  aos.writeFloats(x)
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
