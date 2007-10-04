import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

from ldf import *

True = 1
False = 0

#############################################################################
# parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

n1 = 315
n2 = 315
aniso = 100
small = 0.001
niter = 1000
lof = LocalOrientFilter(8)
lof.setGradientSmoothing(1)

#############################################################################
# functions

def main(args):
  goInterp()
  return

def goInterp():
  x = readImage()
  n1,n2 = len(x[0]),len(x)
  v1,v2,ds = getV(x)
  f = Array.zerobyte(n1,n2)
  y = Array.zerofloat(n1,n2)
  z = Array.zerofloat(n1,n2)
  plot(x,-10,10,gray,"x")
  #for itest in [0,1,2,3]:
  for itest in [0]:
    if itest==0 or itest==2:
      for i2 in [n2/2]:
        for i1 in range(n1):
          f[i2][i1] = 1
          y[i2][i1] = i1-n1/2
      cmin = -n1/2
      cmax = n1/2
      sigma = 1000
      aniso = 1
    else:
      for i2 in [n2/2]:
        for i1 in range(n1):
          f[i2][i1] = 1
          y[i2][i1] = x[i2][i1]
      cmin = -10
      cmax = 10
      sigma = 1
      aniso = 100
    lif = LocalInterpolationFilter(aniso,small,niter)
    if itest<2:
      ds = None
      es = None
    else:
      ds = makeBlock(n1,n2)
      es = makeBlock(n1,n2)
    z = Array.copy(y)
    """
    lsf = LocalSmoothingFilter(sigma,aniso)
    lsf.applyPass(ds,es,v1,y,z)
    plot(y,0.0,0.0,jet)
    plot(z,0.0,0.0,jet)
    return
    """
    lif.applyLinear(sigma,aniso,ds,es,v1,f,z)
    plot(y,cmin,cmax,jet,"y"+str(itest))
    plot(z,cmin,cmax,jet,"z"+str(itest))
    plot(Array.sub(y,z),0,0,jet)
    plot2(x,z,cmin,cmax,"xz"+str(itest))

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def getV(x):
  n1 = len(x[0])
  n2 = len(x)
  v1 = Array.zerofloat(n1,n2)
  v2 = Array.zerofloat(n1,n2)
  el = Array.zerofloat(n1,n2)
  lof.apply(x,None,None,None,v1,v2,None,None,el)
  return v1,v2,el

def makeBlock(n1,n2):
  ds = Array.zerofloat(n1,n2);
  #for i2 in range(n2/3,n2):
  #  for i1 in range(n1):
  #    ds[i2][i1] = 1.0
  for i1 in range(n1):
    for i2 in range(100-(i1*50)/315,n2):
      ds[i2][i1] = 1.0
  return ds
 
#############################################################################
# plot

def plot(f,cmin=0.0,cmax=0.0,cmap=ColorMap.GRAY,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  frame(p,png)

def plot2(f,g,cmin=0,cmax=0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.getGray())
  pv.setClips(-10,10)
  pv = p.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin!=cmax:
    pv.setClips(cmin,cmax)
  pv.setColorModel(ColorMap.getJet(0.3))
  frame(p,png)

fontSize = 24
#width = 500
#height = 520
width = 700
height = 680
widthColorBar = 80

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  #p = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.LEFT_TOP)
  #p.addColorBar()
  #p.setColorBarWidthMinimum(widthColorBar)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(100,6,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
