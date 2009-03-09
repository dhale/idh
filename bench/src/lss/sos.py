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

from lss import *

#############################################################################
# parameters

smNone = LocalSemblanceFilter.Smoothing.NONE
smBoxcar = LocalSemblanceFilter.Smoothing.BOXCAR
smGaussian = LocalSemblanceFilter.Smoothing.GAUSSIAN
smLaplacian = LocalSemblanceFilter.Smoothing.LAPLACIAN

d2U = LocalSemblanceFilter.Direction2.U
d2V = LocalSemblanceFilter.Direction2.V
d2UV = LocalSemblanceFilter.Direction2.UV

plotTitleBarHeight = 23
plotWidthColorBar = 80
plotWidthColorBarTotal = plotWidthColorBar+53
dataClip = 9.5

"""
n1,n2 = 251,357 # Teapot Dome slice vertical
fileName = "tp73.dat"
plotPref = "tpd"
dataDir = "/data/seis/tp/"
dataScale = dataClip/4.0
halfWidth = 2
halfWidth1 = 1*halfWidth
halfWidth2 = 6*halfWidth
sigmaTensor = 8.0
"""

n1,n2 = 500,500 # Atwater channels slice horizontal
fileName = "atwj1s.dat"
plotPref = "atw"
dataDir = "/data/seis/atw/"
dataScale = dataClip/18000.0
halfWidth = 20
halfWidth1 = 1*halfWidth
halfWidth2 = 1*halfWidth
sigmaTensor = 12.0

def setPlotWidthHeight():
  global plotWidth,plotHeight
  plotWidth = 800
  plotHeight = plotWidth*n1/n2
  plotWidth += plotWidthColorBarTotal
  plotHeight += plotTitleBarHeight

setPlotWidthHeight()
plotFontSize = 24
#plotPngDir = "./png/"
plotPngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

small = 0.01
niter = 1000

#############################################################################
# functions

def main(args):
  f = goImage()
  goSmooth2(halfWidth1,d2UV,f)
  #goSemblance2(halfWidth1,halfWidth2,d2U,f)
  #goSemblanceClassic()
  return

def goImage():
  f = readImage(n1,n2,fileName)
  plot(f,-dataClip,dataClip,gray,"f")
  return f

def goSmooth2(hw,d,f):
  t = computeTensors(sigmaTensor,f)
  k1 = 1+4*halfWidth1
  k2 = 1+4*halfWidth1
  f = makeImpulses(n1,n2,k1,k2)
  for sm in [smBoxcar,smGaussian,smLaplacian]:
    lsf = LocalSemblanceFilter(sm,hw,sm,hw)
    g = lsf.smooth1(d,t,f)
    #plot(g,-dataClip,dataClip,gray,"g")
    plot(g,0.0,0.0,gray,"g")

def goSemblance2(hw1,hw2,d,f):
  t = computeTensors(sigmaTensor,f)
  for sm1 in [smBoxcar,smGaussian,smLaplacian]:
    sm2 = sm1
    if d==d2UV:
      sm2 = smNone
    lsf = LocalSemblanceFilter(sm1,hw1,sm2,hw2)
    s = lsf.semblance(d,t,f)
    plot(s,0.0,1.0,gray,"s")

def goSemblanceClassic():
  f = readImage(n1,n2,fileName)
  t = computeTensors(sigmaTensor,f)
  pmax = 10.0
  hw1 = int(sigma2+0.5)
  hw2 = int(sigma1+0.5)
  s = LocalSemblanceFilter.applyForSlopesInWindow(pmax,hw1,hw2,t,f)
  #plot(f,-dataClip,dataClip,gray,"f")
  plot(s,0.0,1.0,gray,"s")

def computeTensors(sigma,f):
  lof = LocalOrientFilter(sigma)
  t = lof.applyForTensors(f)
  return t
 
def makeImpulses(n1,n2,k1,k2):
  m1 = n1/k1
  m2 = n2/k2
  j1 = (n1-(m1-1)*k1)/2
  j2 = (n2-(m2-1)*k2)/2
  f = Array.zerofloat(n1,n2)
  for i2 in range(m2):
    for i1 in range(m1):
      f[j2+i2*k2][j1+i1*k1] = 1.0
  #return smooth(f)
  return f

def smooth(x):
  n1 = len(x[0])
  n2 = len(x)
  t = Array.zerofloat(n1,n2)
  y = Array.zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply0X(x,t)
  rgf.applyX0(t,y)
  return y

def readImage(n1,n2,fileName):
  f = Array.zerofloat(n1,n2)
  ais = ArrayInputStream(dataDir+fileName)
  ais.readFloats(f)
  ais.close()
  return Array.mul(dataScale,f)

def writeImage(f,fileName):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(f)
  aos.close()
 
#############################################################################
# plot

def plot(f,cmin=0.0,cmax=0.0,cmap=ColorMap.GRAY,png=None):
  print "f min =",Array.min(f),"  max =",Array.max(f)
  n1,n2 = len(f[0]),len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(cmap)
  frame(p,png)

def plot2(f,g,cmin=0,cmax=0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.getGray())
  pv.setClips(-10,10)
  pv = p.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin!=cmax:
    pv.setClips(cmin,cmax)
  pv.setColorModel(ColorMap.getJet(0.3))
  frame(p,png)

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  #p = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.LEFT_TOP)
  p.addColorBar()
  p.setColorBarWidthMinimum(plotWidthColorBar)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(plotFontSize)
  frame.setSize(plotWidth,plotHeight)
  frame.setVisible(True)
  if png and plotPngDir:
    frame.paintToPng(200,6,plotPngDir+plotPref+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
