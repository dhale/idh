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

#############################################################################
# parameters

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
sigmaTensor = 8.0
sigmaSmooth = 2.0
sigma1 = 1.0*sigmaSmooth
sigma2 = 3.0*sigmaSmooth
parallel = True
"""

n1,n2 = 500,500 # Atwater channels slice horizontal
fileName = "atwj1s.dat"
plotPref = "atw"
dataDir = "/data/seis/atw/"
dataScale = dataClip/18000.0
sigmaTensor = 12.0
sigmaSmooth = 6.0
sigma1 = 1.0*sigmaSmooth
sigma2 = 1.0*sigmaSmooth
parallel = True

if parallel:
  plotPref += "p"
else:
  plotPref += "x"

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
  goSmooth(parallel)
  goSemblance(parallel)
  return

def goSemblance(parallel):
  f = readImage(n1,n2,fileName)
  t = computeTensors(sigmaTensor,f)
  s = applySemblance(parallel,sigma1,sigma2,t,f)
  plot(f,-dataClip,dataClip,gray,"f")
  plot(s,0.0,1.0,gray,"s")

def goSmooth(parallel):
  ks = int(5*sigmaSmooth)
  k1,k2 = ks,ks
  f = readImage(n1,n2,fileName)
  i = makeImpulses(n1,n2,k1,k2)
  t = computeTensors(sigmaTensor,f)
  g = applySmooth(parallel,sigmaSmooth,t,f)
  h = applySmooth(parallel,sigmaSmooth,t,i)
  plot(f,-dataClip,dataClip,gray,"f")
  plot(g,-dataClip,dataClip,gray,"g")
  plot(h,-Array.max(h),Array.max(h),gray,"h")

def computeTensors(sigma,f):
  lof = LocalOrientFilter(sigma)
  t = lof.applyForTensors(f)
  return t

def applySemblance(parallel,sigma1,sigma2,t,f):
  c1 = 0.5*sigma1*sigma1
  c2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilter(small,niter)
  g = Array.copy(f)
  if parallel:
    t.setEigenvalues(0.0,1.0)
  else:
    t.setEigenvalues(1.0,0.0)
  lsf.apply(t,c1,f,g)
  ff = Array.mul(f,f)
  gg = Array.mul(g,g)
  lsf.apply(t,c1,Array.copy(ff),ff)
  sff = Array.copy(ff)
  sgg = Array.copy(gg)
  if parallel:
    t.setEigenvalues(1.0,0.0)
  else:
    t.setEigenvalues(0.0,1.0)
  lsf.apply(t,c2,ff,sff)
  lsf.apply(t,c2,gg,sgg)
  s = Array.div(sgg,sff)
  print "s min =",Array.min(s),"  max =",Array.max(s)
  s = Array.clip(0.0,1.0,s)
  return s

def applySmooth(parallel,sigma,t,f):
  g = Array.copy(f)
  c = 0.5*sigma*sigma
  lsf = LocalSmoothingFilter(small,niter)
  if parallel:
    t.setEigenvalues(0.01,1.00)
  else:
    t.setEigenvalues(1.00,0.01)
  lsf.apply(t,c,f,g)
  return g
 
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

def makePlanarSemblance():
  f = readImage(n1,n2,n3,fFile)
  g = Array.zerofloat(n1,n2,n3)
  scale1 = 0.5*sigma1*sigma1
  scale2 = 0.5*sigma2*sigma2
  lsf1 = LocalSmoothingFilterX(scale1,small,niter)
  lsf2 = LocalSmoothingFilterX(scale2,small,niter)
  t = readTensors(tensorsFile)
  t.setEigenvalues(0.0,1.0,1.0)
  lsf2.apply(t,f,g)
  writeImage(g,gFile)
  ff = Array.mul(f,f)
  gg = Array.mul(g,g)
  lsf2.apply(t,Array.copy(ff),ff)
  writeImage(ff,ffFile)
  writeImage(gg,ggFile)
  sff = Array.zerofloat(n1,n2,n3)
  sgg = Array.zerofloat(n1,n2,n3)
  t.setEigenvalues(1.0,0.0,0.0)
  lsf1.apply(t,ff,sff)
  lsf1.apply(t,gg,sgg)
  writeImage(sff,sffFile)
  writeImage(sgg,sggFile)
  ps = Array.div(sgg,sff)
  ps = Array.clip(0.0,1.0,ps)
  writeImage(ps,psFile)

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

def smooth(x):
  n1 = len(x[0])
  n2 = len(x)
  t = Array.zerofloat(n1,n2)
  y = Array.zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply0X(x,t)
  rgf.applyX0(t,y)
  return y
 
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
