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

#############################################################################
# parameters

dataDir = "/data/seis/jss/"
pngDir = "png"
#pngDir = None

# samplings
n1 =  850; d1 = 0.0060; f1 = 0.9
n2 = 1280; d2 = 0.0125; f2 = 2.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)

# Which pair of images? (0 = gold, 1 = conv, 2 = simsrc)
index = 1
prefix = "s"+str(index)

# scale factors for different pairs of images (make them easier to compare)
scales = [1.00,6.00,2.73]

# local correlation filter
lmax = 3
lmin = -lmax
lcfSigma = 12.0
lcfType = LocalCorrelationFilter.Type.SIMPLE
lcfWindow = LocalCorrelationFilter.Window.GAUSSIAN
lcf = LocalCorrelationFilter(lcfType,lcfWindow,lcfSigma)

#############################################################################
# functions

def main(args):
  #doPair(2)
  doAllPairs()
  return

def doAllPairs():
  for i in range(3):
    doPair(i)

def doPair(i):
  global index,prefix
  index = i
  prefix = "s"+str(index)
  f,g = doImages()
  fw,gw = doWhiten(f,g)
  u1,u2 = doShifts(fw,gw)

def doImages():
  f = readImage(prefix+"f")
  g = readImage(prefix+"g")
  f,g = doScale(f,g)
  plot(f,50.0,"f")
  plot(g,50.0,"g")
  return f,g

def doScale(f,g):
  fs = Array.mul(scales[index],f)
  gs = Array.mul(scales[index],g)
  return fs,gs

def doWhiten(f,g):
  fw = whiten(f)
  gw = whiten(g)
  plot(fw,5.0,"fw")
  plot(gw,5.0,"gw")
  return fw,gw

def doShifts(f,g):
  u1 = Array.zerofloat(n1,n2)
  u2 = Array.zerofloat(n1,n2)
  du = Array.zerofloat(n1,n2)
  h = Array.copy(g)
  sf = LocalShiftFinder(lcfSigma)
  for iter in range(3):
    sf.find1(lmin,lmax,f,h,du)
    sf.shift1(du,u1,u2,h)
    sf.find2(lmin,lmax,f,h,du)
    sf.shift2(du,u1,u2,h)
  u1 = Array.mul(1000.0*d1,u1)
  u2 = Array.mul(1000.0*d2,u2)
  plotu(u1,5.0,"u1")
  plotu(u2,5.0,"u2")
  return u1,u2

def readImage(fileName):
  fileName = dataDir+fileName+".dat"
  ais = ArrayInputStream(fileName)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def whiten(f):
  fw = Array.zerofloat(n1,n2)
  sf = LocalShiftFinder(lcfSigma)
  sf.whiten(1.0,f,fw)
  return fw

def doScales():
  f0 = readImage("s0f")
  f1 = readImage("s1f")
  f2 = readImage("s2f")
  a01 = Array.sum(Array.mul(f0,f1))/Array.sum(Array.mul(f1,f1))
  a02 = Array.sum(Array.mul(f0,f2))/Array.sum(Array.mul(f2,f2))
  print "a01 =",a01," a02 =",a02

#############################################################################
# plot

def plot(f,clip=0.0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  pv = p.addPixels(s1,s2,f)
  if clip!=0.0:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  clipMin,clipMax = pv.getClipMin(),pv.getClipMax()
  print "plot clip min =",clipMin," max =",clipMax
  frame(p,png)

def plotu(u,clip=0.0,png=None):
  p = panel()
  pv = p.addPixels(s1,s2,u)
  pv.setColorModel(ColorMap.RED_WHITE_BLUE)
  if clip!=0.0:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(1.0,99.0)
  frame(p,png)

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  p.setHLabel("distance (km)")
  p.setVLabel("depth (km)")
  p.addColorBar()
  p.setColorBarWidthMinimum(90)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(36)
  frame.setSize(1100,800)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(300,6,pngDir+"/"+prefix+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

