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

from jss import *

#############################################################################
# parameters

dataDir = "/data/seis/jss/"
#pngDir = "png/jss/"
pngDir = None

# samplings
n1 = 1400; d1 = 0.004; f1 = 0.8
n2 = 1600; d2 = 0.010; f2 = 2.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)

# Which pair of images? (1 = conv, 2 = simsrc)
index = 1
prefix = "s"+str(index)

# scale factors for different pairs of images (make them easier to compare)
#scales = [1.00,6.00,2.73]

# local correlation filter
lmax = 3
lmin = -lmax
lcfSigma1 = 12.0
lcfSigma2 = 12.0 # make this 36.0 to reduce vertical stripes in uxv
lcfWhiten = 6.0
lcfType = LocalCorrelationFilter.Type.SIMPLE
lcfWindow = LocalCorrelationFilter.Window.GAUSSIAN

# Clip for shifts
uclip = 7.0

#############################################################################
# functions

def main(args):
  doPair(1)
  #doAllPairs()
  return

def doAllPairs():
  for i in [1,2,3]:
    doPair(i)

def doPair(i):
  global index,prefix
  index = i
  prefix = "s"+str(index)
  f,g = doImages()
  fw,gw = doWhiten(f,g)
  uza1 = measureShiftsV(f,g)
  uxa,uza = measureShifts(fw,gw)
  if i==1:
    uxv,uzv,ux,uz = deriveShifts(uxa,uza)
    writeImage(uz,prefix+"uz")

def doImages():
  f = readImage(prefix+"f")
  g = readImage(prefix+"g")
  plot(f,20.0,"f")
  plot(g,20.0,"g")
  return f,g

def doScale(f,g):
  fs = Array.mul(scales[index],f)
  gs = Array.mul(scales[index],g)
  return fs,gs

def doWhiten(f,g):
  fw = whiten(f)
  gw = whiten(g)
  plot(fw,2.0,"fw")
  plot(gw,2.0,"gw")
  return fw,gw

def measureShifts(f,g):
  u1 = Array.zerofloat(n1,n2)
  u2 = Array.zerofloat(n1,n2)
  du = Array.zerofloat(n1,n2)
  h = Array.copy(g)
  sf = LocalShiftFinder(lcfSigma1,lcfSigma2)
  for iter in range(3):
    sf.find1(lmin,lmax,f,h,du)
    sf.shift1(du,u1,u2,h)
    sf.find2(lmin,lmax,f,h,du)
    sf.shift2(du,u1,u2,h)
  u1 = Array.mul(1000.0*d1,u1)
  u2 = Array.mul(1000.0*d2,u2)
  plotu(u2,uclip,"uxa")
  plotu(u1,uclip,"uza")
  return u2,u1

def measureShiftsV(f,g):
  u1 = Array.zerofloat(n1,n2)
  sf = LocalShiftFinder(lcfSigma1,lcfSigma2)
  sf.find1(lmin,lmax,f,g,u1)
  u1 = Array.mul(1000.0*d1,u1)
  plotu(u1,uclip,"uza1")
  return u1

def deriveShifts(uxa,uza):
  dz,dx = 1000.0*d1,1000.0*d2 # convert d1 and d2 to m
  nz,nx = n1,n2
  r = 5.0
  uxv = Array.zerofloat(nz,nx)
  uzv = Array.zerofloat(nz,nx)
  ux = Array.zerofloat(nz,nx)
  uz = Array.zerofloat(nz,nx)
  Shifts.shifts(r,dx,dz,uxa,uza,uxv,uzv,ux,uz)
  plotu(uxv,uclip,"uxv")
  plotu(uzv,uclip,"uzv")
  plotu(ux,0.5*uclip,"ux")
  plotu(uz,0.5*uclip,"uz")
  return uxv,uzv,ux,uz

def readImage(fileName):
  fileName = dataDir+fileName+".dat"
  ais = ArrayInputStream(fileName)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(f,fileName):
  outfile = dataDir+fileName+".dat"
  aos = ArrayOutputStream(outfile)
  aos.writeFloats(f)
  aos.close()

def whiten(f):
  fw = Array.zerofloat(n1,n2)
  sf = LocalShiftFinder(lcfWhiten)
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
  #print "plot clip min =",clipMin," max =",clipMax
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
  frame.setSize(1100,770)
  #frame.setSize(2200,1540)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(300,6,pngDir+prefix+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

