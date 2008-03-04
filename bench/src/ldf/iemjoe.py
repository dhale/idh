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

dataDir = "/data/seis/joe/"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

paintBar = None

n1 = 200
n2 = 200
aniso = 10
small = 0.001
niter = 1000
lof = LocalOrientFilter(8)
lof.setGradientSmoothing(1)

def readImage(fileName):
  ais = ArrayInputStream(dataDir+fileName)
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
  ds = el
  return v1,v2,ds

iem = None
pvv = None
vnull = 0.0
v = Array.zerofloat(n1,n2)
x = readImage("xz174.dat")
v1,v2,ds = getV(x)
x = readImage("x174.dat")

#############################################################################
# functions

def main(args):
  global iem,v,pvv
  p = panel()
  #p.addColorBar()
  p.setVLabel("z (samples)")
  p.setHLabel("x (samples)")
  pv = p.addPixels(x)
  pv = p.addPixels(x)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setClips(-0.005,0.005);
  pv.setColorModel(ColorMap.GRAY)
  pv = p.addPixels(v)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.getJet(0.3))
  pvv = pv
  f = frame(p)
  mm = f.getModeManager()
  iem = ImageEditMode(mm,pv,vnull,v)
  makeToolBarAndMenus(f)
  f.setVisible(True)
  return

def interp000():
  interpolate(0)
def interp010():
  interpolate(10)
def interp100():
  interpolate(100)
def interpolate(aniso):
  print "interpolate"
  f = Array.zerobyte(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if v[i2][i1]!=vnull:
        f[i2][i1] = 1
  ldt = LocalDiffusionTensors2(1.0/aniso,1.0,None,None,v1)
  small = 0.001
  niter = 1000
  lif = LocalInterpolationFilterIc(small,niter)
  lif.apply(ldt,f,v)
  pvv.set(v)
def interpolateOld(aniso):
  print "interpolate"
  lif = LocalInterpolationFilter(aniso,small,niter)
  f = Array.zerobyte(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if v[i2][i1]!=vnull:
        f[i2][i1] = 1
  ds = None
  es = None
  lif.applyLinear(ds,es,v1,f,v)
  pvv.set(v)
 
#############################################################################
# User interface

fontSize = 24
width = 800
height = 750
widthColorBar = 80

class CallbackAction(AbstractAction):
  def __init__(self,name,desc,callback):
    AbstractAction.__init__(self)
    self.putValue(self.NAME,name)
    self.putValue(self.SHORT_DESCRIPTION,desc)
    self.callback = callback
  def actionPerformed(self,event):
    self.callback()

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  p.setColorBarWidthMinimum(widthColorBar)
  return p

def frame(panel):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  #frame.setVisible(True)
  return frame

def makeToolBarAndMenus(frame):
    modeMenu = JMenu("Mode")
    modeMenu.setMnemonic('M')
    tzm = frame.getTileZoomMode()
    tzmItem = ModeMenuItem(tzm)
    modeMenu.add(tzmItem)
    iemItem = ModeMenuItem(iem)
    modeMenu.add(iemItem)

    menuBar = JMenuBar()
    menuBar.add(modeMenu)

    toolBar = JToolBar(SwingConstants.VERTICAL)
    toolBar.setRollover(True)
    tzmButton = ModeToggleButton(tzm)
    toolBar.add(tzmButton)
    iemButton = ModeToggleButton(iem)
    iemButton.setText("Edit");
    toolBar.add(iemButton)

    intAction = CallbackAction("Int000","Interpolate",interp000)
    intButton = JButton(intAction)
    toolBar.add(intButton)
    intAction = CallbackAction("Int010","Interpolate",interp010)
    intButton = JButton(intAction)
    toolBar.add(intButton)
    intAction = CallbackAction("Int100","Interpolate",interp100)
    intButton = JButton(intAction)
    toolBar.add(intButton)

    frame.add(toolBar,BorderLayout.WEST)
    frame.setJMenuBar(menuBar)


 
#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
