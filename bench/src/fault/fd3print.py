#############################################################################
# Print figures for fault displacements in 3D images

import sys
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fault import *

gmin,gmax,gint,gatt,glab = -5.5,5.5,2.0,None,"Log amplitude"
t1min,t1max,t1int,t1att = 0.0,15.0,2.0,None
t1lab = "Vertical component of throw (ms)"
flmin,flmax,flint,flatt,fllab = 0.5,1.0,0.05,None,"Fault likelihood"
background = Color(255,255,255) # pure white for print
pngDir = "png/"
#pngDir = None

def main(args):
  #goFigures3d()
  goFigures3f()

def goFigures3f():
  global k1,k2,k3,pngPre
  for subset in ["a","b"]:
    setupForF3dSubset(subset)
    g = readImage("g")
    gs = readImage("gs")
    h = readImage("h")
    t1 = readImage("t1"); t1 = mul(4.0,t1)
    fl = readImage("fl")
    flt = readImage("flt")
    for kkk in kkks:
      k1,k2,k3,pngPre = kkk
      #plot3f(g,None,gmin,gmax,gint,gatt,glab,"g")
      #plot3f(gs,None,gmin,gmax,gint,gatt,glab,"gs")
      plot3f(g,t1,t1min,t1max,t1int,t1att,t1lab,"gt1")
      plot3f(h,t1,t1min,t1max,t1int,t1att,t1lab,"ht1")

def goFigures3d():
  global k1,k2,k3
  k1,k2,k3 = 366,15,96 # good for 3D displays
  s = rimage("tpsz"); plot3d1(s,smin,smax,sint,slog,slab)

def setupForF3dSubset(subset):
  global s1,s2,s3
  global n1,n2,n3
  global dataDir,dataPre,kkks
  dataDir = "/data/seis/f3d/faults/"
  n1,n2,n3 = 462,951,591
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.004,0.000,0.000
  if subset=="a": # deeper normal faults
    dataPre = "a"
    j1,j2,j3 = 240+130,50,100
    m1,m2,m3 = 90,221,220
    kkks = [
      [26,38,119,"t26"]]
      #[26,38,119,"t26"], # ??? 26(39), 38, 119
      #[42,69,162,"t42"]] # ??? 42,     69, 162(124)
  elif subset=="b": # shallower conical faults
    dataPre = "b"
    #j1,j2,j3 = 240,50,100
    #m1,m2,m3 = 120,221,220
    j1,j2,j3 = 240+15,50,100
    m1,m2,m3 = 90,221,220
    kkks = [
      [55,47,146,"t55"],
      [59,52,142,"t59"]]
  elif subset=="s1": # conical fault above normal faults
    dataPre = "s1"
    j1,j2,j3 = 240,50,100
    m1,m2,m3 = 222,221,220
  n1,n2,n3 = m1,m2,m3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

def readImage(fileName):
  if dataPre=="b":
    f = zerofloat(n1+30,n2,n3)
  else:
    f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+dataPre+fileName+".dat")
  ais.readFloats(f)
  ais.close()
  if dataPre=="b":
    f = copy(90,n2,n3,15,0,0,f)
  return f

def plot3f(g,c=None,cmin=0,cmax=0,cint=None,cmap=None,clab=None,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1, 50)
  if c:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar(clab)
    if cint:
      cb.setInterval(cint)
  else:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar("Log amplitude")
    #cb.setInterval(2.0)
  pp.setInterval1(0.1)
  #pp.setInterval2(1.0)
  #pp.setInterval3(1.0)
  if c:
    pv12 = PixelsView(s1,s2,slice12(k3,c))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,slice13(k2,c))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,slice23(k1,c))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    for pv in [pv12,pv13,pv23]:
      #pv.setColorModel(bwrNotch(1.0))
      pv.setColorModel(jetRamp(1.0))
      #pv.setColorModel(jetFill(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(120)
  pf.setFontSizeForPrint(8,240)
  #pf.setFontSizeForPrint(8,504)
  #pf.setSize(1008,800)
  pf.setSize(1008,672)
  pf.setVisible(True)
  if png and pngDir:
    #pf.paintToPng(720,3.3,pngDir+pngPre+dataPre+png+".png")
    pf.paintToPng(720,7.0,pngDir+pngPre+dataPre+png+".png")

def plot3d1(s,cmin,cmax,cint,logType,logLabel,horizons=[]):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.viewCanvas.setBackground(background)
  frame.orbitView.setAzimuth(-65.0)
  if logLabel:
    cbar = addColorBar3d(frame,logLabel,cint)
    ipg.addColorMapListener(cbar)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,cbar)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)

def plot3d2(s,g,cmin,cmax,cint,logType,logLabel,horizons=[],cval=0):
  world = World()
  ipg = addImage2ToWorld(world,s,g)
  ipg.setClips1(smin,smax)
  ipg.setClips2(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  if cval:
    addContourToWorld(world,g,cval)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.orbitView.setAzimuth(-65.0)
  frame.viewCanvas.setBackground(background)
  cbar = addColorBar3d(frame,logLabel,cint)
  ipg.addColorMap2Listener(cbar)

def addColorBar3d(frame,label,cint=None):
  cbar = ColorBar(label)
  cbar.setFont(Font("Arial",Font.PLAIN,48)) # ~ 8*1460/240 for one-column
  cbar.setBackground(background)
  if cint:
    cbar.setInterval(cint)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,1.0/256,256))
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
    """
    if i<96:
      a[i] = 1.0
    elif i<128:
      a[i] = alpha*(128.0-i)/32.0
    elif i<160:
      a[i] = alpha*(i-127.0)/32.0
    else:
      a[i] = 1.0
    """
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
