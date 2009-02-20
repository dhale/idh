import sys
from math import *
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

#True = 1
#False = 0

#############################################################################
# parameters

fontSize = 24
width = 300
height = 800
widthColorBar = 80
dataDir = "/data/seis/oz/"
#pngDir = "./png"
pngDir = None

#############################################################################
# functions

def main(args):
  for i in range(1,6):
    display(i)
  return

def display(index):
  s1,s2,f = readData(index)
  f = tpow(2.0,s1,f)
  f = gpow(0.5,f)
  plot(s1,s2,f)

def fileName(index):
  if index<10:
    name = "oz0"+str(index)+".F"
  else:
    name = "oz"+str(index)+".F"
  return dataDir+name

def samplings(index):
  if index==1:
    n1,d1,f1 = 1275,0.004,0.004
    n2,d2,f2 = 53,0.100584,-2.615184
  elif index==2:
    n1,d1,f1 = 1025,0.004,0.004
    n2,d2,f2 = 127,0.030480,-2.834640
  elif index==3:
    n1,d1,f1 = 1500,0.004,0.004
    n2,d2,f2 = 24,0.103632,0.103632
  elif index==4:
    n1,d1,f1 = 1275,0.004,0.004
    n2,d2,f2 = 52,0.100,-2.550
  elif index==5:
    n1,d1,f1 = 3000,0.002,0.002
    n2,d2,f2 = 51,0.100,-2.500
  else:
    return None,None
  return Sampling(n1,d1,f1),Sampling(n2,d2,f2)

def gpow(p,f):
  return Array.mul(Array.sgn(f),Array.pow(Array.abs(f),p))

def tpow(p,s,f):
  nt,dt,ft = s.count,s.delta,s.first
  nx = len(f)
  g = Array.copy(f)
  for ix in range(nx):
    Array.mul(g[ix],Array.pow(Array.rampfloat(ft,dt,nt),p),g[ix])
  return g

def readData(index):
  s1,s2 = samplings(index)
  n1,n2 = s1.count,s2.count
  ais = ArrayInputStream(fileName(index),ByteOrder.BIG_ENDIAN)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return s1,s2,f

#############################################################################
# plot

def plot(s1,s2,f,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  pv = p.addPixels(s1,s2,f)
  pv.setPercentiles(1.0,99.0)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(p,png)

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  #p.addColorBar()
  #p.setColorBarWidthMinimum(widthColorBar)
  p.setHLabel("offset (km)");
  p.setVLabel("time (s)");
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
