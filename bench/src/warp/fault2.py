#############################################################################
# Dynamic time warping for fault images
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
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

#from fault import *
#from fault.Util import *

#############################################################################

pngDir = "./png/fault"
#pngDir = None

s1,s2 = None,None
label1,label2 = None,None

def main(args):
  #goFaultImages()
  goFaultShifts()

def goFaultShifts():
  global smoothShifts,smoothSigma1,smoothSigma2; 
  maxlag = 20
  uclips = (0,8)
  f,g = makeFaultImages()
  fclips = (-3.0,3.0)
  dw = DynamicWarping(0,maxlag)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.REFLECT)
  dw.setStrainMax(0.25,0.25)
  dw.setErrorSmoothing(2)
  dw.setShiftSmoothing(1.0,2.0)
  u = dw.findShifts(f,g)
  h = dw.applyShifts(u,g)
  #v = copy(u)
  #LocalShiftFinder(maxlag,maxlag/2).find1(0,maxlag,f,g,v)
  uclips = (0.0,7*4)
  plot = plotp
  plot(f,fclips,label="Amplitude",png="faultf")
  plot(g,fclips,label="Amplitude",png="faultg")
  plot(mul(4,u),uclips,label="Shift (ms)",png="faultu")
  #plot(mul(4,v),uclips,label="Shift (ms)",png="faultv")
  plot(h,fclips,label="Amplitude",png="faulth")

def goFaultImages():
  f,g = makeFaultImages()
  clips = (-3,3)
  plot(f,clips,label="Amplitude",png="faultf")
  plot(g,clips,label="Amplitude",png="faultg")

""" 
sampling of 3D image subset used to extract fault images
j1,j2,j3 = 240, 50,100
n1,n2,n3 = 222,221,220
d1,d2,d3 = 0.004,0.025,0.025
f1,f2,f3 = 0.964,0.000,0.000
"""
def makeFaultImages():
  n1,n2 = 222,220
  m1,m2 =  90,216
  j1,j2 = 130,2
  d1,d2 = 0.004,0.025
  f1,f2 = 1.484,2.550
  global s1,s2,label1,label2
  s1 = Sampling(m1,d1,f1)
  s2 = Sampling(m2,d2,f2)
  label1 = "Time (s)"
  label2 = "Distance along fault strike (km)"
  f = readImage("/data/seis/f3d/faults/s1/gfm.dat",n1,n2)
  g = readImage("/data/seis/f3d/faults/s1/gfp.dat",n1,n2)
  f = copy(m1,m2,j1,j2,f)
  g = copy(m1,m2,j1,j2,g)
  return f,g
def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x
 
#############################################################################
# plot

def plotp(f,clips=None,title=None,label=None,png=None):
  n1,n2 = len(f[0]),len(f)
  width,height = 550,325
  cbwm = 60
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(Color.WHITE)
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  pv = sp.addPixels(s1,s2,f)
  if png[-1]!="u":
    gv = sp.addGrid("H-")
    gv.setColor(Color.BLACK)
  if clips:
    pv.setClips(clips[0],clips[1])
  if title:
    sp.setTitle(title)
  if label:
    sp.addColorBar(label)
  sp.setVLabel(label1)
  sp.setHLabel(label2)
  sp.setHInterval(1.0)
  sp.setVInterval(0.1)
  sp.setFontSizeForPrint(8,240)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,3.33333,pngDir+"/p"+png+".png")

def plots(f,clips=None,title=None,label=None,png=None):
  n1,n2 = len(f[0]),len(f)
  width,height = 600,400
  cbwm = 70
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(Color.WHITE)
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  pv = sp.addPixels(s1,s2,f)
  if png[-1]=="u":
    pv.setColorModel(ColorMap.JET)
  if png[-1]!="u":
    gv = sp.addGrid("H-.")
    gv.setColor(Color.YELLOW)
  if clips:
    pv.setClips(clips[0],clips[1])
  if title:
    sp.setTitle(title)
  if label:
    sp.addColorBar(label)
  sp.setVLabel(label1)
  sp.setHLabel(label2)
  sp.setHInterval(1.0)
  sp.setVInterval(0.1)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(768,2.0,pngDir+"/s"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
