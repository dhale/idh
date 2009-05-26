import sys
from org.python.util import PythonObjectInputStream
from math import *
from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *

from lss import *

#############################################################################
# parameters

smNone = LocalSemblanceFilter.Smoothing.NONE
smBoxcar = LocalSemblanceFilter.Smoothing.BOXCAR
smGaussian = LocalSemblanceFilter.Smoothing.GAUSSIAN
smLaplacian = LocalSemblanceFilter.Smoothing.LAPLACIAN

d3U = LocalSemblanceFilter.Direction3.U
d3V = LocalSemblanceFilter.Direction3.V
d3W = LocalSemblanceFilter.Direction3.W
d3UV = LocalSemblanceFilter.Direction3.UV
d3UW = LocalSemblanceFilter.Direction3.UW
d3VW = LocalSemblanceFilter.Direction3.VW

plotTitleBarHeight = 23
plotWidthColorBar = 80
plotWidthColorBarTotal = plotWidthColorBar+53
fClip = 9

def setTpd(): # Teapot Dome slice vertical
  global n1,n2,n3,k1,k2,k3,dataDir,dataPref,fScale
  global halfWidth,halfWidth1,halfWidth2,sigmaTensor
  global plotWidth,plotHeight
  global plotWidth2,plotHeight2
  n1,n2,n3 = 251,161,357
  #k1,k2,k3 = 183,62,138 # good slices
  k1,k2,k3 = 183,73,138 # good slices
  dataDir = "/data/seis/tp/"
  dataPref = "tp3"
  fScale = fClip/4.0
  #halfWidth = 4
  halfWidth = 2
  halfWidth1 = 1*halfWidth
  halfWidth2 = 4*halfWidth
  sigmaTensor = 8.0
  plotWidth = 1040
  plotHeight = 745
  plotWidth = 1035
  plotHeight = 670

def setAtw(): # Atwater channels slice horizontal
  global n1,n2,n3,k1,k2,k3,dataDir,dataPref,fScale
  global halfWidth,halfWidth1,halfWidth2,sigmaTensor
  global plotWidth,plotHeight
  global plotWidth2,plotHeight2
  n1,n2,n3 = 129,500,500
  k1,k2,k3 =  40,260,190 # good slices?
  dataDir = "/data/seis/atw/"
  dataPref = "atw"
  fScale = fClip/15000.0
  halfWidth = 2
  halfWidth1 = 1*halfWidth
  halfWidth2 = 1*halfWidth
  #halfWidth2 = 4*halfWidth
  sigmaTensor = 12.0
  plotWidth = 1040
  plotHeight = 610
  plotWidth2 = 780
  plotHeight2 = 670

def setPlotWidthHeight():
  global plotWidth,plotHeight
  plotWidth = 900
  plotHeight = 700
  #plotHeight = plotWidth*n1/n2
  plotWidth += plotWidthColorBarTotal
  plotHeight += plotTitleBarHeight

plotFontSize = 32
plotPngDir = "./png/sos/"
#plotPngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

#############################################################################
# functions

def main(args):
  setAtw(); goAll()
  #setTpd(); goAll()
  return

def goAll():
  #goImage()
  #goTensors()
  #goSemblanceVW()
  #goSemblanceW()
  #goSemblanceClassic()
  #goSmoothGSW()
  #goPlot()
  goPlot2()

def goPlot2():
  global plotWidth,plotHeight
  plotWidth = plotWidth2
  plotHeight = plotHeight2
  for name in ["svwl2_2","swl2_2"]:
    s = readImage(n1,n2,n3,name)
    s = slice1(40,s)
    p = panel()
    ps = p.addPixels(s)
    ps.setClips(0.0,1.0)
    frame(p,name)

def goImage():
  f = readImage(n1,n2,n3,"f",fScale)
  plot3([f])

def goTensors():
  f = readImage(n1,n2,n3,"f",fScale)
  lof = LocalOrientFilter(sigmaTensor)
  d = lof.applyForTensors(f)
  writeTensors(d,"st12")

def goSemblanceVW():
  hw1 = halfWidth1
  hw2 = halfWidth2
  f = readImage(n1,n2,n3,"f",fScale)
  t = readTensors()
  #for sm1 in [smBoxcar,smGaussian,smLaplacian]:
  for sm1 in [smLaplacian]:
    sm2 = sm1
    lsf = LocalSemblanceFilter(sm1,hw1,sm2,hw2)
    s = lsf.semblance(d3VW,t,f)
    name = "svw"+smstr(sm1)+str(hw1)+"_"+str(hw2)
    writeImage(s,name)
    print "s min =",Array.min(s),"max =",Array.max(s)
    plot3([f,s])

def goSemblanceW():
  hw1 = halfWidth1
  hw2 = halfWidth2
  f = readImage(n1,n2,n3,"f",fScale)
  t = readTensors()
  for sm1 in [smLaplacian]:
    sm2 = sm1
    lsf = LocalSemblanceFilter(sm1,hw1,sm2,hw2)
    s = lsf.semblance(d3W,t,f)
    name = "sw"+smstr(sm1)+str(hw1)+"_"+str(hw2)
    writeImage(s,name)
    print "s min =",Array.min(s),"max =",Array.max(s)
    plot3([f,s])

def goSmoothGSW():
  hw = 20
  hw1 = halfWidth1
  hw2 = halfWidth2
  sm = smLaplacian
  sname = "sw"+smstr(sm)+str(hw1)+"_"+str(hw2)
  gname = "g"+sname
  t = readTensors()
  f = readImage(n1,n2,n3,"f",fScale)
  s = readImage(n1,n2,n3,sname,1.0)
  eu = Array.copy(s)
  ev = Array.copy(s)
  ew = Array.copy(s)
  t.getEigenvalues(eu,ev,ew)
  e1 = Array.div(Array.sub(ev,ew),eu)
  s = Array.mul(s,e1)
  #s = Array.mul(s,s)
  t.setEigenvalues(0.0,0.0,1.0)
  lsf = LocalSmoothingFilter()
  c = hw*(hw+1)/6.0
  g = Array.copy(f)
  lsf.apply(t,c,s,f,g)
  writeImage(s,gname)
  plot3([f,g])

def goSemblanceClassic():
  pmax = 10.0
  hw1 = halfWidth1
  hw2 = halfWidth2
  f = readImage(n1,n2,n3,"f",fScale)
  t = readTensors()
  s = LocalSemblanceFilter.semblanceForSlopes(pmax,hw1,hw2,t,f)
  name = "ssc"+str(hw1)+"_"+str(hw2)
  writeImage(s,name)
  print "s min =",Array.min(s),"max =",Array.max(s)
  plot3([f,s])

def goPlot():
  #for name in ["f","svwl2_8","svwb2_8","svwg2_8","ssc2_8"]:
  #for name in ["f","svwl2_2","swl2_2"]:
  for name in ["f","svwl2_8"]:
    if name=="f":
      f = readImage(n1,n2,n3,name,fScale)
      plotp3(k1,k2,k3,f,-fClip,fClip,name)
    else:
      s = readImage(n1,n2,n3,name)
      plotp3(k1,k2,k3,s,0.0,1.0,name)

def slice1(k1,f):
  return Array.reshape(n2,n3,Array.flatten(Array.copy(1,n2,n3,40,0,0,f)))

def computeTensors(sigma,f):
  f = readImage(n1,n2,n3,"f")
  lof = LocalOrientFilter(sigma)
  d = lof.applyForTensors(f)
  return d

def readImage(n1,n2,n3,fileName,scale=1.0):
  f = Array.zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+dataPref+fileName+".dat")
  ais.readFloats(f)
  ais.close()
  return Array.mul(scale,f)

def writeImage(f,fileName):
  aos = ArrayOutputStream(dataDir+dataPref+fileName+".dat")
  aos.writeFloats(f)
  aos.close()
 
def readTensors():
  tensorsFile = "st"+str(int(sigmaTensor))
  fis = FileInputStream(dataDir+dataPref+tensorsFile+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors
 
def writeTensors(tensors,tensorsFile):
  fos = FileOutputStream(dataDir+dataPref+tensorsFile+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

def smstr(sm):
  if sm==smBoxcar:
    return "b"
  elif sm==smGaussian:
    return "g"
  else:
    return "l"
 
#############################################################################
# plot

def plot3(flist):
  world = World()
  for f in flist:
    ipg = ImagePanelGroup(f)
    ipg.setPercentiles(1.0,99.0)
    #ipg.setClips(0.0,1.0)
    world.addChild(ipg)
    clipMin = ipg.getClipMin()
    clipMax = ipg.getClipMax()
    print "clip min =",clipMin,"max =",clipMax
  frame = TestFrame(world)
  frame.setVisible(True)

def plotp3(k1,k2,k3,f,cmin,cmax,png=None):
  print "plotp3: min =",Array.min(f)," max =",Array.max(f)
  panel = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    Sampling(n1),Sampling(n2),Sampling(n3),f)
  panel.addColorBar()
  panel.setSlice23(k1)
  panel.setSlice13(k2)
  panel.setSlice12(k3)
  panel.setInterval1(100.0)
  panel.setInterval2(100.0)
  panel.setInterval3(100.0)
  panel.setLabel1("time")
  panel.setLabel2("inline")
  panel.setLabel3("crossline")
  panel.setColorBarWidthMinimum(plotWidthColorBar)
  panel.setClips(cmin,cmax)
  panel.setLineColor(Color.BLACK)
  return frame(panel,png)

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  p.addColorBar()
  p.setColorBarWidthMinimum(plotWidthColorBar)
  #p.setHInterval(100)
  #p.setVInterval(100)
  #p.setHLabel(hlabel)
  #p.setVLabel(vlabel)
  #p.setHLabel(" ")
  #p.setVLabel(" ")
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(plotFontSize)
  frame.setSize(plotWidth,plotHeight)
  frame.setVisible(True)
  if png and plotPngDir:
    frame.paintToPng(200,6,plotPngDir+dataPref+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
