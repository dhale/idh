import sys
from math import *

from java.awt import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

from tp import *

def main(args):
  #convertHorizons()
  #plotCurves("gamma",1.0)
  #plotCurves("velocity",40.0)
  plotCurves("density",1.0)
  #plotCurves("porosity",0.0)
  #plotAll()

# Data directories.
resampDir = "/data/seis/tp/resamp/"
horizonsDir = resampDir+"horizons/"
horizonsTextDir = "/data/seis/tp/Transform/horizons/"

# Important files.
imageFile = resampDir+"tp3z.dat"
wdataFile = resampDir+"tp3logs.dat"

# Coordinate sampling.
n1=401; d1=0.004; f1=0.200; s1 = Sampling(n1,d1,f1)
n2=161; d2=0.025; f2=0.000; s2 = Sampling(n2,d2,f2)
n3=357; d3=0.025; f3=0.000; s3 = Sampling(n3,d3,f3)

# Horizon names.
"""
horizonNames = [
  "KF2F2WC",
  "CrowMountainCRMT",
  "FallRiverDKOT",
  "TensleepASand",
  "TensleepBbaseC1Dolo"
  "BasementPC",
]
"""
horizonNames = [
  "CrowMountainCRMT",
]

def zHorizonFile(name):
  return horizonsDir+"tp3z"+name+".dat"

def zHorizonRead(name):
  return Horizon.readBinary(zHorizonFile(name))

def convertHorizons():
  for name in horizonNames:
    h = Horizon.readText(horizonsTextDir+"z"+name+".txt")
    print name," ns =",h.ns," nt=",h.nt
    h.writeBinary(zHorizonFile(name))

def plotCurves(curve,fnull):
  world = World()
  wdata = readWellLogData()
  #image = wdata.rasterizeLogsWith(curve,fnull,s1,s2,s3)
  #wdata.printCounts(image,fnull)
  #printStats(image)
  #ipg = ImagePanelGroup(s1,s2,s3,image)
  #ipg = ImagePanelGroup(Sampling(n1),Sampling(n2),Sampling(n3),image)
  #ipg.setColorModel(ColorMap.JET)
  #world.addChild(ipg)
  image = readImage(imageFile)
  printStats(image)
  ipg = ImagePanelGroup(s1,s2,s3,image)
  world.addChild(ipg)
  addHorizonGroups(world,horizonNames)
  addWellGroups(world,wdata,curve)
  makeFrame(world)

def plotAll():
  world = World()
  wdata = readWellLogData()
  image = readImage(imageFile)
  printStats(image)
  ipg = ImagePanelGroup(s1,s2,s3,image)
  world.addChild(ipg)
  addHorizonGroups(world,horizonNames)
  addWellGroups(world,wdata,"velocity")
  makeFrame(world)

def addHorizonGroups(world,horizonNames):
  for hname in horizonNames:
    h = zHorizonRead(hname)
    tg = makeTriangleGroup(h)
    world.addChild(tg)

def addWellGroups(world,wdata,curve):
  for log in wdata.getLogsWith(curve):
    pg = makePointGroup(log)
    world.addChild(pg)

def makeFrame(world):
  frame = TestFrame(world)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,1.0,3.0)
  view.setScale(2.0)
  view.setAzimuth(10.0)
  frame.setSize(1200,900)
  frame.setVisible(True)
  return frame

def readWellLogData():
  return WellLog.Data.readBinary(wdataFile)

def readImage(fileName):
  x = Array.zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def printStats(image):
  print "image min =",Array.min(image)," max =",Array.max(image)

def makeTriangleGroup(horizon):
  ijk = horizon.getIABC()
  xyz = horizon.getX321()
  tg = TriangleGroup(ijk,xyz)
  tg.setStates(horizonStates)
  return tg;
def makeHorizonStates():
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.CYAN)
  states.add(cs)
  lms = LightModelState()
  lms.setTwoSide(True)
  states.add(lms)
  ms = MaterialState()
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
  ms.setSpecular(Color.WHITE)
  ms.setShininess(100.0)
  states.add(ms)
  return states
horizonStates = makeHorizonStates()

def makePointGroup(log):
  n = log.n
  xyz = Array.zerofloat(3*n)
  Array.copy(n,0,1,log.x3,0,3,xyz)
  Array.copy(n,0,1,log.x2,1,3,xyz)
  Array.copy(n,0,1,log.x1,2,3,xyz)
  #pg = PointGroup(0.020,xyz)
  pg = PointGroup(xyz)
  pg.setStates(logStates)
  return pg
def makeLogStates():
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.YELLOW)
  states.add(cs)
  """
  lms = LightModelState()
  lms.setTwoSide(True)
  states.add(lms)
  ms = MaterialState()
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
  ms.setSpecular(Color.WHITE)
  ms.setShininess(100.0)
  states.add(ms)
  """
  """
  ls = LineState()
  ls.setSmooth(True)
  ls.setWidth(5)
  states.add(ls)
  """
  return states
logStates = makeLogStates()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
