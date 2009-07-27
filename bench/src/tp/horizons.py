"""
Reads and reformats horizons from Teapot Dome.
Author: Dave Hale, Colorado School of Mines
Version: 2009.06.07
"""
from imports import *

#############################################################################
def main(args):
  #setGlobals("t"); makeHorizons() # make binary files for time horizons
  #setGlobals("z"); makeHorizons() # make binary files for depth horizons
  setGlobals("t"); viewHorizons() # view time horizons
  #setGlobals("z"); viewHorizons() # view depth horizons

# Horizon names with colors, ordered by increasing time/depth.
horizonColors = {
  "KF2F2WC":Color.RED,
  "FallRiverDKOT":Color.GREEN,
  "CrowMountainCRMT":Color.BLUE,
  "TensleepASand":Color.CYAN,
  "TensleepBbaseC1Dolo":Color.MAGENTA,
  "BasementPC":Color.YELLOW
}
horizonNames = horizonColors.keys()

def setGlobals(w):
  global what,tssHorizonDir,csmHorizonDir,csmSeismicImage
  global s1,s2,s3,time,depth
  what = w
  tpDir = "/data/seis/tp/"
  tssHorizonDir = tpDir+"tss/horizons/"
  s2 = Sampling(357,0.025,0.000)
  s3 = Sampling(161,0.025,0.000)
  if what=="t":
    csmHorizonDir = tpDir+"csm/horizont/"
    csmSeismicImage = tpDir+"csm/seismict/tpst.dat"
    s1 = Sampling(1501,0.002,0.000)
    time = True; depth = False
  else:
    csmHorizonDir = tpDir+"csm/horizonz/"
    csmSeismicImage = tpDir+"csm/seismicz/tpsz.dat"
    s1 = Sampling(2762,0.002,0.000)
    time = False; depth = True

def makeHorizons():
  for name in horizonNames:
    h = Horizon.readText(tssHorizonDir+what+name+".txt",time)
    h.clip(s2,s3)
    print name," ns =",h.ns," nt=",h.nt
    h.writeBinary(horizonFile(name))

def viewHorizons():
  x = readImage(csmSeismicImage)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  world = World()
  world.addChild(ipg)
  addHorizonGroups(world)
  frame = makeFrame(world)

def horizonFile(name):
  return csmHorizonDir+"tph"+what+name+".dat"

def horizonRead(name):
  return Horizon.readBinary(horizonFile(name))

def addHorizonGroups(world):
  for hname in horizonNames:
    h = horizonRead(hname)
    c = horizonColors[hname]
    tg = makeTriangleGroup(h,c)
    world.addChild(tg)

def makeFrame(world):
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,1.0,3.0)
  view.setScale(2.0)
  view.setAzimuth(-50.0)
  frame.setSize(1200,900)
  frame.setVisible(True)
  return frame

def readImage(filename):
  n1,n2,n3 = s1.count,s2.count,s3.count
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(filename)
  ais.readFloats(x)
  ais.close()
  return x

def makeTriangleGroup(horizon,color):
  ijk = horizon.getIABC()
  xyz = horizon.getX321()
  tg = TriangleGroup(ijk,xyz)
  tg.setColor(color)
  return tg;

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
