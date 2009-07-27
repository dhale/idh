"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

setupForDirectory("subz_401_4_600")
s1,s2,s3 = getSamplings()
seismicFile = "tpsz"
logSet = "deep"

def main(args):
  #logType = "density"
  logType = "velocity"
  #logType = "porosity"
  #logType = "gamma"
  #gridSibson(logType)
  #gridNearest(logType)
  gridBlendedP(logType)
  display(logType)

def simplegFile(logType):
  return "tpg"+logType[0]
def interpgFile(logType):
  return "tpi"+logType[0]
def interppFile(logType):
  return "tpp"+logType[0]
def interptFile(logType):
  return "tpt"+logType[0]

def gridBlendedP(logType):
  gfile = simplegFile(logType)
  p = readImage(gfile)
  bi = BlendedGridder3()
  t = bi.gridNearest(0.0,p)
  pfile = interppFile(logType)
  tfile = interptFile(logType)
  writeImage(pfile,p)
  writeImage(tfile,t)

def getScatteredSamples(logType):
  g = readImage(simplegFile(logType))
  f,x1,x2,x3 = SimpleGridder3.getGriddedSamples(0.0,s1,s2,s3,g)
  return f,x1,x2,x3

def gridNearest(logType):
  f,x1,x2,x3 = getScatteredSamples(logType)
  print "got scattered samples: n =",len(f)
  ni = NearestGridder3(f,x1,x2,x3)
  print "constructed nearest gridder"
  g = ni.grid(s1,s2,s3)
  print "gridding complete: min =",min(g)," max =",max(g)
  gfile = interpgFile(logType)
  writeImage(gfile,g)

def gridSibson(logType):
  f,x1,x2,x3 = getScatteredSamples(logType)
  print "got scattered samples: n =",len(f)
  si = SibsonGridder3(f,x1,x2,x3)
  print "constructed Sibson gridder"
  g = si.grid(s1,s2,s3)
  print "gridding complete: min =",min(g)," max =",max(g)
  gfile = interpgFile(logType)
  writeImage(gfile,g)

def display(logType):
  x = readImage(seismicFile)
  g = readImage(interpgFile(logType))
  world = World()
  addImage2ToWorld(world,x,g)
  addLogsToWorld(world,logSet,logType)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
