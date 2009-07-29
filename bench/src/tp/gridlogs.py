"""
Grids well log curves.
"""
from tputils import *

setupForSubset("subz_401_4_600")
seismicFile = "tpsz"
s1,s2,s3 = getSamplings()
logSet = "deep" # use only deep logs

def main(args):
  #gridAll()
  display("velocity")

def gridAll():
  for logType in ["velocity","density","porosity","gamma"]:
    grid(logType)

def grid(logType):
  f,x1,x2,x3 = readLogSamplesMerged(logSet,logType)
  print logType+": min =",min(f)," max =",max(f)
  sg = SimpleGridder3(f,x1,x2,x3)
  sg.setNullValue(0.0)
  g = sg.grid(s1,s2,s3)
  gfile = "tpg"+logType[0]
  writeImage(gfile,g)

def griddedFile(logType):
  return "tpg"+logType[0]

def display(logType):
  x = readImage(seismicFile)
  g = readImage(griddedFile(logType))
  world = World()
  addImage2ToWorld(world,x,g)
  addLogsToWorld(world,logSet,logType)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
