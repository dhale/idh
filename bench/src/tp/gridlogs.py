"""
Grids well log curves.
"""
from tputils import *

#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
seismicFile = "tpsz"
s1,s2,s3 = getSamplings()
logSet = "deep" # use only deep logs
#smooth = 50 # half-width of Gaussian smoothing filter
smooth = 0 # for no smoothing

def main(args):
  gridAll()
  #displayAll()
  #grid("velocity")
  #display("velocity")
  #display("density")
  #display("porosity")
  #display("gamma")
  
def gridAll():
  for logType in ["velocity","density","porosity","gamma"]:
    grid(logType)

def displayAll():
  for logType in ["velocity","density","porosity","gamma"]:
    display(logType)

def grid(logType):
  fnull = 0.0
  wlg = WellLogGridder(s1,s2,s3,fnull)
  fl,x1l,x2l,x3l = readLogSamples(logSet,logType,smooth)
  for i in range(len(fl)):
    f,x1,x2,x3 = fl[i],x1l[i],x2l[i],x3l[i]
    print logType+": inserting log",i," min =",min(f)," max =",max(f)
    wlg.insertWellLog(f,x1,x2,x3)
  g = wlg.getGriddedValues()
  gfile = griddedFile(logType)
  writeImage(gfile,g)

def gridSimpleBinningDoesNotWorkWell(logType):
  f,x1,x2,x3 = readLogSamplesMerged(logSet,logType)
  print logType+": min =",min(f)," max =",max(f)
  sg = SimpleGridder3(f,x1,x2,x3)
  sg.setNullValue(0.0)
  g = sg.grid(s1,s2,s3)
  #sg.trimEndsOfCurves(g)
  gfile = griddedFile(logType)
  writeImage(gfile,g)

def griddedFile(logType):
    return "tpg"+logType[0]

def display(logType):
  x = readImage(seismicFile)
  g = readImage(griddedFile(logType))
  world = World()
  addImage2ToWorld(world,x,g)
  #addLogsToWorld(world,logSet,logType)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
