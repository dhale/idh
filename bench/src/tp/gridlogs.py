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

def main(args):
  #gridAll()
  #displayAll()
  omit = 0
  #grid("velocity",omit=omit)
  #display("velocity",omit=omit)
  for k in [5,6,7,9,19,24]: #range(29):
    gridOne("velocity",k)
  #grid("velocity")
  #display("velocity")
  #display("density")
  #display("porosity")
  #display("gamma")

def gridOne(logType,k):
  p = readImage("ig6/tpp"+logType[0]+"b")
  q = readImage("ig6/tpq"+logType[0]+"b")
  fl,x1l,x2l,x3l = readLogSamples(logSet,logType)
  wlg = WellLogGridder(s1,s2,s3,0.0)
  fk,x1k,x2k,x3k = fl[k],x1l[k],x2l[k],x3l[k]
  fg,x1g,x2g,x3g = wlg.getGriddedSamples(fk,x1k,x2k,x3k)
  fp = wlg.getGriddedValues(x1k,x2k,x3k,p)
  fq = wlg.getGriddedValues(x1k,x2k,x3k,q)
  print x2g[0],x3g[0]
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(400,900)
  sp.setTitle("Log "+str(k))
  sp.setHLimits(2.0,7.0)
  sp.setVLimits(0.6,2.0)
  sp.setHLabel("Velocity (km/s)")
  sp.setVLabel("Depth (km)")
  pv = sp.addPoints(x1k,fk)
  pv.setLineColor(Color.BLACK)
  pv = sp.addPoints(x1g,fg)
  pv.setLineColor(Color.LIGHT_GRAY)
  pv.setLineWidth(5.0)
  pv = sp.addPoints(x1g,fp)
  pv.setLineColor(Color.RED)
  pv.setLineWidth(5.0)
  pv = sp.addPoints(x1g,fq)
  pv.setLineColor(Color.BLUE)
  pv.setLineWidth(5.0)
  
def gridAll():
  for logType in ["velocity","density","porosity","gamma"]:
    grid(logType)

def displayAll():
  for logType in ["velocity","density","porosity","gamma"]:
    display(logType)

def grid(logType,omit=-1):
  fnull = 0.0
  wlg = WellLogGridder(s1,s2,s3,fnull)
  fl,x1l,x2l,x3l = readLogSamples(logSet,logType)
  for i in range(len(fl)):
    f,x1,x2,x3 = fl[i],x1l[i],x2l[i],x3l[i]
    if i==omit:
      print logType+": omitting log",i," min =",min(f)," max =",max(f)
    else:
      print logType+": inserting log",i," min =",min(f)," max =",max(f)
      wlg.insertWellLog(f,x1,x2,x3)
  g = wlg.getGriddedValues()
  gfile = griddedFile(logType,omit)
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

def griddedFile(logType,omit=-1):
  if omit<0:
    return "tpg"+logType[0]
  elif omit<10:
    return "tpg"+logType[0]+"o0"+str(omit)
  else:
    return "tpg"+logType[0]+"o"+str(omit)

def display(logType,omit=-1):
  x = readImage(seismicFile)
  g = readImage(griddedFile(logType,omit))
  world = World()
  addImage2ToWorld(world,x,g)
  #addLogsToWorld(world,logSet,logType)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
