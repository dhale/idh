"""
Cross-validation (leave-one-out)
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
logSet = "deep" # use only deep logs
logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.6
#logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.0,3.0
#logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
#logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
vomit = [5,6,7,9,19,24] # velocity logs to leave out (omit)
smin,smax = -5.5,5.5
sfile = "tpsz" # seismic image
esfile = "tpets" # eigen-tensors scaled by semblances
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,

def main(args):
  goCrossVal()

def goCrossVal():
  for omit in [5]:
    crossval(omit)

def crossval(omit):
  gridWellLogs(omit)
  gridBlendedP(omit)
  gridBlendedQ(omit)

def gridWellLogs(omit=-1):
  print "gridWellLogs:",logType,"without log",omit
  fnull = 0.0
  wlg = WellLogGridder(s1,s2,s3,fnull)
  fl,x1l,x2l,x3l = readLogSamples(logSet,logType)
  for i in range(len(fl)):
    f,x1,x2,x3 = fl[i],x1l[i],x2l[i],x3l[i]
    if i!=omit:
      wlg.insertWellLog(f,x1,x2,x3)
  g = wlg.getGriddedValues()
  writeImage(griddedFile(omit),g)
  
def gridBlendedP(omit=-1):
  print "gridBlendedP:",logType,"without log",omit
  e = getEigenTensors()
  bi = BlendedGridder3(e)
  p = readImage(griddedFile(omit))
  t = bi.gridNearest(0.0,p)
  writeImage(nearestFile(omit),p)
  writeImage(mintimeFile(omit),t)

def gridBlendedQ(omit=-1):
  print "gridBlendedQ:",logType,"without log",omit
  e = getEigenTensors()
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(griddedFile(omit))
  t = readImage(mintimeFile(omit))
  t = clip(0.0,50.0,t)
  q = copy(p)
  #bg.gridBlended(t,p,q)
  writeImage(blendedFile(omit),q)

def getEigenTensors():
  e = readTensors(esfile)
  return e

def griddedFile(omit=-1):
  if omit<0:
    return "tpg"+logType[0]
  elif omit<10:
    return "tpg"+logType[0]+"o0"+str(omit)
  else:
    return "tpg"+logType[0]+"o"+str(omit)

def nearestFile(omit=-1):
  if omit<0:
    return "tpp"+logType[0]
  elif omit<10:
    return "tpp"+logType[0]+"o0"+str(omit)
  else:
    return "tpp"+logType[0]+"o"+str(omit)

def mintimeFile(omit=-1):
  if omit<0:
    return "tpt"+logType[0]
  elif omit<10:
    return "tpt"+logType[0]+"o0"+str(omit)
  else:
    return "tpt"+logType[0]+"o"+str(omit)

def blendedFile(omit=-1):
  if omit<0:
    return "tpq"+logType[0]
  elif omit<10:
    return "tpq"+logType[0]+"o0"+str(omit)
  else:
    return "tpq"+logType[0]+"o"+str(omit)

#############################################################################
#############################################################################
#############################################################################

def gridOne(logType,k):
  p = readImage("ig6/tpp"+logType[0]+"b")
  q = readImage("tpq"+logType[0]+"b")
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

def display(logType,omit=-1):
  x = readImage(seismicFile)
  g = readImage(griddedFile(omit))
  world = World()
  addImage2ToWorld(world,x,g)
  #addLogsToWorld(world,logSet,logType)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
