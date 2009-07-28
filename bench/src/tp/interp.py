"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

setupForDirectory("subz_401_4_600")
s1,s2,s3 = getSamplings()
logSet = "d" # deep logs only
logType = "v" # velocity
method = "b" # blended

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,
gfile = "tpg"+logType # simple gridding with null for unknown samples
pfile = "tpp"+logType+method # values of nearest known samples
qfile = "tpq"+logType+method # output of blended gridder
tfile = "tpt"+logType+method # times to nearest known samples

def main(args):
  #gridSibson() # very slow!
  #gridNearest()
  gridBlendedP()
  #gridBlendedQ()
  s = readImage(sfile)
  p = readImage(pfile)
  display(s,p)

def gridBlendedP():
  e = getEigenTensors()
  bi = BlendedGridder3(e)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def gridBlendedQ():
  bi = BlendedGridder3()
  p = readImage(pfile)
  t = readImage(tfile)
  t = clip(0.0,100.0,t)
  q = copy(p)
  bi.gridBlended(t,p,q)
  writeImage(qfile,q)

def getEigenTensors():
  e = readTensors(efile)
  s1 = readImage(s1file)
  s2 = readImage(s2file)
  s3 = readImage(s3file)
  e.setEigenvalues(s3,s2,s1)
  return e

def getScatteredSamples():
  g = readImage(gfile)
  f,x1,x2,x3 = SimpleGridder3.getGriddedSamples(0.0,s1,s2,s3,g)
  return f,x1,x2,x3

def gridNearest():
  f,x1,x2,x3 = getScatteredSamples()
  print "got scattered samples: n =",len(f)
  ni = NearestGridder3(f,x1,x2,x3)
  print "constructed nearest gridder"
  g = ni.grid(s1,s2,s3)
  print "gridding complete: min =",min(g)," max =",max(g)
  writeImage(gfile,g)

def gridSibson():
  f,x1,x2,x3 = getScatteredSamples()
  print "got scattered samples: n =",len(f)
  si = SibsonGridder3(f,x1,x2,x3)
  print "constructed Sibson gridder"
  g = si.grid(s1,s2,s3)
  print "gridding complete: min =",min(g)," max =",max(g)
  writeImage(gfile,g)

def display(s,g):
  world = World()
  addImage2ToWorld(world,s,g)
  addLogsToWorld(world,logSet,logType)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
