"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
logSet = "d" # deep logs only
logType = "v" # velocity
method = "b" # blended

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors (structure tensors)
esfile = "tpets" # eigen-tensors scaled by semblances
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,
gfile = "tpg"+logType # simple gridding with null for unknown samples
pfile = "tpp"+logType+method # values of nearest known samples
qfile = "tpq"+logType+method # output of blended gridder
tfile = "tpt"+logType+method # times to nearest known samples

def main(args):
  #gridBlendedP()
  #gridBlendedQ()
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  #t = readImage(tfile); print "t min =",min(t)," max =",max(t)
  q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #display1(s)
  display(s,p,3.0,5.0)
  #display(s,t,0.0,100.0)
  display(s,q,3.0,5.0)

def gridBlendedP():
  e = getEigenTensors(0.01)
  bi = BlendedGridder3(e)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def gridBlendedQ():
  e = getEigenTensors(0.0)
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(pfile)
  t = readImage(tfile)
  t = clip(0.0,10.0,t)
  q = copy(p)
  bg.gridBlended(t,p,q)
  #writeImage(qfile,q)

def getEigenTensors(eps):
  e = readTensors(efile)
  s1 = readImage(s1file); print "s1 min =",min(s1)," max =",max(s1)
  s2 = readImage(s2file); print "s2 min =",min(s2)," max =",max(s2)
  s3 = readImage(s3file); print "s3 min =",min(s3)," max =",max(s3)
  s1 = clip(eps,1.0,s1)
  s2 = clip(eps,1.0,s2)
  s3 = clip(eps,1.0,s3)
  e.setEigenvalues(s3,s2,s1)
  #writeTensors(esfile,e)
  return e

def display1(s):
  world = World()
  ipg = addImageToWorld(world,s)
  addLogsToWorld(world,logSet,logType)
  #addHorizonToWorld(world,"CrowMountainCRMT")
  #addHorizonToWorld(world,"TensleepASand")
  addHorizonToWorld(world,"TensleepBbaseC1Dolo")
  makeFrame(world)

def display(s,g,cmin,cmax):
  world = World()
  ipg = addImage2ToWorld(world,s,g)
  ipg.setClips2(cmin,cmax)
  addLogsToWorld(world,logSet,logType)
  #addHorizonToWorld(world,"CrowMountainCRMT")
  #addHorizonToWorld(world,"TensleepASand")
  addHorizonToWorld(world,"TensleepBbaseC1Dolo")
  makeFrame(world)

#############################################################################

def dumpTensors():
  e = getEigenTensors()
  d = zerofloat(6)
  for i3 in range(14,17):
    for i2 in range(110,113):
      for i1 in range(249,252):
        e.getTensor(i1,i2,i3,d)
        dump(d)

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

def testSpd():
  e = getEigenTensors()
  lsf = LocalSmoothingFilter(0.01,10000)
  t = readImage(tfile)
  t = clip(0.0,10.0,t)
  s = mul(t,t)
  lsf.testSpd(e,0.5,s)

def gridBlended2():
  sdir = "s3_84/"
  s = readSlice3(sdir+"tpsz")
  p = readSlice3(sdir+"tpgd"); cmin,cmax = 2000,2800
  #display2(s)
  display2(s,p,cmin,cmax)
  q = copy(p)
  lof = LocalOrientFilter(8)
  e = lof.applyForTensors(s)
  e.setEigenvalues(0.01,1.0)
  bi = BlendedGridder2(e)
  t = bi.gridNearest(0.0,p)
  t = clip(0,100,t)
  e.setEigenvalues(0.0,1.0)
  bi.setTensors(e)
  bi.setSmoothness(0.5)
  bi.gridBlended(t,p,q)
  #display2(s,p,cmin,cmax)
  display2(s,q,cmin,cmax)

def display2(s,g=None,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.getPlotPanel().setColorBarWidthMinimum(80)
  pv = sp.addPixels(s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g!=None:
    pv = sp.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.getJet(0.3))
    if cmin!=cmax:
      pv.setClips(cmin,cmax)

#############################################################################
run(main)
