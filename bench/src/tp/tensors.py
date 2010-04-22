"""
Computes structure eigen-tensors.
"""
from tputils import *

#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

sfile = "tpsz" # seismic image
mfile = "tpmz" # image mask
efile = "tpet" # eigen-tensors
esfile = "tpets" # eigen-tensors scaled by semblance
slfile = "tps1" # semblance w,uv
spfile = "tps2" # semblance vw,u
sifile = "tps3" # semblance uvw,

def main(args):
  #makeStructureTensors()
  display()

def makeStructureTensors():
  sigma = 8.0
  s = readImage(sfile)
  m = readImage(mfile)
  mask = ZeroMask(m)
  lof = LocalOrientFilter(sigma)
  e = lof.applyForTensors(s)
  mask.apply((1.0,0.0,0.0,1.0,0.0,1.0),e)
  writeTensors(efile,e)

def display():
  s = readImage(sfile)
  #et = readTensors(efile)
  et = readTensors(esfile)
  #eu = zerofloat(n1,n2,n3)
  #ev = zerofloat(n1,n2,n3)
  #ew = zerofloat(n1,n2,n3)
  #et.getEigenvalues(eu,ev,ew)
  #eu = mul(eu,eu)
  #ev = mul(ev,ev)
  #ev = mul(ev,ev)
  #ew = mul(ew,ew)
  #et.setEigenvalues(eu,ev,ew)
  world = World()
  ipg = addImageToWorld(world,s)
  addTensorsInImage(ipg.getImagePanel(Axis.X),et,20)
  addTensorsInImage(ipg.getImagePanel(Axis.Y),et,20)
  addTensorsInImage(ipg.getImagePanel(Axis.Z),et,20)
  makeFrame(world)

def displayOld():
  s = readImage(sfile)
  e = readTensors(efile)
  eu = zerofloat(n1,n2,n3)
  ev = zerofloat(n1,n2,n3)
  ew = zerofloat(n1,n2,n3)
  e.getEigenvalues(eu,ev,ew)
  #e0 = div(ew/eu)         #  isotropy: e0 = (ew   )/eu
  #e1 = div(sub(ev,ew)/eu) # linearity: e1 = (ev-ew)/eu
  e2 = div(sub(eu,ev),eu) # planarity: e2 = (eu-ev)/eu
  world = World()
  addImage2ToWorld(world,s,e2)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
