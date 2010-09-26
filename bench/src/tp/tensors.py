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
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,

def main(args):
  #makeTensors()
  scaleTensors(0.001)
  #display()

def makeTensors():
  sigma = 8.0
  s = readImage(sfile)
  m = readImage(mfile)
  mask = ZeroMask(m)
  lof = LocalOrientFilter(sigma)
  e = lof.applyForTensors(s)
  mask.apply((1.0,0.0,0.0,1.0,0.0,1.0),e)
  writeTensors(efile,e)

def scaleTensors(eps):
  e = readTensors(efile)
  s1 = readImage(s1file); print "s1 min =",min(s1)," max =",max(s1)
  s2 = readImage(s2file); print "s2 min =",min(s2)," max =",max(s2)
  s3 = readImage(s3file); print "s3 min =",min(s3)," max =",max(s3)
  pow(s2,4.0,s2)
  fill(eps,s3)
  s1 = clip(eps,1.0,s1)
  s2 = clip(eps,1.0,s2)
  s3 = clip(eps,1.0,s3)
  e.setEigenvalues(s3,s2,s1)
  writeTensors(esfile,e)

def display():
  k1,k2,k3 = 366,15,96
  s = readImage(sfile)
  et = readTensors(esfile)
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(-5.5,5.5)
  ipg.setSlices(k1,k2,k3)
  addTensorsInImage(ipg.getImagePanel(Axis.X),et,20)
  addTensorsInImage(ipg.getImagePanel(Axis.Y),et,20)
  addTensorsInImage(ipg.getImagePanel(Axis.Z),et,20)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.orbitView.setAzimuth(-65.0)
  background = Color(254,254,254)
  frame.viewCanvas.setBackground(background)

#############################################################################
run(main)
