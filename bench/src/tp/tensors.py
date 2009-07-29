"""
Computes structure eigen-tensors.
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

sigma = 8.0
sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors

def main(args):
  makeStructureTensors()
  #display()

def makeStructureTensors():
  s = readImage(sfile)
  lof = LocalOrientFilter(sigma)
  e = lof.applyForTensors(s)
  writeTensors(efile,e)

def display():
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
