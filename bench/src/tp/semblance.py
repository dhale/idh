"""
Computes semblance images.
"""
from tputils import *

setupForDirectory("subz_401_4_600")
s1,s2,s3 = getSamplings()

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,

def main(args):
  #semblance1()
  semblance2()
  display(sfile,s2file)

def semblance1():
  lsf = LocalSemblanceFilter(2,2)
  e = readTensors(efile)
  s = readImage(sfile)
  s1 = lsf.semblance(LocalSemblanceFilter.Direction3.W,e,s)
  writeImage(s1file,s1)

def semblance2():
  lsf = LocalSemblanceFilter(2,8)
  e = readTensors(efile)
  s = readImage(sfile)
  s2 = lsf.semblance(LocalSemblanceFilter.Direction3.VW,e,s)
  writeImage(s2file,s2)

def display(sfile,tfile):
  s = readImage(sfile)
  t = readImage(tfile)
  world = World()
  addImage2ToWorld(world,s,t)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
