"""
Computes semblance images.
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,

lsf1 = LocalSemblanceFilter(2,2)
lsf2 = LocalSemblanceFilter(2,8)
lsf3 = LocalSemblanceFilter(16,0)

def main(args):
  #semblance1()
  #semblance2()
  #semblance3()
  display(sfile,s1file)
  display(sfile,s2file)
  display(sfile,s3file)

def semblance1():
  e = readTensors(efile)
  s = readImage(sfile)
  s1 = lsf1.semblance(LocalSemblanceFilter.Direction3.W,e,s)
  writeImage(s1file,s1)

def semblance2():
  e = readTensors(efile)
  s = readImage(sfile)
  s2 = lsf2.semblance(LocalSemblanceFilter.Direction3.VW,e,s)
  writeImage(s2file,s2)

def semblance3():
  e = readTensors(efile)
  s = readImage(sfile)
  s3 = lsf3.semblance(LocalSemblanceFilter.Direction3.UVW,e,s)
  writeImage(s3file,s3)

def display(sfile,smfile):
  s = readImage(sfile)
  sm = readImage(smfile)
  print "semblance: min =",min(sm)," max =",max(sm)
  world = World()
  addImage2ToWorld(world,s,sm)
  addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
