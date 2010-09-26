"""
Computes semblance images.
"""
from tputils import *

#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()

sfile = "tpsz" # seismic image
mfile = "tpmz" # mask image
efile = "tpet" # eigen-tensors
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,

lsf1 = LocalSemblanceFilter(2,2)
lsf2 = LocalSemblanceFilter(2,8)
lsf3 = LocalSemblanceFilter(16,0)
#lsf1 = LocalSemblanceFilter(8,2)
#lsf2 = LocalSemblanceFilter(8,2)
#lsf3 = LocalSemblanceFilter(8,0)
mask = ZeroMask(readImage(mfile))

def main(args):
  semblance1()
  semblance2()
  semblance3()
  maskSemblance()
  #display(sfile,s1file)
  #display(sfile,s2file)
  #display(sfile,s3file)

def maskSemblance():
  for smfile in [s1file,s2file,s3file]:
    sm = readImage(smfile)
    if smfile==s3file:
      mask.apply(0.0001,sm)
    else:
      mask.apply(1.00,sm)
    writeImage(smfile,sm)

def semblance1():
  e = readTensors(efile)
  s = readImage(sfile)
  s1 = lsf1.semblance(LocalSemblanceFilter.Direction3.W,e,s)
  mask.apply(1.00,s1)
  writeImage(s1file,s1)

def semblance2():
  e = readTensors(efile)
  s = readImage(sfile)
  s2 = lsf2.semblance(LocalSemblanceFilter.Direction3.VW,e,s)
  mask.apply(1.00,s2)
  writeImage(s2file,s2)

def semblance3():
  e = readTensors(efile)
  s = readImage(sfile)
  s3 = lsf3.semblance(LocalSemblanceFilter.Direction3.UVW,e,s)
  mask.apply(0.01,s3)
  writeImage(s3file,s3)

def display(sfile,smfile):
  s = readImage(sfile)
  sm = readImage(smfile)
  print "semblance: min =",min(sm)," max =",max(sm)
  world = World()
  ipg = addImage2ToWorld(world,s,sm)
  ipg.setClips2(0,1)
  #addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
