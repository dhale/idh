"""
Computes a mask for samples that are zero or near zero.
"""
from tputils import *

#setupForSubset("subt_251_4_500")
#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

sfile = "tpsz" # seismic image
mfile = "tpmz" # mask image
#sfile = "tpst" # seismic image
#mfile = "tpmt" # mask image

def main(args):
  s = readImage(sfile)
  mask = ZeroMask(0.1,10.0,1.0,1.0,s)
  m = mask.getAsFloats()
  writeImage(mfile,m)
  display()

def display():
  s = readImage(sfile)
  m = readImage(mfile)
  world = World()
  addImage2ToWorld(world,s,m)
  makeFrame(world)

#############################################################################
run(main)
