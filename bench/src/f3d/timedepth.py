"""
Converts F3 from time to depth.
NOTE: the time-depth volume produced by this script has artifacts.
These are caused by errors in the file F3_Stacking_Vel.txt.
Author: Dave Hale, Colorado School of Mines
Version: 2012.12.30
"""
from f3utils import *
seismicDir = getSeismicDir()
s1,s2,s3 = getSamplings()

#############################################################################
def main(args):
  getz()

def getz():
  td = TimeDepth(seismicDir+"odt/F3_Stacking_Vel.txt")
  z = td.getZ(s1,s2,s3)
  return
  f = readImage("f3d")
  world = World()
  addImage2ToWorld(world,f,z,cmin=-1.0,cmax=1.0)
  makeFrame(world)

#############################################################################
run(main)
