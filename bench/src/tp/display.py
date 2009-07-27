"""
Displays Teapot Dome data.
"""
from tputils import *

#setupForDirectory("subt_251_4_500")
setupForDirectory("subz_401_4_600")

def main(args):
  world = World()
  x = readImage("tpsz")
  addImageToWorld(world,x)
  #addAllHorizonsToWorld(world)
  #addHorizonToWorld(world,"CrowMountainCRMT")
  addHorizonToWorld(world,"TensleepASand")
  addLogsToWorld(world,"d","p")
  makeFrame(world)

#############################################################################
run(main)
