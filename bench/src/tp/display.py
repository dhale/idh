"""
Displays Teapot Dome data.
"""
from tputils import *

#setupForSubset("subt_251_4_500")
setupForSubset("subz_401_4_600")

def main(args):
  displaySlice("tpsz",ColorMap.GRAY)
  displaySlice("tpgv",ColorMap.JET)
  displaySlice("tpgd",ColorMap.JET)
  displaySlice("tpgp",ColorMap.JET)
  displaySlice("tpgg",ColorMap.JET)
  #displaySubset()

def displaySlice(name,cmap):
  x = readSlice3("s3_84/"+name)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  pv = sp.addPixels(x)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmap!=None:
    pv.setColorModel(cmap)

def displaySubset():
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
