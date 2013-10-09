"""
Displays Teapot Dome data.
"""
from tputils import *

#setupForSubset("subt_251_4_500")
setupForSubset("subz_401_4_600")

def main(args):
  #displaySlices()
  displaySubset()

def displaySubset():
  world = World()
  x = readImage("tpsz")
  addImageToWorld(world,x)
  #addAllHorizonsToWorld(world)
  addHorizonToWorld(world,"CrowMountainCRMT")
  addHorizonToWorld(world,"TensleepASand")
  addHorizonToWorld(world,"TensleepBbaseC1Dolo")
  addLogsToWorld(world,"d","v",cmin=2.5,cmax=4.5)
  makeFrame(world)

def displaySlices():
  displaySlice("tpsz",ColorMap.GRAY)
  displaySlice("tpgv",ColorMap.JET)
  displaySlice("tpgd",ColorMap.JET)
  displaySlice("tpgp",ColorMap.JET)
  displaySlice("tpgg",ColorMap.JET)

def displaySlice(name,cmap):
  x = readSlice3("s3_84/"+name)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  pv = sp.addPixels(x)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmap!=None:
    pv.setColorModel(cmap)

#############################################################################
run(main)
