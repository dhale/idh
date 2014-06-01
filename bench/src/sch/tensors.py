"""
Computes tensors for local smoothing filters.
"""

from schutils import *
setupForSubset("s2")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gfile = "g" # seismic image
dfile = "d" # smoothing tensors

def main(args):
  makeTensors()
  display()

def makeTensors():
  g = readImage(gfile)
  lof = LocalOrientFilter(8.0)
  d = lof.applyForTensors(g)
  d.invertStructure(0.0,2.0,4.0)
  writeTensors(dfile,d)

def display():
  g = readImage(gfile)
  d = readTensors(dfile)
  world = World()
  ipg = addImageToWorld(world,g)
  ipg.setClips(-1.0,1.0)
  addTensorsInImage(ipg.getImagePanel(Axis.X),d,20)
  addTensorsInImage(ipg.getImagePanel(Axis.Y),d,20)
  addTensorsInImage(ipg.getImagePanel(Axis.Z),d,20)
  frame = makeFrame(world)
  frame.orbitView.setAzimuth(-65.0)
  frame.setSize(1460,980)

#############################################################################
run(main)
