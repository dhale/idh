"""
Displays some images.
"""

from mbsutils import *
#setupForSubset("s1")
setupForSubset("s1b")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

def main(args):
  goDisplay()
  #goDisplaySamples()

def goDisplay():
  world = World()
  for image in ["g","gs"]:
    g = readImage(image)
    ipg = addImageToWorld(world,g)
    ipg.setClips(-1.0,1.0)
  frame = makeFrame(world)
  frame.setSize(1200,900)

def goDisplaySamples():
  g = readImage("gs")
  sf = SimpleFrame()
  ipg = sf.addImagePanels(g)
  ipg.setClips(-1.0,1.0)
  sf.orbitView.setScale(3.0)
  sf.setSize(1200,900)

#############################################################################
run(main)


