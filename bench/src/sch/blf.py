"""
Structure-oriented bilateral smoothing filter
"""

from schutils import *
setupForSubset("s2")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gfile = "g" # input seismic image
gsfile = "gs" # output smoothed image
dfile = "d" # smoothing tensors

def main(args):
  #goSmoothingFilter()
  goBilateralFilter()
  display()
 
def goSmoothingFilter():
  d = readTensors(dfile)
  g = readImage(gfile)
  sigmaS = 4.0
  c = 0.5*sigmaS*sigmaS
  lsf = LocalSmoothingFilter()
  gs = zerofloat(n1,n2,n3)
  lsf.apply(d,c,g,gs)
  writeImage(gsfile,gs)
 
def goBilateralFilter():
  d = readTensors(dfile)
  g = readImage(gfile)
  sigmaS = 4.0 # sufficient to attenuate noise
  #sigmaR = 0.3 # leaves too much noise speckle
  sigmaR = 0.7 # smoothing ok?
  #sigmaR = 1.0 # smooths perhaps a bit too much
  #sigmaR = computeSigmaR(g) # ~ 0.165
  print "sigmaS =",sigmaS," sigmaR =",sigmaR
  bf = BilateralFilter(sigmaS,sigmaR)
  gs = zerofloat(n1,n2,n3)
  bf.apply(d,g,gs)
  writeImage(gsfile,gs)

def computeSigmaR(g):
  return 0.5*(Quantiler.estimate(0.75,g)-Quantiler.estimate(0.25,g))

def display():
  g = readImage(gfile)
  gs = readImage(gsfile)
  world = World()
  ipg = addImageToWorld(world,g)
  ipg.setClips(-1.0,1.0)
  ipg = addImageToWorld(world,gs)
  ipg.setClips(-1.0,1.0)
  frame = makeFrame(world)
  #frame.orbitView.setAzimuth(-65.0)
  frame.setSize(1460,980)

#############################################################################
run(main)
