"""
Displays 3D images using the Mines Java Toolkit
"""
from common import *

#############################################################################
# time,trace,line sampling and other info in global variables

n1,n2,n3 = None,None,None
d1,d2,d3 = None,None,None
f1,f2,f3 = None,None,None
s1,s2,s3 = None,None,None
clip = None
dataDir = "/data/seis/" # base directory for all images
def setupForImage(name):
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3,clip,dataDir
  n1,n2,n3 = 1001,1001,321 # numbers of samples
  d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
  f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
  dataDir += "nor/csm/dat/"
  if name is "full":
    clip = 1.0
  elif name is "zoft":
    pass
  elif name is "tofz":
    pass
  else:
    print "do not recognize",name
    sys.exit(0)
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

#############################################################################
# Add displays to this simple function, as necessary.
def goDisplay(what):
  if what=="full":
    setupForImage("full")
    g = readImage("norne2006full")
    showOne(g,clip=clip)
  elif what=="zoft":
    setupForImage("zoft")
    g = readImage("zoft")
    showOne(g,cmin=0,cmax=4)
  elif what=="tofz":
    setupForImage("tofz")
    g = readImage("tofz")
    showOne(g,cmin=0,cmax=4)
  else:
    print "Sorry, do not know how to display ",what
    sys.exit(0)

#############################################################################
# Utilities

def main(args):
  goDisplay(args[1])

def readImage(name):
  fileName = dataDir+name+".dat"
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+name+".dat")
  ais.readFloats(f)
  ais.close()
  return f

def showOne(g,cmin=None,cmax=None,clip=None,cmap=None):
  sf = SimpleFrame()
  ipg = sf.addImagePanels(s1,s2,s3,g)
  if cmin and cmax:
    ipg.setClips(cmin,cmax)
  elif clip:
    ipg.setClips(-clip,clip)
  else:
    ipg.setPercentiles(2,98)
    print "clips: min =",ipg.clipMin," max =",ipg.clipMax
  if cmap:
    ipg.setColorModel(cmap)
  sf.orbitView.setScale(4.0)
  sf.orbitView.setAxesScale(1.0,1.0,0.5)
  sf.setSize(1000,800)

def showTwo(g1,g2,cmin=None,cmax=None,clip=None,cmap=None):
  sf = SimpleFrame()
  for g in [g1,g2]:
    ipg = sf.addImagePanels(s1,s2,s3,g)
    if cmin and cmax:
      ipg.setClips(cmin,cmax)
    elif clip:
      ipg.setClips(-clip,clip)
    if cmap:
      ipg.setColorModel(cmap)
  sf.orbitView.setScale(4.0)
  sf.orbitView.setAxesScale(1.0,1.0,0.5)
  sf.setSize(1000,800)

def showTwoOverlay(g1,g2,
  cmin1=None,cmax1=None,clip1=None,cmap1=None,
  cmin2=None,cmax2=None,clip2=None,cmap2=None):
  ipg = ImagePanelGroup2(s1,s2,s3,g1,g2)
  if cmin1 and cmax1:
    ipg.setClips1(cmin1,cmax1)
  elif clip1:
    ipg.setClips1(-clip1,clip1)
  if cmap1:
    ipg.setColorModel1(cmap1)
  if cmin2 and cmax2:
    ipg.setClips2(cmin2,cmax2)
  elif clip2:
    ipg.setClips2(-clip2,clip2)
  if cmap2:
    ipg.setColorModel2(cmap2)
  sf = SimpleFrame()
  sf.world.addChild(ipg)
  sf.orbitView.setScale(4.0)
  sf.orbitView.setAxesScale(1.0,1.0,0.5)
  sf.setSize(1000,800)

def showThree(g1,g2,g3,cmin=None,cmax=None,clip=None,cmap=None):
  sf = SimpleFrame()
  for g in [g1,g2,g3]:
    ipg = sf.addImagePanels(s1,s2,s3,g)
    if cmin and cmax:
      ipg.setClips(cmin,cmax)
    elif clip:
      ipg.setClips(-clip,clip)
    if cmap:
      ipg.setColorModel(cmap)
  sf.orbitView.setScale(4.0)
  sf.orbitView.setAxesScale(1.0,1.0,0.5)
  sf.setSize(1000,800)

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
    """
    if i<96:
      a[i] = 1.0
    elif i<128:
      a[i] = alpha*(128.0-i)/32.0
    elif i<160:
      a[i] = alpha*(i-127.0)/32.0
    else:
      a[i] = 1.0
    """
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

#############################################################################
run(main)
