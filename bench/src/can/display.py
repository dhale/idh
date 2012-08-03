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
dataDir = "../" # image data assumed to be in the parent directory
def setupForImage(name):
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3,clip,dataDir
  if name[0:3]=="atw":
    if name[3:5]=="00":
      n1,n2,n3 = 851,1932,604 # numbers of samples
    else:
      n1,n2,n3 = 851,1932,603 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 10000.0
    dataDir += "atw/"
  elif name=="can":
    n1,n2,n3 = 1501,4001,3001 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0
    dataDir += "can/"
  elif name=="f3d":
    n1,n2,n3 = 462,951,591 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0
    dataDir += "f3d/"
  elif name=="gbc":
    n1,n2,n3 = 2000,150,145 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0
    dataDir += "gbc/dat/"
  elif name[0:3]=="mbs":
    if name[3:8]=="Small":
      n1,n2,n3 = 501,448,422 # numbers of samples
    else:
      n1,n2,n3 = 501,769,560 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0
    dataDir += "mbs/dat/"
  elif name=="nor":
    n1,n2,n3 = 1001,1001,321 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0
    dataDir += "nor/dat/"
  elif name=="pen":
    n1,n2,n3 = 1501,480,456 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 3000.0
    dataDir += "pen/"
  elif name=="pnz":
    n1,n2,n3 = 751,1001,1001 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0e5
    dataDir += "pnz/dat/"
  elif name=="seam":#NOT IMPLEMENTED
    n1,n2,n3 = 2000,150,145 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 1.0
    dataDir += "gbc/dat/"
  elif name=="tpd":
    #n1,n2,n3 = 2762,188,345 # numbers of samples
    n1,n2,n3 = 1501,357,161 # numbers of samples
    d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals
    f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
    clip = 3.0
    dataDir += "tpd/csm/seismict/"
  else:
    print "do not recognize",name
    sys.exit(0)
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

#############################################################################
# Add displays to this simple function, as necessary.
def goDisplay(what):
  if what[0:3]=="atw":
    if what=="atw":
      print "must specify atw00 or atw01"
      sys.exit(0);
    setupForImage(what)
    g = readImage(what)
    print "atw: min=%f, max=%f"%(min(g),max(g))
    print "clip =",clip
    showOne(g,clip=clip)
  elif what=="can":
    setupForImage("can")
    g = readImage("canNZ")
    showOne(g,clip=clip)
  elif what=="f3d":
    setupForImage("f3d")
    g = readImage("f3draw")
    showOne(g,clip=clip)
  elif what=="gbc":
    setupForImage("gbc")
    g1 = readImage("pp")
    g2 = readImage("ps1")
    g3 = readImage("ps2")
    showThree(g1,g2,g3,clip=clip)
  elif what[0:3]=="mbs":
    if what=="mbs":
      print "must specify mbsSmall or mbsLarge"
      sys.exit(0);
    setupForImage(what)
    if what[3:8]=="Small":
      g = readImage("pstm_fraw_s1")
    else:
      g = readImage("pstm_raw_s1")
    showOne(g,clip-clip)
  elif what=="nor":
    setupForImage("nor")
    g = readImage("norne2006full")
    showOne(g,clip=clip)
  elif what=="pen":
    setupForImage("pen")
    g = readImage("pen")
    showOne(g,clip=clip)
  elif what[0:3]=="pnz":
    if what=="pnz":
      print "must specify pnz00 or pnz01 or ..."
      sys.exit(0)
    setupForImage("pnz")
    g = readImage(what)
    showOne(g,clip=clip)
  elif what=="seam":
    print "not implemented"
  elif what=="tpd":
    setupForImage("tpd")
    g = readImage("tpst")
    showOne(g,clip=clip)
  else:
    print "Sorry, do not know how to display ",what
    sys.exit(0)

#############################################################################
# Utilities

def main(args):
  if len(args)<2:
    print "Usage display.sh <what>, where <what> = f3d, pen, mbs, ..."
    return
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
