"""
Displays 3D images using the Mines Java Toolkit
"""
from common import *

#############################################################################
# time,trace,line sampling and other info in global variables

n1,n2,n3 = 751,1001,501
d1,d2,d3 = 1.0,1.0,1.0 # 0.004,0.0125,0.0125
f1,f2,f3 = 0.0,0.0,0.0
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
clip = 1e5
dataDir = "/data/seis/pnz/" # base directory for all images

#############################################################################

def main(args):
  #goDisplay(args[1])
  #goDisplay("f")
  #goDisplay("gs")
  goDisplay("ggs")
  #goSlice()

def goDisplay(what):
  f = readImage(what)
  print "f: min=%f, max=%f"%(min(f),max(f))
  print "clip =",clip
  showOne(f,clip=clip)

def goSlice():
  f = readImage("ggs")
  m2,m3 = 501,501
  j2,j3 = 165,0
  g = zerofloat(m2,m3)
  for i3 in range(j3,j3+m3):
    for i2 in range(j2,j2+m2):
      g[i3-j3][i2-j2] = f[i3][i2][256]
  writeImage("ggs256",g)
  SimplePlot.asPixels(g)

#############################################################################
# Utilities

def readImage(name):
  fileName = dataDir+name+".dat"
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+name+".dat")
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(name,f):
  fileName = dataDir+name+".dat"
  aos = ArrayOutputStream(dataDir+name+".dat")
  aos.writeFloats(f)
  aos.close()

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
