# Make good subsets for fault processing

from common import *

n1,n2,n3 = 1501,4001,3001 # numbers of samples
d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals (for now)
f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

dataDir = "/data/dhale/can/"
clip = 1.0

def main(args):
  #subset1(True)
  #subset1b(True)
  #subset1c(True)
  pass

def subset1(make):
  """
  s1 = 3.2 GB volume with interesting faults
  use this volume for scanning and thinning
  """
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  j1,j2,j3 = 450,1463,0 # first samples in subset
  m1,m2,m3 = 501,1601,1001 # numbers of samples in subset
  if make:
    x = zerofloat(m1)
    aif = ArrayFile("/data/seis/can/canNZ.dat","r")
    aos = ArrayOutputStream(dataDir+"s1/g.dat")
    for i3 in range(j3,j3+m3):
      for i2 in range(j2,j2+m2):
        aif.seek(4*(j1+n1*(i2+n2*i3)))
        aif.readFloats(x)
        aos.writeFloats(x)
    aos.close()
    aif.close()
  n1,n2,n3 = m1,m2,m3
  d1,d2,d3 = d1,d2,d3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

def subset1b(make):
  """
  s1b = bottom part of subset s1, below chaotic faulting
  used for surface extraction and throw estimation
  make this subset after scanning and thinning
  """
  subset1(False)
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  j1,j2,j3 = 200,0,0 # first samples in subset
  m1,m2,m3 = 301,1601,1001 # numbers of samples in subset
  if make:
    x = zerofloat(n1)
    y = zerofloat(m1)
    for what in ["g","gs","p2","p3","fl","fp","ft","flt","fpt","ftt"]:
      ais = ArrayInputStream(dataDir+"s1/"+what+".dat")
      aos = ArrayOutputStream(dataDir+"s1b/"+what+".dat")
      for i3 in range(n3):
        for i2 in range(n2):
          ais.readFloats(x)
          copy(m1,j1,x,0,y)
          aos.writeFloats(y)
      ais.close()
      aos.close()
  n1,n2,n3 = m1,m2,m3
  d1,d2,d3 = d1,d2,d3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

def subset1c(make):
  """
  s1c = small part of subset s1b, below chaotic faulting
  used for surface extraction and throw estimation
  make this subset after scanning and thinning
  """
  subset1b(False)
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  j1,j2,j3 =   0,350,230 # first samples in subset
  m1,m2,m3 = 301,301,301 # numbers of samples in subset
  if make:
    x = zerofloat(n1)
    y = zerofloat(m1)
    for what in ["g","gs","p2","p3","fl","fp","ft","flt","fpt","ftt"]:
      ais = ArrayInputStream(dataDir+"s1b/"+what+".dat")
      aos = ArrayOutputStream(dataDir+"s1c/"+what+".dat")
      ais.skipBytes(4*n1*n2*j3)
      for i3 in range(j3,j3+m3):
        ais.skipBytes(4*n1*j2)
        for i2 in range(j2,j2+m2):
          ais.readFloats(x)
          copy(m1,j1,x,0,y)
          aos.writeFloats(y)
        ais.skipBytes(4*n1*(n2-j2-m2))
      ais.skipBytes(4*n1*n2*(n3-j3-m3))
      ais.close()
      aos.close()
  n1,n2,n3 = m1,m2,m3
  d1,d2,d3 = d1,d2,d3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

#############################################################################

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
