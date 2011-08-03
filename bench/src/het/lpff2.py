import sys
#from math import *
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from het import *

pngDir = None
#pngDir = "png/het/"

#############################################################################

def main(args):
  #goImage()
  #goFilter()
  goHetMul()

def goHetMul():
  x,s1,s2,clip = getImageF3d(); png = "f3d_"
  #x,s1,s2,clip = getImageTpd(); png = "tpd_"
  plot(x,s1,s2,clip,png=png+"input")
  doHetMul(4,x,s1,s2,clip,png=png+"lpff")

def narrowBandFilter(fn,rn,x):
  nf = NotchFilter(fn,rn)
  y = like(x)
  nf.apply1ForwardReverse(x,y)
  return sub(x,y)
  
def doHetMul(sigma,x,s1,s2,clip,png=None):
  n1,n2 = s1.count,s2.count
  lpff = LocalPeakFrequencyFilter(sigma)
  f = like(x)
  for i2 in range(n2):
    f[i2] = lpff.findPeakFrequencies(x[i2])
  favg = sum(f)/n1/n2
  print "favg =",favg/s1.delta
  f = mul(f,1.0/s1.delta)
  plot(f,s1,s2)
  """
  #f = fillfloat(favg,n1,n2)
  w = mul(2.0*PI,sub(f,0.002))
  s = like(x)
  for i2 in range(n2):
    p = w[i2][0]
    for i1 in range(n1):
      s[i2][i1] = sin(p)
      p += w[i2][i1]
  """
  s = narrowBandFilter(favg,0.8,x)
  plot(s,s1,s2,0.1*clip)
  y = mul(s,x)
  y = mul(2,y)
  plot(y,s1,s2,clip)
  bf = ButterworthFilter(favg,6,ButterworthFilter.Type.LOW_PASS)
  z = like(y)
  bf.apply1Forward(y,z)
  bf.apply1Reverse(z,z)
  plot(z,s1,s2,clip)

def goFilter():
  #x,s1,s2,clip = getImageF3d(); png = "f3d_"
  x,s1,s2,clip = getImageTpd(); png = "tpd_"
  #x,s1,s2,clip = getImageAtw(); png = "atw_"
  plot(x,s1,s2,clip,png=png+"input")
  doFilter(20,x,s1,s2,clip,png=png+"lpff")
  
def doFilter(sigma,x,s1,s2,clip,png=None):
  n1,n2 = s1.count,s2.count
  lpff = LocalPeakFrequencyFilter(sigma)
  f = like(x)
  for i2 in range(n2):
    f[i2] = lpff.findPeakFrequencies(x[i2])
    mul(1.0/s1.delta,f[i2],f[i2])
  plot(f,s1,s2,cbar="Frequency (Hz)",png=png+"f")

def goImage():
  x,s1,s2,clip = getImage()
  plot(x,s1,s2,clip)

def getImage():
  return getImageF3d()
  #return getImageTpd()
  #return getImageSyn()

def getImageF3d():
  n1,n2 = 462,951
  d1,d2 = 0.004,0.025
  f1,f2 = 0.004,0.000
  fileName = "/data/seis/f3d/f3d75.dat"
  x = readImage(fileName,n1,n2)
  subset = True
  if subset:
    j1,j2 = 240,0
    n1,n2 = n1-j1,440
    f1,f2 = f1+j1*d1,f2+j2*d2
    x = copy(n1,n2,j1,j2,x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,5.0

def getImageTpd():
  n1,n2 = 251,357
  d1,d2 = 0.004,0.025
  f1,f2 = 0.500,0.000
  fileName = "/data/seis/tp/csm/oldslices/tp73.dat"
  x = readImage(fileName,n1,n2)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,2.0

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def like(x):
  n1,n2 = len(x[0]),len(x)
  return zerofloat(n1,n2)

#############################################################################
# plot

def plot(f,s1,s2,clip=0.0,t=None,cbar="",limits=None,png=None):
  n1,n2 = len(f[0]),len(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setFontSizeForPrint(8.0,240)
  sp.setSize(1020,730)
  if limits:
    sp.setHInterval(1.0)
    sp.setVInterval(0.1)
  else:
    sp.setHInterval(2.0)
    sp.setVInterval(0.2)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Time (s)")
  if limits:
    sp.setLimits(limits[0],limits[1],limits[2],limits[3])
  if cbar!=None:
    if len(cbar)>0:
      cbar = sp.addColorBar(cbar)
    else:
      cbar = sp.addColorBar()
  sp.plotPanel.setColorBarWidthMinimum(120)
  pv = sp.addPixels(s1,s2,f)
  if clip!=0.0:
    pv.setClips(-clip,clip)
    pv.setColorModel(ColorMap.GRAY)
  else:
    pv.setColorModel(ColorMap.JET)
  if png and pngDir:
    sp.paintToPng(720,3.3,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
