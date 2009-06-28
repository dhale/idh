import sys
from math import *
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

from lcc import *

#############################################################################
# parameters

fontSize = 24
width = 640
height = 505
widthColorBar = 80
dataDir = "/data"
pngDir = None

n1 = 315
n2 = 315
order = 1
sigma = 4

#############################################################################
# functions

def main(args):
  goLattice()
  return

def goLattice():
  x = doImage()
  print "xrms =",rms(x)
  #doForwardInverse14(x,order,sigma)
  #doForwardInverse41(x,order,sigma)
  doTexture(x,order,sigma)

def doImage():
  x = readImage()
  #x = transpose(x)
  plot(x,10.0,"x")
  return x

def doTexture(x,order,sigma):
  y = zerofloat(n1,n2)
  z13 = zerofloat(n1,n2)
  z24 = zerofloat(n1,n2)
  llf = LocalLatticeFilter(sigma)
  cq1 = llf.findQ1(order,x)
  cq2 = llf.findQ2(order,x)
  cq3 = llf.findQ3(order,x)
  cq4 = llf.findQ4(order,x)
  r = sub(randfloat(n1,n2),0.5)
  #r = zerofloat(n1,n2)
  #r[n2/2][n1/3] = 1.0
  llf.applyQ1Inverse(cq1,r,z13)
  llf.applyQ3Inverse(cq3,z13,z13)
  plot(z13,0.0,None)
  llf.applyQ2Inverse(cq2,r,z24)
  llf.applyQ4Inverse(cq4,z24,z24)
  plot(z24,0.0,None)
  z = add(z13,z24)
  plot(z,0.0,None)
  llf.applyQ3Inverse(cq3,r,z13)
  llf.applyQ1Inverse(cq1,z13,z13)
  plot(z13,0.0,None)
  llf.applyQ4Inverse(cq4,r,z24)
  llf.applyQ2Inverse(cq2,z24,z24)
  plot(z24,0.0,None)
  z = add(z13,z24)
  plot(z,0.0,None)

def doForwardInverse14(x,order,sigma):
  y = zerofloat(n1,n2)
  llf = LocalLatticeFilter(sigma)
  cq1 = llf.findQ1(order,x)
  llf.applyQ1Forward(cq1,x,y)
  cq4 = llf.findQ4(order,y)
  llf.applyQ4Forward(cq4,y,y)
  print "yrms =",rms(y)
  plot(y,2.0,"y")
  z = zerofloat(n1,n2)
  llf.applyQ4Inverse(cq4,y,z)
  llf.applyQ1Inverse(cq1,z,z)
  #plot(z,10.0,"z")
  print "max |z-x|:",max(abs(sub(z,x)))
  r = sub(randfloat(n1,n2),0.5)
  s = zerofloat(n1,n2)
  llf.applyQ4Inverse(cq4,r,s)
  llf.applyQ1Inverse(cq1,s,s)
  plot(s,0.0,"s")

def doForwardInverse41(x,order,sigma):
  y = zerofloat(n1,n2)
  llf = LocalLatticeFilter(sigma)
  cq4 = llf.findQ4(order,x)
  llf.applyQ4Forward(cq4,x,y)
  cq1 = llf.findQ1(order,y)
  llf.applyQ1Forward(cq1,y,y)
  print "yrms =",rms(y)
  plot(y,2.0,"y")
  z = zerofloat(n1,n2)
  llf.applyQ1Inverse(cq1,y,z)
  llf.applyQ4Inverse(cq4,z,z)
  #plot(z,10.0,"z")
  print "max |z-x|:",max(abs(sub(z,x)))
  r = sub(randfloat(n1,n2),0.5)
  s = zerofloat(n1,n2)
  llf.applyQ1Inverse(cq1,r,s)
  llf.applyQ4Inverse(cq4,s,s)
  plot(s,0.0,"s")

def rms(x):
  n1 = len(x[0])
  n2 = len(x)
  s = 1.0/(n1*n2)
  return sqrt(s*sum(mul(x,x)))

def plotc(c):
  m = order
  c1 = zerofloat(n1,n2)
  c2 = zerofloat(n1,n2)
  for k in range(m):
    for i2 in range(n2):
      for i1 in range(n1):
        c1[i2][i1] = c[i2][i1*m*2+k*2  ]
        c2[i2][i1] = c[i2][i1*m*2+k*2+1]
    plot(c1,1.0,"c1")
    plot(c2,1.0,"c2")

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def flip2(f):
  n1 = len(f[0])
  n2 = len(f)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    copy(f[n2-1-i2],g[i2])
  return g

#############################################################################
# plot

def plot(f,clip=0.0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  if n1<50 and n2<50:
    s1 = Sampling(n1,1,-(n1-1)/2)
    s2 = Sampling(n2,1,-(n2-1)/2)
  pv = p.addPixels(s1,s2,f)
  if clip!=0.0:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  frame(p,png)

def panel():
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p.addColorBar()
  p.setColorBarWidthMinimum(widthColorBar)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(200,6,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
