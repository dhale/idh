import sys
from math import *

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

dataDir = "/data/seis/atw/"

# John Mathewson's subsets, resampled to 20 x 20 m trace spacing
# atwj1.dat = John's area1 file (i1=40 is a good 2D slice)
# atwj3.dat = John's area3 file
n1= 129; d1=0.0040; f1=0.0000
n2= 500; d2=0.0200; f2=0.0000
n3= 500; d3=0.0200; f3=0.0000

def main(args):
  doPlot3d()
  #doSlice1()
  return

def doPlot3d():
  x = readImage("atwj1.dat")
  plot3d(x)

def doSlice1():
  x = readImage("atwj1.dat")
  y = slice1(40,x)
  writeImage(y,"atwj1s.dat")

def readImage(fileName):
  ais = ArrayInputStream(dataDir+fileName)
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(x,fileName):
  aos = ArrayOutputStream(dataDir+fileName)
  aos.writeFloats(x)
  aos.close()

def slice1(k1,x):
  s = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      s[i3][i2] = x[i3][i2][k1]
  SimplePlot.asPixels(s)
  return s

def plot3d(x):
  print "x min =",min(x)," max =",max(x)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  #s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  clip = 2.5e4
  ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
