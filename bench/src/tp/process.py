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
from edu.mines.jtk.sgl.test import *

dataDir = "/data/seis/tp/"

# Rotated and resampled data.
#n1=1501; d1=0.002; f1=0.000
#n2=161;  d2=0.025; f2=0.000
#n3=357;  d3=0.025; f3=0.000

# After subset and subsampling in time.
n1=251; d1=0.004; f1=0.500
n2=161; d2=0.025; f2=0.000
n3=357; d3=0.025; f3=0.000

def main(args):
  x = readSubset()
  plot3d(x)
  return

def readSubset():
  ais = ArrayInputStream(dataDir+"tp3s.dat")
  x = Array.zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def subset():
  m1,m2,m3 = 251,n2,n3
  j1,j2,j3 = 251,0,0
  x = Array.zerofloat(n1,n2,n3)
  y = Array.zerofloat(2*m1,m2,m3)
  z = Array.zerofloat(m1,m2,m3)
  ais = ArrayInputStream(dataDir+"tp3r.dat")
  ais.readFloats(x)
  ais.close()
  sx = SimpleFloat3(x)
  sx.get123(2*m1,m2,m3,j1,j2,j3,y)
  Array.copy(m1,m2,m3,0,0,0,2,1,1,y,0,0,0,1,1,1,z)
  aos = ArrayOutputStream(dataDir+"tp3s.dat")
  aos.writeFloats(z)
  aos.close()
  plot3d(z)

def plot12(i3):
  ais = ArrayInputStream(dataDir+"tp3.dat")
  x = Array.zerofloat(n1,n2)
  ais.skipBytes(4*n1*n2*i3)
  ais.readFloats(x)
  ais.close()
  print "x min =",Array.min(x)," max =",Array.max(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(x)
  #pv.setPercentiles(1.0,99.0)
  pv.setClips(-1.0e-6,1.0e-6)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)

def plot3d(x):
  print "x min =",Array.min(x)," max =",Array.max(x)
  #s1 = Sampling(n1,d1,f1)
  #s2 = Sampling(n2,d2,f2)
  #s3 = Sampling(n3,d3,f3)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  s1,s2,s3 = Sampling(n1,1.0,0.0),Sampling(n2,1.0,0.0),Sampling(n3,1.0,0.0)
  x3 = SimpleFloat3(x)
  ipg = ImagePanelGroup(s3,s2,s1,x3)
  #clip = 1.0e-6
  #ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
