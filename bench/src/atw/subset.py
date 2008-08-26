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

import atw.Convert as Convert

dataDir = "/data/seis/atw/"

# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header
n1i=851 # number of time samples (1st dimension)
d1i=0.0040 # time sampling interval
f1i=0.0000 # time of first sample
n2i=1932 # number of traces in 2nd (inline) dimension (24.15 km ~ 15 mi)
d2i=0.0125
f2i=0.0000
n3i=1207 # number of traces in 3rd (crossline) dimension (24.14 km ~ 15 mi)
d3i=0.0200
f3i=0.0000

# Subset parameters.
#n1= 851; d1=0.0040; f1=0.0000
#n2= 643; d2=0.0125; f2=0.0000
#n3= 402; d3=0.0200; f3=0.0000
#i1s,i2s,i3s = 0,1000,350 # shallow prograding sequence
i1s,i2s,i3s = 0,800,750
n1= 851; d1=0.0040; f1=0.0000
n2= 641; d2=0.0125; f2=0.0000
n3= 401; d3=0.0200; f3=0.0000

def main(args):
  #plot12(0);
  #plot12(100);
  #testFormat()
  #subset(i1s,i2s,i3s)
  ais = ArrayInputStream(dataDir+"atws.dat")
  y = Array.zerofloat(n1,n2,n3)
  ais.readFloats(y)
  ais.close()
  plot3d(y)
  return

def subset(i1s,i2s,i3s):
  ais = ArrayInputStream(dataDir+"atw.dat")
  aos = ArrayOutputStream(dataDir+"atws.dat")
  x = Array.zerofloat(n1)
  for i3 in range(i3s):
    ais.skipBytes(4*n1i*n2i)
  for i3 in range(n3):
    ais.skipBytes(4*n1i*i2s)
    for i2 in range(n2):
      ais.skipBytes(4*i1s)
      ais.readFloats(x)
      aos.writeFloats(x)
      ais.skipBytes(4*(n1i-i1s-n1))
    ais.skipBytes(4*n1i*(n2i-i2s-n2))
  ais.close()
  aos.close()

def plot12(i3):
  ais = ArrayInputStream(dataDir+"atw.dat")
  x = Array.zerofloat(n1,n2)
  ais.skipBytes(4*n1*n2*i3)
  ais.readFloats(x)
  ais.close()
  print "x min =",Array.min(x)," max =",Array.max(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(x)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  #pv.setPercentiles(1.0,99.0)
  pv.setClips(-3.5e4,3.5e4)

def plot3d(x):
  print "x min =",Array.min(x)," max =",Array.max(x)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  #s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  clip = 2.5e4
  ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def testFormat():
  xi = Array.zeroint(n1i,n2i)
  x1 = Array.zerofloat(n1i,n2i)
  x2 = Array.zerofloat(n1i,n2i)
  infile = dataDir+"atw.dat";
  ais = ArrayInputStream(infile)
  ais.readInts(xi)
  ais.close()
  for i2 in range(n2i):
    Convert.ibmToFloat(xi[i2],x1[i2])
    Convert.ieeeToFloat(xi[i2],x2[i2])
  SimplePlot.asPixels(x1)
  SimplePlot.asPixels(x2)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
