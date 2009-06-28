import sys

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

import joe.Convert as Convert

dataDir = "/data/seis/joe/"

# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header
n1i=1051 # number of time samples (1st dimension)
d1i=0.004 # time sampling interval
f1i=0.000 # time of first sample
n2i=500 # number of traces in 2nd dimension
d2i=0.005
f2i=0.000
n3i=500 # number of traces in 3rd dimension
d3i=0.005
f3i=0.000

k1=750 # index of first sample in subvolume
n1=251; d1=d1i; f1=f1i+k1*d1 # time sampling in subvolume
n2=500; d2=d2i; f2=f2i
n3=500; d3=d3i; f3=f3i

def main(args):
  #testFormat()
  readFormat()
  #readSegy()
  #x = readImage()
  #plot3d(x)
  return

def plot12(i3):
  ais = ArrayInputStream(dataDir+"win34.dat")
  x = zerofloat(n1,n2)
  ais.skipBytes(4*n1*n2*i3)
  ais.readFloats(x)
  ais.close()
  print "x min =",min(x)," max =",max(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(x)
  #pv.setPercentiles(1.0,99.0)
  pv.setClips(-1.0e-6,1.0e-6)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)

def plot3d(x):
  print "x min =",min(x)," max =",max(x)
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  s3 = Sampling(n3,d3,f3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  #clip = 1.0e-6
  #ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def readFormat():
  infile = dataDir+"QuinSy_migPh0time.sgy"
  ais = ArrayInputStream(infile)
  ais.skipBytes(nhead)
# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
  h = zeroshort(nbhed)
  ais.readShorts(h)
  print "format =",h[12]
  ais.close()

def testFormat():
  infile = dataDir+"QuinSy_migPh0time.sgy"
  ais = ArrayInputStream(infile)
  ais.skipBytes(nhead+nbhed)
  ais.skipBytes(nthed)
  xi = zeroint(n1i)
  ais.readInts(xi)
  ais.close()
  x1 = zerofloat(n1i)
  x2 = zerofloat(n1i)
  Convert.ibmToFloat(xi,x1)
  Convert.ieeeToFloat(xi,x2)
  SimplePlot.asPoints(x1)
  SimplePlot.asPoints(x2)
  #dump(x1)
  #dump(x2)

def readSegy():
  infile = dataDir+"QuinSy_migPh0time.sgy"
  outfile = dataDir+"win34.dat"
  ais = ArrayInputStream(infile)
  aos = ArrayOutputStream(outfile)
  ais.skipBytes(nhead+nbhed)
  x = zerofloat(n1i)
  y = zerofloat(n1)
  for i in range(n2i*n3i):
    if i%1000==0:
      print "i =",i
    ais.skipBytes(nthed)
    ais.readFloats(x)
    copy(n1,k1,x,0,y)
    #dump(y)
    #print "y min =",min(y)," max =",max(y)
    aos.writeFloats(y)
  ais.close()
  aos.close()

def readImage():
  ais = ArrayInputStream(dataDir+"win34.dat")
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
