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

import tp.Convert as Convert

dataDir = "/data/seis/tp/"

# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header
n1i=1501 # number of time samples (1st dimension)
d1i=0.002 # time sampling interval
f1i=0.000 # time of first sample
n2i=188 # number of traces in 2nd (inline) dimension
d2i=0.033
f2i=0.000
n3i=345 # number of traces in 3rd (crossline) dimension
d3i=0.033
f3i=0.000

# Rotation and resampling of (x2,x3) coordinates
# Rotation center (in km) and angle (in radians)
# x2Old = x2c + x2New*cos(phi) + x3New*sin(phi)
# x3Old = x3c - x2New*sin(phi) + x3New*cos(phi)
x3c,x2c,phi = 0.935,4.192,-0.485364 # (-27.8093 degrees)
n1=1501; d1=0.002; f1=0.000
n2=161;  d2=0.025; f2=0.000
n3=357;  d3=0.025; f3=0.000

def main(args):
  testFormat()
  #readFormat()
  #readSegy()
  #resample()
  #ais = ArrayInputStream(dataDir+"tp3r.dat")
  #y = Array.zerofloat(n1,n2,n3)
  #ais.readFloats(y)
  #ais.close()
  #plot3d(y)
  #ais = ArrayInputStream(dataDir+"tp3.dat")
  #y = Array.zerofloat(n1i,n2i,n3i)
  #ais.readFloats(y)
  #ais.close()
  #plot3d(y)
  return

def resample():
  cphi = cos(phi)
  sphi = sin(phi)
  x2i = Array.zerofloat(n2,n3)
  x3i = Array.zerofloat(n2,n3)
  for i3 in range(n3):
    x3 = f3+i3*d3
    for i2 in range(n2):
      x2 = f2+i2*d2
      x2i[i3][i2] = x2c+x2*cphi+x3*sphi
      x3i[i3][i2] = x3c-x2*sphi+x3*cphi
  x = Array.zerofloat(n1i,n2i,n3i)
  y = Array.zerofloat(n1,n2,n3)
  x23 = Array.zerofloat(n2i,n3i)
  y23 = Array.zerofloat(n2,n3)
  sx = SimpleFloat3(x)
  sy = SimpleFloat3(y)
  si = SincInterpolator()
  si.setUniformSampling(n2i,d2i,f2i,n3i,d3i,f3i)
  ais = ArrayInputStream(dataDir+"tp3.dat")
  ais.readFloats(x)
  ais.close()
  print "x min/max =",Array.min(x),Array.max(x)
  for i1 in range(n1):
    if i1%100==0: print "i1 =",i1
    sx.get23(n2i,n3i,i1,0,0,x23)
    print "x23 min/max =",Array.min(x23),Array.max(x23)
    si.setUniformSamples(x23)
    for i3 in range(n3):
      for i2 in range(n2):
        y23[i3][i2] = si.interpolate(x2i[i3][i2],x3i[i3][i2])
    sy.set23(n2,n3,i1,0,0,y23)
  aos = ArrayOutputStream(dataDir+"tp3r.dat")
  aos.writeFloats(y)
  aos.close()
  plot3d(y)

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

def readFormat():
  infile = dataDir+"filt_mig.sgy"
  ais = ArrayInputStream(infile)
  ais.skipBytes(nhead)
# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
  h = Array.zeroshort(nbhed)
  ais.readShorts(h)
  print "format =",h[12]
  Array.dump(h)
  ais.close()

def testFormat():
  xi = Array.zeroint(n1i)
  x1 = Array.zerofloat(n1i)
  x2 = Array.zerofloat(n1i)
  infile = dataDir+"filt_mig.sgy"
  ais = ArrayInputStream(infile)
  ais.skipBytes(nhead+nbhed)
  ais.skipBytes(n3i/2*n2i*(nthed+4*n1i))
  ais.skipBytes(n2i/2*(nthed+4*n1i))
  ais.skipBytes(nthed)
  ais.readInts(xi)
  ais.close()
  Convert.ibmToFloat(xi,x1)
  Convert.ieeeToFloat(xi,x2)
  SimplePlot.asPoints(x1)
  SimplePlot.asPoints(x2)
  #Array.dump(xi)
  #Array.dump(x1)
  #Array.dump(x2)

def readSegy():
  infile = dataDir+"filt_mig.sgy"
  outfile = dataDir+"tp3.dat"
  ais = ArrayInputStream(infile)
  aos = ArrayOutputStream(outfile)
  ais.skipBytes(nhead+nbhed)
  x = Array.zerofloat(n1i)
  y = Array.zerofloat(n1i)
  for i in range(n2i*n3i):
    if i%1000==0:
      print "i =",i
    ais.skipBytes(nthed)
    ais.readFloats(x)
    Array.copy(n1i,x,y)
    #Array.dump(y)
    #print "y min =",Array.min(y)," max =",Array.max(y)
    aos.writeFloats(y)
  ais.close()
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
