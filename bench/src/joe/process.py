import sys

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

from lcc import *
from joe import *

datadir = "/data/seis/joe/"

# 
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

# sampling before windowing to get rid of wierd edges
#k1=750 # index of first sample in subvolume
#n1=251; d1=d1i; f1=f1i+k1*d1 # time sampling in subvolume
#n2=500; d2=d2i; f2=f2i
#n3=500; d3=d3i; f3=f3i

# sampling after windowing
#n1,n2,n3 = 251,400,500
#d1,d2,d3 = d1i,d2i,d3i
#f1,f2,f3 = f1i,f2i,f3i

# sampling after windowing
n1,n2,n3 = 200,200,200
d1,d2,d3 = 1.0,1.0,1.0
f1,f2,f3 = 0.0,0.0,0.0

def main(args):
  #window()
  #x = readFile("x.dat",n1,n2,n3)
  #plot3d(x)
  #u2,u3,ep = planes()
  #x = readFile("x.dat",n1,n2,n3)
  #u2 = readFile("u2.dat",n1,n2,n3)
  #ep = readFile("ep.dat",n1,n2,n3)
  #plot3ds((x,u2))
  #plot3ds((x,ep))
  #ldf()
  #y = readFile("y.dat",n1,n2,n3)
  #w1,w2,w3,el = lines()
  #plot3ds((y,el))
  llf()
  return

def ldf():
  x = readFile("x.dat",n1,n2,n3)
  u2 = readFile("xu2.dat",n1,n2,n3)
  u3 = readFile("xu3.dat",n1,n2,n3)
  ldf = LocalDipFilter()
  y = Array.zerofloat(n1,n2,n3)
  #ldf.applyDip(0.05,u2,u3,x,y)
  ldf.applyNotch(0.01,u2,u3,x,y)
  writeFile("y.dat",y)
  plot3ds((x,y))

def llf():
  x = readFile("x.dat",n1,n2,n3)
  y = readFile("y.dat",n1,n2,n3)
  w1 = readFile("yw1.dat",n1,n2,n3)
  w2 = readFile("yw2.dat",n1,n2,n3)
  w3 = readFile("yw3.dat",n1,n2,n3)
  el = readFile("yel.dat",n1,n2,n3)
  print "w1 min/max =",Array.min(w1),Array.max(w1)
  print "w2 min/max =",Array.min(w2),Array.max(w2)
  print "w3 min/max =",Array.min(w3),Array.max(w3)
  llf = LocalLineFilter()
  z = Array.zerofloat(n1,n2,n3)
  llf.applyLine(0.1,w1,w2,w3,y,z)
  #llf.applyNotch(0.01,w1,w2,w3,y,z)
  z = Array.mul(el,Array.sub(y,z))
  writeFile("z.dat",z)
  plot3ds((x,z))

def window():
  x = readFile("win34.dat",251,500,500)
  #j1,j2,j3 =  0, 50, 0
  #m1,m2,m3 = n1,400,n3
  j1,j2,j3 =   0,150,250
  m1,m2,m3 = 200,200,200
  y = Array.copy(m1,m2,m3,j1,j2,j3,x)
  writeFile("x.dat",y)

def planes():
  x = readFile("x.dat",n1,n2,n3)
  lof = LocalOrientFilter(12.0)
  u2 = Array.zerofloat(n1,n2,n3)
  u3 = Array.zerofloat(n1,n2,n3)
  ep = Array.zerofloat(n1,n2,n3)
  lof.apply(x,
    None,None,
    None,u2,u3,
    None,None,None,
    None,None,None,
    None,None,None,
    ep,None)
  writeFile("xu2.dat",u2)
  writeFile("xu3.dat",u3)
  writeFile("xep.dat",ep)
  return u2,u3,ep

def lines():
  x = readFile("y.dat",n1,n2,n3)
  lof = LocalOrientFilter(6.0)
  w1 = Array.zerofloat(n1,n2,n3)
  w2 = Array.zerofloat(n1,n2,n3)
  w3 = Array.zerofloat(n1,n2,n3)
  el = Array.zerofloat(n1,n2,n3)
  lof.apply(x,
    None,None,
    None,None,None,
    None,None,None,
    w1,w2,w3,
    None,None,None,
    None,el)
  writeFile("yw1.dat",w1)
  writeFile("yw2.dat",w2)
  writeFile("yw3.dat",w3)
  writeFile("yel.dat",el)
  return w1,w2,w3,el

def plot12(i3):
  ais = ArrayInputStream(datadir+"win34.dat")
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

def plot3ds(xs):
  world = World()
  for x in xs:
    print "x min =",Array.min(x)," max =",Array.max(x)
    n1 = len(x[0][0])
    n2 = len(x[0])
    n3 = len(x)
    s1 = Sampling(n1)
    s2 = Sampling(n2)
    s3 = Sampling(n3)
    ipg = ImagePanelGroup(s3,s2,s1,SimpleFloat3(x))
    ipg.setPercentiles(1,99)
    world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def plot3d(x,cmin=-0.01,cmax=0.01):
  print "x min =",Array.min(x)," max =",Array.max(x)
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  s3 = Sampling(n3,d3,f3)
  x3 = SimpleFloat3(x)
  ipg = ImagePanelGroup(s3,s2,s1,x3)
  ipg.setClips(cmin,cmax)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def readFile(file,n1,n2,n3):
  ais = ArrayInputStream(datadir+file)
  x = Array.zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def writeFile(file,x):
  aos = ArrayOutputStream(datadir+file)
  aos.writeFloats(x)
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
