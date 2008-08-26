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
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *

from ldf import *

True = 1
False = 0

#############################################################################
# parameters
n1,n2,n3 = 105,105,105
#n1,n2,n3 = 195,195,195
#n1,n2,n3 = 315,315,315
sigma = 16
ffile = "filters.dat"
small = 0.01
niter = 100


#############################################################################
# functions

def main(args):
  #combineFilters()
  #doAmp(0,0)
  x = smooth(makeRandom(n1,n2,n3))
  doDiff(x)

def makeImpulse(n1,n2,n3):
  x = Array.zerofloat(n1,n2,n3)
  x[(n3-1)/2][(n2-1)/2][(n1-1)/2] = 1.0
  return x

def makeRandom(n1,n2,n3):
  r = Random(314159)
  return Array.sub(Array.randfloat(r,n1,n2,n3),0.5)

def smooth(x):
  n1 = len(x[0][0])
  n2 = len(x[0])
  n3 = len(x)
  y = Array.zerofloat(n1,n2,n3)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply000(x,y)
  return y

def makeVectorsConstant(theta,phi,n1,n2,n3):
  uss = UnitSphereSampling(16);
  theta *= pi/180
  phi *= pi/180
  v1 = cos(theta)
  v2 = sin(theta)*sin(phi)
  v3 = sin(theta)*cos(phi)
  kv = uss.getIndex(v3,v2,v1)
  print "v1 =",v1," v2 =",v2," v3 =",v3," kv =",kv
  return Array.fillshort(kv,n1,n2,n3);

def makeVectorsRadial(n1,n2,n3):
  k1,k2,k3 = n1/2,n2/2,n3/2
  uss = UnitSphereSampling(16)
  iv = Array.zeroshort(n1,n2,n3)
  for i3 in range(n3):
    v3 = float(i3-k3)
    for i2 in range(n2):
      v2 = float(i2-k2)
      for i1 in range(n1):
        v1 = 0.1+float(i1-k1)
        sv = 1.0/sqrt(v1*v1+v2*v2+v3*v3)
        u1,u2,u3 = sv*v1,sv*v2,sv*v3
        iv[i3][i2][i1] = uss.getIndex(u3,u2,u1)
  return iv

def doDiff(x):
  ds = None
  iv = makeVectorsRadial(n1,n2,n3)
  #iv = makeVectorsConstant(90,0,n1,n2,n3)
  y = Array.copy(x)
  #LocalDiffusionFilterMp.setFiltersFile(ffile)
  #ldf = LocalDiffusionFilterMp(sigma)
  ldf = LocalDiffusionFilterCg(sigma,small,niter)
  ldf.applyLinearPass(ds,iv,x,y)
  #ldf.applyPlanarPass(ds,iv,x,y)
  plot3d(x)
  plot3d(y)
  

def doAmp(theta,phi):
  #ldf = LocalDiffusionFilterMp(sigma)
  ldf = LocalDiffusionFilterMp(sigma,ffile)
  #ldf = LocalDiffusionFilterCg(sigma,small,niter)
  x = makeImpulse(n1,n2,n3)
  ds = None
  iv = makeVectorsConstant(theta,phi,n1,n2,n3)
  h = Array.zerofloat(n1,n2,n3)
  #ldf.applyInlinePass(ds,iv,x,h)
  ldf.applyNormalPass(ds,iv,x,h)
  ah = frequencyResponse(h)
  plot3d(ah)

def frequencyResponse(x):
  n1 = len(x[0][0])
  n2 = len(x[0])
  n3 = len(x)
  n1 = FftComplex.nfftSmall(n1)
  n2 = FftComplex.nfftSmall(n2)
  n3 = FftComplex.nfftSmall(n3)
  xr = Array.copy(n1,n2,n3,x)
  xi = Array.zerofloat(n1,n2,n3)
  cx = Array.cmplx(xr,xi)
  fft1 = FftComplex(n1)
  fft2 = FftComplex(n2)
  fft3 = FftComplex(n3)
  fft1.complexToComplex1(1,n2,n3,cx,cx)
  fft2.complexToComplex2(1,n1,n3,cx,cx)
  fft3.complexToComplex3(1,n1,n2,cx,cx)
  ax = Array.cabs(cx)
  a = Array.zerofloat(n1,n2,n3)
  j1 = n1/2
  j2 = n2/2
  j3 = n3/2
  Array.copy(n1-j1,n2-j2,n3-j3,0,0,0,ax,j1,j2,j3,a)
  Array.copy(j1,n2-j2,n3-j3,n1-j1,0,0,ax,0,j2,j3,a)
  Array.copy(n1-j1,j2,n3-j3,0,n2-j2,0,ax,j1,0,j3,a)
  Array.copy(n1-j1,n2-j2,j3,0,0,n3-j3,ax,j1,j2,0,a)
  Array.copy(j1,j2,n3-j3,n1-j1,n2-j2,0,ax,0,0,j3,a)
  Array.copy(n1-j1,j2,j3,0,n2-j2,n3-j3,ax,j1,0,0,a)
  Array.copy(j1,n2-j2,j3,n1-j1,0,n3-j3,ax,0,j2,0,a)
  Array.copy(j1,j2,j3,n1-j1,n2-j2,n3-j3,ax,0,0,0,a)
  return a

def combineFilters():
  n2 = 4*8*21*9
  d2 = Array.zerofloat(n2/4)
  af2 = ArrayFile("filters2.dat","r")
  af2.readInt()
  af2.readInt()
  af2.readInt()
  af2.readFloats(d2)
  af2.close()
  n3 = 4*56*481*9
  d3 = Array.zerofloat(n3/4)
  af3 = ArrayFile("filters3.dat","r")
  af3.readInt()
  af3.readFloats(d3)
  af3.close()
  af = ArrayFile("filters.dat","rw")
  af.writeInt(1)
  af.writeInt(2)
  af.writeInt(n2)
  af.writeFloats(d2)
  af.writeInt(3)
  af.writeInt(n3)
  af.writeFloats(d3)
  af.close()

def makeFilters():
  ldf = LocalDiffusionFilterMp(sigma)
  ldf.save(ffile)

#############################################################################
# plot

def plot3d(x):
  print "x min =",Array.min(x)," max =",Array.max(x)
  ipg = ImagePanelGroup(x)
  #ipg.setColorModel(ColorMap.JET)
  ipg.setColorModel(ColorMap.GRAY)
  ipg.setClips(-0.1,0.1)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
