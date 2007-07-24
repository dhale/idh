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
sigma = 16
ffile = "filters.dat"


#############################################################################
# functions

def main(args):
  doAmp(0,0)
  return

def makeImpulse(n1,n2,n3):
  x = Array.zerofloat(n1,n2,n3)
  x[(n3-1)/2][(n2-1)/2][(n1-1)/2] = 1.0
  return x

def makeVectors(dip,azi,n1,n2,n3):
  uss = UnitSphereSampling(16);
  dip *= pi/180.0
  azi *= pi/180.0
  v1 = cos(dip)
  v2 = -sin(dip)*cos(azi)
  v3 = -sin(dip)*sin(azi)
  iv = Array.fillshort(uss.getIndex(v3,v2,v1),n1,n2,n3);
  return iv

def doAmp(dip,azi):
  ldf = LocalDiffusionFilterMp(sigma,ffile)
  x = makeImpulse(n1,n2,n3)
  ds = None
  iv = makeVectors(dip,azi,n1,n2,n3)
  h = Array.zerofloat(n1,n2,n3)
  ldf.applyInlinePass(ds,iv,x,h)
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

#############################################################################
# plot

def plot3d(x):
  print "x min =",Array.min(x)," max =",Array.max(x)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  ipg = ImagePanelGroup(s1,s2,s3,SimpleFloat3(x))
  ipg.setColorModel(ColorMap.JET)
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
