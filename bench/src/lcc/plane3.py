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
from edu.mines.jtk.util.ArrayMath import *

from lcc import *

#############################################################################
# parameters


#############################################################################
# functions

def main(args):
  #doAmp(30,90)
  doSymTest()
  return

def makeImpulse(dip,azi):
  n1,n2,n3 = 105,105,105
  dip *= pi/180.0
  azi *= pi/180.0
  u1c = cos(dip)
  u2c = -sin(dip)*cos(azi)
  u3c = -sin(dip)*sin(azi)
  u1 = fillfloat(u1c,n1,n2,n3)
  u2 = fillfloat(u2c,n1,n2,n3)
  u3 = fillfloat(u3c,n1,n2,n3)
  x = zerofloat(n1,n2,n3)
  x[(n3-1)/2][(n2-1)/2][(n1-1)/2] = 1.0
  return x,u1,u2,u3

def doAmp(dip,azi):
  ldf = LocalDipFilter()
  x,u1,u2,u3 = makeImpulse(dip,azi)
  h = copy(x)
  ldf.applyDip(0.05,u2,u3,x,h)
  ah = frequencyResponse(h)
  plot3d(ah)

def doSymTest():
  n1,n2,n3 = 105,105,105
  ldf = LocalDipFilter()
  u2 = mul(0.5,randfloat(n1,n2,n3))
  u3 = mul(0.5,randfloat(n1,n2,n3))
  x = sub(randfloat(n1,n2,n3),0.5)
  y = sub(randfloat(n1,n2,n3),0.5)
  ax = zerofloat(n1,n2,n3)
  ay = zerofloat(n1,n2,n3)
  ldf.applyForward(0.1,0.1,u2,u3,x,ax)
  ldf.applyForward(0.1,0.1,u2,u3,y,ay)
  yax = sum(mul(y,ax))
  xay = sum(mul(x,ay))
  print "yax =",yax," xay=",xay

def frequencyResponse(x):
  n1 = len(x[0][0])
  n2 = len(x[0])
  n3 = len(x)
  n1 = FftComplex.nfftSmall(n1)
  n2 = FftComplex.nfftSmall(n2)
  n3 = FftComplex.nfftSmall(n3)
  xr = copy(n1,n2,n3,x)
  xi = zerofloat(n1,n2,n3)
  cx = cmplx(xr,xi)
  fft1 = FftComplex(n1)
  fft2 = FftComplex(n2)
  fft3 = FftComplex(n3)
  fft1.complexToComplex1(1,n2,n3,cx,cx)
  fft2.complexToComplex2(1,n1,n3,cx,cx)
  fft3.complexToComplex3(1,n1,n2,cx,cx)
  ax = cabs(cx)
  a = zerofloat(n1,n2,n3)
  j1 = n1/2
  j2 = n2/2
  j3 = n3/2
  copy(n1-j1,n2-j2,n3-j3,0,0,0,ax,j1,j2,j3,a)
  copy(j1,n2-j2,n3-j3,n1-j1,0,0,ax,0,j2,j3,a)
  copy(n1-j1,j2,n3-j3,0,n2-j2,0,ax,j1,0,j3,a)
  copy(n1-j1,n2-j2,j3,0,0,n3-j3,ax,j1,j2,0,a)
  copy(j1,j2,n3-j3,n1-j1,n2-j2,0,ax,0,0,j3,a)
  copy(n1-j1,j2,j3,0,n2-j2,n3-j3,ax,j1,0,0,a)
  copy(j1,n2-j2,j3,n1-j1,0,n3-j3,ax,0,j2,0,a)
  copy(j1,j2,j3,n1-j1,n2-j2,n3-j3,ax,0,0,0,a)
  return a

#############################################################################
# plot

def plot3d(x):
  print "x min =",min(x)," max =",max(x)
  ipg = ImagePanelGroup(x)
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
