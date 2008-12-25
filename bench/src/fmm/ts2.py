#############################################################################
# Tests 2D time solvers.

import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

from fmm import *

#############################################################################
# global parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

#############################################################################
# test functions

def main(args):
  #testConstant()
  testSine()
  #testTsai()
  #testSeismic()
  return
 
def testConstant():
  n1,n2 = 1001,1001
  tensors = ConstantTensors2(n1,n2,20.0,0.010,1.000)
  for ts in [TimeSolver2(n1,n2,tensors),
             AfmmSolver2(n1,n2,tensors)]:
    ts.zeroAt(2*(n1-1)/4,2*(n2-1)/4)
    plot(ts.getTimes(),0,0,prism)
 
def testSine():
  n1,n2 = 601,601
  tensors = SineTensors2(n1,n2,400.0,30.0)
  i1 = [7*n1/8,3*n1/8,5*n1/8,2*n1/8]
  i2 = [1*n2/8,4*n2/8,4*n2/8,6*n2/8]
  for ts in [TimeSolver2(n1,n2,tensors),
             AfmmSolver2(n1,n2,tensors)]:
    #ts.zeroAt(2*(n1-1)/4,2*(n2-1)/4)
    for i in range(len(i1)):
      ts.zeroAt(i1[i],i2[i])
    plot(ts.getTimes(),0,0,prism)
 
def testTsai():
  n1,n2 = 101,101
  tensors = TsaiTensors2(n1,n2)
  for ts in [TimeSolver2(n1,n2,tensors),
             AfmmSolver2(n1,n2,tensors)]:
    #ts.zeroAt(2*(n1-1)/4,2*(n2-1)/4)
    i1 = [3*n1/8,3*n1/8,5*n1/8,1*n1/8]
    i2 = [2*n2/8,5*n2/8,6*n2/8,6*n2/8]
    for i in range(len(i1)):
      ts.zeroAt(i1[i],i2[i])
    plot(ts.getTimes(),0,0,jet)

def testSeismic():
  n1,n2 = 251,357
  x = readFloats(n1,n2,dataDir+"/seis/tp/tp73.dat")
  #x = gpow(x,0.75)
  plot(x,0,0,gray)
  tensors = SeismicTensors2(x,100,3)
  ts = TimeSolver2(n1,n2,tensors)
  #ts.zeroAt(65,180)
  ts.zeroAt(95,180)
  #ts.zeroAt(118,180)
  #ts.zeroAt( 83,180)
  plot(ts.getTimes(),0,0,prism)

#############################################################################
# other functions

def readFloats(n1,n2,fileName):
  x = Array.zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def gpow(x,gamma):
  return Array.mul(Array.sgn(x),Array.pow(Array.abs(x),gamma))

def coherence(x,sigma):
  n1,n2 = len(x[0]),len(x)
  lof1 = LocalOrientFilter(sigma*2)
  lof2 = LocalOrientFilter(sigma*8)
  u11 = Array.zerofloat(n1,n2)
  u21 = Array.zerofloat(n1,n2)
  su1 = Array.zerofloat(n1,n2)
  sv1 = Array.zerofloat(n1,n2)
  u12 = Array.zerofloat(n1,n2)
  u22 = Array.zerofloat(n1,n2)
  su2 = Array.zerofloat(n1,n2)
  sv2 = Array.zerofloat(n1,n2)
  lof1.apply(x,None,u11,u21,None,None,su1,sv1,None)
  lof2.apply(x,None,u12,u22,None,None,su2,sv2,None)
  c = u11;
  for i2 in range(n2):
    for i1 in range(n1):
      u11i = u11[i2][i1]
      u21i = u21[i2][i1]
      su1i = su1[i2][i1]
      sv1i = sv1[i2][i1]
      u12i = u12[i2][i1]
      u22i = u22[i2][i1]
      su2i = su2[i2][i1]
      sv2i = sv2[i2][i1]
      s111 = (su1i-sv1i)*u11i*u11i+sv1i
      s121 = (su1i-sv1i)*u11i*u21i     
      s221 = (su1i-sv1i)*u21i*u21i+sv1i
      s112 = (su2i-sv2i)*u12i*u12i+sv2i
      s122 = (su2i-sv2i)*u12i*u22i     
      s222 = (su2i-sv2i)*u22i*u22i+sv2i
      s113 = s111*s112+s121*s122
      s223 = s121*s122+s221*s222
      t1 = s111+s221
      t2 = s112+s222
      t3 = s113+s223
      t12 = t1*t2
      if t12>0.0:
        c[i2][i1] = t3/t12
      else:
        c[i2][i1] = 0.0
  return c

#############################################################################
# tensors

class SeismicTensors2(EigenTensors2):
  """
  2D tensors computed from a seismic image
  """
  def __init__(self,x,alpha,sigma):
    EigenTensors2.__init__(self,len(x[0]),len(x))
    n1,n2 = len(x[0]),len(x)
    u1 = Array.zerofloat(n1,n2)
    u2 = Array.zerofloat(n1,n2)
    su = Array.zerofloat(n1,n2)
    sv = Array.zerofloat(n1,n2)
    lof = LocalOrientFilter(sigma)
    lof.apply(x,None,u1,u2,None,None,su,sv,None)
    #c = coherence(x,sigma)
    #c = Array.pow(c,8.0)
    #c = Array.div(1.0,Array.sub(1.0,c))
    #c = Array.div(1.0,Array.pow(Array.sub(1.0,c),0.25))
    c = Array.fillfloat(1.0,n1,n2)
    #plot(c,0,0,jet)
    smax = Array.max(su)
    scale = (alpha*alpha-1.0)/smax
    for i2 in range(n2):
      for i1 in range(n1):
        du = c[i2][i1]/(1.0+scale*su[i2][i1])
        dv = c[i2][i1]/(1.0+scale*sv[i2][i1])
        self.setEigenvalues(i1,i2,du,dv)
        self.setEigenvectorU(i1,i2,u1[i2][i1],u2[i2][i1])

class ConstantTensors2(EigenTensors2):
  """
  2D constant tensors
  dip is angle clockwise between eigenvector v and x2 axis (in degrees)
  du is eigenvalue corresponding to eigenvector u
  dv is eigenvalue corresponding to eigenvector v
  """
  def __init__(self,n1,n2,dip,du,dv):
    EigenTensors2.__init__(self,n1,n2)
    a = dip*pi/180.0
    u1 =  cos(a)
    u2 = -sin(a)
    v1 = -u2
    v2 =  u1
    for i2 in range(n2):
      for i1 in range(n1):
        self.setEigenvalues(i1,i2,du,dv)
        self.setEigenvectorU(i1,i2,u1,u2)

class SineTensors2(EigenTensors2):
  """
  2D seismic-like sine waves
  ampmax is max amplitude of sine waves
  dipmax is maximum dip (in degrees) of sine waves
  """
  def __init__(self,n1,n2,ampmax,dipmax):
    EigenTensors2.__init__(self,n1,n2)
    b1 = 9.0*2.0*pi/(n1-1) # 9 cycles vertically
    b2 = 3.0*2.0*pi/(n2-1) # 3 cycles horizontally
    a1 = ampmax
    a2 = atan(dipmax*pi/180.0)/b2
    f = Array.zerofloat(n1,n2)
    for i2 in range(n2):
      for i1 in range(n1):
        s2 = a2*sin(b2*i2)
        f[i2][i1] = a1*cos(b1*(i1+s2))
        e1 = -a1*b1*sin(b1*(i1+s2))
        e2 = e1*a2*b2*cos(b2*i2)
        den = 1.0+e1*e1+e2*e2
        d11 = (1.0+e2*e2)/den
        d22 = (1.0+e1*e1)/den
        d12 = -e1*e2/den
        self.setTensor(i1,i2,(d11,d12,d22))
    plot(f)

class TsaiTensors2(EigenTensors2):
  """
  2D geodetic tensors for Tsai et al.'s test function.
  """
  def __init__(self,n1,n2):
    EigenTensors2.__init__(self,n1,n2)
    a1 = 2.0*pi
    a2 = 2.0*pi
    d1 = 2.0/(n1-1)
    d2 = 2.0/(n2-1)
    f1 = -1.0
    f2 = -1.0
    f = Array.zerofloat(n1,n2)
    for i2 in range(n2):
      x2 = f2+i2*d2
      for i1 in range(n1):
        x1 = f1+i1*d1
        f[i2][i1] = cos(a1*x1)*sin(a2*x2)
        e1 = -a1*sin(a1*x1)*sin(a2*x2)
        e2 =  a2*cos(a1*x1)*cos(a2*x2)
        den =  1.0+e1*e1+e2*e2
        d11 =  1.0-e1*e1/den
        d22 =  1.0-e2*e2/den
        d12 = -e1*e2/den
        self.setTensor(i1,i2,(d11,d12,d22))
    plot(f)

#############################################################################
# plotting

fontSize = 24
width = 920
height = 900

def plot(f,cmin=0.0,cmax=0.0,cmap=jet,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(cmap)
  frame(p,png)

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(100,6,pngDir+"/"+png+".png")
  return frame

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  return p

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
