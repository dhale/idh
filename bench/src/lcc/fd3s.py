#############################################################################
# Fault displacements from 3D images

import sys
from org.python.util import PythonObjectInputStream
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from lcc import FaultFinder3SB
from dnp import LocalSlopeFinder

#############################################################################

def makeFaultFinder():
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 20.0
  return FaultFinder3SB(slopeMax,shiftMax,thetaMax)

def main(args):
  goScan()

def goScan():
  ff = makeFaultFinder()
  s1,s2,s3,f = imageF3d()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  f = slog(f)
  #plot3(f)
  p = ff.findSlopes(f)
  snd = ff.semblanceNumDen(p,f)
  #f = ff.taper(10,f)
  #phi,theta = 90,-8
  #phi,theta = -50,-8
  #phi,theta = 50,-8
  phi,theta = 0.0,0.0
  sp = ff.makePhiSampling(phi,phi)
  st = ff.makeThetaSampling(theta,theta)
  ff.setPhiSampling(sp)
  ff.setThetaSampling(st)
  print "scanning ..."
  c,p,t = ff.faultPhiThetaScan(snd)
  print "c min =",min(c)," max =",max(c)
  plot3(f,c,0,1)

def goSlopes():
  ff = makeFaultFinder()
  s1,s2,s3,f = imageF3d()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  f = slog(f)
  f = ff.taper(10,f)
  plot3(f)
  p2,p3 = ff.findSlopes(f)
  p2 = clip(-1,1,p2)
  p3 = clip(-1,1,p3)
  plot3(f,p2)
  plot3(f,p3)

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))
 
#############################################################################
# data read/write

def readImage(n1,n2,n3,fileName):
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(f,fileName):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(f)
  aos.close()

def imageF3d():
  n1,n2,n3 = 462,951,591
  #j1,j2,j3 = 240,0,0
  #m1,m2,m3 = 222,440,440
  j1,j2,j3 = 240, 50,100
  m1,m2,m3 = 222,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.004,0.000,0.000
  firstTime = False
  if firstTime:
    fileName = "/data/seis/f3d/f3d.dat"
    af = ArrayFile(fileName,"r")
    x = zerofloat(m1,m2,m3)
    for i3 in range(m3):
      for i2 in range(m2):
        af.seek(4*(j1+n1*(i2+j2+n2*(i3+j3))))
        af.readFloats(x[i3][i2])
    af.close()
    #writeImage(x,"/data/seis/f3d/f3ds1.dat")
  n1,n2,n3 = m1,m2,m3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  x = readImage(m1,m2,m3,"/data/seis/f3d/f3ds1.dat")
  return s1,s2,s3,x

#############################################################################
# plotting

def plot3(f,g=None,gmin=None,gmax=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame()
  if g==None:
    ipg = sf.addImagePanels(f)
  else:
    ipg = ImagePanelGroup2(f,g)
    sf.world.addChild(ipg)
  sf.setSize(1200,900)
  sf.orbitView.setScale(3.0)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
