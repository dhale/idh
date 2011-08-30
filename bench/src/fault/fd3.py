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

from fault import FaultFinder3
from dnp import LocalSlopeFinder

#############################################################################

def main(args):
  goScan()
  #goShear()
  #goAlign()
  #goRotate()
  #goSlopes()

def goScan():
  ff = makeFaultFinder()
  s1,s2,s3,f = imageF3d()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  f = slog(f)
  p2,p3 = ff.findSlopes(f)
  f = ff.taper(10,f)
  #phi,theta = 90,-8
  #phi,theta = -50,-8
  phi,theta = 50,-8
  sp = ff.makePhiSampling(phi,phi)
  st = ff.makeThetaSampling(theta,theta)
  ff.setPhiSampling(sp)
  ff.setThetaSampling(st)
  c,p,t = ff.scan(p2,p3,f)
  plot3(f)
  plot3(f,c,0,1)

def goAlign():
  ff = makeFaultFinder()
  s1,s2,s3,f = imageF3d()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  f = slog(f)
  f = ff.taper(10,f)
  plot3(f)
  #p2 = p3 = zerofloat(n1,n2,n3)
  p2,p3 = ff.findSlopes(f)
  phi = 50
  r = FaultFinder.Rotator(phi,n1,n2,n3)
  f = r.rotate(f)
  p2 = r.rotate(p2)
  p3 = r.rotate(p3)
  fm,fp = ff.align(phi,p2,p3,f)
  nullsToZeros(fm)
  nullsToZeros(fp)
  fe = sub(fp,fm)
  print "rms error = ",sqrt(sum(mul(fe,fe))/n1/n2/n3)
  plot3(f)
  plot3(fm)
  plot3(fp)
  plot3(fe)

def goShear():
  ff = makeFaultFinder()
  s1,s2,s3,f = imageF3d()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  f = slog(f)
  shear = -0.13
  g = ff.shear(shear,f)
  h = ff.unshear(shear,g)
  e = sub(h,f)
  print "shear-unshear error max =",max(abs(e))
  plot3(f)
  plot3(g)
  plot3(h)
  plot3(e)

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

def goRotate():
  ff = makeFaultFinder()
  s1,s2,s3,f = imageF3d()
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  f = slog(f)
  f = ff.taper(10,f)
  plot3(f)
  #r = FaultFinder.Rotator(-50.0,n1,n2,n3)
  r = FaultFinder.Rotator(-90.0,n1,n2,n3)
  g = r.rotate(f)
  h = r.unrotate(g)
  plot3(g)
  plot3(h)

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

def makeFaultFinder():
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 20.0
  return FaultFinder3(slopeMax,shiftMax,thetaMax)
 
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

def nullsToZeros(f):
  FaultFinder3.nullsToZeros(f)

def plot3(f,g=None,gmin=None,gmax=None):
  nullsToZeros(f)
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame()
  if g==None:
    ipg = sf.addImagePanels(f)
  else:
    nullsToZeros(g)
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
