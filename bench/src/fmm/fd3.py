#############################################################################
# Fault detection

import sys
from org.python.util import PythonObjectInputStream
from math import *
from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from ldf import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *


#############################################################################
# global parameters

dataDir = "/data/seis/tp/"
#pngDir = "./png"
pngDir = None

n1,n2,n3 = 251,161,357
#n1,n2,n3 = 10,10,10
imageFile = dataDir+"tp3s.dat"
tensorsFile = dataDir+"tpst4.dat"
psemblanceFile = dataDir+"psem4.dat"
psmoothedFile = dataDir+"psmo4.dat"
psubtractFile = dataDir+"psub4.dat"
lsemblanceFile = dataDir+"lsem4.dat"
pvarianceFile = dataDir+"pvar4.dat"
lnpsFile = dataDir+"lnps4.dat"
lnpeFile = dataDir+"lnpe4.dat"
faultFile = dataDir+"tpflt4.dat"
sigma = 4.0
sigma2 = 1.0*sigma
sigma1 = 1.0*sigma
small = 0.01
niter = 100

#############################################################################
# test functions

def main(args):
  #doTensors()
  doPlanarSemblance()
  #doPlanarSmooth()
  #doPlanarSubtract()
  #doPlanarVariance()
  #doLinearSemblance()
  #doLinearNotPlanarSemblance()
  #doLinearNotPlanarEigenvalues()
  #doFaultMask()

#############################################################################
# functions

def doPlanarSemblance():
  f = readImage(n1,n2,n3,imageFile)
  s = semblancePlanar(f)
  writeImage(s,psemblanceFile)
  plot3s([f,s])

def doPlanarSmooth():
  f = readImage(n1,n2,n3,imageFile)
  s = readImage(n1,n2,n3,psemblanceFile)
  print "psemblance min =",Array.min(s)," max =",Array.max(s)
  d = readTensors(tensorsFile)
  g = smoothPlanar(s,d,f)
  writeImage(g,psmoothedFile)
  plot3s([f,g])

def doPlanarSubtract():
  f = readImage(n1,n2,n3,imageFile)
  g = readImage(n1,n2,n3,psmoothedFile)
  h = Array.sub(f,g)
  writeImage(h,psubtractFile)
  plot3s([f,h])

def semblanceIso(f,g):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  rgf = RecursiveGaussianFilter(sigma)
  sf = Array.zerofloat(n1,n2,n3)
  sg = Array.zerofloat(n1,n2,n3)
  rgf.apply000(f,sf)
  rgf.apply000(g,sg)
  Array.sub(f,sf,sf)
  Array.sub(g,sg,sg)
  fg = Array.mul(sf,sg)
  ff = Array.mul(sf,sf)
  gg = Array.mul(sg,sg)
  print "before smoothing"
  print "fg min =",Array.min(fg)," max =",Array.max(fg)
  print "ff min =",Array.min(ff)," max =",Array.max(ff)
  print "gg min =",Array.min(gg)," max =",Array.max(gg)
  rgf.apply000(fg,fg)
  rgf.apply000(ff,ff)
  rgf.apply000(gg,gg)
  print "after smoothing"
  print "fg min =",Array.min(fg)," max =",Array.max(fg)
  print "ff min =",Array.min(ff)," max =",Array.max(ff)
  print "gg min =",Array.min(gg)," max =",Array.max(gg)
  sn = Array.mul(fg,fg)
  sd = Array.mul(ff,gg)
  return Array.div(sn,sd)

def semblancePlanar(f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  scale1 = 0.5*sigma1*sigma1
  scale2 = 0.5*sigma2*sigma2
  lsf1 = LocalSmoothingFilterX(scale1,small,niter)
  lsf2 = LocalSmoothingFilterX(scale2,small,niter)
  g = Array.zerofloat(n1,n2,n3)
  t = readTensors(tensorsFile)
  t.setEigenvalues(0.0,1.0,1.0)
  lsf2.apply(t,f,g)
  h = Array.sub(f,g)
  ff = Array.mul(f,f)
  hh = Array.mul(h,h)
  plot3s([ff,hh])
  sff = Array.zerofloat(n1,n2,n3)
  shh = Array.zerofloat(n1,n2,n3)
  t.setEigenvalues(1.0,0.0,0.0)
  lsf1.apply(t,ff,sff)
  lsf1.apply(t,hh,shh)
  plot3s([sff,shh])
  Array.add(0.0001*Array.max(sff),sff,sff)
  return Array.sub(1.0,Array.div(shh,sff))

def smoothPlanar(s,d,f):
  eu = Array.fillfloat(0.0,n1,n2,n3)
  ev = Array.fillfloat(1.0,n1,n2,n3)
  ew = Array.fillfloat(1.0,n1,n2,n3)
  Array.mul(s,ev,ev)
  Array.mul(s,ew,ew)
  d.setEigenvalues(eu,ev,ew)
  scale2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilterX(scale2,small,niter)
  g = Array.zerofloat(n1,n2,n3)
  lsf.apply(d,f,g)
  return g

def getEigenvalues(d):
  eu = Array.zerofloat(n1,n2,n3)
  ev = Array.zerofloat(n1,n2,n3)
  ew = Array.zerofloat(n1,n2,n3)
  d.getEigenvalues(eu,ev,ew)
  return eu,ev,ew

def doLinearEigenvalues():
  f = readImage(n1,n2,n3,imageFile)
  d = readTensors(tensorsFile)
  eu = Array.zerofloat(n1,n2,n3)
  ev = Array.zerofloat(n1,n2,n3)
  ew = Array.zerofloat(n1,n2,n3)
  d.getEigenvalues(eu,ev,ew)
  es = 0.001*Array.max(eu)
  eu = Array.add(es,eu)
  ep = Array.div(Array.sub(eu,ev),eu)
  el = Array.div(Array.sub(ev,ew),eu)
  ei = Array.div(ew,eu)
  el = Array.mul(el,Array.sub(1.0,ep))
  el = Array.mul(el,Array.sub(1.0,ei))
  writeImage(el,lnpeFile)
  plot3s([f,Array.sub(1.0,el)])

def doLinearSemblance():
  f = readImage(n1,n2,n3,imageFile)
  d = readTensors(tensorsFile)
  ff = Array.zerofloat(n1,n2,n3)
  sn = Array.zerofloat(n1,n2,n3)
  sd = Array.zerofloat(n1,n2,n3)
  d.setEigenvalues(0.0,0.0,1.0)
  scale2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilterX(scale2,small,niter)
  lsf.apply(d,f,sn)
  Array.mul(sn,sn,sn)
  Array.mul(f,f,ff)
  lsf.apply(d,ff,sd)
  ssn = Array.zerofloat(n1,n2,n3)
  ssd = Array.zerofloat(n1,n2,n3)
  d.setEigenvalues(1.0,1.0,0.0)
  scale1 = 0.5*sigma1*sigma1
  lsf = LocalSmoothingFilterX(scale1,small,niter)
  lsf.apply(d,sn,ssn)
  lsf.apply(d,sd,ssd)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply000(ssn,ssn)
  rgf.apply000(ssd,ssd)
  s = Array.div(ssn,ssd)
  plot3s([f,s])
  writeImage(s,lsemblanceFile)

def doFaultMask():
  f = readImage(n1,n2,n3,imageFile)
  s = readImage(n1,n2,n3,semblanceFile)
  g = Array.mul(Array.sub(1.0,s),f)
  plot3s([f,g])
  writeImage(g,faultFile)

def doPlanarVariance():
  f = readImage(n1,n2,n3,imageFile)
  d = readTensors(tensorsFile)
  g = Array.zerofloat(n1,n2,n3)
  vn = Array.zerofloat(n1,n2,n3)
  vd = Array.zerofloat(n1,n2,n3)
  svn = Array.zerofloat(n1,n2,n3)
  svd = Array.zerofloat(n1,n2,n3)
  d.setEigenvalues(0.0,1.0,1.0)
  scale2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilterX(scale2,small,niter)
  lsf.apply(d,f,g)
  Array.sub(f,g,g)
  Array.mul(g,g,vn)
  Array.mul(f,f,vd)
  d.setEigenvalues(1.0,0.0,0.0)
  scale1 = 4.0*0.5*sigma1*sigma1
  lsf = LocalSmoothingFilterX(scale1,small,niter)
  lsf.apply(d,vn,svn)
  lsf.apply(d,vd,svd)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply000(svn,svn)
  rgf.apply000(svd,svd)
  svb = 1.0e-3*Array.max(svd)
  Array.add(svb,svd,svd)
  v = Array.div(svn,svd)
  Array.sub(1.0,v,v)
  plot3s([f,v])
  writeImage(v,pvarianceFile)

def doTensors():
  tensors = computeTensors(imageFile)
  writeTensors(tensors,tensorsFile)
 
def readTensors(tensorsFile):
  fis = FileInputStream(tensorsFile)
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors
 
def writeTensors(tensors,tensorsFile):
  fos = FileOutputStream(tensorsFile)
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

def computeTensors(imageFile):
  f = readImage(n1,n2,n3,imageFile)
  lof = LocalOrientFilter(sigma)
  d = lof.applyForTensors(f)
  return d

def readImage(n1,n2,n3,fileName):
  f = Array.zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(f,fileName):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(f)
  aos.close()

def applyGaussianFilter(sigma,f):
  g = Array.copy(f)
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply000(f,g)
  return g

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

def plot3(f,clip=0.0):
  f = Array.copy(f)
  world = World()
  ipg = ImagePanelGroup(f)
  if clip!=0.0:
    ipg.setClips(-clip,clip)
  else:
    ipg.setPercentiles(0.0,100.0)
  clipMin = ipg.getClipMin()
  clipMax = ipg.getClipMax()
  print "clip min =",clipMin,"max =",clipMax
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def plot3s(flist):
  world = World()
  for f in flist:
    ipg = ImagePanelGroup(f)
    ipg.setPercentiles(1.0,99.0)
    #ipg.setClips(0.0,1.0)
    world.addChild(ipg)
    clipMin = ipg.getClipMin()
    clipMax = ipg.getClipMax()
    print "clip min =",clipMin,"max =",clipMax
  frame = TestFrame(world)
  frame.setVisible(True)

def plot(s1,s2,s3,f,contour):
  tg = TriangleGroup(contour.i,contour.x,contour.u)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.CYAN)
  states.add(cs)
  lms = LightModelState()
  lms.setTwoSide(True)
  states.add(lms)
  ms = MaterialState()
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
  ms.setSpecular(Color.WHITE)
  ms.setShininess(100.0)
  states.add(ms)
  tg.setStates(states);
  ipg = ImagePanelGroup(s1,s2,s3,f)
  ipg.setColorModel(jet)
  world = World()
  world.addChild(tg)
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
