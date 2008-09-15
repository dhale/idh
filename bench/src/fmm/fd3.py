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
varianceFile = dataDir+"tpvar4.dat"
semblanceFile = dataDir+"tpsem4.dat"
sigma = 4.0
sigma2 = sigma
sigma1 = sigma
small = 0.01
niter = 100

#############################################################################
# test functions

def main(args):
  doTensors()
  #doVariance()
  #doSemblance()

#############################################################################
# functions

def doSemblance():
  f = readImage(n1,n2,n3,imageFile)
  d = readTensors(tensorsFile)
  ff = Array.zerofloat(n1,n2,n3)
  sn = Array.zerofloat(n1,n2,n3)
  sd = Array.zerofloat(n1,n2,n3)
  d.setEigenvalues(0.0,1.0,1.0)
  scale2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilter(scale2,small,niter)
  lsf.apply(d,f,sn)
  Array.mul(sn,sn,sn)
  Array.mul(f,f,ff)
  lsf.apply(d,ff,sd)
  ssn = Array.zerofloat(n1,n2,n3)
  ssd = Array.zerofloat(n1,n2,n3)
  d.setEigenvalues(1.0,0.0,0.0)
  scale1 = 0.5*sigma1*sigma1
  lsf = LocalSmoothingFilter(scale1,small,niter)
  lsf.apply(d,sn,ssn)
  lsf.apply(d,sd,ssd)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply000(ssn,ssn)
  rgf.apply000(ssd,ssd)
  s = Array.div(ssn,ssd)
  plot3s([f,s])
  writeImage(s,semblanceFile)

def doVariance():
  f = readImage(n1,n2,n3,imageFile)
  d = readTensors(tensorsFile)
  g = Array.zerofloat(n1,n2,n3)
  vn = Array.zerofloat(n1,n2,n3)
  vd = Array.zerofloat(n1,n2,n3)
  svn = Array.zerofloat(n1,n2,n3)
  svd = Array.zerofloat(n1,n2,n3)
  d.setEigenvalues(0.0,1.0,1.0)
  scale2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilter(scale2,small,niter)
  lsf.apply(d,f,g)
  Array.sub(f,g,g)
  Array.mul(g,g,vn)
  Array.mul(f,f,vd)
  d.setEigenvalues(1.0,0.0,0.0)
  scale1 = 4.0*0.5*sigma1*sigma1
  lsf = LocalSmoothingFilter(scale1,small,niter)
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
  writeImage(v,varianceFile)

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
    ipg.setPercentiles(1.0,99.0)
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
