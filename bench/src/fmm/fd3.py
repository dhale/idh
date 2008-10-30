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
tensorsFile = dataDir+"tp_st8.dat" # structure tensors
fFile = dataDir+"tp_f.dat" # input image
gFile = dataDir+"tp_pg.dat" # planar-smoothed image
ffFile = dataDir+"tp_pff.dat" # square of planar-smoothed image
ggFile = dataDir+"tp_pgg.dat" # planar-smoothed image-squared
sffFile = dataDir+"tp_psff.dat" # numerator of planar semblance
sggFile = dataDir+"tp_psgg.dat" # denominator of planar semblance
psFile = dataDir+"tp_ps.dat" # planar semblance
tsFile = dataDir+"tp_ts.dat" # tensor semblance
sigma = 3.0
sigma2 = 1.0*sigma
sigma1 = 2.0*sigma
small = 0.01
niter = 100

#############################################################################
# test functions

def main(args):
  #makeTensors()
  #makePlanarSemblance()
  #plotPlanarSemblance()
  #plotPlanarSemblanceLimits()
  makeTensorSemblance()
  #plotPlanarTensorSemblance()

#############################################################################
# functions

def plotPlanarTensorSemblance():
  f = readImage(n1,n2,n3,fFile)
  ps = readImage(n1,n2,n3,psFile)
  ts = readImage(n1,n2,n3,tsFile)
  plot3s([f,ps])
  plot3s([f,ts])
  plot3s([ps,ts])

def makeTensorSemblance():
  sigma1,sigma2 = 2.0,8.0
  f = readImage(n1,n2,n3,fFile)
  s = SemblanceFilter.applyTensorTraces(sigma1,sigma2,f)
  writeImage(s,tsFile)
  plot3s([f,s])

def makePlanarSemblance():
  f = readImage(n1,n2,n3,fFile)
  g = Array.zerofloat(n1,n2,n3)
  scale1 = 0.5*sigma1*sigma1
  scale2 = 0.5*sigma2*sigma2
  lsf1 = LocalSmoothingFilterX(scale1,small,niter)
  lsf2 = LocalSmoothingFilterX(scale2,small,niter)
  t = readTensors(tensorsFile)
  t.setEigenvalues(0.0,1.0,1.0)
  lsf2.apply(t,f,g)
  writeImage(g,gFile)
  ff = Array.mul(f,f)
  gg = Array.mul(g,g)
  lsf2.apply(t,Array.copy(ff),ff)
  writeImage(ff,ffFile)
  writeImage(gg,ggFile)
  sff = Array.zerofloat(n1,n2,n3)
  sgg = Array.zerofloat(n1,n2,n3)
  t.setEigenvalues(1.0,0.0,0.0)
  lsf1.apply(t,ff,sff)
  lsf1.apply(t,gg,sgg)
  writeImage(sff,sffFile)
  writeImage(sgg,sggFile)
  ps = Array.div(sgg,sff)
  ps = Array.clip(0.0,1.0,ps)
  writeImage(ps,psFile)

def plotPlanarSemblance():
  f = readImage(n1,n2,n3,fFile)
  g = readImage(n1,n2,n3,gFile)
  ff = readImage(n1,n2,n3,ffFile)
  gg = readImage(n1,n2,n3,ggFile)
  sff = readImage(n1,n2,n3,sffFile)
  sgg = readImage(n1,n2,n3,sggFile)
  ps = readImage(n1,n2,n3,psFile)
  plot3s([f,g])
  plot3s([ff,gg])
  plot3s([sff,sgg])
  plot3s([f,ps])

def plotPlanarSemblanceLimits():
  f = readImage(n1,n2,n3,fFile)
  sff = readImage(n1,n2,n3,sffFile)
  sgg = readImage(n1,n2,n3,sggFile)
  ps = Array.div(sgg,sff)
  pslo = Array.clip(-0.01,0.01,ps)
  pshi = Array.clip( 0.99,1.01,ps)
  print "pslo min =",Array.min(pslo)," max =",Array.max(pslo)
  print "pshi min =",Array.min(pshi)," max =",Array.max(pshi)
  plot3s([f,pslo])
  plot3s([f,pshi])

def makeTensors():
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

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
