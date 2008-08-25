#############################################################################
# Tests MarchingCubes

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

from fmm import *

#############################################################################
# global parameters

dataDir = "/data/seis/tp/"
#pngDir = "./png"
pngDir = None

n1,n2,n3 = 251,161,357
#n1,n2,n3 = 10,10,10

#############################################################################
# test functions

def main(args):
  imageFile = dataDir+"tp3s.dat"
  tensorsFile = dataDir+"et3s211.dat"
  makeTensors(imageFile,tensorsFile)

#############################################################################
# tensors functions

def makeTensors(imageFile,tensorsFile):
  tensors = computeTensors(imageFile)
  writeTensors(tensors,tensorsFile)
 
def readTensors(tensorsFile):
  fis = FileInputStream(ptFile)
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
  alpha,beta,gamma = 2.0,1.0,1.0
  f = readImage(n1,n2,n3,imageFile)
  u1 = Array.zerofloat(n1,n2,n3)
  u2 = Array.zerofloat(n1,n2,n3)
  w1 = Array.zerofloat(n1,n2,n3)
  w2 = Array.zerofloat(n1,n2,n3)
  su = Array.zerofloat(n1,n2,n3)
  sv = Array.zerofloat(n1,n2,n3)
  sw = Array.zerofloat(n1,n2,n3)
  lof = LocalOrientFilter(3.0)
  lof.apply(f,
    None,None,
    u1,u2,None,
    None,None,None,
    w1,w2,None,
    su,sv,sw,
    None,None)
  sumax = Array.max(su)
  sbias = sumax*0.01
  su = Array.add(su,sbias)
  sv = Array.add(sv,sbias)
  sw = Array.add(sw,sbias)
  print "max su =",Array.max(su)," sv =",Array.max(sv)," sw =",Array.max(sw)
  #dc = Array.pow(Array.sub(1.0,coherence(sigma,x)),-gamma)
  dc = Array.fillfloat(1.0,n1,n2,n3)
  du = Array.mul(dc,Array.pow(su,-alpha))
  dv = Array.mul(du,Array.pow(Array.div(sv,su),-beta))
  dw = Array.mul(du,Array.pow(Array.div(sw,su),-beta))
  print "max du =",Array.max(du)," dv =",Array.max(dv)," dw =",Array.max(dw)
  ds = 1.0/Array.max(dw)
  du = Array.mul(ds,du)
  dv = Array.mul(ds,dv)
  dw = Array.mul(ds,dw)
  print "ds =",ds
  print "max du =",Array.max(du)," dv =",Array.max(dv)," dw =",Array.max(dw)
  return EigenTensors3(u1,u2,w1,w2,du,dv,dw,True) 

def readImage(n1,n2,n3,fileName):
  f = Array.zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(f)
  ais.close()
  return f

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

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
  ipg = ImagePanelGroup(s3,s2,s1,SimpleFloat3(f))
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
