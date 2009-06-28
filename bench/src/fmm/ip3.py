#############################################################################
# Image painting in 3D.

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
from edu.mines.jtk.util.ArrayMath import *

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
  #tensorsFile = dataDir+"tp_et211.dat"
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
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  w1 = zerofloat(n1,n2,n3)
  w2 = zerofloat(n1,n2,n3)
  su = zerofloat(n1,n2,n3)
  sv = zerofloat(n1,n2,n3)
  sw = zerofloat(n1,n2,n3)
  lof = LocalOrientFilter(3.0)
  lof.apply(f,
    None,None,
    u1,u2,None,
    None,None,None,
    w1,w2,None,
    su,sv,sw,
    None,None)
  sumax = max(su)
  sbias = sumax*0.01
  su = add(su,sbias)
  sv = add(sv,sbias)
  sw = add(sw,sbias)
  print "max su =",max(su)," sv =",max(sv)," sw =",max(sw)
  #ps = readImage(n1,n2,n3,dataDir+"tp_ps.dat")
  #print "ps min =",min(ps)," max =",max(ps)
  #dc = pow(sub(1.01,ps),-gamma)
  dc = fillfloat(1.0,n1,n2,n3)
  du = mul(dc,pow(su,-alpha))
  dv = mul(du,pow(div(sv,su),-beta))
  dw = mul(du,pow(div(sw,su),-beta))
  print "min du =",min(du)," dv =",min(dv)," dw =",min(dw)
  print "max du =",max(du)," dv =",max(dv)," dw =",max(dw)
  ds = 1.0/max(dw)
  du = mul(ds,du)
  dv = mul(ds,dv)
  dw = mul(ds,dw)
  print "ds =",ds
  print "min du =",min(du)," dv =",min(dv)," dw =",min(dw)
  print "max du =",max(du)," dv =",max(dv)," dw =",max(dw)
  return EigenTensors3(u1,u2,w1,w2,du,dv,dw,True) 

def readImage(n1,n2,n3,fileName):
  f = zerofloat(n1,n2,n3)
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
