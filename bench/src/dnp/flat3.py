import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from dnp import *

seismicDir = "/data/seis/tp/csm/seismict/subt_251_4_500/"
ffile = "tpst"
mfile = "tpmt"
s1 = Sampling(251,0.004,0.500)
s2 = Sampling(357,0.025,0.000)
s3 = Sampling(161,0.025,0.000)
n1,n2,n3 = s1.count,s2.count,s3.count

def main(args):
  #display()
  #slopes()
  normals()

def normals():
  f = readImage(ffile)
  m = readImage(mfile)
  u2 = copy(f)
  u3 = copy(f)
  ep = copy(f)
  Flattener.normals(f,u2,u3,ep)
  mask = ZeroMask(m)
  mask.apply(0.0,u2)
  mask.apply(0.0,u3)
  mask.apply(0.0,ep)
  for g in [u2,u3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
    makeFrame(world)

def slopes():
  f = readImage(ffile)
  m = readImage(mfile)
  p2 = copy(f)
  p3 = copy(f)
  Flattener.slopes(1,1,f,p2,p3)
  mask = ZeroMask(m)
  mask.apply(0.0,p2)
  mask.apply(0.0,p3)
  world = World()
  addImage2ToWorld(world,f,p2)
  makeFrame(world)
  world = World()
  addImage2ToWorld(world,f,p3)
  makeFrame(world)

def display():
  f = readImage(ffile)
  world = World()
  addImageToWorld(world,f)
  makeFrame(world)

#############################################################################
# read/write files

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

def readSlice3(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

from org.python.util import PythonObjectInputStream
def readTensors(name):
  """
  Reads tensors from file with specified basename; e.g., "tpet".
  """
  fis = FileInputStream(seismicDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors
def writeTensors(name,tensors):
  """
  Writes tensors to file with specified basename; e.g., "tpet".
  """
  fos = FileOutputStream(seismicDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

#############################################################################
# graphics

def addImageToWorld(world,image):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  ipg.setColorModel2(ColorMap.getJet(0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(2.0)
  view.setAzimuth(-50.0)
  frame.setSize(1200,900)
  frame.setVisible(True)
  return frame

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
