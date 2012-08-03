"""
Jython utilities for Canning 3D.
Author: Dave Hale, Colorado School of Mines
Version: 2012.08.03
"""
from imports import *

#############################################################################
# Internal constants

_datdir = "/data/dhale/can/"

#############################################################################
# Setup

def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2,s3
  Example: setupForSubset("canf")
  """
  global seismicDir
  global s1,s2,s3
  if name=="s1":
    """ subset good for faults """
    seismicDir = _datdir+"s1/"
    n1,n2,n3 = 501,1601,1001
    d1,d2,d3 = 1.0,1.0,1.0
    f1,f2,f3 = 0.0,0.0,0.0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  else:
    print "unrecognized subset:",name
    System.exit

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

#############################################################################
# read/write files

def readImage(name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
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
  ois.close()
  return tensors
def writeTensors(name,tensors):
  """
  Writes tensors to file with specified basename; e.g., "tpet".
  """
  fos = FileOutputStream(seismicDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  oos.close()

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

def addTensorsInImage(ip,et,esize):
  tp = TensorsPanel(s1,s2,s3,et)
  tp.setEllipsoidSize(esize)
  ip.getFrame().addChild(tp)
  return tp

def getHorizonColor(name):
  return _horizonColors[name]

def makeHorizonTriangles(horizon,color):
  ijk = horizon.getIABC()
  xyz = horizon.getX321()
  tg = TriangleGroup(ijk,xyz)
  tg.setColor(color)
  return tg

def addHorizonToWorld(world,name):
  horizon = readHorizon(name)
  color = getHorizonColor(name)
  tg = makeHorizonTriangles(horizon,color)
  world.addChild(tg)
  return tg

def addAllHorizonsToWorld(world):
  for name in _horizonNames:
    addHorizonToWorld(world,name)

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1250,900)
  frame.setVisible(True)
  return frame

#############################################################################
# other

def addNoise(s,x):
  """ adds s percent random noise to reduce problems with zero traces """
  xdif = max(x)-min(x)
  n1,n2,n3 = s1.count,s2.count,s3.count
  r = randfloat(n1,n2,n3)
  x = add(mul(sub(r,0.5),s*xdif),x)
  return x

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
