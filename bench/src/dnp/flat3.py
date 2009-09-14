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
  #flatten()
  #flattenWithSlopes()
  flattenTest()

n1,n2,n3 = 101,101,101
s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
def flattenTest():
  """Test for t(tau,x) = tau*(1+a*sin(bx)*sin(cy))"""
  x = rampfloat(0,0,0,1,n1,n2,n3)
  y = rampfloat(0,0,1,0,n1,n2,n3)
  t = rampfloat(0,1,0,0,n1,n2,n3)
  smax = 5.0
  a = smax/(n1-1)
  b = 2*PI/(n2-1)
  c = 2*PI/(n3-1)
  bx = mul(b,x)
  cy = mul(c,y)
  cosbx,coscy = cos(bx),cos(cy)
  sinbx,sincy = sin(bx),sin(cy)
  asinbx,asincy = mul(a,sinbx),mul(a,sincy)
  bcosbx,ccoscy = mul(b,cosbx),mul(c,coscy)
  asinbxsincy = mul(asinbx,sincy)
  den = add(1,asinbxsincy)
  p2 = div(mul(t,mul(asinbx,ccoscy)),den)
  p3 = div(mul(t,mul(bcosbx,asincy)),den)
  ep = fillfloat(1,n1,n2,n3)
  fl = FlattenerS(8.0,0.1)
  sf = fl.findShifts(p2,p3,ep) # found shifts
  se = neg(mul(t,asinbxsincy)) # exact shifts
  world = World()
  addImageToWorld(world,sf,jet,-smax,smax)
  addImageToWorld(world,se,jet,-smax,smax)
  makeFrame(world)

def flattenWithSlopes():
  f = readImage("tpst")
  p2 = readImage("tpp2")
  p3 = readImage("tpp3")
  ep = readImage("tpep")
  fl = FlattenerS(8.0,0.1)
  s = fl.findShifts(p2,p3,ep)
  print "s min =",min(s),"max =",max(s)
  #writeImage("tpss",s)
  g = fl.applyShifts(f,s)
  #writeImage("tpsf",g)
  world = World()
  addImage2ToWorld(world,f,g)
  s = clip(-5,5,s)
  addImage2ToWorld(world,g,s)
  makeFrame(world)

def flatten():
  f = readImage(ffile)
  m = readImage(mfile)
  fl = FlattenerS(8.0,0.1)
  s = fl.findShifts(f,m)
  writeImage("tpss",s)
  g = fl.applyShifts(f,s)
  writeImage("tpsf",g)
  world = World()
  addImageToWorld(world,f)
  addImage2ToWorld(world,g,s)
  makeFrame(world)

def slopes():
  f = readImage(ffile)
  m = readImage(mfile)
  p2 = copy(f)
  p3 = copy(f)
  ep = copy(f)
  lsf = LocalSlopeFinder(8.0,5.0)
  lsf.findSlopes(f,p2,p3,ep);
  zm = ZeroMask(m)
  zero = 0.00;
  tiny = 0.01;
  zm.apply(zero,p2);
  zm.apply(zero,p3);
  zm.apply(tiny,ep);
  writeImage("tpp2",p2)
  writeImage("tpp3",p3)
  writeImage("tpep",ep)
  for g in [p2,p3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
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

gray = ColorMap.GRAY
jet = ColorMap.JET

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
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
