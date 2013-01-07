"""
Jython utilities for F3 dataset.
Author: Dave Hale, Colorado School of Mines
Version: 2012.12.30
"""
from imports import *

#############################################################################
# Internal constants

_f3dDir = "/data/seis/f3d/"

#############################################################################
# Setup

seismicDir = _f3dDir
wellLogsDir = _f3dDir
s1 = Sampling(462,0.004,0.004)
s2 = Sampling(951,0.025,0.000)
s3 = Sampling(591,0.025,0.000)
def setupForSubset(name):
  global s1,s2,s3
  if name=="all":
    s1 = Sampling(462,0.004,0.004)
    s2 = Sampling(951,0.025,0.000)
    s3 = Sampling(651,0.025,0.000)

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

#############################################################################
# read/write files

def readImage(name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "f3d"
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
  name: base name of image file; e.g., "f3gp"
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
  Reads tensors from file with specified basename; e.g., "f3et".
  """
  fis = FileInputStream(seismicDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors
def writeTensors(name,tensors):
  """
  Writes tensors to file with specified basename; e.g., "f3et".
  """
  fos = FileOutputStream(seismicDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

def readLogSamples(type,smooth=0):
  """ 
  Reads log curves from the specified set that have the specified type.
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  fileName = wellLogsDir+"f3dwell.dat"
  wdata = WellLog.Data.readBinary(fileName)
  logs = wdata.getLogsWith(type)
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    print log.name
    if smooth: 
      log.smooth(smooth)
    #samples = log.getSamples(type,s1,s2,s3) # limit to seismic grid
    samples = log.getSamples(type,None,None,None) # include all logs
    if samples:
      f,x1,x2,x3 = samples
      fl.append(f)
      x1l.append(x1)
      x2l.append(x2)
      x3l.append(x3)
  return fl,x1l,x2l,x3l

def readLogSamplesMerged(set,type,smooth=0):
  """ 
  Same as readLogSamples, except log sample values for all logs are
  merged into one array, so that this function returns four arrays 
  (f,x1,x2,x3) instead of four lists of arrays.
  """
  fl,x1l,x2l,x3l = readLogSamples(set,type,smooth)
  n = 0
  for f in fl:
    n += len(f)
  f = zerofloat(n)
  x1 = zerofloat(n)
  x2 = zerofloat(n)
  x3 = zerofloat(n)
  n = 0
  for i,fi in enumerate(fl):
    x1i,x2i,x3i = x1l[i],x2l[i],x3l[i]
    ni = len(fi)
    copy(ni,0,1,fi,n,1,f)
    copy(ni,0,1,x1i,n,1,x1)
    copy(ni,0,1,x2i,n,1,x2)
    copy(ni,0,1,x3i,n,1,x3)
    n += ni
  return f,x1,x2,x3

def getWellIntersections(type,x1):
  fileName = wellLogsDir+"welli.dat"
  wdata = WellLog.Data.readBinary(fileName)
  x2,x3 = wdata.getIntersections(type,x1)
  return x2,x3

#############################################################################
# graphics

def addImageToWorld(world,image,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2,cmin=0,cmax=0):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  if cmin<cmax:
    ipg.setClips1(cmin,cmax)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.setAlpha(ColorMap.PRISM,0.7))
  world.addChild(ipg)
  return ipg

def addTensorsInImage(ip,et,esize):
  tp = TensorsPanel(s1,s2,s3,et)
  tp.setEllipsoidSize(esize)
  ip.getFrame().addChild(tp)
  return tp

def addLogsToWorld(world,type,cmin=0,cmax=0,cbar=None,smooth=0):
  samples = readLogSamples(type,smooth)
  print "number of logs =",len(samples[0])
  lg = makeLogPoints(samples,type,cmin,cmax,cbar)
  #lg = makeLogLines(samples,type,cmin,cmax)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.YELLOW)
  states.add(cs)
  lg.setStates(states)
  world.addChild(lg)

def makeLogPoints(samples,type,cmin,cmax,cbar):
  lg = Group()
  fl,x1l,x2l,x3l = samples
  for i,f in enumerate(fl):
    f = fl[i]
    x1 = x1l[i]
    x2 = x2l[i]
    x3 = x3l[i]
    pg = makePointGroup(f,x1,x2,x3,cmin,cmax,cbar)
    lg.addChild(pg)
  return lg

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.JET)
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(4)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

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
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
