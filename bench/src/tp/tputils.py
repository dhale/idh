"""
Jython utilities for Teapot Dome.
Author: Dave Hale, Colorado School of Mines
Version: 2009.07.25
"""
from imports import *

#############################################################################
# Internal constants

_tpDir = "/data/seis/tpd/"
_csmDir = _tpDir+"csm/"
_horizontDir = _csmDir+"horizont/"
_horizonzDir = _csmDir+"horizonz/"
_seismictDir = _csmDir+"seismict/"
_seismiczDir = _csmDir+"seismicz/"
_horizonColors = {
  "KF2F2WC":Color.RED,
  "FallRiverDKOT":Color.GREEN,
  "CrowMountainCRMT":Color.CYAN,
  "TensleepASand":Color.YELLOW,
  "TensleepBbaseC1Dolo":Color.MAGENTA,
  #"BasementPC":Color.BLUE
}
_horizonNames = _horizonColors.keys()

#############################################################################
# Setup

def setupForSubset(name):
  """
  Setup for a specified directory includes:
    time or depth flags
    seismic, horizon and well log directories
    samplings s1,s2,s3
  Example: setupForSubset("subt_251_4_500")
  """
  global _time,_depth,_tz
  global seismicDir,horizonDir,wellLogsDir
  global s1,s2,s3
  subdir = name.split("/")[-1] # remove any base directory names
  tz = subdir[3] # subdirectory name begins with "subt" or "subz"
  if tz=="t":
    _time = True; _depth = False; _tz = "t"
    horizonDir = _horizontDir
    seismicDir = _seismictDir+subdir+"/"
  else:
    _time = False; _depth = True; _tz = "z"
    horizonDir = _horizonzDir
    seismicDir = _seismiczDir+subdir+"/"
  wellLogsDir = _csmDir+"welllogs/"
  parts = subdir.split("_")
  n1 = int(parts[1])
  d1 = float(parts[2])*0.001
  f1 = float(parts[3])*0.001
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(357,0.025,0.000)
  s3 = Sampling(161,0.025,0.000)

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

def readHorizon(name):
  """ 
  Reads a horizon with specified name.
  name: horizon name; e.g., "CrowMountainCRMT" 
  Returns a horizon.
  """
  fileName = horizonDir+"tph"+_tz+name+".dat"
  return Horizon.readBinary(fileName)

def readLogSamples(set,type,smooth=0):
  """ 
  Reads log curves from the specified set that have the specified type.
  set: "s" for shallow, "d" for deep, or "a" for all
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  fileName = wellLogsDir+"tpw"+set[0]+".dat"
  wdata = WellLog.Data.readBinary(fileName)
  logs = wdata.getLogsWith(type)
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    if smooth: 
      log.smooth(smooth)
    samples = log.getSamples(type,s1,s2,s3)
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

def getWellIntersections(set,type,x1):
  fileName = wellLogsDir+"tpw"+set[0]+".dat"
  wdata = WellLog.Data.readBinary(fileName)
  x2,x3 = wdata.getIntersections(type,x1)
  return x2,x3

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

def addLogsToWorld(world,set,type,cmin=0,cmax=0,cbar=None,smooth=0):
  samples = readLogSamples(set,type,smooth)
  #print "number of logs =",len(samples[0])
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
