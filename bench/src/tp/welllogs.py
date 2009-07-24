import sys
from math import *

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from tp import *

# Directories and files for well logs, headers, directional surveys
tpDir = "/data/seis/tp/"
doeWellLogsDir = tpDir+"doe/WellLogs/"
csmWellLogsDir = tpDir+"csm/welllogs/"
doeWellHeaders = doeWellLogsDir+"WellHeaders.txt"
doeDirectionalSurveys = doeWellLogsDir+"DirectionalSurveys.txt"
doeWellLogs = ""
csmWellLogs = ""

# File for 3D seismic depth image to view with wells.
csmSeismiczDir = tpDir+"csm/seismicz/"
csmSeismic = csmSeismiczDir+"tpsz.dat"
s1 = Sampling(2762,0.002,0.000)
s2 = Sampling(357,0.025,0.000)
s3 = Sampling(161,0.025,0.000)

def main(args):
  #makeBinaryWellLogs("test")
  #makeBinaryWellLogs("deep")
  #makeBinaryWellLogs("shallow")
  #makeBinaryWellLogs("all")
  #viewWellCoordinates("deep")
  #viewWellCurves("deep","velocity")
  #viewWellCurves("deep","density")
  #viewWellCurves("deep","gamma")
  #viewWellCurves("deep","porosity")
  #viewWellsWithSeismic("all","velocity")
  #viewWellsWithSeismic("deep","velocity")
  viewWellsWithSeismic("all","gamma")
  #viewWellsWithSeismic("deep","gamma")
  viewElevations("deep")

def setGlobals(what):
  global csmWellLogs,doeWellLogs
  if what=="all":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/"
    csmWellLogs = csmWellLogsDir+"tpwa.dat"
  elif what=="deep":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/Deeper_LAS_files/"
    csmWellLogs = csmWellLogsDir+"tpwd.dat"
  elif what=="shallow":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/Shallow_LAS_files/"
    csmWellLogs = csmWellLogsDir+"tpws.dat"
  elif what=="test":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/test/"
    csmWellLogs = csmWellLogsDir+"tpwt.dat"

# Reads well locations in CSM (x1,x2,x3) coordinates
def readLocations(what):
  setGlobals(what)
  wldata = WellLog.Data.readBinary(csmWellLogs)
  wldata.printInfo()
  nlog = wldata.size()
  x1 = zerofloat(nlog)
  x2 = zerofloat(nlog)
  x3 = zerofloat(nlog)
  ilog = 0
  for log in wldata.getAll():
    map = Coordinates.Map(log.xe,log.yn,log.ze)
    csm = Coordinates.Csm(map)
    x1[ilog] = csm.x1
    x2[ilog] = csm.x2
    x3[ilog] = csm.x3
    ilog += 1
  return x1,x2,x3

# Grids z(x,y) for specified samplings.
def gridData(x,y,z,sx,sy):
  si = SibsonInterpolator2(z,x,y)
  #si.setGradientPower(1.0)
  si.setBounds(sx,sy)
  return si.interpolate(sx,sy)

# Returns estimated errors for interpolated data.
def errors(x,y,z,sx,sy):
  si = SibsonInterpolator2(z,x,y)
  si.setBounds(sx,sy)
  n = len(x)
  v = zerofloat(n)
  for i in range(n):
    v[i] = si.validate(i)
  return sub(z,v)

def viewElevations(what):
  x1,x2,x3 = readLocations(what)
  x,y,z = x2,x3,x1
  sx,sy = s2,s3
  zi = gridData(x,y,z,sx,sy)
  e = errors(x,y,z,sx,sy)
  ei = gridData(x,y,e,sx,sy)
  zi = mul(1000.0,zi)
  ei = mul(1000.0,ei)
  viewPointsPixels(x,y,z,sx,sy,zi,"Elevation (m)")
  viewPointsPixels(x,y,e,sx,sy,ei,"Validation error (m)")

def viewPointsPixels(x,y,z,sx,sy,zs,zlabel):
  sp = SimplePlot()
  sp.setHLabel("X (km)")
  sp.setVLabel("Y (km)")
  pp = sp.getPlotPanel()
  pp.setColorBarWidthMinimum(100)
  pv = sp.addPixels(sx,sy,zs)
  pv.setColorModel(ColorMap.JET)
  pv = sp.addPoints(x,y)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(3)
  sp.addColorBar(zlabel)
  sp.setSize(970,440)

def makeBinaryWellLogs(what):
  setGlobals(what)
  wldata = WellLog.Data(doeWellLogs,doeWellHeaders,doeDirectionalSurveys)
  print "well log data"
  wldata.printInfo()
  wldata.clip(s2,s3)
  print "after clipping"
  wldata.printInfo()
  wldata.writeBinary(csmWellLogs)

def dumpWellHeaders(what):
  setGlobals(what)
  whdata = WellHeader.Data(doeWellHeaders)
  whdata.printInfo()

def dumpDirectionalSurveys(what):
  setGlobals(what)
  dsdata = DirectionalSurvey.Data(doeDirectionalSurveys)
  dsdata.printInfo()

def viewWellCoordinates(what):
  setGlobals(what)
  wldata = WellLog.Data.readBinary(csmWellLogs)
  wldata.printInfo()
  spz1 = SimplePlot()
  sp23 = SimplePlot()
  zlist = []
  x1list = []
  x2list = []
  x3list = []
  for log in wldata.getAll():
    zlist.append(log.z)
    x1list.append(log.dx1dz())
    x2list.append(log.x2)
    x3list.append(log.x3)
  tv = PointsView(zlist,x1list)
  spz1.add(tv)
  tv = PointsView(x2list,x3list)
  sp23.add(tv)

def viewWellCurves(what,curve):
  setGlobals(what)
  wldata = WellLog.Data.readBinary(csmWellLogs)
  zlist,clist = [],[]
  for log in wldata.getLogsWith(curve):
    if curve=="velocity":
      z,c = getNonNullValues(log.x1,log.v)
    elif curve=="density":
      z,c = getNonNullValues(log.x1,log.d)
    elif curve=="gamma":
      z,c = getNonNullValues(log.x1,log.g)
    elif curve=="porosity":
      z,c = getNonNullValues(log.x1,log.p)
    zlist.append(z)
    clist.append(c)
  sp = SimplePlot()
  sp.setHLabel("Depth (km)")
  if curve=="velocity":
    sp.setVLabel("Velocity (km/s)")
  elif curve=="density":
    sp.setVLabel("Density (kg/m3)")
  elif curve=="gamma":
    sp.setVLabel("Gamma Ray")
  elif curve=="porosity":
    sp.setVLabel("Porosity")
  tv = PointsView(zlist,clist)
  sp.add(tv)

def getNonNullValues(z,c):
  zlist = FloatList()
  clist = FloatList()
  cnull = -999.2500
  n = len(z)
  for i in range(n):
    if c[i]!=cnull:
      zlist.add(z[i])
      clist.add(c[i])
  return zlist.trim(),clist.trim()

def viewWellsWithSeismic(what,curve):
  setGlobals(what)
  wdata = WellLog.Data.readBinary(csmWellLogs)
  n1c,n2c,n3c = s1.count,s2.count,s3.count
  ais = ArrayInputStream(csmSeismic)
  x = zerofloat(n1c,n2c,n3c)
  ais.readFloats(x)
  ais.close()
  print "x min =",min(x)," max =",max(x)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  world = World()
  world.addChild(ipg)
  addWellGroups(world,wdata,curve)
  frame = TestFrame(world)
  view = frame.getOrbitView()
  frame.setVisible(True)

def addWellGroups(world,wdata,curve):
  logs = wdata.getLogsWith(curve)
  print "number of logs =",len(logs)
  for log in logs:
    pg = makePointGroup(log)
    world.addChild(pg)

def makePointGroup(log):
  n = log.n
  xyz = zerofloat(3*n)
  copy(n,0,1,log.x3,0,3,xyz)
  copy(n,0,1,log.x2,1,3,xyz)
  copy(n,0,1,log.x1,2,3,xyz)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.YELLOW)
  states.add(cs)
  pg = PointGroup(xyz)
  pg.setStates(states)
  return pg

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
