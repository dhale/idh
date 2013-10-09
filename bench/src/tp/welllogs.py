"""
Reads and reformats well logs from Teapot Dome.
Author: Dave Hale, Colorado School of Mines
Version: 2009.06.07
"""
from imports import *

#############################################################################
def main(args):
  #fixDensityLogs()
  #makeBinaryWellLogs("test")
  #makeBinaryWellLogs("deep")
  #makeBinaryWellLogs("shallow")
  #makeBinaryWellLogs("all")
  #viewWellCoordinates("deep")
  #viewReflectivity(490252305400)
  #viewWellCurves("deep","velocity")
  #viewWellCurves("deep","density")
  #viewWellCurves("deep","gamma")
  #viewWellCurves("deep","porosity")
  #viewWellsWithSeismic("all","velocity")
  #viewWellsWithSeismic("deep","velocity")
  #viewWellsWithSeismic("all","gamma")
  viewWellsWithSeismic("deep","gamma")
  #viewElevations("deep")
  pass

def viewReflectivity(id):
  setGlobals("all")
  wldata = WellLog.Data.readBinary(csmWellLogs)
  log = wldata.get(id)
  n = log.n
  v = log.v
  d = log.d
  z = log.z
  r = zerofloat(n)
  fnull = WellLog.NULL_VALUE
  for i in range(1,n):
    m = i-1
    vi,di = v[i],d[i]
    vm,dm = v[m],d[m]
    if vi!=fnull and di!=fnull and vm!=fnull and dm!=fnull:
      zi,zm = vi*di,vm*dm
      r[i] = (zm-zi)/(zm+zi)
  for i in range(n):
    if z[i]==fnull: z[i] = 0.0
    if v[i]==fnull: v[i] = 2.0
    if d[i]==fnull: d[i] = 2.0
  print "r: min =",min(r)," max =",max(r)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(z)
  sp.setSize(300,800)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(z,v)
  sp.setSize(300,800)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(z,d)
  sp.setSize(300,800)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(z,r)
  sp.setSize(300,800)

# Directories and files for well logs, headers, directional surveys
tpDir = "/data/seis/tpd/"
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

# Fix units for density logs; convert kg/m3 to g/cc.
# This should be unnecessary with latest WellLog.java.
def fixDensityLogs():
  for what in ["deep","shallow","all"]:
    print "fixing",what
    setGlobals(what)
    wldata = WellLog.Data.readBinary(csmWellLogs)
    for log in wldata.getLogsWith("density"):
      n,d = log.n,log.d
      dnull = WellLog.NULL_VALUE
      for i in range(n):
        if d[i]!=dnull:
          d[i] *= 0.001
    wldata.writeBinary(csmWellLogs)

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
  flist,zlist = [],[]
  for log in wldata.getLogsWith(curve):
    #log.despike(3)
    #log.smooth(25)
    f,z,y,x = log.getSamples(curve)
    flist.append(f)
    zlist.append(z)
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
  tv = PointsView(zlist,flist)
  sp.add(tv)

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
  frame = SimpleFrame(world)

def addWellGroups(world,wdata,curve):
  logs = wdata.getLogsWith(curve)
  print "number of logs =",len(logs)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.YELLOW)
  states.add(cs)
  wellsGroup = Group()
  wellsGroup.setStates(states)
  for log in logs:
    f,x1,x2,x3 = log.getSamples(curve,s1,s2,s3)
    pg = makePointGroup(x1,x2,x3)
    wellsGroup.addChild(pg)
  world.addChild(wellsGroup)

def makePointGroup(x1,x2,x3):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  pg = PointGroup(xyz)
  return pg

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
