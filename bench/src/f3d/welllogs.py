"""
Reads and reformats well logs for F3.
Author: Dave Hale, Colorado School of Mines
Version: 2012.12.29
"""
from imports import *

#############################################################################
def main(args):
  #makeBinaryWellLogs()
  viewWellCoordinates()
  #viewWellCurves("velocity")
  #viewWellCurves("density")
  #viewWellCurves("gamma")
  #viewWellCurves("porosity")
  #viewWellsWithSeismic("velocity")
  #viewWellsWithSeismic("gamma")
  pass

# Directories and files for well logs, headers, directional surveys
f3dDir = "/data/seis/f3d/"
odtWellLogsDir = f3dDir+"odt/"
csmWellLogsDir = f3dDir
odtWellLogs = odtWellLogsDir
csmWellLogs = csmWellLogsDir+"f3dwell.dat"

# File for 3D seismic depth image to view with wells.
csmSeismiczDir = f3dDir
csmSeismic = csmSeismiczDir+"f3d.dat"
s1 = Sampling(462,0.004,0.004)
s2 = Sampling(951,0.025,0.000)
s3 = Sampling(591,0.025,0.000)

# Reads well locations in CSM (x1,x2,x3) coordinates
def readLocations():
  wldata = WellLog.Data.readBinary(csmWellLogs)
  wldata.printInfo()
  nlog = wldata.size()
  x1 = zerofloat(nlog)
  x2 = zerofloat(nlog)
  x3 = zerofloat(nlog)
  ilog = 0
  for log in wldata.getAll():
    odt = Coordinates.Odt(log.xe,log.yn,log.ze)
    csm = Coordinates.Csm(odt)
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

def makeBinaryWellLogs():
  wldata = WellLog.Data(odtWellLogs)
  print "well log data"
  wldata.printInfo()
  #wldata.clip(s2,s3)
  #print "after clipping"
  #wldata.printInfo()
  wldata.writeBinary(csmWellLogs)

def viewWellCoordinates():
  wldata = WellLog.Data.readBinary(csmWellLogs)
  wldata.printInfo()
  spz1 = SimplePlot()
  sp23 = SimplePlot()
  zlist = []
  x1list = []
  x2list = []
  x3list = []
  for log in wldata.getAll():
    print log.name
    zlist.append(log.z)
    x1list.append(log.x1)
    x2list.append(log.x2[0])
    x3list.append(log.x3[0])
  tv = PointsView(zlist,x1list)
  spz1.add(tv)
  tv = PointsView(x2list,x3list)
  tv.setLineStyle(PointsView.Line.NONE)
  tv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  sp23.add(tv)

def viewWellCurves(curve):
  wldata = WellLog.Data.readBinary(csmWellLogs)
  flist,zlist = [],[]
  colors = [Color.RED,Color.GREEN,Color.BLUE,Color.MAGENTA]
  icolor = 0
  sp = SimplePlot()
  sp.setHLabel("Time (s)")
  if curve=="velocity":
    sp.setVLabel("Velocity (km/s)")
  elif curve=="density":
    sp.setVLabel("Density (gm/cc)")
  elif curve=="gamma":
    sp.setVLabel("Gamma Ray")
  elif curve=="porosity":
    sp.setVLabel("Porosity")
  for log in wldata.getLogsWith(curve):
    #log.despike(3)
    log.smooth(25)
    f,z,y,x = log.getSamples(curve)
    tv = PointsView(z,f)
    tv.setLineColor(colors[icolor])
    icolor = (icolor+1)%len(colors)
    sp.add(tv)

def viewWellsWithSeismic(curve):
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
