"""
Make figures (for slides and prints) for interpolation results.
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
logSet = "d" # deep logs only
smin,smax,sint,slog,slab = -5.5,5.5,2.0,None,"Amplitude"
vmin,vmax,vint,vlog,vlab = 2.4,5.6,0.5,"v","Velocity (km/s)"
dmin,dmax,dint,dlog,dlab = 2.2,2.8,0.1,"d","Density (g/cc)"
pmin,pmax,pint,plog,plab = 0.0,0.4,0.1,"p","Porosity"
gmin,gmax,gint,glog,glab = 0,200,50,"g","Gamma ray (API units)"
zmin,zmax,zint,zlog,zlab = 5,15,2,None,"Impedance (g/cc x km/s)"
vomit = [9,5,24,7,19,6] # velocity logs to leave out (omit)
vlabel = {9:"A",5:"B",24:"C",7:"D",19:"E",6:"F"} # labels for velocity logs
slides = True # True for slides, False for prints
prints = not slides # True for prints, False for slides
if slides:
  background = Color(253,254,255) # will be made transparent for slides
else:
  background = Color(255,255,255) # pure white for prints
#pngDir = None
pngDir = "png/"

"""
horizons = [
  "CrowMountainCRMT",
  "TensleepASand",
  "TensleepBbaseC1Dolo"]
"""


def main(args):
  #goFigures3d()
  #goFigures3f()
  goFiguresCv()

def goFigures3d():
  global k1,k2,k3
  k1,k2,k3 = 366,15,96 # good for 3D displays
  s = rimage("tpsz"); plot3d1(s,smin,smax,sint,slog,slab)
  #plot3d1(s,vmin,vmax,vint,vlog,vlab)
  #plot3d1(s,vmin,vmax,vint,vlog,vlab,["CrowMountainCRMT"])
  #plot3d1(s,vmin,vmax,vint,vlog,vlab,["TensleepASand"])
  #p = rimage("tppvb"); plot3d2(s,p,vmin,vmax,vint,vlog,vlab)
  #p = rimage("tppdb"); plot3d2(s,p,dmin,dmax,dint,dlog,dlab)
  #p = rimage("tpppb"); plot3d2(s,p,pmin,pmax,pint,plog,plab)
  #p = rimage("tppgb"); plot3d2(s,p,gmin,gmax,gint,glog,glab)
  #q = rimage("tpqvb"); plot3d2(s,q,vmin,vmax,vint,vlog,vlab)
  #plot3d2(s,q,vmin,vmax,vint,vlog,vlab,["CrowMountainCRMT"])
  #plot3d2(s,q,vmin,vmax,vint,vlog,vlab,["TensleepASand"])
  #q = rimage("tpqdb"); plot3d2(s,q,dmin,dmax,dint,dlog,dlab)
  #q = rimage("tpqpb"); plot3d2(s,q,pmin,pmax,pint,plog,plab)
  #q = rimage("tpqgb"); plot3d2(s,q,gmin,gmax,gint,glog,glab)
  #q = getImpedance();  plot3d2(s,q,zmin,zmax,zint,zlog,zlab)

def goFigures3f():
  kkks = [
    [100,246,54,"z10/"], # 1.0 km - bad logs; deep fault
    [225,170,96,"z15/"]] # 1.5 km - low-velocity layers
  #k1,k2,k3,subd = 100,246,54,"z10/" # 1.0 km - bad logs; deep fault
  #k1,k2,k3,subd = 225,170,96,"z15/" # 1.5 km - low-velocity layers
  ##k1,k2,k3 = 228,170,96 # intersect low-velocity layers
  ##k1,k2,k3 = 225,281,54 # depth = 1.5 km
  ##k1,k2,k3 = 110,170,96 # depth = 1.04 km - shallow low-gamma layer
  ##k1,k2,k3 = 75,170,96 # depth = 0.9 km
  global k1,k2,k3,subd
  for kkk in kkks:
    k1,k2,k3,subd = kkk
    s=rimage("tpsz");plot3f(s,None,smin,smax,sint,slog,slab,"tpsz")
    p=rimage("tppvb");plot3f(s,p,vmin,vmax,vint,vlog,vlab,"tppvb")
    q=rimage("tpqvb");plot3f(s,q,vmin,vmax,vint,vlog,vlab,"tpqvb")
    p=rimage("tppdb");plot3f(s,p,dmin,dmax,dint,dlog,dlab,"tppdb")
    q=rimage("tpqdb");plot3f(s,q,dmin,dmax,dint,dlog,dlab,"tpqdb")
    p=rimage("tpppb");plot3f(s,p,pmin,pmax,pint,plog,plab,"tpppb")
    q=rimage("tpqpb");plot3f(s,q,pmin,pmax,pint,plog,plab,"tpqpb")
    p=rimage("tppgb");plot3f(s,p,gmin,gmax,gint,glog,glab,"tppgb")
    q=rimage("tpqgb");plot3f(s,q,gmin,gmax,gint,glog,glab,"tpqgb")

def goFiguresCv():
  #plotWellPoints(1.0)
  #plotWellPoints(1.5)
  #for omit in vomit:
  #  plotLogs(omit)
  plotLogs(6)

def plot3d1(s,cmin,cmax,cint,logType,logLabel,horizons=[]):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.viewCanvas.setBackground(background)
  frame.orbitView.setAzimuth(-65.0)
  if logLabel:
    cbar = addColorBar3d(frame,logLabel,cint)
    ipg.addColorMapListener(cbar)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,cbar)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)

def plot3d2(s,g,cmin,cmax,cint,logType,logLabel,horizons=[],cval=0):
  world = World()
  ipg = addImage2ToWorld(world,s,g)
  ipg.setClips1(smin,smax)
  ipg.setClips2(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  if cval:
    addContourToWorld(world,g,cval)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.orbitView.setAzimuth(-65.0)
  frame.viewCanvas.setBackground(background)
  cbar = addColorBar3d(frame,logLabel,cint)
  ipg.addColorMap2Listener(cbar)

def addContourToWorld(world,g,cval):
  mc = MarchingCubes(s1,s2,s3,g)
  mc.setSwap13(True)
  c = mc.getContour(cval)
  tg = TriangleGroup(c.i,c.x,c.u)
  tg.setColor(Color.RED)
  world.addChild(tg)

def addColorBar3d(frame,label,cint=None):
  cbar = ColorBar(label)
  if slides:
    cbar.setFont(Font("Arial",Font.PLAIN,60)) # determined experimentally
  else:
    cbar.setFont(Font("Arial",Font.PLAIN,48)) # ~ 8*1460/240 for one-column
  cbar.setBackground(background)
  if cint:
    cbar.setInterval(cint)
  if slides:
    cbar.setPreferredSize(Dimension(250,100))
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def plot3f(s,g=None,cmin=0,cmax=0,cint=None,clog=None,clab=None,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Depth (km)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(smin,smax)
  if g:
    pp.setLineColor(Color.BLACK)
    cb = pp.addColorBar(clab)
    if cint:
      cb.setInterval(cint)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(2.0)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if g:
    pv12 = PixelsView(s1,s2,slice12(k3,g))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,g))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,g))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
    if clog:
      x2,x3 = getWellIntersections(logSet,clog,s1.getValue(k1))
      pv23 = PointsView(x2,x3)
      pv23.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
      pv23.setLineStyle(PointsView.Line.NONE)
      pv23.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv23.setMarkSize(6.0)
      pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  if slides:
    pp.setColorBarWidthMinimum(140)
    pf.setFontSizeForSlide(1.0,1.0)
    pf.setSize(1033,811)
  else:
    pp.setColorBarWidthMinimum(125)
    pf.setFontSizeForPrint(8,240)
    pf.setSize(995,790)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(400,3.3,pngDir+subd+png+".png")

def plotLogs(k):
  fw,x1w,x2w,x3w = readLog(wellLog(k))
  fg,x1g,x2g,x3g = readLog(griddedLog(k))
  fq,x1q,x2q,x3q = readLog(blendedLog(k))
  fs,x1s,x2s,x3s = smoothLog(fg),x1g,x2g,x3g
  #fs,x1s,x2s,x3s = fg,x1g,x2g,x3g
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if slides:
    sp.setSize(450,900)
    sp.setFontSizeForSlide(0.5,1.0)
  else:
    sp.setSize(450,900)
    sp.setFontSizeForPrint(8,160)
  #sp.setTitle("Log "+str(k))
  sp.setBackground(background)
  sp.plotPanel.setHInterval(1.0)
  sp.plotPanel.setVInterval(0.2)
  sp.setHLimits(2.800,6.200)
  sp.setVLimits(0.575,2.025)
  sp.setHLabel("Velocity (km/s)")
  sp.setVLabel("Depth (km)")
  if True:
    pv = sp.addPoints(x1w,fw)
    pv.setLineColor(Color.GRAY)
  if True:
    pv = sp.addPoints(x1s,fs)
    pv.setLineColor(Color.BLACK)
    pv.setLineWidth(4.0)
  if True:
    pv = sp.addPoints(x1q,fq)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(4.0)
  if pngDir:
    png = pngDir+"cvl/tpwsqv"+vlabel[k]+".png"
    #png = pngDir+"cvl/tpqv"+vlabel[k]+".png"
    sp.paintToPng(600,3.3,png)

def plotWellPoints(x1):
  x2,x3 = getWellIntersections(logSet,"v",x1)
  pv = PointsView(x2,x3)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(8.0)
  if slides:
    ap = PlotPanel.AxesPlacement.NONE
    pp = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP,ap)
  else:
    pp = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  pp.setHLimits(s2.first,s2.last)
  pp.setVLimits(s3.first,s3.last)
  if prints:
    pp.setHLabel("Crossline (km)")
    pp.setVLabel("Inline (km)")
  pp.addTiledView(pv)
  pf = PlotFrame(pp)
  if slides:
    pf.setSize(600,293)
  else:
    pf.setSize(600,333)
    pf.setFontSizeForPrint(8,240)
  pf.setBackground(background)
  pf.setVisible(True)
  if pngDir:
    png = pngDir+"cvl/tpwp"+str(int(0.5+10*x1))+".png"
    pf.paintToPng(600,3.3,png)

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

def getImpedance():
  d = rimage("tpqdb"); print "d min =",min(d)," max =",max(d)
  v = rimage("tpqvb"); print "v min =",min(v)," max =",max(v)
  i = mul(d,v)
  return i

def rimage(fileName):
  return readImage(fileName)

def smoothLog(f):
  #sigma = 2.0*s1.delta/0.0001524 # 0.0001524 km = 6 inches
  sigma = 2.5
  lpad = int(3.0*sigma)
  n = len(f)
  npad = lpad+n+lpad
  fpad = zerofloat(npad)
  gpad = zerofloat(npad)
  for i in range(lpad):
    fpad[i] = f[0]
  copy(n,0,f,lpad,fpad)
  for i in range(lpad+n,npad):
    fpad[i] = f[-1]
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply0(fpad,gpad)
  g = copy(n,lpad,gpad)
  return g 

def readLog(fileName):
  ais = ArrayInputStream(getSeismicDir()+fileName+".dat")
  n = ais.readInt()
  f = zerofloat(n)
  x1 = zerofloat(n)
  x2 = zerofloat(n)
  x3 = zerofloat(n)
  ais.readFloats(f)
  ais.readFloats(x1)
  ais.readFloats(x2)
  ais.readFloats(x3)
  ais.close()
  return f,x1,x2,x3

def readWellLog(set,type,index):
  fl,x1l,x2l,x3l = readLogSamples(set,type,smooth)
  return fl[index],x1l[index],x2l[index],x3l[index]

def griddedOmit(omit):
  return "tpg"+omitSuffix(omit)
def nearestOmit(omit):
  return "tpp"+omitSuffix(omit)
def mintimeOmit(omit):
  return "tpt"+omitSuffix(omit)
def blendedOmit(omit):
  return "tpq"+omitSuffix(omit)
def wellLog(ilog):
  return "tpw"+ilogSuffix(ilog)
def griddedLog(ilog):
  return "tpg"+ilogSuffix(ilog)
def nearestLog(ilog):
  return "tpp"+ilogSuffix(ilog)
def mintimeLog(ilog):
  return "tpt"+ilogSuffix(ilog)
def blendedLog(ilog):
  return "tpq"+ilogSuffix(ilog)
def omitSuffix(omit): # suffix for image files with a log omitted
  if omit<10:
    return "vo0"+str(omit)
  else:
    return "vo"+str(omit)
def ilogSuffix(ilog): # suffix for log files with specified index
  if ilog<10:
    return "vl0"+str(ilog)
  else:
    return "vl"+str(ilog)

#############################################################################
run(main)
