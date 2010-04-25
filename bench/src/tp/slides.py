"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
logSet = "d" # deep logs only
logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.6
#logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.2,2.8
#logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
#logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
method = "b" # blended
smin,smax = -5.5,5.5

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors (structure tensors)
esfile = "tpets" # eigen-tensors scaled by semblances
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,
gfile = "tpg"+logType # simple gridding with null for unknown samples
pfile = "tpp"+logType+method # values of nearest known samples
qfile = "tpq"+logType+method # output of blended gridder
tfile = "tpt"+logType+method # times to nearest known samples
def setSubDir(dir):
  global sfile,efile,esfile,s1file,s2file,s3file,gfile,pfile,qfile,tfile
  sfile = dir+sfile
  efile = dir+efile
  esfile = dir+esfile
  s1file = dir+s1file
  s2file = dir+s2file
  s3file = dir+s3file
  gfile = dir+gfile
  pfile = dir+pfile
  qfile = dir+qfile
  tfile = dir+tfile
#setSubDir("igfig/")
setSubDir("ig6/")

horizons = ["CrowMountainCRMT"]
"""
horizons = [
  "CrowMountainCRMT",
  "TensleepASand",
  "TensleepBbaseC1Dolo"]
"""

pngDir = "png/"
#pngDir = None

#k1,k2,k3 = 228,170,74 # 2D displays
#k1,k2,k3 = 228,170,106 # 2D displays
#k1,k2,k3 = 228,170,96 # 2D displays
k1,k2,k3 = 366,15,96 # 3D displays

def main(args):
  goFigure3()
  #goFigureS()
  #goImpedance()

def goFigure3():
  global k1,k2,k3
  k1,k2,k3 = 366,15,96
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  #display1(s,True,cmin=vmin,cmax=vmax)
  #display1(s,False)
  #display1(s,False,["CrowMountainCRMT"])
  #display1(s,False,["CrowMountainCRMT","TensleepASand"])
  display1(s,True,["CrowMountainCRMT","TensleepASand"],vmin,vmax)
  #p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #t = readImage(tfile); print "t min =",min(t)," max =",max(t)
  #display(s,p,vmin,vmax,logType)
  #display(s,q,vmin,vmax,logType)
  #display(s,t,0.0,200.0,logType)
  display(s,q,vmin,vmax,logType,["CrowMountainCRMT"])
  display(s,q,vmin,vmax,logType,["TensleepASand"])

def goFigureS():
  global k1,k2,k3
  #k1,k2,k3 = 228,170,96 # intersect low-velocity layers
  k1,k2,k3 = 225,170,96 # intersect low-velocity layers at 1.5 km
  #k1,k2,k3 = 225,281,54 # depth = 1.5 km
  #k1,k2,k3 = 110,170,96 # depth = 1.04 km - shallow low-gamma layer
  #k1,k2,k3 = 100,246,54 # depth = 1.0 km - show bad logs; good deep fault
  #k1,k2,k3 = 75,170,96 # depth = 0.9 km
  #dir = "igfig/"
  dir = "ig6/"
  s = readImage("tpsz")
  #display3(s,None,None,-5.0,5.0,2.0,"Amplitude","tpsz")
  p = readImage(dir+"tppvbii")
  display3(s,p,"v",2.4,5.6,0.5,"Velocity (km/s)","tppvbii")
  q = readImage(dir+"tpqvbii") 
  display3(s,q,"v",2.4,5.6,0.5,"Velocity (km/s)","tpqvbii")
  """
  p = readImage(dir+"tppdb") 
  display3(s,p,"d",2.2,2.8,0.2,"Density (g/cc)","tppdb")
  q = readImage(dir+"tpqdb") 
  display3(s,q,"d",2.2,2.8,0.2,"Density (g/cc)","tpqdb")
  p = readImage(dir+"tpppb") 
  display3(s,p,"p",0.0,0.4,0.1,"Porosity","tpppb")
  q = readImage(dir+"tpqpb") 
  display3(s,q,"p",0.0,0.4,0.1,"Porosity","tpqpb")
  p = readImage(dir+"tppgb") 
  display3(s,p,"g",0,200,50,"Gamma ray (API units)","tppgb")
  q = readImage(dir+"tpqgb") 
  display3(s,q,"g",0,200,50,"Gamma ray (API units)","tpqgb")
  """

def goImpedance():
  global k1,k2,k3,logLabel
  k1,k2,k3 = 366,15,96
  logType = None
  logLabel = "Impedance (g/cc x km/s)"
  d = readImage("ig6/tpqdb"); print "d min =",min(d)," max =",max(d)
  v = readImage("ig6/tpqvb"); print "v min =",min(v)," max =",max(v)
  s = readImage("ig6/tpsz"); print "s min =",min(s)," max =",max(s)
  i = mul(d,v)
  imin,imax = 5,15
  display(s,i,imin,imax,logType)
  #display2S(s,i,imin,imax,logType)

def gridBlendedP():
  e = getEigenTensors(0.01)
  bi = BlendedGridder3(e)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def gridBlendedQ():
  e = getEigenTensors(0.0)
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(pfile)
  t = readImage(tfile)
  t = clip(0.0,10.0,t)
  q = copy(p)
  bg.gridBlended(t,p,q)
  writeImage(qfile,q)

def getEigenTensors(eps):
  e = readTensors(efile)
  s1 = readImage(s1file); print "s1 min =",min(s1)," max =",max(s1)
  s2 = readImage(s2file); print "s2 min =",min(s2)," max =",max(s2)
  s3 = readImage(s3file); print "s3 min =",min(s3)," max =",max(s3)
  s1 = clip(eps,1.0,s1)
  s2 = clip(eps,1.0,s2)
  s3 = clip(eps,1.0,s3)
  e.setEigenvalues(s3,s2,s1)
  #writeTensors(esfile,e)
  return e

background = Color(254,254,254)
def display2S(s,g,cmin,cmax,logType,horizons=[]):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  ipg = addImageToWorld(world,g)
  ipg.setClips(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)

def display(s,g,cmin,cmax,logType,horizons=[]):
  world = World()
  ipg = addImage2ToWorld(world,s,g)
  ipg.setClips1(smin,smax)
  ipg.setClips2(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.orbitView.setAzimuth(-65.0)
  frame.viewCanvas.setBackground(background)
  cbar = addColorBar(frame,logLabel)
  ipg.addColorMap2Listener(cbar)

def display1(s,wells=True,horizons=[],cmin=0,cmax=0):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.viewCanvas.setBackground(background)
  frame.orbitView.setAzimuth(-65.0)
  if wells:
    cbar = addColorBar(frame,logLabel)
    addLogsToWorld(world,logSet,logType,cmin,cmax,cbar)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)

def addColorBar(frame,label):
  cbar = ColorBar(logLabel)
  cbar.setFont(Font("Arial",Font.PLAIN,60)) # size by experimenting
  cbar.setBackground(background)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def display2(s,g=None,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.getPlotPanel().setColorBarWidthMinimum(80)
  pv = sp.addPixels(s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g!=None:
    pv = sp.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.getJet(0.3))
    if cmin!=cmax:
      pv.setClips(cmin,cmax)

def display3(s,g=None,type=None,cmin=0,cmax=0,cint=None,clab=None,png=None):
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
  #pp.setColorBarWidthMinimum(125)
  pp.setColorBarWidthMinimum(140)
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
    if type:
      x2,x3 = getWellIntersections(logSet,type,s1.getValue(k1))
      pv23 = PointsView(x2,x3)
      pv23.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
      pv23.setLineStyle(PointsView.Line.NONE)
      pv23.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv23.setMarkSize(6.0)
      pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  #pf.setFontSizeForPrint(8,225)
  #pf.setSize(995,790)
  pf.setFontSizeForSlide(1.0,1.0)
  pf.setSize(1033,811)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(400,3.3,pngDir+png+".png")

#############################################################################

def dumpTensors():
  e = getEigenTensors()
  d = zerofloat(6)
  for i3 in range(14,17):
    for i2 in range(110,113):
      for i1 in range(249,252):
        e.getTensor(i1,i2,i3,d)
        dump(d)

def getScatteredSamples():
  g = readImage(gfile)
  f,x1,x2,x3 = SimpleGridder3.getGriddedSamples(0.0,s1,s2,s3,g)
  return f,x1,x2,x3

def gridNearest():
  f,x1,x2,x3 = getScatteredSamples()
  print "got scattered samples: n =",len(f)
  ni = NearestGridder3(f,x1,x2,x3)
  print "constructed nearest gridder"
  g = ni.grid(s1,s2,s3)
  print "gridding complete: min =",min(g)," max =",max(g)
  writeImage(gfile,g)

def gridSibson():
  f,x1,x2,x3 = getScatteredSamples()
  print "got scattered samples: n =",len(f)
  si = SibsonGridder3(f,x1,x2,x3)
  print "constructed Sibson gridder"
  g = si.grid(s1,s2,s3)
  print "gridding complete: min =",min(g)," max =",max(g)
  writeImage(gfile,g)

def testSpd():
  e = getEigenTensors()
  lsf = LocalSmoothingFilter(0.01,10000)
  t = readImage(tfile)
  t = clip(0.0,10.0,t)
  s = mul(t,t)
  lsf.testSpd(e,0.5,s)

def gridBlended2():
  sdir = "s3_84/"
  s = readSlice3(sdir+"tpsz")
  p = readSlice3(sdir+"tpgd"); cmin,cmax = 2000,2800
  #display2(s)
  display2(s,p,cmin,cmax)
  q = copy(p)
  lof = LocalOrientFilter(8)
  e = lof.applyForTensors(s)
  e.setEigenvalues(0.01,1.0)
  bi = BlendedGridder2(e)
  t = bi.gridNearest(0.0,p)
  t = clip(0,100,t)
  e.setEigenvalues(0.0,1.0)
  bi.setTensors(e)
  bi.setSmoothness(0.5)
  bi.gridBlended(t,p,q)
  #display2(s,p,cmin,cmax)
  display2(s,q,cmin,cmax)

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

#############################################################################
run(main)
