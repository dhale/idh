"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
method = "b" # blended
logSet = "d" # deep logs only
#logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.6
logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.0,2.8
#logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
#logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
smin,smax = -5.5,5.5
#smooth = 50 # half-width of smoothing filter for logs
smooth = 0 # no smoothing

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

horizons = ["CrowMountainCRMT"]
"""
horizons = [
  "CrowMountainCRMT",
  "TensleepASand",
  "TensleepBbaseC1Dolo"]
"""

#pngDir = "png/"
pngDir = None

#k1,k2,k3 = 228,170,74 # 2D displays
#k1,k2,k3 = 228,170,106 # 2D displays
#k1,k2,k3 = 228,170,96 # 2D displays
k1,k2,k3 = 366,15,96 # 3D displays

def main(args):
  goInterp()
  #goFigures()
  #goImpedance()

def goInterp():
  global k1,k2,k3
  k1,k2,k3 = 366,15,96
  #gridBlendedP()
  #gridBlendedQ()
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  #display1(s,True,cmin=vmin,cmax=vmax)
  #display1(s,False)
  #display1(s,False,["CrowMountainCRMT","TensleepASand"])
  #display1(s,True,["CrowMountainCRMT","TensleepASand"])
  p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #t = readImage(tfile); print "t min =",min(t)," max =",max(t)
  display(s,p,vmin,vmax,logType)
  display(s,q,vmin,vmax,logType)
  #display(s,t,0.0,100.0,logType)
  #display(s,q,vmin,vmax,logType,["CrowMountainCRMT"])
  #display(s,q,vmin,vmax,logType,["TensleepASand"])

def goFigures():
  global k1,k2,k3
  k1,k2,k3 = 228,170,96 # intersect low-velocity layers
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  #p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  #q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #display3(s,None,0.0,0.0,"tpsz")
  #display3(s,p,vmin,vmax,"tppvb")
  #display3(s,q,vmin,vmax,"tpqvb")
  p = readImage("ig6/tppvb"); print "p min =",min(p)," max =",max(p)
  display3(s,p,vmin,vmax,"tppvb")
  p = readImage("tppvo09"); print "p min =",min(p)," max =",max(p)
  display3(s,p,vmin,vmax,"tppvb")

def gridBlendedP():
  e = getEigenTensors()
  bi = BlendedGridder3(e)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def gridBlendedQ():
  e = getEigenTensors()
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(pfile)
  t = readImage(tfile)
  t = clip(0.0,50.0,t)
  q = copy(p)
  bg.gridBlended(t,p,q)
  writeImage(qfile,q)

def getEigenTensors():
  e = readTensors(esfile)
  return e

def goImpedance():
  global k1,k2,k3,logLabel
  k1,k2,k3 = 366,15,96
  logType = None
  logLabel = "Impedance (g/cc x km/s)"
  d = readImage("tpqdb"); print "d min =",min(d)," max =",max(d)
  v = readImage("tpqvb"); print "v min =",min(v)," max =",max(v)
  s = readImage("tpsz"); print "s min =",min(s)," max =",max(s)
  i = mul(d,v)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply2XX(mul(0.5,log(i)),i)
  #imin,imax = min(i),max(i)
  imin,imax = -0.05,0.05
  display2S(s,i,imin,imax,logType)

def display2S(s,g,cmin,cmax,logType,horizons=[]):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  ipg = addImageToWorld(world,g)
  ipg.setClips(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,smooth=smooth)
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
    addLogsToWorld(world,logSet,logType,cmin,cmax,smooth=smooth)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)
  cbar = addColorBar(frame,logLabel)
  ipg.addColorMap2Listener(cbar)

def display1(s,wells=True,horizons=[],cmin=0,cmax=0):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  if wells:
    cbar = addColorBar(frame,logLabel)
    addLogsToWorld(world,logSet,logType,cmin,cmax,cbar,smooth=smooth)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)

def addColorBar(frame,label):
  cbar = ColorBar(logLabel)
  cbar.setFont(cbar.getFont().deriveFont(36.0))
  frame.add(cbar,BorderLayout.EAST)
  #frame.viewCanvas.setBackground(frame.getBackground())
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

def display3(s,g=None,cmin=0,cmax=0,png=None):
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
    cb = pp.addColorBar(logLabel)
    cb.setInterval(1.0)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(5.0)
  pp.setColorBarWidthMinimum(100)
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
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  pf.setSize(996,815)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

#############################################################################

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
