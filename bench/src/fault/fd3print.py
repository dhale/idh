#############################################################################
# Slides for fault displacements in 3D images

from common import *

gmin,gmax,gint,glab = -5.5,5.5,2.0,"Amplitude"
smin,smax,sint,slab = -16.0,16.0,5.0,"Vertical component of throw (ms)"
t1min,t1max,t1int,t1lab = -0.00001,16.0,5.0,"Vertical component of throw (ms)"
flmin,flmax,flint,fllab = 0.5,1.0,0.1,"Fault likelihood"
background = Color.WHITE
pngDir = "png/"
#pngDir = None

def main(args):
  goFigures3()

def goFigures3():
  global s1,s2,s3
  s1,s2,s3 = samplingS1A()
  #s1,s2,s3 = samplingS1B()
  g = readImage("g")
  p2 = readImage("p2")
  p3 = readImage("p3")
  gs = readImage("gs")
  fl = readImage("fl")
  flt = readImage("flt")
  fs = makeFaultSurfer()
  surfs = getSurfs(fs)
  #xyz = surfs[5].sampleFaultDip()
  s = getShifts(gs,fs,surfs,p2,p3)
  t1,t2,t3 = getThrows(fs,surfs)
  s = mul(-4.0,s)
  t1 = mul(4.0,t1)
  print "s: min =",min(s)," max =",max(s)
  print "t1: min =",min(t1)," max =",max(t1)
  #kkkks = [[80,26,38,119]
  #         [80,55,47,146],
  #         [80,59,52,142]]
  #kkkks = [[80,55,47,146]]
  #kkkks = [[80,59,52,142]]
  kkkks = [[80,55,35,135]]
  global k1d,k1f,k2,k3,kp
  plot = plot3d
  #plot = plot3f
  for kkkk in kkkks:
    k1d,k1f,k2f,k3f = kkkk
    #for kp in range(0,35):
    for kp in range(0,1):
      k2 = k2f+kp
      k3 = k3f+kp
      """
      plot(g,a=s,amin=smin,amax=smax,amap=bwrNotch(1.0),
           alab=slab,aint=sint,png="s")
      plot(g,a=t1,amin=t1min,amax=t1max,amap=jetRamp(1.0),
           alab=t1lab,aint=t1int,png="t")
      """
      plot(g,a=t1,amin=t1min,amax=t1max,amap=jetRamp(0.0),
           alab=t1lab,aint=t1int,surfs=surfs,smax=t1max/4,png="ts")
      """
      plot(g,a=fl,amin=flmin,amax=flmax,amap=jetRamp(1.0),
           alab=fllab,aint=flint,png="fl")
      plot(g,a=flt,amin=flmin,amax=flmax,amap=jetRamp(1.0),
           alab=fllab,aint=flint,png="flt")
      #plot(g,a=fl,amin=flmin,amax=flmax,amap=jetRamp(0.0),
      #     alab=fllab,aint=flint,surfs=surfs,png="fls")
      plot(g,png="g")
      """

def makeFaultSurfer():
  fl = readImage("fl")
  fp = readImage("fp")
  ft = readImage("ft")
  fs = FaultSurfer3([fl,fp,ft])
  fs.setThreshold(0.5)
  return fs
def getSurfs(fs):
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,2000)
  return surfs
def getShifts(gs,fs,surfs,p2,p3):
  s = fs.findShifts(20.0,surfs,gs,p2,p3)
  print "s: min =",min(s)," max =",max(s)
  return s
def getThrows(fs,surfs):
  t1,t2,t3 = fs.findThrows(-0.0001,surfs)
  return t1,t2,t3 

def plot3d(g,a=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           xyz=None,surfs=None,smax=None,png=None):
  n1 = len(g[0][0])
  n2 = len(g[0])
  n3 = len(g)
  sf = SimpleFrame()
  sf.setBackground(background)
  if a==None:
    ipg = sf.addImagePanels(g)
    ipg.setClips(gmin,gmax)
  else:
    ipg = ImagePanelGroup2(g,a)
    if amap==None:
      amap = jetFill(0.8)
    ipg.setColorModel2(amap)
    if amin and amax:
      ipg.setClips1(gmin,gmax)
      ipg.setClips2(amin,amax)
    sf.world.addChild(ipg)
  ipg.setSlices(k1d,k2,k3)
  if xyz:
    pg = PointGroup(0.2,xyz)
    ss = StateSet()
    cs = ColorState()
    cs.setColor(Color.YELLOW)
    ss.add(cs)
    pg.setStates(ss)
    #ss = StateSet()
    #ps = PointState()
    #ps.setSize(5.0)
    #ss.add(ps)
    #pg.setStates(ss)
    sf.world.addChild(pg)
  if surfs:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    #if not smax:
    #  ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    for surf in surfs:
      #surf.blocky()
      if smax:
        xyz,uvw,rgb = surf.getXyzUvwRgbShifts(smax)
      else:
        xyz,uvw,rgb = surf.getXyzUvwRgb()
      #qg = QuadGroup(False,xyz,rgb)
      qg = QuadGroup(True,xyz,rgb) #qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
    sf.world.addChild(sg)
  #ipg.setSlices(209,12,18)
  #ipg.setSlices(200,0,0)
  #ipg.setSlices(80,9,13)
  #ipg.setSlices(80,9,209)
  #sf.setSize(1300,1100)
  #sf.setSize(1040,1124)
  if a==None:
    cbar = addColorBar3d(sf,glab,gint)
    ipg.addColorMapListener(cbar)
  else:
    cbar = addColorBar3d(sf,alab,aint)
    ipg.addColorMap2Listener(cbar)
  sf.setSize(1000,800)
  sf.setWorldSphere(n3/2,n2/2,n1/2,0.5*sqrt(n1*n1+n2*n2+n3*n3))
  #sf.orbitView.setAzimuthAndElevation(90,40)
  #sf.orbitView.setAzimuthAndElevation(-49,57)
  #sf.orbitView.setAzimuthAndElevation(60,60)
  #sf.orbitView.setAzimuthAndElevation(-90,80)
  #sf.orbitView.setScale(1.42)
  # good for subset a
  ipg.setSlices(80,38,119) # t26
  #ipg.setSlices(k1d,k2,k3)
  sf.orbitView.setAzimuthAndElevation(-73,51)
  sf.orbitView.setScale(1.35)
  sf.orbitView.setTranslate(Vector3(0.0750,0.0664,0.0441))
  # good for zoom of subset a above
  #sf.orbitView.setAzimuthAndElevation(-73,51)
  #sf.orbitView.setScale(8.806)
  #sf.orbitView.setTranslate(Vector3(-0.0261,0.0664,-0.2864))
  # good for zoom of overlap in subset a
  #ipg.setSlices(72,200,182)
  #sf.orbitView.setAzimuthAndElevation(68.8,32.7)
  #sf.orbitView.setScale(3.71)
  #sf.orbitView.setTranslate(Vector3(0.00957,0.45415,-0.17196))
  # good for subset b
  #ipg.setSlices(k1d,k2,k3)
  #sf.orbitView.setAzimuthAndElevation(-99,67)
  #sf.orbitView.setScale(1.56)
  #sf.orbitView.setTranslate(Vector3(0.0435,0.0550,-0.0157))
  # good for closeup view
  #sf.orbitView.setAzimuthAndElevation(225.72,44.38)
  #sf.orbitView.setScale(14.54)
  #sf.orbitView.setTranslate(Vector3(-0.4886,0.1457,-0.3072))
  sf.viewCanvas.setBackground(sf.getBackground())
  sf.setVisible(True)
  if png and pngDir:
    png = pngDir+sampling()+"t"+str(k1d)+png
    sf.paintToFile(png+".png");
    cbar.paintToPng(cbar.getWidth(),1.0,png+"c.png")

def plot3f(g,a=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1f,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1, 70)
  pp.setClips(gmin,gmax)
  if a:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar(alab)
    if aint:
      cb.setInterval(aint)
  else:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(2.0)
  pp.setInterval1(0.1)
  pp.setInterval2(1.0)
  pp.setInterval3(1.0)
  if a:
    pv12 = PixelsView(s1,s2,slice12(k3,a))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,slice13(k2,a))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,slice23(k1f,a))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(amap)
      if amin!=amax:
        pv.setClips(amin,amax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(170)
  pf.setFontSizeForSlide(1.0,0.8)
  pf.setSize(1000,800)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+sampling()+"f"+str(k1f)+str(kp)+png
    pf.paintToPng(360,7.0,png+".png")

def addColorBar3d(frame,clab,cint=None):
  cbar = ColorBar(clab)
  cbar.setFont(Font("Arial",Font.PLAIN,40))
  cbar.setBackground(background)
  if cint:
    cbar.setInterval(cint)
  cbar.setWidthMinimum(120)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

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
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

def setupForF3dSubset(subset):
  global s1,s2,s3
  global n1,n2,n3
  global dataDir,dataPre,kkks
  dataDir = "/data/seis/f3d/faults/"
  n1,n2,n3 = 462,951,591
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.004,0.000,0.000
  if subset=="a": # deeper normal faults
    dataPre = "a"
    j1,j2,j3 = 240+130,50,100
    m1,m2,m3 = 90,221,220
    kkks = [
      [26,38,119,"t26"]]
      #[26,38,119,"t26"], # ??? 26(39), 38, 119
      #[42,69,162,"t42"]] # ??? 42,     69, 162(124)
  elif subset=="b": # shallower conical faults
    dataPre = "b"
    #j1,j2,j3 = 240,50,100
    #m1,m2,m3 = 120,221,220
    j1,j2,j3 = 240+15,50,100
    m1,m2,m3 = 90,221,220
    kkks = [
      [55,47,146,"t55"],
      [59,52,142,"t59"]]
  elif subset=="s1": # conical fault above normal faults
    dataPre = "s1"
    j1,j2,j3 = 240,50,100
    m1,m2,m3 = 222,221,220
  n1,n2,n3 = m1,m2,m3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
