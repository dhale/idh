"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""

from fakeutils import *
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gsxfile = "gsx" # image after lsf with fault likelihoods
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
epfile  = "ep" # eigenvalue-derived planarity
snfile  = "sn" # semblance numerators
sdfile  = "sd" # semblance denominators
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
frsfile = "frs" # fault relative shifts
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)

# Processing begins here. When experimenting with one part of this demo, we
# can disable other parts that have already written results to files.
def main(args):
  #goFakeData()
  #goSlopes()
  #goScan()
  goSkin()
  #goThin()
  #goSmooth()
  #goSurfing()
  #goDisplay("gx")

def goFakeData():
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  #sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
  nplanar = 0
  conjugate = True
  conical = True
  impedance = False
  wavelet = True
  noise = 0.0
  gx,p2,p3 = FakeData.seismicAndSlopes3d2014A(
      sequence,nplanar,conjugate,conical,impedance,wavelet,noise)
  writeImage(gxfile,gx)
  writeImage(p2kfile,p2)
  writeImage(p3kfile,p3)
  print "gx min =",min(gx)," max =",max(gx)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  gmin,gmax,gmap = -3.0,3.0,ColorMap.GRAY
  if impedance:
    gmin,gmax,gmap = 0.0,1.4,ColorMap.JET
  plot3(gx,cmin=gmin,cmax=gmax,cmap=gmap)
  #plot3(gx,p2,cmap=bwrNotch(1.0))
  #plot3(gx,p3,cmap=bwrNotch(1.0))

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 16.0,1.0,1.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  p2k = readImage(p2kfile)
  p3k = readImage(p3kfile)
  print "p2  min =",min(p2)," max =",max(p2)
  print "p2k min =",min(p2k)," max =",max(p2k)
  print "p3  min =",min(p3)," max =",max(p3)
  print "p3k min =",min(p3k)," max =",max(p3k)
  print "ep min =",min(ep)," max =",max(ep)
  #plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0))
  #plot3(gx,p2k,cmin=-1,cmax=1,cmap=bwrNotch(1.0))
  #plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0))
  #plot3(gx,p3k,cmin=-1,cmax=1,cmap=bwrNotch(1.0))
  #plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0))

def goScan():
  print "goScan ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gxfile)
  gx = FaultScanner.taper(10,0,0,gx);
  sigmaPhi,sigmaTheta = 4,20
  minPhi,maxPhi = 0,360
  minTheta,maxTheta = 65,85
  fsc = FaultScanner(sigmaPhi,sigmaTheta)
  sw = Stopwatch()
  sw.restart()
  fl,fp,ft = fsc.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
  sw.stop()
  print "fl min =",min(fl)," max =",max(fl)
  print "fp min =",min(fp)," max =",max(fp)
  print "ft min =",min(ft)," max =",max(ft)
  print "time =",sw.time()
  writeImage(flfile,fl)
  writeImage(fpfile,fp)
  writeImage(ftfile,ft)
  plot3(gx)
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0))

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
  writeImage(fltfile,flt)
  writeImage(fptfile,fpt)
  writeImage(fttfile,ftt)
  plot3(gx)
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0))
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0))
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0))
  plot3(gx,ftt,cmin=60,cmax=90,cmap=jetFillExceptMin(1.0))

def goSmooth():
  print "goSmooth ..."
  flstop = 0.25
  fsigma = 16.0
  gx = readImage(gxfile)
  fl = readImage(fltfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,fl,gx)
  writeImage(gsxfile,gsx)
  #plot3(gx)
  plot3(gsx)

def goSkin():
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner([fl,fp,ft])
  fs.setGrowLikelihoods(0.2,0.7)
  fs.setMinSkinSize(4000)
  cells = fs.findCells()
  plot3(gx)
  plot3(gx,cells=cells)
  print "total number of cells =",len(cells)
  skins = fs.findSkins(cells)
  print "total number of skins =",len(skins)
  for iskin,skin in enumerate(skins):
    print "number of cells in skin",iskin,"=",skin.size()
    #cells = skin.getCells()
    #plot3(gx,cells=cells)
    plot3(gx,skins=[skin],links=True,curve=False,trace=False)
  plot3(gx,skins=skins,links=False,curve=True,trace=True);

  
def goDisplay(what):
  def show2(g1,g2):
    world = World()
    addImageToWorld(world,g1).setClips(-0.5,0.5)
    addImageToWorld(world,g2).setClips(-0.5,0.5)
    makeFrame(world).setSize(1200,900)
  if what=="gx":
    g = readImage("gx")
    plot3(g)
  elif what=="gs":
    gs = readImage("gs")
    #plot3(gs)
    g = readImage("g0")
    show2(g,gs)
  elif what=="gflt":
    g = readImage("g0")
    fl = readImage("flt")
    fl = pow(fl,0.5) # for display only?
    plot3(g,fl,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  elif what=="gfrs":
    g = readImage("g")
    s = readImage("frs")
    plot3(g,s,cmin=-5.0,cmax=5.0,cmap=bwrNotch(1.0))
  elif what=="gfs1":
    g = readImage("g")
    fs1 = readImage("fs1")
    plot3(g,fs1,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  else:
    print "do not know how to display ",what

#############################################################################
# graphics

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
def hueFill(alpha):
  return ColorMap.getHue(0.0,1.0,alpha)
def hueFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.getHue(0.0,1.0),a)

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,
          xyz=None,cells=None,skins=None,smax=None,
          links=False,curve=False,trace=False):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-3.0,3.0)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-3.0,3.0)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    sf.world.addChild(ipg)
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
  if cells:
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgb(0.5,cells)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    if links:
      size = 0.5 
    for skin in skins:
      if smax:
        xyz,uvw,rgb = skin.getXyzUvwRgbShifts(smax)
      else:
        xyz,uvw,rgb = skin.getCellXyzUvwRgb(size)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
      if curve or trace:
        cell = skin.getCellNearestCentroid()
        if curve:
          xyz = cell.getFaultCurveXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
        if trace:
          xyz = cell.getFaultTraceXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
      if links:
        xyz = skin.getCellLinksXyz()
        lg = LineGroup(xyz)
        sg.addChild(lg)
    sf.world.addChild(sg)
  #ipg.setSlices(95,21,51)
  ipg.setSlices(95,5,95)
  sf.setSize(700,700)
  vc = sf.getViewCanvas()
  ov = sf.getOrbitView()
  vc.setBackground(Color.WHITE)
  ov.setAzimuthAndElevation(-55.0,25.0)
  ov.setTranslate(Vector3(-0.0677,-0.0421,-0.0445))
  ov.setScale(2.2)

#############################################################################
run(main)
