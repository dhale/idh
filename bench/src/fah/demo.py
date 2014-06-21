"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""

from fakeutils import *
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

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

def main(args):
  #goFakeData()
  #goSlopes()
  #goScan()
  goThin()
  #goSmooth()
  #goSurfing()
  #goDisplay("gx")

def goFakeData():
  folding = True
  faulting = True
  impedance = False
  wavelet = True
  noise = 0.0
  gx,p2,p3 = FakeData.seismicAndSlopes3d2014A(
      folding,faulting,impedance,wavelet,noise)
  writeImage(gxfile,gx)
  writeImage(p2kfile,p2)
  writeImage(p3kfile,p3)
  print "gx min =",min(gx)," max =",max(gx)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  #plot3(gx)
  #plot3(gx,p2,cmap=bwrNotch(1.0))
  #plot3(gx,p3,cmap=bwrNotch(1.0))

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 16.0,1.0,1.0,5.0
  p2,p3,ep = FaultScanner3.slopes(sigma1,sigma2,sigma3,pmax,gx)
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
  #plot3(gx,p2, cmin=-2,cmax=2,cmap=bwrNotch(1.0))
  #plot3(gx,p2k,cmin=-2,cmax=2,cmap=bwrNotch(1.0))
  #plot3(gx,p3, cmin=-2,cmax=2,cmap=bwrNotch(1.0))
  #plot3(gx,p3k,cmin=-2,cmax=2,cmap=bwrNotch(1.0))
  #plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0))

def goScan():
  print "goScan ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gxfile)
  gx = FaultScanner3.taper(10,0,0,gx);
  sigmaPhi,sigmaTheta = 4,20
  minPhi,maxPhi = -90,90
  minTheta,maxTheta = -20,20
  fsc = FaultScanner3(sigmaPhi,sigmaTheta)
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
  plot3(gx,fl,cmin=0,cmax=1,cmap=jetRamp(1.0))

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  plot3(gx,fl,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  fd = FaultScanner3.fd([fl,fp,ft])
  fd = abs(fd)
  plot3(gx,fd,cmap=jetRamp(1.0))
  #plot3(gx,fd,cmin=-0.1,cmax=0.1,cmap=bwrFill(0.3))
  return
  fl,fp,ft = FaultScanner3.thin([fl,fp,ft])
  writeImage(fltfile,fl)
  writeImage(fptfile,fp)
  writeImage(fttfile,ft)
  plot3(gx,fl,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  #plot3(gx,fp,cmin=-90,cmax=90,cmap=jetRamp(1.0))

def goSmooth():
  print "goSmooth ..."
  flstop = 0.2
  fsigma = 16.0
  gx = readImage(gxfile)
  fl = readImage(fltfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gsx = FaultScanner3.smooth(flstop,fsigma,p2,p3,fl,gx)
  writeImage(gsxfile,gsx)
  plot3(gx)
  plot3(gsx)
  p2 = readImage(p2kfile)
  p3 = readImage(p3kfile)
  gsx = FaultScanner3.smooth(flstop,fsigma,p2,p3,fl,gx)
  plot3(gsx)
  
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
  elif what=="gft1":
    g = readImage("g")
    ft1 = readImage("ft1")
    plot3(g,ft1,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  else:
    print "do not know how to display ",what

def goSurfing():
  g = readImage(gfile)
  gs = readImage(gsfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fs = FaultSurfer3([fl,fp,ft])
  fs.setThreshold(0.4)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,4000)
  plot3(g,surfs=surfs)
  s = fs.findShifts(20.0,surfs,gs,p2,p3)
  print "s: min =",min(s)," max =",max(s)
  writeImage(frsfile,s)
  t1,t2,t3 = fs.findThrows(-0.012345,surfs)
  writeImage(ft1file,t1)
  writeImage(ft2file,t2)
  writeImage(ft3file,t3)
  #plot3(g,surfs=surfs,smax=5.0)
  plot3(g,s,-2.0,2.0,cmap=bwrNotch(1.0))
  #plot3(g,surfs=surfs)
  #plot3(g,surfs=surfs,smax=-3.75)
  #plot3(g,surfs=surfs,smax= 3.75)
  #plot3(g,s,-10,10,cmap=bwrFill(0.7))
  #plot3(g,t1,-10.0,10.0,cmap=bwrFill(0.7))
  #plot3(g,t2,-0.50,0.50,cmap=bwrFill(0.7))
  #plot3(g,t3,-0.50,0.50,cmap=bwrFill(0.7))
  #plot3(g)

#############################################################################
# graphics

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,
          xyz=None,surfs=None,smax=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    ipg.setClips(-3.0,3.0)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-3.0,3.0)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if cmin and cmax:
      ipg.setClips2(cmin,cmax)
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
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
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
  ipg.setSlices(100,51,24);
  sf.setSize(700,700)
  vc = sf.getViewCanvas()
  ov = sf.getOrbitView()
  vc.setBackground(Color.WHITE)
  ov.setAzimuthAndElevation(40.0,25.0)
  ov.setScale(2.2)
  ov.setTranslate(Vector3(-0.0102,-0.0508,0.0395))

#############################################################################
run(main)
