"""
Fault processing
"""

from schutils import *
setupForSubset("s2")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gfile = "g" # input image (maybe after bilateral filtering)
gsfile = "gs" # image after lsf with fault likelihoods
p2file = "p2" # inline slopes
p3file = "p3" # crossline slopes
epfile = "ep" # eigenvalue-based planarity
snfile = "sn" # semblance numerators
sdfile = "sd" # semblance denominators
flfile = "fl" # fault likelihoods
fpfile = "fp" # fault strikes (phi)
ftfile = "ft" # fault dips (theta)
fltfile = "flt" # thinned fault likelihoods
fptfile = "fpt" # thinned fault strikes (phi)
fttfile = "ftt" # thinned fault dips (theta)
frsfile = "frs" # fault relative shifts
ft1file = "ft1" # 1st component of fault throws
ft2file = "ft2" # 2nd component of fault throws
ft3file = "ft3" # 3rd component of fault throws

def main(args):
  #goDisplay("g")
  #goSlopes()
  #goSemblance()
  #goScan()
  #goThin()
  #goSmooth()
  #goSurfing()
  #goDisplay("g")
  goDisplay("gs")
  #goDisplay("gflt")
  #goDisplay("gfrs")
  #goDisplay("gft1")

def goDisplay(what):
  def show2(g1,g2):
    world = World()
    addImageToWorld(world,g1).setClips(-0.5,0.5)
    addImageToWorld(world,g2).setClips(-0.5,0.5)
    makeFrame(world).setSize(1200,900)
  if what=="g":
    g = readImage("g")
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

def goSmooth():
  smoothed = False
  g = readImage(gfile)
  if smoothed:
    gs = readImage(gsfile)
  else:
    g = readImage(gfile)
    fl = readImage(fltfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    gs = FaultScanner3.smooth(16.0,p2,p3,fl,g)
    writeImage(gsfile,gs)
    #plot3(g,fl,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
    #plot3(gs,fl,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  world = World()
  addImageToWorld(world,g).setClips(-1,1)
  addImageToWorld(world,gs).setClips(-1,1)
  makeFrame(world).setSize(1200,900)

def goThin():
  thinned = False
  if not thinned:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    rgf = RecursiveGaussianFilter(1)
    rgf.applyXX0(fl,fl)
    rgf.applyX0X(fl,fl)
    #ref = RecursiveExponentialFilter(1.0)
    #ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    #ref.apply2(fl,fl); ref.apply2(fl,fl)
    #ref.apply3(fl,fl); ref.apply3(fl,fl)
    fl,fp,ft = FaultScanner3.thin([fl,fp,ft])
    writeImage(fltfile,fl)
    writeImage(fptfile,fp)
    writeImage(fttfile,ft)
  else:
    fl = readImage(fltfile)
    fp = readImage(fptfile)
    ft = readImage(fttfile)
  g = readImage(gfile)
  plot3(g,fl,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0))
  plot3(g,fp,cmin=-90,cmax=90,cmap=jetRamp(1.0))
  #plot3(g,ft,cmin=-30,cmax=30)

def goScan():
  print "goScan ..."
  sn = readImage(snfile)
  sd = readImage(sdfile)
  sigmaPhi,sigmaTheta = 8,40 # sampling intervals are relatively small
  #sigmaPhi,sigmaTheta = 4,20
  minPhi,maxPhi = -90,90
  minTheta,maxTheta = -25,25
  #minPhi,maxPhi = 40,40
  #minTheta,maxTheta = -10,-10
  fsc = FaultScanner3(sigmaPhi,sigmaTheta,[sn,sd])
  print "scanning ..."
  sw = Stopwatch()
  sw.restart()
  fl,fp,ft = fsc.scan(minPhi,maxPhi,minTheta,maxTheta)
  sw.stop()
  print "time =",sw.time(),", fl min =",min(fl)," max =",max(fl)
  writeImage(flfile,fl)
  writeImage(fpfile,fp)
  writeImage(ftfile,ft)
  g = readImage(gfile)
  plot3(g)
  plot3(g,fl,cmin=0,cmax=1,cmap=jetRamp(1.0))
  if minPhi<maxPhi:
    plot3(g,fp,minPhi,maxPhi)
  if minTheta<maxTheta:
    plot3(g,ft,minTheta,maxTheta)
  print "done"

def goSemblance():
  print "goSemblance"
  g = readImage(gfile)
  #print "applying log gain ..."
  #g = slog(g)
  print "reading slopes ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  print "computing semblance num/den ..."
  fse = FaultSemblance()
  g = fse.taper(10,g)
  sn,sd = fse.semblanceNumDen(p2,p3,g)
  print "writing semblance num/den ..."
  writeImage(snfile,sn)
  writeImage(sdfile,sd)
  print "done"

def goSlopes():
  g = readImage(gfile)
  #g = slog(g)
  #plot3(g)
  sigma1,sigma2,pmax = 16.0,2.0,5.0
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf.findSlopes(g,p2,p3,ep)
  print "p2: min =",min(p2)," max =",max(p2)
  print "p3: min =",min(p3)," max =",max(p3)
  print "ep: min =",min(ep)," max =",max(ep)
  plot3(g,p2,cmin=-1,cmax=1,cmap=bwrNotch(1.0))
  plot3(g,p3,cmin=-1,cmax=1,cmap=bwrNotch(1.0))
  plot3(g,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0))
  writeImage(p2file,p2)
  writeImage(p3file,p3)

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

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
    """
    if i<96:
      a[i] = 1.0
    elif i<128:
      a[i] = alpha*(128.0-i)/32.0
    elif i<160:
      a[i] = alpha*(i-127.0)/32.0
    else:
      a[i] = 1.0
    """
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,h=None,
          xyz=None,surfs=None,smax=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  sf.setBackground(Color.WHITE)
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    ipg.setClips(-0.5,0.5)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-0.5,0.5)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if cmin and cmax:
      ipg.setClips2(cmin,cmax)
    sf.world.addChild(ipg)
  if h!=None:
    ipg = sf.addImagePanels(s1,s2,s3,h)
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
  sf.orbitView.setScale(2.0)
  sf.viewCanvas.setBackground(sf.getBackground())

#############################################################################
run(main)
