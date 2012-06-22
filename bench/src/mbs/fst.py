"""
Fault processing
"""

from mbsutils import *
setupForSubset("s1")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gfile = "gs" # input image (after smoothing with soblf)
gsfile = "gsf" # image after lsf with fault likelihoods
p2file = "p2" # inline slopes
p3file = "p3" # crossline slopes
snfile = "sn" # semblance numerators
sdfile = "sd" # semblance denominators
flfile = "fl" # fault likelihoods
fpfile = "fp" # fault strikes (phi)
ftfile = "ft" # fault dips (theta)
fltfile = "flt" # thinned fault likelihoods
fptfile = "fpt" # thinned fault strikes (phi)
fttfile = "ftt" # thinned fault dips (theta)

def main(args):
  goSlopes()
  #goSemblance()
  #goScan()
  #goThin()
  #goSmooth()
  #goSurfing()

def goSurfing():
  g = readImage(gfile)
  gs = readImage(gsfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSurfer3([fl,fp,ft])
  fs.setThreshold(0.6)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,4000)
  plot3(g,surfs=surfs)
  s = fs.findShifts(20.0,surfs,gs)
  print "s: min =",min(s)," max =",max(s)
  plot3(g,surfs=surfs,smax=5.0)
  plot3(g,s,-5.0,5.0,gmap=bwrNotch(1.0))
  #t1,t2,t3 = fs.findThrows(-0.12345,surfs)
  #plot3(g,surfs=surfs)
  #plot3(g,surfs=surfs,smax=-3.75)
  #plot3(g,surfs=surfs,smax= 3.75)
  #plot3(g,s,-10,10,gmap=bwrFill(0.7))
  #plot3(g,t1,-10.0,10.0,gmap=bwrFill(0.7))
  #plot3(g,t2,-0.50,0.50,gmap=bwrFill(0.7))
  #plot3(g,t3,-0.50,0.50,gmap=bwrFill(0.7))
  #plot3(g)

def goSmooth():
  g = readImage(gfile)
  #g = slog(g)
  fl = readImage(fltfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gs = FaultScanner3.smooth(16.0,p2,p3,fl,g)
  writeImage(gsfile,gs)
  plot3(g,fl,gmin=0.0,gmax=1.0,gmap=jetRamp(1.0))
  plot3(gs,fl,gmin=0.0,gmax=1.0,gmap=jetRamp(1.0))
  plot3(g)
  plot3(gs)

def goThin():
  thinned = True
  if not thinned:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    #rgf = RecursiveGaussianFilter(1)
    #rgf.applyXX0(fl,fl)
    #rgf.applyX0X(fl,fl)
    ref = RecursiveExponentialFilter(1.0)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply2(fl,fl); ref.apply2(fl,fl)
    ref.apply3(fl,fl); ref.apply3(fl,fl)
    fl,fp,ft = FaultScanner3.thin([fl,fp,ft])
    writeImage(fltfile,fl)
    writeImage(fptfile,fp)
    writeImage(fttfile,ft)
  else:
    fl = readImage(fltfile)
    fp = readImage(fptfile)
    ft = readImage(fttfile)
  g = readImage(gfile)
  #g = slog(g)
  plot3(g,fl,gmin=0.0,gmax=1.0,gmap=jetRamp(1.0))
  #plot3(g,fp,-90,90)
  #plot3(g,ft,-30,30)

def goScan():
  sn = readImage(snfile)
  sd = readImage(sdfile)
  sigmaPhi,sigmaTheta = 4,20
  fsc = FaultScanner3(sigmaPhi,sigmaTheta,[sn,sd])
  print "scanning ..."
  sw = Stopwatch()
  sw.restart()
  fl,fp,ft = fsc.scan(-90,90,-30,30)
  #fl,fp,ft = fsc.scan(40,50,-30,30)
  sw.stop()
  print "time =",sw.time(),", fl min =",min(fl)," max =",max(fl)
  writeImage(flfile,fl)
  writeImage(fpfile,fp)
  writeImage(ftfile,ft)
  g = readImage(gfile)
  #g = slog(g)
  plot3(g,fl,0,1)
  plot3(g,fp,-90,90)
  plot3(g,ft,-30,30)

def goSemblance():
  g = readImage(gfile)
  g = slog(g)
  print "reading slopes ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  print "computing semblance num/den"
  fse = FaultSemblance()
  g = fse.taper(10,g)
  sn,sd = fse.semblanceNumDen(p2,p3,g)
  print "writing semblance num/den"
  writeImage(snfile,sn)
  writeImage(sdfile,sd)

def goSlopes():
  g = readImage(gfile)
  #g = slog(g)
  #plot3(g)
  sigma1,sigma2,pmax = 16.0,4.0,5.0
  lsf = LocalSlopeFinder(sigma1,pmax)
  lsf.setSigma2(sigma2)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  lsf.findSlopes(g,p2,p3,None)
  p2 = clip(-1,1,p2)
  p3 = clip(-1,1,p3)
  print "p2: min =",min(p2)," max =",max(p2)
  print "p3: min =",min(p3)," max =",max(p3)
  plot3(g,p2,gmap=bwrNotch(1.0))
  plot3(g,p3,gmap=bwrNotch(1.0))
  #writeImage(p2file,p2)
  #writeImage(p3file,p3)

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

def plot3(f,g=None,gmin=None,gmax=None,gmap=None,h=None,
          xyz=None,surfs=None,smax=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame()
  sf.setBackground(Color.WHITE)
  if g==None:
    ipg = sf.addImagePanels(f)
    ipg.setClips(-1.0,1.0)
  else:
    ipg = ImagePanelGroup2(f,g)
    ipg.setClips1(-1.0,1.0)
    if gmap==None:
      gmap = jetFill(0.8)
    ipg.setColorModel2(gmap)
    if gmin and gmax:
      ipg.setClips2(gmin,gmax)
    sf.world.addChild(ipg)
  if h!=None:
    ipg = sf.addImagePanels(h)
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
