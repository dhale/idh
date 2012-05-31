#############################################################################
# Fault displacements from 3D images

import sys
from org.python.util import PythonObjectInputStream
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fault import *

#############################################################################

def main(args):
  #goS1()
<<<<<<< HEAD
  #goS1A()
=======
  goS1A()
>>>>>>> scripting
  goS1B()
def goS1():
  samplingS1()
  #goSurfingFake()
  #goBenchmarkSmoothers()
  #goSlopes()
  #goAlign()
  #goSemblance()
  #goScan()
  #goThin()
  #goSmooth()
  #goPartsAB()
def goS1A():
  samplingS1A()
  goSurfing()
  #goUnfault()
def goS1B():
  samplingS1B()
  goSurfing()
  #goUnfault()

def goSurfing():
  g = readImage("g"); g = slog(g)
  #h = readImage("h"); h = slog(h)
  gs = readImage("gs")
  fl = readImage("fl")
  fp = readImage("fp")
  ft = readImage("ft")
  fs = FaultSurfer3([fl,fp,ft])
  fs.setThreshold(0.5)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  #surfs = fs.getSurfsWithSize(surfs,2000)
  surfs = fs.getSurfsWithSize(surfs,2000)
  s = fs.findShifts(20.0,surfs,gs)
  print "s: min =",min(s)," max =",max(s)
  #t1,t2,t3 = fs.findThrows(-0.12345,surfs)
  #plot3(g,surfs=surfs)
  #plot3(g,surfs=surfs,smax=-3.75)
  #plot3(g,surfs=surfs,smax= 3.75)
  plot3(g,surfs=surfs,smax= 3.75)
  #plot3(h,surfs=surfs,smax= 3.75)
  #plot3(g,s,-10,10,gmap=bwrFill(0.7))
  #plot3(g,t1,-10.0,10.0,gmap=bwrFill(0.7))
  #plot3(g,t2,-0.50,0.50,gmap=bwrFill(0.7))
  #plot3(g,t3,-0.50,0.50,gmap=bwrFill(0.7))
  #plot3(g,s,-5,5,gmap=bwrNotch(1.0))
  #plot3(g,s,-5,5,gmap=bwrNotch(1.0))
  plot3(g,s,-1,1,gmap=bwrNotch(1.0))
  plot3(g)

def goUnfault():
  g = readImage("g")
  # Simon did the unfaulting
  h = readImage("h")
  plot3(g,h=h)

def goSurfingFake():
  #n1,n2,n3 = 101,102,103
  n1,n2,n3 = 51,52,53
  #n1,n2,n3 = 41,42,43
  #n1,n2,n3 = 11,12,13
  #g = sub(randfloat(n1,n2,n3),0.5)
  g = zerofloat(n1,n2,n3)
  f,p,t = Util.fakeSpheresFpt(n1,n2,n3)
  print "f max =",max(f)
  fs = FaultSurfer3([f,p,t])
  fs.setThreshold(0.5)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,100)
  #xyz = fs.sampleFaultDip(surfs[0])
  #s = fs.findShifts(surfs)
  #plot3(g,s,surfs=surfs)
  plot3(g,f,0,1,surfs=surfs)
  #sp = SimplePlot()
  #pv = sp.addPixels(fs.slice1(n1/2,s))
  #pv.setColorModel(ColorMap.JET)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)

def goSmooth():
  doSmooth = False
  if doSmooth:
    g = readImage("g")
    fl = readImage("flt")
    p2 = readImage("p2")
    p3 = readImage("p3")
    gs = FaultScanner3.smooth(16.0,p2,p3,fl,g)
    writeImage(gs,"gs")
  g = readImage("g"); g = slog(g)
  gs = readImage("gs"); gs = slog(gs)
  plot3(g)
  plot3(gs)

def goThin():
  doThin = False
  if doThin:
    f = readImage("fl")
    p = readImage("fp")
    t = readImage("ft")
    #rgf = RecursiveGaussianFilter(1)
    #rgf.applyXX0(f,f)
    #rgf.applyX0X(f,f)
    ref = RecursiveExponentialFilter(1.0)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply2(f,f); ref.apply2(f,f)
    ref.apply3(f,f); ref.apply3(f,f)
    f,p,t = FaultScanner3.thin([f,p,t])
    writeImage(f,"flt")
    writeImage(p,"fpt")
    writeImage(t,"ftt")
  g = readImage("g"); g = slog(g)
  f = readImage("flt")
  p = readImage("fpt")
  t = readImage("ftt")
  plot3(g,f,0,1)
  plot3(g,p,-90,90)
  plot3(g,t,-15,15)

def goScan():
  doScan = False
  if doScan:
    g = readImage("g"); g = slog(g)
    fse = FaultSemblance()
    #g = fse.taper(10,g)
    sn = readImage("sn")
    sd = readImage("sd")
    sigmaPhi,sigmaTheta = 4,20
    fsc = FaultScanner3(sigmaPhi,sigmaTheta,[sn,sd])
    print "scanning ..."
    sw = Stopwatch()
    sw.restart()
    fl,fp,ft = fsc.scan(-90,90,-15,15)
    sw.stop()
    print "time =",sw.time(),", fl min =",min(fl)," max =",max(fl)
    writeImage(fl,"fl")
    writeImage(fp,"fp")
    writeImage(ft,"ft")
  g = readImage("g"); g = slog(g)
  fl = readImage("fl")
  fp = readImage("fp")
  ft = readImage("ft")
  plot3(g,fl,0,1)
  plot3(g,fp,-90,90)
  plot3(g,ft,-15,15)

def goSemblance():
  g = readImage("g"); g = slog(g)
  print "reading slopes ..."
  p2 = readImage("p2")
  p3 = readImage("p3")
  print "computing semblance num/den"
  fse = FaultSemblance()
  g = fse.taper(10,g)
  sn0,sd0 = fse.semblanceNumDen(p2,p3,g)
  print "writing semblance num/den"
  writeImage(sn0,"sn")
  writeImage(sd0,"sd")
  return
  print "semblances for different vertical smoothings:"
  for sigma in [0,2,4,8]:
    ref = RecursiveExponentialFilter(sigma)
    sn = copy(sn0)
    sd = copy(sd0)
    ref.apply1(sn,sn)
    ref.apply1(sd,sd)
    s = fse.semblanceFromNumDen(sn,sd)
    print "sigma =",sigma," s min =",min(s)," max =",max(s)
    plot3(g,s,0,1)

def goAlign():
  g = readImage("g"); g = slog(g)
  fse = FaultSemblance()
  ref = RecursiveExponentialFilter(4)
  p2 = readImage("p2") # semblance with estimated slopes
  p3 = readImage("p3")
  sn,sd = fse.semblanceNumDen(p2,p3,g)
  ref.apply1(sn,sn)
  ref.apply1(sd,sd)
  s = fse.semblanceFromNumDen(sn,sd)
  plot3(g,s,0,1)
  p2 = zerofloat(n1,n2,n3) # semblance with zero slopes
  p3 = zerofloat(n1,n2,n3)
  sn,sd = fse.semblanceNumDen(p2,p3,g)
  ref.apply1(sn,sn)
  ref.apply1(sd,sd)
  s = fse.semblanceFromNumDen(sn,sd)
  plot3(g,s,0,1)

def goSlopes():
  g = readImage("g"); g = slog(g)
  plot3(g)
  fse = FaultSemblance()
  p2,p3 = fse.slopes(g)
  p2 = clip(-1,1,p2)
  p3 = clip(-1,1,p3)
  plot3(g,p2)
  plot3(g,p3)
  writeImage(p2,"p2")
  writeImage(p3,"p3")

def goBenchmarkSmoothers():
  g = readImage("g"); g = slog(g)
  fse = FaultSemblance()
  #g = fse.taper(10,g)
  sn = readImage("sn")
  sd = readImage("sd")
  smoothers = [
    (FaultScanner3.Smoother.ROTATE_AND_SHEAR,"ras"),
    (FaultScanner3.Smoother.FFT_GAUSSIAN,"fft")
  ]
  sw = Stopwatch()
  for smoother,name in smoothers:
    print "smoother =",name
    fsc = FaultScanner3(4,20,[sn,sd],smoother)
    print "  scanning ..."
    sw.restart()
    f,p,t = fsc.scan(90,90,-8,8)
    sw.stop()
    print "  time =",sw.time(),", f min =",min(f)," max =",max(f)
    plot3(g,f,0,1)

def goPartsAB():
  samplingS1A()
  n1a,n2a,n3a = n1,n2,n3
  j1a,j2a,j3a = 130,0,0
  samplingS1B()
  n1b,n2b,n3b = n1,n2,n3
  j1b,j2b,j3b = 15,0,0
  samplingS1()
  global dataSub
  for fileName in ["fl","fp","ft","flt","g","gs"]:
    dataSub = "s1/"
    x = readImage(fileName)
    x = copy(n1a,n2a,n3a,j1a,j2a,j3a,x)
    dataSub = "s1a/"
    writeImage(x,fileName)
    dataSub = "s1/"
    x = readImage(fileName)
    x = copy(n1b,n2b,n3b,j1b,j2b,j3b,x)
    dataSub = "s1b/"
    writeImage(x,fileName)

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))
 
#############################################################################
# data read/write

dataDir,dataSub = "",""
def readImage(fileName):
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+dataSub+fileName+".dat")
  ais.readFloats(f)
  ais.close()
  print "readImage: min =",min(f)," max =",max(f)
  return f
def writeImage(f,fileName):
  aos = ArrayOutputStream(dataDir+dataSub+fileName+".dat")
  aos.writeFloats(f)
  aos.close()
def samplingS1():
  global n1,n2,n3
  global s1,s2,s3
  global dataDir,dataSub
  dataDir = "/data/seis/f3d/faults/"
  dataSub = "s1/"
  n1,n2,n3 = 222,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.964,1.250,2.500
  s1,s2,s3 = samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3)
def samplingS1A():
  global n1,n2,n3
  global s1,s2,s3
  global dataDir,dataSub
  dataDir = "/data/seis/f3d/faults/"
  dataSub = "s1a/"
  n1,n2,n3 = 90,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 1.484,1.250,2.500
  s1,s2,s3 = samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3)
def samplingS1B():
  global n1,n2,n3
  global s1,s2,s3
  global dataDir,dataSub
  dataDir = "/data/seis/f3d/faults/"
  dataSub = "s1b/"
  n1,n2,n3 = 90,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 1.024,1.250,2.500
  s1,s2,s3 = samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3)
def samplings(n1,d1,f1,n2,d2,f2,n3,d3,f3):
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  print "s1:",s1.count,s1.delta,s1.first
  print "s2:",s2.count,s2.delta,s2.first
  print "s3:",s3.count,s3.delta,s3.first
  return s1,s2,s3

#############################################################################
# plotting

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,fillfloat(0.0,1.0/256,256))
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
  sf.setBackground(Color(255,255,255))
  if g==None:
    ipg = sf.addImagePanels(f)
  else:
    ipg = ImagePanelGroup2(f,g)
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
  #ipg.setSlices(209,12,18)
  #ipg.setSlices(200,0,0)
  #ipg.setSlices(80,9,13)
  #ipg.setSlices(80,9,209)
  #sf.setSize(1300,1100)
  sf.setSize(1040,1124)
  sf.setWorldSphere(n3/2,n2/2,n1/2,0.5*sqrt(n1*n1+n2*n2+n3*n3))
  #sf.orbitView.setAzimuthAndElevation(90,40)
  #sf.orbitView.setAzimuthAndElevation(-49,57)
  #sf.orbitView.setAzimuthAndElevation(60,60)
  #sf.orbitView.setAzimuthAndElevation(-90,80)
  #sf.orbitView.setScale(1.42)
  # good for subset a
  ipg.setSlices(80,38,119)
  sf.orbitView.setAzimuthAndElevation(-73,51)
  sf.orbitView.setScale(1.35)
  sf.orbitView.setTranslate(Vector3(0.0587,0.0626,0.0068))
  sf.viewCanvas.setBackground(sf.getBackground())
  # good for subset b
  #ipg.setSlices(80,47,146) # t55
  #ipg.setSlices(80,52,142) # t59
  #sf.orbitView.setAzimuthAndElevation(-99,67)
  #sf.orbitView.setScale(1.56)
  #sf.orbitView.setTranslate(Vector3(0.0435,0.0550,-0.0157))
  # good for closeup view
  #sf.orbitView.setAzimuthAndElevation(225.72,44.38)
  #sf.orbitView.setScale(14.54)
  #sf.orbitView.setTranslate(Vector3(-0.4886,0.1457,-0.3072))
  #sf.viewCanvas.setBackground(sf.getBackground())

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
