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
  #goBenchmarkSmoothers()
  #goSlopes()
  #goAlign()
  #goSemblance()
  #goScan()
  #goThin()
  #goSmooth()
  #goShifts()
  #goSurfingFake()
  #goSurfing()
  #goUnfault()
  #goSubsets()

def go120to90():
  global dataDir,dataPre
  dataDir = "/data/seis/f3d/faults/"
  dataPre = ""
  def subset(fileName):
    n1,n2,n3 = 120,221,220 # shallow incoherent
    s1,s2,s3 = Sampling(n1,1.0,0.0),Sampling(n2,1.0,0.0),Sampling(n3,1.0,0.0)
    x = readImage(n1,n2,n3,fileName)
    x = copy(90,n2,n3,15,0,0,x)
    writeImage(x,fileName)
  subset("bfl")
  subset("bflt")
  subset("bfp")
  subset("bft")
  subset("bg")
  subset("bgs")
  #subset("bh")
  subset("bt1")
  subset("bt2")
  subset("bt3")

def goSurfing():
  global dataDir,dataPre
  dataDir = "/data/seis/f3d/faults/"
  dataPre = "b"
  n1,n2,n3 = 90,221,220 # deeper coherent
  s1,s2,s3 = Sampling(n1,1.0,0.0),Sampling(n2,1.0,0.0),Sampling(n3,1.0,0.0)
  g = readImage(n1,n2,n3,"g")
  #h = readImage(n1,n2,n3,"h")
  s1,s2,s3,gignore = imageF3d()
  gs = readImage(n1,n2,n3,"gs")
  fl = readImage(n1,n2,n3,"fl")
  fp = readImage(n1,n2,n3,"fp")
  ft = readImage(n1,n2,n3,"ft")
  print "s1:",s1.count,s1.delta,s1.first
  print "s2:",s2.count,s2.delta,s2.first
  print "s3:",s3.count,s3.delta,s3.first
  fs = FaultSurfer3([fl,fp,ft])
  fs.setThreshold(0.5)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,2000)
  fl = fp = ft = None
  p2 = readImage(n1,n2,n3,"p2")
  p3 = readImage(n1,n2,n3,"p3")
  s = fs.findShifts(20.0,surfs,gs,p2,p3)
  print "s: min =",min(s)," max =",max(s)
  #t1,t2,t3 = fs.findThrows(-0.12345,surfs)
  #plot3(g,surfs=surfs)
  #plot3(g,surfs=surfs,smax=-3.75)
  plot3(g,surfs=surfs,smax= 3.75)
  #plot3(h,surfs=surfs,smax= 3.75)
  #plot3(g,s,-10,10,gmap=bwrFill(0.7))
  #plot3(g,t1,-10.0,10.0,gmap=bwrFill(0.7))
  #plot3(g,t2,-0.50,0.50,gmap=bwrFill(0.7))
  #plot3(g,t3,-0.50,0.50,gmap=bwrFill(0.7))
  plot3(g,s,-5,5,gmap=bwrNotch(1.0))
  plot3(g)

def goSubsets():
  def subset(s1,s2,s3,g):
    n1,j1 = 90,130 # deeper coherent
    #n1,j1 = 120,0 # shallow incoherent
    s1 = Sampling(n1,s1.delta,s1.first+j1*s1.delta)
    g = copy(n1,s2.count,s3.count,j1,0,0,g)
    return s1,s2,s3,g
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  gs = readImage(n1,n2,n3,"gs")
  fl = readImage(n1,n2,n3,"fl")
  fp = readImage(n1,n2,n3,"fp")
  ft = readImage(n1,n2,n3,"ft")
  flt = readImage(n1,n2,n3,"flt")
  s1,s2,s3,g  = subset(s1,s2,s3,g)
  s1,s2,s3,gs = subset(s1,s2,s3,gs)
  s1,s2,s3,fl = subset(s1,s2,s3,fl)
  s1,s2,s3,fp = subset(s1,s2,s3,fp)
  s1,s2,s3,ft = subset(s1,s2,s3,ft)
  s1,s2,s3,flt = subset(s1,s2,s3,flt)
  global dataDir,dataPre
  dataDir = "/data/seis/f3d/faults/"
  dataPre = "a"
  writeImage(gs,"gs")
  writeImage(fl,"fl")
  writeImage(fp,"fp")
  writeImage(ft,"ft")
  writeImage(flt,"flt")

def goUnfault():
  global dataDir,dataPre
  dataDir = "/data/seis/f3d/faults/"
  dataPre = ""
  n1,n2,n3 = 90,221,220
  #n1,n2,n3 = 120,221,220
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  s3 = Sampling(n3,1.0,0.0)
  g = readImage(n1,n2,n3,"ag")
  h = readImage(n1,n2,n3,"ah")
  #g = readImage(n1,n2,n3,"bg")
  #h = readImage(n1,n2,n3,"bh")
  plot3(g,h=h)
  #g = slog(g)
  #q1 = zerofloat(n1,n2,n3)
  #t1 = readImage(n1,n2,n3,"at1")
  #plot3(g,t1,-10.0,10.0,gmap=bwrFill(0.7))
  #st1,sx1,sx2,sx3 = SimpleGridder3.getGriddedSamples(-0.12345,s1,s2,s3,t1)
  #bg = BlendedGridder3()
  #d1 = bg.gridNearest(-0.12345,t1)
  #bg.gridBlended(d1,t1,q1)
  #sg = SibsonGridder3(st1,sx1,sx2,sx3)
  #q1 = sg.grid(s1,s2,s3)
  #plot3(g,q1,-10.0,10.0,gmap=bwrFill(0.7))

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

def goShifts():
  s1,s2,s3,g = imageF3d()
  f = readImage(n1,n2,n3,"fl")
  p = readImage(n1,n2,n3,"fp")
  t = readImage(n1,n2,n3,"ft")
  #f,p,t = FaultScanner3.fakeFpt(s1.count,s2.count,s3.count)
  #plot3(sub(1.0,f),None,0,1)
  g = slog(g)
  #g = readImage(n1,n2,n3,"gs")
  #plot3(sub(1.0,f),None,0,1)
  #plot3(g,f,0,1)
  #plot3(g,p)
  faults = FaultScanner3.findFaults([f,p,t],1000)
  print "faults: count =",faults.count
  f = faults.getLikelihoods()
  plot3(g,f,0,1)
  #s = faults.getShifts()
  #plot3(g,s)
  return

def goSmooth():
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  fl = readImage(n1,n2,n3,"flt")
  #p2 = readImage(n1,n2,n3,"p2")
  #p3 = readImage(n1,n2,n3,"p3")
  #gs = FaultScanner3.smooth(16.0,p2,p3,fl,g)
  #writeImage(gs,"gs")
  gs = readImage(n1,n2,n3,"gs")
  #plot3(g,fl,0,1)
  plot3(g)
  plot3(gs)
  plot3(sub(1.0,pow(fl,4)))

def goThin():
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  f = readImage(n1,n2,n3,"fl")
  p = readImage(n1,n2,n3,"fp")
  t = readImage(n1,n2,n3,"ft")
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
  plot3(g,f,0,1)
  plot3(g,p,-90,90)
  plot3(g,t,-15,15)

def goScan():
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  fse = FaultSemblance()
  #g = fse.taper(10,g)
  sn = readImage(n1,n2,n3,"sn")
  sd = readImage(n1,n2,n3,"sd")
  sigmaPhi,sigmaTheta = 4,20
  fsc = FaultScanner3(sigmaPhi,sigmaTheta,[sn,sd])
  print "scanning ..."
  sw = Stopwatch()
  sw.restart()
  f,p,t = fsc.scan(-90,90,-15,15)
  sw.stop()
  print "time =",sw.time(),", f min =",min(f)," max =",max(f)
  writeImage(f,"fl")
  writeImage(p,"fp")
  writeImage(t,"ft")
  plot3(g,f,0,1)
  plot3(g,p,-90,90)
  plot3(g,t,-15,15)

def goSemblance():
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  print "reading slopes ..."
  p2 = readImage(n1,n2,n3,"p2")
  p3 = readImage(n1,n2,n3,"p3")
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
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  n1,n2,n3 = len(g[0][0]),len(g[0]),len(g)
  fse = FaultSemblance()
  ref = RecursiveExponentialFilter(4)
  p2 = readImage(n1,n2,n3,"p2") # semblance with estimated slopes
  p3 = readImage(n1,n2,n3,"p3")
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
  s1,s2,s3,g = imageF3d()
  g = slog(g)
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
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  fse = FaultSemblance()
  #g = fse.taper(10,g)
  sn = readImage(n1,n2,n3,"sn")
  sd = readImage(n1,n2,n3,"sd")
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

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))
 
#############################################################################
# data read/write

dataDir,dataPre = "",""

def readImage(n1,n2,n3,fileName):
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+dataPre+fileName+".dat")
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(f,fileName):
  aos = ArrayOutputStream(dataDir+dataPre+fileName+".dat")
  aos.writeFloats(f)
  aos.close()

def imageF3d():
  global dataDir,dataPre
  dataDir = "/data/seis/f3d/faults/"
  dataPre = "s1"
  n1,n2,n3 = 462,951,591
  #j1,j2,j3 = 240,0,0
  #m1,m2,m3 = 222,440,440
  j1,j2,j3 = 240, 50,100
  m1,m2,m3 = 222,221,220
  d1,d2,d3 = 0.004,0.025,0.025
  f1,f2,f3 = 0.004,0.000,0.000
  firstTime = False
  if firstTime:
    af = ArrayFile(dataDir+"f3d.dat","r")
    x = zerofloat(m1,m2,m3)
    for i3 in range(m3):
      for i2 in range(m2):
        af.seek(4*(j1+n1*(i2+j2+n2*(i3+j3))))
        af.readFloats(x[i3][i2])
    af.close()
    writeImage(x,"g")
  n1,n2,n3 = m1,m2,m3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  x = readImage(n1,n2,n3,"g")
  return s1,s2,s3,x

#############################################################################
# plotting

def plot2(f,x23p=None,x23l=None,fmin=None,fmax=None,
          label=None,title=None,png=None):
  n1 = len(f[0])
  n2 = len(f)
  panel = panel2()
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  if title:
    panel.setTitle(title)
  panel.setColorBarWidthMinimum(100)
  panel.setLimits(0,0,n1-1,n2-1)
  #panel.setLimits(0,0,100,100)
  #panel.setLimits(100,160,130,190)
  pv = panel.addPixels(f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  if fmin==None: fmin = min(f)
  if fmax==None: fmax = max(f)
  pv.setClips(fmin,fmax)
  if x23p:
    x2p,x3p = x23p
    pv = panel.addPoints(x2p,x3p)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(5.0)
  if x23l:
    x2l,x3l = x23l
    pv = panel.addPoints(x2l,x3l)
    pv.setLineWidth(2.0)
  frame2(panel,png)

def panel2():
  #panel = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)
  return panel

def frame2(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.8)
  #frame.setSize(1290,777)
  #frame.setSize(1490,977)
  frame.setSize(980,890)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+"/"+png+".png")
  return frame

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
