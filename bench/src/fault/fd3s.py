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
  #goTrack()
  #goShifts()
  goSurfing()
  #goSurfingFake()

def goSurfingFake():
  n1,n2,n3 = 101,102,103
  #n1,n2,n3 = 51,52,53
  #n1,n2,n3 = 41,42,43
  #n1,n2,n3 = 11,12,13
  g = sub(randfloat(n1,n2,n3),0.5)
  f,p,t = Util.fakeFpt(n1,n2,n3)
  fs = FaultSurfer3([f,p,t])
  fs.setThreshold(0.5)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,1000)
  plot3(g,f,0,1,surfs=surfs)

def goSurfing():
  def subset(s1,s2,s3,g):
    n1,j1 = 70,130 # deeper coherent
    #n1,j1 = 130,0 # shallow incoherent
    s1 = Sampling(n1,s1.delta,s1.first+j1*s1.delta)
    g = copy(n1,s2.count,s3.count,j1,0,0,g)
    return s1,s2,s3,g
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  f = readImage(s1,s2,s3,"fl")
  p = readImage(s1,s2,s3,"fp")
  t = readImage(s1,s2,s3,"ft")
  s1,s2,s3,g = subset(s1,s2,s3,g)
  s1,s2,s3,f = subset(s1,s2,s3,f)
  s1,s2,s3,p = subset(s1,s2,s3,p)
  s1,s2,s3,t = subset(s1,s2,s3,t)
  fs = FaultSurfer3([f,p,t])
  fs.setThreshold(0.5)
  quads = fs.findQuads()
  quads = fs.linkQuads(quads)
  surfs = fs.findSurfs(quads)
  surfs = fs.getSurfsWithSize(surfs,1000)
  plot3(g,surfs=surfs)
  #plot3(g)

def goShifts():
  s1,s2,s3,g = imageF3d()
  f = readImage(s1,s2,s3,"fl")
  p = readImage(s1,s2,s3,"fp")
  t = readImage(s1,s2,s3,"ft")
  #f,p,t = FaultScanner3.fakeFpt(s1.count,s2.count,s3.count)
  #plot3(sub(1.0,f),None,0,1)
  g = slog(g)
  #g = readImage(s1,s2,s3,"gs")
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

def goTrack():
  s1,s2,s3,g = imageF3d()
  f = readImage(s1,s2,s3,"flt")
  g = slog(g)
  #g = readImage(s1,s2,s3,"gs")
  #plot3(sub(1.0,f),None,0,1)
  #return
  p = readImage(s1,s2,s3,"fp")
  t = readImage(s1,s2,s3,"ft")
  ft = FaultTracker3([f,p,t])
  #i1,i2,i3 = 189,184,111
  #i1,i2,i3 = 165,163,121
  #i1,i2,i3 = 182,146,108
  i1,i2,i3 = 209,99,29
  #i1,i2,i3 = 195,54,37
  ft.setThreshold(0.1)
  kf = ft.track(i1,i2,i3)
  print "number of points =",len(kf[0])
  xyz = ft.xyz(kf)
  plot3(g,f,0,1,xyz)
  return
  gfm,gfp = sampleFault(kf,g)
  writeImage(gfm,"gfm")
  writeImage(gfp,"gfp")
  SimplePlot.asPixels(gfm)
  SimplePlot.asPixels(gfp)

def sampleFault(kf,g):
  n1,n2,n3 = len(g[0][0]),len(g[0]),len(g)
  k1,k2,k3 = kf[0],kf[1],kf[2]
  nk = len(k1)
  #j2 = fillint(n2,n1,n3)
  j2 = fillint(-1,n1,n3)
  for ik in range(nk):
    i1,i2,i3 = k1[ik],k2[ik],k3[ik]
    #if i2>j2[i3][i1]:
    if i2>j2[i3][i1]:
      j2[i3][i1] = i2
  gfm = zerofloat(n1,n3)
  gfp = zerofloat(n1,n3)
  for i3 in range(n3):
    for i1 in range(n1):
      i2 = j2[i3][i1]
      #if i2<n2:
      if i2>=0:
        i2m = max(0,i2-2)
        i2p = min(n2-1,i2+2)
        gfm[i3][i1] = g[i3][i2m][i1]
        gfp[i3][i1] = g[i3][i2p][i1]
  return gfm,gfp

def goSmooth():
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  fl = readImage(s1,s2,s3,"flt")
  #p2 = readImage(s1,s2,s3,"p2")
  #p3 = readImage(s1,s2,s3,"p3")
  #gs = FaultScanner3.smooth(16.0,p2,p3,fl,g)
  #writeImage(gs,"gs")
  gs = readImage(s1,s2,s3,"gs")
  #plot3(g,fl,0,1)
  plot3(g)
  plot3(gs)
  plot3(sub(1.0,pow(fl,4)))

def goThin():
  s1,s2,s3,g = imageF3d()
  g = slog(g)
  f = readImage(s1,s2,s3,"fl")
  p = readImage(s1,s2,s3,"fp")
  t = readImage(s1,s2,s3,"ft")
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
  sn = readImage(s1,s2,s3,"sn")
  sd = readImage(s1,s2,s3,"sd")
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
  p2 = readImage(s1,s2,s3,"p2")
  p3 = readImage(s1,s2,s3,"p3")
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
  p2 = readImage(s1,s2,s3,"p2") # semblance with estimated slopes
  p3 = readImage(s1,s2,s3,"p3")
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
  sn = readImage(s1,s2,s3,"sn")
  sd = readImage(s1,s2,s3,"sd")
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

def readImage(s1,s2,s3,fileName):
  n1,n2,n3 = s1.count,s2.count,s3.count
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
  x = readImage(s1,s2,s3,"g")
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

def plot3(f,g=None,gmin=None,gmax=None,xyz=None,surfs=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame()
  if g==None:
    ipg = sf.addImagePanels(f)
  else:
    ipg = ImagePanelGroup2(f,g)
    if gmin and gmax:
      ipg.setClips2(gmin,gmax)
    if gmin==0.0 and gmax==1.0:
      updateColorModel2(ipg,0.8)
    sf.world.addChild(ipg)
  if xyz:
    pg = PointGroup(xyz)
    ss = StateSet()
    ps = PointState()
    ps.setSize(3.0)
    ss.add(ps)
    pg.setStates(ss)
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
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    for surf in surfs:
      xyz,uvw,rgb = surf.getXyzUvwRgb()
      #qg = QuadGroup(False,xyz,rgb)
      qg = QuadGroup(True,xyz,rgb)
      #qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
    sf.world.addChild(sg)
  #ipg.setSlices(209,12,18)
  ipg.setSlices(200,0,0)
  sf.setSize(1200,1100)
  sf.setWorldSphere(n3/2,n2/2,n1/2,0.5*sqrt(n1*n1+n2*n2+n3*n3))
  sf.orbitView.setScale(1.5)

def updateColorModel2(ipg,alpha):
  n = 256
  r = zerobyte(n)
  g = zerobyte(n)
  b = zerobyte(n)
  a = zerobyte(n)
  icm = ipg.getColorModel2()
  icm.getReds(r)
  icm.getGreens(g)
  icm.getBlues(b)
  for i in range(n):
    ai = int(255.0*alpha*i/n)
    if ai>127:
      ai -= 256
    a[i] = ai
  #if alpha<1.0:
    #r[n/2] = r[n/2-1] = -1
    #g[n/2] = g[n/2-1] = -1
    #b[n/2] = b[n/2-1] = -1
    #a[n/2  ] = a[n/2-1] = 0
    #a[n/2+1] = a[n/2-2] = 0
    #a[n/2+2] = a[n/2-3] = 0
    #a[0] = a[1] = 0
  icm = IndexColorModel(8,n,r,g,b,a)
  ipg.setColorModel2(icm)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
