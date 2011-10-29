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
  goTrack()

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

def plot3(f,g=None,gmin=None,gmax=None,xyz=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame()
  if g==None:
    ipg = sf.addImagePanels(f)
  else:
    ipg = ImagePanelGroup2(f,g)
    sf.world.addChild(ipg)
  if xyz!=None:
    pg = PointGroup(xyz)
    ps = PointState()
    ps.setSize(3.0)
    ss = StateSet()
    ss.add(ps)
    pg.setStates(ss)
    sf.world.addChild(pg)
  ipg.setSlices(209,12,18)
  sf.setSize(1200,1100)
  sf.orbitView.setScale(2.5)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
