#############################################################################
# Dynamic time warping for 3D images

from imports import *
from util import FakeData

#############################################################################

#pngDir = "./png"
pngDir = None

#seed = abs(Random().nextInt()/1000000)
seed = 588
print "seed =",seed

smax = 12 # max shift for synthetic test
nrms = 0.0 # rms noise/signal ratio
strainMax1 = 0.25 # not less than smax*2*pi/n1
strainMax2 = 0.20 # not less than smax*pi/n2
strainMax3 = 0.20 # not less than smax*pi/n3
n1,n2,n3 = 201,201,201 # numbers of samples
d1,d2,d3 = 0.004,0.025,0.025 # sampling intervals
f1,f2,f3 = 0.000,0.000,0.000 # first samples
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
label1,label2,label3 = "Time (s)","Inline (km)","Crossline (km)"
dataDir = "/data/seis/fake/"

def main(args):
  #goFakeImages()
  goFakeShifts()

def goFakeImages():
  f,g,s = makeFakeImages(smax,nrms)
  print "s: min =",min(s),"max =",max(s)
  writeImage(dataDir+"fakef.dat",f)
  writeImage(dataDir+"fakeg.dat",g)
  writeImage(dataDir+"fakes.dat",s)
  #show(f)
  #show(g)
  #show(s)

def goFakeShifts():
  f = readImage(dataDir+"fakef.dat",n1,n2,n3)
  g = readImage(dataDir+"fakeg.dat",n1,n2,n3)
  s = readImage(dataDir+"fakes.dat",n1,n2,n3)
  #f,g,s = makeFakeImages(smax,nrms)
  print "s: min =",min(s),"max =",max(s)
  esmooth = 2
  usmooth = 1.0
  mlag = 4+smax
  dw = DynamicWarping(-mlag,mlag)
  dw.setStrainMax(strainMax1,strainMax2,strainMax3)
  dw.setErrorSmoothing(esmooth)
  dw.setShiftSmoothing(usmooth)
  dw.setWindowSizeAndOverlap(51,51,0.5,0.5)
  u = dw.findShifts(f,g)
  h = dw.applyShifts(u,g)
  print "s: min =",min(s),"max =",max(s)
  print "u: min =",min(u),"max =",max(u)
  show(s)
  show(u)
  show(g)
  show(f)
  show(h)

def makeFakeImages(smax,nrms):
  f = FakeData.seismic3d2010A(n1,n2,n3,20.0,10.0,30.0,0.5,0.0);
  w = Warp3.sinusoid(0.5*smax,0.0,0.0,0.5*smax,0.0,0.0,n1,n2,n3)
  #w = Warp3.constant(smax,0.0,0.0,n1,n2,n3)
  g = w.warp1(f)
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  s = w.u1x()
  return f,g,s

#############################################################################
# utilities

def readImage(fileName,n1,n2,n3):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(x)
  aos.close()

def addRickerWavelet(fpeak,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      Conv.conv(nh,-ih,h,n1,0,f[i3][i2],n1,0,g[i3][i2])
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,f,seed=0):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  g = mul(2.0,sub(randfloat(r,n1,n2,n3),0.5))
  g = addRickerWavelet(0.125,g) # bandlimited in time
  rgf = RecursiveGaussianFilter(1.0)
  rgf.applyX0X(g,g)
  rgf.applyXX0(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2/n3)
  grms = sqrt(sum(mul(g,g))/n1/n2/n3)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

def show(f):
  frame = SimpleFrame()
  ip = frame.addImagePanels(f)
  ip.setSlices(n1-1,n2/2,n3/2)
  frame.getOrbitView().setScale(2.0)
  frame.setSize(900,900)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
