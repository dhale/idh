import sys
from math import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from ldf import BilateralFilter

gauss = BilateralFilter.Type.GAUSS
huber = BilateralFilter.Type.HUBER
tukey = BilateralFilter.Type.TUKEY

#############################################################################

def main(args):
  goFilter()
  goSmooth()
  #goSemblance()
  #goNormalize()

def getImage():
  n1,n2 = 462,951
  fileName = "/data/seis/f3d/f3d75.dat"
  return readImage(fileName,n1,n2),6.0

def goFilter():
  x,xclip = getImage()
  plot(x,xclip)
  t = imageTensors(2.0,0.0001,x)
  #c = semblance(2.0,t,x)
  #c = pow(c,8.0)
  #plot(c,1.0)
  #eu = mul(0.0001,c)
  #ev = mul(1.0000,c)
  #t.setEigenvalues(eu,ev)
  s = localScale(3.0,x)
  x = div(x,s)
  #plot(x,2.0)
  y = bilateralFilter(30.0,1.0,t,x)
  #plot(y,2.0)
  y = mul(s,y)
  plot(y,xclip)

def goSmooth():
  x,xclip = getImage()
  plot(x,xclip)
  t = imageTensors(2.0,0.0001,x)
  c = semblance(2.0,t,x)
  c = pow(c,8.0)
  #plot(c,1.0)
  eu = mul(0.0001,c)
  ev = mul(1.0000,c)
  t.setEigenvalues(eu,ev)
  y = smooth(30.0,t,x)
  plot(y,xclip)

def goSemblance():
  x,xclip = getImage()
  plot(x,xclip)
  t = imageTensors(2.0,0.001,x)
  u,s,y = normalize(80.0,t,x)
  s = semblance(2.0,t,y)
  plot(s,1.0)

def goNormalize():
  x,xclip = getImage()
  t = imageTensors(2.0,0.001,x)
  u,s,y = normalize(80.0,t,x)
  print "x: min =",min(x)," max =",max(x)
  print "u: min =",min(u)," max =",max(u)
  print "s: min =",min(s)," max =",max(s)
  print "y: min =",min(y)," max =",max(y)
  plot(x,xclip)
  plot(u,xclip)
  plot(s,xclip)
  plot(y,2.0)

def goEdges():
  n1,n2 = 251,357
  clip = 4.5
  fileName = "/data/seis/tp/csm/oldslices/tp73.dat"
  x = readImage(fileName,n1,n2)
  t = makeImageTensors(x)
  plot(x,clip)
  sigma1 = 4.0; c1 = 0.5*sigma1*sigma1
  sigma2 = 8.0; c2 = 0.5*sigma2*sigma2
  sigma3 = 2.0; c3 = 0.5*sigma3*sigma3
  lsf = LocalSmoothingFilter()
  rgf = RecursiveGaussianFilter(sigma3)
  y1 = copy(x)
  y2 = copy(x)
  lsf.apply(t,c1,x,y1)
  lsf.apply(t,c2,x,y2)
  #plot(y1)
  #plot(y2)
  y = sub(y2,y1)
  plot(y)
  z = mul(x,y)
  plot(z)
  rgf.apply00(z,z)
  plot(z)

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def imageTensors(sigma,eu,x):
  n1,n2 = len(x[0]),len(x)
  lof = LocalOrientFilter(2.0*sigma,2.0)
  lof.setGradientSmoothing(sigma)
  t = lof.applyForTensors(x)
  t.setEigenvalues(eu,1.000)
  return t

def bilateralFilter(sigmaS,sigmaX,t,x):
  y = like(x)
  bf = BilateralFilter(sigmaS,sigmaX)
  bf.apply(t,x,y)
  return y

def smoothS(x):
  y = like(x)
  lsf = LocalSmoothingFilter()
  lsf.applySmoothS(x,y)
  return y

def smooth(sigma,t,x):
  z = copy(x)
  if t==None:
    rgf = RecursiveGaussianFilter(sigma)
    rgf.apply00(x,z)
  else:
    ldk = LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D71)
    lsf = LocalSmoothingFilter(0.001,1000,ldk)
    y = copy(x)
    #lsf.applySmoothS(y,y)
    #lsf.applySmoothL(kmax,y,y)
    lsf.apply(t,0.5*sigma*sigma,y,z)
  return z

def localScale(sigma,x):
  s = mul(x,x)
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply00(copy(s),s)
  smax = max(s)
  smin = 1.0e-6*smax
  clip(smin,smax,s,s)
  sqrt(s,s)
  return s

def localMuSigma(sigma,t,x):
  u = smooth(sigma,t,x) # u = <x>
  zero(u)
  s = sub(x,u) # s = x-u
  mul(s,s,s) # s = (x-u)^2
  s = smooth(sigma,t,s) # s = <(x-<x>)^2>
  s = smooth(2.0,None,s)
  clip(0,max(s),s,s) # 0 <= s
  sqrt(s,s) # s = sqrt(<(x-<x>)^2>)
  return u,s

def clipSmall(small,s):
  smax = max(s)
  smin = small*smax
  return clip(smin,smax,s)

def normalize(sigma,t,x):
  u,s = localMuSigma(sigma,t,x)
  s = clipSmall(0.0001,s)
  return u,s,div(x,s)
  #return u,s,div(sub(x,u),s)

def powerGain(x,p):
  return mul(sgn(x),pow(abs(x),p))

def semblance(sigma,t,s):
  lsf = LocalSemblanceFilter(int(sigma),2*int(sigma))
  return lsf.semblance(LocalSemblanceFilter.Direction2.V,t,s)

def like(x):
  n1,n2 = len(x[0]),len(x)
  return zerofloat(n1,n2)

#############################################################################
# plot

pngDir = "./png"
#pngDir = None

def plot(f,clip=0.0,png=None):
  n1,n2 = len(f[0]),len(f)
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(s1,s2,f)
  if clip!=0.0:
    if clip==1.0:
      pv.setClips(0.0,clip)
    else:
      pv.setClips(-clip,clip)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #if clip==1.0:
  #  pv.setColorModel(ColorMap.JET)
  #sp.setFontSizeForSlide(1.0,1.0)
  sp.setSize(1192,863)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(100,6,pngDir+"/"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
