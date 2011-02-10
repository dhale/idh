import sys
#from math import *
from java.awt import *
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

#pngDir = None
pngDir = "png/blf/"

#############################################################################

def main(args):
  #goImage()
  goFilter()
  #goFilterImpulses()
  #goFilterRandom()

def goFilter():
  #lim = (2.1,1.4,7.9,1.85)
  lim = (2.3,0.98,7.7,1.42)
  #x,s1,s2,clip = getImageF3d(); png = "f3d_"
  x,s1,s2,clip = getImageTpd(); png = "tpd_"
  plot(x,s1,s2,clip,png=png+"input")
  #plot(x,s1,s2,clip,limits=lim,png=png+"inputz")
  sigmaS = 16.0
  sigmaR = qqd(x); print "sigmaS =",sigmaS," sigmaR =",sigmaR
  doFilter(sigmaS,sigmaR,False,False,x,s1,s2,clip,png=png+"blf")
  doFilter(sigmaS,sigmaR,True,False,x,s1,s2,clip,png=png+"blft")
  #doFilter(sigmaS,sigmaR,True,False,x,s1,s2,clip,limits=lim,png=png+"blftz")
  doFilter(sigmaS,sigmaR,True,True,x,s1,s2,clip,png=png+"blftc")
  doFilter(sigmaS,100.0*sigmaR,True,False,x,s1,s2,clip,png=png+"lsft")
  doFilter(sigmaS,100.0*sigmaR,True,True,x,s1,s2,clip,png=png+"lsftc")

def doFilter(sigmaS,sigmaR,useT,useC,x,s1,s2,clip,limits=None,png=None):
  if useT:
    t = imageTensors(2.0,x)
    if useC:
      t.scale(coherence(2.0,t,x))
  else:
    t = None
  y = bilateralFilter(sigmaS,sigmaR,t,x)
  plot(y,s1,s2,clip,limits=limits,png=png)
  plot(sub(x,y),s1,s2,clip*0.251,limits=limits,png=png+"d")

def goFilterImpulses():
  xa,s1,s2,clip = getImage()
  xb = makeImpulses(12,len(xa[0]),len(xa))
  t = imageTensors(2.0,xa)
  #t.scale(coherence(2.0,t,xa))
  plot(xa,s1,s2,1.3)
  #plot(xb,s1,s2,0.9)
  for sigmaR in [0.1,1.0,100.0]:
  #for sigmaR in [100.0]:
    y = like(xa)
    bf = BilateralFilter(30.0,sigmaR)
    bf.setType(BilateralFilter.Type.TUKEY)
    bf.applyAB(t,xa,xb,y)
    y = smoothS(y)
    #plot(y,s1,s2,0.5*max(y))
    plot(y,s1,s2,0.1)

def goFilterRandom():
  xa,s1,s2,clip = getImage()
  plot(xa,s1,s2,clip)
  n1,n2 = len(xa[0]),len(xa)
  xb = makeRandom(n1,n2)
  t = imageTensors(2.0,xa)
  #t.scale(coherence(2.0,t,xa))
  #plot(xa,s1,s2,1.3)
  #plot(xb,s1,s2,0.5)
  xqqd = qqd(x)
  for sigmaR in [xqqd,100*xqqd]:
    y = like(xa)
    bf = BilateralFilter(30.0,sigmaR)
    bf.setType(BilateralFilter.Type.TUKEY)
    bf.applyAB(t,xa,xb,y)
    plot(y,s1,s2,0.1)
    bf.apply(t,xa,y)
    plot(y,s1,s2,clip)

def goImage():
  x,s1,s2,clip = getImage()
  plot(x,s1,s2,clip)

def getImage():
  return getImageF3d()
  #return getImageTpd()
  #return getImageSyn()

def getImageF3d():
  n1,n2 = 462,951
  d1,d2 = 0.004,0.025
  f1,f2 = 0.004,0.000
  fileName = "/data/seis/f3d/f3d75.dat"
  x = readImage(fileName,n1,n2)
  subset = True
  if subset:
    j1,j2 = 240,0
    n1,n2 = n1-j1,440
    f1,f2 = f1+j1*d1,f2+j2*d2
    x = copy(n1,n2,j1,j2,x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,5.0

def getImageTpd():
  n1,n2 = 251,357
  d1,d2 = 0.004,0.025
  f1,f2 = 0.500,0.000
  fileName = "/data/seis/tp/csm/oldslices/tp73.dat"
  x = readImage(fileName,n1,n2)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,2.0

def getImageSyn():
  n1,n2 = 251,357
  d1,d2 = 1.0,1.0
  f1,f2 = 0.0,0.0
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  x = zerofloat(n1,n2)
  for i2 in range(n2):
    if i2<n2/2:
      x[i2] = sin(rampfloat(0.00,0.1,n1))
    else:
      x[i2] = sin(rampfloat(3.14,0.1,n1))
    x[i2] = mul(1.4,x[i2])
  return x,s1,s2,1.4

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def imageTensors(sigma,x):
  """ 
  Returns tensors for guiding filters along features in specified image.
  """
  n1,n2 = len(x[0]),len(x)
  lof = LocalOrientFilter(2.0*sigma,2.0)
  lof.setGradientSmoothing(sigma)
  t = lof.applyForTensors(x) # structure tensors
  t.invertStructure(0.0,4.0) # inverted with ev = 1.0, eu = small
  #eu = fillfloat(0.0001,n1,n2)
  #ev = fillfloat(1.0000,n1,n2)
  #t.setEigenvalues(eu,ev)
  return t

def bilateralFilter(sigmaS,sigmaX,t,x):
  y = like(x)
  bf = BilateralFilter(sigmaS,sigmaX)
  if t:
    bf.apply(t,x,y)
  else:
    bf.apply(x,y)
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
    lsf = LocalSmoothingFilter(0.001,int(10*sigma))
    y = copy(x)
    #lsf.applySmoothS(y,y)
    #lsf.applySmoothL(kmax,y,y)
    lsf.apply(t,0.5*sigma*sigma,y,z)
  return z

def coherence(sigma,t,x):
  s = semblance(sigma,t,x) # structure-oriented semblance s
  s = pow(s,16.0) # make smaller sembances in [0,1] much smaller
  return s

def semblance(sigma,t,s):
  lsf = LocalSemblanceFilter(int(sigma),4*int(sigma))
  return lsf.semblance(LocalSemblanceFilter.Direction2.V,t,s)

def qqd(x):
  return 0.5*(Quantiler.estimate(0.75,x)-Quantiler.estimate(0.25,x))

def like(x):
  n1,n2 = len(x[0]),len(x)
  return zerofloat(n1,n2)

def makeImpulses(ni,n1,n2):
  x = zerofloat(n1,n2)
  ns = max(n1/ni,n2/ni)
  m1 = (n1-1)/ns
  m2 = (n2-1)/ns
  j1 = (n1-1-(m1-1)*ns)/2
  j2 = (n2-1-(m2-1)*ns)/2
  for i2 in range(j2,n2,ns):
    for i1 in range(j1,n1,ns):
      x[i2][i1] = 1.0
  return x

def makeRandom(n1,n2):
  x = mul(2.0,sub(randfloat(n1,n2),0.5))
  return smooth(1.0,None,x)

#############################################################################
# plot

def plot(f,s1,s2,clip=0.0,t=None,cbar="",limits=None,png=None):
  n1,n2 = len(f[0]),len(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setFontSizeForPrint(8.0,240)
  sp.setSize(1020,750)
  #sp.setSize(1090,760)
  #sp.setSize(1090,800)
  if limits:
    sp.setHInterval(1.0)
    sp.setVInterval(0.1)
  else:
    sp.setHInterval(2.0)
    sp.setVInterval(0.2)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Time (s)")
  if limits:
    sp.setLimits(limits[0],limits[1],limits[2],limits[3])
  if cbar!=None:
    if len(cbar)>0:
      cbar = sp.addColorBar(cbar)
    else:
      cbar = sp.addColorBar()
    if clip==2.0:
      cbar.setInterval(1.0)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(s1,s2,f)
  if clip!=0.0:
    if clip==1.0:
      pv.setClips(0.0,clip)
    else:
      pv.setClips(-clip,clip)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if clip==1.0:
    #pv.setColorModel(ColorMap.JET)
    pv.setColorModel(ColorMap.GRAY)
  else:
    #pv.setColorModel(ColorMap.GRAY_YELLOW_RED)
    pv.setColorModel(ColorMap.GRAY)
  if t:
    tv = TensorsView(s1,s2,t)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(30)
    #tv.setScale(3)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if png and pngDir:
    sp.paintToPng(720,3.3,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
