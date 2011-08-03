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

from het import *

pngDir = None
#pngDir = "png/het/"

#############################################################################

def main(args):
  goShifts()
  #goPeakF()
  #goImage()
  #goHetMul()
  #goAmpPhase()

def goPeakF():
  x,s1,s2,clip = getImage()
  plot(x,s1,s2,clip)
  lbf = LocalBurgFilter()
  sigma1 = 16
  for sigma2 in [0,1,2,4,8,16]:
    c = lbf.getCoefficients(sigma1,sigma2,2,x)
    pf = lbf.getPeakFrequencies(c)
    pf = mul(pf,1/s1.delta)
    plot(pf,s1,s2)

def goShifts():
  x,s1,s2,xclip = getImage()
  plot(x,s1,s2,xclip)
  n1,n2 = s1.count,s2.count
  lbf = LocalBurgFilter()
  sigma1,sigma2 = 16,2
  c = lbf.getCoefficients(sigma1,sigma2,2,x)
  pf = lbf.getPeakFrequencies(c)
  pf = sum(pf)/n1/n2
  print "peak frequency average =",pf/s1.delta,"Hz"
  ref = RecursiveExponentialFilter(16,0)
  #ref.setExtrapolation(RecursiveExponentialFilter.Extrapolation.CONSTANT)
  pfList = [pf-0.04,pf,pf+0.04]
  for pf in pfList:
    p = rampfloat(0.0,2.0*PI*pf,0.0,n1,n2)
    c = cos(p)
    s = sin(p)
    cx = mul(c,x)
    sx = mul(s,x)
    #plot(cx,s1,s2,xclip)
    #plot(sx,s1,s2,xclip)
    ref.apply(cx,cx)
    ref.apply(sx,sx)
    #plot(cx,s1,s2,0.2*xclip)
    #plot(sx,s1,s2,0.2*xclip)
    cx = cmplx(cx,sx)
    #ax = carg(cx)
    #plot(ax,s1,s2)
    cy = like(cx)
    lag2 = 5
    ccopy(n1,n2-lag2,0,lag2,cx,0,0,cy)
    for i2 in range(n2-lag2,n2):
      ccopy(cx[n2-1],cy[i2])
    cxy = cmul(cx,cconj(cy))
    axy = carg(cxy)
    plot(axy,s1,s2)
  return
  lag2 = 1
  #plot(y,s1,s2,xclip)
  cx = mul(c,x)
  sx = mul(s,x)
  plot(cx,s1,s2,xclip)
  ref = RecursiveExponentialFilter(32,0.0)
  ref.apply(cx,cx)
  ref.apply(sx,sx)
  plot(cx,s1,s2,0.2*xclip)
  cx = cmplx(cx,sx)
  cy = like(cx)
  ccopy(n1,n2-lag2,0,lag2,cx,0,0,cy)
  for i2 in range(n2-lag2,n2):
    ccopy(cx[n2-1],cy[i2])
  cxy = cmul(cx,cconj(cy))
  axy = carg(cxy)
  plot(axy,s1,s2)
  axy = unwrap2(axy)
  plot(axy,s1,s2)
  return
  ax = carg(cx)
  plot(ax,s1,s2)
  ax = unwrap2(ax)
  plot(ax,s1,s2)

def unwrap1(p):
  n = len(p)
  q = copy(p)
  twopi = 2.0*PI
  for i in range(1,n):
    dp = p[i]-p[i-1]
    absdp = abs(dp)
    if (abs(dp+twopi)<absdp):
      dp += twopi
    elif (abs(dp-twopi)<absdp):
      dp -= twopi
    q[i] = q[i-1]+dp
  return q
def unwrap2(p):
  n = len(p)
  q = like(p)
  for i in range(n):
    q[i] = unwrap1(p[i])
  return q

def goAmpPhase():
  x,s1,s2,clip = getImage()
  plot(x,s1,s2,clip)
  a,p = ampPhase2(x)
  c,s = cos(p),sin(p)
  plot(a,s1,s2,clip)
  plot(p,s1,s2)
  plot(c,s1,s2,1.0)
  plot(s,s1,s2,1.0)

def goHetMul():
  x,s1,s2,clip = getImage()
  plot(x,s1,s2,clip)
  fpeak = getPeakFrequency(x)
  print "peak frequency =",fpeak/s1.delta
  n1,n2 = s1.count,s2.count
  shift = -0.00
  #for makeMod in [makeCosineModulator,makeBandModulator,makeBurgModulator]:
  for makeMod in [makeBurgModulator]:
    #g = makeMod(shift,x)
    g = x
    plot(g,s1,s2,0.5*clip)
    y = mul(g,x)
    plot(y,s1,s2,0.5*clip)
    z = smoothingFilter(2.0,y)
    #z = lowPassFilter(fpeak,y)
    plot(z,s1,s2,0.5*clip)

def getPeakFrequency(x):
  n1,n2 = len(x[0]),len(x)
  f = getPeakFrequencies(n1/4.0,x)
  return sum(f)/n1/n2

def getPeakFrequencies(sigma,x):
  n1,n2 = len(x[0]),len(x)
  lbf = LocalBurgFilter()
  f = zerofloat(n1,n2)
  for i2 in range(n2):
    c = lbf.getCoefficients(sigma,2,x[i2])
    f[i2] = lbf.getPeakFrequencies(c)
  return f

def makeCosineModulator(shift,x):
  n1,n2 = len(x[0]),len(x)
  fpeak = getPeakFrequency(x)
  fpeak += shift
  phase = rampfloat(0.0,2.0*PI*fpeak,0.0,n1,n2)
  return cos(phase)

def makeBandModulator(shift,x):
  fpeak = getPeakFrequency(x)
  fpeak += shift
  y = like(x)
  nf = NotchFilter(fpeak,0.9)
  nf.apply1ForwardReverse(x,y)
  return sub(x,y)

def makeBurgModulator(shift,x):
  n1,n2 = len(x[0]),len(x)
  y = like(x)
  lbf = LocalBurgFilter()
  for i2 in range(n2):
    c = lbf.getCoefficients(16,2,x[i2])
    c = lbf.zeroPeakFrequencies(c)
    c = lbf.shiftPeakFrequencies(shift,c)
    t = lbf.applyNotch(0.2,c,x[i2])
    t = reverse(t)
    t = lbf.applyNotch(0.2,c,t)
    t = reverse(t)
    y[i2] = t
  return sub(x,y)

def squareRootGain(x):
  return mul(sgn(x),sqrt(abs(x)))

def smoothingFilter(sigma,x):
  rgf = RecursiveGaussianFilter(sigma)
  y = like(x)
  rgf.apply0X(x,y)
  return y

def lowPassFilter(fpass,x):
  y = like(x)
  bf = ButterworthFilter(fpass,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(x,y)
  return y

def goFilter():
  #x,s1,s2,clip = getImageF3d(); png = "f3d_"
  x,s1,s2,clip = getImageTpd(); png = "tpd_"
  plot(x,s1,s2,clip,png=png+"input")
  doFilter(20,x,s1,s2,clip,png=png+"lbf")
  
def doFilter(sigma,x,s1,s2,clip,png=None):
  n1,n2 = s1.count,s2.count
  lbf = LocalPeakFrequencyFilter(sigma)
  f = like(x)
  for i2 in range(n2):
    f[i2] = lbf.findPeakFrequencies(x[i2])
    mul(1.0/s1.delta,f[i2],f[i2])
  plot(f,s1,s2,cbar="Frequency (Hz)",png=png+"f")

def goImage():
  x,s1,s2,clip = getImage()
  plot(x,s1,s2,clip)

def getImage():
  return getImageF3d()
  #return getImageTpd()

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

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def like(x):
  n1,n2 = len(x[0]),len(x)
  return zerofloat(n1,n2)

hbf = HilbertTransformFilter()
def ampPhase1(x):
  n = len(x)
  y = zerofloat(n)
  hbf.apply(n,x,y)
  return sqrt(add(mul(x,x),mul(y,y))),atan21(y,x)
def ampPhase2(x):
  n = len(x)
  a,p = like(x),like(x)
  for i in range(n):
    a[i],p[i] = ampPhase1(x[i])
  return a,p
def atan21(y,x):
  n = len(x)
  a = zerofloat(n)
  for i in range(n):
    a[i] = atan2(y[i],x[i])
  return a


#############################################################################
# plot

def plot(f,s1,s2,clip=0.0,t=None,cbar="",limits=None,png=None):
  n1,n2 = len(f[0]),len(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setFontSizeForPrint(8.0,240)
  sp.setSize(1020,730)
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
  sp.plotPanel.setColorBarWidthMinimum(120)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if clip!=0.0:
    pv.setClips(-clip,clip)
    pv.setColorModel(ColorMap.GRAY)
  else:
    pv.setColorModel(ColorMap.JET)
  if png and pngDir:
    sp.paintToPng(720,3.3,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
