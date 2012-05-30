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
  goSimple()

def goSimple():
  st,x,y = makeXY()
  SimplePlot.asPoints(st,x)
  SimplePlot.asPoints(st,y)
  fpeak = getPeakFrequency(x)
  print "peak frequency =",fpeak/st.delta
  nt = st.count
  fft = Fft(nt)
  fft.setCenter(True)
  sf = fft.getFrequencySampling1()
  ax = cabs(fft.applyForward(x))
  ay = cabs(fft.applyForward(y))
  SimplePlot.asPoints(sf,ax)
  SimplePlot.asPoints(sf,ay)
  shift = 0.05
  gx = makeCosineModulator(shift,x)
  gy = makeCosineModulator(shift,y)
  SimplePlot.asPoints(st,gx)
  SimplePlot.asPoints(st,gy)
  xgx = mul(x,gx)
  ygy = mul(y,gy)
  SimplePlot.asPoints(st,xgx)
  SimplePlot.asPoints(st,ygy)
  axgx = cabs(fft.applyForward(xgx))
  aygy = cabs(fft.applyForward(ygy))
  SimplePlot.asPoints(sf,axgx)
  SimplePlot.asPoints(sf,aygy)
  xgx = lowPassFilter(fpeak,xgx)
  ygy = lowPassFilter(fpeak,ygy)
  #xgx = smoothingFilter(25.0,xgx)
  #ygy = smoothingFilter(25.0,ygy)
  SimplePlot.asPoints(st,xgx)
  SimplePlot.asPoints(st,ygy)

def addRickerWavelet(fpeak,x):
  n = x.length
  y = like(x)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih);
  Conv.conv(nh,-ih,h,n,0,x,n,0,y)
  return y

def ricker(fpeak,t):
  x = PI*fpeak*time
  xx = x*x
  return (1.0-2.0*xx)*exp(-xx)

def makeXY():
  nt,dt,ft = 251,0.004,0.0
  st = Sampling(nt,dt,ft)
  fpeak = 0.1/dt
  wpeak = 2.0*PI*fpeak
  sigma = 5*dt
  sigmas = sigma*sigma
  s = 10*dt
  t = rampfloat(ft,dt,nt)
  p = mul(wpeak,t) # phase
  t = sub(t,5.0*s)
  x = mul(cos(p),exp(neg(mul(0.5/sigmas,mul(t,t)))))
  t = sub(t,5.0*s)
  y = mul(cos(p),exp(neg(mul(0.5/sigmas,mul(t,t)))))
  t = sub(t,5.0*s)
  z = mul(cos(p),exp(neg(mul(0.5/sigmas,mul(t,t)))))
  x = add(x,y)
  y = add(y,z)
  return st,x,y

def envelope(x):
  n = len(x)
  y = like(x)
  hbf = HilbertTransformFilter()
  hbf.apply(n,x,y)
  return sqrt(add(mul(x,x),mul(y,y)))

def getPeakFrequency(x):
  n = len(x)
  f = getPeakFrequencies(n/4.0,x)
  return sum(f)/n

def getPeakFrequencies(sigma,x):
  n = len(x)
  lbf = LocalBurgFilter()
  c = lbf.getCoefficients(sigma,2,x)
  f = lbf.getPeakFrequencies(c)
  return f

def makeCosineModulator(shift,x):
  n = len(x)
  fpeak = getPeakFrequency(x)
  fpeak += shift
  phase = rampfloat(0.0,2.0*PI*fpeak,n)
  return cos(phase)

def makeBandModulator(shift,x):
  fpeak = getPeakFrequency(x)
  fpeak += shift
  y = like(x)
  nf = NotchFilter(fpeak,0.9)
  nf.applyForwardReverse(x,y)
  return sub(x,y)

def makeBurgModulator(shift,x):
  n = len(x)
  y = like(x)
  lbf = LocalBurgFilter()
  c = lbf.getCoefficients(16,2,x)
  c = lbf.zeroPeakFrequencies(c)
  c = lbf.shiftPeakFrequencies(shift,c)
  t = lbf.applyNotch(0.9,c,x)
  t = reverse(t)
  t = lbf.applyNotch(0.9,c,t)
  y = reverse(t)
  return sub(x,y)

def squareRootGain(x):
  return mul(sgn(x),sqrt(abs(x)))

def smoothingFilter(sigma,x):
  y = like(x)
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply0(x,y)
  return y

def lowPassFilter(fpass,x):
  y = like(x)
  bf = ButterworthFilter(fpass,6,ButterworthFilter.Type.LOW_PASS)
  bf.applyForwardReverse(x,y)
  return y

def like(x):
  return zerofloat(len(x))

#############################################################################
# plot

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
