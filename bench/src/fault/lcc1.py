#############################################################################
# Local cross-correlation for 1D sequences

import sys
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
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fault import Warp1
from fault.Util import *

#############################################################################

#pngDir = "./png"
pngDir = None

seed = 152 # c<0 & disc for many scales
seed = abs(Random().nextInt()/1000000)
seed = 1648 # disc for many scales
seed = 481 # OK for smaller scales, disc for larger scales
seed = 258 # shows slope constraints for constant shift = 16
print "seed =",seed

def main(args):
  goCorrelate()

def goCorrelate():
  n = 501
  fpeak = 0.125
  shift = 2.0/fpeak
  sigma = 1.0/fpeak
  mlag = 2*int(shift)
  nlag = 1+2*mlag
  nrms = 0.50
  w = Warp1.constant(shift,n)
  #w = Warp1.sinusoid(shift,n)
  s = getTrueShifts(n,w)
  f = makeRandomEvents(n,seed=seed); 
  #f = makeCosine(fpeak,n)
  g = w.warp(f)
  f = addRickerWavelet(fpeak,f)
  g = addRickerWavelet(fpeak,g)
  f = addNoise(nrms,fpeak,f,seed=10*seed+1)
  g = addNoise(nrms,fpeak,g,seed=10*seed+2)
  ff,gg = f,g
  #for scale in [1,2,4,8,16,32,64,128]:
  scale = 1.0
  ssigma = scale*sigma
  c = correlate(ssigma,-mlag,mlag,ff,gg)
  #c = cfake(n,-mlag,mlag)
  u = findShifts(ssigma,-mlag,mlag,ff,gg)
  v = findShiftsSmooth(0.0,16*sigma,-mlag,mlag,c)
  plot(ff,gg,c,s,u)
  plot(ff,gg,c,s,v)

def cfake(n1,lmin,lmax):
  nl = 1+lmax-lmin
  lpeak = 4.0
  c = rampfloat(lmin-lpeak,0.01,1.00,n1,nl)
  c = sub(1.0,mul(c,c))
  return c

def envelope(f):
  n1 = len(f)
  g = zerofloat(n1)
  htf = HilbertTransformFilter()
  htf.apply(n1,f,g)
  g = sqrt(add(mul(f,f),mul(g,g)))
  return g

def correlate(sigma,lmin,lmax,f,g):
  n1 = len(f)
  nl = 1+lmax-lmin
  lcf = makeLcf(sigma)
  lcf.setInputs(f,g)
  c = zerofloat(n1,nl)
  for il in range(nl):
    lag = il+lmin
    lcf.correlate(lag,c[il])
    lcf.normalize(lag,c[il])
  return c

def findShifts(sigma,lmin,lmax,f,g):
  n1 = len(f)
  lsf = makeLsf(sigma)
  u = zerofloat(n1)
  lsf.find1(lmin,lmax,f,g,u)
  return u

def getTrueShifts(n,w):
  d = zerofloat(n)
  for i in range(n):
    d[i] = w.ux(i)
  return d

def makeLcf(sigma):
  lcfType = LocalCorrelationFilter.Type.SIMPLE
  lcfWindow = LocalCorrelationFilter.Window.GAUSSIAN
  return LocalCorrelationFilter(lcfType,lcfWindow,sigma)

def makeLsf(sigma):
  return LocalShiftFinder(sigma)

def makeCosine(freq,n):
  return cos(mul(2.0*PI*freq,rampfloat(0.0,1.0,n)))

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

def plot(f,g,c,s,u):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,100)
  panel.setVLimits(1,slag.first,slag.last)
  fv = panel.addPoints(0,0,s1,f)
  gv = panel.addPoints(0,0,s1,g)
  gv.setLineColor(Color.RED)
  cv = panel.addPixels(1,0,s1,slag,c)
  #cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #cv.setClips(-1.0,1.0)
  cv.setColorModel(ColorMap.JET)
  sv = panel.addPoints(1,0,s)
  sv.setLineColor(Color.WHITE)
  sv.setLineStyle(PointsView.Line.DOT)
  sv.setLineWidth(3)
  uv = panel.addPoints(1,0,u)
  uv.setLineColor(Color.WHITE)
  uv.setLineWidth(3)
  panel.setHLabel("sample")
  panel.setVLabel(0,"f & g")
  panel.setVLabel(1,"lag")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
