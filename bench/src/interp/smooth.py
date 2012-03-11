#############################################################################
# Tests smoothing to match a specified inter-percentile range (IPR).
# Dave Hale, 2012.02.07, Colorado School of Mines
#
# I assume data x = signal s + measurement errors e, where the 
# errors e are uncorrelated with a known distribution (or at 
# least a known variance) and zero mean. I make no assumptions
# about statistics of the signal s. Among other things, this
# means that I do not know the signal-to-noise ratio.
#
# The data x are input to Laplacian smoothing with output
# y = smooth(sigmaS,x) ~ s, for some smoothing parameter 
# sigmaS. The big question: how to best determine sigmaS? 
#
# One answer is to smooth to match an IPR or some other 
# statistic of the noise distribution. This answer leads to 
# the following questions for which this program is relevant. 
# Some of these questions can also be addressed analytically.
#
# Why not instead use Wiener filtering?
# What smoothing parameter sigmaS is best for 
#   constant signal s?
#   a linear signal s?
# How should the best smoothing sigmaS vary with 
#   the number of samples in x?
#   frequency of the signal in s?
#   variance of the errors in e?
# Why match IPR instead of variance?
# How to measure which IPR is best?
# Which IPR is best for errors with
#   a known Gaussian distribution?
#   a known Laplace distribution?
#   an unknown (or incorrectly assumed) distribution?
#   outliers (occasional huge errors)
#############################################################################

from java.awt import *
from java.lang import *
from java.util import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.lapack import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.opt import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

def main(args):
  goTest()

def goTest():
  n = 1000 # number of samples in sequences x, s, e, ...
  ntrial = 100 # = 1 for one trial with plots of s, x, and y = smoothed(x)
  sigmaE = 0.1 # half-width of distribution for errors
  plist = [5,10,15,20,25,30,35] # list of lower percentile of IPRs to test
  flist = [0.0001,0.001,0.01] # list of frequencies to test
  for f in flist:
    print "frequency =",f,"cycles/sample"
    s = signal(n,f)
    sigmaSList = {}
    for p in plist:
      sigmaSList[p] = []
    for itrial in range(ntrial):
      e = noiseGaussian(n,sigmaE)
      #e = noiseLaplace(n,sigmaE)
      x = add(s,e)
      for p in plist:
        sigmaS = findSigmaS(p,x,sigmaE)
        sigmaSList[p].append(sigmaS)
        if ntrial==1:
          y = smooth(sigmaS,x)
          sp = SimplePlot()
          sp.setSize(900,500)
          sp.setTitle("p="+str(p)+"  sigmaS="+str(sigmaS))
          sp.addPoints(x)
          sp.addPoints(s).setLineColor(Color.BLUE)
          sp.addPoints(y).setLineColor(Color.RED)
    if ntrial>1: # statistics can be computed if more than one trial
      for p in plist:
        s = sigmaSList[p]
        sum1 = 0.0
        sum2 = 0.0
        for i in range(ntrial):
          sum1 += s[i]
          sum2 += s[i]*s[i]
        mu = sum1/ntrial
        sd = sqrt(sum2/ntrial-mu*mu) # valid only for ntrial >> 1
        print "  p =",p,"  mu =",mu,"  sd =",sd


# Returns y = smoothed(x) for specified smoothing half-width sigmaS
def smooth(sigmaS,x):
  y = zerofloat(len(x))
  ref = RecursiveExponentialFilter(sigmaS)
  ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE)
  ref.apply(x,y)
  return y

# Binary search for the sigmaS that matches the IPR for percentile p
def findSigmaS(p,x,sigmaE):
  n = len(x)
  func = Func(p,x,sigmaE)
  bzf = BrentZeroFinder(func)
  sigmasMax = 8.0
  while func.evaluate(sigmasMax)<0.0:
    sigmasMax *= 2.0
    if sigmasMax>n/2:
      return sigmasMax
  tolerance = 0.001*sigmasMax
  return bzf.findZero(0.0,sigmasMax,tolerance)

# Inverse CDF for Gaussian distribution
icdfGaussian = {
  5:1.6449, # 5 percent of values are greater than 1.6449
  10:1.2816,
  15:1.0364,
  20:0.8416,
  25:0.6745,
  30:0.5244,
  35:0.3853,
  40:0.2533,
  45:0.1257
}

# Inverse CDF for Laplace distribution
def icdfLap(p):
  if p<0.5:
    return  log(1.0+2.0*(p-0.5))
  else:
    return -log(1.0-2.0*(p-0.5))
icdfLaplace = {
  5:icdfLap(0.95),
  10:icdfLap(0.90),
  15:icdfLap(0.85),
  20:icdfLap(0.80),
  25:icdfLap(0.75),
  30:icdfLap(0.70),
  35:icdfLap(0.65),
  40:icdfLap(0.60),
  45:icdfLap(0.55),
}
 
# Function is < 0 if sigmaS is too small, > 0 if sigmaS too large
class Func(BrentZeroFinder.Function):
  def __init__(self,p,x,sigmaE):
    self.p = p
    self.x = x
    self.t = sigmaE*icdfGaussian[p]
    #self.t = sigmaE*icdfLaplace[p]
  def evaluate(self,sigmaS):
    p,x,t = self.p,self.x,self.t
    y = smooth(sigmaS,x)
    n = len(y)
    m = 0
    for i in range(n):
      if abs(x[i]-y[i])>t:
        m += 1
    return m-n*2.0*0.01*p

# Signal is a cosine with specified frequency.
def signal(n,freq):
  return cos(rampfloat(0.0,2.0*PI*freq,n))

# Noise distribution is Gaussian with specified std dev.
def noiseGaussian(n,sigma):
  r = Random()
  e = zerofloat(n)
  for i in range(n):
    e[i] = sigma*r.nextGaussian()
  return e

# Noise distribution is Laplace with specified std dev.
def noiseLaplace(n,sigma):
  r = Random()
  e = zerofloat(n)
  for i in range(n):
    u = r.nextDouble()-0.5
    if u<0.0:
      e[i] = -sigma*log(1.0+2.0*u)
    else:
      e[i] =  sigma*log(1.0-2.0*u)
  return e
  
#############################################################################
# Run the function main on the Swing thread
import sys
from javax.swing import *
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
if __name__ == "__main__":
  run(main)
