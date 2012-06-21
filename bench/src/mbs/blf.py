"""
Applies bilateral filter to mbs data.

Author: Dave Hale, Colorado School of Mines
Version: 2012.06.18
"""
from imports import *

"""
Subset of PstmLarge image:
i1min,i1max = 150, 650,  n1 = 501
i2min,i2max = 490,1258,  n2 = 769
i3min,i3max = 358, 917,  n3 = 560
"""
global n1,n2,n3
n1,n2,n3 = 501,769,560
d1,d2,d3 = 0.002,0.016764,0.016764 # 2 ms, 55 ft, 55 ft
f1,f2,f3 = 0.300,0.000000,0.000000
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
datdir = "/data/seis/mbs/dat/"

def main(args):
  goScale()
  
def goScale():
  datfile = datdir+"pstm_raw_s1.dat"
  f = readImage(datfile)
  mul(0.0001,f,f)
  show3d(f,clip=1.0)
  
def show3d(f,clip=None):
  print "show3d: f min =",min(f)," max =",max(f)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

def readImage(datfile):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datfile)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datfile,x):
  aos = ArrayOutputStream(datfile)
  aos.writeFloats(x)
  aos.close()

def structureTensors(sigma1,sigma23,f):
  lof = LocalOrientFilter(sigma1,sigma23)
  #lof.setGradientSmoothing(sigma)
  t = lof.applyForTensors(f)
  return t

def diffusionTensors(sigma,x):
  t = structureTensors(sigma,x) # structure tensors
  t.invertStructure(0.0,2.0,4.0) # invert with ew = 1, ev small, eu smaller
  return t

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
