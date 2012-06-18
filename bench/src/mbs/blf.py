"""
Applies bilateral filter to mbs data.

Author: Dave Hale, Colorado School of Mines
Version: 2012.06.18
"""
from imports import *

"""
Geometry (corner points)
i2a

iline,xline =1422

"""

"""
Subset of PstmLarge image
i1min,i1max = 150, 650 # n1 = 501 (time = 0.3 - 1.3 s)
i2min,i2max = 358, 917 # n2 = 560 (iline = 358 -  917)
i3min,i3max = 490,1258 # n3 = 763 (xline = 490 - 1258)
"""
global n1,n2,n3
n1,n2,n3 = 501,560,763
d1,d2,d3 = 0.002,0.016764,0.016764 # 2 ms, 55 ft, 55 ft
f1,f2,f3 = 0.300,0.000000,0.000000
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
datDir = "/data/seis/mbs/dat/"

def main(args):
  goScale()
  
def goScale():
  datFile = datDir+"pstm_raw_s1.dat"
  x = readImage(datFile)
  #mul(10000.0,x,x)
  #writeImage(datFile,x)
  show3(x,clip=1.0)

def show3(x,clip=0.0):
  print "show3: min =",min(x)," max =",max(x)
  frame = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  frame.setSize(1200,1000)
  #view.setAxesScale(1.0,1.0,4.0)
  ipg = frame.addImagePanels(s1,s2,s3,x)
  if clip>0.0:
    ipg.setClips(-clip,clip)
  view = frame.orbitView
  view.setAzimuthAndElevation(40.0,30.0)
  view.setScale(3.0)
  view.setAxesScale(1.0,1.0,8.0)

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

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
