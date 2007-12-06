import sys
from math import *
from java.awt import *
from java.lang import *
from javax.swing import *
from lcc import *
from sw import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

True = 1
False = 0

#n1 = 1501  
n1 = 301  
d1 = 0.004
f1 = 0.0
s1 = Sampling(n1,d1,f1)

n2 = 623  
d2 = 0.025
f2 = 0
s2 = Sampling(n2,d2,f2)

n3 = 367
d3 = 0.025
f3 = 0
s3 = Sampling(n3,d3,f3)

datadir = "/data/seis/sw/all/"
datbdir = "/datb/seis/sw/all/"
datcdir = "/datb/seis/sw/all/"

##############################################################################
# Read/write

def readFloats3(file):
  f = Array.zerofloat(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readFloats(f)
  af.close()
  return f

def writeFloats3(file,f):
  af = ArrayFile(datadir+file,"rw")
  af.writeFloats(f)
  af.close()
  return f

def readBytes3(file):
  b = Array.zerobyte(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readBytes(b)
  af.close()
  return b

def writeBytes3(file,lag):
  af = ArrayFile(datadir+file,"rw")
  af.writeBytes(lag)
  af.close()

def readFloats12(file,i3):
  f = readFloats3(file)
  return slice12(f,i3)

def readFloats13(file,i2):
  f = readFloats3(file)
  return slice13(f,i2)

def readFloats23(file,i1):
  f = readFloats3(file)
  return slice23(f,i1)

def slice12(f,i3):
  f12 = Array.copy(f[i3])
  return f12

def slice13(f,i2):
  f13 = Array.zerofloat(n1,n3)
  for i3 in range(n3):
    Array.copy(n1,f[i3][i2],f13[i3])
  return f13

def slice23(f,i1):
  f23 = Array.zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      f23[i3][i2] = f[i3][i2][i1];
  return f23

def copyFile(x,y):
  print "copyFile:",x,"to",y
  x = ArrayFile(x,"r")
  y = ArrayFile(y,"rw")
  n = 65536
  b = Array.zerobyte(n)
  size = x.length()
  while size>0:
    m = min(n,size)
    x.readBytes(b,0,m)
    y.writeBytes(b,0,m)
    size -= m
  x.close()
  y.close()

def zeroFile(n1,n2,n3,x):
  print "zeroFile:",x
  x = ArrayFile(x,"rw")
  n = 65536
  b = Array.zerobyte(n)
  size = 4*n1*n2*n3
  while size>0:
    m = min(n,size)
    x.writeBytes(b,0,m)
    size -= m
  x.close()

#############################################################################
# Research

class WhitenFilter(FileFloat3Chunks.Filter):
  def __init__(self,sigma1,sigma2,sigma3):
    self.sf = ShiftFinder(sigma1,sigma2,sigma3)
  def apply(self,i1,i2,i3,x,y):
    x,y = x[0],y[0]
    n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
    xmin,xmax = Array.min(x),Array.max(x)
    tiny = 1.0e-5*xmax
    print "WhitenFilter.apply"
    print "  min =",xmin," max =",xmax
    print "  i1 =",i1," i2 =",i2," i3 =",i3
    print "  n1 =",n1," n2 =",n2," n3 =",n3
    r = Array.randfloat(n1,n2,n3)
    Array.mul(tiny,r,r)
    Array.add(r,x,x)
    r = None
    self.sf.whiten(x,y)

def whitenArrayFile(afx,afy):
  sigma1,sigma2,sigma3 = 12,6,6
  wf = WhitenFilter(sigma1,sigma2,sigma3)
  mc = 30000000 # 30 Mfloats
  #mc = 120000000 # 120 Mfloats
  l1,l2,l3 = 3*sigma1,3*sigma2,3*sigma3
  ff3c = FileFloat3Chunks(mc,n1,l1,l1,n2,l2,l2,n3,l3,l3)
  print "whitenArrayFile: chunk size =",ff3c.getChunkSize()
  ff3c.apply(wf,[afx],[afy])

def whiten():
  fs = ["s02","s04"]
  fw = ["w02","w04"]
  #fs = ["s04"]
  #fw = ["w04"]
  for i in range(len(fs)):
    sname = datadir+fs[i]+".dat"
    wname = datadir+fw[i]+".dat"
    print "whitening",sname,"..."
    afs = ArrayFile(sname,"r")
    afw = ArrayFile(wname,"rw")
    whitenArrayFile(afs,afw)
    afs.close()
    afw.close()
    print "... done"

class ShiftFinderFilter(FileFloat3Chunks.Filter):
  def __init__(self,lmin,lmax,sigma1,sigma2,sigma3):
    self.sf = ShiftFinder(sigma1,sigma2,sigma3)
    self.lmin = lmin
    self.lmax = lmax
  def apply(self,i1,i2,i3,x,y):
    n1,n2,n3 = len(x[0][0][0]),len(x[0][0]),len(x[0])
    print "ShiftFinderFilter.apply:"
    print "  i1 =",i1," i2 =",i2," i3 =",i3
    print "  n1 =",n1," n2 =",n2," n3 =",n3
    f,g,u1,u2,u3 = x[0],x[1],x[2],x[3],x[4] # input and output arrays
    print "  f min =",Array.min(f)," max =",Array.max(f)
    print "  g min =",Array.min(g)," max =",Array.max(g)
    n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
    du = Array.zerofloat(n1,n2,n3)
    print "  shift1"
    self.sf.find1(self.lmin,self.lmax,f,g,du)
    print "    du min =",Array.min(du)," max =",Array.max(du)
    self.sf.shift1(du,u1,u2,u3,g)
    print "    complete"
    print "  shift3"
    self.sf.find3(self.lmin,self.lmax,f,g,du)
    print "    du min =",Array.min(du)," max =",Array.max(du)
    self.sf.shift3(du,u1,u2,u3,g)
    print "    complete"
    print "  shift2"
    self.sf.find2(self.lmin,self.lmax,f,g,du)
    print "    du min =",Array.min(du)," max =",Array.max(du)
    self.sf.shift2(du,u1,u2,u3,g)
    print "    complete"
 
def findShifts():
  nshift = 1
  lmin,lmax = -2,2
  sigma1,sigma2,sigma3 = 12,12,12
  mc = 30000000 # 30 Mfloats
  l1,l2,l3 = 3*sigma1,3*sigma2,3*sigma3
  ff3c = FileFloat3Chunks(mc,n1,l1,l1,n2,l2,l2,n3,l3,l3)
  sff = ShiftFinderFilter(lmin,lmax,sigma1,sigma2,sigma3)
  print "findShifts: chunk size =",ff3c.getChunkSize()
  fname = datadir+"w02.dat"
  gname = datadir+"w04.dat"
  zname = datadir+"zeros.dat"
  fa = ArrayFile(fname,"r")
  ga = ArrayFile(gname,"r")
  zeroFile(n1,n2,n3,zname)
  u1a = ArrayFile(zname,"r")
  u2a = ArrayFile(zname,"r")
  u3a = ArrayFile(zname,"r")
  for ishift in range(1,nshift+1):
    print "findShifts: iteration =",ishift
    gb  = ArrayFile(datadir+"gs"+str(ishift)+".dat","rw")
    u1b = ArrayFile(datadir+"u1s"+str(ishift)+".dat","rw")
    u2b = ArrayFile(datadir+"u2s"+str(ishift)+".dat","rw")
    u3b = ArrayFile(datadir+"u3s"+str(ishift)+".dat","rw")
    af = [fa,ga,u1a,u2a,u3a]
    bf = [   gb,u1b,u2b,u3b]
    ip = [    1,  2,  3,  4]
    ff3c.apply(sff,af,bf,ip)
    ga.close()
    u1a.close()
    u2a.close()
    u3a.close()
    ga,u1a,u2a,u3a = gb,u1b,u2b,u3b
  fa.close()
  ga.close()
  u1a.close()
  u2a.close()
  u3a.close()

def main(args):
  #whiten()
  findShifts()
  return

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
