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

n1 = 1501  
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

#############################################################################
# Research

class WhitenFilter(FileFloat3Chunks.Filter):
  def apply(self,i1,i2,i3,x,y):
    x,y = x[0],y[0]
    n1 = len(x[0][0])
    n2 = len(x[0])
    n3 = len(x)
    print "apply"
    print "  i1 =",i1," i2 =",i2," i3 =",i3
    print "  n1 =",n1," n2 =",n2," n3 =",n3
    print "  min =",Array.min(x)," max =",Array.max(x)

def whiten():
  sigma1,sigma2,sigma3 = 12.0,6.0,6.0
  l1,l2,l3 = int(3*sigma1),int(3*sigma2),int(3*sigma3)
  wf = WhitenFilter()
  #mc = 50000000 # 50 Mfloats
  mc = 100000000 # 100 Mfloats
  ff3c = FileFloat3Chunks(mc,n1,l1,l1,n2,l2,l2,n3,l3,l3)
  print "whiten: chunk size =",ff3c.getChunkSize()
  afx = [ArrayFile(datadir+"s02.dat","r")]
  afy = [ArrayFile(datadir+"w02.dat","rw")]
  ff3c.apply(wf,afx,afy)

def main(args):
  whiten()
  return

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
