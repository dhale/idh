import sys
from math import *

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

dataDir = "/data/mri/"
def samplingHead():
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  n1=256; d1=1.000; f1=0.000; s1 = Sampling(n1,d1,f1)
  n2=256; d2=1.000; f2=0.000; s2 = Sampling(n2,d2,f2)
  n3=109; d3=1.000; f3=0.000; s3 = Sampling(n3,d3,f3)
def samplingKnee():
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  n1=256; d1=1.000; f1=0.000; s1 = Sampling(n1,d1,f1)
  n2=256; d2=1.000; f2=0.000; s2 = Sampling(n2,d2,f2)
  n3=127; d3=1.000; f3=0.000; s3 = Sampling(n3,d3,f3)
def samplingCtHead():
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  n1=256; d1=1.000; f1=0.000; s1 = Sampling(n1,d1,f1)
  n2=256; d2=1.000; f2=0.000; s2 = Sampling(n2,d2,f2)
  n3=113; d3=2.000; f3=0.000; s3 = Sampling(n3,d3,f3)
def samplingBrain():
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  n1=256; d1=1.000; f1=0.000; s1 = Sampling(n1,d1,f1)
  n2=256; d2=1.000; f2=0.000; s2 = Sampling(n2,d2,f2)
  n3=109; d3=1.000; f3=0.000; s3 = Sampling(n3,d3,f3)

def main(args):
  #convertImages()
  #x = readImage(samplingHead,"head.dat")
  x = readImage(samplingBrain,"brain.dat")
  #x = readImage(samplingKnee,"knee.dat")
  #x = readImage(samplingCtHead,"cthead.dat")
  plot3d(x)
  return

def readImage(sampling,ffile):
  sampling()
  ais = ArrayInputStream(dataDir+ffile)
  f = zerofloat(n1,n2,n3)
  ais.readFloats(f)
  ais.close()
  return f

def plot3d(x):
  print "x min =",min(x)," max =",max(x)
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  s3 = Sampling(n3,d3,f3)
  sf = SimpleFrame()
  ip = sf.addImagePanels(s1,s2,s3,x)

def convertImages():
  samplingHead()
  convertImage("head.du2","head.dat")
  samplingCtHead()
  convertImage("cthead.du2","cthead.dat")
  samplingKnee()
  convertImage("knee.du2","knee.dat")
  samplingBrain()
  convertImage("brain.du2","brain.dat")

def convertImage(sfile,ffile):
  bo = ByteOrder.LITTLE_ENDIAN
  ais = ArrayInputStream(dataDir+sfile,bo)
  s = zeroshort(n1,n2,n3)
  ais.readShorts(s)
  ais.close()
  f = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        f[i3][i2][i1] = s[i3][i2][i1]
  aos = ArrayOutputStream(dataDir+ffile)
  aos.writeFloats(f)
  aos.close()
  print "wrote",4*n1*n2*n3,"bytes to",ffile

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
