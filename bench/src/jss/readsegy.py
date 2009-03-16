import sys

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

import jss.Convert as Convert

dataDir = "/data/seis/jss/"

# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header

# input segy sampling
"""
n1i=1000 # number of time samples (1st dimension)
d1i=0.006 # 0.006 km = 6 m
f1i=0.000 
n2i=1600 # number of traces in 2nd dimension
d2i=0.0125 # 0.0125 km = 12.5 m
f2i=0.000
k1=150 # index of first sample in windowed data
"""
n1i=1600 # number of time samples (1st dimension)
d1i=0.004 # 0.004 km = 4 m
f1i=0.000 
n2i=2000 # number of traces in 2nd dimension
d2i=0.010 # 0.010 km = 10 m
f2i=0.000

# output data sampling (omit water layer/bottom and sides)
k1=200 # index of first sample in windowed data
n1=n1i-k1
d1=d1i
f1=k1*d1
k2=200 # index of first trace in windowed data
n2=n2i-2*k2
d2=d2i
f2=k2*d2
print "n1 =",n1,"d1 =",d1," f1=",f1
print "n2 =",n2,"d2 =",d2," f2=",f2

sgyFiles = [
  "Ref_Conv.segy", "Ref_SimSrc.segy",
  "Mon_Conv.segy", "Mon_SimSrc.segy"] 

datFiles = [
  "s1f.dat", "s2f.dat",
  "s1g.dat", "s2g.dat"]

def main(args):
  #testFormat()
  #readFormat()
  convert()

def convert():
  nfiles = len(sgyFiles)
  for i in range(nfiles):
    print datFiles[i]
    convertSegy(sgyFiles[i],datFiles[i])
    plot(readImage(datFiles[i]))
  return

def plot(x):
  xmin,xmax = Array.min(x),Array.max(x)
  print "x min =",xmin," max =",xmax
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(756,791)
  sp.setFontSize(36)
  sp.setHLabel("distance (km)")
  sp.setVLabel("depth (km)")
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  pv = sp.addPixels(s1,s2,x)
  pv.setPercentiles(2.0,98.0)
  #pv.setClips(-10,10)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  print "clips: min =",pv.clipMin," max =",pv.clipMax

def readFormat():
  infile = dataDir+"Ref_Gold.segy"
  ais = ArrayInputStream(infile)
  ais.skipBytes(nhead)
# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
  h = Array.zeroshort(nbhed)
  ais.readShorts(h)
  print "format =",h[12]
  ais.close()

def testFormat():
  infile = dataDir+"Ref_Gold.segy"
  ais = ArrayInputStream(infile)
  ais.skipBytes(nhead+nbhed)
  ais.skipBytes(nthed)
  xi = Array.zeroint(n1)
  ais.readInts(xi)
  ais.close()
  x1 = Array.zerofloat(n1)
  x2 = Array.zerofloat(n1)
  Convert.ibmToFloat(xi,x1)
  Convert.ieeeToFloat(xi,x2)
  SimplePlot.asPoints(x1)
  SimplePlot.asPoints(x2)
  #Array.dump(x1)
  #Array.dump(x2)

def convertSegy(infile,outfile):
  infile = dataDir+infile;
  outfile = dataDir+outfile
  ais = ArrayInputStream(infile)
  aos = ArrayOutputStream(outfile)
  ais.skipBytes(nhead+nbhed)
  xi = Array.zeroint(n1i)
  x = Array.zerofloat(n1i)
  y = Array.zerofloat(n1)
  ais.skipBytes(k2*(nthed+4*n1i))
  for i2 in range(n2):
    ais.skipBytes(nthed)
    ais.readInts(xi)
    Convert.ibmToFloat(xi,x)
    Array.copy(n1,k1,x,0,y)
    aos.writeFloats(y)
  ais.close()
  aos.close()

def readImage(file):
  ais = ArrayInputStream(dataDir+file)
  x = Array.zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
