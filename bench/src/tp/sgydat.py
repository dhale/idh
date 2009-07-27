"""
Converts sgy files (SEG-Y format) to dat files (3D arrays of floats)
For Teapot Dome data, this means converting IBM floats to IEEE floats 
and removing all SEG-Y headers.

Note that in Transform's depth images the first trace is missing from
each line. Specifically, in Transform's depth images, I observed that
first trace with indices (  0,  0) has (x,y) = (789048,938850)
 last trace with indices (186,344) has (x,y) = (808604,977165)
The (x,y) coordinates tell us that it is the first (not the last)
trace in each line that is missing. For each line, this conversion 
script recreates the missing trace by duplicating the first trace read.
In this way, this script writes 188 traces per line to all dat files.

Author: Dave Hale, Colorado School of Mines
Version: 2009.06.07
"""
from imports import *

#############################################################################
def main(args):
  #setGlobals("tz"); readFormat()
  #setGlobals("st"); testFormat()
  #process("st")
  #process("sz")
  #process("tz")

def process(what):
  setGlobals(what)
  convert()
  display()

nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header
sgyFile,datFile = "",""
n1,n2,n3 = 0,0,0
missingTrace = False
def setGlobals(what):
  global sgyFile,datFile,n1,n2,n3,missingTrace
  tpDir = "/data/seis/tp/"
  if what=="st": # seismic time image
    doeDir = tpDir+"doe/3D_Seismic/"
    sgyFile = doeDir+"filt_mig.sgy"
    datFile = doeDir+"tpstAll.dat"
    n1,n2,n3 = 1501,188,345
    missingTrace = False
  elif what=="sz": # seismic depth image
    tssDir = tpDir+"tss/"
    sgyFile = tssDir+"tpszAll.sgy"
    datFile = tssDir+"tpszAll.dat"
    n1,n2,n3 = 2762,188,345
    missingTrace = True
  elif what=="tz": # time(depth) image
    tssDir = tpDir+"tss/"
    sgyFile = tssDir+"tptzAll.sgy"
    datFile = tssDir+"tptzAll.dat"
    n1,n2,n3 = 2762,188,345
    missingTrace = True
  print "sgyFile =",sgyFile
  print "datFile =",datFile
  print "n1 =",n1," n2 =",n2," n3 =",n3

def convert():
  print "converting",(n2*n3),"traces"
  ais = ArrayInputStream(sgyFile)
  aos = ArrayOutputStream(datFile)
  ais.skipBytes(nhead+nbhed)
  x = zeroint(n1)
  y = zerofloat(n1)
  n = 0
  if missingTrace:
    j2 = 1
  else:
    j2 = 0
  for i3 in range(n3):
    for i2 in range(j2,n2):
      ais.skipBytes(nthed)
      ais.readInts(x)
      IbmIeee.ibmToFloat(x,y)
      if i2==1 and missingTrace:
        aos.writeFloats(y) # duplicate first trace, if missing
        n += 1
      aos.writeFloats(y)
      n += 1
      if n%1000==0:
        print "converted",n,"traces"
  ais.close()
  aos.close()

def display():
  ais = ArrayInputStream(datFile)
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  print "x min =",min(x)," max =",max(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def readFormat():
  ais = ArrayInputStream(sgyFile)
  ais.skipBytes(nhead)
# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
  h = zeroshort(nbhed)
  ais.readShorts(h)
  ais.close()
  print "dump of binary header as shorts"
  dump(h)
  print "current sampling interval in usec =",h[8]
  print "original sampling interval in usec =",h[9]
  print "number of samples per trace =",h[10]
  print "original number of samples per trace =",h[11]
  format = h[12]
  if format==1:
    print "format = 1 = IBM floating point"
  elif format==5:
    print "format = 5 = IEEE floating point"
  else:
    print "format =",format,"is unknown!"

def testFormat():
  xi = zeroint(n1)
  x1 = zerofloat(n1)
  x2 = zerofloat(n1)
  ais = ArrayInputStream(sgyFile)
  ais.skipBytes(nhead+nbhed)
  ais.skipBytes(n3/2*n2*(nthed+4*n1))
  ais.skipBytes(n2/2*(nthed+4*n1))
  ais.skipBytes(nthed)
  ais.readInts(xi)
  ais.close()
  IbmIeee.ibmToFloat(xi,x1)
  IbmIeee.ieeeToFloat(xi,x2)
  sp = SimplePlot.asPoints(x1); sp.setTitle("Assuming IBM format")
  sp = SimplePlot.asPoints(x2); sp.setTitle("Assuming IEEE format")

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
