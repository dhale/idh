"""
Converts sgy files (SEG-Y format) to dat files (3D arrays of floats)
Removes all SEG-Y headers and, if necessary, converts IBM floats to 
IEEE floats. The byte order for floats in dat files is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Version: 2010.11.05
"""
from imports import *

nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header

#############################################################################
def main(args):
  goParihaka()

def goParihaka():
  """
INLINES  : 1665 - 5599 (INC 1)   CROSSLINES : 2050 - 14026 (INC 2)
TIME: 1501 samples  (6 sec, interval 4ms)
  n1 = 1501
  nbytes = 69,479,807,644
  ntrace = (nbytes-nhead-nbhed)/(240+4*n1)
  """
  n1 = 1501
  datdir = "/data/seis/nz/par/"
  sgyfile = datdir+"Parihaka3d_full.sgy"
  datfile = datdir+"par.dat"
  #readFormat(sgyfile)
  dumpTraceHeaders(sgyfile,n1)
  #testFormat(n1,n2,n3,sgyfile)
  #convert(n1,n2,n3,sgyfile,datfile)
  #display3d(n1,n2,n3,datfile,2.0e4)

def dumpTraceHeaders(sgyfile,n1):
  ais = ArrayInputStream(sgyfile)
  ais.skipBytes(nhead)
  ais.skipBytes(nbhed)
  h = zeroint(nthed/4)
  for i in range(2):
    ais.readInts(h)
    dump(h)
    print "i =",i," iline =",h[47]," xline =",h[48]
    ais.skipBytes(4*n1)
  ais.close()

def goParihakaSubsets():
  #n1,n2,n3 = 751,1501,701
  n1,n2,n3 = 751,1501,601
  datdir = "/data/seis/nz/par/"
  sgyfile = datdir+"parihaka_subset3.sgy"
  datfile = datdir+"par3.dat"
  check(n1,n2,n3,sgyfile)
  #readFormat(sgyfile)
  #testFormat(n1,n2,n3,sgyfile)
  #convert(n1,n2,n3,sgyfile,datfile)
  #display3d(n1,n2,n3,datfile,2.0e4)

def checkSubset(n1,n2,n3,sgyfile):
  print "checking",(n2*n3),"traces"
  nh = nthed/2
  nx = n1
  h = zeroshort(nh)
  x = zeroint(nx)
  ais = ArrayInputStream(sgyfile)
  ais.skipBytes(nhead+nbhed)
  for i3 in range(n3):
    for i2 in range(n2):
      ais.readShorts(h) # read trace header
      ais.readInts(x) # read trace samples
      ns = h[57] # number of samples in this trace
      if ns!=n1:
        print "*** i2 =",i2,"i3 =",i3,"ns =",ns
        print "*** trace header as shorts:"
        dump(h)
        return
    print "checked: i3 =",i3
  ais.close()

def convert(n1,n2,n3,sgyfile,datfile):
  print "converting",(n2*n3),"traces"
  ais = ArrayInputStream(sgyfile)
  aos = ArrayOutputStream(datfile)
  ais.skipBytes(nhead+nbhed)
  x = zeroint(n1)
  y = zerofloat(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      ais.skipBytes(nthed) # skip trace header
      ais.readInts(x) # read trace samples
      IbmIeee.ibmToFloat(x,y) # convert
      aos.writeFloats(y) # write trace samples
    print "converted: i3 =",i3
  ais.close()
  aos.close()

def display3d(n1,n2,n3,datfile,clip=0):
  ais = ArrayInputStream(datfile)
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  print "x min =",min(x)," max =",max(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,x)
  if clip>0.0:
    ipg.setClips(-clip,clip)
  frame.setVisible(True)

def readFormat(sgyfile):
  ais = ArrayInputStream(sgyfile)
  ais.skipBytes(nhead)
# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
  h = zeroshort(nbhed/2)
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

def testFormat(n1,n2,n3,sgyfile):
  xi = zeroint(n1)
  x1 = zerofloat(n1)
  x2 = zerofloat(n1)
  ais = ArrayInputStream(sgyfile)
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
