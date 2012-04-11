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
global n1,n2,n3

#############################################################################
def main(args):
  #goParihaka()
  goF3d()

def goF3d():
  """
  TIME: 462 samples (1.848 sec, interval 4 ms)
  nbytes = 699,003,060
  ntrace = (nbytes-nhead-nbhed)/(240+2*n1)
  """
  fmt = 3 # 2-byte shorts in [-32767,32767]
  global n1,n2,n3
  n1,n2,n3 = 462,951,591
  nbytes = 699003060
  datdir = "/data/seis/f3d/"
  sgyfile = datdir+"f3draw.sgy"
  mapfile = datdir+"f3dmap.dat"
  datfile = datdir+"f3draw.dat"
  displayF3d(datfile)
  #displayF3dTrace(datfile)
  #makeSubsetsF3d(datdir)
  #displaySubsetF3d(datfile,0,500,500,751,501,501)
  #bigSubsetF3d(n1,sgyfile,datfile)
  #readFormat(sgyfile)
  #makeMap(sgyfile,mapfile,nbytes,fmt,n1)
  #dumpTraceHeaders(sgyfile,fmt,n1)
  #convert3(n1,n2,n3,sgyfile,datfile)

def displayF3d(datfile):
  x = readImage(datfile,n1,n2,n3)
  xmax = 30767.0
  display3d(sub(x,clip(-xmax,xmax,x)),1.0)
  #mul(0.001,x,x)
  #display3d(x,10.0)
  #display3d(slog(x),3.0)


def displayF3dTrace(datfile):
  x = readImage(datfile,n1,n2,n3)
  mul(0.001,x,x)
  x = x[170][252]
  SimplePlot.asPoints(x)
  SimplePlot.asPoints(spow(0.5,x))
  SimplePlot.asPoints(slog(x))

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

def goParihaka():
  """
  INLINES  : 1665 - 5599 (INC 1)   CROSSLINES : 2050 - 14026 (INC 2)
  TIME: 1501 samples  (6 sec, interval 4ms)
  nbytes = 69,479,807,644
  ntrace = (nbytes-nhead-nbhed)/(240+4*n1)
  """
  fmt = 1
  n1 = 1501
  nbytes = 69479807644
  datdir = "/data/seis/nz/par/"
  sgyfile = datdir+"Parihaka3d_raw.sgy"
  mapfile = datdir+"map.dat"
  datfile = datdir+"par11.dat"
  #displayParihaka(datfile)
  #makeSubsetsParihaka(datdir)
  #displaySubsetParihaka(datfile,0,500,500,751,501,501)
  #bigSubsetParihaka(n1,sgyfile,datfile)
  readFormat(sgyfile)
  #makeMap(sgyfile,mapfile,nbytes,n1)
  #dumpTraceHeaders(sgyfile,fmt,n1)
  #testFormat(n1,n2,n3,sgyfile)
  #convert(n1,n2,n3,sgyfile,datfile)

def displayParihaka(datfile):
  n1,n2,n3 = 751,1001,1001
  x = readImage(datfile,n1,n2,n3)
  display3d(x,1.0e5)

def makeSubsetsParihaka(datdir):
  m1,m2,m3 = 751,1001,1001
  x = readParihaka(datdir+"parBig.dat",0,0,0,m1,m2,m3)
  writeImage(datdir+"par00.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,1000,0,m1,m2,m3)
  writeImage(datdir+"par01.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,0,1000,m1,m2,m3)
  writeImage(datdir+"par10.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,1000,1000,m1,m2,m3)
  writeImage(datdir+"par11.dat",x)
  display3d(x,1.0e5)

def displaySubsetParihaka(datfile,j1,j2,j3,m1,m2,m3):
  x = readSubsetParihaka(datfile,j1,j2,j3,m1,m2,m3)
  display3d(x,1.0e5)

def readImage(datfile,n1,n2,n3):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datfile)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datfile,x):
  aos = ArrayOutputStream(datfile)
  aos.writeFloats(x)
  aos.close()

def readParihaka(datfile,j1,j2,j3,m1,m2,m3):
  n1,n2,n3 = 1501,2001,2001
  m1 = min(m1,n1)
  m2 = min(m2,n2)
  m3 = min(m3,n3)
  j1 = min(j1,n1-m1)
  j2 = min(j2,n2-m2)
  j3 = min(j3,n3-m3)
  x = zerofloat(m1,m2,m3)
  ais = ArrayInputStream(datfile)
  for i3 in range(j3):
    ais.skipBytes(4*n1*n2)
  ais.skipBytes(4*(j1+n1*j2))
  for i3 in range(m3):
    if i3%10==0: 
      print "i3 =",i3
    for i2 in range(m2):
      ais.readFloats(x[i3][i2])
      ais.skipBytes(4*(n1-m1))
    ais.skipBytes(4*n1*(n2-m2))
  return x

def bigSubsetParihaka(n1,sgyfile,datfile):
  """ big subset 1501x2001x2001 ~ 24 GB
  i1min,i1max =    0, 1500 # time samples
  i2min,i2max = 6500,10500 # increment by 2
  i3min,i3max = 2100, 4100 # increment by 1
  """
  """ A: small subset 501x501x501 ~ 500 MB
  i1min,i1max =    0, 500 # time samples
  i2min,i2max = 7500,8500 # increment by 2
  i3min,i3max = 2500,3000 # increment by 1
  """
  """ B: subset 751x1001x1001 ~ 500 MB
  i1min,i1max =    0, 750 # time samples
  i2min,i2max = 7000,9000 # increment by 2
  i3min,i3max = 3000,4000 # increment by 1
  """
  i1min,i1max =    0, 1500 # time samples
  i2min,i2max = 6500,10500 # increment by 2
  i3min,i3max = 2100, 4100 # increment by 1
  m1 = 1+i1max-i1min
  m2 = 1+(i2max-i2min)/2
  m3 = 1+i3max-i3min
  m23 = m2*m3
  ais = ArrayInputStream(sgyfile)
  aos = ArrayOutputStream(datfile)
  ais.skipBytes(nhead)
  ais.skipBytes(nbhed)
  h = zeroint(nthed/4)
  x = zeroint(n1)
  y = zeroint(m1)
  z = zerofloat(m1)
  nread = 0
  i23 = 0
  while i23<m23:
    nread += 1
    ais.readInts(h)
    i3,i2 = h[49],h[50]
    if i2min<=i2 and i2<=i2max and i3min<=i3 and i3<=i3max:
      i23 += 1
      if i2==i2min:
        print "nread =",nread," i3min =",i3min," i3 =",i3," i3max =",i3max
      ais.readInts(x) # read trace samples
      copy(m1,x,y) # keep only m1 samples
      IbmIeee.ibmToFloat(y,z) # ibm to ieee
      aos.writeFloats(z) # write trace samples
    else:
      ais.skipBytes(4*n1)
  ais.close()
  aos.close()

def makeMap(sgyfile,mapfile,nbytes,fmt,n1):
  if fmt==3:
    bps = 2
  else:
    bps = 4
  ntrace = (nbytes-nhead-nbhed)/(240+bps*n1)
  af = ArrayFile(sgyfile,"r")
  af.skipBytes(nhead)
  af.skipBytes(nbhed)
  h = zeroint(nthed/4)
  #m2,m3 = 8000,16000 # Parihaka
  m2,m3 = 800,1300 # F3D
  m = zerofloat(m2,m3)
  i2min =  Integer.MAX_VALUE
  i2max = -Integer.MAX_VALUE
  i3min =  Integer.MAX_VALUE
  i3max = -Integer.MAX_VALUE
  for i in range(ntrace/10):
    if i%100==0:
      print "i =",i
      print "i2:  min =",i2min," max =",i2max
      print "i3:  min =",i3min," max =",i3max
    af.readInts(h)
    #i2,i3 = h[49],h[50] # Parihaka
    i2,i3 = h[47],h[48] # F3D
    if 0<=i2 and i2<m2 and 0<=i3 and i3<m3:
      m[i3][i2] = 1.0
    #print "i =",i," i2 =",i2," i3 =",i3
    if i2<i2min: i2min = i2
    if i2>i2max: i2max = i2
    if i3<i3min: i3min = i3
    if i3>i3max: i3max = i3
    af.skipBytes(bps*n1)
    af.seek(af.filePointer+100*(nthed+bps*n1))
  af.close()
  print "i2:  min =",i2min," max =",i2max
  print "i3:  min =",i3min," max =",i3max
  sp = SimplePlot()
  pv = sp.addPixels(m)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  aos = ArrayOutputStream(mapfile)
  aos.writeFloats(m)
  aos.close()


def dumpTraceHeaders(sgyfile,fmt,n1):
  af = ArrayFile(sgyfile,"r")
  af.skipBytes(nhead)
  af.skipBytes(nbhed)
  hi = zeroint(nthed/4)
  hs = zeroshort(nthed/2)
  for i in range(5):
    fp = af.getFilePointer()
    af.readInts(hi)
    af.seek(fp)
    af.readShorts(hs)
    print "ensemble number =",hi[5]
    print "trace in ensemble =",hi[6]
    print "coord scale factor =",hs[35]
    print "x,y coord =",hi[45],hi[46]
    print "iline,xline =",hi[47],hi[48]
    print "iline,xline =",hi[49],hi[50]
    #dump(hi)
    #dump(hs)
    if fmt==3:
      af.skipBytes(2*n1)
    else:
      af.skipBytes(4*n1)
  af.close()

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

def convert3(n1,n2,n3,sgyfile,datfile):
  print "converting",(n2*n3),"traces"
  ais = ArrayInputStream(sgyfile)
  aos = ArrayOutputStream(datfile)
  ais.skipBytes(nhead+nbhed)
  x = zeroshort(n1)
  y = zerofloat(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      ais.skipBytes(nthed) # skip trace header
      ais.readShorts(x) # read trace samples
      IbmIeee.shortToFloat(x,y) # convert
      aos.writeFloats(y) # write trace samples
    print "converted: i3 =",i3
  ais.close()
  aos.close()

def display3d(x,clip=0):
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  print "x min =",min(x)," max =",max(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,x)
  if clip>0.0:
    ipg.setClips(-clip,clip)
  frame.setVisible(True)

def displayFile3d(datfile,n1,n2,n3,clip=0):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datfile)
  ais.readFloats(x)
  ais.close()
  display3d(x,clip)

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
  elif format==3:
    print "format = 3 = 2-byte two's complement integer"
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
