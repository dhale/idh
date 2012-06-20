"""
Converts sgy image files (SEG-Y format) to dat files (3D arrays of floats)
Removes all SEG-Y headers and, if necessary, converts the data format to
IEEE floats. The byte order for floats in dat files is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Version: 2012.06.18
"""
from imports import *
global n1,n2,n3

#############################################################################
def main(args):
  #goMbs()
  goNorne()
  #goSino()
  #goF3d()
  #goParihaka()

def goMbs():
  """
  ***************************************************************************
  PstmSmall (pstm_fraw.sgy):
  number of bytes = 959341200
  number of traces = 153740
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  n1 =  1500 (number of samples per trace)
  n2 =   448 (number of traces in inline direction)
  n3 =   422 (number of traces in crossline direction)
  d1 = 0.002000 (time sampling interval, in s)
  d2 = 0.016764 (inline sampling interval, in km)
  d3 = 0.016764 (crossline sampling interval, in km)
  i2min =   601, i2max =  1048 (inline index bounds)
  i3min =  1001, i3max =  1422 (crossline index bounds)
  xmin = 823.725048, xmax = 831.574867 (x coordinate bounds, in km)
  ymin = 245.803217, ymax = 255.588516 (y coordinate bounds, in km)
  grid azimuth =  58.12 degrees
  ***************************************************************************
  PstmLarge (pstm_raw_cut.sgy):
  number of bytes = 4647376320
  number of traces = 795783
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  n1 =  1400 (number of samples per trace)
  n2 =  1119 (number of traces in inline direction)
  n3 =  1189 (number of traces in crossline direction)
  d1 = 0.002000 (time sampling interval, in s)
  d2 = 0.016764 (inline sampling interval, in km)
  d3 = 0.016764 (crossline sampling interval, in km)
  i2min =   350, i2max =  1468 (inline index bounds)
  i3min =   234, i3max =  1422 (crossline index bounds)
  xmin = 822.376308, xmax = 842.769562 (x coordinate bounds, in km)
  ymin = 235.175146, ymax = 255.588516 (y coordinate bounds, in km)
  grid azimuth =  58.12 degrees
  good subset:
    i1min,i1max = 150, 650,  m1 = 501
    i2min,i2max = 490,1258,  m2 = 763
    i3min,i3max = 358, 917,  m3 = 560
  """
  imageType = "PstmSmall" # which image to process
  firstLook = True # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # plots the image
  basedir = "/data/seis/mbs/"
  if imageType=="PstmSmall":
    sgyfile = basedir+"PstmSmall/Marathon20070228/pstm_fraw.sgy"
    datfile = basedir+"dat/pstm_fraw_s1.dat"
    i1min,i1max,i2min,i2max,i3min,i3max = 150,650,601,1048,1001,1422
  elif imageType=="PstmLarge":
    sgyfile = basedir+"PstmLarge/pstm_raw_cut.sgy"
    datfile = basedir+"dat/pstm_raw_s1.dat"
    i1min,i1max,i2min,i2max,i3min,i3max = 150,650,490,1258,358,917
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    si.writeFloats(datfile,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=10000.0)

def goNorne():
  """
  ***************************************************************************
  norne4d_2006-full.sgy:
  number of bytes = 1363689924
  number of traces = 321321
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  n1 =  1001 (number of samples per trace)
  n2 =  1001 (number of traces in inline direction)
  n3 =   321 (number of traces in crossline direction)
  d1 = 0.004000 (time sampling interval, in s)
  d2 = 0.012501 (inline sampling interval, in km)
  d3 = 0.012500 (crossline sampling interval, in km)
  i2min =  1300, i2max =  2300 (inline index bounds)
  i3min =   970, i3max =  1290 (crossline index bounds)
  xmin =  453.320000, xmax =  464.634000 (x coordinate bounds, in km)
  ymin = 7317.354880, ymax = 7329.340160 (y coordinate bounds, in km)
  grid azimuth =  41.80 degrees
  grid reference point:
    i2ref =  1300, i3ref =   970, x =  453.320000, y = 7320.021120
  grid corner points:
    i2min =  1300, i3min =   970, x =  453.320000, y = 7320.021120
    i2max =  2300, i3min =   970, x =  461.652000, y = 7329.340160
    i2min =  1300, i3max =  1290, x =  456.302000, y = 7317.354880
    i2max =  2300, i3max =  1290, x =  464.634000, y = 7326.673920
  ***************************************************************************
  Full Norne corner coordinates (Knut, 16/7-09)
  line	trace	X	Y
  970	2300	461652	7329340.16
  970	1300	453320	7320021.12
  1290	1300	456302	7317355
  1290	2300	464634	7326674
  ***************************************************************************
  E-Segment
  n1,nt = 1001, dt = 0.0040, ft = 0.0 (s)
  n2,ny =  401, dy = 0.0125, fy = 0.0 (km)
  n3,nx =  101, dx = 0.0125, fx = 0.0 (km)
  ***************************************************************************
  Full
  n1,nt = 1001, dt = 0.0040, ft = 0.0 (s)
  n2,ny = 1001, dy = 0.0125, fy = 0.0 (km)
  n3,nx =  321, dx = 0.0125, fx = 0.0 (km)
  """
  firstLook = True # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # plots the image
  basedir = "/data/seis/norne/"
  sgyfile = basedir+"sgy/norne4d_2006-full.sgy"
  datfile = basedir+"dat/full2006.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,1000,1300,2300,970,1290
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    si.writeFloats(datfile,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1000.0)

def goSino():
  """
  Two 2D (x- and z-component) images
  nt = 2001, dt = 0.004 s
  nx =  721, dx = 0.030 km?
  samples [0:1250] in z component ~ samples [0:2000] in x component
  """
  basedir = "/data/seis/sino/"
  for component in ["x","z"]:
    sgyfile = basedir+"260_"+component+"_201-921_stack.segy"
    datfile = basedir+component+"260.dat"
    i1min,i1max,i2min,i2max = 0,2000,201,921
    n1,n2 = 1+i1max-i1min,1+i2max-i2min
    si = SegyImage(sgyfile,ByteOrder.LITTLE_ENDIAN) # non-standard byte order!
    #si.printSummaryInfo()
    #si.printBinaryHeader()
    #si.printTraceHeader(0)
    #si.printTraceHeader(1)
    #plotIbmIeeeFloats(si)
    si.setD2(0.015) # a guess, from looking at group coordinates in headers
    if component=="x":
      si.setFormat(5) # formats appear to be IEEE for x and IBM for z!???
    si.printAllInfo()
    si.writeFloats(datfile,i1min,i1max,i2min,i2max,0,0)
    si.close()
    f = readImage(datfile,n1,n2)
    gain(500,f)
    if component=="z":
      stretch(1.6,f)
    show2d(f,title=component+" component")

def show2d(f,clip=None,title=None):
  print "show2d: f min =",min(f)," max =",max(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  if clip:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(2,98)
  if title:
    sp.setTitle(title)
  sp.setSize(600,1100)
  
def show3d(f,clip=None):
  print "show3d: f min =",min(f)," max =",max(f)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)

def plot23(si):
  i2 = si.getI2sAsFloats()
  i3 = si.getI3sAsFloats()
  sp = SimplePlot()
  sp.setHLabel("inline sample index i2")
  sp.setVLabel("crossline sample index i3")
  pv = sp.addPoints(i2,i3)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(i2,i3)
  sp.setSize(w,h)

def plotXY(si):
  x = si.getXs()
  y = si.getYs()
  sp = SimplePlot()
  sp.setHLabel("x coordinate (km)")
  sp.setVLabel("y coordinate (km)")
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(x,y)
  sp.setSize(w,h)

def goodWidthHeight(x,y):
  xmin,xmax = min(x),max(x)
  ymin,ymax = min(y),max(y)
  w,h = 1000,1000
  if (xmax-xmin)>(ymax-ymin):
    h = int(h*(ymax-ymin)/(xmax-xmin))
  else:
    w = int(w*(xmax-xmin)/(ymax-ymin))
  return w,h

def plotIbmIeeeFloats(si):
  ntrace = si.countTraces()
  fmt = si.getFormat()
  si.setFormat(1) # IBM floats
  fibm = si.getTrace(ntrace/2)
  si.setFormat(5) # IEEE floats
  fieee = si.getTrace(ntrace/2)
  si.setFormat(fmt)
  pp = PlotPanel(2,1)
  pp.setTitle("IBM (top) versus IEEE (bottom)")
  pp.setHLabel(0,"Sample index")
  pp.setVLabel(0,"IBM amplitude")
  pp.setVLabel(1,"IEEE amplitude")
  pp.addPoints(0,0,fibm)
  pp.addPoints(1,0,fieee)
  pf = PlotFrame(pp)
  pf.setSize(1000,800)
  pf.setVisible(True)

def lowpass(f3db,f):
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def stretch(c,f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0.0,1.0/c,n1)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  g = zerofloat(n1)
  for i2 in range(n2):
    si.setUniformSamples(f[i2])
    si.interpolate(n1,1.0/c,0.0,g)
    copy(g,f[i2])

def plot2(x,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(350,1100)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(x)
  pv.setPercentiles(2,98)
  
def display(datfile,clip=0.0):
  x = readImage(datfile,n1,n2,n3)
  print "x min =",min(x)," max =",max(x)
  frame = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  ipg = frame.addImagePanels(s1,s2,s3,x)
  if clip>0.0:
    ipg.setClips(-clip,clip)
  
def goF3d():
  """
  TIME: 462 samples (1.848 sec, interval 4 ms)
  nbytes = 699,003,060
  ntrace = (nbytes-nhead-nbhed)/(240+2*n1)
  """
  fmt = 3 # 2-byte shorts in [-32767,32767]
  global n1,n2,n3
  n1,n2,n3 = 462,951,591
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
  #makeMap(sgyfile,mapfile,fmt,n1)
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
  datdir = "/data/seis/nz/par/"
  sgyfile = datdir+"Parihaka3d_raw.sgy"
  mapfile = datdir+"map.dat"
  datfile = datdir+"par11.dat"
  #displayParihaka(datfile)
  #makeSubsetsParihaka(datdir)
  #displaySubsetParihaka(datfile,0,500,500,751,501,501)
  #bigSubsetParihaka(n1,sgyfile,datfile)
  readFormat(sgyfile)
  #makeMap(sgyfile,mapfile,n1)
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

def readImage(datfile,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
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
  ais = ArrayInputStream(sgyfile,bo)
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

def getTraceHeaderInfo(sgyfile,bo=ByteOrder.BIG_ENDIAN):
  fmt,bps,ntrace,nt,dt = getBinaryHeaderInfo(sgyfile,bo)
  hi = zeroint(240/4) # 240-byte trace header as 4-byte ints
  hs = zeroshort(240/2) # 240-byte trace header as 2-byte ints
  af = ArrayFile(sgyfile,"r")
  af.skipBytes(3200) # skip text file header
  af.skipBytes(400) # skip binary file header

  i2min =  Integer.MAX_VALUE
  i2max = -Integer.MAX_VALUE
  i3min =  Integer.MAX_VALUE
  i3max = -Integer.MAX_VALUE
  xmin =  Float.MAX_VALUE
  xmax = -Float.MAX_VALUE
  ymin =  Float.MAX_VALUE
  ymax = -Float.MAX_VALUE
  for itrace in range(ntrace):
    fp = af.getFilePointer()
    af.readInts(hi)
    af.seek(fp)
    af.readShorts(hs)

    print "ensemble number =",hi[5]
    print "trace in ensemble =",hi[6]
    print "coord scale factor =",hs[35]
    print "x,y coord =",hi[45],hi[46]
    print "iline,xline =",hi[47],hi[48]
    #print "iline,xline =",hi[49],hi[50]
    #dump(hi)
    #dump(hs)
    if fmt==3:
      af.skipBytes(2*n1)
    else:
      af.skipBytes(4*n1)
  af.close()

  


def makeMap(sgyfile,mapfile,fmt,n1):
  if fmt==3:
    bps = 2
  else:
    bps = 4
  nbytes = File(sgyfile).length
  ntrace = (nbytes-nhead-nbhed)/(240+bps*n1)
  af = ArrayFile(sgyfile,"r")
  af.skipBytes(nhead)
  af.skipBytes(nbhed)
  h = zeroint(nthed/4)
  #m2,m3 = 8000,16000 # Parihaka
  #m2,m3 = 800,1300 # F3D
  #m2,m3 = 800,1300 # Norne
  m2,m3 = 1422,1043 # Mbs
  m = zerofloat(m2,m3)
  i2min =  Integer.MAX_VALUE
  i2max = -Integer.MAX_VALUE
  i3min =  Integer.MAX_VALUE
  i3max = -Integer.MAX_VALUE
  for i in range(ntrace):
    #if i%100==0:
    #  print "i =",i
    #  print "i2:  min =",i2min," max =",i2max
    #  print "i3:  min =",i3min," max =",i3max
    af.readInts(h)
    #i2,i3 = h[49],h[50] # Parihaka
    i2,i3 = h[47],h[48] # F3D, Norne, Mbs
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
  if mapfile:
    aos = ArrayOutputStream(mapfile)
    aos.writeFloats(m)
    aos.close()


def dumpTraceHeaders(sgyfile,fmt,n1,ndump=5):
  af = ArrayFile(sgyfile,"r")
  af.skipBytes(nhead)
  af.skipBytes(nbhed)
  hi = zeroint(nthed/4)
  hs = zeroshort(nthed/2)
  for i in range(ndump):
    fp = af.getFilePointer()
    af.readInts(hi)
    af.seek(fp)
    af.readShorts(hs)
    print "ensemble number =",hi[5]
    print "trace in ensemble =",hi[6]
    print "coord scale factor =",hs[35]
    print "x,y coord =",hi[45],hi[46]
    print "iline,xline =",hi[47],hi[48]
    #print "iline,xline =",hi[49],hi[50]
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
  ais = ArrayInputStream(sgyfile,bo)
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
  ais = ArrayInputStream(sgyfile,bo)
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
    print "converted: i3 =",i3," n2 =",n2
  ais.close()
  aos.close()

def convert3(n1,n2,n3,sgyfile,datfile):
  print "converting",(n2*n3),"traces"
  ais = ArrayInputStream(sgyfile,bo)
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
  ais = ArrayInputStream(sgyfile,bo)
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
  ais = ArrayInputStream(sgyfile,bo)
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
