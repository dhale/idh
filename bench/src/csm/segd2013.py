"""
Reads segd files from CSM field camp, assuming these are in Sercel's 
SEG-D Rev 1 format, and writes a dat file containing a 3D array of 
floats. The byte order for floats in the dat file is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Version: 2013.10.15
"""
from imports import *

s1 = Sampling(3001,0.001,0.000) # time sampling
s2 = Sampling(309,0.010,0.0) # station sampling
s3 = Sampling(215,1.0,1003) # first shotpoint is 1003
#s2 = Sampling(277,0.015,-0.030) # channel sampling
#s3 = Sampling(1,1.0,1001.0) # shotpoint station sampling A
segdDir = "/data/seis/csm/fc2013/segd/139/"

#############################################################################
def main(args):
  readAndPlotSegd()

def readAndPlotSegd():
  segdList = File(segdDir).listFiles() # list of segd files
  nshot = len(segdList)-3 # ignore first 3 shots
  #s3 = Sampling(nshot,1,1003) # first shotpoint is 1003
  #print segdList
  i3 = 40
  for segdFile in segdList[i3:i3+1]:
    print segdFile
    sl,sp,rpf,rpl,f = readSegd(segdFile)
    sl,sp = int(sl),int(sp)
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    s2 = Sampling(len(f),0.010,(rpf-sp)*0.010)
    gain2(f)
    g = ppf2(f)
    gb = ppf2(f,True)
    #gain2(g)
    #gain2(gb)
    plot(s1,s2,f,title="Shot "+str(sp)+": raw")
    plot(s1,s2,g,title="Shot "+str(sp)+": ppf")
    plot(s1,s2,gb,title="Shot "+str(sp)+": ppfb")
    plotAmp(s1,s2,f,title="Shot "+str(sp)+": raw")
    plotAmp(s1,s2,g,title="Shot "+str(sp)+": ppf")
    plotAmp(s1,s2,gb,title="Shot "+str(sp)+": ppfb")

def ppf2(f,b=False):
  n1,n2 = len(f[0]),len(f)
  m1,m2 = 40,40
  f1=[0.000, 0.120, 0.120, 0.000]
  f2=[0.000,-0.350, 0.350, 0.000]
  g = copy(f)
  if not b:
    ppf = PolygonPassFilter(n1,n2,m1,m2,f1,f2)
  else:
    ppf = PolygonPassFilterB(n1,n2,f1,f2)
  ppf.apply(f,g)
  return g

def tpow2(f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0,0.002,0.0,n1,n2) # time
  t = pow(t,2.0)
  return mul(t,f)

def tpow3(f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  t = rampfloat(s1.first,s1.delta,0.0,n1,n2) # time
  mul(t,t,t) # time squared
  for f3 in f:
    mul(t,f3,f3)

def gain2(f):
  ref = RecursiveExponentialFilter(100.0)
  for f2 in f: # for all traces, ...
    if max(abs(f2))>0.0:
      g = mul(f2,f2) # square the trace amplitudes
      ref.apply1(g,g) # smooth
      div(f2,sqrt(g),f2) # normalizes

def gain3(f):
  ref = RecursiveExponentialFilter(40.0)
  for f3 in f:
    if max(abs(f3))>0.0:
      g = mul(f3,f3)
      ref.apply1(g,g)
      div(f3,sqrt(g),f3)

def lowpass2(f):
  f3db = 25.0*0.002
  #f3db = 35.0*0.002
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def lowpass3(f):
  bf = ButterworthFilter(0.05,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def plot(s1,s2,f,title=None):
  print "plot f: min =",min(f),"max =",max(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(750,1000)
  sp.setSize(350,700)
  sp.setVLabel("Time (s)")
  if s2.delta==1.0:
    sp.setHLabel("Station")
  else:
    sp.setHLabel("Offset (km)")
  sp.setVLimits(0.0,1.0)
  #sp.setHLimits(s2.first,s2.first+1.5)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(s1,s2,f)
  #pv.setColorModel(ColorMap.BLUE_WHITE_RED)
  pv.setPercentiles(02,98)
  #pv.setClips(-2.5,2.5)

def plotAmp(s1,s2,f,title=None):
  #s1 = Sampling(s1.count,1.0,0.0)
  #s2 = Sampling(s2.count,1.0,0.0)
  fft = Fft(s1,s2)
  fft.setCenter(True)
  sf1 = fft.getFrequencySampling1()
  sf2 = fft.getFrequencySampling2()
  a = cabs(fft.applyForward(f))
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(700,700)
  sp.setVLabel("Frequency (Hz)")
  sp.setHLabel("Wavenumber (cycles/km)")
  #sp.setVLimits(0.0,60.0)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(sf1,sf2,a)
  pv.setColorModel(ColorMap.JET)
  pv.setPercentiles(1,99)
  #pv.setClips(-2.5,2.5)

def readData(s1,s2,s3,fileName,bo=ByteOrder.LITTLE_ENDIAN):
  n1,n2,n3 = s1.count,s2.count,s3.count
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,bo)
  ais.readFloats(f)
  ais.close()
  return f

def writeData(flist,fileName,bo=ByteOrder.LITTLE_ENDIAN):
  n3 = len(flist)
  print "writing",n3," shot records to",fileName
  aos = ArrayOutputStream(fileName,bo)
  for f in flist:
    aos.writeFloats(f)
  aos.close()

def readSegd(segdFile):
  n1,n2 = 3001,309 # number of samples, number of traces
  gh = zerobyte(32) # general header
  th = zerobyte(20) # trace header
  the = zerobyte(32) # trace header extension
  csh = zerobyte(32) # channel set header
  ais = ArrayInputStream(segdFile,ByteOrder.BIG_ENDIAN)
  ais.readBytes(gh) # general header 1
  fn = bcd2(gh,0) # file number
  ais.readBytes(gh) # general header 2
  ais.readBytes(gh) # general header 3
  sln = bin5(gh,3) # source line number
  spn = bin5(gh,8) # source point number
  print "file =",segdFile
  #print "fn = ",fn," sln =",sln," spn =",spn
  cns = 0 # channel set number for seismic traces
  nct = 0 # total number of channels, including aux channels
  for ics in range(16): # for each channel set header, ...
    ais.readBytes(csh) # read channel set header
    cn = csh[1] # channel set number
    ct = (csh[10]>>4)&0xf # channel type (in high 4 bits)
    nc = bcd2(csh,8) # number of channels
    if nc>0: # if we have channels of this type, ...
      #print "cn =",cn," nc =",nc," ct =",ct
      if ct==1: # if seismic, ...
        cns = cn # remember channel set number for seismic
        ncs = nc # remember number of seismic channels
      nct += nc # count total number of channels
  #print "nct =",nct,"cns =",cns
  ais.skipBytes(1024) # skip extended header
  ais.skipBytes(1024) # skip external header
  rpf = 1
  rpl = 1
  f = None
  for ict in range(nct): # for all channels (including aux channels)
    ais.readBytes(th) # trace header
    cn = th[3] # channel set number
    ic = bcd2(th,4) # channel (trace) number
    ais.readBytes(the) # trace header extension 1
    rln = bin3(the,0) # receiver line number
    rpn = bin3(the,3) # receiver point number
    n1 = bin3(the,7) # number of samples
    #print "ic =",ic," rln =",rln," rpn =",rpn," n1 =",n1
    if ic==1:
      rpf = rpn
    elif ic==n2:
      rpl = rpn
    ais.skipBytes(6*len(the)) # skip trace header extensions 2-7
    if cn==cns: # if seismic channel, ...
      #print "ic =",ic," rln =",rln," rpn =",rpn
      if not f:
        f = zerofloat(n1,n2) # the traces
      ais.readFloats(f[ic-1]) # get the seismic trace
    else:
      ais.skipBytes(4*n1) # skip the auxiliary trace
  ais.close()
  f = mul(1.0e-14,f) # scale values to approximate range [-10,10]
  return sln,spn,rpf,rpl,f

def bcd2(b,k):
  """ Returns binary-coded decimal integer from bytes k,k+1 in b."""
  return (1000*((b[k  ]>>4)&0xf)+100*(b[k  ]&0xf)+
            10*((b[k+1]>>4)&0xf)+  1*(b[k+1]&0xf))

def displayLine10(vib):
  if vib=="A":
    f = readData(s1,s2,s3a,shotDir+"shota.dat")
  elif vib=="B":
    f = readData(s1,s2,s3b,shotDir+"shotb.dat")
  lowpass3(f)
  tpow3(f)
  gain3(f)
  sf = SimpleFrame()
  ip = sf.addImagePanels(f)
  ip.setPercentiles(2,98)

def bin3(b,k):
  """ Returns binary integer from bytes k,k+1,k+2 in b."""
  b0 = b[k  ]
  b1 = b[k+1]
  b2 = b[k+2] 
  if b0<0: b0 += 256
  if b1<0: b1 += 256
  if b2<0: b2 += 256
  return (b0<<16)|(b1<<8)|(b2)

def bin5(b,k):
  """ Returns binary integer from bytes k,k+1,...,k+4 in b."""
  b0 = b[k  ]
  b1 = b[k+1]
  b2 = b[k+2] 
  b3 = b[k+3] 
  b4 = b[k+4] 
  if b0<0: b0 += 256
  if b1<0: b1 += 256
  if b2<0: b2 += 256
  if b3<0: b3 += 256
  if b4<0: b4 += 256
  return b0*65536.0+b1*256.0+b2+b3/256.0+b4/65536.0

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
