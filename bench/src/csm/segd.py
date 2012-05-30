"""
Reads segd files from CSM field camp, assuming these are in Sercel's 
SEG-D Rev 1 format, and writes a dat file containing a 3D array of 
floats. The byte order for floats in the dat file is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Version: 2012.05.09
"""
from imports import *

"""
Line10, 2011
station spacing is 30 m
line was shot east to west
Vibrator A is 3 stations to the east of channel 1
Vibrator B is 81 stations to the west of channel 120
"""
"""
Line10, 2012
station spacing is 15 m
line was shot SW to NE
shots from 1003 - 1217
Vibrator A is 3 stations to the east of channel 1
Vibrator B is 81 stations to the west of channel 120
"""
s1 = Sampling(4001,0.002,0.000) # time sampling
s2 = Sampling(342,1,954) # station sampling, sweep 1
s3 = Sampling(215,1.0,1003) # first shotpoint is 1003
#s2 = Sampling(277,0.015,-0.030) # channel sampling
#s3 = Sampling(1,1.0,1001.0) # shotpoint station sampling A
#shotDir = "/data/seis/csm/fc2012/"
#segdDir = "/data/seis/csm/fc2012/segd/test139/"
shotDir = "/data/seis/csm/fc2012/line141s10/"
segdDir = "/data/seis/csm/fc2012/segd/line141s10/"

#############################################################################
def main(args):
  #readLine141Segd()
  displayLine141()
  #displayLine140S1()
  #readLine140Segd()
  #readTestSegd()

def readLine141Segd():
  #global s3
  segdList = File(segdDir).listFiles() # list of segd files
  nshot = len(segdList)-3 # ignore first 3 shots
  #s3 = Sampling(nshot,1,1003) # first shotpoint is 1003
  g = zerofloat(s1.count,s2.count,s3.count)
  #print segdList
  for segdFile in segdList[3:]:
    print segdFile
    sl,sp,rpf,rpl,f = readSegd(segdFile)
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    i3 = int(sp-s3.first)
    zero(f[42]) # no geophone string at station 996
    copy(f,g[i3])
    #lowpass2(f)
    #tpow2(f)
    #gain2(f)
    #plot(s1,s2,f,title="Shot "+str(sp))
    #plotAmp(s1,s2,f,title="Shot "+str(sp))
  writeData(g,shotDir+"shotsp.dat")

def displayLine141():
  f = readData(s1,s2,s3,shotDir+"shotsp.dat")
  lowpass3(f)
  tpow3(f)
  gain3(f)
  sf = SimpleFrame()
  ip = sf.addImagePanels(f)
  #ip.setClips(-2.5,2.5)

def displayLine140S1():
  f = readData(s1,s2,s3,shotDir+"shotsp.dat")
  lowpass3(f)
  tpow3(f)
  gain3(f)
  sf = SimpleFrame()
  ip = sf.addImagePanels(f)
  #ip.setClips(-2.5,2.5)

def readLine140Segd():
  segdList = File(segdDir).listFiles() # list of segd files
  nshot = len(segdList)
  g = zerofloat(s1.count,s2.count,nshot)
  #print segdList
  ishot = 0
  for segdFile in segdList[:]:
    print segdFile
    sl,sp,rpf,rpl,f = readSegd(segdFile)
    #s2 = Sampling(len(f),0.015,-0.030) # offset sampling
    #s2 = Sampling(len(f),1,954) # station sampling
    #tpow2(f)
    #lowpass2(f)
    #gain2(f)
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    copy(f,g[ishot])
    ishot += 1
    #plot(s1,s2,f,title="Shot "+str(ishot))
    #plotAmp(s1,s2,f,title="Test "+str(ishot))
  #sf = SimpleFrame()
  #ip = sf.addImagePanels(g)
  #ip.setPercentiles(2,98)
  writeData(g,shotDir+"shots.dat")
  writeData(g,shotDir+"shotsp.dat",bo=ByteOrder.LITTLE_ENDIAN)

def readTestSegd():
  segdList = File(segdDir).listFiles() # list of segd files
  #print segdList
  itest = 0
  for segdFile in segdList:
    print segdFile
    sl,sp,rpf,rpl,f = readSegd(segdFile)
    s1 = Sampling(len(f[0]),0.002,0.000) # time sampling
    #s2 = Sampling(len(f),0.015,-0.030) # offset sampling
    s2 = Sampling(len(f),1,1001) # station sampling
    tpow2(f)
    lowpass2(f)
    gain2(f)
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    itest += 1
    plot(s1,s2,f,title="Test "+str(itest))
    #plotAmp(s1,s2,f,title="Test "+str(itest))

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

def tpow2(f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0,0.002,0.0,n1,n2) # time
  mul(t,t,t) # time squared
  return mul(t,f)

def tpow3(f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  t = rampfloat(s1.first,s1.delta,0.0,n1,n2) # time
  mul(t,t,t) # time squared
  for f3 in f:
    mul(t,f3,f3)

def gain2(f):
  ref = RecursiveExponentialFilter(40.0)
  for f2 in f:
    if max(abs(f2))>0.0:
      g = mul(f2,f2)
      ref.apply1(g,g)
      div(f2,sqrt(g),f2)

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
  sp.setSize(900,900)
  sp.setVLabel("Time (s)")
  if s2.delta==1.0:
    sp.setHLabel("Station")
  else:
    sp.setHLabel("Offset (km)")
  sp.setVLimits(0.0,8.0)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(s1,s2,f)
  #pv.setColorModel(ColorMap.BLUE_WHITE_RED)
  pv.setPercentiles(1,99)
  #pv.setClips(-2.5,2.5)

def plotAmp(s1,s2,f,title=None):
  fft = Fft(s1)
  sf = fft.getFrequencySampling1()
  ff = zerofloat(sf.count,s2.count)
  for i2 in range(s2.count):
    ff[i2] = cabs(fft.applyForward(f[i2]))
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(750,1000)
  sp.setSize(900,900)
  sp.setVLabel("Frequency (Hz)")
  if s2.delta==1.0:
    sp.setHLabel("Station")
  else:
    sp.setHLabel("Offset (km)")
  sp.setVLimits(0.0,120.0)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(sf,s2,ff)
  pv.setColorModel(ColorMap.JET)
  pv.setPercentiles(1,99)
  #pv.setClips(-2.5,2.5)

def readSegd(segdFile):
  #n1,n2 = 4001,230 # number of samples, number of traces
  #n1,n2 = 4001,277 # number of samples, number of traces (1 sweep)
  n1,n2 = 4001,342 # number of samples, number of traces (5 sweeps)
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
  #print "file =",segdFile
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

def readLine10Segd():
  segdList = File(segdDir).listFiles() # list of segd files
  n1,n2 = s1.count,s2.count
  fzeros = zerofloat(n1,n2)
  spfa = int(s3a.first) # first shot point for vib A
  spfb = int(s3b.first) # first shot point for vib B
  faList,fbList = [],[] # lists of shot records
  spa,spb = spfa-1,spfb-1 # shots last appended
  for segdFile in segdList:
    sl,sp,rpf,rpl,f = readSegd(segdFile)
    #print segdFile
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    if sl==10: # vibrator A
      if sp==spa: # if same station as before, count the last one
        faList.pop()
      for i in range(spa+1,sp): # if necessary, insert zero records
        faList.append(fzeros)
        print "a: zero ",i
      print "a: append ",sp
      spa = sp
      faList.append(f)
    elif sl==20: # vibrator B
      if sp==spb:
        fbList.pop()
      for i in range(spb+1,sp):
        fbList.append(fzeros)
        print "b: zero ",i
      print "b: append ",sp
      spb = sp
      fbList.append(f)
  na = len(faList)
  nb = len(fbList)
  print "na =",na," nb =",nb
  writeData(faList,shotDir+"shota.dat")
  writeData(fbList,shotDir+"shotb.dat")

def readSegd2011(segdFile):
  n1,n2 = 3001,120 # number of samples, number of traces
  f = zerofloat(n1,n2) # the traces
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
  #print "file =",segdFile
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
  for ict in range(nct): # for all channels (including aux channels)
    ais.readBytes(th) # trace header
    cn = th[3] # channel set number
    ic = bcd2(th,4) # channel (trace) number
    ais.readBytes(the) # trace header extension 1
    rln = bin3(the,0) # receiver line number
    rpn = bin3(the,3) # receiver point number
    if ic==1:
      rpf = rpn
    elif ic==120:
      rpl= rpn
    ais.skipBytes(6*len(the)) # skip trace header extensions 2-7
    if cn==cns: # if seismic channel, ...
      #print "ic =",ic," rln =",rln," rpn =",rpn
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
