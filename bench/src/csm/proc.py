"""
Processing.
Author: Dave Hale, Colorado School of Mines
Version: 2012.05.20
"""
from imports import *

s1 = Sampling(4001,0.002,0.000) # time sampling
s2 = Sampling(342,1,954) # station sampling, sweep 1
s3 = Sampling(105,1,1003) # first shotpoint is 1003
n1,n2,n3 = s1.count,s2.count,s3.count
#shotDir = "/data/seis/csm/fc2012/"
#segdDir = "/data/seis/csm/fc2012/segd/test139/"
shotDir = "/data/seis/csm/fc2012/line141s10/"
segdDir = "/data/seis/csm/fc2012/segd/line141s10/"

#############################################################################
def main(args):
  #process()
  display()

def process():
  f = readData(s1,s2,s3,shotDir+"shots.dat")
  lowpass3(f)
  tpow3(f)
  #g = copy(f)
  #gain3(g)
  #for i3 in range(n3):
  #  plot(s1,s2,g[i3],title="Shot "+str(s3.getValue(i3)))
  gain3(f,hw=4000.0)
  muteAirwave(f)
  #killNearOffsets(f)
  taperEdges(f)
  killGroundRoll(f)
  writeData(f,shotDir+"shotsq.dat",bo=ByteOrder.LITTLE_ENDIAN)
  #gain3(f)
  #for i3 in range(n3):
  #  plot(s1,s2,f[i3],title="Shot "+str(s3.getValue(i3)))

def display():
  f = readData(s1,s2,s3,shotDir+"shotsq.dat",bo=ByteOrder.LITTLE_ENDIAN)
  #lowpass3(f)
  #tpow3(f)
  #gain3(f)
  sf = SimpleFrame()
  ip = sf.addImagePanels(f)
  ip.setPercentiles(1,99)
  #ip.setClips(-2.5,2.5)

def taperEdges(f):
  t1 = 50
  h = fillfloat(1.0,s1.count,s2.count)
  for i2 in range(n2):
    for i1 in range(0,t1+t1):
      h[i2][i1] = max(0.0,float(i1-t1)/t1)
    for i1 in range(n1-t1-t1,n1):
      h[i2][i1] = max(0.0,float(n1-t1-i1)/t1)
  for i3 in range(n3):
    mul(h,f[i3],f[i3])

def killNearOffsets(f):
  lkill = 3
  for i3 in range(n3):
    i2 = int(s3.getValue(i3)-s2.first)
    i2min = i2-lkill
    i2max = i2+lkill
    for i2 in range(i2min,i2max+1):
      #scale = max(0.0,1.0-sin(0.5*PI*(i2-i2min)/lkill))
      scale = 0.0
      mul(scale,f[i3][i2],f[i3][i2])

def killGroundRoll(f):
  v = 1.100 # velocity
  p = 1/v # slope
  ci = CubicInterpolator(
    CubicInterpolator.Method.LINEAR,4,
    [p-0.4,p-0.3,p+0.3,p+0.4],[1,0,0,1])
  s2 = Sampling(n2,0.015,0.0)
  fft = Fft(s1,s2)
  fft.setComplex(False)
  fft.setCenter2(True)
  fft.setPadding1(200)
  fft.setPadding2(100)
  sw = fft.getFrequencySampling1()
  sk = fft.getFrequencySampling2()
  nw,nk = sw.count,sk.count
  h = fillfloat(1.0,nw,nk)
  for ik in range(nk):
    k = sk.getValue(ik)
    for iw in range(nw):
      w = sw.getValue(iw)
      if w!=0.0:
        h[ik][iw] = min(1.0,ci.interpolate(abs(k/w)))
  h = cmplx(h,zerofloat(nw,nk))
  for i3 in range(n3):
    g = copy(f[i3])
    g = fft.applyForward(g)
    cmul(h,g,g)
    g = fft.applyInverse(g)
    copy(g,f[i3])

def muteAirwave(f):
  vair = 0.330
  lmute = 0.2/s1.delta
  nmute = 1+2*lmute
  for i3 in range(n3):
    for i2 in range(n2):
      f32 = f[i3][i2]
      offset = (s2.getValue(i2)-s3.getValue(i3))*0.015
      imute = int(abs(offset)/vair/s1.delta)
      i1min = max(0,imute-lmute)
      i1max = min(n1-1,imute+lmute)
      for i1 in range(i1min,i1max+1):
        f32[i1] = 0.0

def readData(s1,s2,s3,fileName,bo=ByteOrder.BIG_ENDIAN):
  n1,n2,n3 = s1.count,s2.count,s3.count
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,bo)
  ais.readFloats(f)
  ais.close()
  return f

def writeData(flist,fileName,bo=ByteOrder.BIG_ENDIAN):
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

def gain3(f,hw=40.0):
  ref = RecursiveExponentialFilter(hw)
  for f3 in f:
    if max(abs(f3))>0.0:
      g = mul(f3,f3)
      ref.apply1(g,g)
      div(f3,add(0.0001,sqrt(g)),f3)

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
  sp.setVLimits(0.0,4.0)
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

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
