"""
Processing.
Author: Dave Hale, Colorado School of Mines
Version: 2012.05.20
---
Receiver stations: 954 - 1295 ( 954 <=> 0.000)
Source stations:  1003 - ???? (1003 <=> 7.350)
"""
from imports import *

s1 = Sampling(4001,0.002,0.000) # time sampling
s2 = Sampling(342,0.015,0.000) # receiver sampling (first group at 954)
s3 = Sampling(215,0.015,0.735) # shot sampling (first shot at 1003)
#s3 = Sampling(1,0.015,0.735)
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
#shotDir = "/data/seis/csm/fc2012/"
#segdDir = "/data/seis/csm/fc2012/segd/test139/"
shotDir = "/data/seis/csm/fc2012/line141s10/"
segdDir = "/data/seis/csm/fc2012/segd/line141s10/"

#############################################################################
def main(args):
  #process()
  display()

def process():
  f = readData(shotDir+"shotsp.dat")
  #lowpass(35.0,f)
  tpow(f)
  balance(f)
  #g = copy(f)
  #for i3 in range(n3):
  #  plot(g[i3],title="Shot at "+str(s3.getValue(i3)))
  muteAirWave(f)
  taperEdges(f)
  removeSlowWaves(f)
  #muteFirstBreak(f)
  #balance(f)
  #for i3 in range(n3):
  #  plot(f[i3],title="Shot at "+str(s3.getValue(i3)))
  writeData(f,shotDir+"shotsq.dat")

def display():
  f = readData(shotDir+"shotsq.dat")
  sf = SimpleFrame()
  ip = sf.addImagePanels(f)
  ip.setPercentiles(1,99)
  #ip.setClips(-2.5,2.5)

def balance(f):
  mf = MedianFinder(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      ma = mf.findMedian(abs(f[i3][i2]))
      if ma==0.0:
        ma = 0.00001
      div(f[i3][i2],ma,f[i3][i2])

def taperEdges(f):
  t1 = 50
  h = fillfloat(1.0,n1,n2)
  for i2 in range(n2):
    for i1 in range(0,t1+t1):
      h[i2][i1] = max(0.0,float(i1-t1)/t1)
    for i1 in range(n1-t1-t1,n1):
      h[i2][i1] = max(0.0,float(n1-t1-i1)/t1)
  for i3 in range(n3):
    mul(h,f[i3],f[i3])

def muteAirWave(f):
  vel = 0.330 # km/s
  lmute = 0.2/d1
  nmute = 1+2*lmute
  for i3 in range(n3):
    for i2 in range(n2):
      f32 = f[i3][i2]
      offset = s2.getValue(i2)-s3.getValue(i3)
      imute = s1.indexOfNearest(abs(offset)/vel)
      i1min = max(0,imute-lmute)
      i1max = min(n1-1,imute+lmute)
      for i1 in range(i1min,i1max+1):
        f32[i1] = 0.0

def muteFirstBreak(f):
  vel = 4.000 # km/s
  kmute = s1.indexOfNearest(0.3)
  for i3 in range(n3):
    for i2 in range(n2):
      f32 = f[i3][i2]
      offset = s2.getValue(i2)-s3.getValue(i3)
      imute = s1.indexOfNearest(abs(offset)/vel)
      for i1 in range(0,kmute+imute):
        f32[i1] = 0.0

def muteNearOffsets(f):
  lkill = 3
  for i3 in range(n3):
    i2 = s2.indexOfNearest(s3.getValue(i3))
    i2min = i2-lkill
    i2max = i2+lkill
    for i2 in range(i2min,i2max+1):
      #scale = max(0.0,1.0-sin(0.5*PI*(i2-i2min)/lkill))
      scale = 0.0
      mul(scale,f[i3][i2],f[i3][i2])

"""
refracted shear?
shot 116
321-93: 0.456 s
155-102: 0.795 km
vel = 1.75
"""

def removeSlowWaves(f):
  #vgr = 1.1 # ground-roll velocity
  vgr = 0.1 # ground-roll velocity
  vrs = 2.3 # refracted shear wave?
  slopeFilter(1.0/vrs,1.0/vgr,f)

def slopeFilter(pmin,pmax,f):
  ci = CubicInterpolator(
    CubicInterpolator.Method.LINEAR,4,
    [pmin-0.1,pmin,pmax,pmax+0.1],[1,0,0,1])
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

def readData(fileName,bo=ByteOrder.LITTLE_ENDIAN):
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,bo)
  ais.readFloats(f)
  ais.close()
  return f

def writeData(f,fileName,bo=ByteOrder.LITTLE_ENDIAN):
  n3 = len(f)
  print "writing",n3," shot records to",fileName
  aos = ArrayOutputStream(fileName,bo)
  for i3 in range(n3):
    print "  writing i3 =",i3
    aos.writeFloats(f[i3])
  print "  closing ..."
  aos.close()
  print "  done"

def tpow(f):
  t = rampfloat(f1,d1,0.0,n1,n2) # time
  mul(t,t,t) # time squared
  for f3 in f:
    mul(t,f3,f3)

def gain(f,hw=40.0):
  ref = RecursiveExponentialFilter(hw)
  for f3 in f:
    if max(abs(f3))>0.0:
      g = mul(f3,f3)
      ref.apply1(g,g)
      div(f3,add(0.0001,sqrt(g)),f3)

def lowpass(f3db,f):
  bf = ButterworthFilter(f3db*d1,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def plot(f,title=None):
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

def plotAmp(f,title=None):
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
