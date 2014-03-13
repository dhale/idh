#############################################################################
# Demo wavelet estimation from warping.
# This version minimizes || f - HSLAg ||, using iterative alternating
# updates to H and A, where HA need not equal the identity operator I.
# (1) find h such that f ~ SLAg * h (where h is a shaping filter)
# (2) find a such that f ~ HSLGa (solving non-Toeplitz system)

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import WaveletWarpingHA

#############################################################################

#pngDir = "./png/sinow/"
pngDir = None

def main(args):
  #goSimpleTest()
  goSino()

def goSino():
  na,ka = 41,-20 # sampling for inverse A of wavelet in PS image
  nh,kh = 41,-20 # sampling for wavelet H in PP image
  nt,dt,ft = 501,0.004,0.000 # used for plotting only
  nx,dx,fx = 721,0.015,0.000
  sa = Sampling(na,dt,ka*dt)
  sh = Sampling(nh,dt,kh*dt)
  st = Sampling(nt,dt,ft)
  sx = Sampling(nx,dx,fx)
  tmin,tmax = 100,400 # PP time window
  sfac = 1.000 # stabilization factor
  wha = 0.000 # weight for HA = I terms
  ha = True # True to find h, then a, then h, ...; False, for a, h, a, ...
  niter = 0
  f,g,u = getSinoImages() # PP image, PS image, and warping u(t,x)
  ww = WaveletWarpingHA()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  ww.setWeightHA(0.0)
  slg = ww.applyS(u,ww.applyL(u,g)) # PS warping without wavelets
  e1 = ww.rms(sub(f,slg))
  for ha in [True,False]:
    suffix = ""
    if ha:
      ag = zerofloat(na); ag[-ka] = 1.0 # initial inverse a in g
      hf = ww.getWaveletH(nh,kh,na,ka,ag,u,f,g) # wavelet h in f
      for jiter in range(niter):
        ag = ww.getInverseA(na,ka,nh,kh,hf,u,f,g) # inverse a in g
        hf = ww.getWaveletH(nh,kh,na,ka,ag,u,f,g) # wavelet h in f
    else:
      hf = zerofloat(nh); hf[-kh] = 1.0 # initial wavelet h in f
      ag = ww.getWaveletH(na,ka,nh,kh,hf,u,f,g) # inverse a in g
      for jiter in range(niter):
        hf = ww.getWaveletH(nh,kh,na,ka,ag,u,f,g) # wavelet h in f
        ag = ww.getInverseA(na,ka,nh,kh,hf,u,f,g) # inverse a in g
    hg = ww.getWaveletH(na,ka,ag,nh,kh) # wavelet in g
    hslag = ww.applyHSLA(na,ka,ag,nh,kh,hf,u,g) # PS warping with wavelets
    ew = ww.rms(sub(f,hslag))
    print "  e1 =",e1," ew =",ew
    rhslag = ww.rms(hslag); hslag = mul(hslag,1.0/rhslag)
    print "rhslag =",rhslag
    plotImage(st,sx,hslag,zoom=True,png="psww"+suffix)
    plotWaveletsPpPs(sh,hf,hg,png="wavelets"+suffix)
  rf = ww.rms(f); f = mul(f,1.0/rf)
  rslg = ww.rms(slg); slg = mul(slg,1.0/rslg)
  print "rf =",rf," rslg =",rslg
  plotImage(st,sx,f,zoom=True,png="pp")
  plotImage(st,sx,slg,zoom=True,png="psw1")

def goSimpleTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  #na,ka = 11,-5 # sampling for inverse wavelet A
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.00
  wha = 10.0
  st = Sampling(nt,dt,ft)
  for mp in [True,False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for r in [0.5,2.0]: # 0.5 for stretch; 2.0 for squeeze
      aw = zerofloat(na); aw[-ka] = 1.0
      hw = zerofloat(nh); hw[-kh] = 1.0
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)
      u = rampfloat(0.0,r,nt)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      ww = WaveletWarpingHA()
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ww.setWeightHA(wha)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      #dump(ak)
      for iter in range(500):
        aw = ww.getInverseA(na,ka,nh,kh,hw,u,f,g) # estimated inverse
        hw = ww.getWaveletH(nh,kh,na,ka,aw,u,f,g) # estimated wavelet
        #dump(aw)
      #hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
      sg = ww.applyS(u,g)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      hslag = ww.applyH(nh,kh,hw,slag)
      normalizeMax(hw)
      normalizeMax(hk)
      title = "r = "+str(r)
      #plotSequences(st,[f,g],labels=["f","g"],title=title)
      #plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title)

def getSinoImages():
  dataDir = "/data/seis/sino/warp/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  u = readImage(dataDir+"shifts.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  gain(100,f)
  gain(100,g)
  return f,g,u

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def normalizeMax(f):
  return mul(1.0/max(abs(f)),f,f)

def normalizeRms(f):
  mul(1.0/rms(f),f,f)
def rms(f):
  return sqrt(sum(mul(f,f))/len(f)/len(f[0]))

def makeImpulses(r,nt,ni):
  p = zerofloat(nt)
  q = zerofloat(nt)
  tmax = nt-1
  if r<=1.0:
    ts = rampfloat(tmax/(ni+1),tmax/(ni+1),ni)
  else:
    ts = rampfloat(tmax/(ni+1)/r,tmax/(ni+1)/r,ni)
  si = SincInterp.fromErrorAndFrequency(0.01,0.45)
  rj = -1.0
  for ji in range(ni):
    tj = ts[ji]
    tp = ts[ji]
    tq = r*tp
    if r<=1.0:
      tq += (1.0-r)*tmax
    #print "tp =",tp," tq =",tq
    rj = -rj
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
  return p,q

def addWavelet(fpeak,decay,p,mp=False):
  w = 2.0*PI*fpeak
  if not mp:
    decay *= 2.0
    w -= 2.0*PI*0.04
  r = exp(-decay)
  a1,a2 = -2.0*r*cos(w),r*r
  #print "a =",[1,a1,a2]
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(p,t)
  if not mp:
    w = 2.0*PI*(fpeak+0.04)
    poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
    zeros = []
    gain = 1.0
    rcf = RecursiveCascadeFilter(poles,zeros,gain)
    rcf.applyReverse(t,x)
  else:
    copy(t,x)
  conv(2,0,[1.0,-0.95],len(x),0,copy(x),len(x),0,x) # attenuate DC
  return x

def getWavelet(fpeak,decay,nh,kh,mp=False):
  x = zerofloat(nh)
  x[-kh] = 1.0
  return addWavelet(fpeak,decay,x,mp)

def plotSequence(x,xmax=None,title=None):
  sp = SimplePlot.asPoints(x)
  if xmax==None:
    xmax = max(abs(max(x)),abs(min(x)))
    xmax *= 1.05
  sp.setVLimits(-xmax,xmax)
  if title:
    sp.setTitle(title)

def plotSequences(st,xs,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(ix,0,st,xi)
    if labels:
      pp.setVLabel(ix,labels[ix])
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setVisible(True)
  if title:
    pp.setTitle(title)

def plotWavelets(st,hs,hmax=None,title=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  cs = [Color.BLACK,Color.GRAY,Color.LIGHT_GRAY]
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if hs[ih]:
      #pv = sp.addPoints(st,hs[ih])
      #pv.setLineStyle(ls[ih])
      sv = sp.addSequence(st,hs[ih])
      sv.setColor(cs[ih])
      #pv.setLineWidth(2)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  sp.setVLabel("Amplitude (normalized)")
  sp.setHLabel("Time (s)")
  if title:
    sp.setTitle(title)

def plotImage(st,sx,f,zoom=False,png=None):
  wpt = 240
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setFontSizeForPrint(8,wpt)
  sp.setSize(350,450)
  pv = sp.addPixels(st,sx,f)
  pv.setClips(-2.0,2.0)
  if zoom:
    sp.setVLimits(0.400,1.600) # zoom consistent with tmin,tmax = 100,400
  sp.setHLabel("Distance (km)")
  sp.setVLabel("PP time (s)")
  if pngDir and png:
    sp.paintToPng(360,wpt/72.0,pngDir+png+".png")

def plotWaveletsPpPs(st,h1,h2,png=None):
  wpt = 240
  pp = PlotPanel(2,1)
  h1 = mul(h1,1.0/max(abs(h1)))
  h2 = mul(h2,1.0/max(abs(h2)))
  sv1 = pp.addSequence(0,0,st,h1)
  sv2 = pp.addSequence(1,0,st,h2)
  pp.setVLimits(0,-0.65,1.05)
  pp.setVLimits(1,-0.65,1.05)
  pp.setVLabel(0,"PP wavelet")
  pp.setVLabel(1,"PS wavelet")
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(350,450)
  pf.setFontSizeForPrint(8,wpt)
  pf.setVisible(True)
  if pngDir and png:
    pf.paintToPng(720,wpt/72.0,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
