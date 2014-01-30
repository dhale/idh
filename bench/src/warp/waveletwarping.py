#############################################################################
# Demo wavelet estimation from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import WaveletWarping

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  goSimpleTest()

def goSimpleTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.1,0.10 # peak frequency and decay for wavelet
  na,ka = 11,-5 # sampling for inverse wavelet A
  nh,kh = 161,-80 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  for zp in [False,True]: # for minimum-phase and zero-phase wavelets, ...
    hk = getWavelet(freq,decay,nh,kh,zp) # known wavelet
    for r in [0.5,2.0]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(freq,decay,p,zp)
      g = addWavelet(freq,decay,q,zp)
      u = rampfloat(0.0,r,nt)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      ww = WaveletWarping()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
      hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
      dump(ak)
      dump(aw)
      sg = ww.applyS(u,g)
      af = ww.applyA(na,ka,aw,f)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      baf = ww.applyB(af)
      bslag = ww.applyB(slag)
      hslag = ww.applyH(nh,kh,hw,slag)
      nhw = normalize(hw)
      nhk = normalize(hk)
      title = "r = "+str(r)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
      plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
      plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
      plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title)

def normalize(h):
  return div(h,max(max(h),-min(h)))

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
    print "tp =",tp," tq =",tq
    rj = -rj
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
  return p,q

def addWavelet(fpeak,decay,p,zp=False):
  if zp:
    decay *= 2.0
  r = exp(-decay)
  w = 2.0*PI*fpeak
  a1,a2 = -2.0*r*cos(w),r*r
  #print "a =",[1,a1,a2]
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(p,t)
  if zp:
    rcf.applyReverse(t,x)
  else:
    copy(t,x)
  return x

def getWavelet(fpeak,decay,nh,kh,zp=False):
  x = zerofloat(nh)
  x[-kh] = 1.0
  return addWavelet(fpeak,decay,x,zp)

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
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(ls[ih])
      #pv.setLineWidth(2)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  sp.setVLabel("Amplitude (normalized)")
  sp.setHLabel("Time (s)")
  if title:
    sp.setTitle(title)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
