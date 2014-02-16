#############################################################################
# Demo estimation of two wavelets from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import WaveletWarping2

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  goSimpleTest()

def goSimpleTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  fpeakf,decayf = 0.08,0.10 # peak frequency and decay for wavelet in f
  fpeakg,decayg = 0.05,0.08 # peak frequency and decay for wavelet in g
  #nab,kab = 81,-20 # sampling for inverse wavelets A and B
  nab,kab = 7,-3 # sampling for inverse wavelets A and B
  ncd,kcd = 181,-90 # sampling for wavelets C and D
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  for mp in [False]: # True, for minimum-phase; False for other
    ck = getWavelet(fpeakf,decayf,ncd,kcd,mp) # known wavelet c in f
    dk = getWavelet(fpeakg,decayg,ncd,kcd,mp) # known wavelet d in g
    for r in [0.5]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(fpeakf,decayf,p,mp)
      g = addWavelet(fpeakg,decayg,q,mp)
      u = rampfloat(0.0,r,nt)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      ww = WaveletWarping2()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ak,bk = ww.getWaveletCD(ncd,kcd,[ck,dk],nab,kab) # known inverses
      ak0 = ak[-kab]; ak = div(ak,ak0); bk = div(bk,ak0)
      aw,bw = ww.getInverseAB(nab,kab,u,f,g) # estimated inverse wavelets
      cw,dw = ww.getWaveletCD(nab,kab,[aw,bw],ncd,kcd) # estimated wavelets
      print "ak & aw:"; dump(ak); dump(aw)
      print "bk & bw:"; dump(bk); dump(bw)
      sg = ww.applyS(u,g)
      af = ww.applyA(nab,kab,aw,f)
      bg = ww.applyB(nab,kab,bw,g)
      lbg = ww.applyL(u,bg) # lowpass, if squeezing
      slbg = ww.applyS(u,lbg)
      waf = ww.applyW(af)
      wslbg = ww.applyW(slbg)
      cslbg = ww.applyC(ncd,kcd,cw,slbg)
      ncw = normalize(cw)
      ndw = normalize(dw)
      nck = normalize(ck)
      ndk = normalize(dk)
      title = "r = "+str(r)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      #plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      #plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[af,bg],labels=["Af","Bg"],title=title)
      plotSequences(st,[af,lbg],labels=["Af","LBg"],title=title)
      plotSequences(st,[af,slbg],labels=["Af","SLBg"],title=title)
      plotSequences(st,[waf,wslbg],labels=["WAf","WSLBg"],title=title)
      plotSequences(st,[f,cslbg],labels=["f","CSLBg"],title=title)
      plotWavelets(Sampling(ncd,dt,kcd*dt),[ncw,nck],title=title)
      plotWavelets(Sampling(ncd,dt,kcd*dt),[ndw,ndk],title=title)

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
  #conv(2,0,[1.0,-0.95],len(x),0,copy(x),len(x),0,x) # attenuate DC
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
