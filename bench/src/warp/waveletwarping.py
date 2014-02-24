#############################################################################
# Demo wavelet estimation from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import WaveletWarping

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  #goSimpleTest()
  goSino()

def goSino():
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  nt,dt,ft = 501,0.004,0.000 # used for plotting only
  nx,dx,fx = 721,0.015,0.000
  f,g,u = getSinoImages()
  st = Sampling(nt,dt,ft)
  sx = Sampling(nx,dx,fx)
  tmin,tmax = 88,363
  fmin,fmax = 0.0,0.5
  sfac = 1.000
  ww = WaveletWarping()
  ww.setFrequencyRange(fmin,fmax)
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
  hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
  nhw = normalize(hw)
  af = ww.applyA(na,ka,aw,f)
  ag = ww.applyA(na,ka,aw,g)
  slag = ww.applyS(u,ww.applyL(u,ag))
  slg = ww.applyS(u,ww.applyL(u,g))
  cmin,cmax = -0.5,0.5
  #plotImage(st,sx,haslg,fmin=cmin,fmax=cmax,zoom=True)
  #plotImage(st,sx,hslag,fmin=cmin,fmax=cmax,zoom=True)
  plotImage(st,sx,f,fmin=cmin,fmax=cmax,zoom=True)
  plotImage(st,sx,slg,fmin=cmin,fmax=cmax,zoom=True)
  plotImage(st,sx,af,fmin=cmin,fmax=cmax,zoom=True)
  plotImage(st,sx,slag,fmin=cmin,fmax=cmax,zoom=True)
  f0 = ww.rms(f)
  f1 = ww.rms(af)
  g0 = ww.rms(slg)
  g1 = ww.rms(slag)
  e0 = ww.rms(sub(f,slg))
  e1 = ww.rms(sub(af,slag))
  print "f0 =",f0," f1 =",f1
  print "g0 =",g0," g1 =",g1
  print "e0 =",e0," e1 =",e1
  #plotWavelets(Sampling(nh,dt,kh*dt),[nhw])

def goSimpleTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  for mp in [True,False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for r in [0.5,2.0]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)
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
      #dump(ak)
      #dump(aw)
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
      #plotSequences(st,[f,g],labels=["f","g"],title=title)
      #plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      #plotSequences(st,[f,g],labels=["f","g"],title=title)
      #plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
      #plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
      #plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
      #plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
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

def plotImage(st,sx,f,fmin=None,fmax=None,pp=True,zoom=False,png=None):
  wpt = 252 # width (in points) of plot
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(st,sx,f)
  if fmax:
    if fmin:
      pv.setClips(fmin,fmax)
    else:
      pv.setClips(-fmax,fmax)
  sp.setFontSizeForPrint(8,wpt)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("PP time (s)")
  if zoom:
    if pp:
      sp.setVLimits(0.352,1.452)
      #sp.setLimits(290,90,410,305)
    else:
      sp.setLimits(290,290,410,560)
  sp.setSize(430,490)
  if pngDir and png:
    sp.paintToPng(720,wpt/72.0,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
