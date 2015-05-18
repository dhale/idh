#############################################################################
# Demo wavelet estimation from warping.
# This version minimizes || f - HSAg || such that a0 = 1.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import WaveletWarpingAH
from warp import WarpingFilter

#############################################################################

#pngDir = "./png/sinow/"
pngDir = None

def main(args):
  #goTest()
  #goSinoImages()
  goSinoTraces()

def goTest():
  nt,ni = 981,500 # number of time samples; number of impulses
  ffreq,fdecay = 0.08,0.05 # peak frequency and decay for wavelet
  gfreq,gdecay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 11,-5 # sampling for inverse wavelet A
  nb,kb = 11,-5 # sampling for inverse wavelet B
  nc,kc = 181,-90 # sampling for wavelet C
  nd,kd = 181,-90 # sampling for wavelet D
  dt,ft = 0.004,0.000 # used for plotting only
  st = Sampling(nt,dt,ft)
  for mp in [False]: # True, for minimum-phase; False for mixed-phase
    ck = getWavelet(ffreq,fdecay,nc,kc,mp) # known wavelet C
    dk = getWavelet(gfreq,gdecay,nc,kc,mp) # known wavelet D
    normalizeMax(ck)
    normalizeMax(dk)
    #for r0,r1 in [(3.0,1.5),(2.0,2.0)]: # <1 for stretch; >1 for squeeze
    for r0,r1 in [(3.0,1.5)]: # <1 for stretch; >1 for squeeze
      title = "r0 = "+str(r0)+"  r1 = "+str(r1)
      u,p,q = logupq(r0,r1,nt,ni)
      f = addWavelet(ffreq,fdecay,p,mp)
      g = addWavelet(gfreq,fdecay,q,mp)
      ww = WaveletWarpingAH()
      ww.setTimeRange(0,nt/2)
      ww.setStabilityFactor(0.0001)
      ww.setMaxIterations(40)
      #bw,cw = ww.getAHSimple(nb,kb,nc,kc,u,f,g);
      bw,cw = ww.getAHNewton(nb,kb,nc,kc,u,f,g);
      #bw,cw = ww.getAHGaussNewton(nb,kb,nc,kc,u,f,g);
      print "rmse =",ww.rms(sub(f,ww.applyHSA(nb,kb,bw,nc,kc,cw,u,g)))
      dump(bw);
      dw = ww.getWavelet(nd,kd,nb,kb,bw)
      sg = ww.warp(u,g)
      bg = ww.convolve(nb,kb,bw,g)
      sbg = ww.warp(u,bg)
      csbg = ww.convolve(nc,kc,cw,sbg)
      normalizeMax(cw)
      normalizeMax(dw)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,csbg],labels=["f","CSBg"],title=title)
      plotWavelets(Sampling(nc,dt,kc*dt),[cw,ck],title=title+" C")
      plotWavelets(Sampling(nd,dt,kd*dt),[dw,dk],title=title+" D")

def goSinoTraces():
  na,ka = 41,-20 # sampling for inverse A
  nh,kh = 41,-20 # sampling for wavelet H
  nt,dt,ft = 501,0.004,0.000 # used for plotting only
  #nx,dx,fx = 721,0.015,0.000
  sa = Sampling(na,dt,ka*dt)
  sh = Sampling(nh,dt,kh*dt)
  st = Sampling(nt,dt,ft)
  f,g,u = getSinoTraces(16) # PP trace f(t), PS trace g(u), and warping u(t)
  ww = WaveletWarpingAH()
  ww.setTimeRange(100,500) # PP time window
  ww.setStabilityFactor(0.001) # stability factor
  ww.setMaxIterations(100) # max Gauss-Newton iterations
  a,h = ww.getAHGaussNewton(na,ka,nh,kh,u,f,g)
  sg = ww.warp(u,g)
  ag = ww.convolve(na,ka,a,g)
  sag = ww.warp(u,ag)
  hsag = ww.convolve(nh,kh,h,sag)
  dump(a)
  g = copy(nt,g)
  sg = copy(nt,sg)
  sag = copy(nt,sag)
  hsag = copy(nt,hsag)
  fmhsag= sub(f,hsag)
  plotWavelets(sa,[a],title="A")
  plotWavelets(sh,[h],title="H")
  plotTraces(st,[f,hsag,fmhsag],labels=["f","HSAg","f-HSAg"])
def getSinoTraces(m2=1):
  dataDir = "/data/seis/sino/warp/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  u = readImage(dataDir+"shifts.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  i2 = 300 # first trace in stack
  fs = copy(f[i2])
  gs = copy(g[i2])
  us = copy(u[i2])
  for j2 in range(1,m2):
    add(f[i2+j2],fs,fs)
    add(g[i2+j2],gs,gs)
    add(u[i2+j2],us,us)
  mul(1.0/m2,fs,fs)
  mul(1.0/m2,gs,gs)
  mul(1.0/m2,us,us)
  gain(100,fs)
  gain(100,gs)
  return fs,gs,us

def goSinoImages():
  na,ka = 21,-10 # sampling for inverse A
  nh,kh = 81,-40 # sampling for wavelet H
  nt,dt,ft = 501,0.004,0.000 # used for plotting only
  #nx,dx,fx = 721,0.015,0.000
  nx,dx,fx = 101,0.015,0.000
  sa = Sampling(na,dt,ka*dt)
  sh = Sampling(nh,dt,kh*dt)
  st = Sampling(nt,dt,ft)
  sx = Sampling(nx,dx,fx)
  tmin,tmax = 100,500 # PP time window
  f,g,u = getSinoImages() # PP image, PS image, and warping u(t,x)
  ww = WaveletWarpingAH()
  ww.setTimeRange(100,400) # PP time window
  ww.setStabilityFactor(0.0001) # stability factor
  ww.setMaxIterations(20) # max Gauss-Newton iterations
  a,h = ww.getAHGaussNewton(na,ka,nh,kh,u,f,g)
  dump(a)
  sg = ww.warp(u,g)
  ag = ww.convolve(na,ka,a,g)
  sag = ww.warp(u,ag)
  hsag = ww.convolve(nh,kh,h,sag)
  normalizeRms(f)
  normalizeRms(sg)
  normalizeRms(hsag)
  plotImage(st,sx,f,zoom=True,png="f")
  plotWavelets(sa,[a],title="A")
  plotWavelets(sh,[h],title="H")
  plotImage(st,sx,sg,zoom=True,png="Sg")
  plotImage(st,sx,hsag,zoom=True,png="HSAg")
def getSinoImages():
  dataDir = "/data/seis/sino/warp/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  u = readImage(dataDir+"shifts.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  f = copy(n1f,101,0,300,f)
  g = copy(n1g,101,0,300,g)
  u = copy(n1f,101,0,300,u)
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
  maxf = abs(max(f))
  minf = abs(min(f))
  if minf<maxf:
    scale = 1.0/maxf
  else:
    scale = -1.0/minf
  return mul(scale,f,f)

def normalizePow(f):
  scale = 1.0/sqrt(sum(mul(f,f)))
  mul(scale,f,f)

def normalizeRms(f):
  mul(1.0/rms(f),f,f)
def normalizeRms1(f):
  mul(1.0/rms1(f),f,f)
def rms(f):
  return sqrt(sum(mul(f,f))/len(f)/len(f[0]))
def rms1(f):
  return sqrt(sum(mul(f,f))/len(f))

def logupq(r0,r1,nt,ni):
  rand = Random(3)
  umax = nt-1
  us = mul(umax,randfloat(rand,ni));
  rs = pow(mul(sub(randfloat(rand,ni),0.5),2.0),5.0)
  if r0==r1:
    u = rampfloat(0.0,r0,nt)
    def t(u):
      return u/r0
  else:
    a = umax/log(r0/r1)
    b = r0/a
    u = mul(a,log(rampfloat(1.0,b,nt)))
    def t(u):
      return (exp(u/a)-1.0)/b
  p = zerofloat(nt)
  q = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
  for ji in range(ni):
    ui = us[ji]
    ri = rs[ji]
    ti = t(ui)
    si.accumulate(ti,ri,nt,1.0,0.0,p)
    si.accumulate(ui,ri,nt,1.0,0.0,q)
  return u,p,q

def logupqOld(r0,r1,nt,ni):
  umax = nt-1
  us = rampfloat(umax/(ni+1),umax/(ni+1),ni)
  if r0==r1:
    u = rampfloat(0.0,r0,nt)
    def t(u):
      return u/r0
  else:
    a = umax/log(r0/r1)
    b = r0/a
    u = mul(a,log(rampfloat(1.0,b,nt)))
    def t(u):
      return (exp(u/a)-1.0)/b
  p = zerofloat(nt)
  q = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
  for ji in range(ni):
    ui = us[ji]
    ti = t(ui)
    ri = 1.0
    si.accumulate(ti,ri,nt,1.0,0.0,p)
    si.accumulate(ui,ri,nt,1.0,0.0,q)
  return u,p,q

def expupq(r0,r1,nt,ni):
  tmax = nt-1
  b = log(r1/r0)/tmax
  a = r0/b
  u = sub(tmax,mul(a,sub(exp(rampfloat(0.0,b,nt)),1.0)))
  ts = rampfloat(tmax/(ni+1),tmax/(ni+1),ni)
  p = zerofloat(nt)
  q = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
  rj = -1.0
  for ji in range(ni):
    tj = ts[ji]
    tp = ts[ji]
    tq = tmax-a*(exp(b*tp)-1.0)
    #print "tp =",tp," tq =",tq
    rj = -rj
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
  return u,p,q

def nmoupq(r0,r1,nt,ni):
  tmax = nt-1
  a = tmax*tmax*(1.0/(r1*r1)-1.0)
  u = sqrt(add(pow(rampfloat(0,1.0,nt),2.0),a))
  ts = rampfloat(tmax/(ni+1),tmax/(ni+1),ni)
  p = zerofloat(nt)
  q = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
  rj = -1.0
  for ji in range(ni):
    tj = ts[ji]
    tp = ts[ji]
    tq = sqrt(tp*tp+a)
    #print "tp =",tp," tq =",tq
    rj = -rj
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
  return u,p,q

def makeImpulses(r,nt,ni):
  p = zerofloat(nt)
  q = zerofloat(nt)
  tmax = nt-1
  if r<=1.0:
    ts = rampfloat(tmax/(ni+1),tmax/(ni+1),ni)
  else:
    ts = rampfloat(tmax/(ni+1)/r,tmax/(ni+1)/r,ni)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
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

def plotTraces(st,xs,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(ix,0,st,xi)
    pp.setVLimits(ix,-3.0,3.0)
    if labels:
      pp.setVLabel(ix,labels[ix])
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(1400,700)
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
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(ls[ih])
      pv.setLineWidth(2)
      #sv = sp.addSequence(st,hs[ih])
      #sv.setColor(cs[ih])
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
  #sp.setFontSizeForPrint(8,wpt)
  #sp.setSize(650,850)
  sp.setSize(200,850)
  pv = sp.addPixels(st,sx,f)
  pv.setClips(-2.0,2.0)
  if zoom:
    #sp.setVLimits(0.400,1.600) # zoom consistent with tmin,tmax = 100,400
    sp.setVLimits(0.400,2.000) # zoom consistent with tmin,tmax = 100,500
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
  pf.setSize(650,850)
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
