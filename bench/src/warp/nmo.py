#############################################################################
# Demo wavelet estimation from NMO stretch

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import NormalMoveout,WaveletNmo

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  goEstimateWaveletFromGather(args[1])

def goEstimateWaveletFromGather(name):
  """ Estimates wavelet from a gather sampled in time and offset """
  print name
  if name == "syn1": # Synthetic with one reflector
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    nref,vnmo = 1,2.0 # number of reflectors and NMO velocity
    tran,tbed = False,False
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.00
    texp,tbal = 0.00,0
    tmin,tmax,perc = 0.75,1.75,100.0
    zp = True # zero-phase?
    if zp:
      na,ka = 5,-2
      nh,kh = 251,-125
      decay *= 4.0
    else:
      na,ka = 11,0
      nh,kh = 151,-25
    hsyn = getArWavelet(freq,decay,st,nh,kh,zp)
  elif name == "synr": # Synthetic with random reflectors
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    tran,tbed = True,False
    nref,vnmo = 40,2.0 # number of reflectors and NMO velocity
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.00
    texp,tbal = 0.00,0
    tmin,tmax,perc = 0.0,1.75,100.0
    zp = False # zero-phase?
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 151,-25 # sampling for wavelet h
    hsyn = getArWavelet(freq,decay,st,nh,kh)
  elif name == "synt": # Synthetic with random thin beds
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    tran,tbed = True,True
    nref,vnmo = 40,2.0 # number of reflectors and NMO velocity
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.0001
    texp,tbal = 0.00,0
    tmin,tmax,perc = 0.15,1.75,100.0
    zp = True # zero-phase?
    if zp:
      na,ka = 5,-2
      nh,kh = 251,-125
      decay *= 4.0
    else:
      na,ka = 11,0
      nh,kh = 151,-25
    hsyn = getArWavelet(freq,decay,st,nh,kh,zp)
  elif name == "oz01": # Vibroseis
    st = Sampling(1275,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(53,0.100584,-2.615184); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 3.00 # NMO velocity
    fmin,fmax,sfac = 5.0,80.0,1.01
    texp,tbal = 0.00,100
    tmin,tmax,perc = 1.5,2.5,99
    na,ka = 21,0 # sampling for inverse wavelet a
    nh,kh = 51,-25 # sampling for wavelet h
    hsyn = None
  elif name == "oz04": # Vibroseis
    st = Sampling(1275,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(52,0.1,-2.55); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 3.00 # NMO velocity
    fmin,fmax,sfac = 5.0,80.0,1.01
    texp,tbal = 0.00,100
    tmin,tmax,perc = 0.0,5.0,99
    na,ka = 21,0 # sampling for inverse wavelet a
    nh,kh = 51,-25 # sampling for wavelet h
    hsyn = None
  elif name == "oz16": # Airgun
    st = Sampling(1325,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(  48,0.025,0.233); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 1.95 # NMO velocity
    fmin,fmax,sfac = 5.0,50.0,1.001
    texp,tbal = 0.00,100
    tmin,tmax,perc = 0.8,2.3,98
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 201,-50 # sampling for wavelet h
    hsyn = None
  elif name == "oz30": # Airgun
    st = Sampling(2175,0.004,0.00400); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(  96,0.025,0.23075); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 1.55 # NMO velocity
    fmin,fmax,sfac = 5.0,100.0,1.00
    texp,tbal = 0.00,100
    #tmin,tmax,perc = 1.2,2.2,98
    #tmin,tmax,perc = 2.5,3.5,98
    tmin,tmax,perc = 1.2,3.2,98
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 151,-25 # sampling for wavelet h
    hsyn = None
  f = zerofloat(nt,nx)
  if name[0:3]=="syn":
    p = makeCmpReflections(vnmo,nref,st,sx,random=tran,thinBeds=tbed)
    f = addArWavelet(freq,decay,st,sx,p,zp)
  else:
    ais = ArrayInputStream("/data/seis/oz/"+name+".F")
    ais.readFloats(f)
    ais.close()
  f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  vnmo = fillfloat(vnmo,nt)
  nmo = NormalMoveout()
  #nmo.setStretchMax(9.0)
  wn = WaveletNmo(nmo)
  wn.setFrequencyRange(fmin*dt,fmax*dt)
  wn.setTimeRange(int(tmin/dt),int(tmax/dt))
  wn.setStabilityFactor(sfac)
  e = nmo.apply(st,sx,vnmo,f)
  apef = wn.getInverseAPef(na,ka,f)
  anmo = wn.getInverseANmo(na,ka,st,sx,vnmo,f)
  hpef = wn.getWaveletH(na,ka,apef,nh,kh);
  hnmo = wn.getWaveletH(na,ka,anmo,nh,kh);
  plotWavelets(Sampling(nh,st.delta,kh*st.delta),[hnmo,hpef,hsyn],
               title="estimated wavelets")
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  alist = [apef,anmo]
  hlist = [hpef,hnmo]
  tlist = ["PEF","NMO"]
  for ia in range(0,len(alist)):
    a = alist[ia]
    h = hlist[ia]
    t = tlist[ia]
    nah = na+nh
    kah = ka+kh
    ah = zerofloat(nah)
    conv(na,ka,a,nh,kh,h,nah,kah,ah)
    g = wn.applyHNmoA(na,ka,a,nh,kh,h,st,sx,vnmo,f)
    epef = wn.getVariancePef(na,ka,a,f)
    enmo = wn.getVarianceNmo(na,ka,a,st,sx,vnmo,f)
    enor = wn.getNormalizedVarianceNmo(na,ka,a,st,sx,vnmo,f)
    print t+": epef =",epef," enmo =",enmo," enor =",enor
    print " a ="; dump(a)
    plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,
      title=t+": improved NMO")

def stackError(g):
  nt,nx = len(g[0]),len(g)
  s = copy(g[0])
  for ix in range(1,nx):
    s = add(s,g[ix])
  s = div(s,nx)
  g = copy(g)
  for ix in range(nx):
    g[ix] = sub(g[ix],s)
  return g

def rms(g):
  nt,nx = len(g[0]),len(g)
  return sqrt(sum(mul(g,g))/nt/nx)

def tpow(power,st,f):
  """Applies t^power gain."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  tp = pow(rampfloat(ft,dt,nt),power) # sampled times raised to power
  g = zerofloat(nt,nx)
  for ix in range(nx):
    mul(tp,f[ix],g[ix])
  return g

def balance(sigma,f):
  f = add(max(f)*0.00001,f)
  ff = mul(f,f)
  RecursiveExponentialFilter(sigma).apply1(ff,ff)
  return div(f,sqrt(ff))

def normalize(h):
  return div(h,max(max(h),-min(h)))

def makeCmpReflections(vel,nref,st,sx,random=False,thinBeds=False):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  p = zerofloat(nt,nx)
  if random:
    rand = Random(21) #14, 3141, 6, 17
    ts = add(ft,mul((nt-1)*dt,randfloat(rand,nref)))
    rs = sub(mul(2.0,randfloat(rand,nref)),1.0)
  else:
    ts = rampfloat(nt*dt/(nref+1),nt*dt/(nref+1),nref)
    rs = fillfloat(1.0,nref)
  if thinBeds:
    tsc = copy(ts)
    rsc = copy(rs)
    nref *= 2
    ts = zerofloat(nref)
    rs = zerofloat(nref)
    for iref in range(0,nref,2):
      ts[iref] = tsc[iref/2]
      ts[iref+1] = ts[iref]+5.0*dt
      rs[iref] = rsc[iref/2]
      rs[iref+1] = -rs[iref]
  si = SincInterp.fromErrorAndFrequency(0.01,0.45)
  for jx in range(nx):
    xj = sx.getValue(jx)
    cj = (xj*xj)/(vel*vel)
    for jr in range(nref):
      tj = ts[jr]
      rj = rs[jr]
      tj = sqrt(tj*tj+cj)
      si.accumulate(tj,rj,nt,dt,ft,p[jx])
  return p

def getArWavelet(fpeak,decay,st,nh,kh,zp=False):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  if zp:
    gain *= sqrt(1.0+a1*a1+a2*a2)
  x = zerofloat(nh)
  t = zerofloat(nh)
  h = zerofloat(nh)
  x[-kh] = 1.0
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(x,t)
  if zp:
    rcf.applyReverse(t,h)
  else:
    copy(t,h)
  return h

def addArWavelet(fpeak,decay,st,sx,p,zp=False):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  if zp:
    gain *= sqrt(1.0+a1*a1+a2*a2)
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.apply1Forward(p,t)
  if zp:
    rcf.apply1Reverse(t,x)
  else:
    copy(t,x)
  return x

def plotGather(st,sx,p,tmin=None,tmax=None,perc=None,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.setHLabel("Offset (km)")
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  if tmin!=None and tmax!=None:
    sp.setVLimits(tmin,tmax)
  pv = sp.addPixels(st,sx,p)
  if perc:
    pv.setPercentiles(100-perc,perc)

def plotSequence(s,x,xmax=None,title=None):
  sp = SimplePlot.asPoints(s,x)
  if xmax==None:
    xmax = max(abs(max(x)),abs(min(x)))
    xmax *= 1.05
  sp.setVLimits(-xmax,xmax)
  if title:
    sp.setTitle(title)

def plotWavelets(st,hs,hmax=None,title=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(ls[ih])
      pv.setLineWidth(2)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  if title:
    sp.setTitle(title)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
