#############################################################################
# Demo nmo stretch

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import WaveletNmo,WarpedWavelet,ShapingFilter

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  #goEstimateWaveletFromGather(args[1])
  goEstimateWaveletFromVigGather()

def goEstimateWaveletFromVigGather():
  """ Estimates wavelet from Viking Graben gather """
  st = Sampling(1500,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
  #sx = Sampling(60,0.01524,-.97902); nx,dx,fx = sx.count,sx.delta,sx.first
  sx = Sampling(60,0.050,-3.212); nx,dx,fx = sx.count,sx.delta,sx.first
  ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
  vs=1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735
  #vnmo = CubicInterpolator(ts,vs).interpolate(rampfloat(ft,dt,nt))
  vnmo = 1.550
  smax = 100
  fmin,fmax,sfac = 0.0,50.0,1.00
  texp,tbal = 2.0,100
  tmin,tmax,perc = 1.0,2.0,98.0
  na,ka = 11,0 # sampling for inverse wavelet a
  nh,kh = 301,-50 # sampling for wavelet h
  f = zerofloat(nt,nx)
  ais = ArrayInputStream("/data/seis/vig/scdp800to809.dat")
  ais.skipBytes(4*1500*60*0)
  ais.readFloats(f)
  ais.close()
  f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)
  wn = WaveletNmo(st,sx,smax,vnmo)
  wn.setFrequencyRange(fmin,fmax)
  wn.setStabilityFactor(sfac)
  e = wn.applyNmo(f)
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  ad = wn.getInverseAPef(na,ka,f)
  ai = wn.getInverseA(na,ka,f)
  ta = 0
  for a in [ad,ai]:
    ta += 1
    h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
    nah = na+nh
    kah = ka+kh
    ah = zerofloat(nah)
    conv(na,ka,a,nh,kh,h,nah,kah,ah)
    g = wn.applyHNmoA(na,ka,a,nh,kh,h,f)
    #g = wn.applyBNmoA(na,ka,a,f)
    epef = wn.getVariancePef(na,ka,a,f)
    evar = wn.getVariance(na,ka,a,f)
    enor = wn.getNormalizedVariance(na,ka,a,f)
    print str(ta)+": epef =",epef," evar =",evar," enor =",enor
    print " a ="; dump(a)
    plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,
      title=str(ta)+": improved NMO")
    #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,
    #  title=str(ta)+": stack error")
    plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),
                 title=str(ta)+": estimated wavelet")
    #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),
    #             title="inverse wavelet")
    #plotSequence(Sampling(nah,st.delta,kah*st.delta),normalize(ah),
    #             title="unit impulse")
  #d = wn.getDifferenceGathers(na,ka,f)
  #for ia in range(na):
  #  plotGather(st,sx,d[ia],
  #             tmin=tmin,tmax=tmax,perc=perc,title="lag="+str(ka+ia))

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
    na,ka = 3,0 # sampling for inverse wavelet a
    nh,kh = 301,-50 # sampling for wavelet h
  elif name == "synr": # Synthetic with random reflectors
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    tran,tbed = True,False
    nref,vnmo = 40,2.0 # number of reflectors and NMO velocity
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.00
    texp,tbal = 0.00,0
    tmin,tmax,perc = 0.0,1.75,100.0
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 301,-50 # sampling for wavelet h
  elif name == "synt": # Synthetic with random thin beds
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    tran,tbed = True,True
    nref,vnmo = 40,2.0 # number of reflectors and NMO velocity
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.00
    texp,tbal = 0.00,0
    tmin,tmax,perc = 0.0,1.75,100.0
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 301,-50 # sampling for wavelet h
  elif name == "oz01": # Vibroseis
    st = Sampling(1275,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(53,0.100584,-2.615184); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 3.00 # NMO velocity
    fmin,fmax,sfac = 5.0,80.0,1.01
    texp,tbal = 0.00,100
    tmin,tmax,perc = 1.5,2.5,99
    na,ka = 21,0 # sampling for inverse wavelet a
    nh,kh = 51,-25 # sampling for wavelet h
  elif name == "oz04": # Vibroseis
    st = Sampling(1275,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(52,0.1,-2.55); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 3.00 # NMO velocity
    fmin,fmax,sfac = 5.0,80.0,1.01
    texp,tbal = 0.00,100
    tmin,tmax,perc = 0.0,5.0,99
    na,ka = 21,0 # sampling for inverse wavelet a
    nh,kh = 51,-25 # sampling for wavelet h
  elif name == "oz16": # Airgun
    st = Sampling(1325,0.004,0.004); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(  48,0.025,0.233); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 1.95 # NMO velocity
    fmin,fmax,sfac = 5.0,50.0,1.001
    texp,tbal = 0.00,100
    tmin,tmax,perc = 0.8,2.3,98
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 201,-50 # sampling for wavelet h
  elif name == "oz30": # Airgun
    st = Sampling(2175,0.004,0.00400); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(  96,0.025,0.23075); nx,dx,fx = sx.count,sx.delta,sx.first
    vnmo = 1.55 # NMO velocity
    fmin,fmax,sfac = 5.0,50.0,1.00
    texp,tbal = 0.00,100
    #tmin,tmax,perc = 2.5,3.0,98
    tmin,tmax,perc = 1.2,3.2,98
    #tmin,tmax,perc = 1.2,4.2,99.5
    na,ka = 11,0 # sampling for inverse wavelet a
    nh,kh = 201,-50 # sampling for wavelet h
  f = zerofloat(nt,nx)
  if name[0:3]=="syn":
    p = makeCmpReflections(vnmo,nref,st,sx,random=tran,thinBeds=tbed)
    f = addArWavelet(freq,decay,st,sx,p)
  else:
    ais = ArrayInputStream("/data/seis/oz/"+name+".F")
    ais.readFloats(f)
    ais.close()
  f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)
  smax = 2.0
  wn = WaveletNmo(st,sx,smax,vnmo)
  wn.setFrequencyRange(fmin,fmax)
  wn.setStabilityFactor(sfac)
  e = wn.applyNmo(f)
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  ad = wn.getInverseAPef(na,ka,f)
  ai = wn.getInverseA(na,ka,f)
  ta = 0
  for a in [ad,ai]:
    ta += 1
    h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
    #zero(h); h[-kh] = 1.0
    nah = na+nh
    kah = ka+kh
    ah = zerofloat(nah)
    conv(na,ka,a,nh,kh,h,nah,kah,ah)
    #g = wn.applyHNmoA(na,ka,a,nh,kh,h,f)
    g = wn.applyBNmoA(na,ka,a,f)
    epef = wn.getVariancePef(na,ka,a,f)
    evar = wn.getVariance(na,ka,a,f)
    enor = wn.getNormalizedVariance(na,ka,a,f)
    print str(ta)+": epef =",epef," evar =",evar," enor =",enor
    print " a ="; dump(a)
    plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,
      title=str(ta)+": improved NMO")
    #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,
    #  title=str(ta)+": stack error")
    plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),
                 title=str(ta)+": estimated wavelet")
    #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),
    #             title="inverse wavelet")
    #plotSequence(Sampling(nah,st.delta,kah*st.delta),normalize(ah),
    #             title="unit impulse")
  """
  d = wn.getDifferenceGathers(na,ka,f)
  for ia in range(na):
    plotGather(st,sx,d[ia],
               tmin=tmin,tmax=tmax,perc=perc,title="lag="+str(ka+ia))
  """

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

def addArWavelet(fpeak,decay,st,sx,p):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.apply1Forward(p,x)
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

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
