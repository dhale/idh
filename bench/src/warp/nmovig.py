#############################################################################
# Demo wavelet estimation via NMO for Viking Graben CMP gathers

from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import NormalMoveout,WaveletNmo

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  #goNmoVigGather()
  goEstimateWaveletFromVigGather()

def goNmoVigGather():
  """ Applies NMO correction to Viking Graben gather """
  st,sx,f = readVigGather(800)
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
  vs=1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735
  vnmoP = CubicInterpolator(ts,vs).interpolate(rampfloat(ft,dt,nt))
  vnmoM = fillfloat(1.475,nt)
  vnmo = vnmoM
  texp,tbal,smax = 2.0,100,9.00
  tmin,tmax,perc = 0.4,2.4,98.0
  #f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)
  nmo = NormalMoveout()
  nmo.setStretchMax(smax)
  g = nmo.apply(st,sx,vnmo,f)
  h = nmo.stackAndReplicate(g)
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,perc=perc,title="input");
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="output");
  plotGather(st,sx,h,tmin=tmin,tmax=tmax,perc=perc,title="stack");

def goEstimateWaveletFromVigGather():
  """ Estimates wavelet from Viking Graben gather """
  st,sx,f = readVigGather(800)
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
  vs=1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735
  vnmoP = CubicInterpolator(ts,vs).interpolate(rampfloat(ft,dt,nt))
  vnmoM = fillfloat(1.500,nt) # water velocity, for multiples
  vnmo = vnmoM
  na,ka = 16,0 # sampling for inverse wavelet a
  nh,kh = 151,-25 # sampling for wavelet h
  fmin,fmax,sfac = 5.0,75.0,1.00
  texp,tbal,smax = 2.0,100,9.00
  tmin,tmax,perc = 0.4,1.9,98.0
  #f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)
  nmo = NormalMoveout()
  nmo.setStretchMax(smax)
  wn = WaveletNmo(nmo)
  wn.setTimeRange(int(tmin/dt),int(tmax/dt))
  wn.setFrequencyRange(fmin*dt,fmax*dt)
  wn.setStabilityFactor(sfac)
  e = nmo.apply(st,sx,vnmo,f)
  apef = wn.getInverseAPef(na,ka,f)
  anmo = wn.getInverseANmo(na,ka,st,sx,vnmo,f)
  hpef = wn.getWaveletH(na,ka,apef,nh,kh);
  hnmo = wn.getWaveletH(na,ka,anmo,nh,kh);
  plotWavelets(Sampling(nh,st.delta,kh*st.delta),[hnmo,hpef],
               title="estimated wavelets")
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,perc=perc,title="CMP gather")
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
    print tlist[ia]+": epef =",epef," enmo =",enmo," enor =",enor
    print " a ="; dump(a)
    plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,
      title=t+": improved NMO")
    #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,
    #  title=t+": stack error")
    #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),
    #             title="inverse wavelet")
    #plotSequence(Sampling(nah,st.delta,kah*st.delta),normalize(ah),
    #             title="unit impulse")
  #d = wn.getDifferenceGathers(na,ka,st,sx,vnmo,f)
  #for ia in range(na):
  #  plotGather(st,sx,d[ia],
  #             tmin=tmin,tmax=tmax,perc=perc,title="lag="+str(ka+ia))

def readVigGather(icdp):
  st = Sampling(1500,0.004,0.004)
  if icdp%2==0:
    sx = Sampling(60,0.050,0.262) # sampling for even cdps
  else:
    sx = Sampling(60,0.050,0.287) # sampling for odd cdps
  nt,dt,ft = st.count,st.delta,st.first
  nx,dx,fx = sx.count,sx.delta,sx.first
  f = zerofloat(nt,nx)
  ais = ArrayInputStream("/data/seis/vig/scdp800to809.dat")
  ais.skipBytes(4*nt*nx*(icdp-800))
  ais.readFloats(f)
  ais.close()
  for ix in range(nx/2):
    fi = f[ix]
    f[ix] = f[nx-1-ix]
    f[nx-1-ix] = fi
  return st,sx,f

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
  for ih in range(nh):
    pv = sp.addPoints(st,hs[ih])
    pv.setLineStyle(ls[ih])
    pv.setLineWidth(2)
  if hmax==None:
    hmax = max(abs(max(hs)),abs(min(hs)))
    hmax *= 1.05
  sp.setVLimits(-hmax,hmax)
  if title:
    sp.setTitle(title)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
