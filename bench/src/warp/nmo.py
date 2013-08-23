#############################################################################
# Demo nmo stretch

from imports import *

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  go()

def go():
  n = 501
  st = Sampling(501,0.004,0.0)
  sx = Sampling(201,0.010,0.0)
  vel = 2.0
  nr = 50
  p = makeCmpReflections(vel,nr,st,sx)
  hp = addArWavelet(30,0.1,st,sx,p)
  shp = applyNmo(vel,st,sx,hp)
  sp = applyNmo(vel,st,sx,p)
  hsp = addArWavelet(30,0.1,st,sx,sp)
  plotCmpGather(st,sx,p)
  plotCmpGather(st,sx,hp)
  plotCmpGather(st,sx,shp)
  plotCmpGather(st,sx,hsp)

def applyNmo(vel,st,sx,p):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  q = copy(p)
  si = SincInterp.fromErrorAndFrequency(0.01,0.4)
  for jx in range(nx):
    xj = sx.getValue(jx)
    cj = (xj*xj)/(vel*vel)
    ti = sqrt(add(pow(rampfloat(ft,dt,nt),2.0),cj))
    si.interpolate(nt,dt,ft,p[jx],nt,ti,q[jx])
  return q

def plotCmpGather(st,sx,p):
  sp = SimplePlot.asPixels(st,sx,p)
  sp.setSize(400,750)

def makeCmpReflections(vel,nr,st,sx):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  p = zerofloat(nt,nx)
  ts = add(ft,mul((nt-1)*dt,randfloat(nr)))
  rs = sub(mul(2.0,randfloat(nr)),1.0)
  si = SincInterp.fromErrorAndFrequency(0.01,0.45)
  for jx in range(nx):
    xj = sx.getValue(jx)
    cj = (xj*xj)/(vel*vel)
    for jr in range(nr):
      tj = ts[jr]
      rj = rs[jr]
      tj = sqrt(tj*tj+cj)
      si.accumulate(tj,rj,nt,dt,ft,p[jx])
  return p

def addArWavelet(fpeak,decay,st,sx,p):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.apply1Forward(p,x)
  return x

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),9.0)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
