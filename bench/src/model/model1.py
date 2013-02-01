#############################################################################
# Synthetic seismograms for vertically propagating plane in 1D model

import sys
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from model import *

#############################################################################

def main(args):
  goDemo2()

def goDemo1():
  sm = SeismicModel1D()
  sm.setSourceType(SeismicModel1D.SourceType.LAND_VIBROSEIS)
  sm.setSensorType(SeismicModel1D.SensorType.GEOPHONE)
  sm.setSurfaceReflectionCoefficient(0.0)
  sm.addLayer(0.0,2.0,2.0,1.0e6)
  sm.addLayer(1.0,6.0,2.0,1.0e6)
  sm.addSource(0.0,1.0)
  sm.addSensor(0.0)
  fpeak = 25.0
  sm.setRickerWavelet(fpeak)
  #sm.dumpLayers()
  nt = 1101
  dt = 0.004
  fref = 0.5/dt
  st = Sampling(nt,dt,0.0)
  s = sm.makeSeismograms(nt,dt,fref)[0]
  subtractRickerWavelet(nt,dt,fpeak,s)
  print "min =",min(s)," max =",max(s)
  SimplePlot.asPoints(st,s)

def goDemo2():
  nt,dt = 1251,0.004
  nz,dz = 10000,0.00039
  v = add(1.987,mul(0.01,randfloat(nz)))
  p = add(2.123,mul(0.01,randfloat(nz)))
  fpeak = 25.0
  st = Sampling(nt,dt,0.0)
  s1 = model1Seismogram(nz,dz,v,p,nt,dt,fpeak); print "s1 done"
  s2 = simpleSeismogram(nz,dz,v,p,nt,dt,fpeak); print "s2 done"
  ds =  sub(s1,s2)
  smin = min(min(s1),min(s2))
  smax = max(max(s1),max(s2))
  plot(st,s1,smin,smax)
  plot(st,s2,smin,smax)
  plot(st,ds,smin,smax)

def plot(st,s,smin=None,smax=None):
  sp = SimplePlot()
  pv = sp.addPoints(st,s)
  if smin and smax:
    sp.setVLimits(smin,smax)
  sp.setSize(1200,500)

def model1Seismogram(nz,dz,v,p,nt,dt,fpeak):
  sm = SeismicModel1D()
  for iz in range(nz):
    sm.addLayer(iz*dz,v[iz],p[iz],1.0e6)
  sm.setSurfaceReflectionCoefficient(0.0)
  sm.setSourceType(SeismicModel1D.SourceType.LAND_VIBROSEIS)
  sm.setSensorType(SeismicModel1D.SensorType.GEOPHONE)
  sm.addSource(0.0,1.0)
  sm.addSensor(0.0)
  sm.setDecay(0.1)
  sm.setOversample(2)
  sm.setRickerWavelet(fpeak)
  #sm.dumpLayers()
  s = sm.makeSeismograms(nt,dt,0.5/dt)[0]
  subtractRickerWavelet(nt,dt,fpeak,s)
  return s

def subtractRickerWavelet(nt,dt,fpeak,s):
  w = getRickerWavelet(fpeak,dt)
  hw = (len(w)-1)/2
  for it in range(min(nt,hw)):
    s[it] -= w[hw+it] # subtract direct arrival at time zero

def simpleSeismogram(nz,dz,v,p,nt,dt,fpeak):
  scale = PI*fpeak*dt
  h = int(10.0/scale)
  s = zerofloat(nt)
  zp = v[0]*p[0]
  t = 0.0
  for iz in range(nz-1):
    t += 2.0*dz/v[iz]
    zm = zp
    zp = v[iz+1]*p[iz+1]
    r = (zm-zp)/(zm+zp)
    itlo = max(0,int(t/dt-h))
    ithi = min(nt-1,int(t/dt+h))
    for it in range(itlo,ithi):
      ti = it*dt
      s[it] += r*ricker(fpeak,ti-t)
  return s

def makeSequences():
  n = 500
  fpeak = 0.125
  shift = 2.0/fpeak
  #w = Warp1Function.constant(shift,n)
  w = WarpFunction1.sinusoid(shift,n)
  #f = makeCosine(fpeak,n)
  f = makeRandomEvents(n,seed=seed); 
  g = w.warp(f)
  f = addRickerWavelet(fpeak,f)
  g = addRickerWavelet(fpeak,g)
  s = zerofloat(n)
  for i in range(n):
    s[i] = w.ux(i)
  return f,g,s

def makeCosine(freq,n):
  return cos(mul(2.0*PI*freq,rampfloat(0.0,1.0,n)))

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def convolveWithRickerWavelet(fpeak,dt,f):
  w = getRickerWavelet(fpeak,dt)
  nw = len(w)
  kw = -(nw-1)/2
  nt = len(f)
  g = zerofloat(nt)
  Conv.conv(nw,kw,w,nt,0,f,nt,0,g)
  return g

def getRickerWavelet(fpeak,dt):
  scale = PI*fpeak*dt
  i0 = int(10.0/scale)
  nt = 1+2*i0
  w = zerofloat(nt)
  for it in range(nt):
    x = scale*(it-i0)
    w[it] = (1.0-2.0*x*x)*exp(-x*x)
  return w
def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
