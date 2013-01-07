#############################################################################
# Demo seismic data processing with Oz Yilmaz's common-source gathers.

import sys
from java.awt import *
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

# Directories for binary data files and png images
dataDir = "/data/seis/oz/"
pngDir = None #"./png"

# Maps data file name to metadata
metaMap = {
  "oz01":{
    "area":"South Texas",
    "source":"Vibroseis",
    "st":Sampling(1275,0.004,0.004),
    "sx":Sampling(53,0.100584,-2.615184)},
  "oz02":{
    "area":"West Texas",
    "source":"Vibroseis",
    "st":Sampling(1025,0.004,0.004),
    "sx":Sampling(127,0.030480,-2.834640)},
  "oz03":{
    "area":"Louisiana",
    "source":"Dynamite",
    "st":Sampling(1500,0.004,0.004),
    "sx":Sampling(24,0.103632,0.103632)},
  "oz04":{
    "area":"Turkey",
    "source":"Vibroseis",
    "st":Sampling(1275,0.004,0.004),
    "sx":Sampling(52,0.1,-2.55)},
  "oz05":{
    "area":"South America",
    "source":"Dynamite",
    "st":Sampling(3000,0.002,0.002),
    "sx":Sampling(51,0.1,-2.5)},
}

#############################################################################
# Data processing

def main(args):
  #for name in ["oz01"]: #metaMap:
  for name in metaMap:
    process(name)
  return

def process(name):
  st = metaMap[name]["st"]
  sx = metaMap[name]["sx"]
  source = metaMap[name]["source"]
  title = name+": "+source
  f = read(name)
  f = tpow(2.0,st,f)
  plot(name,st,sx,f,name)
  #f = nmo(3.0,st,sx,f)
  #plot(name+": before",st,sx,f,name)
  #f = align(f)
  #plot(name+": after",st,sx,f,name)
  #g,t = alignWithTimes(st,f)
  #plot(name+": before",st,sx,f,name)
  #plot(name+": after",st,sx,g,name)

def read(name):
  """Reads gather with specified name."""
  nt = metaMap[name]["st"].count
  nx = metaMap[name]["sx"].count
  fileName = dataDir+name+".F" # suffix F implies floats
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  f = zerofloat(nt,nx)
  ais.readFloats(f)
  ais.close()
  return f

def tpow(power,st,f):
  """Applies t^power gain."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  tp = pow(rampfloat(ft,dt,nt),power) # sampled times raised to power
  g = zerofloat(nt,nx)
  for ix in range(nx):
    mul(tp,f[ix],g[ix])
  return g

def gpow(power,f):
  """Applies sgn(f)*|f|^power gain."""
  return mul(sgn(f),pow(abs(f),power))

def nmo(v,st,sx,f):
  """Applies NMO correction for specified velocity v."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  t = rampfloat(ft,dt,nt) # sampled times
  tt = mul(t,t) # sampled times squared
  si = SincInterp() # high-fidelity sinc interpolation
  g = zerofloat(nt,nx) # the output gather
  for ix in range(nx): # loop over all traces in gather
    x = sx.getValue(ix) # source-receiver offset x
    t = sqrt(add(tt,(x*x)/(v*v))) # sqrt(t*t+(x*x)/(v*v))
    si.interpolate(nt,dt,ft,f[ix],nt,t,g[ix]) # interpolated trace
  return g

def stack(f):
  nt,nx = len(f[0]),len(f)
  g = zerofloat(nt)
  for ix in range(nx):
    add(f[ix],g,g)
  return g

def align(f):
  """Scary alignment warps traces to flatten reflections."""
  nt,nx = len(f[0]),len(f)
  s = stack(f)
  g = zerofloat(nt,nx)
  du = zerofloat(nt)
  ut = zerofloat(nt)
  sigma = 20 # Gaussian window half-width sigma = 20 samples
  lsf = LocalShiftFinder(sigma)
  for ix in range(nx):
    lsf.find1(-10,10,s,f[ix],du)
    copy(f[ix],g[ix])
    lsf.shift1(du,ut,g[ix])
  return g

def alignWithTimes(st,f):
  """Scary alignment warps traces to flatten reflections."""
  nt,nx = len(f[0]),len(f)
  dt,ft = st.delta,st.first
  s = stack(f)
  g = zerofloat(nt,nx)
  t = zerofloat(nt,nx)
  tau = rampfloat(ft,dt,nt)
  du = zerofloat(nt)
  ut = zerofloat(nt)
  sigma = 20 # Gaussian window half-width sigma = 20 samples
  lsf = LocalShiftFinder(sigma)
  for ix in range(nx):
    lsf.find1(-20,20,s,f[ix],du)
    copy(f[ix],g[ix])
    lsf.shift1(du,ut,g[ix])
    mul(dt,du,du)
    sub(tau,du,t[ix])
  return g,t

def plot(title,st,sx,f,png=None):
  f = gpow(0.5,f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,900)
  sp.setVLabel("Time (s)")
  sp.setHLabel("Offset (km)")
  sp.setTitle(title)
  pv = sp.addPixels(st,sx,f)
  pv.setPercentiles(1.0,99.0)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if png and pngDir:
    sp.paintToPng(100,6,pngDir+"/"+png+".png")

#############################################################################
# Run everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
