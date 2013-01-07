import os,sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from wavelet import *

def main(args):
  #demo(makeImpulse)
  #demo(makeSweep)
  demo(readTrace)

def demo(makeInput):
  st,sf,x = makeInput()
  plotTrace(st,x)
  nt = len(x)
  yr = zerofloat(nt)
  yi = zerofloat(nt)
  mt = MorletTransform(st,sf)
  slogf = mt.getLogFrequencySampling()
  y = mt.apply(x)
  ya = mt.abs(y)
  yp = mt.arg(y)
  plotAmplitude(st,slogf,ya)
  plotPhase(st,slogf,yp)
  plotReal(st,slogf,y)
  plotImag(st,slogf,y)

def makeImpulse():
  nt = 251
  x = zerofloat(nt)
  x[nt/2] = 1.0
  st = Sampling(nt,1.0,0.0) # time is in samples
  nf,fmin,fmax = 101,0.025,0.40 # fmin,fmax are in cycles/sample
  sf = MorletTransform.frequencySampling(nf,fmin,fmax)
  return st,sf,x

def makeSweep():
  nt,fmin,fmax = 251,0.025,0.4
  a = fmin
  b = 0.5*(fmax-fmin)/(nt-1)
  t = rampfloat(0.0,1.0,nt)
  p = mul(2.0*PI,add(mul(fmin,t),mul(b,mul(t,t))))
  x = sin(p)
  st = Sampling(nt,1.0,0.0) # time is in samples
  nf,fmin,fmax = 101,fmin,fmax # fmin,fmax are in cycles/sample
  sf = MorletTransform.frequencySampling(nf,fmin,fmax)
  return st,sf,x

def readTrace():
  nt,nx = 251,357
  p = zerofloat(nt,nx)
  ais = ArrayInputStream(getDataDir()+"tp73.dat")
  ais.readFloats(p)
  ais.close()
  st = Sampling(nt,0.004,0.5)
  nf,fmin,fmax = 101,10.0,100.0
  sf = MorletTransform.frequencySampling(nf,fmin,fmax)
  return st,sf,p[nx/2]

def plotTrace(st,x):
  sp = SimplePlot.asSequence(st,x)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude")
  sp.setTitle("Input trace")
  sp.paintToPng(600,3.33,"png/trace.png")

def plotReal(st,slogf,y):
  sp = SimplePlot()
  pv = sp.addPixels(st,slogf,y[0])
  pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  pv.setColorModel(ColorMap.GRAY)
  sp.addColorBar("Amplitude")
  sp.setHLabel("Time (s)")
  sp.setVLabel("log10[ Frequency (Hz) ]")
  sp.setTitle("Morlet transform: real part")
  sp.plotPanel.setColorBarWidthMinimum(100)
  sp.paintToPng(600,3.33,"png/real.png")

def plotImag(st,slogf,y):
  sp = SimplePlot()
  pv = sp.addPixels(st,slogf,y[1])
  pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  pv.setColorModel(ColorMap.GRAY)
  sp.addColorBar("Amplitude")
  sp.setHLabel("Time (s)")
  sp.setVLabel("log10[ Frequency (Hz) ]")
  sp.setTitle("Morlet transform: imaginary part")
  sp.plotPanel.setColorBarWidthMinimum(100)
  sp.paintToPng(600,3.33,"png/imaginary.png")

def plotAmplitude(st,slogf,y):
  sp = SimplePlot()
  pv = sp.addPixels(st,slogf,y)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  pv.setColorModel(ColorMap.JET)
  sp.addColorBar("Amplitude")
  sp.setHLabel("Time (s)")
  sp.setVLabel("log10[ Frequency (Hz) ]")
  sp.setTitle("Morlet transform: amplitude")
  sp.plotPanel.setColorBarWidthMinimum(100)
  sp.paintToPng(600,3.33,"png/amplitude.png")

def plotPhase(st,slogf,y):
  y = mul(360.0,y)
  sp = SimplePlot()
  pv = sp.addPixels(st,slogf,y)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  pv.setColorModel(ColorMap.getHue(0.0,1.0))
  sp.addColorBar("Phase (degrees)")
  sp.setHLabel("Time (s)")
  sp.setVLabel("log10[ Frequency (Hz) ]")
  sp.setTitle("Morlet transform: phase")
  sp.plotPanel.setColorBarWidthMinimum(100)
  sp.paintToPng(600,3.33,"png/phase.png")

#############################################################################
# utilities

def getDataDir():
  scriptDir = sys.path[0]
  baseIndex = scriptDir.find("jtk"+os.sep+"src")
  if baseIndex<0:
    baseIndex = scriptDir.find("idh"+os.sep+"bench")
  if baseIndex<0:
    return None
  dataDir = scriptDir[:baseIndex]+"jtk"+os.sep+"data"+os.sep
  return dataDir

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
