#############################################################################
# Demonstrates smoothing with a recursive two-sided exponential filter.

import os,sys
from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

#############################################################################
izv = RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE
izs = RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE

pngDir = None
pngDir = "./png/"

def main(args):
  #goImpulseResponses()
  goSmooth()
  goNormalize()

def goImpulseResponses():
  n,d,f = 41,1.0,-20.0
  m = 10
  sigma = sqrt(m*(m+1)/3.0)
  print "m =",m," a =",sigmaToA(sigma)," sigma =",sigma
  s = Sampling(n,d,f)
  x = zerofloat(n)
  e = zerofloat(n)
  g = zerofloat(n)
  r = zerofloat(n)
  x[(n-1)/2] = 1.0
  ref = RecursiveExponentialFilter(sigma)
  rgf = RecursiveGaussianFilter(sigma)
  rrf = RecursiveRectangleFilter(-m,m)
  ref.apply(x,e)
  rgf.apply0(x,g)
  rrf.apply(x,r)
  sp = SimplePlot()
  svr = sp.addSequence(s,r)
  svr.setColor(Color.RED)
  svg = sp.addSequence(s,g)
  svg.setColor(Color.BLUE)
  sve = sp.addSequence(s,e)
  sve.setColor(Color.BLACK)
  #sp.setTitle("Impulse responses of smoothing filters")
  #sp.setFontSizeForPrint(12,504)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(800,600)
  if pngDir:
    sp.paintToPng(300,7.0,pngDir+"smoothimp.png")

def goNormalize():
  """ 2500 samples, sampling interval = 2 ms, from Gary Ohlhoeft """
  """ y = x-<x>, z = y/sqrt(<y y>) """
  st,x = readSequence("SpringDashpot.txt")
  sigma = 100.0
  ref = RecursiveExponentialFilter(sigma)
  print "sigma =",sigma," a =",sigmaToA(sigma)
  ref.setEdges(izs)
  y = zerofloat(st.count)
  ref.apply(x,y)
  y = sub(y,x)
  yy = mul(y,y)
  ref.apply(yy,yy)
  z = sqrt(yy)
  v = div(y,z)
  sp = SimplePlot()
  pv = sp.addPoints(st,x)
  pv.setLineColor(Color.RED)
  pv = sp.addPoints(st,y)
  sp.setHLimits(-0.05,5.05)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Displacement (cm)")
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(800,500)
  if pngDir:
    sp.paintToPng(300,7.0,pngDir+"debias.png")
  sp = SimplePlot()
  pv = sp.addPoints(st,y)
  pv = sp.addPoints(st,z)
  pv.setLineColor(Color.RED)
  pv.setLineWidth(3)
  sp.setHLimits(-0.05,5.05)
  sp.setVLimits(-0.21,0.21)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Displacement (cm)")
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(800,500)
  if pngDir:
    sp.paintToPng(300,7.0,pngDir+"amplitude.png")
  sp = SimplePlot()
  pv = sp.addPoints(st,v)
  sp.setHLimits(-0.05,5.05)
  sp.setVLimits(-3.10,3.10)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Displacement (normalized)")
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(800,500)
  if pngDir:
    sp.paintToPng(300,7.0,pngDir+"normalize.png")

def goSmooth():
  """ 2500 samples, sampling interval = 2 ms, from Gary Ohlhoeft """
  st,x = readSequence("SpringDashpot.txt")
  sp = SimplePlot()
  pvx = sp.addPoints(st,x)
  pvx.setLineColor(Color.RED)
  info = {
    izv:("smoothizv","Zero values assumed off the ends"),
    izs:("smoothizs","Smoothing noisy data"),
  }
  sp.setHLimits(-0.05,5.05)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Displacement (cm)")
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(800,500)
  if pngDir:
    sp.paintToPng(300,7.0,pngDir+"smoothraw.png")
  sigma = 20.0
  ref = RecursiveExponentialFilter(sigma)
  print "sigma =",sigma," a =",sigmaToA(sigma)
  for edges in info:
    ref.setEdges(edges)
    y = zerofloat(st.count)
    ref.apply(x,y)
    sp = SimplePlot()
    pvx = sp.addPoints(st,x)
    pvy = sp.addPoints(st,y)
    pvx.setLineColor(Color.RED)
    pvy.setLineWidth(3)
    sp.setHLimits(-0.05,5.05)
    sp.setHLabel("Time (s)")
    sp.setVLabel("Displacement (cm)")
    #sp.setTitle(info[edges][1])
    #sp.setFontSizeForPrint(12,504)
    sp.setFontSizeForSlide(1.0,0.9)
    sp.setSize(800,500)
    if pngDir:
      sp.paintToPng(300,7.0,pngDir+info[edges][0]+".png")

def sigmaToA(sigma):
  ss = sigma*sigma
  return (1.0+ss-sqrt(1.0+2.0*ss))/ss

def readSequence(fileName):
  fileName = getDataDir()+fileName
  s = Scanner(BufferedReader(FileReader(fileName)))
  s.findWithinHorizon("nx1=",0)
  nt = s.nextInt()
  s.findWithinHorizon("dx1=",0)
  dt = s.nextDouble()
  s.findWithinHorizon("fx1=",0)
  ft = s.nextDouble()
  s.findWithinHorizon("f=",0)
  s.nextLine()
  x = zerofloat(nt)
  for it in range(nt):
    x[it] = s.nextFloat()
    s.nextLine()
  s.close()
  st = Sampling(nt,dt,ft)
  return st,x

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

