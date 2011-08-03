import sys
from math import *
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

from lss import *

#############################################################################
# parameters

smNone = LocalSemblanceFilterX.Smoothing.NONE
smBoxcar = LocalSemblanceFilterX.Smoothing.BOXCAR
smGaussian = LocalSemblanceFilterX.Smoothing.GAUSSIAN
smLaplacian = LocalSemblanceFilterX.Smoothing.LAPLACIAN

plotWidth = 800
plotHeight = 600
plotPngDir = "./png/"
#plotPngDir = None

#############################################################################
# functions

def main(args):
  goSmoothingFilters()

def goSmoothingFilters():
  n1 = 43
  hw = 10
  h1s,a1s = [],[]
  sms = [smLaplacian,smGaussian,smBoxcar]
  nms = ["l","g","b"]
  for i in range(len(sms)):
    lsf = LocalSemblanceFilterX(sms[i],hw,smNone,hw)
    h1 = lsf.smooth1(makeImpulse(n1))
    a1 = spectrum(h1)
    h1s.append(h1)
    a1s.append(a1)
    plotSequence(h1,"h1"+nms[i])
  plotSequences(h1s,"h1s")
  plotSpectra(a1s,"a1s")

def makeImpulse(n1):
  f = zerofloat(n1)
  f[n1/2] = 1.0
  return f

def spectrum(h1):
  n1 = len(h1)
  nfft = FftReal.nfftSmall(n1*100)
  fft = FftReal(nfft)
  npad = nfft+2
  hpad = zerofloat(npad)
  copy(h1,hpad)
  fft.realToComplex(-1,hpad,hpad)
  a1  = cabs(hpad)
  return a1
 
#############################################################################
# plot

def plotSequence(f,png=None):
  n1 = len(f)
  s1 = Sampling(n1,1.0,-(n1-1)/2.0)
  p = panel()
  p.setHLabel("sample index")
  p.setHInterval(10.0)
  p.setVLimits(-0.002,0.119)
  sv = p.addSequence(s1,f)
  frame(p,png)

def plotSequences(fs,png=None):
  nf = len(fs)
  n1 = len(fs[0])
  s1 = Sampling(n1,1.0,-(n1-1)/2.0)
  cs = [Color.BLUE,Color.RED,Color.BLACK]
  p = panel()
  p.setHLabel("sample index")
  p.setHInterval(10.0)
  p.setVLimits(-0.002,0.119)
  for i in range(nf):
    sv = p.addSequence(s1,fs[i])
    sv.setColor(cs[i])
  frame(p,png)

def plotSpectra(fs,png=None):
  ns = len(fs)
  nk = len(fs[0])
  dk = 0.5/(nk-1)
  fk = 0.0
  sk = Sampling(nk,dk,fk)
  cs = [Color.BLUE,Color.RED,Color.BLACK]
  p = panel()
  p.setHLimits(0.0,0.1)
  p.setVLimits(0.0,1.0)
  p.setHLabel("frequency (cycles/sample)")
  for i in range(ns):
    pv = p.addPoints(sk,fs[i])
    pv.setLineWidth(5)
    pv.setLineColor(cs[i])
  frame(p,png)

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSizeForSlide(1.0,0.9)
  frame.setSize(plotWidth,plotHeight)
  frame.setVisible(True)
  if png and plotPngDir:
    frame.paintToPng(720,3.3,plotPngDir+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
