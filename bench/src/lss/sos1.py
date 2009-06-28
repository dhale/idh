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

smNone = LocalSemblanceFilter.Smoothing.NONE
smBoxcar = LocalSemblanceFilter.Smoothing.BOXCAR
smGaussian = LocalSemblanceFilter.Smoothing.GAUSSIAN
smLaplacian = LocalSemblanceFilter.Smoothing.LAPLACIAN

plotWidth = 800
plotHeight = 600
plotFontSize = 28
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
  for sm in [smLaplacian,smGaussian,smBoxcar]:
    lsf = LocalSemblanceFilter(sm,hw,smNone,hw)
    h1 = lsf.smooth1(makeImpulse(n1))
    a1 = spectrum(h1)
    h1s.append(h1)
    a1s.append(a1)
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

def plotSequences(fs,png=None):
  nf = len(fs)
  n1 = len(fs[0])
  s1 = Sampling(n1,1.0,-(n1-1)/2.0)
  cs = [Color.BLUE,Color.RED,Color.BLACK]
  p = panel()
  #p.setHLabel("sample index")
  p.setHLabel(" ")
  p.setHInterval(10.0)
  for i in range(nf):
    sv = p.addSequence(s1,fs[i])
    sv.setColor(cs[i])
  frame(p,png)

def plotSequencesX(fs,png=None):
  nf = len(fs)
  n1 = len(fs[0])
  s1 = Sampling(n1,1.0,-(n1-1)/2.0)
  cs = [Color.BLUE,Color.RED,Color.BLACK]
  ms = [PointsView.Mark.FILLED_CIRCLE,
        PointsView.Mark.FILLED_CIRCLE,
        PointsView.Mark.FILLED_CIRCLE]
  p = panel()
  #p.setHLabel("sample index")
  p.setHLabel(" ")
  for i in range(nf):
    pv = p.addPoints(s1,fs[i])
    pv.setMarkColor(cs[i])
    pv.setMarkStyle(ms[i])
    pv.setMarkSize(12)
    pv.setLineStyle(PointsView.Line.NONE)
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
  #p.setHLabel("frequency (cycles/sample)")
  p.setHLabel(" ")
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
  frame.setFontSize(plotFontSize)
  frame.setSize(plotWidth,plotHeight)
  frame.setVisible(True)
  if png and plotPngDir:
    frame.paintToPng(200,6,plotPngDir+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
