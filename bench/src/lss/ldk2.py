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

from lss import *

D22 = LocalDiffusionKernel.Stencil.D22
D24 = LocalDiffusionKernel.Stencil.D24
D33 = LocalDiffusionKernel.Stencil.D33
D71 = LocalDiffusionKernel.Stencil.D71

#############################################################################

def main(args):
  #benchStencils()
  showStencils()
  #showDerivativeApproximations()

def benchStencils():
  stencils = [D22,D24,D33,D71]
  n1,n2 = 1000,1000
  s1,s2 = Sampling(n1),Sampling(n2)
  x = randfloat(n1,n2)
  h = zerofloat(n1,n2)
  d = makeConstantTensors(n1,n2,1.0,22)
  for stencil in stencils:
    ldk = LocalDiffusionKernel(stencil)
    sw = Stopwatch()
    sw.restart()
    while sw.time()<1.0: # warmup
      zero(h)
      ldk.apply(d,x,h)
    n = 0
    sw.restart()
    while sw.time()<2.0:
      zero(h)
      ldk.apply(d,x,h)
      n += 1
    sw.stop()
    print "stencil =",stencil," rate =",int(n/sw.time()),"megasamples/s"
    plot(s1,s2,h)

def showStencils():
  stencils = [D22,D24,D33,D71]
  n1,n2 = 13,13
  s1,s2 = Sampling(n1,1,-(n1-1)/2),Sampling(n2,1,-(n2-1)/2)
  x = zerofloat(n1,n2); x[n2/2][n1/2] = 1.0
  d = makeConstantTensors(n1,n2,1.0,22)
  nh = len(stencils)
  h = zerofloat(n1,n2,nh)
  for ih,stencil in enumerate(stencils):
    ldk = LocalDiffusionKernel(stencil)
    ldk.apply(d,x,h[ih])
    plot(s1,s2,h[ih])
  for ih in range(nh):
    plotA(s1,s2,h[ih])

def showDerivativeApproximations():
  for m in [1,2,3,4,6]:
    nh = m+1+m
    kw = KaiserWindow.fromErrorAndLength(0.002,nh)
    print "nh =",nh," m =",m
    sh = Sampling(nh,1,-m)
    h = zerofloat(nh)
    for i in range(1,m+1):
      wi = kw.evaluate(i)
      h[m+i] = wi*cos(PI*i)/i
      h[m-i] = -h[m+i]
    dump(h)
    SimplePlot.asSequence(sh,h)
    fft = Fft(sh)
    fft.setPadding(500)
    sk = fft.getFrequencySampling1()
    hk = fft.applyForward(h)
    ak = div(cabs(hk),2*PI)
    SimplePlot.asPoints(sk,ak)

def makeConstantTensors(n1,n2,aniso,angle):
  angle *= PI/180
  ca = cos(angle)
  sa = sin(angle)
  u1 = fillfloat(ca,n1,n2)
  u2 = fillfloat(sa,n1,n2)
  au = fillfloat(1-aniso,n1,n2)
  av = fillfloat(1,n1,n2)
  return EigenTensors2(u1,u2,au,av)
 
#############################################################################
# plot

def plot(s1,s2,h):
  p = panel()
  p.addColorBar("Amplitude")
  p.setColorBarWidthMinimum(140)
  pv = p.addPixels(s1,s2,h)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #pv.setClips(-0.001,0.001)
  frame(p,810,634)

def plotA(s1,s2,h):
  fft = Fft(s1,s2)
  fft.setCenter(True)
  fft.setPadding(200)
  sk1 = fft.getFrequencySampling1()
  sk2 = fft.getFrequencySampling2()
  hk = fft.applyForward(h)
  ak = cabs(hk)
  print "ak: min =",min(ak)," max =",max(ak)
  p = panel()
  p.addColorBar("Amplitude")
  p.setColorBarWidthMinimum(140)
  pv = p.addPixels(sk1,sk2,ak)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(0.0,2.5)
  cv = p.addContours(sk1,sk2,ak)
  cv.setLineColor(Color.BLACK)
  cv.setContours(Sampling(8,0.1,0.1))
  frame(p,820,634)
  #SimplePlot.asPoints(sk1,ak[sk2.count/2])

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  #p = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.LEFT_TOP)
  #p.setHInterval(100)
  #p.setVInterval(100)
  #p.setHLabel(hlabel)
  #p.setVLabel(vlabel)
  #p.setHLabel(" ")
  #p.setVLabel(" ")
  return p

def frame(panel,width,height,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and plotPngDir:
    frame.paintToPng(200,6,plotPngDir+plotPref+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
