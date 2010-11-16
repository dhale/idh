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

#############################################################################
# parameters

def j1(x):
  ax = abs(x)
  if ax<8.0:
    xx = x*x
    num = x*(72362614232.0 + xx*(-7895059235.0 +
                             xx*(242396853.1 +
                             xx*(-2972611.439 +
                             xx*(15704.48260+
                             xx*(-30.16036606))))))
    den = 144725228442.0 + xx*(2300535178.0 +
                           xx*(18583304.74 +
                           xx*(99447.43394 +
                           xx*(376.9991397 +
                           xx))))
    y = num/den
  else:
    z = 8.0/ax
    zz = z*z
    t1 = 1.0 + zz*(0.183105e-2 +
               zz*(-0.3516396496e-4 +
               zz*(0.2457520174e-5 +
               zz*(-0.240337019e-6))))
    t2 = 0.04687499995 + zz*(-0.2002690873e-3 +
                         zz*(0.8449199096e-5 +
                         zz*(-0.88228987e-6 +
                         zz*0.105787412e-6)))
    am = ax-2.356194491
    y = sqrt(0.636619772/ax)*(cos(am)*t1-z*sin(am)*t2)
    if x<0.0:
      y = -y
  return y

def h1(x):
  if x==0.0:
    return 1.0
  else:
    pix = PI*x
    return sin(pix)/(pix)

def h2(r):
  if r==0.0:
    return PI/4.0
  else:
    pir = PI*r
    return j1(pir)/(2.0*r)

def h3(r):
  if r==0.0:
    return PI/6.0
  else:
    pir = PI*r
    return (sin(pir)-pir*cos(pir))/(2.0*r*pir*pir)

def showIdealLowpassFilters():
  nr,dr,fr = 1001,0.01,0.0
  r = rampfloat(fr,dr,nr)
  f1 = zerofloat(nr)
  f2 = zerofloat(nr)
  f3 = zerofloat(nr)
  for i in range(nr):
    f1[i] = h1(r[i])
    f2[i] = h2(r[i])
    f3[i] = h3(r[i])
  sp = SimplePlot()
  pv = sp.addPoints(r,f1); pv.setLineColor(Color.BLACK)
  pv = sp.addPoints(r,f2); pv.setLineColor(Color.RED)
  pv = sp.addPoints(r,f3); pv.setLineColor(Color.BLUE)

def showH2():
  k1max,k2max = 0.4,0.2
  #kw1 = KaiserWindow.fromErrorAndWidth(0.001,1.0-2.0*k1max)
  #kw2 = KaiserWindow.fromErrorAndWidth(0.001,1.0-2.0*k2max)
  #kw1 = KaiserWindow.fromErrorAndWidth(0.001,0.5-k1max)
  #kw2 = KaiserWindow.fromErrorAndWidth(0.001,0.5-k2max)
  #kw1 = KaiserWindow.fromErrorAndWidth(0.001,0.05)
  #kw2 = KaiserWindow.fromErrorAndWidth(0.001,0.05)
  kw1 = KaiserWindow.fromErrorAndLength(0.001,19)
  kw2 = KaiserWindow.fromErrorAndLength(0.001,19)
  n1,n2 = int(kw1.getLength()),int(kw2.getLength())
  print "n1 =",n1," n2 =",n2
  n1,n2 = n1/2*2+1,n2/2*2+1
  print "n1 =",n1," n2 =",n2
  d1,d2 = 1.0,1.0
  f1,f2 = -(n1-1)/2,-(n2-1)/2
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  h = zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = s2.getValue(i2)
    y2 = 2.0*k2max*x2
    w2 = 2.0*k2max*kw2.evaluate(x2)
    for i1 in range(n1):
      x1 = s1.getValue(i1)
      y1 = 2.0*k1max*x1
      w1 = 2.0*k1max*kw1.evaluate(x1)
      r = sqrt(y1*y1+y2*y2)
      h[i2][i1] = w1*w2*h2(r)
  #h = mul(h,1.0/sum(h))
  #h = ((0.0625,0.1250,0.0625),(0.1250,0.2500,0.1250),(0.0625,0.1250,0.0625))
  #h = ((-0.5,0,-0.5),(0,2,0),(-0.5,0,-0.5))
  #h = ((0,0.125,0),(0.125,0.500,0.125),(0,0.125,0))
  print "h: sum =",sum(h)
  dump(copy(1+n1/2,1+n2/2,h))
  plotH(s1,s2,h)
  plotA(s1,s2,h)

def main(args):
  showH2()
 
#############################################################################
# plot

def plotH(s1,s2,h):
  p = panel()
  p.addColorBar("Amplitude")
  p.setColorBarWidthMinimum(140)
  pv = p.addPixels(s1,s2,h)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
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
  cv = p.addContours(sk1,sk2,ak)
  cv.setLineColor(Color.BLACK)
  cv.setContours(Sampling(8,0.1,0.1))
  frame(p,820,634)
  SimplePlot.asPoints(sk1,ak[sk2.count/2])

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
