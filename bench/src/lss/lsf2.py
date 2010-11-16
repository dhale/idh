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

def main(args):
  #testDerivative()
  #testConstant()
  testCircle()
  #testChirp()

def testDerivative():
  for m in [1,2,3,4,6]:
    nh = m+1+m
    kw = KaiserWindow.fromErrorAndLength(0.005,nh)
    #kw = KaiserWindow.fromErrorAndLength(0.002,nh)
    #kw = KaiserWindow.fromErrorAndWidth(0.005,0.20)
    #m = int(kw.getLength()+1)/2
    #nh = m+1+m
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

def testConstant():
  n1,n2 = 101,101
  s1,s2 = Sampling(n1),Sampling(n2)
  f = zerofloat(n1,n2); f[n2/2][n1/2] = 1.0
  #f[   0][   0] = 1.0
  #f[n2-1][   0] = 1.0
  #f[   0][n1-1] = 1.0
  #f[n2-1][n1-1] = 1.0
  f = applyBandPassFilter(f)
  #plot(s1,s2,f)
  lsf = LocalSmoothingFilterX(0.01,1000)
  for angle in [20]: #[0,15,30,45]:
    for method in [33,22,24,71,91]:
      lsf.METHOD = method
      t = makeConstantTensors(n1,n2,angle)
      g = zerofloat(n1,n2)
      #lsf.applyLhs(t,1.0,None,f,g); sub(g,f,g)
      lsf.apply(t,5.0,f,g)
      #plot(s1,s2,g)
      plotA(s1,s2,g)

def testChirp():
  #kmax = 0.5/sqrt(2.0)
  kmax = 0.5
  n1,n2 = 251,251
  s1,s2 = Sampling(n1),Sampling(n2)
  f = makeChirp(kmax,n1,n2)
  plot(s1,s2,f)
  #plotA(s1,s2,f)
  f = applyBandPassFilter(f)
  plot(s1,s2,f)
  #plotA(s1,s2,f)
  t = makeCircleTensors(n1,n2)
  g = zerofloat(n1,n2)
  lsf = LocalSmoothingFilterX(0.01,1000)
  for method in [22,24,33,91]:
    copy(f,g)
    lsf.METHOD = method
    lsf.apply(t,200.0,f,g)
    plot(s1,s2,g)
    #plotA(s1,s2,g)

def testCircle():
  n1,n2 = 401,401
  kr,kt = 20,40
  s1,s2 = Sampling(n1),Sampling(n2)
  f = makeImpulsesOnCircles(n1,n2,kr,kt)
  t = makeCircleTensors(n1,n2)
  f = applyBandPassFilter(f)
  plot(s1,s2,f)
  lsf = LocalSmoothingFilterX(0.01,1000)
  g = zerofloat(n1,n2)
  for method in [33,712]:
    lsf.METHOD = method
    n = 0
    sw = Stopwatch()
    sw.start()
    while sw.time()<=10.0:
      copy(f,g)
      lsf.apply(t,200.0,f,g)
      n += 1
    sw.stop()
    print "n =",n," rate =",n/sw.time()
    plot(s1,s2,g)

def applyBandPassFilter(f):
  bpf = BandPassFilter(0.00,0.45,0.10,0.01)
  bpf.setExtrapolation(BandPassFilter.Extrapolation.ZERO_SLOPE)
  g = zerofloat(len(f[0]),len(f))
  bpf.apply(f,g)
  return g

def makeImpulsesOnGrid(n1,n2,k1,k2):
  m1 = n1/k1
  m2 = n2/k2
  j1 = (n1-(m1-1)*k1)/2
  j2 = (n2-(m2-1)*k2)/2
  f = zerofloat(n1,n2)
  for i2 in range(m2):
    for i1 in range(m1):
      f[j2+i2*k2][j1+i1*k1] = 1.0
  return f

def makeImpulsesOnCircles(n1,n2,kr,kt):
  f = zerofloat(n1,n2)
  nr = int(sqrt(n1*n1+n2*n2))
  mr = nr/kr
  jr = (nr-(mr-1)*kr)/2
  for ir in range(jr,nr,kr):
    r = ir
    nt = int(0.5*PI*r)
    mt = nt/kt
    jt = (nt-(mt-1)*kt)/2
    for it in range(jt,nt,kt):
      t = 0.5*PI*it/nt
      i1 = int(r*cos(t)+0.5)
      i2 = int(r*sin(t)+0.5)
      if i1<n1 and i2<n2:
        f[i2][i1] = 1.0
  return f

def makeChirp(kmax,n1,n2):
  n = max(n1,n2)
  s = kmax*PI*n
  x1 = rampfloat(0.0,1.0/n,n1)
  x2 = rampfloat(0.0,1.0/n,n2)
  c = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      rr = x1[i1]*x1[i1]+x2[i2]*x2[i2]
      c[i2][i1] = cos(s*rr)
  return c

def makeCircleTensors(n1,n2):
  #c1,c2 = (n1-1.0)/2.0,(n2-1.0)/2.0
  c1,c2 = 0,0
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = i2-c2
    for i1 in range(n1):
      x1 = i1-c1
      xs = x1*x1+x2*x2
      if xs!=0.0:
        xs = 1.0/sqrt(xs)
      else:
        xs = 1.0
      u1[i2][i1] = xs*x1
      u2[i2][i1] = xs*x2
      au[i2][i1] = 0.0
      av[i2][i1] = 1.0
  #plot(Sampling(n1),Sampling(n2),u1)
  #plot(Sampling(n1),Sampling(n2),u2)
  return EigenTensors2(u1,u2,au,av)

def makeConstantTensors(n1,n2,angle):
  angle *= PI/180
  ca = cos(angle)
  sa = sin(angle)
  u1 = fillfloat(ca,n1,n2)
  u2 = fillfloat(sa,n1,n2)
  au = fillfloat(1,n1,n2)
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
  #pv.setClips(-1.0,1.0)
  pv.setClips(-0.025,0.025)
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
  pv.setClips(0.0,1.0)
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
