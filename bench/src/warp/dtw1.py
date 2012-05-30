#############################################################################
# Dynamic time warping for 1D sequences

from imports import *

#############################################################################

#pngDir = "./png"
pngDir = None

seed = 99 
seed = 954 
seed = 2127 
seed = abs(Random().nextInt()/1000000)
seed = 877 
print "seed =",seed

nrms = 0.5
strainMax = 0.2

def main(args):
  #goCostMatrix()
  #goDistanceMatrix()
  #goDistance()
  #goAccumulate()
  goFindShifts()

def goFindShifts():
  ml = 33
  f,g,s = makeSequences()
  dw = DynamicWarping(-ml,ml)
  dw.setStrainMax(strainMax)
  dw.setErrorSmoothing(0)
  e = dw.computeErrors(f,g)
  u = dw.findShifts(f,g)
  dw.smoothShifts(u,u)
  plot(f,g,etran(e),s,u)
  dw.setErrorSmoothing(1)
  e = dw.smoothErrors(e)
  u = dw.findShifts(f,g)
  dw.smoothShifts(u,u)
  plot(f,g,etran(e),s,u)

def goAccumulate():
  ml = 33
  f,g,s = makeSequences()
  dw = DynamicWarping(-ml,ml)
  dw.setStrainMax(strainMax)
  e = dw.computeErrors(f,g)
  #e = normalize(e)
  d = dw.accumulateForward(e)
  #u = dw.backtrackReverse(d,e)
  u1 = dw.backtrackReverse(d,e)
  #u = dw.smooth(u)
  es = dw.smoothErrors(e)
  #es = normalize(es)
  ds = dw.accumulateForward(es)
  #u = dw.backtrackReverse(ds,es)
  u11 = dw.backtrackReverse(ds,es)
  print "max|u11-u1| = ",max(abs(sub(u11,u1)))
  #u = dw.smooth(u)
  print " esumu =",dw.sumErrors(e,u1), "  esums =",dw.sumErrors(e,s);
  print "essumu =",dw.sumErrors(es,u11)," essums =",dw.sumErrors(es,s);
  plot(f,g,etran(e))
  plot(f,g,etran(d))
  plot(f,g,etran(es))
  plot(f,g,etran(e),s,u1)
  plot(f,g,etran(es),s,u11)

def smooth(u):
  v = copy(u)
  rgf = RecursiveGaussianFilter(4)
  rgf.apply0(u,v)
  return v

def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

def etran(e):
  e = normalize(e)
  return transpose(pow(e,0.25))

def goDistance():
  ref = RecursiveExponentialFilter(8.0)
  ml = 33
  nl = 1+2*ml
  f,g,s = makeSequences()
  c = dtwCost(-ml,ml,f,g)
  #for cl in c:
  #  ref.apply(cl,cl)
  for limit in [1,2,3]:
    d = dtwDistance(c,limit)
    u = dtwShifts(-ml,ml,d)
    ref.apply(u,u)
    #plot(f,g,pow(c,0.5),s,u)
    plot(f,g,pow(d,0.5),s,u)

def makeSequences():
  n = 501
  fpeak = 0.125
  shift = 2.0/fpeak
  #w = Warp1.constant(shift,n)
  w = Warp1.sinusoid(shift,n)
  #f = makeCosine(fpeak,n)
  f = makeRandomEvents(n,seed=seed); 
  g = w.warp(f)
  f = addRickerWavelet(fpeak,f)
  g = addRickerWavelet(fpeak,g)
  f = addNoise(nrms,fpeak,f,seed=10*seed+1)
  g = addNoise(nrms,fpeak,g,seed=10*seed+2)
  s = zerofloat(n)
  for i in range(n):
    s[i] = w.ux(i)
  return f,g,s

def makeCostMatrix(f,g):
  n = len(f)
  c = zerofloat(n,n)
  for i in range(n):
    sub(f,g[i],c[i])
    mul(c[i],c[i],c[i])
  return c

def goCostMatrix():
  f,g,s = makeSequences()
  c = makeCostMatrix(f,g)
  plotWithMatrix(f,g,sqrt(c))

def goDistanceMatrix():
  f,g,s = makeSequences()
  c = makeCostMatrix(f,g)
  d = dtwDistanceMatrix(c)
  uf,ug = dtwWarpMatrix(d)
  #ref = RecursiveExponentialFilter(4.0)
  #ref.apply1(uf,uf)
  #ref.apply1(ug,ug)
  plotWithMatrix(f,g,pow(c,0.25),uf,ug)
  plotWithMatrix(f,g,pow(d,0.25),uf,ug)

def makeCosine(freq,n):
  return cos(mul(2.0*PI*freq,rampfloat(0.0,1.0,n)))

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

def plotWithMatrix(f,g,c,uf,ug):
  n = len(f)
  panel = PlotPanel(2,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.mosaic.setWidthElastic(0,25)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setWidthElastic(1,100)
  panel.mosaic.setHeightElastic(1,100)
  fv = panel.addPoints(0,1,f)
  gv = panel.addPoints(0,1,g)
  gv.setLineColor(Color.RED)
  fv.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  gv.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  gv = panel.addPoints(1,0,g)
  fv = panel.addPoints(1,0,f)
  fv.setLineColor(Color.RED)
  cv = panel.addPixels(1,1,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #cv.setClips(-1.0,1.0)
  cv.setColorModel(ColorMap.GRAY)
  uv = panel.addPoints(1,1,ug,uf)
  uv.setLineColor(Color.RED)
  dv = panel.addPoints(1,1,(0.0,n-1.0),(0.0,n-1.0))
  dv.setLineColor(Color.WHITE)
  panel.setHLabel(1,"sample")
  panel.setVLabel(1,"sample")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,1000)
  frame.setVisible(True)

def plot(f,g,c,s=None,u=None,clip=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,100)
  panel.setVLimits(1,slag.first,slag.last)
  fv = panel.addPoints(0,0,s1,f)
  gv = panel.addPoints(0,0,s1,g)
  gv.setLineColor(Color.RED)
  cv = panel.addPixels(1,0,s1,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if clip:
    cv.setClips(0.0,clip)
  cv.setColorModel(ColorMap.JET)
  if s:
    sv = panel.addPoints(1,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(1,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.setHLabel("sample")
  panel.setVLabel(0,"f & g")
  panel.setVLabel(1,"lag")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
