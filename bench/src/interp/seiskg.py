#############################################################################
# Seismic-guided simulation and kriging.

from shared import *

s1,s2 = None,None
n1,d1,f1 = None,None,None
n2,d2,f2 = None,None,None

def main(args):
  makePnzSamplings()
  #goBilateral()
  #goTensors()
  #goSimulate()
  goKriging()

def goTensors():
  g,s1,s2 = readPnzImage()
  et = makeStructureTensors(s1,s2,g)
  plot(None,None,None,s1,s2,g,et=et)

def goSimulate():
  g,s1,s2 = readPnzImage()
  et = makeStructureTensors(s1,s2,g)
  #et = makeIdentityTensors(s1,s2)
  p,px,x1,x2 = simulatePorosity(s1,s2,et)
  pmin,pmax = min(p),max(p)
  s1s = Sampling(1+s1.count/5,5*s1.delta,s1.first)
  s2s = Sampling(1+s2.count/5,5*s2.delta,s2.first)
  ps = gridSimple(px,x1,x2,s1s,s2s)
  plot(None,None,None,s1,s2,g,cmap=gray,clab="Seismic amplitude")
  plot(None,None,None,s1,s2,p,cmap=jet,clab="Porosity")
  #plot(px,x1,x2,s1,s2,p,cmin=pmin,cmax=pmax,cmap=jet,clab="Porosity")
  plot(None,None,None,s1s,s2s,ps,cmin=pmin,cmax=pmax,cmap=tjet,clab="Porosity")

def goKriging():
  g,s1,s2 = readPnzImage()
  #plot(None,None,None,s1,s2,g,cmap=gray,clab="Seismic amplitude")
  et = makeStructureTensors(s1,s2,g)
  #et = makeIdentityTensors(s1,s2)
  p,px,x1,x2 = simulatePorosity(s1,s2,et)
  pmin,pmax = min(p),max(p)
  plot(None,x1,x2,s1,s2,p,mv=True,cmin=pmin,cmax=pmax,cmap=jet,clab="Porosity")
  #s1s = Sampling(1+s1.count/5,5*s1.delta,s1.first)
  #s2s = Sampling(1+s2.count/5,5*s2.delta,s2.first)
  s1s,s2s = s1,s2
  ps = gridSimple(px,x1,x2,s1s,s2s)
  #plot(None,None,None,s1s,s2s,ps,cmin=pmin,cmax=pmax,cmap=tjet,clab="Porosity")
  #et = makeIdentityTensors(s1,s2) # use wrong tensors!
  for pa in [False,True]:
    pk = gridKriging(px,x1,x2,s1,s2,et,pa)
    print "pa =",pa
    printMaxError(px,x1,x2,s1,s2,pk)
    plot(None,None,None,s1,s2,pk,cmin=pmin,cmax=pmax,cmap=jet,clab="Porosity")
    print "rms error =",rmsError(pk,p)

sigmaM,shapeM,rangeM = 1.0,1.0,1.0
sigmaD = 0.00*sigmaM

def simulatePorosity(s1,s2,et):
  sc = SmoothCovariance(sigmaM,shapeM,rangeM,2)
  r = makeRandomImage()
  p = copy(r)
  sc.applyHalf(s1,s2,et,p)
  p = add(0.25,mul(0.02,p))
  px,x1,x2 = randomSamples(200,s1,s2,p)
  return p,px,x1,x2

def gridKriging(f,x1,x2,s1,s2,et,pa):
  cm = SmoothCovariance(sigmaM,shapeM,rangeM,2)
  #cm.testSpd(s1.count,s2.count,et)
  kg = KrigingGridder2(f,x1,x2)
  kg.setModelCovariance(cm);
  kg.setDataError(sigmaD)
  kg.setPolyTrend(0)
  kg.setPaciorek(pa)
  if pa:
    kg.setTensors(et)
  else:
    kg.setTensors(et)
  return kg.grid(s1,s2)

def rmsError(x,y):
  n1,n2 = len(x[0]),len(x)
  z = sub(x,y)
  return sqrt(sum(mul(z,z))/n1/n2)

def printMaxError(px,x1,x2,s1,s2,p):
  emax = -1.0
  i1max,i2max = -1,-1
  n = len(px)
  for i in range(n):
    i1 = s1.indexOfNearest(x1[i])
    i2 = s2.indexOfNearest(x2[i])
    e = abs(p[i2][i1]-px[i])
    if e>emax: 
      emax = e
      i1max = i1
      i2max = i2
  print "emax =",emax
  #print "i1max =",i1max," i2max =",i2max
  #print "x1max =",s1.getValue(i1max)," x2max =",s2.getValue(i2max)

seed = 885
#seed = -1917
#seed = 2046
#r = Random()
#seed = r.nextInt()/1000000
print "seed =",seed
def makeRandomImage():
  r = Random(seed)
  f = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      f[i2][i1] = r.nextGaussian()
  RecursiveGaussianFilter(1.0).apply00(f,f)
  return f
def randomSamples(n,s1,s2,p):
  r = Random(seed)
  n1,n2 = s1.count,s2.count
  f = zerofloat(n)
  x1 = zerofloat(n)
  x2 = zerofloat(n)
  mark = zeroint(n1,n2)
  for i in range(n):
    marked = False
    while not marked:
      i1 = r.nextInt(n1)
      i2 = r.nextInt(n2)
      ok = True
      m = 5
      for j2 in range(max(0,i2-m),min(n2,i2+m+1)):
        for j1 in range(max(0,i1-m),min(n1,i1+m+1)):
          if mark[j2][j1]>0:
            ok = False
      if ok:
        f[i] = p[i2][i1]
        x1[i] = s1.getValue(i1)
        x2[i] = s2.getValue(i2)
        mark[i2][i1] = 1
        marked = True
  return f,x1,x2

def makeIdentityTensors(s1,s2):
  n1,n2 = s1.count,s2.count;
  au = fillfloat(1.0,n1,n2)
  av = fillfloat(1.0,n1,n2)
  u1 = fillfloat(1.0,n1,n2)
  u2 = zerofloat(n1,n2)
  return EigenTensors2(u1,u2,au,av)

def makeStructureTensors(s1,s2,f):
  lof = LocalOrientFilter(4.0)
  lof.setGradientSmoothing(1.0)
  s = lof.applyForTensors(f)
  s.invertStructure(0.5,3.0)
  return s

def makePnzSamplings():
  global s1,n1,d1,f1,s2,n2,d2,f2
  s1 = Sampling(251,0.0125,0.0)
  s2 = Sampling(251,0.0125,0.0)
  n1,d1,f1 = s1.count,s1.delta,s1.first
  n2,d2,f2 = s2.count,s2.delta,s2.first

def readPnzImage():
  global s1,s2
  makePnzSamplings()
  f = readImage("/data/seis/pnz/ggs256",501,501)
  mul(1e-5,f,f)
  m1,m2 = s1.count,s2.count
  j1,j2 = 250,0
  f = copy(m1,m2,j1,j2,f)
  s1 = Sampling(m1,s1.delta,0.0)
  s2 = Sampling(m2,s2.delta,0.0)
  return f,s1,s2

def gridSimple(f,x1,x2,s1,s2):
  return SimpleGridder2(f,x1,x2).grid(s1,s2);

def gridBlended(f,x1,x2,s1,s2,smooth=0.5,et=None):
  bg = BlendedGridder2(f,x1,x2)
  bg.setSmoothness(smooth)
  if et:
    bg.setTensors(et)
    #d = zerofloat(3); et.getTensor(0,0,d); bg.setTensor(d[0],d[1],d[2])
  return bg.grid(s1,s2);

def goBilateral():
  sigmaS,sigmaR = 4.0,0.25
  f,s1,s2 = readPnzImage()
  g = applyBilateralFilter(sigmaS,sigmaR,f)
  plot(None,None,None,s1,s2,f,mv=False)
  plot(None,None,None,s1,s2,g,mv=False)

def applyBilateralFilter(sigmaS,sigmaR,f):
  g = copy(f)
  bf = BilateralFilter(sigmaS,sigmaR)
  bf.apply(f,g)
  return g

#############################################################################
# Plotting

pngDir = None # directory to use for png files
#pngDir = "../../png/interp" # directory to use for png files

def plot(f,x1,x2,s1,s2,g,png=None,
         cmin=0,cmax=0,cmap=None,clab=None,mv=True,et=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setSize(1000,810)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  sp.setLimits(s1.first,s2.first,s1.last,s2.last)
  if not clab:
    clab = "Amplitude"
  sp.addColorBar(clab)
  sp.plotPanel.setColorBarWidthMinimum(200)
  #sp.plotPanel.setVInterval(100)
  pv = sp.addPixels(s1,s2,g)
  if cmin==cmax:
    cmin,cmax = min(g),max(g)
  pv.setClips(cmin,cmax)
  if cmap:
    pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if mv and x1 and x2:
    mv = sp.addPoints(x1,x2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkColor(Color.WHITE)
    mv.setMarkSize(5)
  if et:
    tv = TensorsView(s1,s2,et)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2.0)
    tv.setScale(2.0)
    sp.add(tv)
  if png and pngDir: sp.paintToPng(360,3.33,pngDir+png+".png")

def makeTransparentColorModel():
  a = fillfloat(1.0,256); a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a);

gray = ColorMap.GRAY
jet = ColorMap.JET
tjet = makeTransparentColorModel()

#############################################################################
# Data

def readImage(fileName,n1,n2):
  ais = ArrayInputStream(fileName+".dat")
  x = zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

#############################################################################
if __name__ == "__main__":
  run(main)
