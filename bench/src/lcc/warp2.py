import sys
from java.awt import *
from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from lcc import *

#############################################################################
# parameters

fontSize = 24
#width = 640
#height = 505
#widthColorBar = 80
width = 600
height = 622
widthColorBar = 0
dataDir = "/data"
#pngDir = "."
pngDir = None

n1 = 315
n2 = 315

d1max = 3.0
d2max = 3.0
disp = Displacement2.gaussian(d1max,d2max,n1,n2)

lmax = 7
lmin = -lmax

lcfSigma = 12.0
lcfType = LocalCorrelationFilter.Type.SIMPLE
lcfWindow = LocalCorrelationFilter.Window.GAUSSIAN
lcf = LocalCorrelationFilter(lcfType,lcfWindow,lcfSigma)

#############################################################################
# functions

def main(args):
  #goImages()
  #goLcc()
  #goLagSearch()
  #goSequentialShifts()
  #doExact()
  goErrors()
  return

def goLcc():
  f,g = doImages()
  doLcc(f,g,False,False)
  doLcc(f,g,True,False)
  doLcc(f,g,True,True)

def goImages():
  f,g = doImages()
  preprocess(f,g,True,False)
  preprocess(f,g,True,True)

def goSequentialShifts():
  f,g = doImages()
  doSequentialShifts(f,g,True,True,True)
  #doSequentialShifts(f,g,True,True,False)

def goLagSearch():
  f,g = doImages()
  doLagSearch(f,g,False,False)
  doLagSearch(f,g,True,False)
  doLagSearch(f,g,True,True)

def doExact():
  e1 = disp.u1x()
  e2 = disp.u2x()
  plotu(e1,d1max,"e1")
  plotu(e2,d2max,"e2")

def doImages():
  f = readImage()
  g = warpImage(f)
  plot(f,0.0,"f")
  plot(g,0.0,"g")
  return f,g

"""
3*3*5*7
15*21
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
  1   3   5   7   9    11    13    15    17    19
correlation tiles now sampled for
  20  50  80  110  140 170   200   230   260   290 <= xcor tiles
  20          110            200               290 <= axes labels
was
  15  45  75  105  135 165   195   225   255   285
"""
def doLcc(f,g,whiten,smooth,tail=""):
  f,g,suffix = preprocess(f,g,whiten,smooth)
  l1 = lmax
  l2 = lmax
  m1 = 1+2*l1
  m2 = 1+2*l2
  lcf.setInputs(f,g)
  c = zerofloat(10*m1,10*m2)
  t = zerofloat(n1,n2)
  #k1 = ( 38, 98,113,150, 98,172)
  #k2 = (172,216,230,105,132,227)
  #nk = len(k1)
  #ck = zerofloat(m1,m2,nk)
  for lag2 in range(-l2,l2+1):
    for lag1 in range(-l1,l1+1):
      lcf.correlate(lag1,lag2,t)
      lcf.normalize(lag1,lag2,t)
      #copy(n1/m1,n2/m2,l1,l2,m1,m2,t,l1+lag1,l2+lag2,m1,m2,c)
      copy(10,10,20,20,2*m1,2*m2,t,l1+lag1,l2+lag2,m1,m2,c)
      #for k in range(nk):
      #  ck[k][l2+lag2][l1+lag1] = t[k2[k]][k1[k]]
  #c = copy(10*m1,10*m2,m1,m2,c)
  plotLccAllTiles(c,"lcca"+suffix+tail)
  c = copy(m1,m2,3*m1,3*m2,c)
  plotLccOneTile(c,"lcc1"+suffix+tail)
  #for k in range(nk):
  #  plot(ck[k],0.0,"lcc"+suffix+"_"+str(k1[k])+"_"+str(k2[k]))

def getExactU():
  u1 = disp.u1x()
  u2 = disp.u2x()
  #plotu(u1,d1max)
  #plotu(u2,d2max)
  return u1,u2

def getQuadraticFitU(f,g,whiten,smooth):
  f,g,suffix = preprocess(f,g,whiten,smooth)
  lcf.setInputs(f,g)
  l1 = zerobyte(n1,n2)
  l2 = zerobyte(n1,n2)
  lcf.findMaxLags(lmin,lmax,lmin,lmax,l1,l2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lcf.refineLags(l1,l2,u1,u2)
  #plotu(u1,d1max)
  #plotu(u2,d2max)
  return u1,u2

def getCyclicSearchU(f,g,whiten,smooth,interp=True):
  f,g,suffix = preprocess(f,g,whiten,smooth)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  du = zerofloat(n1,n2)
  h = copy(g)
  sf = ShiftFinder(lcfSigma)
  sf.setInterpolateDisplacements(interp)
  for iter in range(4):
    sf.find1(lmin,lmax,f,h,du)
    sf.shift1(du,u1,u2,h)
    sf.find2(lmin,lmax,f,h,du)
    sf.shift2(du,u1,u2,h)
    plotu(u1,d1max)
    plotu(u2,d2max)
  return u1,u2

def goErrors():
  f,g = doImages()
  u1e,u2e = getExactU()
  u1p,u2p = getQuadraticFitU(f,g,False,False)
  f,g,suffix = preprocess(f,g,True,True)
  u1q,u2q = getQuadraticFitU(f,g,False,False)
  u1c,u2c = getCyclicSearchU(f,g,False,False)
  rms1p,rms2p = rmsError(u1p,u1e),rmsError(u2p,u2e)
  max1p,max2p = maxError(u1p,u1e),maxError(u2p,u2e)
  rms1q,rms2q = rmsError(u1q,u1e),rmsError(u2q,u2e)
  max1q,max2q = maxError(u1q,u1e),maxError(u2q,u2e)
  rms1c,rms2c = rmsError(u1c,u1e),rmsError(u2c,u2e)
  max1c,max2c = maxError(u1c,u1e),maxError(u2c,u2e)
  rms1r,rms2r = rms1c/rms1q,rms2c/rms2q
  max1r,max2r = max1c/max1q,max2c/max2q
  print "rms error p",rms1p,rms2p
  print "rms error q",rms1q,rms2q
  print "rms error c",rms1c,rms2c
  print "rms error c/q",rms1r,rms2r
  print "max error p",max1p,max2p
  print "max error q",max1q,max2q
  print "max error c",max1c,max2c
  print "max error c/q",max1r,max2r
  """
  nbin = 10
  nsum1 = zeroint(nbin)
  nsum2 = zeroint(nbin)
  rms1q,rms2q = zerofloat(nbin),zerofloat(nbin)
  rms1c,rms2c = zerofloat(nbin),zerofloat(nbin)
  sum1q,sum2q = zerofloat(nbin),zerofloat(nbin)
  sum1c,sum2c = zerofloat(nbin),zerofloat(nbin)
  max1q,max2q = zerofloat(nbin),zerofloat(nbin)
  max1c,max2c = zerofloat(nbin),zerofloat(nbin)
  for i2 in range(n2):
    for i1 in range(n1):
      u1i = u1e[i2][i1]
      u2i = u2e[i2][i1]
      e1q = Math.abs(u1q[i2][i1]-u1i)
      e2q = Math.abs(u2q[i2][i1]-u2i)
      e1c = Math.abs(u1c[i2][i1]-u1i)
      e2c = Math.abs(u2c[i2][i1]-u2i)
      j1 = int(nbin*Math.abs(u1i-int(u1i))+0.5)%nbin
      j2 = int(nbin*Math.abs(u2i-int(u2i))+0.5)%nbin
      nsum1[j1] += 1
      nsum2[j2] += 1
      sum1q[j1] += e1q*e1q
      sum2q[j2] += e2q*e2q
      sum1c[j1] += e1c*e1c
      sum2c[j2] += e2c*e2c
      max1q[j1] = Math.max(max1q[j1],e1q)
      max2q[j2] = Math.max(max2q[j2],e2q)
      max1c[j1] = Math.max(max1c[j1],e1c)
      max2c[j2] = Math.max(max2c[j2],e2c)
  for j2 in range(nbin):
    rms2q[j2] = Math.sqrt(sum2q[j2]/nsum2[j2])
    rms2c[j2] = Math.sqrt(sum2c[j2]/nsum2[j2])
  for j1 in range(nbin):
    rms1q[j1] = Math.sqrt(sum1q[j1]/nsum1[j1])
    rms1c[j1] = Math.sqrt(sum1c[j1]/nsum1[j1])
  rms1r = div(rms1c,rms1q)
  rms2r = div(rms2c,rms2q)
  max1r = div(max1c,max1q)
  max2r = div(max2c,max2q)
  print "counts"; dump(nsum1); dump(nsum2)
  print "rms error q"; dump(rms1q); dump(rms2q)
  print "rms error c"; dump(rms1c); dump(rms2c)
  print "rms error c/q"; dump(rms1r); dump(rms2r)
  print "max error q"; dump(max1q); dump(max2q)
  print "max error c"; dump(max1c); dump(max2c)
  print "max error c/q"; dump(max1r); dump(max2r)
  fracs = rampfloat(0.0,1.0/nbin,nbin)
  sp = SimplePlot()
  sp.setSize(691,702)
  sp.setVLimits(0.0,0.16)
  sp.setHLimits(-0.05,0.95)
  sp.setVLabel("rms error in shift (samples)")
  sp.setHLabel("distance to nearest sample (samples)")
  prms1q = sp.addPoints(fracs,rms1q)
  prms1q.setLineStyle(PointsView.Line.NONE)
  prms1q.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  prms1q.setMarkSize(24.0)
  prms2q = sp.addPoints(fracs,rms2q)
  prms2q.setLineStyle(PointsView.Line.NONE)
  prms2q.setMarkStyle(PointsView.Mark.HOLLOW_SQUARE)
  prms2q.setMarkSize(24.0)
  prms1c = sp.addPoints(fracs,rms1c)
  prms1c.setLineStyle(PointsView.Line.NONE)
  prms1c.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  prms1c.setMarkSize(24.0)
  prms2c = sp.addPoints(fracs,rms2c)
  prms2c.setLineStyle(PointsView.Line.NONE)
  prms2c.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
  prms2c.setMarkSize(24.0)
  sp = SimplePlot()
  sp.setSize(691,702)
  sp.setVLimits(0.0,1.0)
  sp.setHLimits(-0.05,0.95)
  sp.setVLabel("max error in shift (samples)")
  sp.setHLabel("distance to nearest sample (samples)")
  pmax1q = sp.addPoints(fracs,max1q)
  pmax1q.setLineStyle(PointsView.Line.NONE)
  pmax1q.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  pmax1q.setMarkSize(24.0)
  pmax2q = sp.addPoints(fracs,max2q)
  pmax2q.setLineStyle(PointsView.Line.NONE)
  pmax2q.setMarkStyle(PointsView.Mark.HOLLOW_SQUARE)
  pmax2q.setMarkSize(24.0)
  pmax1c = sp.addPoints(fracs,max1c)
  pmax1c.setLineStyle(PointsView.Line.NONE)
  pmax1c.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pmax1c.setMarkSize(24.0)
  pmax2c = sp.addPoints(fracs,max2c)
  pmax2c.setLineStyle(PointsView.Line.NONE)
  pmax2c.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
  pmax2c.setMarkSize(24.0)
  """

def maxError(y,x):
  return max(abs(sub(y,x)))

def rmsError(y,x):
  n1,n2 = len(x[0]),len(x)
  e = sub(y,x)
  return Math.sqrt(sum(mul(e,e))/(n1*n2))

def doLagSearch(f,g,whiten,smooth):
  fsave = copy(f)
  f,g,suffix = preprocess(f,g,whiten,smooth)
  lcf.setInputs(f,g)
  l1 = zerobyte(n1,n2)
  l2 = zerobyte(n1,n2)
  lcf.findMaxLags(lmin,lmax,lmin,lmax,l1,l2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lcf.refineLags(l1,l2,u1,u2)
  plotu(u1,d1max,"u1"+suffix+"ls")
  plotu(u2,d2max,"u2"+suffix+"ls")
  return u1,u2
 
def doSequentialShifts(f,g,whiten,smooth,interp=True):
  fsave = copy(f)
  f,g,suffix = preprocess(f,g,whiten,smooth)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  du = zerofloat(n1,n2)
  h = copy(g)
  sf = ShiftFinder(lcfSigma)
  sf.setInterpolateDisplacements(interp)
  for iter in range(4):
    sf.find1(lmin,lmax,f,h,du)
    print "1: du min =",min(du),"max =",max(du)
    sf.shift1(du,u1,u2,h)
    print "1: u1 min =",min(u1),"max =",max(u1)
    doLcc(f,h,False,False,"_1"+str(iter))
    sf.find2(lmin,lmax,f,h,du)
    print "2: du min =",min(du),"max =",max(du)
    sf.shift2(du,u1,u2,h)
    doLcc(f,h,False,False,"_2"+str(iter))
    print "2: u2 min =",min(u2),"max =",max(u2)
    plotu(u1,d1max,"u1"+suffix+"ss"+str(iter))
    plotu(u2,d2max,"u2"+suffix+"ss"+str(iter))
  #plotu(u1,d1max,"u1"+suffix+"ss")
  #plotu(u2,d2max,"u2"+suffix+"ss")
  return u1,u2

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def warpImage(f):
  return disp.warp(f)

def preprocess(f,g,whiten,smooth):
  if not whiten:
    return f,g,""
  sf = ShiftFinder(lcfSigma)
  fw = copy(f)
  gw = copy(g)
  sigma = 0.0
  suffix = "w"
  if smooth:
    sigma = 1.0
    suffix += "s"
  sf.whiten(sigma,f,fw)
  sf.whiten(sigma,g,gw)
  plot(fw,0.0,"f"+suffix)
  plot(gw,0.0,"g"+suffix)
  return fw,gw,suffix

#############################################################################
# plot

def plot(f,clip=0.0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  if n1<50 and n2<50:
    s1 = Sampling(n1,1,-(n1-1)/2)
    s2 = Sampling(n2,1,-(n2-1)/2)
  pv = p.addPixels(s1,s2,f)
  if clip!=0.0:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  clipMin,clipMax = pv.getClipMin(),pv.getClipMax()
  print "clip min =",clipMin," max =",clipMax
  frame(p,png)

def plotLccAllTiles(c,png=None):
  p = panel()
  pv = p.addPixels(c)
  pv.setClips(-1.0,1.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  n1,n2 = len(c[0]),len(c)
  l1,l2 = lmax,lmax
  k1,k2 = 1+2*l1,1+2*l2
  m1,m2 = 1+n1/k1,1+n2/k2
  x1 = zerofloat(2,m2)
  x2 = zerofloat(2,m2)
  for i2 in range(m2):
    x1[i2][0] = -0.5
    x1[i2][1] = n1-0.5
    x2[i2][0] = i2*k2-0.5
    x2[i2][1] = x2[i2][0]
  p.addPoints(x1,x2).setLineWidth(1)
  x1 = zerofloat(2,m1)
  x2 = zerofloat(2,m1)
  for i1 in range(m1):
    x1[i1][0] = i1*k1-0.5
    x1[i1][1] = x1[i1][0]
    x2[i1][0] = -0.5
    x2[i1][1] = n2-0.5
  p.addPoints(x1,x2).setLineWidth(1)
  frame(p,png)

def plotLccOneTile(c,png=None):
  p = panel()
  pv = p.addPixels(c)
  pv.setClips(-1.0,1.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #n1,n2 = len(c[0]),len(c)
  #l1,l2 = lmax,lmax
  #k1,k2 = 1+2*l1,1+2*l2
  #m1,m2 = 1+n1/k1,1+n2/k2
  #x1 = (-0.5,n1-0.5)
  #x2 = (float(lmax),float(lmax))
  #pv = p.addPoints(x1,x2)
  #pv.setLineColor(Color.BLUE)
  #pv.setLineWidth(3)
  #x1 = (float(lmax),float(lmax))
  #x2 = (-0.5,n2-0.5)
  #pv = p.addPoints(x1,x2)
  #pv.setLineColor(Color.BLUE)
  #pv.setLineWidth(3)
  frame(p,png)

def plotu(u,clip=0.0,png=None):
  p = panel()
  pv = p.addPixels(u)
  pv.setColorModel(ColorMap.JET)
  if clip!=0.0:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(1.0,99.0)
  frame(p,png)


def makev(u1,u2,clipv=0.0,lvec=0):
  n1 = len(u1[0])
  n2 = len(u1)
  if lvec==0:
    lvec = min(n1,n2)/21-2
  lvec = 1+(lvec/2)*2
  if clipv==0.0:
    clipv = max(max(u1),max(u2))
  scale = lvec/clipv
  print "scale =",scale
  nv1 = n1/lvec
  nv2 = n2/lvec
  kv1 = int(float(n1-1)/float(nv1-1))
  kv2 = int(float(n2-1)/float(nv2-1))
  jv1 = (n1-1-(nv1-1)*kv1)/2
  jv2 = (n2-1-(nv2-1)*kv2)/2
  v1 = zerofloat(2,nv1*nv2)
  v2 = zerofloat(2,nv1*nv2)
  iv = 0
  for iv2 in range(nv2):
    i2 = jv2+iv2*kv2
    for iv1 in range(nv1):
      i1 = jv1+iv1*kv1
      v1[iv][0] = i1
      v1[iv][1] = i1+u1[i2][i1]*scale
      v2[iv][0] = i2
      v2[iv][1] = i2+u2[i2][i1]*scale
      iv = iv+1
  return v1,v2

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  if widthColorBar>0:
    p.addColorBar()
    p.setColorBarWidthMinimum(widthColorBar)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(100,6,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
