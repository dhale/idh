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
pngDir = "."
#pngDir = None

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
  doExact()
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

def doLcc(f,g,whiten,smooth):
  f,g,suffix = preprocess(f,g,whiten,smooth)
  l1 = lmax
  l2 = lmax
  m1 = 1+2*l1
  m2 = 1+2*l2
  lcf.setInputs(f,g)
  c = Array.zerofloat(n1,n2)
  t = Array.zerofloat(n1,n2)
  #k1 = ( 38, 98,113,150, 98,172)
  #k2 = (172,216,230,105,132,227)
  #nk = len(k1)
  #ck = Array.zerofloat(m1,m2,nk)
  for lag2 in range(-l2,l2+1):
    for lag1 in range(-l1,l1+1):
      lcf.correlate(lag1,lag2,t)
      lcf.normalize(lag1,lag2,t)
      Array.copy(n1/m1,n2/m2,l1,l2,m1,m2,t,l1+lag1,l2+lag2,m1,m2,c)
      #for k in range(nk):
      #  ck[k][l2+lag2][l1+lag1] = t[k2[k]][k1[k]]
  c = Array.copy(10*m1,10*m2,m1,m2,c)
  plotLccAllTiles(c,"lcca"+suffix)
  c = Array.copy(m1,m2,6*m1,7*m2,c)
  plotLccOneTile(c,"lcc1"+suffix)
  #for k in range(nk):
  #  plot(ck[k],0.0,"lcc"+suffix+"_"+str(k1[k])+"_"+str(k2[k]))

def doLagSearch(f,g,whiten,smooth):
  fsave = Array.copy(f)
  f,g,suffix = preprocess(f,g,whiten,smooth)
  lcf.setInputs(f,g)
  l1 = Array.zerobyte(n1,n2)
  l2 = Array.zerobyte(n1,n2)
  lcf.findMaxLags(lmin,lmax,lmin,lmax,l1,l2)
  u1 = Array.zerofloat(n1,n2)
  u2 = Array.zerofloat(n1,n2)
  lcf.refineLags(l1,l2,u1,u2)
  plotu(u1,d1max,"u1"+suffix+"ls")
  plotu(u2,d2max,"u2"+suffix+"ls")
  return u1,u2
 
def doSequentialShifts(f,g,whiten,smooth,interp=True):
  fsave = Array.copy(f)
  f,g,suffix = preprocess(f,g,whiten,smooth)
  u1 = Array.zerofloat(n1,n2)
  u2 = Array.zerofloat(n1,n2)
  du = Array.zerofloat(n1,n2)
  h = Array.copy(g)
  sf = ShiftFinder(lcfSigma)
  sf.setInterpolateDisplacements(interp)
  for iter in range(4):
    sf.find1(lmin,lmax,f,h,du)
    print "1: du min =",Array.min(du),"max =",Array.max(du)
    sf.shift1(du,u1,u2,h)
    print "1: u1 min =",Array.min(u1),"max =",Array.max(u1)
    sf.find2(lmin,lmax,f,h,du)
    print "2: du min =",Array.min(du),"max =",Array.max(du)
    sf.shift2(du,u1,u2,h)
    print "2: u2 min =",Array.min(u2),"max =",Array.max(u2)
    plotu(u1,d1max,"u1"+suffix+"ss"+str(iter))
    plotu(u2,d2max,"u2"+suffix+"ss"+str(iter))
  #plotu(u1,d1max,"u1"+suffix+"ss")
  #plotu(u2,d2max,"u2"+suffix+"ss")
  return u1,u2

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def warpImage(f):
  return disp.warp(f)

def preprocess(f,g,whiten,smooth):
  if not whiten:
    return f,g,""
  sf = ShiftFinder(lcfSigma)
  fw = Array.copy(f)
  gw = Array.copy(g)
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
  x1 = Array.zerofloat(2,m2)
  x2 = Array.zerofloat(2,m2)
  for i2 in range(m2):
    x1[i2][0] = -0.5
    x1[i2][1] = n1-0.5
    x2[i2][0] = i2*k2-0.5
    x2[i2][1] = x2[i2][0]
  p.addPoints(x1,x2).setLineWidth(2)
  x1 = Array.zerofloat(2,m1)
  x2 = Array.zerofloat(2,m1)
  for i1 in range(m1):
    x1[i1][0] = i1*k1-0.5
    x1[i1][1] = x1[i1][0]
    x2[i1][0] = -0.5
    x2[i1][1] = n2-0.5
  p.addPoints(x1,x2).setLineWidth(2)
  frame(p,png)

def plotLccOneTile(c,png=None):
  p = panel()
  pv = p.addPixels(c)
  pv.setClips(-1.0,1.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  n1,n2 = len(c[0]),len(c)
  l1,l2 = lmax,lmax
  k1,k2 = 1+2*l1,1+2*l2
  m1,m2 = 1+n1/k1,1+n2/k2
  x1 = (-0.5,n1-0.5)
  x2 = (float(lmax),float(lmax))
  pv = p.addPoints(x1,x2)
  pv.setLineColor(Color.BLUE)
  pv.setLineWidth(3)
  x1 = (float(lmax),float(lmax))
  x2 = (-0.5,n2-0.5)
  pv = p.addPoints(x1,x2)
  pv.setLineColor(Color.BLUE)
  pv.setLineWidth(3)
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
    clipv = max(Array.max(u1),Array.max(u2))
  scale = lvec/clipv
  print "scale =",scale
  nv1 = n1/lvec
  nv2 = n2/lvec
  kv1 = int(float(n1-1)/float(nv1-1))
  kv2 = int(float(n2-1)/float(nv2-1))
  jv1 = (n1-1-(nv1-1)*kv1)/2
  jv2 = (n2-1-(nv2-1)*kv2)/2
  v1 = Array.zerofloat(2,nv1*nv2)
  v2 = Array.zerofloat(2,nv1*nv2)
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
