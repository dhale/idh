import sys
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
width = 640
height = 505
widthColorBar = 80
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
  goLcc()
  #goLagSearch()
  #goSequentialShifts()
  return

def goLcc():
  f,g = doImages()
  doLcc(f,g,False,False)
  #doLcc(f,g,True,False)
  doLcc(f,g,True,True)

def goImages():
  f,g = doImages()
  preprocess(f,g,True,False)
  preprocess(f,g,True,True)

def goSequentialShifts():
  f,g = doImages()
  doSequentialShifts(f,g,True,True)

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
  f,g = doImages()
  f,g,suffix = preprocess(f,g,whiten,smooth)
  l1 = lmax
  l2 = lmax
  m1 = 1+2*l1
  m2 = 1+2*l2
  lcf.setInputs(f,g)
  c = Array.zerofloat(n1,n2)
  t = Array.zerofloat(n1,n2)
  k1 = ( 38, 98,113,150, 98,172)
  k2 = (172,216,230,105,132,227)
  nk = len(k1)
  ck = Array.zerofloat(m1,m2,nk)
  for lag2 in range(-l2,l2+1):
    for lag1 in range(-l1,l1+1):
      lcf.correlate(lag1,lag2,t)
      lcf.normalize(lag1,lag2,t)
      Array.copy(n1/m1,n2/m2,l1,l2,m1,m2,t,l1+lag1,l2+lag2,m1,m2,c)
      for k in range(nk):
        ck[k][l2+lag2][l1+lag1] = t[k2[k]][k1[k]]
  """
  for i1 in range(n1):
    for i2 in range(0,n2,m2):
      c[i2][i1] = -1.0
    for i2 in range(n2-1,n2,m2):
      c[i2][i1] = -1.0
  for i2 in range(n2):
    for i1 in range(0,n1,m1):
      c[i2][i1] = -1.0
    for i1 in range(n1-1,n1,m1):
      c[i2][i1] = -1.0
  """
  plot(c,1.0,"lcc"+suffix)
  #for k in range(nk):
  #  plot(ck[k],0.0,"lcc"+suffix+"_"+str(k1[k])+"_"+str(k2[k]))

def doLagSearch(f,g,whiten,smooth):
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
 
def doSequentialShifts(f,g,whiten,smooth):
  f,g,suffix = preprocess(f,g,whiten,smooth)
  u1 = Array.zerofloat(n1,n2)
  u2 = Array.zerofloat(n1,n2)
  du = Array.zerofloat(n1,n2)
  ha = Array.copy(g)
  hb = Array.copy(g)
  sf = ShiftFinder(lcfSigma)
  for iter in range(4):
    sf.find1(lmin,lmax,f,ha,du)
    print "1: du min =",Array.min(du),"max =",Array.max(du)
    sf.shift1(du,ha,hb,u1,u2)
    print "1: u1 min =",Array.min(u1),"max =",Array.max(u1)
    ht = ha; ha = hb; hb = ht
    sf.find2(lmin,lmax,f,ha,du)
    print "2: du min =",Array.min(du),"max =",Array.max(du)
    sf.shift2(du,ha,hb,u1,u2)
    print "2: u2 min =",Array.min(u2),"max =",Array.max(u2)
    ht = ha; ha = hb; hb = ht
    plotu(u1,d1max,"u1"+suffix+"ss"+str(iter))
    plotu(u2,d2max,"u2"+suffix+"ss"+str(iter))
  plotu(u1,d1max,"u1"+suffix+"ss")
  plotu(u2,d2max,"u2"+suffix+"ss")
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

def plotu(u,clip=0.0,png=None):
    p = panel()
    pv = p.addPixels(u)
    pv.setColorModel(ColorMap.JET)
    if clip!=0.0:
      pv.setClips(-clip,clip)
    else:
      pv.setPercentiles(1.0,99.0)
    frame(p,png)

def panel():
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
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
    frame.paintToPng(200,6,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
