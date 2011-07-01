#############################################################################
# Fault displacements from 2D images

import sys
from java.awt import *
from java.awt.image import *
from java.io import *
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

from lcc import FaultFinder2SB
from dnp import LocalSlopeFinder

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  #goShifts()
  goThin()
  #goScan()

def getFaultFinder():
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 20.0
  return FaultFinder2SB(slopeMax,shiftMax,thetaMax)

def getImage():
  #return imageSyn()
  return imageF3d()
  #return imageTpd()

def goShifts():
  s1,s2,f = getImage()
  n1,n2 = len(f[0]),len(f)
  #plot2(s1,s2,f,title="input")
  f = slog(f)
  plot2(s1,s2,f,title="log input")
  ff = getFaultFinder()
  p = ff.findSlopes(f)
  #plot2(s1,s2,f,p,title="slopes")
  f = ff.taper(10,f)
  snd = ff.semblanceNumDen(p,f)
  ct = ff.faultThetaScan(snd)
  ct = ff.faultThetaThin(ct)
  plot2(s1,s2,f,ct[0],gmin=0,gmax=1.0,title="faults")
  #plot2(s1,s2,f,ct[1],title="thetas")
  g = ff.smooth(16.0,p,ct[0],f);
  plot2(s1,s2,g,title="smoothed")
  u = ff.findShifts(ct,g)
  plot2(s1,s2,g,u,gmin=-3,gmax=3,title="shifts")

def goThin():
  s1,s2,f = getImage()
  n1,n2 = len(f[0]),len(f)
  #plot2(s1,s2,f,title="input")
  f = slog(f)
  plot2(s1,s2,f,title="log input")
  ff = getFaultFinder()
  p = ff.findSlopes(f)
  #plot2(s1,s2,f,p,title="slopes")
  #f = ff.taper(10,f)
  snd = ff.semblanceNumDen(p,f)
  ct = ff.faultThetaScan(snd)
  plot2(s1,s2,f,ct[0],gmin=0,gmax=1.0,title="c scan")
  plot2(s1,s2,f,ct[1],title="t scan")
  ct = ff.faultThetaThin(ct)
  plot2(s1,s2,f,ct[0],gmin=0,gmax=1.0,title="c thin")
  plot2(s1,s2,f,ct[1],title="t thin")
  g = ff.smooth(16.0,p,ct[0],f);
  plot2(s1,s2,g,title="smoothed")
  #plot2(s1,s2,sexp(g),title="exp smoothed")

def goScan():
  ff = getFaultFinder()
  s1,s2,f = getImage()
  n1,n2 = len(f[0]),len(f)
  #plot2(s1,s2,f,title="input")
  f = slog(f)
  plot2(s1,s2,f,title="log input")
  p = ff.findSlopes(f)
  #plot2(s1,s2,f,p,title="slopes")
  f = ff.taper(10,f)
  snd = ff.semblanceNumDen(p,f)
  st = Sampling(41,1.00,-20.0)
  #st = Sampling(1,1.00,0.0)
  nt,dt,ft = st.count,st.delta,st.first
  for jt in range(nt):
    t = st.getValue(jt)
    s = ff.semblance(t,snd)
    c = sub(1.0,pow(s,8))
    title = "theta = %6.3f" % t
    plot2(s1,s2,f,c,gmin=0,gmax=1.0,title=title)

def findShifts(sigma,min1,max1,lag2,f):
  n1,n2 = len(f[0]),len(f)
  u = zerofloat(n1,n2)
  lsf = LocalShiftFinder(sigma)
  lsf.setSmoothShifts(True)
  #lsf.setSmoothShifts(False)
  for i2 in range(n2):
    i2m = max(i2-lag2,0)
    i2p = min(i2+lag2,n2-1)
    lsf.find1(min1,max1,f[i2m],f[i2p],u[i2])
  return u

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))
 
#############################################################################
# data read/write

def readImage(n1,n2,fileName):
  f = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(f,fileName):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(f)
  aos.close()

def imageTpd():
  #ft,fx = 0.500,0.000
  #dt,dx = 0.004,0.025
  f1,f2 = 0.0,0.0
  d1,d2 = 1.0,1.0
  n1,n2 = 251,357
  fileName = "/data/seis/tp/csm/oldslices/tp73.dat"
  x = readImage(n1,n2,fileName)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return s1,s2,x

def imageF3d():
  n1,n2 = 462,951
  d1,d2 = 0.004,0.025
  f1,f2 = 0.004,0.000
  fileName = "/data/seis/f3d/f3d75.dat"
  x = readImage(n1,n2,fileName)
  subset = True
  if subset:
    j1,j2 = 240,0
    n1,n2 = n1-j1,440
    f1,f2 = f1+j1*d1,f2+j2*d2
    x = copy(n1,n2,j1,j2,x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return s1,s2,x

def imageSyn():
  s1,s2,f = imageTpd()
  n1,n2 = s1.count,s2.count
  g = zerofloat(n1,n2)
  f1 = copy(f[n2/8])
  f2 = shiftRamp(f1)
  for i2 in range(n2):
    for i1 in range(n1):
      if i2<n2/4:
        copy(f1,g[i2])
      elif i2<n2/3:
        copy(f2,g[i2])
      elif i2<2*n2/3:
        copy(f1,g[i2])
      else:
        copy(f2,g[i2])
  #zero(g[1*n2/3])
  #zero(g[2*n2/3])
  #r = randomNoise(1.0,n1,n2)
  #g = add(g,r)
  #rgf = RecursiveGaussianFilter(2.0)
  #rgf.applyX0(g,g)
  return s1,s2,g

def shiftRamp(f):
  n = len(f)
  g = copy(f)
  t = rampfloat(0.0,1.0-8.0/(n-1),n)
  si = SincInterpolator()
  si.setUniform(n,1.0,0.0,f)
  si.interpolate(n,t,g)
  return g

def randomNoise(a,n1,n2):
  ran = Random(3)
  r = mul(2.0*a,sub(randfloat(ran,n1,n2),0.5))
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply00(r,r)
  return r
 
#############################################################################
# plotting

def plot2(s1,s2,f,g=None,gmin=None,gmax=None,
          label=None,title=None,png=None):
  n1 = len(f[0])
  n2 = len(f)
  panel = panel2()
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  if title:
    panel.setTitle(title)
  panel.setColorBarWidthMinimum(180)
  #pv = panel.addPixels(s1,s2,f)
  pv = panel.addPixels(f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  #pv.setClips(-4.5,4.5)
  if g:
    alpha = 0.5
    #pv = panel.addPixels(s1,s2,g)
    pv = panel.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    if gmin==None: gmin = min(g)
    if gmax==None: gmax = max(g)
    pv.setClips(gmin,gmax)
    pv.setColorModel(ColorMap.getJet(alpha))
    #updateColorModel(pv,0.9)
  frame2(panel,png)

def updateColorModel(pv,alpha):
    n = 256
    r = zerobyte(n)
    g = zerobyte(n)
    b = zerobyte(n)
    a = zerobyte(n)
    icm = pv.getColorModel()
    icm.getReds(r)
    icm.getGreens(g)
    icm.getBlues(b)
    ia = int(255.0*alpha)
    if ia>128:
      ia -= 256
    for i in range(n):
      a[i] = ia
    if alpha<1.0:
      #r[n/2] = r[n/2-1] = -1
      #g[n/2] = g[n/2-1] = -1
      #b[n/2] = b[n/2-1] = -1
      a[n/2  ] = a[n/2-1] = 0
      a[n/2+1] = a[n/2-2] = 0
      a[n/2+2] = a[n/2-3] = 0
    icm = IndexColorModel(8,n,r,g,b,a)
    pv.setColorModel(icm)

def panel2():
  #panel = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)
  return panel

def frame2(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  frame.setFontSizeForSlide(1.0,0.8)
  #frame.setSize(1290,777)
  #frame.setSize(1490,977)
  frame.setSize(1490,815)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

"""
Outline:
find slopes p
use p to align minus-plus images fm and fp
for all thetas
  shear fm and fp
  compute correlation coefficients c
  unshear c
  compute fault likelihood d
  remember dmax and corresponding theta

smooth dmax laterally and pick peaks
for all fault curves
  gather samples on both sides of fault
  cross-correlate to find displacements
"""
