#############################################################################
# Fault displacements

import sys
from org.python.util import PythonObjectInputStream
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

from lcc import FaultFinder
from dnp import LocalSlopeFinder

"""
Outline:
find slopes p
use p to align minus-plus images fm and fp
for all thetas
  shear fm and fp
  compute correlation coefficients c
  unshear c
  remember cmax and corresponding theta
smooth max c laterally and pick correlation peaks
for all fault curves
  gather samples on both sides of fault
  cross-correlate to find displacements
"""

#############################################################################
# global parameters

dataDir = "/data/seis/tp/csm/oldslices/"
#pngDir = "./png"
pngDir = None

def main(args):
  #goShifts()
  goThetas()
  #goScan()
  #goAlign()
  #goShear()
  #goShiftsX()
  #goFilter()

def goShifts():
  #s1,s2,f = imageSyn()
  #s1,s2,f = imageF3d()
  s1,s2,f = imageTeapot()
  n1,n2 = len(f[0]),len(f)
  plot2(s1,s2,f,title="input")
  #f = spow(0.5,f)
  f = slog(f)
  plot2(s1,s2,f,title="log input")
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 20.0
  ff = FaultFinder(slopeMax,shiftMax,thetaMax)
  p = ff.findSlopes(f)
  #plot2(s1,s2,f,p,title="slopes")
  ct = ff.findThetas(p,f)
  plot2(s1,s2,f,ct[0],title="faults")
  plot2(s1,s2,f,ct[1],title="thetas")
  g = ff.smooth(16.0,p,ct[0],f);
  plot2(s1,s2,g,title="smoothed")
  u = ff.findShifts(ct,g)
  plot2(s1,s2,g,u,gmin=-3,gmax=3,title="shifts")

def goThetas():
  s1,s2,f = imageSyn()
  #s1,s2,f = imageF3d()
  #s1,s2,f = imageTeapot()
  n1,n2 = len(f[0]),len(f)
  plot2(s1,s2,f,title="input")
  f = slog(f)
  plot2(s1,s2,f,title="log input")
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 20.0
  ff = FaultFinder(slopeMax,shiftMax,thetaMax)
  p = ff.findSlopes(f)
  #plot2(s1,s2,f,p,title="slopes")
  fmp = ff.align(1,p,f)
  ct = ff.ctScan(fmp)
  plot2(s1,s2,f,ct[0],title="c scan")
  plot2(s1,s2,f,ct[1],title="t scan")
  ct = ff.ctThin(ct)
  plot2(s1,s2,f,ct[0],title="c thin")
  plot2(s1,s2,f,ct[1],title="t thin")
  g = ff.smooth(16.0,p,ct[0],f);
  plot2(s1,s2,g,title="smoothed")
  plot2(s1,s2,sexp(g),title="exp smoothed")

def goScan():
  s1,s2,f = imageSyn()
  #s1,s2,f = imageF3d()
  #s1,s2,f = imageTeapot()
  n1,n2 = len(f[0]),len(f)
  plot2(s1,s2,f,title="input")
  #f = spow(0.3,f)
  #f = slog(f)
  #plot2(s1,s2,f,title="log input")
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 20.0
  ff = FaultFinder(slopeMax,shiftMax,thetaMax)
  p = ff.findSlopes(f)
  #plot2(s1,s2,f,p,title="slopes")
  fmp = ff.align(1,p,f)
  #st = Sampling(41,1.00,-20.0)
  st = Sampling(1,1.00,0.0)
  nt,dt,ft = st.count,st.delta,st.first
  for jt in range(nt):
    t = st.getValue(jt)
    c = ff.xcor(t,fmp)
    c = sub(1.0,c)
    title = "theta = %6.3f" % t
    plot2(s1,s2,f,c,title=title)

def goShear():
  #s1,s2,f = imageF3d()
  s1,s2,f = imageTeapot()
  n1,n2 = len(f[0]),len(f)
  f = slog(f)
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 45.0
  ff = FaultFinder(slopeMax,shiftMax,thetaMax)
  shear = 0.13
  g = ff.shear(shear,f)
  h = ff.unshear(shear,g)
  e = sub(h,f)
  print "shear-unshear error max =",max(abs(e))
  plot2(s1,s2,f,title="input")
  plot2(s1,s2,g,title="sheared")
  plot2(s1,s2,h,title="unsheared")
  plot2(s1,s2,e,title="error")

def goAlign():
  s1,s2,f = imageSyn()
  #s1,s2,f = imageF3d()
  #s1,s2,f = imageTeapot()
  n1,n2 = len(f[0]),len(f)
  f = slog(f)
  slopeMax = 5.0
  shiftMax = 15.0
  thetaMax = 45.0
  ff = FaultFinder(slopeMax,shiftMax,thetaMax)
  p = ff.findSlopes(f)
  fm,fp = ff.align(1,p,f)
  fe = sub(fm,fp)
  plot2(s1,s2,f,title="input")
  plot2(s1,s2,fm,title="minus")
  plot2(s1,s2,fp,title="plus")
  plot2(s1,s2,fe,title="error")

def goShiftsX():
  #s1,s2,f = imageSyn()
  #s1,s2,f = imageTeapot()
  s1,s2,f = imageF3d()
  f = taper(f)
  n1,n2 = len(f[0]),len(f)
  ff = FaultFinder(0.3,15,10,7)
  fs = ff.shear(-0.14,f)
  f = copy(n1,n2,fs)
  plot2(s1,s2,f)
  u = zerofloat(n1,n2)
  c = zerofloat(n1,n2)
  d = zerofloat(n1,n2)
  mlag = 15
  nlag = 1+2*mlag
  sigma = 3.0*mlag
  for k in [2]:
    u = findShifts(sigma,-mlag,mlag,k,f)
    #w = pow(c,2.0) 
    #us = smooth(0.9,w,u)
    #plot2(s1,s2,f,us)
    plot2(s1,s2,f,u)
    #k2 = 240 # for imageSyn with k=3
    #k2 = 145 # teapot and f3d
    k2 = 192 # sheared f3d
    fm = f[k2-k]
    fp = f[k2+k]
    lcf = LocalCorrelationFilter(
      LocalCorrelationFilter.Type.SIMPLE,
      LocalCorrelationFilter.Window.GAUSSIAN,
      sigma)
    lcf.setInputs(fm,fp)
    c = zerofloat(n1,nlag)
    for ilag in range(nlag):
      lag = ilag-mlag
      lcf.correlate(lag,c[ilag])
      lcf.normalize(lag,c[ilag])
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    sp.setSize(250,780)
    sp.addPoints(fm).setLineColor(Color.RED)
    sp.addPoints(fp).setLineColor(Color.BLUE)
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    sp.addColorBar()
    sp.setSize(500,780)
    sl = Sampling(nlag,1.0,-mlag)
    pv = sp.addPixels(Sampling(n1),sl,c)
    pv.setColorModel(ColorMap.JET)
    pv = sp.addPoints(u[k2])
    pv.setLineWidth(5)

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

def taper(f):
  f = copy(f)
  n1,n2 = len(f[0]),len(f)
  m = 5
  t = fillfloat(1.0,n1)
  for i in range(m):
    ti = 0.5*(1.0-cos(PI*(i+1)/m))
    t[i] = ti
    t[n1-1-i] = ti
  for i2 in range(n2):
    mul(t,f[i2],f[i2])
  return f

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

def smooth(a,w,f):
  aw = mul(a,w)
  bw = sub(1.0,aw)
  cw = sub(1.0,bw)
  fm = smoothM(a,w,f)
  fp = smoothP(a,w,f)
  g = div(sub(add(fm,fp),mul(bw,f)),cw)
  return add(smoothM(a,w,f),smoothP(a,w,f))

def smoothM(a,w,f):
  n1,n2 = len(f[0]),len(f)
  g = copy(f)
  aw = zerofloat(n1)
  bw = zerofloat(n1)
  for i2 in range(1,n2):
    mul(a,w[i2],aw)
    sub(1.0,aw,bw)
    add(mul(bw,g[i2-1]),mul(aw,f[i2]),g[i2]) 
  return g

def smoothP(a,w,f):
  n1,n2 = len(f[0]),len(f)
  g = copy(f)
  aw = zerofloat(n1)
  bw = zerofloat(n1)
  for i2 in range(n2-2,-1,-1):
    mul(a,w[i2],aw)
    sub(1.0,aw,bw)
    add(mul(bw,g[i2+1]),mul(aw,f[i2]),g[i2]) 
  return g

def goFilter():
  #s1,s2,f = imageTeapot()
  s1,s2,f = imageSyn()
  n1,n2 = s1.count,s2.count
  a = 0.95
  fm = filterM(a,f)
  fp = filterP(a,f)
  plot2(s1,s2,f)
  plot2(s1,s2,fm)
  plot2(s1,s2,fp)

def smoothM(a,w,f):
  n1,n2 = len(f[0]),len(f)
  g = copy(f)
  aw = zerofloat(n1)
  bw = zerofloat(n1)
  for i2 in range(1,n2):
    mul(a,w[i2],aw)
    sub(1.0,aw,bw)
    add(mul(bw,g[i2-1]),mul(aw,f[i2]),g[i2]) 
  return g

def smoothP(a,w,f):
  n1,n2 = len(f[0]),len(f)
  g = copy(f)
  aw = zerofloat(n1)
  bw = zerofloat(n1)
  for i2 in range(n2-2,-1,-1):
    mul(a,w[i2],aw)
    sub(1.0,aw,bw)
    add(mul(bw,g[i2+1]),mul(aw,f[i2]),g[i2]) 
  return g

def filterM(a,f):
  n2 = len(f)
  g = copy(f)
  for i2 in range(1,n2):
    add(mul(a,g[i2-1]),mul(1.0-a,f[i2]),g[i2]) 
  return g

def filterP(a,f):
  n2 = len(f)
  g = copy(f)
  for i2 in range(n2-2,-1,-1):
    add(mul(a,g[i2+1]),mul(1.0-a,f[i2]),g[i2]) 
  return g
 
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

def imageTeapot():
  #ft,fx = 0.500,0.000
  #dt,dx = 0.004,0.025
  ft,fx = 0.0,0.0
  dt,dx = 1.0,1.0
  nt,nx = 251,357
  image = zerofloat(nt,nx)
  st,sx = Sampling(nt,dt,ft),Sampling(nx,dx,fx)
  ais = ArrayInputStream(dataDir+"tp73.dat")
  ais.readFloats(image)
  ais.close()
  return st,sx,image

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
  s1,s2,f = imageTeapot()
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
    if not gmin: gmin = min(g)
    if not gmax: gmax = max(g)
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
