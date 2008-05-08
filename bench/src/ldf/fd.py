import sys
from math import *
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

from ldf import *

True = 1
False = 0

#############################################################################
# parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

n1 = 315
n2 = 315
aniso = 100
small = 0.00001
niter = 1000
lof = LocalOrientFilter(8)
lof.setGradientSmoothing(1)

#############################################################################
# functions

def main(args):
  goStencils()
  #goGaussian()
  #goTestAni()
  return

def frequencyResponseFd(rs,ldt):
  x = Array.zerofloat(3,3)
  x[1][1] = 1.0
  y = Array.zerofloat(3,3)
  ldk = LocalDiffusionKernel(rs)
  ldk.apply(ldt,x,y)
  return frequencyResponse(y)

def frequencyResponseEx(ldt):
  d = Array.zerofloat(3)
  ldt.getTensor(0,0,d)
  d11,d12,d22 = d[0],d[1],d[2]
  nk1,nk2 = 315,315
  dk1,dk2 = 2*pi/nk1,2*pi/nk2
  fk1,fk2 = -pi+dk1/2,-pi+dk2/2
  a = Array.zerofloat(nk1,nk2)
  for ik2 in range(nk2):
    k2 = fk2+ik2*dk2
    for ik1 in range(nk1):
      k1 = fk1+ik1*dk1
      a[ik2][ik1] = d11*k1*k1+2.0*d12*k1*k2+d22*k2*k2
  return a

def errorSample(e):
  n1,n2 = len(e[0]),len(e)
  k1,k2 = (n1-1)/2,(n2-1)/2
  mj = 1
  for jj in [[29,0],[21,20],[0,29],[-20,21]]:
    j1,j2 = k1+mj*jj[0],k2+mj*jj[1]
    ej = abs(e[j2][j1])
    print "e =",ej

def testAni(rs,v1):
  v1 = Array.fillfloat(v1,3,3)
  ldt = LocalDiffusionTensors2(0.0,1.0,None,None,v1)
  fa = frequencyResponseFd(rs,ldt)
  fb = frequencyResponseEx(ldt)
  plotf(fa)
  #plotf(fb)
  e = Array.sub(fb,fa)
  #plotf(e)
  errorSample(e)

def goTestAni():
  for v1 in [1.0,1/sqrt(2.0)]:
    for rs in [0,1/12.0]:
      print "v1 =",v1,"rs =",rs
      testAni(rs,v1)

def goStencils():
  rs0 = 0.0/12.0
  rs1 = 1.0/12.0
  rs2 = 2.0/12.0
  rs3 = 3.0/12.0
  iso = False
  for theta in [90,-45,45]:
    for rs in [rs0,rs1,rs2,rs3]:
      #stencilOne(iso,rs,theta)
      stencilAvg(iso,rs,theta)
  iso = True
  theta = 0
  #for rs in [rs0,rs1,rs2,rs3]:
  #  stencilAvg(iso,rs,theta)

def stencilOne(iso,rs,theta):
  print "one: iso =",iso,"rs =",rs,"theta =",theta
  png = "one"+str(iso)+str(int(rs*12))+str(theta)
  if iso:
    a = 1.0
    b = 0.0
    c = 1.0
  else:
    theta *= pi/180.0
    v1 = sin(theta)
    v2 = cos(theta)
    a = v1*v1
    b = v1*v2
    c = v2*v2
  r = (1+sqrt(1-4*rs))/2
  s = (1-sqrt(1-4*rs))/2
  g1 = [[-s,s],[-r,r]]
  g2 = [[-s,-r],[s,r]]
  g1g1 = Array.zerofloat(3,3)
  g1g2 = Array.zerofloat(3,3)
  g2g1 = Array.zerofloat(3,3)
  g2g2 = Array.zerofloat(3,3)
  h = Array.zerofloat(3,3)
  Conv.xcor(2,2,0,0,g1,2,2,0,0,g1,3,3,-1,-1,g1g1)
  Conv.xcor(2,2,0,0,g1,2,2,0,0,g2,3,3,-1,-1,g1g2)
  Conv.xcor(2,2,0,0,g2,2,2,0,0,g1,3,3,-1,-1,g2g1)
  Conv.xcor(2,2,0,0,g2,2,2,0,0,g2,3,3,-1,-1,g2g2)
  for i2 in range(3):
    for i1 in range(3):
      h[i2][i1] = a*g1g1[i2][i1]+b*(g1g2[i2][i1]+g2g1[i2][i1])+c*g2g2[i2][i1]
  Array.dump(Array.transpose(h))
  plotf(frequencyResponse(h),png)

def stencilAvg(iso,rs,theta):
  print "avg: iso =",iso,"rs =",rs,"theta =",theta
  png = "avg"+str(iso)+str(int(rs*12))+str(theta)
  theta *= pi/180.0
  n1,n2 = 3,3
  x = Array.zerofloat(n1,n2)
  y = Array.zerofloat(n1,n2);
  x[n2/2][n1/2] = 1.0
  v1 = Array.fillfloat(sin(theta),n1,n2)
  if iso:
    ldt = LocalDiffusionTensors2(1.0,0.0,None,None,v1)
  else:
    ldt = LocalDiffusionTensors2(0.0,1.0,None,None,v1)
  ldk = LocalDiffusionKernel(rs)
  ldk.apply(ldt,x,y)
  Array.dump(Array.transpose(y))
  plotf(frequencyResponse(y),png)

def gaussianResponse(nx,sigma,rs,ldt):
  t = Array.zerofloat(3)
  ldt.getTensor(0,0,t)
  a,b,c = t[0],t[1],t[2]
  dx = 8.0*sigma/(nx-1)
  fx = -(nx-1)/2.0*dx
  f = Array.zerofloat(nx,nx)
  g = Array.zerofloat(nx,nx)
  h = Array.zerofloat(nx,nx)
  oss = 1.0/(sigma*sigma)
  for i2 in range(nx):
    x2 = fx+i2*dx
    for i1 in range(nx):
      x1 = fx+i1*dx
      g00 = exp(-0.5*oss*(x1*x1+x2*x2))
      g20 = oss*(1.0-x1*x1*oss)*g00
      g11 =     (   -x1*x2*oss)*g00
      g02 = oss*(1.0-x2*x2*oss)*g00
      f[i2][i1] = g00
      g[i2][i1] = a*g20+2.0*b*g11+c*g02
  ldk = LocalDiffusionKernel(rs)
  ldk.apply(ldt,f,h)
  Array.mul(1.0/(dx*dx),h,h)
  e = Array.sub(h,g)
  for i2 in range(nx):
    for i1 in range(nx):
      if i1==0 or i1==nx-1 or i2==0 or i2==nx-1:
        e[i2][i1] = 0.0
  #print "g min =",Array.min(g)," max =",Array.max(g)
  #print "h min =",Array.min(h)," max =",Array.max(h)
  #print "e abs =",Array.max(Array.abs(e))
  print "e rms =",rmsError(h,g)
  #plot(g,0,0,jet)
  #plot(h,0,0,jet)
  plot(e,0,0,jet)
def goGaussian():
  rs0 = 0.0/12.0
  rs1 = 1.0/12.0
  rs3 = 3.0/12.0
  nx = 129
  sigma = 1.0
  theta = pi/4
  for rs in [rs0,rs1,rs3]:
    print "rs =",rs
    for theta in [0*pi/16,1*pi/16,2*pi/16,3*pi/16,4*pi/16]:
      v1 = Array.fillfloat(sin(theta),nx,nx)
      ldt = LocalDiffusionTensors2(0.0,1.0,None,None,v1)
      gaussianResponse(nx,sigma,rs,ldt)

def rmsError(a,b):
  n1,n2 = len(a[0]),len(a)
  s = 0.0
  for i2 in range(1,n2-1):
    for i1 in range(1,n1-1):
      d = a[i2][i1]-b[i2][i1]
      s += d*d
  return sqrt(s/(n2-2)/(n1-2))

def frequencyResponse(x):
  n1,n2 = 315,315
  n1 = FftComplex.nfftSmall(n1)
  n2 = FftComplex.nfftSmall(n2)
  xr = Array.zerofloat(n1,n2)
  xi = Array.zerofloat(n1,n2)
  Array.copy(len(x[0]),len(x),x,xr)
  cx = Array.cmplx(xr,xi)
  fft1 = FftComplex(n1)
  fft2 = FftComplex(n2)
  fft1.complexToComplex1(1,n2,cx,cx)
  fft2.complexToComplex2(1,n1,cx,cx)
  ax = Array.cabs(cx)
  a = Array.zerofloat(n1,n2)
  j1 = n1/2
  j2 = n2/2
  Array.copy(n1-j1,n2-j2,0,0,ax,j1,j2,a)
  Array.copy(j1,j2,n1-j1,n2-j2,ax,0,0,a)
  Array.copy(n1-j1,j2,0,n2-j2,ax,j1,0,a)
  Array.copy(j1,n2-j2,n1-j1,0,ax,0,j2,a)
  return a
 
#############################################################################
# Pythagorean triples
"""
( 3, 4, 5)	( 5, 12, 13)	( 7, 24, 25)	( 8, 15, 17)
( 9, 40, 41)	(11, 60, 61)	(12, 35, 37)	(13, 84, 85)
(16, 63, 65)	(20, 21, 29)	(28, 45, 53)	(33, 56, 65)
(36, 77, 85)	(39, 80, 89)	(48, 55, 73)	(65, 72, 97)
"""
 
#############################################################################
# plot

def plot(f,cmin=0.0,cmax=0.0,cmap=ColorMap.GRAY,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(0.0,100.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  frame(p,png)

def plot2(f,g,cmin=0,cmax=0,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  pv = p.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.getGray())
  pv.setClips(-10,10)
  pv = p.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin!=cmax:
    pv.setClips(cmin,cmax)
  pv.setColorModel(ColorMap.getJet(0.3))
  frame(p,png)

def plotf(f,png=None):
  n1 = len(f[0])
  n2 = len(f)
  p = panel()
  pv = p.addPixels(f)
  pv.setClips(0.0,4.0)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(p,png)

fontSize = 24
width = 600
height = 622
widthColorBar = 80

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  #p = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.LEFT_TOP)
  #p.addColorBar()
  #p.setColorBarWidthMinimum(widthColorBar)
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
