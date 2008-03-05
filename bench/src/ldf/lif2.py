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
pngDir = "./png"
#pngDir = None

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
  #goInterp()
  #goStencils()
  goGaussian()
  return

def goInterp():
  x = readImage()
  #x = flip2(x)
  n1,n2 = len(x[0]),len(x)
  v1,v2,ds = getV(x)
  #theta = 0.15*pi
  #v1 = Array.fillfloat(sin(theta),n1,n2)
  #v2 = Array.fillfloat(cos(theta),n1,n2)
  f = Array.zerobyte(n1,n2)
  y = Array.zerofloat(n1,n2)
  z = Array.zerofloat(n1,n2)
  #plot(x,-10,10,gray,"x")
  #for itest in [0,1,2,3]:
  for itest in [0,1,2,3]:
    Array.zero(f);
    Array.zero(y);
    Array.zero(z);
    if itest<2:
      d0 = None
      d1 = None
      i2s = [n2/2]
    else:
      d0 = makeFault(n1,n2)
      d1 = makeFault(n1,n2)
      i2s = [5,n2/2]
    if itest==0 or itest==2:
      for i2 in i2s:
        for i1 in range(n1):
          f[i2][i1] = 1
          y[i2][i1] = float(i1)/n1
      cmin =  0.0
      cmax =  1.0
      ldt = LocalDiffusionTensors2(0.1,1.0,d0,d1,v1)
    else:
      for i2 in i2s:
        for i1 in range(n1):
          f[i2][i1] = 1
        #y[i2] = smoothSgn(x[i2])
        y[i2] = Array.copy(x[i2])
      cmin = -10
      cmax = 10
      ldt = LocalDiffusionTensors2(0.01,1.0,d0,d1,v1)
    Array.copy(y,z)
    small = 0.001
    niter = 1000
    lif = LocalInterpolationFilterIc(small,niter)
    lif.apply(ldt,f,z)
    #plot(y,cmin,cmax,jet,"y"+str(itest))
    plot(z,cmin,cmax,jet,"z"+str(itest))
    plot2(x,z,cmin,cmax,"xz"+str(itest))

def smoothSgn(x):
  n = len(x)
  y = Array.zerofloat(n)
  for i in range(n):
    if x[i]<0:
      y[i] = -10
    else:
      y[i] = 10
  lgf = RecursiveGaussianFilter(2.0)
  lgf.apply0(y,y)
  return y

def goInterpOld():
  x = readImage()
  n1,n2 = len(x[0]),len(x)
  v1,v2,ds = getV(x)
  #theta = 0.1*pi
  #v1 = Array.fillfloat(sin(theta),n1,n2)
  #v2 = Array.fillfloat(cos(theta),n1,n2)
  #ds = Array.fillfloat(1.0,n1,n2)
  f = Array.zerobyte(n1,n2)
  y = Array.zerofloat(n1,n2)
  z = Array.zerofloat(n1,n2)
  #plot(x,-10,10,gray,"x")
  #for itest in [0,1,2,3]:
  for itest in [0]:
    if itest==0 or itest==2:
      for i2 in [n2/2]:
        for i1 in range(n1):
          f[i2][i1] = 1
          y[i2][i1] = sin(2*pi*i1/n1) 
          #y[i2][i1] = 2.0*i1/n1-1.0
      cmin = -1.0
      cmax =  1.0
      sigma = 0.1
      aniso = 1000
    else:
      for i2 in [n2/2]:
        for i1 in range(n1):
          f[i2][i1] = 1
          y[i2][i1] = x[i2][i1]
      cmin = -10
      cmax = 10
      sigma = 1
      aniso = 100
    lif = LocalInterpolationFilter(aniso,small,niter)
    if itest<2:
      ds = None
      es = None
    else:
      ds = makeBlock(n1,n2)
      es = makeBlock(n1,n2)
    z1 = Array.copy(y)
    z2 = Array.copy(y)
    lif.apply(ds,v1,f,z1)
    lif.applyLinear(ds,es,v1,f,z2)
    z2[0][n1-1] = 1.0
    #lif.applyLinear(sigma,aniso,ds,es,v1,f,z)
    #plot(y,cmin,cmax,jet,"y"+str(itest))
    plot(z1,cmin,cmax,jet,"z1"+str(itest))
    plot(z2,cmin,cmax,jet,"z2"+str(itest))
    #plot(Array.sub(z,y),cmin,cmax,jet)
    #plot2(x,z,cmin,cmax,"xz"+str(itest))

def readImage():
  fileName = dataDir+"/seis/vg/junks.dat"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  f = Array.zerofloat(n1,n2)
  ais.readFloats(f)
  ais.close()
  return f

def getV(x):
  n1 = len(x[0])
  n2 = len(x)
  v1 = Array.zerofloat(n1,n2)
  v2 = Array.zerofloat(n1,n2)
  el = Array.zerofloat(n1,n2)
  lof.apply(x,None,None,None,v1,v2,None,None,el)
  return v1,v2,el

def makeFault(n1,n2):
  ds = Array.fillfloat(1.0,n1,n2);
  #for i2 in range(n2/3,n2):
  #  for i1 in range(n1):
  #    ds[i2][i1] = 1.0
  for i1 in range(n1):
    i2 = 100-(i1*50)/n1
    ds[i2  ][i1] = 0.001
    ds[i2+1][i1] = 0.001
  return ds

def makeBlock(n1,n2):
  ds = Array.fillfloat(0.001,n1,n2);
  #for i2 in range(n2/3,n2):
  #  for i1 in range(n1):
  #    ds[i2][i1] = 1.0
  for i1 in range(n1):
    for i2 in range(100-(i1*50)/315,n2):
      ds[i2][i1] = 1.0
  return ds

def flip2(f):
  n1 = len(f[0])
  n2 = len(f)
  g = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    Array.copy(f[n2-1-i2],g[i2])
  return g

def goStencils():
  rs0 = 0.0/12.0
  rs1 = 1.0/12.0
  rs3 = 3.0/12.0
  iso = False
  for theta in [90,-45,45]:
    for rs in [rs0,rs1,rs3]:
      stencilOne(iso,rs,theta)
      stencilAvg(iso,rs,theta)
  iso = True
  theta = 0
  for rs in [rs0,rs1,rs3]:
    stencilAvg(iso,rs,theta)

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
  print "e abs =",Array.max(Array.abs(e))
  #plot(g,0,0,jet)
  #plot(h,0,0,jet)
  plot(e,0,0,jet)
def goGaussian():
  rs0 = 0.0/12.0
  rs1 = 1.0/12.0
  rs3 = 3.0/12.0
  nx = 65
  sigma = 1.0
  theta = pi/4
  for rs in [rs0,rs1,rs3]:
    for theta in [0*pi/16,1*pi/16,2*pi/16,3*pi/16,4*pi/16]:
      v1 = Array.fillfloat(sin(theta),nx,nx)
      ldt = LocalDiffusionTensors2(0.0,1.0,None,None,v1)
      gaussianResponse(nx,sigma,rs,ldt)

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
#width = 500
#height = 520
width = 700
height = 680
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
