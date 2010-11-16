#############################################################################
# Blended-neighbor gridding in 1D

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
from edu.mines.jtk.la import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fmm import *

#############################################################################
# global parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM


#############################################################################
# test functions

def main(args):
  goInterp()

def goInterp():
  fnull = -999.0
  s,f = gridFritschCarlton(fnull)
  #s,f = gridRamp(fnull)
  t,p = nearestNeighbor(f,fnull)
  q = blendedNeighbor(0.50,t,p)
  plot(fnull,s,f,p,q)
  #q = blendedNeighbor(0.50,t,q)
  #plot(fnull,s,f,p,q)
  #q = radialBasis(fnull,s,f)
  #plot(fnull,s,f,p,q)

def gridRamp(fnull):
  #x = [0.0,0.1,0.3,0.7,0.8,0.85,1.0]
  #y = [0.0,0.1,0.3,0.7,0.8,0.85,1.0]
  #x = [-1.0,0.0,1.0]
  #y = [-1.0,0.0,1.0]
  x = [0.0,1.0]
  y = [0.0,1.0]
  xmin,xmax = x[0],x[-1]
  fx = x[0]
  dx = 0.01
  nx = 1+int((xmax-xmin)/dx+0.5)
  sx = Sampling(nx,dx,fx)
  f = fillfloat(fnull,nx)
  n = len(x)
  for i in range(n):
    ix = sx.indexOfNearest(x[i])
    f[ix] = y[i]
    #if i!=5:
    #  f[ix] = 0.0
  return sx,f

def gridFritschCarlton(fnull):
  x = [7.99,8.09,8.19,8.7,9.2,10.0,12.0,15.0,20.0]
  y = [0.000000,2.76429e-5,4.37498e-2,0.169183,0.469428,
       0.943740,0.998636,0.999919,0.999994]
  xmin,xmax = x[0],x[-1]
  fx = xmin
  dx = 0.01
  nx = 1+int((xmax-xmin)/dx+0.5)
  sx = Sampling(nx,dx,fx)
  f = fillfloat(fnull,nx)
  n = len(x)
  for i in range(n):
    ix = sx.indexOfNearest(x[i])
    f[ix] = y[i]
    #if i!=5:
    #  f[ix] = 0.0
  return sx,f
 
def nearestNeighbor(f,fnull):
  n = len(f)
  t = zerofloat(n)
  p = zerofloat(n)
  for i in range(n):
    ilo,ihi = i,i
    while ilo>=0 and f[ilo]==fnull:
      ilo -= 1
    while ihi<n and f[ihi]==fnull:
      ihi += 1
    if ilo>=0 and (i-ilo<=ihi-i or ihi==n):
      t[i] = i-ilo
      p[i] = f[ilo]
    else:
      t[i] = ihi-i
      p[i] = f[ihi]
  return t,p

def blendedNeighborX(s,t,p):
  n = len(t)
  q = zerofloat(n)
  a = zerofloat(n)
  c = zerofloat(n)
  b = fillfloat(1.0,n)
  t = clip(0.0,1.0e6,t)
  stt = s*(t[0]*t[0]+t[1]*t[1])*0.5
  b[0] += stt; c[0] -= stt
  for i in range(1,n-1):
    a[i] -= stt; b[i] += stt
    stt = s*(t[i]*t[i]+t[i+1]*t[i+1])*0.5
    b[i] += stt; c[i] -= stt
  a[n-1] -= stt; b[n-1] += stt
  m = TridiagonalFMatrix(n,a,b,c)
  m.solve(p,q)
  return q

def blendedNeighbor(s,t,p):
  n = len(p)
  q = zerofloat(n)
  a = zerofloat(n)
  b = zerofloat(n)
  c = zerofloat(n)
  b[0] = 1.0+s*t[0]*t[0]
  c[0] =    -s*t[0]*t[1]
  for i in range(1,n-1):
    a[i] =        -s*t[i]*t[i-1]
    b[i] = 1.0+2.0*s*t[i]*t[i]
    c[i] =        -s*t[i]*t[i+1]
  a[n-1] =    -s*t[n-1]*t[n-2]
  b[n-1] = 1.0+s*t[n-1]*t[n-1]
  m = TridiagonalFMatrix(n,a,b,c)
  m.solve(p,q)
  return q

def radialBasis(fnull,s,f):
  def rbf(x):
    if x<0:
      return -x*x*x
    else:
      return x*x*x
  x,y = getKnownSamples(fnull,s,f)
  m = len(x)
  a = DMatrix(m,m)
  b = DMatrix(m,1)
  for j in range(m):
    for k in range(m):
      a.set(j,k,rbf(x[j]-x[k]))
    b.set(j,0,y[j])
  lud = DMatrixLud(a)
  v = lud.solve(b)
  w = zerofloat(m)
  for j in range(m):
    w[j] = v.get(j,0)
  n = s.count
  q = zerofloat(n)
  for i in range(n):
    xi = s.getValue(i)
    for j in range(m):
      q[i] += w[j]*rbf(xi-x[j])
  return q

def getKnownSamples(fnull,s,f):
  n = s.count
  m = 0
  for i in range(n):
    if f[i]!=fnull:
      m += 1
  x = zerofloat(m)
  y = zerofloat(m)
  m = 0
  for i in range(n):
    if f[i]!=fnull:
      x[m] = s.getValue(i)
      y[m] = f[i]
      m += 1
  return x,y

def plot(fnull,s,f,p,q):
  x,y = getKnownSamples(fnull,s,f)
  sp = SimplePlot()
  #sp.setVLimits(-0.1,1.1)
  pv = sp.addPoints(x,y)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  pv = sp.addPoints(s,p)
  pv.setLineStyle(PointsView.Line.DASH)
  pv = sp.addPoints(s,q)

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(fontSize)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(100,6,pngDir+"/"+png+".png")
  return frame

def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  return p

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
