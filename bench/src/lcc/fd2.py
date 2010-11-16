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

#############################################################################
# global parameters

dataDir = "/data/seis/tp/csm/oldslices/"
#pngDir = "./png"
pngDir = None

def main(args):
  #goShifts()
  goShifts2()
  #goFilter()

def goShifts2():
  s1,s2,f = imageSynth()
  n1,n2 = len(f[0]),len(f)
  u = zerofloat(n1,n2)
  c = zerofloat(n1,n2)
  d = zerofloat(n1,n2)
  plot2Teapot(s1,s2,f)
  for k in range(1,6):
    u,c,d = findShifts(30,-10,10,k,f)
    w = pow(c,2.0) 
    us = smooth(0.9,w,u)
    #plot2Teapot(s1,s2,f,c)
    #plot2Teapot(s1,s2,f,u)
    plot2Teapot(s1,s2,f,us)

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

def findShifts(sigma,min1,max1,lag2,f):
  n1,n2 = len(f[0]),len(f)
  u = zerofloat(n1,n2)
  c = zerofloat(n1,n2)
  d = zerofloat(n1,n2)
  lsf = LocalShiftFinder(sigma)
  for i2 in range(n2):
    i2m = max(i2-lag2,0)
    i2p = min(i2+lag2,n2-1)
    lsf.find1(min1,max1,f[i2m],f[i2p],u[i2],c[i2],d[i2])
  return u,c,d

def goFilter():
  #s1,s2,f = imageTeapot()
  s1,s2,f = imageSynth()
  n1,n2 = s1.count,s2.count
  a = 0.95
  fm = filterM(a,f)
  fp = filterP(a,f)
  plot2Teapot(s1,s2,f)
  plot2Teapot(s1,s2,fm)
  plot2Teapot(s1,s2,fp)

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

def shiftRamp(f):
  n = len(f)
  g = copy(f)
  t = rampfloat(0.0,1.0-8.0/(n-1),n)
  si = SincInterpolator()
  si.setUniform(n,1.0,0.0,f)
  si.interpolate(n,t,g)
  return g

def imageSynth():
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
  r = randomNoise(4.0,n1,n2)
  g = add(g,r)
  rgf = RecursiveGaussianFilter(2.0)
  rgf.applyX0(g,g)
  return s1,s2,g

def randomNoise(a,n1,n2):
  r = mul(2.0*a,sub(randfloat(n1,n2),0.5))
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply00(r,r)
  return r

def goShifts():
  #s1,s2,f = imageTeapot()
  s1,s2,f = imageSynth()
  plot2Teapot(s1,s2,f)
  n1,n2 = s1.count,s2.count
  u = zerofloat(n1,n2)
  c = zerofloat(n1,n2)
  d = zerofloat(n1,n2)
  sigma = 30
  min1 = -15
  max1 =  15
  lsf = LocalShiftFinder(sigma)
  for i2 in range(1,n2-1):
    lsf.find1(min1,max1,f[i2],f[i2+1],u[i2],c[i2],d[i2])
  #u = mul(0.5,u)
  copy(u[1],u[0])
  copy(c[1],c[0])
  copy(d[1],d[0])
  copy(u[n2-2],u[n2-1])
  copy(c[n2-2],c[n2-1])
  copy(d[n2-2],d[n2-1])
  #plot2Teapot(s1,s2,f,c)
  #plot2Teapot(s1,s2,f,d)
  w = pow(d,2.0)
  plot2Teapot(s1,s2,f,w)
  plot2Teapot(s1,s2,f,u)
  #ul,vl = shiftsFilterLR(0.5,w,u)
  #ur,vr = shiftsFilterRL(0.5,w,u)
  #plot2Teapot(s1,s2,f,ul)
  #plot2Teapot(s1,s2,f,ur)
  us = shiftsFilter(0.5,0.9,w,u)
  print "us: min =",min(us)," max =",max(us)
  plot2Teapot(s1,s2,f,us)
  #ul,vl = shiftsFilterLR(0.9,w,u)
  #plot2Teapot(s1,s2,f,ul)
  #plot2Teapot(s1,s2,f,sqrt(vl))

"""
x x x . x x x
"""
def shiftsFilter(t,a,w,u):
  ul,vl = shiftsFilterLR(a,w,u)
  ur,vr = shiftsFilterRL(a,w,u)
  n1,n2 = len(u[0]),len(u)
  tt = t*t
  us = zerofloat(n1,n2)
  for i2 in range(1,n2-1):
    for i1 in range(n1):
      dl = u[i2][i1]-ul[i2-1][i1]
      dr = u[i2][i1]-ur[i2+1][i1]
      if dl*dl>tt*vl[i2][i1] and dr*dr>tt*vr[i2][i1]:
        if dl*dl>0.04 and dr*dr>0.04:
          us[i2][i1] = u[i2][i1]
  return us

def shiftsFilterLR(a,w,u):
  b = 1.0-a
  n1,n2 = len(u[0]),len(u)
  ws = copy(w[0]) # running sum of weights
  us = mul(w[0],u[0]) # running sum numerator for mean shift
  vs = zerofloat(n1) # running sum numerator for variance
  uu = div(us,ws) # mean shift = us/ws
  vv = div(vs,ws) # variance = vs/ws, initially zero
  ul = zerofloat(n1,n2) # shift smoothed
  vl = zerofloat(n1,n2) # variance smoothed
  for i2 in range(1,n2):
    for i1 in range(n1):
      wi = w[i2][i1]
      ws[i1] = a*ws[i1]+b*wi
      ui = u[i2][i1]
      us[i1] = a*us[i1]+b*wi*ui
      uu[i1] = us[i1]/ws[i1]
      vi = (ui-uu[i1])*(ui-uu[i1])
      vs[i1] = a*vs[i1]+b*wi*vi
      vv[i1] = vs[i1]/ws[i1]
      ul[i2][i1] = uu[i1]
      vl[i2][i1] = vv[i1]
  return ul,vl

def shiftsFilterRL(a,w,u):
  b = 1.0-a
  n1,n2 = len(u[0]),len(u)
  ws = copy(w[n2-1]) # running sum of weights
  us = mul(w[n2-1],u[n2-1]) # running sum numerator for mean shift
  vs = zerofloat(n1) # running sum numerator for variance
  uu = div(us,ws) # mean shift = us/ws
  vv = div(vs,ws) # variance = vs/ws, initially zero
  ur = zerofloat(n1,n2) # shift smoothed
  vr = zerofloat(n1,n2) # variance smoothed
  for i2 in range(n2-2,-1,-1):
    for i1 in range(n1):
      wi = w[i2][i1]
      ws[i1] = a*ws[i1]+b*wi
      ui = u[i2][i1]
      us[i1] = a*us[i1]+b*wi*ui
      uu[i1] = us[i1]/ws[i1]
      vi = (ui-uu[i1])*(ui-uu[i1])
      vs[i1] = a*vs[i1]+b*wi*vi
      vv[i1] = vs[i1]/ws[i1]
      ur[i2][i1] = uu[i1]
      vr[i2][i1] = vv[i1]
  return ur,vr
  
def medianFilter(m,w,u):
  n1,n2 = len(u[0]),len(u)
  v = zerofloat(n1,n2)
  m2p1 = m*2+1
  mf = MedianFinder(m2p1)
  ww = zerofloat(m2p1)
  uu = zerofloat(m2p1)
  for i2 in range(n2):
    i2m = max(i2-m,0)
    i2p = min(i2m+m+m,n2-1)
    i2m = min(i2m,i2p-m-m)
    for i1 in range(n1):
      k2 = 0
      for j2 in range(i2m,i2p+1):
        ww[k2] = w[j2][i1]
        uu[k2] = u[j2][i1]
        k2 += 1
      v[i2][i1] = mf.findMedian(ww,uu)
  return v
 
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
 
#############################################################################
# plotting

def plot2Teapot(s1,s2,f,g=None,label=None,png=None):
  n1 = len(f[0])
  n2 = len(f)
  panel = panel2Teapot()
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(180)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-4.5,4.5)
  if g:
    alpha = 0.5
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    gclip = max(abs(min(g)),abs(max(g)))
    pv.setClips(-gclip,gclip)
    pv.setColorModel(ColorMap.getJet(alpha))
    #updateColorModel(pv,alpha)
  frame2Teapot(panel,png)

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
      r[n/2] = r[n/2-1] = -1
      g[n/2] = g[n/2-1] = -1
      b[n/2] = b[n/2-1] = -1
      a[n/2] = a[n/2-1] = 0
    icm = IndexColorModel(8,n,r,g,b,a)
    pv.setColorModel(icm)

def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  return panel

def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(1290,777)
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
