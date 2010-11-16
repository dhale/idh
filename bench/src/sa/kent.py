# Process data from Kent Broadhead
info = """
(from datafortest.m, provided by Kent)
load final krig xsoft ysoft xwell ywell zwell A0 cokrig 

figure
subplot(2,2,1),imagesc(xsoft,ysoft,krig,[.13,.2]),set(gca,'ydir','n'),colorbar,title('Kriged Wells')
hold on,plot(xwell,ywell,'+')
subplot(2,2,2),imagesc(xsoft,ysoft,A0),set(gca,'ydir','n'),colorbar,title('Soft Data')
subplot(2,2,4),imagesc(xsoft,ysoft,cokrig,[.13,.2]),set(gca,'ydir','n'),colorbar,title('Cokriging')
hold on,plot(xwell,ywell,'+')

% Dave, the above is for display purpose so that you can correlate 
back to the Powerpoint I gave you.
%The only actual data you need is the well points (xwell, ywell, zwell), 
the soft data A0 (attribute from seismic) and the soft data xy's 
(xsoft, ysoft). 
%Once you compute your interpolation, pass it back to me (I can 
read a segy) and I will incoporte it into a data integration with 
amplitudes (like cokriging) and make various comparison figures 
again kriging, cokriging, etc. 
"""

from math import *
import sys
from math import *
from java.awt import *
from java.lang import *
from java.io import *
from java.util import *
from java.nio import *
from javax.swing import *

from com.jmatio.io import *
from com.jmatio.types import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

dataDir = "/data/sa/kent/"
nx,dx,fx = 50,450.000,195000.0
ny,dy,fy = 50,469.388,2758500.0 # using average dy, real dy = 400 or 500
n1,d1,f1 = 51,1.0,0.0
n2,d2,f2 = 51,1.0,0.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)

def main(args):
  xw,yw,zw,s = readMatlab()
  plot(xw,yw,s,"Amplitude","Seismic image")
  p = grid(xw,yw,zw)
  plot(xw,yw,p,"Porosity","Image-ignorant")
  writeImage(dataDir+"pii.txt",p)
  p = grid(xw,yw,zw,s)
  plot(xw,yw,p,"Porosity","Image-guided")
  writeImage(dataDir+"pig.txt",p)

def grid(x,y,f,s=None):
  bg = BlendedGridder2(f,x,y)
  bg.setSmoothness(0.7)
  if s:
    d = makeImageTensors(s)
    bg.setTensors(d)
  q = bg.grid(s1,s2)
  return q

def plot(x,y,s,slabel="",title=""):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setTitle(title)
  sp.setHLabel("x")
  sp.setVLabel("y")
  sp.plotPanel.setColorBarWidthMinimum(100)
  sp.setSize(934,830)
  sp.addColorBar(slabel)
  pv = sp.addPixels(s)
  if slabel=="Porosity":
    pv.setClips(0.13,0.20)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if x and y:
    pv = sp.addPoints(x,y)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setLineWidth(3)
    pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
    pv.setMarkColor(Color.WHITE)

def readMatlab():
  mfr = MatFileReader(dataDir+"final.mat")
  mfh = mfr.getMatFileHeader()
  print mfh
  mc = mfr.getContent()
  for e in mc:
    print e
  xw = matfloat1(mc.get("xwell").getArray())
  yw = matfloat1(mc.get("ywell").getArray())
  zw = matfloat1(mc.get("zwell").getArray())
  xw = mul(1.0/dx,sub(xw,fx))
  yw = mul(1.0/dy,sub(yw,fy))
  s = matfloat2(mc.get("A0").getArray())
  #s = matfloat2(mc.get("krig").getArray())
  """
  xs = mc.get("xsoft").getArray()
  print "fx =",xs[0][0]
  print "dx =",xs[ 50][0]-xs[  0][0]
  print "dx =",xs[100][0]-xs[ 50][0]
  print "dx =",xs[150][0]-xs[100][0]
  ys = mc.get("ysoft").getArray()
  print "fy =",ys[0][0]
  dela = 0.0
  for iy in range(1,ny):
    deli = ys[iy][0]-ys[iy-1][0]
    dela += deli
    print "dy =",deli
  dela /= ny-1
  print "dela =",dela
  """
  return xw,yw,zw,s

def writeImage(fileName,f):
  fos = FileOutputStream(fileName)
  pw = PrintWriter(fos)
  for iy in range(ny):
    for ix in range(nx):
      pw.println(f[iy][ix])
  pw.close()

def matfloat1(x):
  x = flatten(x)
  n = len(x)
  y = zerofloat(n)
  for i in range(n):
    y[i] = x[i]
  return y

def matfloat2(f):
  g = zerofloat(n1,n2)
  for ix in range(nx):
    for iy in range(ny):
      g[ix][iy] = f[ix][iy]
  for i2 in range(n2):
    g[i2][n1-1] = g[i2][n1-2]
  for i1 in range(n1):
    g[n2-1][i1] = g[n2-2][i1]
  return g

def makeImageTensorsX(s):
  """ 
  Returns guiding tensors interpolation along features in specified image.
  """
  n1,n2 = len(s[0]),len(s)
  lof = LocalOrientFilter(6)
  t = lof.applyForTensors(s)
  su = zerofloat(n1,n2)
  sv = zerofloat(n1,n2)
  t.getEigenvalues(su,sv)
  e1 = div(sub(su,sv),su)
  e2 = div(sv,su)
  e2 = pow(e2,4)
  eu = e2
  ev = add(e1,e2)
  t.setEigenvalues(eu,ev)
  plot(None,None,e2,"Isotropy","Isotropy")
  return t

def makeImageTensors(s):
  """ 
  Returns guiding tensors interpolation along features in specified image.
  """
  n1,n2 = len(s[0]),len(s)
  lof = LocalOrientFilter(6)
  t = lof.applyForTensors(s)
  lsf = LocalSemblanceFilter(8,8)
  eu = lsf.semblance(LocalSemblanceFilter.Direction2.UV,t,s)
  ev = lsf.semblance(LocalSemblanceFilter.Direction2.V ,t,s)
  #eu = pow(eu,2.0)
  eu = clip(0.0001,1.0,eu)
  ev = clip(0.0001,1.0,ev)
  t.setEigenvalues(eu,ev)
  plot(None,None,eu,"Eigenvalue","Eigenvalue eu")
  plot(None,None,ev,"Eigenvalue","Eigenvalue ev")
  return t

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

