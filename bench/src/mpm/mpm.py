#############################################################################
# Smoothing of seismicity data provided by Morgan P. Moschetti.

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

dataDir = "/data/earth/mpm/"
pngDir = "./"
#pngDir = None

#############################################################################

def main(args):
  goSmooth()
  #goDisplay()

def goDisplay():
  sx,sy,s = readSeis("vals_10a_ceus.txt")
  plotSeis(sx,sy,s)

def goSmooth():
  sx,sy,s = readSeis("vals_10a_ceus.txt")
  t = getTensors(s)
  r = smoothAni(5.0,t,s)
  plotSeis(sx,sy,s,t,clip=12.0,png="s")
  plotSeis(sx,sy,r,t,clip=2.0,png="r")

#############################################################################
# Smoothing

def smoothIso(sigma,s):
  rgf = RecursiveGaussianFilter(sigma)
  r = copy(s)
  rgf.apply00(s,r)
  return r

def smoothAni(sigma,t,s):
  r = copy(s)
  c = 0.5*sigma*sigma
  lsf = LocalSmoothingFilter()
  lsf.apply(t,c,s,r)
  lsf.applySmoothS(r,r)
  return r

def getTensors(s):
  p0,p1 = 0.1,4.0
  sigmaS = 8.0
  sigmaG = sigmaS/4.0
  nx,ny = len(s[0]),len(s)
  lof = LocalOrientFilter(sigmaS)
  lof.setGradientSmoothing(sigmaG)
  s = add(s,mul(0.01,randfloat(nx,ny))) # tiny numbers instead of zeros
  t = lof.applyForTensors(s)
  t.invertStructure(p0,p1)
  return t

#############################################################################
# Plotting

def plotSeis(sx,sy,s,t=None,clip=None,png=None):
  sp = SimplePlot()
  sp.setSize(1200,700)
  sp.setVLabel("Latitude")
  sp.setHLabel("Longitude")
  sp.addColorBar("Seismicity")
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(sx,sy,s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  if clip:
    pv.setClips(0.0,clip)
  if t:
    tv = TensorsView(sx,sy,t)
    tv.setLineColor(Color.WHITE)
    sp.add(tv)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,3.3,pngDir+png+".png")

#############################################################################
# Input/output

def readSeis(fileName):
  nx,ny = 501,255
  dx,dy = 0.1,0.1
  fx,fy = -115.0,24.6
  sx,sy = Sampling(nx,dx,fx),Sampling(ny,dy,fy)
  scanner = Scanner(FileInputStream(dataDir+fileName))
  s = zerofloat(ny,nx)
  for ix in range(nx):
    for iy in range(ny):
      s[ix][iy] = scanner.nextFloat()
  s = transpose(s) # make x the fast dimension
  scanner.close()
  return sx,sy,s

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
