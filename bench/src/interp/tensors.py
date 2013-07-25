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
from edu.mines.jtk.util.ArrayMath import *

def main(args):
  goTensors()
  #goGradients()

def goGradients():
  f,s1,s2 = readPnzImage()
  plotPnz(f,s1,s2,png="pnz")
  rgf = RecursiveGaussianFilter(1.0)
  g1 = copy(f); rgf.apply1X(f,g1); rgf.applyX0(g1,g1)
  g2 = copy(f); rgf.applyX1(f,g2); rgf.apply0X(g2,g2)
  plotPnz(g1,s1,s2,png="pnzeg1")
  plotPnz(g2,s1,s2,png="pnzeg2")
  g11 = mul(g1,g1)
  g12 = mul(g1,g2)
  g22 = mul(g2,g2)
  ggmax = max(max(g11),max(g12),max(g22))
  cmin,cmax = -0.1*ggmax,0.1*ggmax
  plotPnz(g11,s1,s2,cmin=cmin,cmax=cmax,png="pnzeg11")
  plotPnz(g12,s1,s2,cmin=cmin,cmax=cmax,png="pnzeg12")
  plotPnz(g22,s1,s2,cmin=cmin,cmax=cmax,png="pnzeg22")
  rgf = RecursiveGaussianFilter(4.0)
  rgf.apply00(g11,g11)
  rgf.apply00(g12,g12)
  rgf.apply00(g22,g22)
  ggmax = max(max(g11),max(g12),max(g22))
  cmin,cmax = -0.2*ggmax,0.2*ggmax
  plotPnz(g11,s1,s2,cmin=cmin,cmax=cmax,png="pnzes11")
  plotPnz(g12,s1,s2,cmin=cmin,cmax=cmax,png="pnzes12")
  plotPnz(g22,s1,s2,cmin=cmin,cmax=cmax,png="pnzes22")

def goTensors():
  #reads = [readPnzImage,readTpImage]
  #plots = [plotPnz,plotTp]
  #prefs = ["pnze","tpe"]
  reads = [readPnzImage]
  plots = [plotPnz]
  prefs = ["pnze"]
  for i,read in enumerate(reads):
    plot = plots[i]
    pref = prefs[i]
    g,s1,s2 = read()
    plot(g,s1,s2,png=pref)
    lof = LocalOrientFilter(4.0)
    s = lof.applyForTensors(g)
    d00 = EigenTensors2(s); d00.invertStructure(0.0,0.0)
    d01 = EigenTensors2(s); d01.invertStructure(0.0,1.0)
    d02 = EigenTensors2(s); d02.invertStructure(0.0,2.0)
    d04 = EigenTensors2(s); d04.invertStructure(0.0,4.0)
    d11 = EigenTensors2(s); d11.invertStructure(1.0,1.0)
    d12 = EigenTensors2(s); d12.invertStructure(1.0,2.0)
    d14 = EigenTensors2(s); d14.invertStructure(1.0,4.0)
    plot(g,s1,s2,png=pref)
    plot(g,s1,s2,d00,dscale=1,png=pref+"00")
    plot(g,s1,s2,d01,dscale=1,png=pref+"01")
    plot(g,s1,s2,d02,dscale=1,png=pref+"02")
    plot(g,s1,s2,d04,dscale=1,png=pref+"04")
    plot(g,s1,s2,d11,dscale=2,png=pref+"11")
    plot(g,s1,s2,d12,dscale=2,png=pref+"12")
    plot(g,s1,s2,d14,dscale=2,png=pref+"14")

#############################################################################
# plotting

#pngDir = "../../png/" # where to put PNG images of plots
pngDir = None # for no PNG images

backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent

def plotPnz(g,s1,s2,d=None,dscale=1,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  sp.setHInterval(2.0)
  sp.setVInterval(2.0)
  #sp.setFontSizeForPrint(8,240)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(710,750)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(1,99)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(20)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir and png:
    sp.paintToPng(360,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")

def plotTp(g,s1,s2,d=None,dscale=1,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  sp.setHInterval(2.0)
  sp.setVInterval(0.2)
  #sp.setFontSizeForPrint(8,240)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(910,670)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(1,99)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(20)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir and png:
    sp.paintToPng(360,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")

#############################################################################
# data input/output

def readPnzImage():
  s1 = Sampling(501,0.0125,0.0)
  s2 = Sampling(501,0.0125,0.0)
  g = readImage("/data/seis/pnz/ggs256",s1,s2)
  return g,s1,s2

def readTpImage():
  s1 = Sampling(251,0.004,0.500)
  s2 = Sampling(357,0.025,0.000)
  g = readImage("/Users/dhale/Home/git/jtk/data/tp73",s1,s2)
  return g,s1,s2

def readImage(fileName,s1,s2):
  n1,n2 = s1.count,s2.count
  ais = ArrayInputStream(fileName+".dat")
  x = zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(fileName+".dat")
  aos.writeFloats(x)
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
