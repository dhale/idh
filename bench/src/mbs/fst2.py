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

from fault import FaultScanner2,FaultSemblance
from dnp import LocalSlopeFinder

#############################################################################

#pngDir = "./png"
pngDir = None

sigmaTheta = 40
#smoother = FaultScanner2.Smoother.FFT
smoother = FaultScanner2.Smoother.SHEAR

def main(args):
  #goSlopes()
  #goAlign()
  #goSemblance()
  #goScan()
  #goThin()
  goShifts()

def getImage():
  return imageMbs()

def goShifts():
  s1,s2,g = getImage()
  g = slog(mul(2.0,g))
  plot2(s1,s2,g,title="log input")
  fse = FaultSemblance()
  g = fse.taper(10,g)
  p = fse.slopes(g)
  #p = getSlopes(g)
  sn,sd = fse.semblanceNumDen(p,g)
  fsc = FaultScanner2(sigmaTheta,[sn,sd],smoother)
  f,t = fsc.scan(-25,25)
  plot2(s1,s2,g,f,gmin=0,gmax=1,title="fault likelihood")
  #ff,tt = fsc.thin([f,t])
  #plot2(s1,s2,g,ff,gmin=0,gmax=1,title="fault likelihood")
  shiftMin,shiftMax = -20,20
  faults = fsc.findFaults([f,t],shiftMax-shiftMin);
  ff = faults.getLikelihoods()
  plot2(s1,s2,g,ff,gmin=0,gmax=1,title="fault likelihood")
  g = fsc.smooth(16,p,ff,g)
  plot2(s1,s2,g,ff,gmin=0,gmax=1,title="input smoothed")
  p = fse.slopes(g)
  #p = getSlopes(g)
  faults.findShifts(g,p,shiftMin,shiftMax)
  faults.clean()
  s = faults.getShifts()
  print "s min =",min(s)," max =",max(s)
  plot2(s1,s2,g,s,gmin=-12,gmax=12,title="fault throws")

def goThin():
  s1,s2,g = getImage()
  plot2(s1,s2,g,title="input")
  g = slog(mul(2.0,g))
  plot2(s1,s2,g,title="log input")
  fse = FaultSemblance()
  g = fse.taper(10,g)
  for iter in range(1):
    #p = fse.slopes(g)
    p = getSlopes(g)
    #p = zerofloat(len(p[0]),len(p))
    sn,sd = fse.semblanceNumDen(p,g)
    fsc = FaultScanner2(sigmaTheta,[sn,sd],smoother)
    f,t = fsc.scan(-25,25)
    plot2(s1,s2,g,f,gmin=0,gmax=1,title="fault likelihood")
    #plot2(s1,s2,g,t,title="fault dip (degrees)")
    ft,tt = fsc.thin([f,t])
    plot2(s1,s2,g,ft,gmin=0,gmax=1,title="fault likelihood thinned")
    #plot2(s1,s2,g,tt,title="fault dip (degrees) thinned")
    g = fsc.smooth(16,p,ft,g)
    plot2(s1,s2,g,title="input smoothed")

def goScan():
  s1,s2,g = getImage()
  g = slog(g)
  fse = FaultSemblance()
  g = fse.taper(10,g)
  #p = fse.slopes(g)
  p = getSlopes(g)
  sn,sd = fse.semblanceNumDen(p,g)
  fsc = FaultScanner2(sigmaTheta,[sn,sd],FaultScanner2.Smoother.FFT)
  st = Sampling(25,2.0,-25.0)
  for theta in st.values:
    f = fsc.likelihood(theta)
    plot2(s1,s2,g,f,gmin=0.2,gmax=0.7,gmap=jetr,
      title="theta = "+str(int(theta)))
  tmin,tmax = st.first,st.last
  f,t = fsc.scan(tmin,tmax)
  plot2(s1,s2,g,f,gmin=0.2,gmax=0.7,gmap=jetr,title="fault likelihood")
  #plot2(s1,s2,g,t,gmin=tmin,gmax=tmax,title="fault dip (degrees)")

def goSemblance():
  s1,s2,g = getImage()
  g = slog(g)
  fse = FaultSemblance()
  g = fse.taper(10,g)
  #p = fse.slopes(g)
  p = getSlopes(g)
  sn0,sd0 = fse.semblanceNumDen(p,g)
  print "semblances for different vertical smoothings:"
  for sigma in [0,2,4,8]:
    ref = RecursiveExponentialFilter(sigma)
    sn = copy(sn0)
    sd = copy(sd0)
    ref.apply1(sn,sn)
    ref.apply1(sd,sd)
    s = fse.semblanceFromNumDen(sn,sd)
    print "sigma =",sigma," s min =",min(s)," max =",max(s)
    title = "semblance: sigma = "+str(sigma)
    plot2(s1,s2,g,s,gmin=0,gmax=1,title=title)

def goAlign():
  s1,s2,g = getImage()
  g = slog(g)
  n1,n2 = len(g[0]),len(g)
  fse = FaultSemblance()
  #p = fse.slopes(g)
  p = getSlopes(g)
  ref = RecursiveExponentialFilter(4)
  sn,sd = fse.semblanceNumDen(p,g)
  ref.apply1(sn,sn)
  ref.apply1(sd,sd)
  s = fse.semblanceFromNumDen(sn,sd)
  plot2(s1,s2,g,s,gmin=0,gmax=1,title="semblance with alignment")
  p = zerofloat(n1,n2) # semblance with zero slopes
  sn,sd = fse.semblanceNumDen(p,g)
  ref.apply1(sn,sn)
  ref.apply1(sd,sd)
  s = fse.semblanceFromNumDen(sn,sd)
  plot2(s1,s2,g,s,gmin=0,gmax=1,title="semblance without alignment")

def goSlopes():
  s1,s2,g = getImage()
  n1,n2 = s1.count,s2.count
  plot2(s1,s2,g,title="input")
  #g = slog(mul(3.0,g))
  #plot2(s1,s2,g,title="log input")
  p = getSlopes(g)
  plot2(s1,s2,g,p,gmin=-1.5,gmax=1.5,gmap=bwrf,title="slopes")

def getSlopes(g):
  n1,n2 = len(g[0]),len(g)
  lsf = LocalSlopeFinder(8.0,2.0,10.0)
  p = zerofloat(n1,n2)
  lsf.findSlopes(g,p)
  return p

###
def xfindShifts(sigma,min1,max1,lag2,f):
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

def imageMbs():
  #ft,fx = 0.300,0.000
  #dt,dx = 0.002,0.01674
  f1,f2 = 0.0,0.0
  d1,d2 = 1.0,1.0
  n1,n2 = 501,560
  fileName = "/data/seis/mbs/dat/s2/gs667.dat"
  x = readImage(n1,n2,fileName)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return s1,s2,x
 
#############################################################################
# plotting

def plot2(s1,s2,f,g=None,gmin=None,gmax=None,gmap=None,
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
  panel.setColorBarWidthMinimum(200)
  #panel.setHLabel("Inline (km)")
  #panel.setVLabel("Time (s)")
  pv = panel.addPixels(s1,s2,f)
  #pv = panel.addPixels(f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  if max(f)>5:
    pv.setClips(-22,22)
  else:
    pv.setClips(-2.1,2.1)
  if g:
    alpha = 0.5
    pv = panel.addPixels(s1,s2,g)
    #pv = panel.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    if gmin==None: gmin = min(g)
    if gmax==None: gmax = max(g)
    if gmap==None: gmap = jetf
    pv.setClips(gmin,gmax)
    pv.setColorModel(gmap)
  frame2(panel,png)

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
jetf = jetFill(0.5)
jetr = jetRamp(1.0)
bwrf = bwrFill(0.5)
bwrn = bwrNotch(1.0)

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
  #frame.setFontSizeForSlide(1.0,0.8)
  #frame.setSize(1290,777)
  #frame.setSize(1490,977)
  frame.setSize(1290,815)
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
