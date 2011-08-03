import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fd import *

True = 1
False = 0

#############################################################################
# parameters

fontSize = 24
width = 800
height = 600
pngDir = "./png"
#pngDir = None

#############################################################################
# functions

def main(args):
  test1()
  return

def plot(w,png):
  nw = len(w)
  nx = len(w[0])
  colors = [Color.RED,Color.GREEN,Color.BLUE]
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.setVLimits(-30.0,30.0)
  dx = 1.0
  fx = -0.5*(nx-1)*dx
  sx = Sampling(nx,dx,fx)
  for iw in range(nw):
    pv = pp.addPoints(sx,w[iw])
    pv.setLineColor(colors[iw%3])
    pv.setLineWidth(3)
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setFontSizeForSlide(1.0,0.9)
  pf.setSize(width,height)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(720,3.3,pngDir+"/"+png+".png")
 
def makeModel(dl,dr):
  nx = 801
  d = zerofloat(nx)
  v = zerofloat(nx)
  for ix in range(nx):
    v[ix] = 1.0
    if ix<3*nx/4:
      d[ix] = dl
    else:
      d[ix] = dr
  return d,v

def test1():
  nmodel = 4
  lmodel = ["05","03","02","01"]
  dllist = [1.0,1.0,1.0,1.0] # densities left
  drlist = [0.5,0.3,0.2,0.1] # densities right
  for imodel in range(nmodel):
    dl = dllist[imodel]
    dr = drlist[imodel]
    d,v = makeModel(dl,dr)
    nx = len(d)
    rc = (dr-dl)/(dr+dl)
    tc = 1.0+rc
    print "rc =",rc,"  tc =",tc
    dt = 0.5
    dx = 1.0
    nt = 5
    mt = 125
    methods = [
      Wave1.Method.SYMMETRIC,
      Wave1.Method.PRODUCT1,
      Wave1.Method.PRODUCT2
    ]
    nmethod = len(methods)
    f = zerofloat(nx,nmethod,nt)
    for imethod in range(nmethod):
      wave1 = Wave1(methods[imethod],dt,dx,d,v)
      for it in range(nt):
        fi = wave1.step(mt)
        ft = copy(nx/2,fi)
        copy(fi,f[it][imethod])
        fmax = max(ft)
        #print "fmax =",fmax,"  fr =",fmax*rc,"  ft =",fmax*tc
    prefix = "wave1"+lmodel[imodel]
    if imodel==0:
      for it in range(nt):
        png = prefix+str(it)
        plot(f[it],png)
    else:
      it = nt-1
      png = prefix+str(it)
      plot(f[it],png)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
