from java.awt import *
from java.lang import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *
  
def main(args):
  goTest()

def goTest():
  n1,n2 = 75,50
  ta = makeRandomTensors(n1,n2)
  tb = makeRandomTensors(n1,n2)
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  ta.getEigenvalues(au,av)
  pvu = PixelsView(au)
  pvu.setColorModel(ColorMap.HUE)
  pvu.setInterpolation(PixelsView.Interpolation.NEAREST)
  tva = TensorsView(ta)
  tvb = TensorsView(tb)
  tva.setLineColor(Color.RED)
  tvb.setLineColor(Color.BLUE)
  tva.setLineWidth(3)
  tvb.setLineWidth(3)
  tva.setEllipsesDisplayed(10)
  tvb.setEllipsesDisplayed(10)
  sp = SimplePlot()
  tile = sp.plotPanel.getTile(0,0)
  tile.addTiledView(pvu)
  tile.addTiledView(tva)
  tile.addTiledView(tvb)
  sp.addColorBar()

def makeRandomTensors(n1,n2):
  r = mul(10.0,randfloat(n1,n2))
  u1,u2 = cos(r),sin(r)
  au = randfloat(n1,n2)
  av = randfloat(n1,n2)
  return EigenTensors2(u1,u2,au,av)
  
#############################################################################
# Run the function main on the Swing thread
import sys
from javax.swing import *
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
