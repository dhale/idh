#############################################################################
# Henry Pinkard's test

from java.lang import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *

def main(args):
  xs = [-200.0,-200.0, 200.0, 200.0,  0.0]
  ys = [-200.0, 200.0,-200.0, 200.0,  0.0]
  zs = [  10.0,  10.0,  10.0,  10.0,  0.0]
  basis = RadialInterpolator2.Biharmonic()
  ri = RadialInterpolator2(basis,zs,xs,ys)
  rg = RadialGridder2(basis,zs,xs,ys)
  sg = SplinesGridder2(zs,xs,ys)
  nx,ny = 100,100

  sx1 = Sampling(nx,200.0/(nx-1),-100.0)
  sy1 = Sampling(ny,200.0/(ny-1),-100.0)
  z1 = rg.grid(sx1,sy1)
  #z1 = sg.grid(sx1,sy1)

  sx2 = Sampling(nx,300.0/(nx-1),-150.0)
  sy2 = Sampling(ny,300.0/(ny-1),-150.0)
  z2 = rg.grid(sx2,sy2)
  #z2 = sg.grid(sx2,sy2)

  zd = sub(z1,z2)
  print min(z1),min(z2),min(zd)
  print max(z1),max(z2),max(zd)

  sp = SimplePlot()
  sp.addPixels(sx2,sy2,z2).setClips(0.0,4.0)
  sp.addPixels(sx1,sy1,z1).setClips(0.0,4.0)
  
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
