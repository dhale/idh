#############################################################################
# BlendedGridder test

from java.awt import *
from java.lang import *
from java.util import *
from javax.swing import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

#############################################################################

def main(args):
  Parallel.setParallel(True)
  n1,n2,ns = 101,502,101
  random = Random(01234)
  f = randfloat(random,ns)
  x1 = mul(randfloat(random,ns),n1-1)
  x2 = mul(randfloat(random,ns),n2-1)
  et = LocalOrientFilter(1.0).applyForTensors(randfloat(random,n1,n2))
  bg = BlendedGridder2(et,f,x1,x2)
  #bg = BlendedGridder2(f,x1,x2)
  q = bg.grid(Sampling(n1),Sampling(n2))
  print sum(q)

#############################################################################
# Do everything on Swing thread.
import sys
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
