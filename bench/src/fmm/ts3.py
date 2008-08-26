#############################################################################
# Tests 3D time solvers.

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
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *

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
  testConstant()
  return
 
def testConstant():
  n1,n2,n3 = 101,101,101
  d11 = 1.000;
  d12 = 0.900; d22 = 1.000;
  d13 = 0.900; d23 = 0.900; d33 = 1.000 
  tensors = ConstantTensors3(n1,n2,n3,d11,d12,d13,d22,d23,d33)
  ts = TimeSolver3(n1,n2,n3,tensors)
  sw = Stopwatch()
  sw.start()
  ts.zeroAt(2*(n1-1)/4,2*(n2-1)/4,2*(n3-1)/4)
  sw.stop()
  print "time =",sw.time()
  plot(ts.getTimes(),prism)

#############################################################################
# tensors

class ConstantTensors3(EigenTensors3):
  """
  3D constant tensors
  """
  def __init__(self,n1,n2,n3,d11,d12,d13,d22,d23,d33):
    EigenTensors3.__init__(self,n1,n2,n3,False)
    for i3 in range(n3):
      for i2 in range(n2):
        for i1 in range(n1):
          self.setTensor(i1,i2,i3,d11,d12,d13,d22,d23,d33)

#############################################################################
# plotting

def plot(t,cm):
  ipg = ImagePanelGroup(t)
  ipg.setColorModel(cm)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
