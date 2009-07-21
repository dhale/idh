#############################################################################
# Tests 3D time marker.

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
from edu.mines.jtk.util.ArrayMath import *

from fmm import *

#############################################################################
# global parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

def main(args):
  testConstant()
  return
 
def testConstant():
  n1,n2,n3 = 101,101,101
  d11 = 1.000;
  d12 = 0.900; d22 = 1.000;
  d13 = 0.900; d23 = 0.900; d33 = 1.000 
  tensors = ConstantTensors3(n1,n2,n3,d11,d12,d13,d22,d23,d33)
  i1 = [3*(n1-1)/8,5*(n1-1)/8]
  i2 = [5*(n2-1)/8,3*(n2-1)/8]
  i3 = [3*(n3-1)/8,5*(n3-1)/8]
  testMarker(n1,n2,n3,i1,i2,i3,tensors)

def testMarker(n1,n2,n3,i1,i2,i3,tensors):
  times = fillfloat(FLT_MAX,n1,n2,n3)
  marks = zeroint(n1,n2,n3)
  for k in range(len(i1)):
    k1,k2,k3 = i1[k],i2[k],i3[k]
    times[k3][k2][k1] = 0.0
    marks[k3][k2][k1] = 1+k
  tm = TimeMarker3(n1,n2,n3,tensors)
  sw = Stopwatch()
  sw.start()
  tm.apply(times,marks)
  sw.stop()
  print "elapsed time =",sw.time()
  print "times: min =",min(times),"max =",max(times)
  print "marks: min =",min(marks),"max =",max(marks)
  marks = floatsFromInts(marks)
  plot(times,prism)
  plot(marks,jet)

#############################################################################
# other functions

def floatsFromInts(i):
  n1 = len(i[0][0])
  n2 = len(i[0])
  n3 = len(i)
  f = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        f[i3][i2][i1] = i[i3][i2][i1]
  return f

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
