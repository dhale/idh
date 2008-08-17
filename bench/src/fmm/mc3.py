#############################################################################
# Tests MarchingCubes

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
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *

from fmm import *

#############################################################################
# global parameters

dataDir = "/data"
#pngDir = "./png"
pngDir = None

#############################################################################
# test functions

def main(args):
  testSphere()
  return
 
def testSphere():
  n1,n2,n3 = 131,131,131
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  f = Array.zerofloat(n1,n2,n3)
  x1 = (n1-1)/2.0
  x2 = (n2-1)/2.0
  x3 = (n3-1)/2.0
  for i3 in range(n3):
    d3 = i3-x3;
    for i2 in range(n2):
      d2 = i2-x2;
      for i1 in range(n1):
        d1 = i1-x1;
        f[i3][i2][i1] = sqrt(d1*d1+d2*d2+d3*d3)
  c = 0.5*Array.max(f)
  mc = MarchingCubes(s1,s2,s3,f)
  contour = mc.getContour(c)
  plot(s1,s2,s3,f,contour)

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

def plot(s1,s2,s3,f,contour):
  tg = TriangleGroup(contour.i,contour.x,contour.u)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.CYAN)
  states.add(cs)
  #lms = LightModelState()
  #lms.setTwoSide(True)
  #states.add(lms)
  ms = MaterialState()
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
  ms.setSpecular(Color.WHITE)
  ms.setShininess(100.0)
  states.add(ms)
  tg.setStates(states);
  ipg = ImagePanelGroup(s1,s2,s3,SimpleFloat3(f))
  ipg.setColorModel(jet)
  world = World()
  world.addChild(tg)
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
