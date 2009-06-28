import sys
from math import *

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

dataDir = "/data/as/heart/"

def main(args):
  doPlot3d()
  return

def doPlot3d():
  x = readImage("image00.v")
  #x = readImage("image03.v")
  plot3d(x)

def readImage(fileName):
  ais = ArrayInputStream(dataDir+fileName,ByteOrder.LITTLE_ENDIAN)
  global n1,n2,n3
  n1 = ais.readInt()
  n2 = ais.readInt()
  n3 = ais.readInt()
  print "readImage: n1 =",n1," n2 =",n2," n3 =",n3
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(x,fileName):
  aos = ArrayOutputStream(dataDir+fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeInt(len(x[0][0]))
  aos.writeInt(len(x[0]))
  aos.writeInt(len(x))
  aos.writeFloats(x)
  aos.close()

def plot3d(x):
  print "x min =",min(x)," max =",max(x)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  #clip = 2.5e4
  #ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  view = frame.getOrbitView()
  view.setAxesOrientation(View.AxesOrientation.XUP_YRIGHT_ZOUT)
  frame.setVisible(True)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
