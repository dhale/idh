import sys
from math import *

from java.lang import *
from java.io import *
from java.nio import *
from javax.swing import *
from org.python.util import PythonObjectInputStream

from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

dataDir = "/data/as/heart/"
imageFile = dataDir+"image00.v"
tensorsFile = dataDir+"tensors00.v"
fFile = dataDir+"image00.v" # input image
gFile = dataDir+"i00_lg.v" # linear-smoothed image
ffFile = dataDir+"i00_lff.v" # square of linear-smoothed image
ggFile = dataDir+"i00_lgg.v" # linear-smoothed image-squared
sffFile = dataDir+"i00_lsff.v" # numerator of linear semblance
sggFile = dataDir+"i00_lsgg.v" # denominator of linear semblance
lsFile = dataDir+"i00_ls.v" # linear semblance
sigma = 8.0
sigma1 = 3.0
sigma2 = 3.0
small = 0.01
niter = 100

def main(args):
  #doTranspose()
  #doPlot3d()
  #makeTensors()
  makeLinearSemblance()
  return

def doTranspose():
  x = readImage("image00.v")
  y = transpose13(x)
  writeImage(y,"image00t.v")

def doPlot3d():
  x = readImage("image00t.v")
  #x = readImage("image03t.v")
  plot3d(x)

def readImage(fileName):
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  global n1,n2,n3
  n1 = ais.readInt()
  n2 = ais.readInt()
  n3 = ais.readInt()
  print "n1 =",n1," n2 =",n2," n3 =",n3
  x = Array.zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(x,fileName):
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeInt(len(x[0][0]))
  aos.writeInt(len(x[0]))
  aos.writeInt(len(x))
  aos.writeFloats(x)
  aos.close()

def makeLinearSemblance():
  f = readImage(fFile)
  g = Array.zerofloat(n1,n2,n3)
  scale1 = 0.5*sigma1*sigma1
  scale2 = 0.5*sigma2*sigma2
  lsf = LocalSmoothingFilter(small,niter)
  t = readTensors(tensorsFile)
  t.setEigenvalues(0.0,0.0,1.0)
  lsf.apply(t,scale1,f,g)
  writeImage(g,gFile)
  ff = Array.mul(f,f)
  gg = Array.mul(g,g)
  lsf.apply(t,scale1,Array.copy(ff),ff)
  writeImage(ff,ffFile)
  writeImage(gg,ggFile)
  sff = Array.zerofloat(n1,n2,n3)
  sgg = Array.zerofloat(n1,n2,n3)
  t.setEigenvalues(1.0,1.0,0.0)
  lsf.apply(t,scale2,ff,sff)
  lsf.apply(t,scale2,gg,sgg)
  writeImage(sff,sffFile)
  writeImage(sgg,sggFile)
  ls = Array.div(sgg,sff)
  ls = Array.clip(0.0,1.0,ls)
  writeImage(ls,lsFile)

def makeTensors():
  tensors = computeTensors(imageFile)
  writeTensors(tensors,tensorsFile)
 
def readTensors(tensorsFile):
  fis = FileInputStream(tensorsFile)
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors
 
def writeTensors(tensors,tensorsFile):
  fos = FileOutputStream(tensorsFile)
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

def computeTensors(imageFile):
  f = readImage(imageFile)
  lof = LocalOrientFilter(sigma)
  d = lof.applyForTensors(f)
  return d

def plot3d(x):
  print "x min =",Array.min(x)," max =",Array.max(x)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  #clip = 2.5e4
  #ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  #view = frame.getOrbitView()
  #view.setAxesOrientation(View.AxesOrientation.XRIGHT_YUP_ZOUT)
  frame.setVisible(True)

def transpose13(x):
  n1 = len(x[0][0])
  n2 = len(x[0])
  n3 = len(x)
  y = Array.zerofloat(n3,n2,n1)
  xf = Array.flatten(x) # n1*n2*n3
  xi = Array.reshape(n1,n2*n3,xf)
  for i23 in range(n2*n3):
    xi[i23] = Array.reverse(xi[i23])
  xt = Array.transpose(xi) # n2*n3 * n1
  for i1 in range(n1):
    y[i1] = Array.transpose(Array.reshape(n2,n3,xt[i1]))
  return y

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
