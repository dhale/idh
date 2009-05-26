import sys
from math import *

from java.awt import *
from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

def main(args):
  #convertToDepth()
  #plot3d("tp3r.dat")
  #plot3d("../Transform/tp3z.dat")
  plot3d("tp3z.dat")
  return

True = 1
False = 0

dataDir = "/data/seis/tp/resamp/"
nt=401; dt=0.002; ft=0.500
nz=401; dz=0.004; fz=0.770
n1=401; d1=0.004; f1=0.770
n2=161; d2=0.025; f2=0.000
n3=357; d3=0.025; f3=0.000

def ftToKm(ft):
  return ft*0.0003048

def getTimeDepth():
  table = [
    (0.000,    0), (0.132,  453), (0.135,  467), (0.138,  479), (0.176,  682),
    (0.526, 2669), (0.558, 2869), (0.563, 2906), (0.580, 3004), (0.600, 3092),
    (0.643, 3302), (0.699, 3613), (0.725, 3754), (0.734, 3805), (0.742, 3851),
    (0.743, 3858), (0.744, 3863), (0.745, 3872), (0.767, 3981), (0.780, 4064),
    (0.815, 4266), (0.835, 4372), (0.842, 4413), (0.851, 4471), (0.854, 4494),
    (0.873, 4634), (0.973, 5423), (0.983, 5505), (0.993, 5576), (0.993, 5580),
    (1.002, 5649), (1.016, 5754), (1.300, 8000)]
  tlist,zlist = [],[]
  for t,z in table:
    tlist.append(t)
    zlist.append(ftToKm(z))
  return tlist,zlist

def convertToDepth():
  tlist,zlist = getTimeDepth()
  n = len(tlist)
  ci = CubicInterpolator(CubicInterpolator.Method.LINEAR,n,zlist,tlist)
  z = Array.rampfloat(fz,dz,nz)
  t = Array.zerofloat(nz)
  ci.interpolate(nz,z,t)
  p = Array.zerofloat(nt)
  q = Array.zerofloat(nz)
  si = SincInterpolator()
  si.setUniformSampling(nt,dt,ft)
  si.setUniformSamples(p)
  tfile = "tp3r.dat"
  zfile = "tp3z.dat"
  ais = ArrayInputStream(dataDir+tfile)
  aos = ArrayOutputStream(dataDir+zfile)
  for i3 in range(n3):
    for i2 in range(n2):
      ais.readFloats(p)
      si.interpolate(nz,t,q)
      aos.writeFloats(q)
  ais.close()
  aos.close()

def plot3d(fileName):
  x = Array.zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+fileName)
  ais.readFloats(x)
  ais.close()
  print "x min =",Array.min(x)," max =",Array.max(x)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  #clip = 1.0e-6
  #ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
