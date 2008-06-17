import sys
from math import *

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

from ldf import *

True = 1
False = 0

dataDir = "/data/seis/tp/"

# Rotated and resampled data.
#n1=1501; d1=0.002; f1=0.000
#n2=161;  d2=0.025; f2=0.000
#n3=357;  d3=0.025; f3=0.000

# After subset and subsampling in time.
n1=251; d1=0.004; f1=0.500
n2=161; d2=0.025; f2=0.000
n3=357; d3=0.025; f3=0.000

cgray = ColorMap.GRAY
cjet = ColorMap.JET

def main(args):
  #goSubset()
  #goOrient()
  #doPlot3d("d2.dat")
  #goPaint()
  #doPlot3dList(["pp.dat","ph.dat","tp3s.dat"],[cjet,cjet,cgray])
  #goPaintSubsample(10)
  #goPaint()
  #doPlot3dList(["pp10.dat","tp3s.dat"],[cjet,cgray])
  #goPaint()
  #goPaintWells(20)
  #goPaint()
  #doPlot3dList(["pp20w2d.dat","pp.dat"],[cjet,cjet])
  #doPlot3dList(["tp3s.dat"],[cgray])
  goSubset2d()

def goSubset2d():
  x = readFloats("tp3s.dat",n1,n2,n3)
  s = Array.copy(n1,1,n3,0,73,0,x);
  s = Array.flatten(s);
  s = Array.reshape(n1,n3,s);
  SimplePlot.asPixels(s);
  writeFloats("tp73.dat",s)
  
def goPaint():
  dt3 = makeDiffusionTensors3(0.0,1.0,0.001)
  p = readFloats("pw20.dat",n1,n2,n3)
  f = makePaintFlags(p)
  lif = LocalInterpolationFilterIc(0.0001,1000)
  lif.apply(dt3,f,p)
  writeFloats("pp20w2d.dat",p)
  plot3d(p,cjet)

def doPlot3dList(fileNameList,cmapList):
  n = len(fileNameList)
  xList = []
  for i in range(n):
    xList.append(readFloats(fileNameList[i],n1,n2,n3))
  plot3dList(xList,cmapList)

def doPlot3d(fileName,cmap):
  doPlot3dList([fileName],[cmap])

def goOrient():
  x = readFloats("tp3s.dat",n1,n2,n3)
  d1,d2,u2,u3,w1,w2 = orient(x)
  writeFloats("d1.dat",d1);
  writeFloats("d2.dat",d2);
  writeFloats("u2.dat",u2);
  writeFloats("u3.dat",u3);
  writeFloats("w1.dat",w1);
  writeFloats("w2.dat",w2);
  plot3d(d2)

def readFloats(fileName,n1,n2,n3):
  x = Array.zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+fileName)
  ais.readFloats(x)
  ais.close()
  return x

def writeFloats(fileName,x):
  aos = ArrayOutputStream(dataDir+fileName)
  aos.writeFloats(x)
  aos.close()

def makePaintFlags(p):
  n1,n2,n3 = len(p[0][0]),len(p[0]),len(p)
  f = Array.zerobyte(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        if p[i3][i2][i1]!=0.0: f[i3][i2][i1] = 1 
  return f

def goPaintSubsample(m):
  p = readFloats("ph.dat",n1,n2,n3)
  ps = Array.zerofloat(n1,n2,n3)
  for i3 in range(0,n3,m):
    for i2 in range(0,n2,m):
      for i1 in range(n1):
        ps[i3][i2][i1] = p[i3][i2][i1]
  plot3d(ps,cjet)
  writeFloats("ph"+str(m)+".dat",ps)

def goPaintWells(m):
  p = readFloats("pp.dat",n1,n2,n3)
  ps = Array.zerofloat(n1,n2,n3)
  for i3 in range(0,n3,m):
    for i2 in range(0,n2,m):
      for i1 in range(n1):
        ps[i3][i2][i1] = p[i3][i2][i1]
  plot3d(ps,cjet)
  writeFloats("pw"+str(m)+".dat",ps)

def makeDiffusionTensors3(s1,s2,s3):
  d1 = readFloats("d1.dat",n1,n2,n3)
  d2 = readFloats("d2.dat",n1,n2,n3)
  u2 = readFloats("u2.dat",n1,n2,n3)
  u3 = readFloats("u3.dat",n1,n2,n3)
  w1 = readFloats("w1.dat",n1,n2,n3)
  w2 = readFloats("w2.dat",n1,n2,n3)
  return DiffusionTensors3(n1,n2,n3,s1,s2,s3,d1,d2,u2,u3,w1,w2)

def orient(x):
  lof = LocalOrientFilter(8,4,4)
  d1 = Array.zerofloat(n1,n2,n3)
  d2 = Array.zerofloat(n1,n2,n3)
  u2 = Array.zerofloat(n1,n2,n3)
  u3 = Array.zerofloat(n1,n2,n3)
  w1 = Array.zerofloat(n1,n2,n3)
  w2 = Array.zerofloat(n1,n2,n3)
  lof.apply(x,None,None,
    None,u2,u3,
    None,None,None,
    w1,w2,None,
    None,None,None,
    d2,d1)
  print "x min =",Array.min(x)," max =",Array.max(x)
  print "d1 min =",Array.min(d1)," max =",Array.max(d1)
  print "d2 min =",Array.min(d2)," max =",Array.max(d2)
  print "u2 min =",Array.min(u2)," max =",Array.max(u2)
  print "u3 min =",Array.min(u3)," max =",Array.max(u3)
  print "w1 min =",Array.min(w1)," max =",Array.max(w1)
  print "w2 min =",Array.min(w2)," max =",Array.max(w2)
  return d1,d2,u2,u3,w1,w2

def goSubset():
  m1,m2,m3 = 251,n2,n3
  j1,j2,j3 = 251,0,0
  x = Array.zerofloat(n1,n2,n3)
  y = Array.zerofloat(2*m1,m2,m3)
  z = Array.zerofloat(m1,m2,m3)
  ais = ArrayInputStream(dataDir+"tp3r.dat")
  ais.readFloats(x)
  ais.close()
  sx = SimpleFloat3(x)
  sx.get123(2*m1,m2,m3,j1,j2,j3,y)
  Array.copy(m1,m2,m3,0,0,0,2,1,1,y,0,0,0,1,1,1,z)
  aos = ArrayOutputStream(dataDir+"tp3s.dat")
  aos.writeFloats(z)
  aos.close()
  plot3d(z)

def plot12(i3):
  ais = ArrayInputStream(dataDir+"tp3.dat")
  x = Array.zerofloat(n1,n2)
  ais.skipBytes(4*n1*n2*i3)
  ais.readFloats(x)
  ais.close()
  print "x min =",Array.min(x)," max =",Array.max(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(x)
  #pv.setPercentiles(1.0,99.0)
  pv.setClips(-1.0e-6,1.0e-6)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)

def plot3dList(xl,cml):
  world = World()
  for i in range(len(xl)):
    x = xl[i]
    cm = cml[i]
    print "i =",i," x min =",Array.min(x)," max =",Array.max(x)
    #s1 = Sampling(n1,d1,f1)
    #s2 = Sampling(n2,d2,f2)
    #s3 = Sampling(n3,d3,f3)
    n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
    s1,s2,s3 = Sampling(n1,1.0,0.0),Sampling(n2,1.0,0.0),Sampling(n3,1.0,0.0)
    x3 = SimpleFloat3(x)
    ipg = ImagePanelGroup(s3,s2,s1,SimpleFloat3(x))
    ipg.setColorModel(cm)
    world.addChild(ipg)
  frame = TestFrame(world)
  frame.orbitView.setScale(2.0)
  frame.setVisible(True)

def plot3d(x,cmap=cjet):
  plot3dList([x],[cmap])

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
