"""
Displays F3 Demo data.
"""
from common import *

#############################################################################
clip = 5.0
n1,n2,n3 = 462,951,591
d1,d2,d3 = 0.004,0.025,0.025
f1,f2,f3 = 0.004,0.000,0.000
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

#############################################################################
def main(args):
  datfile = "f3d.dat"
  display3d(datfile)
  #displaySlice3(datfile)

def display3d(datfile):
  x = readImage(datfile,n1,n2,n3)
  print "x min =",min(x)," max =",max(x)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,x)
  ipg.setClips(-clip,clip)
  view = frame.getOrbitView()
  view.setScale(3.0)
  view.setAxesScale(1.0,1.0,8.0)
  frame.setSize(1200,900)
  frame.setVisible(True)

def displaySlice3(datfile):
  k3 = 75
  x = readImage(datfile,n1,n2,n3)
  y = x[k3]
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(y)
  pv.setClips(-clip,clip)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
