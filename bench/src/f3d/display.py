"""
Displays F3 data.
"""
from f3utils import *

#############################################################################
setupForDataSet("seta")
s1,s2,s3 = getSamplings()
clip = 1.0

#############################################################################
def main(args):
  display3d()
  #displaySlice3()

def display3d():
  x = readImage("gs8")
  print "x min =",min(x)," max =",max(x)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,x)
  ipg.setClips(-clip,clip)
  view = frame.getOrbitView()
  view.setScale(3.0)
  view.setAxesScale(1.0,1.0,8.0)
  frame.setSize(1200,900)
  frame.setVisible(True)

def displaySlice3():
  n1,n2 = s1.count,s2.count
  x = readImage2(getF3dSlice3Name("g",75),n1,n2)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(x)
  pv.setClips(-clip,clip)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
