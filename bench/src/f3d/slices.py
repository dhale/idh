"""
Slices of F3 Demo data.
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
  k3 = 75 # index of slice
  imagefile = "f3d.dat"
  slicefile = "f3d"+str(k3)+".dat"
  #slice3(75,imagefile,slicefile)
  plot(slicefile)

def slice3(k3,imagefile,slicefile):
  x = readImage(imagefile,n1,n2,n3)
  writeImage(slicefile,x[k3])

def plot(datfile):
  x = readImage(datfile,n1,n2)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setFontSizeForSlide(1.0,1.0)
  sp.setSize(1200,900)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Time (s)")
  pv = sp.addPixels(s1,s2,x)
  pv.setClips(-clip,clip)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
