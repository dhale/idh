"""
Slices of F3 Demo data.
"""
from f3utils import *

#############################################################################
setupForDataSet("seta")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
clip = 1.0

#############################################################################
def main(args):
  #slice1("g",309)
  #slice3("g",75)
  plotSlice1("g",309)
  plotSlice3("g",75)

def slice1(name,k1):
  x = readImage(name)
  x = reshape(n2,n3,flatten(copy(1,n2,n3,k1,0,0,x)))
  writeImage2(getF3dSlice1Name(k1),x)

def slice3(name,k3):
  x = readImage(name)[k3]
  writeImage2(getF3dSlice3Name(k3),x)

def plotSlice1(name,k1):
  x = readImage2(getF3dSlice1Name(name,k1),n2,n3)
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setSize(1100,800)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  pv = sp.addPixels(s2,s3,x)
  pv.setClips(-clip,clip)

def plotSlice3(name,k3):
  x = readImage2(getF3dSlice3Name(name,k3),n1,n2)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(1100,800)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Time (s)")
  pv = sp.addPixels(s1,s2,x)
  pv.setClips(-clip,clip)

#############################################################################
run(main)
