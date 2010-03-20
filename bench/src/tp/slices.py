"""
Makes subdirectories with slices of seismic time or depth images.
For example, the directory with name "s3_84" contains a constant-i3
slice, where i3 = 84.
"""
from tputils import *
#setupForSubset("subz_401_4_600")
setupForSubset("subt_251_4_500")
seismicDir = getSeismicDir()

#############################################################################
def main(args):
  #makeSlice3Z(96)
  makeSlice3T(84)
  makeSlice3T(73)

def makeSlice3T(i3):
  subDir = "s3_"+str(i3)+"/"
  File(seismicDir+subDir).mkdir()
  for name in ["tpst"]:
    x = readImage(name)
    writeImage(subDir+name,x[i3])
    display(x[i3])

def makeSlice3Z(i3):
  subDir = "s3_"+str(i3)+"/"
  File(seismicDir+subDir).mkdir()
  for name in ["tpsz","tpgv","tpgd","tpgg","tpgp"]:
    x = readImage(name)
    writeImage(subDir+name,x[i3])

def display(s,g=None,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.getPlotPanel().setColorBarWidthMinimum(80)
  pv = sp.addPixels(s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g!=None:
    pv = sp.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.getJet(0.3))
    if cmin!=cmax:
      pv.setClips(cmin,cmax)

#############################################################################
run(main)
