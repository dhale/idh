import sys
from java.lang import *
from java.util import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *

def main(args):
  # Data from Ehsan Naeini, Ikon Science
  n1,n2 = 101,701
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  s = readImage("data/seis.dat",n1,n2)
  w = readImage("data/well.dat",n1,n2)

  # The tensor field used to guide the interpolation. The seismic image has
  # much less vertical resolution than the well data, so we set the
  # eigenvalues eu to be tiny, which implies uniformly low correlation
  # perpendicular to seismic reflectors. This means that we are using the
  # seismic image to estimate only the dips of reflectors and the lateral
  # correlations along those reflectors.
  lof = LocalOrientFilter(4.0,1.0)
  tensors = lof.applyForTensors(s)
  tensors.invertStructure(1.0,1.0) # normalized inverse; max eigenvalue = 1
  eu = zerofloat(n1,n2) # eigenvalues for eigenvectors u (across reflectors)
  ev = zerofloat(n1,n2) # eigenvalues for eigenvectors v (along reflectors)
  tensors.getEigenvalues(eu,ev)
  fill(0.0001,eu) # correlation across reflectors = 0.01 = sqrt(eu)
  tensors.setEigenvalues(eu,ev)

  # Blended gridder interpolation. When correlation between well logs is low,
  # say, when wells are in different fault blocks, we may want to increase the
  # smoothness parameter. For smoothness > 0.5, derivatives of the interpolant
  # q along the reflectors (that is, in the direction of eigenvectors v) will
  # be zero at the wells, making the locations of wells less apparent in the
  # interpolant.
  fnull = -9999.0
  bg = BlendedGridder2(tensors)
  #bg.setSmoothness(0.5) # with default smoothness, can see the wells
  bg.setSmoothness(4.0) # increased smoothness, especially near wells
  p = copy(w) # p = nearest neighbor interpolant
  q = zerofloat(n1,n2) # q = blended neighbor interpolant
  t = bg.gridNearest(fnull,p) # t = time map (non-euclidean distances)
  bg.gridBlended(t,p,q)

  # Plot everything.
  plot2(s1,s2,s,cmap=gray,perc=99,title="Seismic image")
  plot2(s1,s2,t,cmap=jet,title="Time map")
  plot2(s1,s2,p,cmap=jet,title="Nearest neighbor interpolation")
  plot2(s1,s2,q,cmap=jet,title="Blended neighbor interpolation")
 
def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot2(s1,s2,f,cmap=gray,perc=100,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  pv.setPercentiles(100-perc,perc)

#############################################################################
# Run everything on Swing thread.
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
