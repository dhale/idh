"""
Make amplitude spectra for all traces in images.
"""

from mbsutils import *
setupForSubset("s1")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gfile = "g" # raw image obtained from SEG-Y file
gsfile = "gs" # image after smoothing with soblf

def main(args):
  goAmplitude()

def goAmplitude():
  g = readImage(gfile)
  gs = readImage(gsfile)
  ga = amplitudeSpectra(g)
  gsa = amplitudeSpectra(gs)
  sf = SimpleFrame()
  ipg = sf.addImagePanels(ga)
  ipg.setColorModel(ColorMap.JET)
  ipg = sf.addImagePanels(gsa)
  ipg.setColorModel(ColorMap.JET)
  sf.orbitView.scale = 2.0
  sf.orbitView.setAxesScales(1.0,1.0,2.0)
  sf.setSize(1500,1000)

def amplitudeSpectra(f):
  fft = Fft(n1)
  a = zerofloat(0,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      a[i3][i2] = cabs(fft.applyForward(f[i3][i2]))
  return a

#############################################################################
run(main)


