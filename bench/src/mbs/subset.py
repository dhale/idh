"""
Make small subsets for research.
"""

from mbsutils import *
setupForSubset("s1")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gfile = "g" # raw image obtained from SEG-Y file
gsfile = "gs" # image after smoothing with soblf

def main(args):
  goSmall3dWithFault()
  #goSliceAcrossFault()

def goSmall3dWithFault():
  """
  i1=[0:279],i2=[407:768],i3=[162:495] copied from subset 1
  in this subset i2=260 corresponds to xline=1157
  """
  g = readImage(gfile)
  m1,m2,m3 = 280,362,334
  j1,j2,j3 =   0,407,162
  h = copy(m1,m2,m3,j1,j2,j3,g)
  writeImage("sb/g",h)
  sf = SimpleFrame()
  sf.addImagePanels(h).setClips(-1,1)

def goSliceAcrossFault():
  """i3=667 copied from subset 1"""
  g = readImage(gfile)
  h = zerofloat(n1,n3)
  for i3 in range(0,n3):
    copy(g[i3][667],h[i3])
  aos = ArrayOutputStream("sa/g667.dat")
  aos.writeFloats(h)
  aos.close()
  SimplePlot.asPixels(h)

#############################################################################
run(main)
