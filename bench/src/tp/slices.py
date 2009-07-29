"""
Makes subdirectories with slices of seismic time or depth images.
For example, the directory with name "s3_84" contains a constant-i3
slice, where i3 = 84.
"""
from tputils import *
setupForSubset("subz_401_4_600")
seismicDir = getSeismicDir()

#############################################################################
def main(args):
  makeSlice3(96)

def makeSlice3(i3):
  subDir = "s3_"+str(i3)+"/"
  File(seismicDir+subDir).mkdir()
  for name in ["tpsz","tpgv","tpgd","tpgg","tpgp"]:
    x = readImage(name)
    writeImage(subDir+name,x[i3])

#############################################################################
run(main)
