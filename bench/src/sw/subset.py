# Subsets the large Shearwater data volume for research processing.
# The time sampling interval is 0.004 s.
# Times sampled in the original volume are [0.0,6.0] s (1501 samples).
# Times sampled in the subset volume are [3.2,4.8] s (301 samples).

import sys
from java.lang import *
from edu.mines.jtk.io import *
from edu.mines.jtk.util import *

n1a = 1501 # = 6.0 s (number of time samples in original volume)
n1s = 601  # = 2.4 s (number of time samples in subset)
j1s = 600  # = 2.4 s (time of first sample in subset)

n2 = 623 # number of samples in inline direction 
n3 = 367 # number of samples in crossline direction

dataall = "/data/seis/sw/all/"
datasub = "/data/seis/sw/sub/"

##############################################################################

def subset(file):
  afile = dataall+file
  sfile = datasub+file
  print "subset:",afile,"to",sfile
  afa = ArrayFile(afile,"r")
  afs = ArrayFile(sfile,"rw")
  fa = Array.zerofloat(n1a)
  fs = Array.zerofloat(n1s)
  for i3 in range(n3):
    for i2 in range(n2):
      afa.readFloats(fa)
      Array.copy(n1s,j1s,fa,0,fs)
      afs.writeFloats(fs)
  afa.close()
  afs.close()

def subsetAll():
  subset("s02.dat")
  subset("s04.dat")
  subset("w02.dat")
  subset("w04.dat")
  subset("u1s3.dat")
  subset("u2s3.dat")
  subset("u3s3.dat")
subsetAll()
