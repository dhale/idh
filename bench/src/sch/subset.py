# Make good subsets for fault processing

from schutils import *

# Sampling for subset s2
n1,n2,n3 = 901,376,547 # numbers of samples
d1,d2,d3 = 1.0,1.0,1.0 # sampling intervals (for now)
f1,f2,f3 = 0.0,0.0,0.0 # values for first samples
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
clip = 1.0

def main(args):
  subset2a(True)

def subset2a(make):
  """ subset s2a of subset s2, with low-pass filtering """
  global n1,n2,n3,d1,d2,d3,f1,f2,f3,s1,s2,s3
  j1,j2,j3 = 150,26,0 # first samples in subset
  m1,m2,m3 = 401,350,547 # numbers of samples in subset
  bpf = BandPassFilter(0.0,0.20,0.05,0.01) # frequencies < half Nyquist
  for name in ["g0","gx"]: # without (g0) and with (gx) bilateral filter
    infile = "/data/seis/sch/dat/s2/"+name+".dat"
    outfile = "/data/seis/sch/dat/s2a/"+name+".dat"
    if make:
      g = zerofloat(m1)
      aif = ArrayFile(infile,"r")
      aos = ArrayOutputStream(outfile)
      for i3 in range(j3,j3+m3):
        for i2 in range(j2,j2+m2):
          aif.seek(4*(j1+n1*(i2+n2*i3)))
          aif.readFloats(g)
          bpf.apply(g,g)
          aos.writeFloats(g)
      aos.close()
      aif.close()
  n1,n2,n3 = m1,m2,m3
  d1,d2,d3 = d1,d2,d3
  f1,f2,f3 = f1+j1*d1,f2+j2*d2,f3+j3*d3
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  for name in ["g0","gx"]:
    outfile = "/data/seis/sch/dat/s2a/"+name+".dat"
    ais = ArrayInputStream(outfile)
    g = zerofloat(n1,n2,n3)
    ais.readFloats(g)
    ais.close()
    showOne(g)

#############################################################################

def showOne(g,cmin=None,cmax=None,clip=None,cmap=None):
  sf = SimpleFrame()
  ipg = sf.addImagePanels(s1,s2,s3,g)
  if cmin and cmax:
    ipg.setClips(cmin,cmax)
  elif clip:
    ipg.setClips(-clip,clip)
  else:
    ipg.setPercentiles(2,98)
    print "clips: min =",ipg.clipMin," max =",ipg.clipMax
  if cmap:
    ipg.setColorModel(cmap)
  sf.orbitView.setScale(3.0)
  #sf.orbitView.setAxesScale(1.0,1.0,0.5)
  sf.setSize(1000,800)

#############################################################################
run(main)
