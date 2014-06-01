"""
Resamples 3D seismic images from DOE coordinates to CSM coordinates.
Resampling includes translation and rotation, in addition to changing
spatial sampling intervals. The translation and rotation better aligns
the (x2,x3) axes with the predominant strike and dip of geologic 
structures apparent in the images.

Author: Dave Hale, Colorado School of Mines
Version: 2009.06.07
"""
from imports import *

#############################################################################
def main(args):
  process("st")
  #process("sz")
  #process("tz")

def process(what):
  setGlobals(what)
  #resample()
  display()

def setGlobals(what):
  global doeFile,csmFile
  global s1d,s2d,s3d # DOE sampling
  global s1c,s2c,s3c # CSM sampling
  tpDir = "/data/seis/tp/"
  if what=="st": # seismic time image
    s1d = Sampling(1501,0.002,0.000)
    s1c = Sampling(1501,0.002,0.000)
    doeFile = tpDir+"doe/3D_Seismic/tpstAll.dat"
    csmFile = tpDir+"csm/seismict/tpst.dat"
  elif what=="sz": # seismic depth image
    s1d = Sampling(2762,0.002,0.000)
    s1c = Sampling(2762,0.002,0.000)
    doeFile = tpDir+"tss/tpszAll.dat"
    csmFile = tpDir+"csm/seismicz/tpsz.dat"
  elif what=="tz": # time depth image
    doeFile = tpDir+"tss/tptzAll.dat"
    csmFile = tpDir+"csm/seismicz/tptz.dat"
  s2d = Sampling(188,0.033528,0.000) # 0.033528 km = 110 ft
  s3d = Sampling(345,0.033528,0.000)
  s2c = Sampling(357,0.025,0.000)
  s3c = Sampling(161,0.025,0.000)
  print "doeFile =",doeFile
  print "csmFile =",csmFile

def resample():
  n1d,n2d,n3d = s1d.count,s2d.count,s3d.count
  d1d,d2d,d3d = s1d.delta,s2d.delta,s3d.delta
  f1d,f2d,f3d = s1d.first,s2d.first,s3d.first
  n1c,n2c,n3c = s1c.count,s2c.count,s3c.count
  d1c,d2c,d3c = s1c.delta,s2c.delta,s3c.delta
  f1c,f2c,f3c = s1c.first,s2c.first,s3c.first
  x2d = zerofloat(n2c,n3c) # DOE coordinates at
  x3d = zerofloat(n2c,n3c) # which to interpolate
  for i3 in range(n3c):
    x3c = f3c+i3*d3c
    for i2 in range(n2c):
      x2c = f2c+i2*d2c
      csm = Coordinates.Csm(x2c,x3c)
      doe = Coordinates.Doe(csm)
      x2d[i3][i2] = doe.x2
      x3d[i3][i2] = doe.x3
  x = zerofloat(n1d,n2d,n3d) # input array
  y = zerofloat(n1c,n2c,n3c) # output array
  x23 = zerofloat(n2d,n3d) # 23 slice of input array
  y23 = zerofloat(n2c,n3c) # 23 slice of output array
  sx = SimpleFloat3(x)
  sy = SimpleFloat3(y)
  si = SincInterpolator()
  ais = ArrayInputStream(doeFile)
  ais.readFloats(x)
  ais.close()
  for i1 in range(n1c):
    if i1%100==0: print "i1 =",i1
    sx.get23(n2d,n3d,i1,0,0,x23)
    for i3 in range(n3c):
      for i2 in range(n2c):
        y23[i3][i2] = si.interpolate(
          n2d,d2d,f2d,n3d,d3d,f3d,x23,x2d[i3][i2],x3d[i3][i2])
    sy.set23(n2c,n3c,i1,0,0,y23)
  aos = ArrayOutputStream(csmFile)
  aos.writeFloats(y)
  aos.close()

def display():
  n1c,n2c,n3c = s1c.count,s2c.count,s3c.count
  ais = ArrayInputStream(csmFile)
  x = zerofloat(n1c,n2c,n3c)
  ais.readFloats(x)
  ais.close()
  print "x min =",min(x)," max =",max(x)
  SimpleFrame.asImagePanels(s1c,s2c,s3c,x)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
