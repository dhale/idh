from imports import *

datdir = "/data/seis/f3d/"

def readImage(datfile,n1,n2,n3=-1):
  ais = ArrayInputStream(datdir+datfile)
  if n3<0:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datfile,x):
  ais = ArrayOutputStream(datdir+datfile)
  ais.writeFloats(x)
  ais.close()
