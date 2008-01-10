import sys
from math import *
from java.awt import *
from java.lang import *
from javax.swing import *
from sw import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

True = 1
False = 0

#n1,d1,f1 = 1501,0.004,0.0
#n1,d1,f1 = 301,0.004,3.6
n1,d1,f1 = 601,0.004,2.4
n2,d2,f2 = 623,0.025,0.0
n3,d3,f3 = 367,0.025,0.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)

dx,dy,dt = d3,d2,d1
sx,sy,st = s3,s2,s1
shifts = ShiftsXY(sx,sy,st)
r = 5.0
v0 = Array.rampfloat(3.5,0.5*dt,n1)
#v0 = Array.rampfloat(4.0,0.5*dt,n1)

datadir = "/data/seis/sw/sub/"

##############################################################################
# Read/write

def readFloats3(file):
  f = Array.zerofloat(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readFloats(f)
  af.close()
  return f

def writeFloats3(file,f):
  af = ArrayFile(datadir+file,"rw")
  af.writeFloats(f)
  af.close()
  return f

#############################################################################
# Research

def computeDdxdt():
  print " computeDdxdt:",
  deltax = readFloats3("u3s3.dat")
  Array.mul(dx,deltax,deltax) # delta x in samples
  ddxdt = shifts.dddt(deltax); # delta x in m
  writeFloats3("c3s3.dat",ddxdt)
  print "min =",Array.min(ddxdt)," max =",Array.max(ddxdt)

def estimateDdxdt():
  print "estimateDdxdt:",
  deltat = readFloats3("u1s3.dat") # delta t in samples
  Array.mul(dt,deltat,deltat) # delta t in s
  ddxdt = shifts.ddxdt(r,v0,deltat)
  writeFloats3("e3s3.dat",ddxdt)
  print "min =",Array.min(ddxdt)," max =",Array.max(ddxdt)

def estimateDeltaX():
  print "estimateDeltaX:",
  deltat = readFloats3("u1s3.dat") # delta t in samples
  Array.mul(dt,deltat,deltat) # delta t in s
  deltax = shifts.getDeltaX(r,v0,deltat) # delta x in km
  Array.mul(1.0/dx,deltax,deltax) # delta x in samples
  writeFloats3("e3s3.dat",deltax)
  print "min =",Array.min(deltax)," max =",Array.max(deltax)

def estimateDeltaY():
  print "estimateDeltaY:",
  deltat = readFloats3("u1s3.dat") # delta t in samples
  Array.mul(dt,deltat,deltat) # delta t in s
  deltay = shifts.getDeltaY(r,v0,deltat) # delta y in km
  Array.mul(1.0/dy,deltay,deltay) # delta y in samples
  writeFloats3("e2s3.dat",deltay)
  print "min =",Array.min(deltay)," max =",Array.max(deltay)

def main(args):
  #computeDdxdt()
  #estimateDdxdt()
  estimateDeltaX()
  estimateDeltaY()
  return

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
