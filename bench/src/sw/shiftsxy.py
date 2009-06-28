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
from edu.mines.jtk.util.ArrayMath import *

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

datadir = "/data/seis/sw/sub24/"

##############################################################################
# Velocity function of time v0(t)

def makev0():
  global v0
  # v0(t) = 3.5 + (t-2.4)*0.5 (a simple linear function)
  # v0 = rampfloat(3.5,0.5*dt,n1)
  # v0 = v0(t) = z'(t), for (z,t) in (m,ms) from checkshot survey
  # depths in ft from checkshot survey
  zl= [0.00,  3168.00,  3568.00,  4018.00,  4418.00,
    5168.00,  5618.00,  6018.00,  6418.00,  7168.00,
    7568.00,  8018.00,  8768.00,  9218.00,  9618.00,
   10028.00, 10718.00, 11168.00, 11618.00, 12368.00,
   12968.00, 13408.00, 13918.00, 14368.00, 14768.00,
   15018.00, 15222.00, 15518.00, 15721.50, 16078.17,
   16468.00, 16768.00, 17068.00, 17368.00, 17618.00,
   17908.00]
  # two-way times in ms from checkshot survey
  tl= [0.00,  1051.00,  1174.80,  1317.00,  1444.80,
    1682.40,  1831.20,  1962.00,  2089.00,  2317.20,
    2437.20,  2567.20,  2776.60,  2896.80,  3001.00,
    3099.40,  3233.00,  3302.00,  3361.00,  3445.60,
    3512.80,  3567.40,  3637.60,  3699.80,  3754.40,
    3788.20,  3814.40,  3854.80,  3883.36,  3965.37,
    4051.25,  4113.00,  4157.80,  4216.40,  4264.00,
    4310.00]
  nt = len(tl)
  z = zerofloat(nt)
  t = zerofloat(nt)
  copy(zl,z)
  copy(tl,t)
  z = mul(0.001*0.3048,z) # depths in km
  t = mul(0.001,t) # two-way times in s
  #dump(z)
  #dump(t)
  ci = CubicInterpolator(CubicInterpolator.Method.LINEAR,nt,t,z)
  v0 = zerofloat(n1)
  for i1 in range(n1):
    t1 = (f1+i1*d1)
    v0[i1] = 2.0*ci.interpolate1(t1)
def plotv0():
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,800)
  sp.addPoints(s1,v0)
  sp.setHLabel("velocity (km/s)")
  sp.setVLabel("time (s)")
  sp.paintToPng(300,3,"vt.png")

##############################################################################
# Read/write

def readFloats3(file):
  f = zerofloat(n1,n2,n3)
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
  mul(dx,deltax,deltax) # delta x in samples
  ddxdt = shifts.dddt(deltax); # delta x in m
  writeFloats3("c3s3.dat",ddxdt)
  print "min =",min(ddxdt)," max =",max(ddxdt)

def estimateDdxdt():
  print "estimateDdxdt:",
  deltat = readFloats3("u1s3.dat") # delta t in samples
  mul(dt,deltat,deltat) # delta t in s
  ddxdt = shifts.ddxdt(r,v0,deltat)
  writeFloats3("e3s3.dat",ddxdt)
  print "min =",min(ddxdt)," max =",max(ddxdt)

def estimateDeltaX():
  print "estimateDeltaX:",
  deltat = readFloats3("u1s3.dat") # delta t in samples
  mul(dt,deltat,deltat) # delta t in s
  deltax = shifts.getDeltaX(r,v0,deltat) # delta x in km
  mul(1.0/dx,deltax,deltax) # delta x in samples
  writeFloats3("e3s3.dat",deltax)
  print "min =",min(deltax)," max =",max(deltax)

def estimateDeltaY():
  print "estimateDeltaY:",
  deltat = readFloats3("u1s3.dat") # delta t in samples
  mul(dt,deltat,deltat) # delta t in s
  deltay = shifts.getDeltaY(r,v0,deltat) # delta y in km
  mul(1.0/dy,deltay,deltay) # delta y in samples
  writeFloats3("e2s3.dat",deltay)
  print "min =",min(deltay)," max =",max(deltay)

def main(args):
  makev0()
  plotv0()
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
