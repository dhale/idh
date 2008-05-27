import sys
from math import *

from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

def main(args):
  #convertHorizons()
  paintHorizons()
  return

True = 1
False = 0

dataDir = "/data/seis/tp/"

# floating point format code should be in bytes 3225-6
# 1 for IBM floating point, 5 for IEEE floating point
nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header
n1i=1501 # number of time samples (1st dimension)
d1i=0.002 # time sampling interval
f1i=0.000 # time of first sample
n2i=188 # number of traces in 2nd (inline) dimension
d2i=0.033
f2i=0.000
n3i=345 # number of traces in 3rd (crossline) dimension
d3i=0.033
f3i=0.000

# Eastings and northings (x and y coordinates) of survey rectangle.
# Obtained from teapot_3d_load.doc.
xll,yll = 788937.0,938846.0 # lower left <=> (i2,i3) = (0,0)
xlr,ylr = 809502.0,939334.0 # lower right <=> (i2,i3) = (n2i-1,0)
xur,yur = 808604.0,977163.0 # lower left <=> (i2,i3) = (n2i-1,n3i-1)
xul,yul = 788039.0,976675.0 # lower left <=> (i2,i3) = (0,n3i-1)

# Scale factors from survey units (km) to eastings and northings.
# The scale factors for x2 and x3 should be almost identical.
x2s = sqrt((xlr-xll)*(xlr-xll)+(ylr-yll)*(ylr-yll))/((n2i-1)*d2i)
x3s = sqrt((xul-xll)*(xul-xll)+(yul-yll)*(yul-yll))/((n3i-1)*d3i)
#print "x2s =",x2s," x3s =",x3s

# Rotation angle of survey rectangle. This is the angle between the 
# x2 (inline) and xe (easting) axes, measured counter-clockwise.
phia = atan2(ylr-yll,xlr-xll)
cosa = cos(phia)
sina = sin(phia)
#print "phia =",phia

# Resampling of (x2,x3) coordinates to trim off excess zeros.
# This is a translation and rotation. The rotation center (in km) 
# and angle (in radians) were determined visually.
# x2Old = x2c + x2New*cos(phi) + x3New*sin(phi)
# x3Old = x3c - x2New*sin(phi) + x3New*cos(phi)
x2c,x3c = 4.192,0.935 # rotation center (km)
phir = -0.485364 # rotation angle (radians, equals -27.8093 degrees)
cosr = cos(phir)
sinr = sin(phir)

# Resampling of x2 and x3 after translation and rotation.
# Spatial sampling intervals become 25 m instead of 33 m.
# Also subsample time.
n1=251; d1=0.004; f1=0.500
n2=161; d2=0.025; f2=0.000
n3=357; d3=0.025; f3=0.000

# From easting-northing (xe,yn) to survey (x2,x3) coordinates.
def toSurvey(xe,yn):
  delx = xe-xll
  dely = yn-yll
  x2 = ( cosa*delx+sina*dely)/x2s
  x3 = (-sina*delx+cosa*dely)/x3s
  return x2,x3

# From survey (x2,x3) coordinates to easting-northing (xe,yn).
def fromSurvey(x2,x3):
  s2 = x2*x2s
  s3 = x3*x3s
  xe = xll+cosa*s2-sina*s3
  yn = yll+sina*s3-cosa*s3
  return xe,yn

# From survey (x2,x3) coordinates to resampled (x2,x3) coordinates
def toResampled(x2,x3):
  del2 = x2-x2c
  del3 = x3-x3c
  r2 = cosr*del2-sinr*del3
  r3 = sinr*del2+cosr*del3
  return r2,r3

# From resampled (x2,x3) coordinates to survey (x2,x3) coordinates
def fromResampled(r2,r3):
  x2 = x2c+cosr*r2+sinr*r3
  x3 = x3c-sinr*r2+cosr*r3
  return x2,x3

#print "toSurvey(xll,yll)",toSurvey(xll,yll)
#print "toSurvey(xlr,ylr)",toSurvey(xlr,ylr)
#print "toSurvey(xur,yur)",toSurvey(xur,yur)
#print "toSurvey(xul,yul)",toSurvey(xul,yul)

def readHorizonLists():
  input = open(dataDir+"CD/3DHorizons.xyz","r")
  lines = input.readlines()
  input.close()
  slast = ""
  index = -1
  slist = []
  xlist = []
  ylist = []
  tlist = []
  lines = lines[1:] # skip headings
  for line in lines:
    a = line.strip().split()
    x,y,s,t = a[2:6]
    s = s.lower().replace("/","")
    if s!=slast:
      slast = s
      slist.append(s)
      xlist.append([])
      ylist.append([])
      tlist.append([])
      index += 1
    xlist[index].append(x)
    ylist[index].append(y)
    tlist[index].append(t)
  n = index+1
  for i in range(n):
    print "horizon",slist[i],"has",len(tlist[i]),"xyt"
  return slist,xlist,ylist,tlist

def convertXY(xlist,ylist,tlist):
  x2list,x3list,tnlist = [],[],[]
  n = len(xlist)
  for i in range(n):
    xs,ys,ts = xlist[i],ylist[i],tlist[i]
    x2i,x3i,ti = [],[],[]
    ns = len(xs)
    for j in range(ns):
      x = float(xs[j])
      y = float(ys[j])
      t = float(ts[j])
      x2,x3 = toSurvey(x,y)
      x2,x3 = toResampled(x2,x3)
      i1 = int((t-f1)/d1+0.5)
      i2 = int((x2-f2)/d2+0.5)
      i3 = int((x3-f3)/d3+0.5)
      if 0<=i1 and i1<n1 and 0<=i2 and i2<n2 and 0<=i3 and i3<n3:
        x2i.append(x2)
        x3i.append(x3)
        ti.append(t)
    x2list.append(x2i)
    x3list.append(x3i)
    tnlist.append(ti)
  return x2list,x3list,tnlist

def writeHorizonLists(slist,x2list,x3list,tlist):
  n = len(slist)
  for i in range(n):
    print "writing",slist[i],"with",len(tlist[i]),"xyt"
    s,x2s,x3s,ts = slist[i],x2list[i],x3list[i],tlist[i]
    ns = len(ts)
    aos = ArrayOutputStream(dataDir+"h"+str(i)+s+".dat")
    aos.writeInt(ns)
    for j in range(ns):
      aos.writeFloat(x2s[j])
      aos.writeFloat(x3s[j])
      aos.writeFloat(ts[j])
    aos.close()

def convertHorizons():
  s,x,y,t = readHorizonLists()
  x2,x3,t = convertXY(x,y,t)
  writeHorizonLists(s,x2,x3,t)

def horizonFilenames():
  list = []
  list.append("h8carlile.dat")
  list.append("h0kf2.dat")
  list.append("h1fallriver.dat")
  list.append("h2lakotamorrison.dat")
  list.append("h3crowmountain.dat")
  list.append("h4redpeak.dat")
  list.append("h5tensleep.dat")
  list.append("h6tensleepbbase.dat")
  list.append("h7basement.dat")
  return list

def readHorizon(fileName):
  ais = ArrayInputStream(dataDir+fileName)
  n = ais.readInt()
  a = Array.zerofloat(3,n)
  ais.readFloats(a)
  ais.close()
  x2 = Array.zerofloat(n)
  x3 = Array.zerofloat(n)
  t = Array.zerofloat(n)
  for i in range(n):
    x2[i] = a[i][0]
    x3[i] = a[i][1]
    t[i] = a[i][2]
  return x2,x3,t

def paintHorizon(index,p,x2,x3,t):
  n = len(t)
  for i in range(n):
    x1i = t[i]
    x2i = x2[i]
    x3i = x3[i]
    i1 = int((x1i-f1)/d1+0.5)
    i2 = int((x2i-f2)/d2+0.5)
    i3 = int((x3i-f3)/d3+0.5)
    p[i3][i2][i1] = index

def paintHorizons():
  p = Array.zerofloat(n1,n2,n3)
  hfiles = horizonFilenames()
  nfile = len(hfiles)
  for ifile in range(nfile):
    x2,x3,t = readHorizon(hfiles[ifile])
    index = ifile+4
    paintHorizon(index,p,x2,x3,t)
  aos = ArrayOutputStream(dataDir+"ph.dat")
  aos.writeFloats(p)
  aos.close()

def plot3d(x):
  #s1 = Sampling(n1,d1,f1)
  #s2 = Sampling(n2,d2,f2)
  #s3 = Sampling(n3,d3,f3)
  n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
  s1,s2,s3 = Sampling(n1,1.0,0.0),Sampling(n2,1.0,0.0),Sampling(n3,1.0,0.0)
  x3 = SimpleFloat3(x)
  ipg = ImagePanelGroup(s3,s2,s1,x3)
  ipg.setColorModel(ColorMap.JET)
  #clip = 1.0e-6
  #ipg.setClips(-clip,clip)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
