import sys
from math import *

from java.awt import *
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

from tp import *

def main(args):
  #plotTimeDepth()
  plotCoordinates()
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
d2i=0.033528 # 110 ft = 33.52800 m
f2i=0.000
n3i=345 # number of traces in 3rd (crossline) dimension
d3i=0.033528 # 110 ft = 33.52800 m
f3i=0.000

# Map coordinates (xe=easting,yn=northing) of survey rectangle.
# Coordinates measured in ft, obtained from teapot_3d_load.doc.
xll,yll = 788937.0,938846.0 # lower left  <=> (i2,i3) = (  0,  0)
xlr,ylr = 809502.0,939334.0 # lower right <=> (i2,i3) = (187,  0)
xur,yur = 808604.0,977163.0 # upper right <=> (i2,i3) = (187,344)
xul,yul = 788039.0,976675.0 # upper left  <=> (i2,i3) = (  0,344)

# Map coordinates from trace headers in Transform's depth image.
# In this image, the first trace of each original line is absent.
#xll,yll = 789048.0,938850.0 # lower left  <=> (i2,i3) = (  0,  0)
#xlr,ylr = 809502.0,939336.0 # lower right <=> (i2,i3) = (186,  0)
#xur,yur = 808604.0,977165.0 # upper right <=> (i2,i3) = (186,344)
#xul,yul = 788150.0,976679.0 # upper left  <=> (i2,i3) = (  0,344)

# Convert coordinates from ft to km.
def ftToKm(ft):
  return ft*0.0003048
xll,yll = ftToKm(xll),ftToKm(yll)
xlr,ylr = ftToKm(xlr),ftToKm(ylr)
xur,yur = ftToKm(xur),ftToKm(yur)
xul,yul = ftToKm(xul),ftToKm(yul)

# Scale factors from survey units (km) to eastings and northings.
# These scale factors should equal one, and are currently unused,
# because eastings and northings measured in ft have already been
# converted to km.
#x2s = sqrt((xlr-xll)*(xlr-xll)+(ylr-yll)*(ylr-yll))/((n2i-1)*d2i)
#x3s = sqrt((xul-xll)*(xul-xll)+(yul-yll)*(yul-yll))/((n3i-1)*d3i)
#print "x2s =",x2s," x3s =",x3s

# Rotation angle of survey rectangle. This is the angle between the 
# x2 (inline) and xe (easting) axes, measured counter-clockwise.
phia = atan2(ylr-yll,xlr-xll)
cosa = cos(phia)
sina = sin(phia)
#print "phia =",phia*180/pi,"degrees"

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
n1=401; d1=0.002; f1=0.500
n2=161; d2=0.025; f2=0.000
n3=357; d3=0.025; f3=0.000

# From map (xe,yn) to survey (x2,x3)
def toSurveyFromMap(xe,yn):
  xe -= xll
  yn -= yll
  x2 =  cosa*xe+sina*yn
  x3 = -sina*xe+cosa*yn
  return x2,x3

# From survey (x2,x3) to map (xe,yn)
def toMapFromSurvey(x2,x3):
  xe = cosa*x2-sina*x3
  yn = sina*x2+cosa*x3
  xe += xll
  yn += yll
  return xe,yn

# From survey (x2,x3) to resampled (x2,x3)
def toResampledFromSurvey(x2,x3):
  x2 -= x2c
  x3 -= x3c
  r2 = cosr*x2-sinr*x3
  r3 = sinr*x2+cosr*x3
  return r2,r3

# From resampled (x2,x3) to survey (x2,x3)
def toSurveyFromResampled(r2,r3):
  x2 =  cosr*r2+sinr*r3
  x3 = -sinr*r2+cosr*r3
  x2 += x2c
  x3 += x3c
  return x2,x3

# From resampled (x2,x3) to map (xe,yn)
def toMapFromResampled(x2,x3):
  x2,x3 = toSurveyFromResampled(x2,x3)
  return toMapFromSurvey(x2,x3)

# Reads the field boundary polygon
def readFieldBoundary():
  input = open(dataDir+"fieldBoundary.txt","r")
  lines = input.readlines()
  input.close()
  xlist = []
  ylist = []
  lines = lines[2:] # skip headings
  for line in lines:
    a = line.strip().split()
    x,y = a[0:2]
    x,y = ftToKm(float(x)),ftToKm(float(y))
    xlist.append(x)
    ylist.append(y)
  return xlist,ylist

def getWellColor(name):
  map = {
    "25-1-X-14"  :Color.MAGENTA,
    "67-1-TpX-10":Color.BLUE,
    "35-1-ShX-10":Color.CYAN,
    "48-X-28"    :Color.RED,
  }
  return map[name]

# Maps well names to well coordinates
def getWellCoordinates(name):
  map = {
    "25-1-X-14"  :(804965.5,948057.8),
    "67-1-TpX-10":(802582.6,951665.6),
    "35-1-ShX-10":(800350.0,953112.9),
    "48-X-28"    :(795788.4,966905.1)
  }
  xe,yn = map[name]
  xe = ftToKm(xe)
  yn = ftToKm(yn)
  return xe,yn

# Maps well names to datum elevations
def getWellElevations(name):
  map = {
    "25-1-X-14"  :5225.0,
    "67-1-TpX-10":5205.7,
    "35-1-ShX-10":5238.0,
    "48-X-28"    :5114.6
  }
  return ftToKm(map[name])

# Maps well names to time-depth lists
def getWellTimeDepth(name):
  map = { 
    "25-1-X-14":[
    (0.241, 2070), (0.372, 2779), (0.441, 3204), (0.479, 3430), (0.527, 3702),
    (0.570, 3936), (0.598, 4083), (0.613, 4165), (0.628, 4257), (0.701, 4717),
    (0.716, 4811), (0.796, 5438), (0.832, 5739), (0.838, 5792), (0.848, 5871),
    (0.850, 5890)],
    "67-1-TpX-10":[
    (0.230, 1790), (0.350, 2442), (0.415, 2852), (0.454, 3090), (0.502, 3364),
    (0.545, 3598), (0.572, 3740), (0.586, 3821), (0.618, 4028), (0.669, 4350),
    (0.682, 4432), (0.760, 5042), (0.795, 5346), (0.801, 5394), (0.804, 5425),
    (0.815, 5512)],
    "35-1-ShX-10":[
    (0.171,  663), (0.182,  724), (0.196,  786), (0.207,  847), (0.219,  909),
    (0.231,  970), (0.243, 1031), (0.254, 1094), (0.267, 1155), (0.280, 1217),
    (0.293, 1278), (0.303, 1340), (0.316, 1402), (0.327, 1463), (0.339, 1525),
    (0.350, 1586), (0.361, 1648), (0.373, 1709), (0.385, 1771), (0.397, 1832),
    (0.408, 1894), (0.420, 1956), (0.434, 2017), (0.444, 2079), (0.456, 2140),
    (0.467, 2202), (0.479, 2263), (0.492, 2325)],
    "48-X-28":[
    (0.000,    0), (0.132,  453), (0.135,  467), (0.138,  479), (0.176,  682),
    (0.526, 2669), (0.558, 2869), (0.563, 2906), (0.580, 3004), (0.600, 3092),
    (0.643, 3302), (0.699, 3613), (0.725, 3754), (0.734, 3805), (0.742, 3851),
    (0.743, 3858), (0.744, 3863), (0.745, 3872), (0.767, 3981), (0.780, 4064),
    (0.815, 4266), (0.835, 4372), (0.842, 4413), (0.851, 4471), (0.854, 4494),
    (0.873, 4634), (0.973, 5423), (0.983, 5505), (0.993, 5576), (0.993, 5580),
    (1.002, 5649), (1.016, 5754), (1.300, 8000)]
  }
  table = map[name]
  tlist,zlist = [],[]
  for t,z in table:
    tlist.append(t)
    zlist.append(ftToKm(z))
  return tlist,zlist

def plotTimeDepth():
  sp = SimplePlot()
  sp.setHLabel("time (s)")
  sp.setVLabel("depth (km)")
  for wellName in ["25-1-X-14","67-1-TpX-10","48-X-28","35-1-ShX-10"]:
    tlist,zlist = getWellTimeDepth(wellName)
    pv = sp.addPoints(tlist,zlist)
    pv.setLineColor(getWellColor(wellName))
    pv.setMarkColor(getWellColor(wellName))
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

def plotCoordinates():
  xe = [xll,xlr,xur,xul,xll]
  yn = [yll,ylr,yur,yul,yll]
  sp = SimplePlot()
  sp.setSize(622,875)
  sp.setHLabel("Easting (km)")
  sp.setVLabel("Northing (km)")
  pv = sp.addPoints(xe,yn)
  pv.setLineColor(Color.BLUE)
  e2 = f2+(n2-1)*d2
  e3 = f3+(n3-1)*d3
  x2 = [f2,e2,e2,f2,f2]
  x3 = [f3,f3,e3,e3,f3]
  xe,yn = [],[]
  for i in range(5):
    #xei,yni = toMapFromResampled(x2[i],x3[i])
    m = Coordinates.Map(Coordinates.Resampled(x2[i],x3[i]))
    xei,yni = ftToKm(m.xe),ftToKm(m.yn)
    xe.append(xei)
    yn.append(yni)
  pv = sp.addPoints(xe,yn)
  pv.setLineColor(Color.RED)
  xe,yn = readFieldBoundary()
  pv = sp.addPoints(xe,yn)
  pv.setLineColor(Color.GREEN)
  for wellName in ["25-1-X-14","67-1-TpX-10","35-1-ShX-10","48-X-28"]:
    x,y = getWellCoordinates(wellName)
    pv = sp.addPoints([x],[y])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(getWellColor(wellName))

#print "toSurveyFromMap(xll,yll)",toSurveyFromMap(xll,yll)
#print "toSurveyFromMap(xlr,ylr)",toSurveyFromMap(xlr,ylr)
#print "toSurveyFromMap(xur,yur)",toSurveyFromMap(xur,yur)
#print "toSurveyFromMap(xul,yul)",toSurveyFromMap(xul,yul)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
