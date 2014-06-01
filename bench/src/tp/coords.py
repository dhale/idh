"""
Shows DOE, Seismic and CSM coordinate rectangles for Teapot Dome data.
Author: Dave Hale, Colorado School of Mines
Version: 2009.07.25
"""
from imports import *

def main(args):
  display()

# Files
tpDir = "/data/seis/tp/"
stFile = tpDir+"doe/3D_Seismic/tpstAll.dat" # seismic time image
fbFile = tpDir+"doe/CD/fieldBoundary.txt" # boundary polygon

# Map coordinates (xe=easting,yn=northing) of DOE seismic rectangle.
# Coordinates measured in ft, obtained from teapot_3d_load.doc.
xll,yll = 788937.0,938846.0 # lower left  <=> (i2,i3) = (  0,  0)
xlr,ylr = 809502.0,939334.0 # lower right <=> (i2,i3) = (187,  0)
xur,yur = 808604.0,977163.0 # upper right <=> (i2,i3) = (187,344)
xul,yul = 788039.0,976675.0 # upper left  <=> (i2,i3) = (  0,344)

# Convert map coordinates from ft to km.
def ftToKm(ft):
  return ft*0.0003048
def kmToFt(km):
  return km*3280.84
xll,yll = ftToKm(xll),ftToKm(yll)
xlr,ylr = ftToKm(xlr),ftToKm(ylr)
xur,yur = ftToKm(xur),ftToKm(yur)
xul,yul = ftToKm(xul),ftToKm(yul)

# CSM sampling of 3D seismic time images
s1c = Sampling(1501,0.002,0.000)
s2c = Sampling(357,0.025,0.000)
s3c = Sampling(161,0.025,0.000)

# DOE sampling of 3D seismic time images
s1d = Sampling(1501,0.002,0.000)
s2d = Sampling(188,0.033528,0.000) # 0.033528 km = 110 ft
s3d = Sampling(345,0.033528,0.000)

# Reads the field boundary polygon
def readFieldBoundary():
  input = open(fbFile,"r")
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

def display():
  sp = SimplePlot()
  #sp.setSize(623,875)
  #sp.setFontSizeForPrint(8,240)
  sp.setSize(621,826)
  sp.setFontSizeForSlide(1.0,1.0)
  sp.setHLabel("Easting (km)")
  sp.setVLabel("Northing (km)")
  sp.plotPanel.setHInterval(2)
  sp.plotPanel.setVInterval(2)
  # DOE seismic image
  s2,s3,ss = makeSeismicSlice(475)
  pv = sp.addPixels(s2,s3,ss)
  # field boundary
  xe,yn = readFieldBoundary()
  pv = sp.addPoints(xe,yn)
  pv.setLineColor(Color.GREEN)
  pv.setLineWidth(3.0)
  pv.setLineStyle(PointsView.Line.DOT)
  #sp.paintToPng(600,3.33,"coords1.png")
  # DOE seismic rectangle
  xe = [xll,xlr,xur,xul,xll]
  yn = [yll,ylr,yur,yul,yll]
  pv = sp.addPoints(xe,yn)
  pv.setLineColor(Color.BLUE)
  pv.setLineWidth(5.0)
  #pv.setLineWidth(3.0)
  #pv.setLineStyle(PointsView.Line.DASH)
  #sp.paintToPng(600,3.33,"coords2.png")
  # CSM seismic rectangle
  n2,d2,f2 = s2c.count,s2c.delta,s2c.first
  n3,d3,f3 = s3c.count,s3c.delta,s3c.first
  e2 = f2+(n2-1)*d2
  e3 = f3+(n3-1)*d3
  x2 = [f2,e2,e2,f2,f2]
  x3 = [f3,f3,e3,e3,f3]
  xe,yn = [],[]
  for i in range(5):
    m = Coordinates.Map(Coordinates.Csm(x2[i],x3[i]))
    xei,yni = ftToKm(m.xe),ftToKm(m.yn)
    xe.append(xei)
    yn.append(yni)
  pv = sp.addPoints(xe,yn)
  pv.setLineColor(Color.RED)
  pv.setLineWidth(5.0)
  sp.paintToPng(600,3.33,"coords3.png")

def makeSeismicSlice(i1):
  # DOE seismic image sampling
  n1d,d1d,f1d = s1d.count,s1d.delta,s1d.first
  n2d,d2d,f2d = s2d.count,s2d.delta,s2d.first
  n3d,d3d,f3d = s3d.count,s3d.delta,s3d.first
  ais = ArrayInputStream(stFile)
  x1 = zerofloat(n1d)
  x23 = zerofloat(n2d,n3d)
  for i3 in range(n3d):
    for i2 in range(n2d):
      ais.readFloats(x1)
      x23[i3][i2] = x1[i1]
  ais.close()
  xmin,xmax = min(x23),max(x23)
  print "xmin =",xmin," xmax =",xmax
  f2m,e2m,d2m = 240.0,248.2,0.025 # chosen by eye to fill the
  f3m,e3m,d3m = 286.0,298.0,0.025 # map-coordinate rectangle 
  n2m = 1+int((e2m-f2m)/d2m+0.5)
  n3m = 1+int((e3m-f3m)/d3m+0.5)
  s2m,s3m = Sampling(n2m,d2m,f2m),Sampling(n3m,d3m,f3m)
  y23 = zerofloat(n2m,n3m)
  si = SincInterpolator()
  e2d,e3d = f2d+(n2d-1)*d2d,f3d+(n3d-1)*d3d
  for i3 in range(n3m):
    x3m = f3m+i3*d3m
    x3m = kmToFt(x3m)
    for i2 in range(n2m):
      x2m = f2m+i2*d2m
      x2m = kmToFt(x2m)
      doe = Coordinates.Doe(Coordinates.Map(x2m,x3m))
      x2d,x3d = doe.x2,doe.x3
      if x2d<f2d or x2d>e2d or x3d<f3d or x3d>e3d:
        y23[i3][i2] = xmax
      else:
        y23[i3][i2] = si.interpolate(n2d,d2d,f2d,n3d,d3d,f3d,x23,x2d,x3d)
  ymin,ymax = min(y23),max(y23)
  print "ymin =",ymin," ymax =",ymax
  return s2m,s3m,y23

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
