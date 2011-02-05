#############################################################################
# Reads, processes and plots elevation data from the GTOPO30 database: 
# http://eros.usgs.gov/#/Find_Data/Products_and_Data_Available/gtopo30_info 
#
# To use this program, download at least one tile from the GTOPO30 database.
# The files for the tile should be in directory with a name like "e100n40". 
# Modify the topoDir directory in the data input/output section below.
#############################################################################
from shared import *

def main(args):
  goTaiwanShow()
  #goTaiwanTopo()
  #goTaiwanSlope()

def goTaiwanShow():
  s1,s2 = getTaiwanSamplings()
  z = readImage("TaiwanTopo",s1,s2)
  s = readImage("TaiwanSlope",s1,s2)
  print "Taiwan:"
  print "  nlon =",s1.count," dlon =",s1.delta," flon =",s1.first
  print "  nlat =",s2.count," dlat =",s2.delta," flat =",s2.first
  z = log10(add(0.0001,z))
  s = log10(add(0.0001,s))
  plot(z,s1,s2,
       cbar="Log10[Elevation (m)]",cmin=0,cmax=3.5,
       contours=0,htic=1.0,vtic=1.0,width=825,height=900)
  plot(s,s1,s2,
       cbar="Log10[Slope (m/m)]",cmin=0,cmax=0,
       contours=0,htic=1.0,vtic=1.0,width=825,height=900)

def goTaiwanTopo():
  s1,s2 = getTaiwanSamplings()
  n1,n2 = s1.count,s2.count
  f1,f2 = s1.first,s2.first
  z,s1,s2 = readTile("e100n40") # elev z, lon sampling s1, lat sampling s2
  plot(z,s1,s2,width=670,height=830)
  j1,j2 = s1.indexOfNearest(f1),s2.indexOfNearest(f2)
  z,s1,s2 = subsetImage(z,s1,s2,n1,n2,j1,j2)
  writeImage("TaiwanTopo",z)
  print "Taiwan:"
  print "  nlon =",s1.count," dlon =",s1.delta," flon =",s1.first
  print "  nlat =",s2.count," dlat =",s2.delta," flat =",s2.first
  z = log10(add(0.0001,z))
  plot(z,s1,s2,contours=True,htic=1.0,vtic=1.0,width=670,height=900)
  plot(z,s1,s2,
       cbar="Log10[Elevation (m)]",cmin=0,cmax=3.5,
       contours=10,htic=1.0,vtic=1.0,width=825,height=900)

def goTaiwanSlope():
  s1,s2 = getTaiwanSamplings()
  z = readImage("TaiwanTopo",s1,s2)
  m = makeMask(z)
  s = computeSlopes(z,s1,s2)
  s = applyMask(s,m)
  writeImage("TaiwanSlope",s)
  print "Taiwan:"
  print "  elevation (m): min =",min(z)," max =",max(z)
  print "    slope (m/m): min =",min(s)," max =",max(s)
  z = log10(add(0.1,z)) # warp to see lowlands
  s = log10(add(0.0001,s)) # warp to see lowlands
  plot(z,s1,s2,cbar="Log10[elevation (m)]",
       htic=1.0,vtic=1.0,width=840,height=910)
  plot(s,s1,s2,cbar="Log10[slope (m/m)]",
       htic=1.0,vtic=1.0,width=840,height=910)

#############################################################################
# processing

def computeSlopes(z,s1,s2):
  """
  Returns slopes in m/m for specified elevations z.
  s1 and s2 are samplings of longitude and latitude.
  m is an optional land mask of 1.0 for land and 0.0 for water.
  """
  n1,n2 = s1.count,s2.count # number of samples for lon,lat
  d1,d2 = s1.delta,s2.delta # sampling intervals for lon,lat
  d1,d2 = toRadians(d1),toRadians(d2) # in radians
  r = 6378137.0 # equitorial radius of earth, in m; assume a sphere
  s = zerofloat(n1,n2) # array of slopes, initially zero
  for i2 in range(1,n2-1): # for all latitudes, ...
    lat = toRadians(s2.getValue(i2)) # latitude, in radians
    r1 = 0.5/(d1*r*cos(lat)) # scale factor for derivative 1
    r2 = 0.5/(d2*r         ) # scale factor for derivative 2
    for i1 in range(1,n1-1): # for all longitudes, ...
      z1 = r1*(z[i2][i1+1]-z[i2][i1-1]) # derivative 1 in m/m
      z2 = r2*(z[i2+1][i1]-z[i2-1][i1]) # derivative 2 in m/m
      s[i2][i1] = sqrt(z1*z1+z2*z2) # slope = magnitude of gradient
  s = add(FLT_MIN,s) # ensure no zero slopes; zeros are water
  return s

def makeMask(z):
  """Returns a land mask, 1.0 for land, 0.0 for water."""
  znull = min(z)
  n1,n2 = len(z[0]),len(z)
  m = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if z[i2][i1]>znull:
        m[i2][i1] = 1.0
  return m

def applyMask(z,m):
  """Applies the specified mask, so that values for water are 0.0."""
  return mul(z,m)

#############################################################################
# plotting

#pngDir = "png/vs30/" # where to put PNG images of plots
pngDir = None # for no PNG images

def plot(z,s1,s2,
         cbar=None,cmin=0,cmax=0,png=None,contours=0,
         htic=0,vtic=0,width=0,height=0):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setHLabel("Longitude (degrees)")
  sp.setVLabel("Latitude (degrees)")
  sp.plotPanel.setColorBarWidthMinimum(110)
  if htic>0: 
    sp.setHInterval(htic)
  if vtic>0: 
    sp.setVInterval(vtic)
  sp.setFontSizeForPrint(8,240)
  if width==0 and height==0:
    width = 900*s1.count/s2.count
    height = 900
  elif width==0:
    width = height*s1.count/s2.count
  elif height==0:
    height = width*s2.count/s1.count
  sp.setSize(width,height)
  if cbar:
    sp.addColorBar(cbar)
  pv = sp.addPixels(s1,s2,z)
  pv.setColorModel(ColorMap.JET)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if contours>0:
    cv = sp.addContours(s1,s2,z)
    cv.setContours(contours)
    cv.setLineColor(Color.BLACK)
  if pngDir and png:
    sp.paintToPng(600,3,pngDir+png+".png") # 600 dpi, 3 inches wide

#############################################################################
# data input/output

topoDir = "/data/earth/gtopo30/" # contains the GTOPO30 directories
imageDir = "/data/earth/vs30/" # directory in which to read/write images

def getTaiwanSamplings():
  """
  Returns samplings of longitude and latitude for the Taiwan subset.
  This subset is consistent with the map of Vs30 on the USGS website.
  """
  nlon,nlat = 312,418
  dlon,dlat = 0.5/60,0.5/60 # 0.0083333333333 degrees = 30 seconds
  flon,flat = 119.60416666,21.84583333
  return Sampling(nlon,dlon,flon),Sampling(nlat,dlat,flat)

def readTile(tileName):
  """
  Reads the gtopo30 tile with specified name; e.g., "e100n40".
  Returns the 2D array of elevations, lon sampling, and lat sampling.
  The first lat and lon correspond to the southwest grid corner.
  """
  fileRoot = topoDir+tileName+"/"+str.upper(tileName)
  fileName = fileRoot+".HDR" # the header file
  info = {} # a dictionary for information read from header file
  s = Scanner(FileInputStream(fileName))
  while s.hasNextLine():
    data = s.nextLine().split()
    info[data[0]] = data[1] # info[key] = value
  s.close()
  nlon = int(info["NCOLS"]) # number of lons
  nlat = int(info["NROWS"]) # number of lats
  dlon = float(info["XDIM"]) # lon sampling interval
  dlat = float(info["YDIM"]) # lat sampling interval
  flon = float(info["ULXMAP"]) # first lon sampled
  flat = float(info["ULYMAP"]) # first lat sampled
  flat -= (nlat-1)*dlat # because we will flip image vertically
  slon = Sampling(nlon,dlon,flon) # lon sampling
  slat = Sampling(nlat,dlat,flat) # lat sampling
  z = zeroshort(nlon,nlat) # elevations are shorts (16-bit integers)
  fileName = fileRoot+".DEM" # the elevation file
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  ais.readShorts(z)
  ais.close()
  z = Util.floatsFromShorts(z) # converts to floats, replaces -9999 with 0.
  z = flip2(z) # flip array so first sample is in lower-left corner
  return z,slon,slat

def subsetImage(z,s1,s2,n1,n2,j1,j2):
  """
  Returns a subset of the specified array z and samplings s1 and s2.
  n1 and n2 are output numbers of samples in 1st and 2nd dimensions.
  j1 and j2 are corresponding indices of first samples output.
  Returns a new array z with samplings s1 and s2 for the subset.
  """
  d1,d2 = s1.delta,s2.delta
  f1,f2 = s1.first,s2.first
  s1 = Sampling(n1,d1,f1+j1*d1)
  s2 = Sampling(n2,d2,f2+j2*d2)
  z = copy(n1,n2,j1,j2,z)
  return z,s1,s2

def readImage(fileName,s1,s2):
  """Reads and returns an image with specified samplings from a file."""
  n1,n2 = s1.count,s2.count
  z = zerofloat(n1,n2)
  ais = ArrayInputStream(imageDir+fileName+".dat")
  ais.readFloats(z)
  ais.close()
  return z

def writeImage(fileName,z):
  """Writes the specified image to a file."""
  aos = ArrayOutputStream(imageDir+fileName+".dat")
  aos.writeFloats(z)
  aos.close()

def flip1(z):
  """Returns a copy of the specified image, with 1st dimension flipped."""
  z = copy(z)
  n2 = len(z)
  for i2 in range(n2):
    z[i2] = reverse(z[i2])
  return z

def flip2(z):
  """Returns a copy of the specified image, with 2nd dimension flipped."""
  z = copy(z)
  n2 = len(z)
  for i2 in range(n2/2):
    j2 = n2-1-i2
    zi2 = z[i2]
    z[i2] = z[j2]
    z[j2] = zi2
  return z

#############################################################################
run(main)
