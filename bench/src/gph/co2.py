#############################################################################
# Demo interpolation on earth's surface using atmospheric CO2 measurements

from shared import *

def main(args):
  printStats()
  #goCo2Table(2009,2010,12)
  #goPlotCo2Table2()
  #printCo2Table()
  #goPlotEarth2()
  #goPlotEarth3()
  #goPlotBoard3()
  #goDemo()

def goCo2Table(fy,ly,ky):
  eimage = readWaterImage()
  names,lats,lons,co2s = readCo2Table()
  sy = yearsSamplingCo2Table()
  ny = sy.count
  iyf = max(sy.indexOfNearest(fy)-1,0)
  iyl = min(sy.indexOfNearest(ly)+1,ny-1)
  while sy.getValue(iyf)<fy:
    iyf += 1
  while sy.getValue(iyl)>ly:
    iyl -= 1
  for iy in range(iyf,iyl+1,ky):
    y = sy.getValue(iy)
    laty,lony,co2y = getGoodData(iy,lats,lons,co2s)
    if not laty:
      continue
    for scale in [1.0,4.0]:
      f,u1,u2 = co2y,lony,laty
      g,s1,s2 = gridBlended(scale,f,u1,u2)
      d = makeTensors(scale,s1,s2)
      title = strMonthYear(y)
      plot2(eimage,f,u1,u2,g,s1,s2,d,
            mv=True,cv=True,tv=False,
            year=y,title=title,png="co2U")
      if iy+ky>iyl:
        cmin,cmax = clipsForYear(y)
        cmap = ColorMap(cmin,cmax,ColorMap.JET)
        cimage = SampledImage.fromFloats(s1,s2,g,cmap)
        cimage.setAlpha(0.8)
        plot3(eimage,cimage)

def goPlotEarth2():
  image = readEarthImage()
  plot2(image,png="earthU")

def goPlotEarth3():
  eimage = readEarthImage()
  plot3(eimage,None)

def goPlotBoard3():
  eimage = makeBoardImage()
  plot3(eimage,None)

def goDemo():
  frame = SimpleFrame(AxesOrientation.XOUT_YRIGHT_ZUP)
  world = frame.getWorld()
  group = makeDemoEarthGroup()
  world.addChild(group)
  view = frame.getOrbitView()
  view.setAzimuth(-90.0)
  view.setWorldSphere(world.getBoundingSphere(True))

def makeDemoEarthGroup():
  model = OblateModel.forWGS84()
  group = OblateImageGroup(model)
  image = readEarthImage()
  #image = makeGrayImage(image)
  group.addImage(image,0.999)
  nlon,nlat = 21,41
  slon = Sampling(nlon,1.0,-120.0)
  slat = Sampling(nlat,1.0,  20.0)
  r = fillbyte(255,nlon,nlat)
  g = zerobyte(nlon,nlat)
  b = zerobyte(nlon,nlat)
  a = fillbyte(100,nlon,nlat)
  image = SampledImage(slon,slat,r,g,b,a)
  group.addImage(image,1.000)
  group.addGrid(Color.BLACK,1,1.001)
  glon = Sampling(11,2.0,-120.0)
  glat = Sampling(11,4.0,  20.0)
  group.addGrid(glon,glat,Color.WHITE,1,1.002)
  return group

#############################################################################
# Plot

pngDir = None
#pngDir = "png/"
backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent
#cmin,cmax = 374.0,379.0 # clips for fields CO2 plots
#cmin,cmax = 320.0,400.0 # clips for ESRL table CO2 plots
#cmin,cmax = 378.0,402.0 # clips for ESRL table CO2 plots, 2009

def plot2(image,
          f=None,u1=None,u2=None,
          g=None,s1=None,s2=None,d=None,
          mv=False,cv=False,tv=False,
          year=0,axes=False,title=None,png=None):
  ppo = PlotPanel.Orientation.X1RIGHT_X2UP
  if axes: ppa = PlotPanel.AxesPlacement.LEFT_BOTTOM
  else: ppa = PlotPanel.AxesPlacement.NONE
  pp = PlotPanel(1,1,ppo,ppa)
  pp.setHLimits(-181,181)
  pp.setVLimits( -91, 91)
  tile = pp.getTile(0,0)
  slon,slat = image.samplingX,image.samplingY
  rgba = image.getFloatsRGBA(False)
  pv = PixelsView(slon,slat,rgba)
  pv.setClips(0.0,1.0)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  tile.addTiledView(pv)
  gwidth = 3
  if cv and g and s1 and s2:
    if year!=0:
      cmin,cmax = clipsForYear(year)
    else:
      cmin,cmax = min(g),max(g)
    gwidth = 1
    cv = PixelsView(s1,s2,g)
    cv.setInterpolation(PixelsView.Interpolation.NEAREST)
    #cmin,cmax = min(f),max(f)
    #cmin -= 1.1*(cmax-cmin)/256 # make null values transparent
    #cv.setColorModel(makeTransparentColorModel(0.8))
    cv.setColorModel(ColorMap.getJet(0.8))
    cv.setClips(cmin,cmax)
    tile.addTiledView(cv)
    cb = pp.addColorBar("CO2 (ppm)")
    cb.setInterval(5)
  if tv and d and s1 and s2:
    gwidth = 1
    tv = TensorsView(s1,s2,d)
    tv.setScale(22.0)
    tv.setLineColor(Color.BLACK)
    tv.setLineWidth(2)
    tv.setEllipsesDisplayed(11)
    tile.addTiledView(tv)
  if mv and f and u1 and u2:
    gwidth = 1
    mv = PointsView(u1,u2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkColor(Color.WHITE)
    mv.setMarkSize(8)
    tile.addTiledView(mv)
  gv = gridLines2(Sampling(13,30,-180),Sampling(7,30,-90),Color.BLACK,gwidth)
  tile.addTiledView(gv)
  pf = PlotFrame(pp)
  pf.setBackground(backgroundColor)
  pf.setFontSizeForPrint(8,240)
  if axes:
    pp.setHLabel("Longitude (degrees)")
    pp.setVLabel("Latitude (degrees)")
    pp.setVInterval(60)
    pp.setHInterval(60)
    pf.setSize(1085,615)
  else:
    if title and cv:
      pp.setTitle(title)
      pf.setSize(1740,890)
    elif cv:
      pf.setSize(1740,770)
    else: 
      pf.setSize(1078,562)
  pf.setVisible(True)
  if png and pngDir: 
    pf.paintToPng(600,3,pngDir+png+".png")

def plot3(eimage,cimage=None):
  model = OblateModel.forWGS84()
  group = OblateImageGroup(model)
  group.addImage(eimage,0.999)
  gwidth = 3
  if cimage:
    group.addImage(cimage,1.000)
    gwidth = 1
  group.addGrid(Color.BLACK,gwidth,1.001)
  frame = SimpleFrame(AxesOrientation.XOUT_YRIGHT_ZUP)
  frame.viewCanvas.setBackground(backgroundColor)
  frame.setSize(1000,1000)
  world = frame.getWorld()
  world.addChild(group)
  view = frame.getOrbitView()
  #view.setAzimuth(-90.0) # good for US
  view.setAzimuth(45.0) # good for showing axes
  view.setWorldSphere(world.getBoundingSphere(True))

def gridLines2(slon,slat,color,width):
  nlon,nlat = slon.count,slat.count
  flon,flat = slon.first,slat.first
  llon,llat = slon.last,slat.last
  x = zerofloat(2,nlon+nlat)
  y = zerofloat(2,nlon+nlat)
  k = 0
  for ilon in range(nlon):
    loni = slon.getValue(ilon)
    x[k][0],x[k][1] = loni,loni
    y[k][0],y[k][1] = flat,llat
    k += 1
  for ilat in range(nlat):
    lati = slat.getValue(ilat)
    x[k][0],x[k][1] = flon,llon
    y[k][0],y[k][1] = lati,lati
    k += 1
  pv = PointsView(x,y)
  pv.setLineColor(color)
  pv.setLineWidth(width)
  return pv

def makeGrayImage(image):
  f = image.getFloatsGray(False)
  sx = image.getSamplingX()
  sy = image.getSamplingY()
  cmap = ColorMap(0.0,1.0,ColorMap.GRAY)
  return SampledImage.fromFloats(sx,sy,f,cmap)

def makeTransparentColorModel(alpha):
  cm = ColorMap.getJet(alpha)
  r = zerobyte(256)
  g = zerobyte(256)
  b = zerobyte(256)
  a = zerobyte(256)
  for i in range(256):
    r[i] = byteFromInt(cm.getRed(i))
    g[i] = byteFromInt(cm.getGreen(i))
    b[i] = byteFromInt(cm.getBlue(i))
    a[i] = byteFromInt(cm.getAlpha(i))
  a[0] = 0
  return IndexColorModel(8,256,r,g,b,a)
def byteFromInt(i):
  if i>127: i -= 256
  return i

#############################################################################
# Test (checkerboard) and background image data

imageDir = "/Users/dhale/Home/box/idh/trunk/bench/src/gph/resources/"

def makeBoardImage(): # checkerboard image for testing
  nlon,nlat = 12,7 # note 12, not 13; need not sample both -180 and 180
  slon = Sampling(nlon,30.0,-180.0)
  slat = Sampling(nlat,30.0, -90.0)
  r = zerobyte(nlon,nlat)
  g = zerobyte(nlon,nlat)
  b = zerobyte(nlon,nlat)
  a = zerobyte(nlon,nlat)
  for ilat in range(nlat):
    for ilon in range(nlon):
      r[ilat][ilon] = -1*((ilat+ilon)%2)
      g[ilat][ilon] = 0
      b[ilat][ilon] = -1*((ilat+ilon+1)%2)
      a[ilat][ilon] = -1
  return SampledImage(slon,slat,r,g,b,a)

def readEarthImage():
  slon = Sampling(2048,360.0/2047,-180.0)
  slat = Sampling(1024,180.0/1023, -90.0)
  image = SampledImage.fromFile(slon,slat,imageDir+"earth.png")
  return image

def readWaterImage():
  slon = Sampling(2048,360.0/2047,-180.0)
  slat = Sampling(1024,180.0/1023, -90.0)
  image = SampledImage.fromFile(slon,slat,imageDir+"water.png")
  sx = image.getSamplingX()
  sy = image.getSamplingY()
  grays = image.getFloatsGray(False)
  grays = clip(0.7,1.0,grays)
  cmap = ColorMap(0.0,1.0,ColorMap.GRAY)
  return SampledImage.fromFloats(sx,sy,grays,cmap)

#############################################################################
# CO2 surface flask data from NOAA ESRL GMD

dataDir = "/Users/dhale/Home/box/idh/trunk/bench/src/gph/co2/"

def readCo2Table():
  s = Scanner(FileInputStream(dataDir+"co2table.txt"))
  names,lats,lons,co2s = [],[],[],[]
  while s.hasNextLine():
    line = s.nextLine()
    if line[0]=="#":
      continue
    data = line.split()
    if len(data)==0:
      continue
    elif data[0]=="station":
      names = data[1:]
      ns = len(names)
    elif data[0]=="latitude":
      for js in range(ns):
        lats.append(float(data[js+1]))
    elif data[0]=="longitude":
      for js in range(ns):
        lons.append(float(data[js+1]))
    else:
      co2y = []
      for js in range(ns):
        co2y.append(float(data[js+1]))
      co2s.append(co2y)
  s.close()
  return names,lats,lons,co2s

def getGoodData(iy,lats,lons,co2s):
  laty,lony,co2y = [],[],[]
  co2i = co2s[iy]
  ns = len(co2i)
  for js in range(ns):
    if co2i[js]>0.0:
      laty.append(lats[js])
      lony.append(lons[js])
      co2y.append(co2i[js])
  if len(laty)>0:
    return laty,lony,co2y
  else:
    return None,None,None

co2TableHeader = """###
# Table of monthly averages of atmospheric CO2 concentations (ppm mol ratios) 
# measured in individual surface flasks at stations around the world. 
# These data are freely available from ftp://ftp.cmdl.noaa.gov/cgg.
#
# For simplicity the data format is a simple tab-delimited table,
# with leading space characters included to align table columns.
# The first row contains "station", then station codes, such as MLO.
# The second row contains "latitude", then station latitudes, in degrees.
# The third row contains "longitude", then station longitudes, in degrees.
# The remaining rows contain a decimal year followed by co2 ppm values, 
# one value per station. Missing values are zero.
#
# The number of samples (rows with year and co2 values) is 600. 
# The sampling interval is 1.0/12.0 (monthly).
# The year of first sample is 1960.04 (middle of January).
# The year of last sample is 2009.96 (middle of December).
# For some early years, values for all stations are missing.
# This table includes only the 84 stations with at least one 
# measured co2 value in the range of years sampled.
#
# Compiled by Dave Hale, Colorado School of Mines, 2001.01.02
###"""

def yearsSamplingCo2Table():
  dy = 1.0/12.0 # sampling interval, in years
  fy = 1960+ 0.5*dy # first year
  ly = 2009+11.5*dy # last year
  ny = int((ly-fy)/dy+1.5) # number of years
  return Sampling(ny,dy,fy)

def readCo2Yearly():
  fileName = "co2yearly.txt"
  years,co2s = [],[]
  s = Scanner(FileInputStream(dataDir+fileName))
  while s.hasNextLine():
    line = s.nextLine()
    if line[0]!="#":
      data = line.split()
      years.append(float(data[0]))
      co2s.append(float(data[1]))
  s.close()
  return years,co2s

def printCo2Table():
  print co2TableHeader
  sy = yearsSamplingCo2Table()
  ny = sy.count
  sll = readCo2Stations() # map: station name -> lat,lon
  scs = {} # map: station name -> array of co2s
  for name in sll:
    co2s = readCo2s(name,sy)
    if co2s:
      scs[name] = co2s
  names = scs.keys() # names of stations with co2s
  names.sort() # sorted alphabetically
  ns = len(names)
  printf("%16s","station") # station names (currently 84)
  for js in range(ns):
    name = names[js]
    printf("\t%8s",name)
  printf("\n")
  printf("%16s","latitude") # latitudes
  for js in range(ns):
    name = names[js]
    lat,lon = sll[name]
    printf("\t%8.2f",lat)
  printf("\n")
  printf("%16s","longitude") # longitudes
  for js in range(ns):
    name = names[js]
    lat,lon = sll[name]
    printf("\t%8.2f",lon)
  printf("\n")
  for jy in range(ny): # years and co2s
    yj = sy.getValue(jy)
    printf("%16.2f",yj)
    for js in range(ns):
      name = names[js]
      co2 = scs[name][jy]
      printf("\t%8.2f",co2)
    printf("\n")

def readCo2Stations():
  fileName = "co2stations.txt"
  sll = {}
  s = Scanner(FileInputStream(dataDir+fileName))
  s.useDelimiter("\t")
  s.nextLine(); s.nextLine(); s.nextLine()
  while s.hasNextLine():
    line = s.nextLine()
    data = line.split("\t")
    station = data[0]
    temp = String(data[0])
    if temp.endsWith(" *"):
      station = station[:-2]
    lat = float(data[2])
    lon = float(data[3])
    sll[station] = (lat,lon)
  s.close()
  return sll

def readCo2s(station,sy):
  ny,dy,fy,ly = sy.count,sy.delta,sy.first,sy.last
  co2s = zerofloat(ny)
  station = str.upper(station)
  fileName = "month/"+str.lower(station)
  if station[:3]=="POC" or station[:3]=="SCS":
    fileName += "_01D1_mm.co2"
    station = station[3:]
  else:
    fileName += "_01D0_mm.co2"
  nco2 = 0
  file = File(dataDir+fileName)
  if file.exists():
    s = Scanner(FileInputStream(file))
    while s.hasNextLine():
      line = s.nextLine()
      data = line.split()
      if len(data)<4:
        continue
      name = data[0]
      if name==station:
        year = float(data[1])
        month = float(data[2])
        y = year+(month-0.5)*dy
        if y<fy or y>ly:
          print "station:",station," year",y,"is out of range"
          continue
        iy = sy.indexOfNearest(y)
        co2s[iy] = float(data[3])
        nco2 += 1
    s.close()
  if nco2>0: return co2s
  else: return None

#############################################################################
# Other stuff

def strMonthYear(y):
  mmap = {1:"Jan", 2:"Feb", 3:"Mar", 4:"Apr", 5:"May", 6:"Jun",
          7:"Jul", 8:"Aug", 9:"Sep",10:"Oct",11:"Nov",12:"Dec"}
  year = int(y)
  month = mmap[int(1.0+(y-year)*12.0)]
  return month+", "+str(year)

co2YearlyAverage = None
def clipsForYear(y):
  global co2YearlyAverage
  if not co2YearlyAverage:
    years,co2s = readCo2Yearly()
    n = len(years)
    method = CubicInterpolator.Method.LINEAR
    co2YearlyAverage = CubicInterpolator(method,n,years,co2s)
  co2 = co2YearlyAverage.interpolate(y)
  cmin = co2-10.0
  cmax = co2+10.0
  return cmin,cmax

def makeTensors(scale,slon,slat):
  nlon,nlat = slon.count,slat.count
  au = zerofloat(nlon,nlat)
  av = fillfloat(1.0,nlon,nlat)
  u1 = fillfloat(1.0,nlon,nlat)
  u2 = fillfloat(0.0,nlon,nlat)
  for ilat in range(nlat):
    lati = slat.getValue(ilat)
    cosl = cos(lati*PI/180)
    auc = scale*scale/(cosl*cosl)
    for ilon in range(nlon):
      au[ilat][ilon] = auc
  return EigenTensors2(u1,u2,au,av)

def gridSamplings():
  #n1,n2 = 360,180
  n1,n2 = 180,90
  d1,d2 = 360.0/n1,180.0/n2
  f1,f2 = -180.0+0.5*d1,-90.0+0.5*d2
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  return s1,s2

def gridSimple(f,x1,x2):
  s1,s2 = gridSamplings()
  sg = SimpleGridder2(f,x1,x2)
  g = sg.grid(s1,s2)
  return g,s1,s2

def gridBlended(scale,f,x1,x2):
  fp,x1p,x2p = padLongitude(f,x1,x2)
  s1,s2 = gridSamplings()
  n1,f1 = s1.count,s1.first
  n1p,f1p = 2*n1,f1-180.0
  s1p = Sampling(n1p,d1,f1p)
  n1p,n2,n1p = s1p.count,s2.count
  d = makeTensors(scale,s1p,s2)
  bg = BlendedGridder2(d,fp,x1p,x2p)
  #bg.setBlending(False)
  gp = bg.grid(s1p,s2)
  g = copy(n1,n2,n1/2,0,gp)
  return g,s1,s2

def padLongitude(f,x1,x2):
  n = len(f)
  for i in range(n):
    if x1[i]<0.0:
      x1.append(x1[i]+360.0)
    else:
      x1.append(x1[i]-360.0)
    f.append(f[i])
    x2.append(x2[i])
  return f,x1,x2

def wrapLon(f,sx,sy):
  nx,ny = sx.count,sy.count
  g = zerofloat(nx+1,ny)
  copy(nx,ny,f,g)
  copy(1,ny,0,0,f,nx,0,g)
  sx = Sampling(nx+1,sx.delta,sx.first)
  return g,sx,sy

def printf(fmt,*varargs):
  sys.stdout.write(fmt % varargs)

def printStats():
  names,lats,lons,co2s = readCo2Table()
  sy = yearsSamplingCo2Table()
  ny = sy.count
  slat = Sampling(180,1.0, -79.5)
  slon = Sampling(360,1.0,-179.5)
  nlat,nlon = slat.count,slon.count
  alat,alon = zerofloat(nlat),zerofloat(nlon)
  clat,clon = zerofloat(nlat),zerofloat(nlon)
  for iy in range(ny-1,ny):
    laty,lony,co2y = getGoodData(iy,lats,lons,co2s)
    if not laty:
      continue
    ns = len(laty)
    for js in range(ns):
      latj = laty[js]
      lonj = lony[js]
      co2j = co2y[js]
      jlat = slat.indexOfNearest(latj)
      jlon = slon.indexOfNearest(lonj)
      alat[jlat] += co2j
      clat[jlat] += 1.0
      alon[jlon] += co2j
      clon[jlon] += 1.0
  SimplePlot.asSequence(slat,alat)
  SimplePlot.asSequence(slon,alon)
  mlat = sum(alat)/sum(clat)
  mlon = sum(alon)/sum(clon)
  print "mlat =",mlat," mlon =",mlon
  vlat = 0.0
  for jlat in range(nlat):
    if clat[jlat]>0:
      dlat = alat[jlat]/clat[jlat]-mlat
      vlat += dlat*dlat
  vlat /= nlat
  vlon = 0.0
  for jlon in range(nlon):
    if clon[jlon]>0:
      dlon = alon[jlon]/clon[jlon]-mlon
      vlon += dlon*dlon
  vlon /= nlon
  print "vlat =",vlat," vlon =",vlon


      

"""
count co2(lat) for 10-degree lat bins
count co2(lon) for 10-degree lon bins
compute avglat = sum(co2(lat))/sum(lat)
compute varlat = sum((co2(lat)-avglat)^2)/sum(lat)
"""


  
  
  

#############################################################################
if __name__ == "__main__":
  run(main)
