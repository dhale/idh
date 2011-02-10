#############################################################################
# Demo interpolation on earth's surface using atmospheric CO2 measurements

from shared import *

pngDir = None
#pngDir = "png/co2/"

def main(args):
  #goCrossValidate(2009.28)
  goAnalyzeErrors()
  #goVariogram(1971.95,1972,1)
  #goVariogram(2009.28,2010,1)
  #goVariogram(2009.28,2010,12,random=True)
  #goVariogram(2009.28,2010,12,random=False)
  #goVariogram(1970.28,2010,12)
  #goGridding(2009.28,2010,12,sg=True,p3=True)
  #goGridding(2009.28,2010,12,sg=False,p3=True)
  #goGridding(2009.28,2010,12,sg=False,p3=False)
  #goGridding(2009.28,2010,12,p3=True)
  #goGridding(1970,2010,1) # Jan-Dec
  #goGridding(1970.04,2010,12) # Jan
  #goGridding(1970.12,2010,12) # Feb
  #goGridding(1970.20,2010,12) # Mar
  #goGridding(1970.29,2010,12) # Apr
  #goGridding(1970.38,2010,12) # May
  #goGridding(1970.46,2010,12) # Jun
  #goGridding(1970.54,2010,12) # Jul
  #goGridding(1970.62,2010,12) # Aug
  #goGridding(1970.70,2010,12) # Sep
  #goGridding(1970.78,2010,12) # Oct
  #goGridding(1970.87,2010,12) # Nov
  #goGridding(1970.95,2010,12) # Dec
  #goPlotCo2Table2()
  #printCo2Table()
  #goPlotEarth2()
  #goPlotEarth3()
  #goTensors()
  #goPlotBoard3()
  #goDemo()

def goTensors():
  eimage = readEarthImage()
  s1,s2 = gridSamplings()
  for scale in [1,4]:
    d = makeTensors(scale,s1,s2)
    plot2(eimage,None,None,None,None,s1,s2,d,
          mv=False,cv=False,tv=True,png="earthte"+str(scale))

def goVariogram(fy,ly,ky,random=False):
  names,lats,lons,co2s = readCo2Table()
  sy = yearsSamplingCo2Table()
  for iy in yearIndicesCo2Table(fy,ly,ky):
    laty,lony,co2y = getGoodData(iy,lats,lons,co2s)
    n = len(laty)
    if random:
      co2y = permute(co2y)
    x,y,s = [],[],[]
    for i in range(n):
      lati,loni,co2i = laty[i],lony[i],co2y[i]
      for j in range(n):
        if j==i:
          continue
        latj,lonj,co2j = laty[j],lony[j],co2y[j]
        d,b = getDistanceAndBearingGreatCircle(lati,loni,latj,lonj)
        b = toRadians(b)
        ss = (co2i-co2j)**2
        x.append(d*sin(b))
        y.append(d*cos(b))
        s.append(ss)
    if random:
      title = "Random permutation"
    else:
      title = monthYearString(sy.getValue(iy))
    png = "co2V"
    if random:
      png += "r"
    plotVariogram(s,x,y,title=title,png=png)
  
def goCrossValidate(y):
  eimage = readWaterImage()
  names,lats,lons,co2s = readCo2Table()
  sy = yearsSamplingCo2Table()
  iy = sy.indexOfNearest(y)
  laty,lony,co2y = getGoodData(iy,lats,lons,co2s)
  for scale in [8.0,16.0]:
    pw = PrintWriter(FileOutputStream("Co2e"+str(int(scale))+".txt"))
    f,u1,u2 = co2y,lony,laty
    n = len(f)
    fm,u1m,u2m = zerofloat(n-1),zerofloat(n-1),zerofloat(n-1)
    for i in range(n):
      fi,u1i,u2i = f[i],u1[i],u2[i]
      copy(i, f, fm); copy(n-i-1,i+1, f,i, fm)
      copy(i,u1,u1m); copy(n-i-1,i+1,u1,i,u1m)
      copy(i,u2,u2m); copy(n-i-1,i+1,u2,i,u2m)
      g,s1,s2 = gridBlended(scale,fm,u1m,u2m)
      i1 = s1.indexOfNearest(u1i)
      i2 = s2.indexOfNearest(u2i)
      e = g[i2][i1]-fi
      print "i =",i,"u1 =",u1i,"u2 =",u2i,"e =",e
      pw.printf("%10.4f %10.4f %10.2f %10.2f%n",(u1i,u2i,fi,e))
      pw.flush()
      #plot2(eimage,f,u1,u2,g,s1,s2,None,mv=True,cv=True,tv=False,year=y)
    pw.close()

def printErrorStats():
  print "     emin   eq25   eq50   eq75   emax   eavg   erms   eaad   emad"
  for name in ["1","2","4","8","16"]:
    fileName = "Co2e"+name
    e,f,x1,x2 = readErrors(fileName)
    n = len(f)
    kmin = 0
    kq25 = int(0.25*n+0.5)
    kq50 = int(0.50*n+0.5)
    kq75 = int(0.75*n+0.5)
    kmax = n-1
    quickSort(e)
    emin = e[kmin]
    eq25 = e[kq25]
    eq50 = e[kq50]
    eq75 = e[kq75]
    emax = e[kmax]
    eavg = sum(e)/n
    erms = sqrt(sum(mul(e,e))/n)
    e = abs(e)
    quickSort(e)
    eaad = sum(e)/n
    emad = e[kq50]
    fmt = name+": %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f"
    print fmt % (emin,eq25,eq50,eq75,emax,eavg,erms,eaad,emad)

def goAnalyzeErrors():
  printErrorStats()
  e1,f,x1,x2 = readErrors("Co2e1")
  e2,f,x1,x2 = readErrors("Co2e4")
  n = len(f)
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,240)
  sp.setHLabel("CO2 error (ppm)")
  sp.setVLabel("Frequency")
  sp.setHInterval(2)
  #sp.setVLimits(0,0.35)
  #sp.setHLimits(-450,450)
  sp.addGrid("HV--")
  styles = [PointsView.Line.SOLID,PointsView.Line.DOT]
  colors = [Color.BLACK,Color.RED]
  widths = [5,5]
  errors = [e2,e1]
  for i in range(len(errors)):
    e = errors[i]
    h = Histogram(e,-10.0,10.0,21)
    #h = Histogram(e)
    sh = h.getBinSampling()
    #print "h # bins =",sh.count
    sv = sp.addPoints(sh,h.getDensities())
    sv.setLineColor(colors[i])
    sv.setLineWidth(widths[i])
    sv.setLineStyle(styles[i])
  if pngDir:
    sp.paintToPng(720,3.3,pngDir+"Co2Error14.png")

def goGridding(fy,ly,ky,sg=False,tv=False,p3=False):
  if sg:
    eimage = readEarthImage()
  else:
    eimage = readWaterImage()
  names,lats,lons,co2s = readCo2Table()
  sy = yearsSamplingCo2Table()
  for iy in yearIndicesCo2Table(fy,ly,ky):
    y = sy.getValue(iy)
    laty,lony,co2y = getGoodData(iy,lats,lons,co2s)
    if not laty:
      continue
    #for scale in [1.0,4.0]:
    for scale in [4.0]:
      f,u1,u2 = co2y,lony,laty
      if sg:
        g,s1,s2 = gridSimple(f,u1,u2)
      else:
        g,s1,s2 = gridBlended(scale,f,u1,u2)
      d = makeTensors(scale,s1,s2)
      #title = monthYearString(y)
      title = None
      pngName = "co2Us"+str(int(scale))+"y"+yearIndexString(iy)
      print pngName," # values =",len(f)
      mv = not sg
      plot2(eimage,f,u1,u2,g,s1,s2,d,
            mv=mv,cv=True,tv=tv,
            year=y,title=title,png=pngName)
      if p3:
        cmin,cmax = clipsForYear(y)
        cmin += (cmax-cmin)*1.1/256 # make null (zero) values transparent
        if sg:
          cmod = makeTransparentColorModel(1.0)
        else:
          cmod = makeTransparentColorModel(0.7)
        cmap = ColorMap(cmin,cmax,cmod)
        cimage = SampledImage.fromFloats(s1,s2,g,cmap)
        #cimage.setAlpha(0.8)
        #plot3(eimage,cimage,az= 45)
        #plot3(eimage,cimage,az=-120)
        #plot3(eimage,cimage,az=  60)
        plot3(eimage,cimage,az=-100)
        #plot3(eimage,cimage,az=  80)
        """
        if sg:
          plot3(eimage,cimage,az=-100)
          plot3(eimage,cimage,az=  80)
        else:
          plot3(eimage,cimage,lats=u2,lons=u1,az=-100)
          plot3(eimage,cimage,lats=u2,lons=u1,az=  80)
        """
      if sg:
        break

def goPlotEarth2():
  image = readEarthImage()
  plot2(image,png="earthU")

def goPlotEarth3():
  eimage = readEarthImage()
  plot3(eimage,az=-90)
  plot3(eimage,az= 45)

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
    cv.setClips(cmin,cmax)
    cmin += 1.1*(cmax-cmin)/256 # make null (zero) values transparent
    if mv:
      cv.setColorModel(makeTransparentColorModel(0.7))
    else:
      cv.setColorModel(makeTransparentColorModel(1.0))
    cv.setColorModel(ColorMap.getJet(0.8))
    tile.addTiledView(cv)
    cb = pp.addColorBar("CO2 (ppm)")
    cb.setInterval(5)
  gv = gridLines2(Sampling(13,30,-180),Sampling(7,30,-90),Color.BLACK,gwidth)
  tile.addTiledView(gv)
  if tv and d and s1 and s2:
    gwidth = 1
    tv = TensorsView(s1,s2,d)
    tv.setScale(22.0)
    tv.setLineColor(Color.BLACK)
    tv.setLineWidth(2)
    tv.setEllipsesDisplayed(Sampling(11,30,-150),Sampling(5,30,-60))
    tile.addTiledView(tv)
  if mv and f and u1 and u2:
    gwidth = 1
    mv = PointsView(u1,u2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkColor(Color.WHITE)
    mv.setMarkSize(8)
    tile.addTiledView(mv)
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
    pf.paintToPng(720,3.3,pngDir+png+".png")
  #runFunc(disposePlotFrame,pf)

def disposePlotFrame(pf):
  print "closing ",pf[0]
  pf[0].setVisible(False)
  pf[0].dispose()

def plot3(eimage,cimage=None,lats=None,lons=None,az=45):
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
  frame.setSize(1100,1100)
  world = frame.getWorld()
  world.addChild(group)
  bs = world.getBoundingSphere(True)
  if lats:
    xyz = xyzFromLatsLons(model,lats,lons)
    xyz = mul(1.002,xyz)
    rgb = fillfloat(1.0,3*len(xyz))
    group = PointGroup(80,xyz,rgb)
    world.addChild(group)
  view = frame.getOrbitView()
  view.setProjection(OrbitView.Projection.ORTHOGRAPHIC)
  view.setAzimuth(az)
  if az==-120 or az==60:
    view.setElevation(0)
  view.setWorldSphere(bs)

def xyzFromLatsLons(model,lats,lons):
  n = len(lats)
  xyz = zerofloat(3*n)
  j = 0
  for i in range(n):
    lati = lats[i]
    loni = lons[i]
    xyzi = model.xyz(lati,loni,0.0)
    xyz[j] = xyzi[0]; j += 1
    xyz[j] = xyzi[1]; j += 1
    xyz[j] = xyzi[2]; j += 1
  return xyz

def plotVariogram(s,x,y,title=None,png=None):
  #sp = SimplePlot()
  #pv = sp.addPoints(x,y)
  #pv.setLineStyle(PointsView.Line.NONE)
  #pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  #pv.setMarkSize(2)
  nx = ny = 37 # for 10-degree increments
  #nx = ny = 73 # for 5-degree increments
  sx = sy = Sampling(nx,360/(nx-1),-180)
  #gridder = BlendedGridder2(s,x,y)
  gridder = SimpleGridder2(s,x,y)
  sg = gridder.grid(sx,sy)
  sg = sqrt(sg)
  #sg = clip(max(sg)/255,max(sg),sg)
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,240)
  if title:
    sp.setTitle(title)
    sp.setSize(800,718)
  else:
    sp.setSize(800,660)
  sp.addColorBar("Standard deviation (ppm)")
  #sp.setHLimits(-180.0,180.0)
  #sp.setVLimits(-180.0,180.0)
  sp.setHLimits(-125.0,125.0)
  sp.setVLimits(-125.0,125.0)
  #sp.setHLimits(-90.0,90.0)
  #sp.setVLimits(-90.0,90.0)
  sp.setHLabel("Geodesic east-west distance (degrees)")
  sp.setVLabel("Geodesic north-south distance (degrees)")
  pv = sp.addPixels(sx,sy,sg)
  #pv.setClips(19/255.0,19)
  pv.setClips(13/255.0,13)
  pv.setColorModel(makeTransparentColorModel(1.0))
  #pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  x,y = makeEllipse(120,30)
  pv = sp.addPoints(x,y)
  pv.setLineColor(Color.RED)
  pv.setLineWidth(5)
  if png and pngDir: 
    sp.paintToPng(720,3.3,pngDir+png+".png")

def makeEllipse(a,b):
  nt = 101
  dt = 2.0*PI/(nt-1)
  ft = 0.0
  x,y = [],[]
  for it in range(nt):
    t = ft+it*dt
    x.append(a*cos(t))
    y.append(b*sin(t))
  return x,y

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

def yearIndicesCo2Table(fy,ly,ky):
  sy = yearsSamplingCo2Table()
  ny = sy.count
  iyf = max(sy.indexOfNearest(fy)-1,0)
  iyl = min(sy.indexOfNearest(ly)+1,ny-1)
  while sy.getValue(iyf)<fy:
    iyf += 1
  while sy.getValue(iyl)>ly:
    iyl -= 1
  return range(iyf,iyl+1,ky)

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

def readErrors(fileName):
  s = Scanner(FileInputStream(fileName+".txt"))
  e,f,x1,x2 = [],[],[],[]
  while s.hasNextLine():
    line = s.nextLine()
    data = line.split()
    x1.append(float(data[0]))
    x2.append(float(data[1]))
    f.append(float(data[2]))
    e.append(float(data[3]))
  s.close()
  e = floatsFromList(e)
  f = floatsFromList(f)
  x1 = floatsFromList(x1)
  x2 = floatsFromList(x2)
  return e,f,x1,x2
def floatsFromList(a):
  b = zerofloat(len(a))
  copy(len(a),a,b)
  return b

#############################################################################
# Other stuff

def monthYearString(y):
  mmap = {1:"Jan", 2:"Feb", 3:"Mar", 4:"Apr", 5:"May", 6:"Jun",
          7:"Jul", 8:"Aug", 9:"Sep",10:"Oct",11:"Nov",12:"Dec"}
  year = int(y)
  month = mmap[int(1.0+(y-year)*12.0)]
  return month+", "+str(year)

def yearIndexString(iy):
  if iy<10:
    return "0"+str(iy)
  else:
    return str(iy)

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

def permute(x):
  n = len(x)
  y = zerofloat(n)
  r = Random(314159)
  for i in range(n):
    j = i+r.nextInt(n-i)
    y[i] = x[j]
  return y

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
  n1,d1,f1 = s1.count,s1.delta,s1.first
  n1p,f1p = 2*n1,f1-180.0
  s1p = Sampling(n1p,d1,f1p)
  n1p,n2 = s1p.count,s2.count
  d = makeTensors(scale,s1p,s2)
  bg = BlendedGridder2(d,fp,x1p,x2p)
  bg.setBlending(False)
  gp = bg.grid(s1p,s2)
  g = copy(n1,n2,n1/2,0,gp)
  return g,s1,s2

def padLongitude(f,x1,x2):
  f,x1,x2 = f[:],x1[:],x2[:]
  n = len(f)
  for i in range(n):
    if x1[i]<0.0:
      x1.append(x1[i]+360.0)
    else:
      x1.append(x1[i]-360.0)
    f.append(f[i])
    x2.append(x2[i])
  f = floatsFromList(f)
  x1 = floatsFromList(x1)
  x2 = floatsFromList(x2)
  return f,x1,x2

def wrapLon(f,sx,sy):
  nx,ny = sx.count,sy.count
  g = zerofloat(nx+1,ny)
  copy(nx,ny,f,g)
  copy(1,ny,0,0,f,nx,0,g)
  sx = Sampling(nx+1,sx.delta,sx.first)
  return g,sx,sy

"""
def testDB():
  #getDistanceAndBearingGreatCircle(0,-45,0,45)
  #getDistanceAndBearingGreatCirclM(0,-45,0,45)
  #getDistanceAndBearingGreatCircle(30,-45,30,45)
  #getDistanceAndBearingGreatCirclM(30,-45,30,45)
  getDistanceAndBearingGreatCircle(30,-45,60,45)
  getDistanceAndBearingGreatCirclM(30,-45,60,45)
def main(args):
  testDB()
"""

def getDistanceAndBearingGreatCircle(lat1,lon1,lat2,lon2):
  d12,b12 = getDistanceAndBearingGC(lat1,lon1,lat2,lon2)
  d21,b21 = getDistanceAndBearingGC(lat2,lon2,lat1,lon1)
  if b21<0.0:
    b21 += 180.0
  else:
    b21 -= 180.0
  d = 0.5*(d12+d21)
  b = 0.5*(b12+b21)
  #print "b12 =",b12," b21 =",b21," b =",b
  #print "d =",d," b =",b
  return d,b

def getDistanceAndBearingGC(lat1,lon1,lat2,lon2):
  lat1,lon1 = toRadians(lat1),toRadians(lon1)
  lat2,lon2 = toRadians(lat2),toRadians(lon2)
  dlon = lon2-lon1
  d = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))
  y = sin(dlon)*cos(lat2)
  x = cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon)
  b = atan2(y,x)
  return toDegrees(d),toDegrees(b)

def getDistanceAndBearingGreatCirclM(lat1,lon1,lat2,lon2):
  latm,lonm = getMidpoint(lat1,lon1,lat2,lon2)
  d,b = getDistanceAndBearingGC(latm,lonm,lat2,lon2)
  d *= 2.0
  #print "d =",d," b =",b
  return d,b

def getMidpoint(lat1,lon1,lat2,lon2):
  lat1,lon1 = toRadians(lat1),toRadians(lon1)
  lat2,lon2 = toRadians(lat2),toRadians(lon2)
  dlon = lon2-lon1
  bx = cos(lat2)*cos(dlon)
  by = cos(lat2)*sin(dlon)
  latm = atan2(sin(lat1)+sin(lat2),sqrt((cos(lat1)+bx)*(cos(lat1)+bx)+by*by)) 
  lonm = lon1+atan2(by,cos(lat1)+bx)
  return toDegrees(latm),toDegrees(lonm)

def getDistanceAndBearingRhumbLine(lat1,lon1,lat2,lon2):
  lat1,lon1 = toRadians(lat1),toRadians(lon1)
  lat2,lon2 = toRadians(lat2),toRadians(lon2)
  dlat,dlon = lat2-lat1,lon2-lon1
  dphi = log(tan(0.5*(lat2+0.5*PI))/tan(0.5*(lat1+0.5*PI)))
  if dphi!=0.0:
    q = dlat/dphi
  else:
    q = cos(lat1)
  if dlon>PI:
    dlon -= 2.0*PI
  elif dlon<-PI:
    dlon += 2.0*PI
  d = sqrt(dlat*dlat+q*q*dlon*dlon) # rhumb-line distance
  b = toDegrees(atan2(dlon,dphi)) # so bearing is constant
  return d,b

def printf(fmt,*varargs):
  sys.stdout.write(fmt % varargs)

#############################################################################
if __name__ == "__main__":
  run(main)
