#############################################################################
# Demo interpolation on earth's surface using data from the World Ocean Atlas
"""
Get mask from objectively analyzed data
Construct tensors
Get sample means and standard deviations
Make random subset of sample means (with new mask)
Subtract annual mean from sample means
Fit using tensor-guided kriging
"""

from shared import *
import gph.Woa as woa
from interp import *

pngDir = None
#pngDir = "../../png/woa/"
slon = woa.getLonSampling()
slat = woa.getLatSampling()

def main(args):
  #goShowData()
  goKriging()

def goKriging():
  namet = "s00an1" # trend
  namea = "s10an1" # NODC gridded
  namem = "s10mn1" # sample mean
  named = "s10sd1" # sample std dev
  namee = "s10se1" # sample std err
  namen = "s10dd1" # number of samples
  namek = "s10kg1" # Kriging gridded
  depth = 100
  level = depthLevel(depth)
  ft = readDataFile(namet,level)
  fa = readDataFile(namea,level)
  fm = readDataFile(namem,level)
  fd = readDataFile(named,level)
  fe = readDataFile(namee,level)
  fn = readDataFile(namen,level)
  fd = woa.sdfix(fn,fd)
  fe = woa.sefix(fn,fe)
  print "median of sd =",woa.nmed(1,10,fn,fd)
  print "median of se =",woa.nmed(1,10,fn,fe)
  livet = woa.mask(ft,fieldNull(namet))
  livea = woa.mask(fa,fieldNull(namea))
  livem = woa.mask(fm,fieldNull(namem))
  lives = woa.subset(Random(1),10,livem)
  #plotField3(namet,depth,ft,live=livet)
  plotField3(namen,depth,fn,live=livem)
  plotField3(namem,depth,fm,live=livem)
  plotField3(namet,depth,ft,live=livet)
  plotField3(namem,depth,sub(fm,ft),live=livem)
  plotField3(namea,depth,sub(fa,ft),live=livea)
  #plotField3(named,depth,fd,live=livem)
  #plotField3(namee,depth,fe,live=livem)
  #plotField3(named,depth,fd,live=lives)
  plotField3(namem,depth,sub(fm,ft),live=lives)
  return
  fr = sub(fm,ft)
  fs,lons,lats = woa.scattered(lives,fr)
  ds,lons,lats = woa.scattered(lives,fd)
  #ds = add(ds,???)
  tensors = woa.makeTensors(livea)
  sc = SmoothCovariance(1.0,1.0,5.0,2)
  nlon,nlat = slon.count,slat.count
  sc.testSpd(nlon,nlat,tensors)
  print "sc"
  kg = KrigingGridder2(tensors,fs,lons,lats)
  kg.setModelCovariance(sc)
  #kg.setPaciorek(True)
  kg.setDataError(ds)
  #kg = BlendedGridder2(tensors,fs,lons,lats)
  print "kg"
  gr = kg.grid(slon,slat)
  plotField3(namek,depth,gr,live=livea)
  gm = add(gr,ft)
  print "gr"
  #plotTrend(trend)
  #plotField3(namem,depth,fr,live=lives)
  #plotField3(namek,depth,gr,live=livea)
  #plotField3(namem,depth,fm,live=lives)
  #plotField3(named,depth,dm,live=lives)
  #plotField3(namek,depth,gm,live=livea)
  #plotField3(namea,depth,fa,live=livea)

def goShowData():
  depth = 100
  #plotFile3("d00dl1",0)
  #plotFile3("b00bn1",depth)
  #plotFile3("s10an1",depth)
  #plotFile3("s10dd1",depth)
  #plotFile3("s10sd1",depth)
  #plotFile3("s10mn1",depth)
  #plotFile3("t00an1",depth)
  plotFile3("t10dd1",depth)
  plotFile3("t10sd1",depth)
  plotFile3("t10se1",depth)
  #plotFile3("t00mn1",depth)
  #plotFile3("t10oa1",depth)
  #plotFile3()

def plotTrend(trend):
  nlat = slat.count
  tlat = zerofloat(nlat)
  for ilat in range(nlat):
    tlat[ilat] = trend.interpolate(slat.getValue(ilat))
  SimplePlot.asPoints(slat,tlat)
def findTrend(m,f):
  nlon,nlat = slon.count,slat.count
  favg = zerofloat(nlat)
  lavg = zerofloat(nlat)
  navg = zeroint(nlat)
  n = 0
  for ilat in range(nlat):
    fsum = 0.0
    nsum = 0
    for ilon in range(nlon):
      if m[ilat][ilon]:
        fsum += f[ilat][ilon]
        nsum += 1
    if nsum>0:
      favg[ilat] = fsum/nsum
      navg[ilat] = nsum
      n += 1
  klat = 0
  for ilat in range(nlat):
    if navg[ilat]>0:
      favg[klat] = favg[ilat]
      lavg[klat] = slat.getValue(ilat)
      klat += 1
  favg = copy(n,favg)
  lavg = copy(n,lavg)
  trend = CubicInterpolator(CubicInterpolator.Method.LINEAR,n,lavg,favg)
  return trend
def subTrend(trend,m,f):
  nlon,nlat = slon.count,slat.count
  g = zerofloat(nlon,nlat)
  for ilat in range(nlat):
    tlat = trend.interpolate(slat.getValue(ilat))
    for ilon in range(nlon):
      if m[ilat][ilon]:
        g[ilat][ilon] = f[ilat][ilon]-tlat
  return g
def addTrend(trend,m,f):
  nlon,nlat = slon.count,slat.count
  g = zerofloat(nlon,nlat)
  for ilat in range(nlat):
    tlat = trend.interpolate(slat.getValue(ilat))
    for ilon in range(nlon):
      if m[ilat][ilon]:
        g[ilat][ilon] = f[ilat][ilon]+tlat
  return g
 
def plotFile3(name=None,depth=0):
  if name:
    f = readDataFile(name,depthLevel(depth))
  plotField3(name,depth,f)

def plotField3(name,depth,f,live=None):
  eimage = readEarthImage()
  #eimage = makeGrayImage(eimage)
  fimage,cmap,cbar,title = None,None,None,None
  if name and f:
    if not live:
      fnull = fieldNull(name)
      live = woa.mask(f,fnull)
    fmin = woa.min(live,f)
    fmax = woa.max(live,f)
    print name,depth,": live =",woa.count(live)," fmin =",fmin," fmax =",fmax
    cmin,cmax = fieldClips(name,live,f)
    cmap = colorMap(name,cmin,cmax)
    fimage = SampledImage.fromFloats(slon,slat,live,f,cmap)
    cbar = valueLabel(name)
    title = makeTitle(name,depth)
  plot3(eimage,fimage,cmap=cmap,cbar=cbar,title=title)

#############################################################################
# WOA data utilities

valueNamesMap = {
  "b":"basin number",
  "d":"depth levels",
  "n":"nitrate",
  "s":"salinity",
  "t":"temperature"
}
valueUnitsMap = {
  "n":"micromole/l",
  "s":"PSU",
  "t":"deg C"
}
fieldNamesMap = {
  "bn":"basin number",
  "dl":"depth level",
  "an":"NODC gridded",
  "mn":"sample mean",
  "dd":"number of samples",
  "sd":"standard deviation",
  "se":"standard error",
  "oa":"mean minus gridded",
  "gp":"number of nearby samples",
  "kg":"kriging gridded"
}
fieldUnitsMap = {
  "bn":"basin number",
  "dl":"depth level",
  "an":None,
  "mn":None,
  "dd":"number of samples",
  "sd":None,
  "se":None,
  "oa":None,
  "gp":"number of nearby values",
  "kg":None,
}
fieldNullsMap = {
  "bn":-100.0,
  "dl":1.0,
  "an":-99.9999,
  "mn":-99.9999,
  "dd":-100.0,
  "sd":-99.9999,
  "se":-99.9999,
  "oa":-99.9999,
  "gp":-100.0,
  "kg":-99.9999
}
monthNamesMap = {
  "00":"Annual",
  "01":"January",
  "02":"February",
  "03":"March",
  "04":"April",
  "05":"May",
  "06":"June",
  "07":"July",
  "08":"August",
  "09":"September",
  "10":"October",
  "11":"November",
  "12":"December"
}
depthLevelsMap = {
  0:1,10:2,20:3,30:4,50:5,75:6,100:7,125:8,150:9,200:10,
  250:11,300:12,400:13,500:14,600:15,700:16,800:17,900:18,
  1000:19,1100:20,1200:21,1300:22,1400:23,1500:24
}
def colorMap(name,fmin,fmax):
  cs = (
    Color.RED,Color.YELLOW,Color.GREEN,
    Color.CYAN,Color.BLUE,Color.MAGENTA,
  )
  nc = len(cs)
  if fieldCode(name)=="bn":
    ran = Random(6)
    n = 1+int(fmax-fmin)
    r = zerofloat(256)
    g = zerofloat(256)
    b = zerofloat(256)
    for i in range(256):
      k = int(i*(n-1.0)/255.0+0.5)
      k = k%nc
      s = 1.0-0.5*i/256
      r[i] = s*cs[k].red/255.0
      g[i] = s*cs[k].green/255.0
      b[i] = s*cs[k].blue/255.0
    return ColorMap(fmin,fmax,r,g,b)
  else:
    return ColorMap(fmin,fmax,ColorMap.JET)
def fieldCode(name):
  return name[3:5]
def fieldNull(name):
  return fieldNullsMap[fieldCode(name)]
def fieldName(name):
  return fieldNamesMap[fieldCode(name)]
def fieldUnits(name):
  return fieldUnitsMap[fieldCode(name)]
def fieldClips(name,m,f):
  code = fieldCode(name)
  if code=="bn":
    return 1,58
  elif code=="dl":
    return 2,40
  else:
    return woa.clips(2,98,m,f)
def valueCode(name):
  return name[0:1]
def valueName(name):
  return valueNamesMap[valueCode(name)]
def valueLabel(name):
  label = fieldUnits(name)
  if label:
    label = label.capitalize()
  else:
    units = valueUnits(name)
    field = fieldName(name)
    if field[0:4]!="NODC":
      field = field.capitalize()
    label = field+" ("+units+")"
  return label
def valueUnits(name):
  return valueUnitsMap[valueCode(name)]
def monthCode(name):
  return name[1:3]
def monthName(name):
  return monthNamesMap[monthCode(name)]
def depthLevel(depth):
  return depthLevelsMap[depth]
def makeTitle(name,depth):
  value = valueName(name)
  month = monthName(name)
  field = fieldName(name)
  if field[0:4]!="NODC":
    field = field.lower()
  title = month+" "+value+" at depth "+str(depth)+" m: "+field
  return title

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

#############################################################################
# Plot

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
  gwidth = 1
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
    #cv.setColorModel(ColorMap.getJet(0.8))
    tile.addTiledView(cv)
    cb = pp.addColorBar("CO2 (ppm)")
    cb.setInterval(5)
  gv = gridLines2(Sampling(13,30,-180),Sampling(7,30,-90),Color.BLACK,gwidth)
  tile.addTiledView(gv)
  if tv and d and s1 and s2:
    gwidth = 1
    tv = TensorsView(s1,s2,d)
    tv.setScale(14.0)
    tv.setLineColor(Color.RED)
    tv.setLineWidth(3)
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
      #pf.setSize(1740,770) # print
      pf.setSize(1740,750) # slide
    else: 
      pf.setSize(1078,562)
  #pf.setFontSizeForPrint(8,240)
  pf.setFontSizeForSlide(1.0,0.9)
  pf.setVisible(True)
  if png and pngDir: 
    pf.paintToPng(720,3.3,pngDir+png+".png")
  #runFunc(disposePlotFrame,pf)

def disposePlotFrame(pf):
  print "closing ",pf[0]
  pf[0].setVisible(False)
  pf[0].dispose()

def plot3(eimage,cimage=None,lats=None,lons=None,
          cmap=None,cbar=None,cint=None,title=None,az=-90):
  model = OblateModel.forWGS84()
  group = OblateImageGroup(model)
  group.addImage(eimage,0.999)
  gwidth = 3
  if cimage:
    group.addImage(cimage,1.000)
    gwidth = 1
  group.addGrid(Color.BLACK,gwidth,1.001)
  frame = SimpleFrame(AxesOrientation.XOUT_YRIGHT_ZUP)
  if title:
    frame.setTitle(title)
  if cmap:
    if cbar:
      cbar = ColorBar(cbar)
    else:
      cbar = ColorBar("")
    cbar.setFont(Font("Arial",Font.PLAIN,48))
    if cint:
      cbar.setInterval(cint)
    cbar.setWidthMinimum(180)
    frame.add(cbar,BorderLayout.EAST)
    cmap.addListener(cbar)
    frame.setSize(1000,750)
  else:
    frame.setSize(790,750)
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
  #view.setProjection(OrbitView.Projection.ORTHOGRAPHIC)
  view.setAzimuth(az)
  if az==-120 or az==60:
    view.setElevation(0)
  view.setWorldSphere(bs)

def addColorBar3d(frame,cmap,clab,cint=None):
  cbar = ColorBar(clab)
  cbar.setFont(Font("Arial",Font.PLAIN,64))
  cbar.setBackground(background)
  if cint:
    cbar.setInterval(cint)
  cbar.setWidthMinimum(220)
  frame.add(cbar,BorderLayout.EAST)
  cmap.addListener(cbar)

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
  #sp.setFontSizeForPrint(8,240)
  sp.setFontSizeForSlide(1.0,0.9)
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

imageDir = "/Users/dhale/Home/git/idh/bench/src/gph/resources/"

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

earthImage = None
def readEarthImage():
  global earthImage
  if not earthImage:
    slon = Sampling(2048,360.0/2047,-180.0)
    slat = Sampling(1024,180.0/1023, -90.0)
    earthImage = SampledImage.fromFile(slon,slat,imageDir+"earth.png")
  return earthImage

waterImage = None
def readWaterImage():
  global waterImage
  if not waterImage:
    slon = Sampling(2048,360.0/2047,-180.0)
    slat = Sampling(1024,180.0/1023, -90.0)
    image = SampledImage.fromFile(slon,slat,imageDir+"water.png")
    sx = image.getSamplingX()
    sy = image.getSamplingY()
    grays = image.getFloatsGray(False)
    grays = clip(0.7,1.0,grays)
    cmap = ColorMap(0.0,1.0,ColorMap.GRAY)
    waterImage = SampledImage.fromFloats(sx,sy,grays,cmap)
  return waterImage

#############################################################################
# WOA data

def readDataFile(fileName,level):
  dataDir = "/data/earth/woa/2009/"
  return Woa.readData(dataDir+fileName,level)

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

def shuffle(x):
  r = Random(314159)
  n = len(x)
  y = x[:]
  y[0] = x[0]
  for i in range(1,n):
    j = r.nextInt(i)
    y[i] = y[j]
    y[j] = x[i]
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
  #n1,n2 = 180,90
  #n1,n2 = 120,60
  n1,n2 = 90,45
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
  bg.setBlending(True)
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
