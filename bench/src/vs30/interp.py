from shared import *

pngDir = None # for no PNG images
#pngDir = "png/vs30/" # where to put PNG images of plots

def main(args):
  #goCrossPlot()
  goAnalyzeErrors()
  #goCrossValidate()
  #goInterp()

def goCrossPlot():
  s1,s2 = getTaiwanSamplings()
  g = readImage("TaiwanSlope",s1,s2)
  f,x1,x2 = readVs30s("TaiwanVs30s",s1,s2)
  g = safeLog10(g) # log(Vs30) ~linear function of log(slope)
  f = safeLog10(f)
  n = len(f)
  h = zerofloat(n)
  for i in range(n):
    i1 = s1.indexOfNearest(x1[i])
    i2 = s2.indexOfNearest(x2[i])
    h[i] = g[i2][i1]
  sp = SimplePlot()
  sp.setSize(840,730)
  sp.setFontSizeForPrint(8,240)
  sp.setHLabel("Log10[slope (m/m)]")
  sp.setVLabel("Log10[Vs30 (m/s)]")
  sp.setHLimits(-3.7,0.2)
  sp.setVLimits( 1.9,3.3)
  #sp.setHLimits(0,0.1)
  #sp.setVLimits(0,810)
  sp.addGrid("HV--")
  h1,f1 = -3.7,AFIT-BFIT*3.7
  h2,f2 =  0.2,AFIT+BFIT*0.2
  pv = sp.addPoints([h1,h2],[f1,f2])
  pv.setLineColor(Color.GRAY)
  pv.setLineWidth(5)
  pv.setLineStyle(PointsView.Line.DASH)
  pv = sp.addPoints(h,f)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(8)
  if pngDir:
    sp.paintToPng(720,3.3,pngDir+"/vs30sVsSlope.png")

def goInterp():
  s1,s2 = getTaiwanSamplings()
  v,x1,x2 = readVs30s("TaiwanVs30s",s1,s2)
  vmin,vmax = min(v),max(v)
  s = readImage("TaiwanSlope",s1,s2) # slope
  mask = makeMask(s)
  logs = safeLog10(s)
  plot(v,x1,x2,logs,s1,s2,cbar="Log10[slope (m/m)]",png="vs30s")
  vmod = getVs30Modeled(logs)
  vmod = clip(vmin,vmax,vmod)
  vmod = maskImage(vmod,mask)
  plot(v,x1,x2,vmod,s1,s2,
       cmin=0,cmax=810,cbar="Vs30 from slope (m/s)",png="vs30vm")
  #tests = [(0,0,1),(0,1,1),(0,2,1),(1,1,300),(1,2,300)]
  tests = [(0,0,1),(0,2,1),(1,1,300)]
  #tests = [(0,0,1)]
  ntest = len(tests)
  for i in range(ntest):
    p0,p1,es = tests[i]
    name = "vs30bg"+str(p0)+str(p1)
    d = makeTensors(logs,mask,p0,p1)
    dm = EigenTensors2(d); maskTensors(dm,mask)
    plot(None,x1,x2,logs,s1,s2,dm,es=es,
      cbar="Log10[slope (m/m)]",png=name+"d")
    r = subtractVs30Modeled(v,x1,x2,logs,s1,s2)
    gridder = BlendedGridder2(d,r,x1,x2)
    for blending in [True]:
      if blending: png = name+"q"
      else:        png = name+"p"
      gridder.setBlending(blending)
      q = gridder.grid(s1,s2)
      #plot(v,x1,x2,q,s1,s2,cmin=-405,cmax=405,
      #  cbar="Vs30 interpolated residual (m/s)",png=png)
      q = addVs30Modeled(q,logs)
      q = clip(vmin,vmax,q)
      q = maskImage(q,mask)
      plot(v,x1,x2,q,s1,s2,cmin=0.0,cmax=810.0,
        cbar="Vs30 interpolated (m/s)",png=png)

def goCrossValidate():
  s1,s2 = getTaiwanSamplings()
  v,x1,x2 = readVs30s("TaiwanVs30s",s1,s2)
  vmin,vmax = min(v),max(v)
  s = readImage("TaiwanSlope",s1,s2) # slope
  mask = makeMask(s)
  logs = safeLog10(s)
  vmod = getVs30Modeled(logs)
  vmod = clip(vmin,vmax,vmod)
  vmod = maskImage(vmod,mask)
  r = subtractVs30Modeled(v,x1,x2,logs,s1,s2)
  n = len(v)
  d = makeTensors(logs,mask,1.0,1.0)
  gridder = BlendedGridder2(d)
  gridder.setBlending(True)
  rm,x1m,x2m = zerofloat(n-1),zerofloat(n-1),zerofloat(n-1)
  pw = PrintWriter(FileOutputStream("TaiwanVs30e.txt"))
  for i in range(n):
    ri,vi,x1i,x2i = r[i],v[i],x1[i],x2[i]
    i1 = s1.indexOfNearest(x1i)
    i2 = s2.indexOfNearest(x2i)
    copy(i, r, rm); copy(n-i-1,i+1, r,i, rm)
    copy(i,x1,x1m); copy(n-i-1,i+1,x1,i,x1m)
    copy(i,x2,x2m); copy(n-i-1,i+1,x2,i,x2m)
    gridder.setScattered(rm,x1m,x2m)
    q = gridder.grid(s1,s2)
    q = addVs30Modeled(q,logs)
    q = clip(vmin,vmax,q)
    q = maskImage(q,mask)
    e = q[i2][i1]-vi
    print "i =",i,"x1 =",x1i,"x2 =",x2i,"e =",e
    #plot(rm,x1m,x2m,q,s1,s2,cmin=0.0,cmax=810.0,cbar="Vs30 (m/s)")
    pw.printf("%10.4f %10.4f %10.2f %10.2f%n",(x1i,x2i,vi,e))
    pw.flush()
  pw.close()

def printErrorStats():
  print "        emin   eq25   eq50   eq75   emax   eavg   erms   eaad   emad"
  for name in ["00d","01d","02d","11d","12d","00e","02e","11e"]:
    fileName = "TaiwanVs30ErrorBg"+name
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
  e1,f,x1,x2 = readErrors("TaiwanVs30ErrorBg11d")
  e2,f,x1,x2 = readErrors("TaiwanVs30ErrorBg00d")
  n = len(f)
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,240)
  sp.setHLabel("Vs30 (m/s)")
  sp.setVLabel("Frequency")
  #sp.setVLimits(0,0.35)
  sp.setHLimits(-450,450)
  sp.addGrid("HV--")
  styles = [PointsView.Line.SOLID,PointsView.Line.DASH]
  colors = [Color.BLACK,Color.RED]
  widths = [5,5]
  errors = [e1,e2]
  for i in range(len(errors)):
    e = errors[i]
    #h = Histogram(e,-500.0,500.0,41)
    h = Histogram(e,-1000.0,1000.0,81)
    sh = h.getBinSampling()
    #print "h # bins =",sh.count
    sv = sp.addPoints(sh,h.getDensities())
    sv.setLineColor(colors[i])
    sv.setLineWidth(widths[i])
    sv.setLineStyle(styles[i])
    #gridder = NearestGridder2(e,x1,x2)
    #q = gridder.grid(s1,s2)
    #q = maskImage(q,mask)
    #plot(e,x1,x2,q,s1,s2,cmin=-500.0,cmax=500.0,cbar="delta Vs30 (m/s)")
  if pngDir:
    sp.paintToPng(720,3.3,pngDir+"vs30sErrorsBg.png")

#############################################################################
# plotting

def plot(f,x1,x2,g,s1,s2,d=None,es=1,cbar=None,cmin=0,cmax=0,png=None):
  sp = SimplePlot()
  sp.setHLabel("Longitude (degrees)")
  sp.setVLabel("Latitude (degrees)")
  sp.setHInterval(1.0)
  sp.setVInterval(1.0)
  sp.setFontSizeForPrint(8,240)
  sp.setSize(865,930)
  if cbar:
    sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(120)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setLineColor(Color.BLACK)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(30)
    tv.setScale(es)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if f and x1 and x2:
    mv = sp.addPoints(x1,x2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkColor(Color.WHITE)
    mv.setMarkSize(6)
  if pngDir and png:
    sp.paintToPng(600,3,pngDir+png+".png")

#############################################################################
# processing

"""
Logs: x = 10^log10(x), ln(x) = ln(10)*log10(x), log10(x) = ln(x)/ln(10)
Eyeball fit of line in crossplot of log10(vs30) versus log10(slope):
log10(v) ~ 2.10 + (0.7/3.0)*(log10(s)+3.5) 
         ~ 2.92 + 0.23*log10(s)
"""
def exp10(x):
  return exp(mul(log(10.0),x))
def safeLog10(x):
  return log10(add(0.0001,x))
def safeExp10(x):
  return sub(exp(mul(log(10.0),x)),0.0001)
AFIT,BFIT = 3.00,0.26
def getVs30Modeled(logs):
  return safeExp10(add(AFIT,mul(BFIT,logs)))
def addVs30Modeled(q,logs):
  return add(q,getVs30Modeled(logs))
def subtractVs30Modeled(v,x1,x2,logs,s1,s2):
  n = len(v)
  logv = zerofloat(n)
  for i in range(n):
    i1 = s1.indexOfNearest(x1[i])
    i2 = s2.indexOfNearest(x2[i])
    logv[i] = AFIT+BFIT*logs[i2][i1]
  return sub(v,safeExp10(logv))
def addLogVs30Modeled(logv,logs):
  return add(logv,add(AFIT,mul(BFIT,logs)))
def subtractLogVs30Modeled(logv,x1,x2,logs,s1,s2):
  n = len(logv)
  r = zerofloat(n)
  for i in range(n):
    i1 = s1.indexOfNearest(x1[i])
    i2 = s2.indexOfNearest(x2[i])
    r[i] = logv[i]-(AFIT+BFIT*logs[i2][i1])
  return r

def extrapSlopes(g,mask):
  n1,n2 = len(g[0]),len(g)
  s1,s2 = Sampling(n1),Sampling(n2)
  f,x1,x2 = [],[],[] # lists of samples on coast
  for i2 in range(1,n2-1):
    for i1 in range(1,n1-1):
      if mask[i2][i1]!=0:
        if (mask[i2-1][i1-1]==0 or
            mask[i2-1][i1  ]==0 or
            mask[i2-1][i1+1]==0 or
            mask[i2  ][i1-1]==0 or
            mask[i2  ][i1+1]==0 or
            mask[i2+1][i1-1]==0 or
            mask[i2+1][i1  ]==0 or
            mask[i2+1][i1+1]==0):
          f.append(g[i2][i1])
          x1.append(float(i1))
          x2.append(float(i2))
  bg = BlendedGridder2(f,x1,x2)
  h = bg.grid(s1,s2)
  g = add(mul(g,mask),mul(h,sub(1.0,mask)))
  return g

def makeTensors(g,mask,p0,p1):
  #g = extrapSlopes(g,mask)
  lof = LocalOrientFilter(4.0)
  d = lof.applyForTensors(g)
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av)
  aw = sub(1.0,mask) # one on water, zero on land
  ab = mul(max(au),aw) # big on water, zero on land
  au = add(mul(au,mask),ab)
  av = add(mul(av,mask),ab)
  d.setEigenvalues(au,av)
  d.invertStructure(p0,p1)
  return d

def makeTensorsGeodetic(e,s1,s2,mask):
  lof = LocalOrientFilter(0.0)
  d = lof.applyForTensors(e)
  n1,n2 = s1.count,s2.count
  d1,d2 = s1.delta,s2.delta
  # m/m = m/sample * sample/deg * deg/m
  #     = m/sample * sample/deg * deg/rad * rad/m
  #     = m/sample * 1.0/(deg/sample*rad/deg*m/rad)
  scale = pow(d1*PI/180.0*6378137.0,-2.0) # for gradient-squared
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av)
  au = add(1.0,mul(scale,au))
  av = add(1.0,mul(scale,av))
  aw = sub(1.0,mask) # one on water, zero on land
  ab = mul(max(au),aw) # big on water, zero on land
  au = add(mul(au,mask),ab)
  av = add(mul(av,mask),ab)
  d.setEigenvalues(au,av)
  d.invert()
  return d

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

def maskImage(z,mask):
  """Applies the specified mask, so that values for water are 0.0."""
  return mul(z,mask)

def maskTensors(d,mask):
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av)
  au = mul(mask,au)
  av = mul(mask,av)
  d.setEigenvalues(au,av)

def median(a):
  n = len(a)
  mf = MedianFinder(n)
  return mf.findMedian(a)

def smoothBoundary(z):
  n1,n2 = len(z[0]),len(z)
  y = copy(z)
  zmin = min(z)
  for i2 in range(1,n2-1):
    for i1 in range(1,n1-1):
      zmm = z[i2-1][i1-1]; zm0 = z[i2-1][i1  ]; zmp = z[i2-1][i1+1]
      z0m = z[i2  ][i1-1]; z00 = z[i2  ][i1  ]; z0p = z[i2  ][i1+1]
      zpm = z[i2+1][i1-1]; zp0 = z[i2+1][i1  ]; zpp = z[i2+1][i1+1]
      if z00==zmin:
        s = 0.0; c = 0.0
        if zmm>zmin:
          s += zmm; c += 1.0
        if zm0>zmin:
          s += zm0; c += 1.0
        if zmp>zmin:
          s += zmp; c += 1.0
        if z0m>zmin:
          s += z0m; c += 1.0
        if z0p>zmin:
          s += z0p; c += 1.0
        if zpm>zmin:
          s += zpm; c += 1.0
        if zp0>zmin:
          s += zp0; c += 1.0
        if zpp>zmin:
          s += zpp; c += 1.0
        if c>0.0:
          y[i2][i1] = s/c
  return y

#############################################################################
# data input/output

dataDir = "/data/earth/vs30/"

def getTaiwanSamplings():
  """
  Returns samplings of longitude and latitude for the Taiwan subset.
  This subset is consistent with the map of Vs30 on the USGS website.
  """
  nlon,nlat = 312,418
  dlon,dlat = 0.5/60,0.5/60 # 0.0083333333333 degrees = 30 seconds
  flon,flat = 119.60416666,21.84583333
  return Sampling(nlon,dlon,flon),Sampling(nlat,dlat,flat)

def readVs30s(fileName,s1,s2):
  """Reads the file and returns lists of measured Vs30, lon, and lat."""
  x1min,x1max = s1.first,s1.last
  x2min,x2max = s2.first,s2.last
  s = Scanner(FileInputStream(dataDir+fileName+".txt"))
  f,x1,x2 = [],[],[]
  while s.hasNextLine():
    line = s.nextLine()
    data = line.split()
    x1i = float(data[0])
    x2i = float(data[1])
    fi = float(data[2])
    if x1min<=x1i and x1i<=x1max and x2min<=x2i and x2i<=x2max:
      x1.append(x1i)
      x2.append(x2i)
      f.append(fi)
  s.close()
  f = floatsFromList(f)
  x1 = floatsFromList(x1)
  x2 = floatsFromList(x2)
  print "read",len(f),"scattered values"
  sg = SimpleGridder2(f,x1,x2) # we have multiple values for some locations
  g = sg.grid(s1,s2) # so we bin the data, averaging values within each bin
  f,x1,x2 = sg.getGriddedSamples(0.0,s1,s2,g)
  print "have",len(f),"values after binning"
  return f,x1,x2
def floatsFromList(a):
  b = zerofloat(len(a))
  copy(a,b)
  return b

def readImage(fileName,s1,s2):
  """Reads and returns an image with specified samplings from a file."""
  n1,n2 = s1.count,s2.count
  z = zerofloat(n1,n2)
  ais = ArrayInputStream(dataDir+fileName+".dat")
  ais.readFloats(z)
  ais.close()
  return z

def writeImage(fileName,z):
  """Writes the specified image to a file."""
  aos = ArrayOutputStream(dataDir+fileName+".dat")
  aos.writeFloats(z)
  aos.close()

def readErrors(fileName):
  s = Scanner(FileInputStream(dataDir+"errors/"+fileName+".txt"))
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

#############################################################################
run(main)
