from shared import *

def main(args):
  goInterp()

def goInterp():
  s1,s2 = getTaiwanSamplings()
  g = readImage("TaiwanSlope",s1,s2)
  f,x1,x2 = readVs30s("TaiwanVs30s",s1,s2)
  mask = makeMask(g)
  h = log10(add(0.0001,g)) # log(Vs30) ~linear function of log(slope)
  d = makeTensors(h,mask)
  dmask = makeTensors(h,mask)
  maskTensors(dmask,mask)
  plot(f,x1,x2,h,s1,s2,cbar="Log10[slope (m/m)]",png="vs30s")
  plot(None,x1,x2,h,s1,s2,dmask,cbar="Log10[slope (m/m)]",png="vs30se")
  f = log(f)
  for tensors in [d,None]:
    if tensors:
      gridder = BlendedGridder2(tensors,f,x1,x2)
      namet = "vs30t"
    else:
      gridder = BlendedGridder2(f,x1,x2)
      namet = "vs30i"
    for blending in [False,True]:
      if blending: 
        name = namet+"q"
      else: 
        name = namet+"p"
      gridder.setBlending(blending)
      #gridder.setSmoothness(1.0)
      q = gridder.grid(s1,s2)
      q = exp(q)
      q = maskImage(q,mask)
      plot(f,x1,x2,q,s1,s2,cmin=0.0,cmax=760.0,cbar="Vs30 (m/s)",png=name)

#############################################################################
# plotting

#pngDir = "png/vs30/" # where to put PNG images of plots
pngDir = None # for no PNG images

def plot(f,x1,x2,g,s1,s2,d=None,cbar=None,cmin=0,cmax=0,png=None):
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
    tv.setScale(300)
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

def makeTensors(g,mask):
  lof = LocalOrientFilter(4.0)
  d = lof.applyForTensors(g)
  d.invertStructure(1.0,2.0)
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av) # eigenvalues are gradients squared; au >= av
  au = mul(au,mask) # au on land, zero on water
  av = mul(av,mask) # av on land, zero on water
  am = sub(1.0,mask) # zero on land, one on water
  au = add(au,am) # merge tensors for land and water
  av = add(av,am) # merge tensors for land and water
  d.setEigenvalues(au,av)
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

#############################################################################
run(main)
