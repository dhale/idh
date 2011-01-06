from shared import *

pngDir = None # for no PNG files
#pngDir = "png/vs30/"

def main(args):
  goGeology()
  #goSlopes()

def goGeology():
  n1,n2 = 351,456
  s1 = Sampling(n1,0.01,118.93-360.0)
  s2 = Sampling(n2,0.01,21.24)
  f,x1,x2 = readScattered("vs30MeasuredTaiwan.txt",s1,s2)
  g = readImage("vs30GeologyTaiwan.dat",n1,n2)
  cleanGeologyImage(g)
  mask = makeMask(g)
  d = makeTensors(mask,g)
  dmask = makeTensors(mask,g)
  maskTensors(mask,dmask)
  g = mul(mask,g)
  plot(f,x1,x2,g,s1,s2,png="vs30g")
  plot(None,None,None,g,s1,s2,dmask,png="vs30ge")
  return
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
      q = gridder.grid(s1,s2)
      q = mul(mask,q)
      print name
      plot(f,x1,x2,q,s1,s2,png=name)

def cleanGeologyImage(g):
  n1,n2 = len(g[0]),len(g)
  for i2 in range(n2):
    for i1 in range(n1):
      gi = g[i2][i1]
      bad = gi<149.0 or gi>151.0
      bad = bad and gi<269.0 or gi>271.0
      bad = bad and gi<559.0 or gi>561.0
      bad = bad and gi<1129.0 or gi>1131.0
      if bad:
        g[i2][i1] = 0.0

def goSlopes():
  n1,n2 = 310,490
  s1 = Sampling(n1,0.00833333333,119.7-360.0)
  s2 = Sampling(n2,0.00833333333,21.75)
  f,x1,x2 = readScattered("vs30MeasuredTaiwan.txt",s1,s2)
  g = readImage("vs30SlopesTaiwan.dat",n1,n2)
  mask = makeMask(g)
  d = makeTensors(mask,g)
  dmask = makeTensors(mask,g)
  maskTensors(mask,dmask)
  g = mul(mask,g)
  plot(f,x1,x2,g,s1,s2,png="vs30s")
  plot(None,None,None,g,s1,s2,dmask,png="vs30se")
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
      q = gridder.grid(s1,s2)
      q = mul(mask,q)
      print name
      plot(f,x1,x2,q,s1,s2,png=name)

def makeTensors(mask,f):
  lof = LocalOrientFilter(3.0)
  d = lof.applyForTensors(f)
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av)
  au = div(av,au) # isotropy
  av = sub(1.0,au) # linearity
  am = mul(sub(1.0,mask),fillfloat(1.0,n1,n2))
  au = add(mul(mask,au),am)
  av = add(mul(mask,av),am)
  au = pow(au,2.0)
  d.setEigenvalues(au,av)
  return d

def maskTensors(mask,d):
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av)
  au = mul(mask,au)
  av = mul(mask,av)
  d.setEigenvalues(au,av)

def plot(f,x1,x2,g,s1,s2,d=None,png=None):
  sp = SimplePlot()
  sp.setHLabel("Longitude (degrees)")
  sp.setVLabel("Latitude (degrees)")
  sp.plotPanel.setHInterval(1.0)
  sp.plotPanel.setVInterval(1.0)
  sp.setFontSizeForPrint(8,240)
  #sp.setSize(670,1000) # without colorbar
  #sp.setSize(866,1044) # for slopes
  sp.setSize(865,920) # for geology
  sp.addColorBar("Vs30 (m/s)")
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setClips(0.0,760.0)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setLineColor(Color.WHITE)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(30)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if f and x1 and x2:
    mv = sp.addPoints(x1,x2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkColor(Color.WHITE)
    mv.setMarkSize(4)
  if pngDir and png:
    sp.paintToPng(600,3,pngDir+png+".png")

#############################################################################
dataDir = "/data/earth/vs30/"

def makeMask(g):
  gnull = min(g)
  n1,n2 = len(g[0]),len(g)
  m = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if g[i2][i1]!=gnull:
        m[i2][i1] = 1.0
  return m

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(dataDir+fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(x)
  ais.close()
  #x = flip1(x)
  x = flip2(x)
  return x

def readScattered(fileName,s1,s2):
  x1min,x1max = s1.first,s1.last
  x2min,x2max = s2.first,s2.last
  s = Scanner(FileInputStream(dataDir+fileName))
  f,x1,x2 = [],[],[]
  while s.hasNextLine():
    line = s.nextLine()
    data = line.split()
    x1i = float(data[0])-360.0
    x2i = float(data[1])
    fi = float(data[2])
    if x1min<=x1i and x1i<=x1max and x2min<=x2i and x2i<=x2max:
      x1.append(x1i)
      x2.append(x2i)
      f.append(fi)
  s.close()
  return f,x1,x2

def flip1(x):
  x = copy(x)
  n2 = len(x)
  for i2 in range(n2):
    x[i2] = reverse(x[i2])
  return x

def flip2(x):
  x = copy(x)
  n2 = len(x)
  for i2 in range(n2/2):
    j2 = n2-1-i2
    xi2 = x[i2]
    x[i2] = x[j2]
    x[j2] = xi2
  return x

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
run(main)
