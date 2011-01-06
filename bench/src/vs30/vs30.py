from shared import *

def main(args):
  goFirst()

def goFirst():
  s1 = Sampling(310,0.00833333333,119.7-360.0)
  s2 = Sampling(490,0.00833333333,21.75)
  n1,n2 = s1.count,s2.count
  g = readImage("vs30SlopesTaiwan.dat",n1,n2)
  mask = makeMask(g)
  f,x1,x2 = readScattered("vs30MeasuredTaiwan.txt",s1,s2)
  d = makeTensors(g)
  v = gridBlended(d,f,x1,x2,s1,s2)
  v = mul(mask,v)
  sp = SimplePlot()
  sp.setHLabel("Longitude (degrees)")
  sp.setVLabel("Latitude (degrees)")
  sp.setFontSizeForPrint(8,240)
  sp.setSize(866,1044)
  sp.addColorBar("Vs30 (m/s)")
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.getGray(0.45,0.55))
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv = sp.addPixels(s1,s2,v)
  pv.setColorModel(makeTransparentColorModel(0.8))
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setClips(10.0,760.0)
  """
  tv = TensorsView(s1,s2,d)
  tv.setLineColor(Color.RED)
  tv.setLineWidth(3)
  tv.setEllipsesDisplayed(30)
  tile = sp.plotPanel.getTile(0,0)
  tile.addTiledView(tv)
  """
  mv = sp.addPoints(x1,x2)
  mv.setLineStyle(PointsView.Line.NONE)
  mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  mv.setMarkColor(Color.WHITE)
  mv.setMarkSize(4)
  sp.paintToPng(600,3,"junk.png")


#############################################################################
dataDir = "/data/earth/vs30/"

def gridBlended(d,f,x1,x2,s1,s2):
  #bg = BlendedGridder2(f,x1,x2)
  bg = BlendedGridder2(d,f,x1,x2)
  #bg.setBlending(False)
  #bg.setSmoothness(1.0)
  g = bg.grid(s1,s2)
  lsf = LocalSmoothingFilter()
  lsf.applySmoothS(g,g)
  return g

def makeMask(g):
  gnull = min(g)
  n1,n2 = len(g[0]),len(g)
  m = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if g[i2][i1]!=gnull:
        m[i2][i1] = 1.0
  return m

def makeTensors(g):
  mask = makeMask(g)
  g = sub(g,min(g))
  g = sin(mul(2.0*PI/max(g),g))
  lof = LocalOrientFilter(4.0)
  d = lof.applyForTensors(g)
  #invertEigenvalues(mask,d)
  lsf = LocalSemblanceFilter(2,2)
  s1 = lsf.semblance(LocalSemblanceFilter.Direction2.V,d,g)
  s2 = lsf.semblance(LocalSemblanceFilter.Direction2.UV,d,g)
  s2 = pow(s2,3.0)
  s1 = mul(mask,s1)
  s2 = mul(mask,s2)
  s1 = clip(0.001,1.0,s1)
  s2 = clip(0.001,1.0,s2)
  d.setEigenvalues(s2,s1)
  return d

def invertEigenvalues(mask,d):
  n1,n2 = d.n1,d.n2
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  d.getEigenvalues(au,av)
  for i2 in range(n2):
    for i1 in range(n1):
      if mask[i2][i1]>0.0:
        au[i2][i1] = 1.0/au[i2][i1]
        av[i2][i1] = 1.0/av[i2][i1]
      else:
        au[i2][i1] = 1.0
        av[i2][i1] = 1.0
  aumax = max(au)
  avmax = max(av)
  amax = max(aumax,avmax)
  au = mul(1.0/amax,au)
  av = mul(1.0/amax,av)
  au = pow(au,2.0)
  av = pow(av,2.0)
  au = clip(0.001,1.0,au)
  av = clip(0.001,1.0,av)
  d.setEigenvalues(au,av)

def subLocalMean(x):
  y = copy(x)
  rgf = RecursiveGaussianFilter(20.0)
  rgf.apply00(x,y)
  y = sub(x,y)
  return y

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
