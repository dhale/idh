#############################################################################
# Demo interpolation of Wolfcamp Aquifer potentiometric levels

from shared import *

#s1 = Sampling(271,1.0,-150.0) # Easting (miles)
#s2 = Sampling(201,1.0,0.0) # Northing (miles)
s1 = Sampling(431,1.0,-240.0) # Easting (km)
s2 = Sampling(291,1.0,  10.0) # Northing (km)
s1a = Sampling(87,5.0,-240.0)
s2a = Sampling(59,5.0,  10.0)
s1b = Sampling(44,10.0,-240.0)
s2b = Sampling(30,10.0,  10.0)

def main(args):
  #demoVariogram()
  demoMethods()

def demoVariogram():
  f,x1,x2 = readWolfcamp()
  ft = fitTrend(f,x1,x2)
  f = subTrend(ft,f,x1,x2)
  hs,vs,hb,vb,hf,vf = makeVariogram(f,x1,x2)
  plotVariogram(hs,vs,None,None,None,None,png="vars")
  plotVariogram(hs,vs,hb,vb,None,None,png="varsb")
  plotVariogram(hs,vs,hb,vb,hf,vf,png="varsbf")

def makeVariogram(f,x1,x2):
  n = len(f)
  nn = n*(n-1)/2
  hs = zerofloat(nn)
  vs = zerofloat(nn)
  k = 0
  for i in range(n):
    for j in range(i+1,n):
      df = f[i]-f[j]
      d1 = x1[i]-x1[j]
      d2 = x2[i]-x2[j]
      hs[k] = sqrt(d1*d1+d2*d2)
      vs[k] = 0.5*df*df
      k += 1
  hmax = max(hs)
  nh = n
  dh = hmax/(nh-1)
  fh = 0.0
  sh = Sampling(nh,dh,fh)
  cb = fillfloat(0.001,nh)
  vb = zerofloat(nh)
  for k in range(nn):
    j = sh.indexOfNearest(hs[k])
    cb[j] += 1.0
    vb[j] += vs[k]
  hb = rampfloat(fh,dh,nh)
  vb = div(vb,cb)
  hf = copy(hb)
  vf = mul(0.004,sub(1.0,exp(neg(div(hf,40.0)))))
  return hs,vs,hb,vb,hf,vf

def demoMethods():
  f,x1,x2 = readWolfcamp()
  gba = gridSimple(f,x1,x2,s1a,s2a);
  gsk = gridKriging(f,x1,x2,s1,s2,sigmaD=0.0,shapeM=0.5,rangeM=80.0)
  gb0 = gridBlended(f,x1,x2,s1,s2,0.5);
  """
  gbb = gridSimple(f,x1,x2,s1b,s2b);
  gnn = gridNearest(f,x1,x2,s1,s2);
  gsk = gridKriging(f,x1,x2,s1,s2,sigmaD=0.0,shapeM=1.0,rangeM=40.0)
  gsh = gridShepard(f,x1,x2,s1,s2);
  gb0 = gridBlended(f,x1,x2,s1,s2,0.5);
  gb1 = gridBlended(f,x1,x2,s1,s2,1.0);
  gs0 = gridSibson(f,x1,x2,s1,s2,False);
  gs1 = gridSibson(f,x1,x2,s1,s2,True);
  """
  plot(f,x1,x2,s1a,s2a,gba,"gba",False)
  plot(f,x1,x2,s1,s2,gsk,"gsk")
  plot(f,x1,x2,s1,s2,gb0,"gb0")
  """
  plot(f,x1,x2,s1b,s2b,gbb,"gbb",False)
  plot(f,x1,x2,s1,s2,gnn,"gnn",False)
  plot(f,x1,x2,s1,s2,gsk,"gsk")
  plot(f,x1,x2,s1,s2,gsh,"gsh")
  plot(f,x1,x2,s1,s2,gb0,"gb0")
  plot(f,x1,x2,s1,s2,gb1,"gb1")
  plot(f,x1,x2,s1,s2,gs0,"gs0")
  plot(f,x1,x2,s1,s2,gs1,"gs1")
  """
  plot3(f,x1,x2,s1,s2,gsk)
  plot3(f,x1,x2,s1,s2,gb0)
  """
  plot3(f,x1,x2,s1,s2,gnn)
  plot3(f,x1,x2,s1,s2,gsk)
  plot3(f,x1,x2,s1,s2,gsh)
  plot3(f,x1,x2,s1,s2,gb0)
  plot3(f,x1,x2,s1,s2,gb1)
  plot3(f,x1,x2,s1,s2,gs0)
  plot3(f,x1,x2,s1,s2,gs1)
  """

def gridSimple(f,x1,x2,s1,s2):
  return SimpleGridder2(f,x1,x2).grid(s1,s2);

def gridBlended(f,x1,x2,s1,s2,smooth=0.5):
  bg = BlendedGridder2(f,x1,x2)
  bg.setSmoothness(smooth)
  return bg.grid(s1,s2)
  """
  sg = SimpleGridder2(f,x1,x2)
  pnull = -FLT_MAX
  tnull = -FLT_MAX
  sg.setNullValue(pnull)
  p = sg.grid(s1,s2);
  n1,n2 = s1.count,s2.count
  t = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if p[i2][i1]==pnull:
        t[i2][i1] = tnull
      else:
        t[i2][i1] = 0.0
  bg.gridNearest(t,p)
  q = zerofloat(n1,n2)
  bg.gridBlended(t,p,q)
  bg.gridBlended(t,copy(q),q)
  return q
  """

def gridKriging(f,x1,x2,s1,s2,sigmaD=0.0,shapeM=1.0,rangeM=40.0):
  kg = KrigingGridder2(f,x1,x2)
  kg.setDataError(sigmaD)
  kg.setModelCovariance(Matern(shapeM,1.0,rangeM));
  kg.setPolyTrend(1);
  return kg.grid(s1,s2)

def gridNearest(f,x1,x2,s1,s2):
  return NearestGridder2(f,x1,x2).grid(s1,s2)

def gridSibson(f,x1,x2,s1,s2,smooth=False):
  sg = SibsonGridder2(f,x1,x2)
  sg.setSmooth(smooth)
  return sg.grid(s1,s2)

def gridDiscreteSibson(f,x1,x2,s1,s2,nsmooth=0):
  sg = DiscreteSibsonGridder2(f,x1,x2)
  sg.setSmooth(nsmooth)
  return sg.grid(s1,s2)

#############################################################################
# Plotting

_pngDir = None # directory to use for png files
#_pngDir = "png/" # directory to use for png files
 
def plotVariogram(hs,vs,hb=None,vb=None,hf=None,vf=None,png=None):
  sp = SimplePlot()
  sp.setSize(1083,750)
  #sp.setFontSizeForSlide(1.0,1.0)
  sp.setFontSize(41)
  sp.setVLimits(0,0.015)
  sp.setHLimits(0,400.0)
  sp.setHLabel("Distance h (km)")
  sp.setVLabel("Semi-variance (km squared)")
  pv = sp.addPoints(hs,vs)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(2.0)
  if hb and vb:
    pv = sp.addPoints(hb,vb)
    pv.setLineStyle(PointsView.Line.NONE)
    #pv.setLineWidth(3.0)
    #pv.setLineColor(Color.BLUE)
    pv.setMarkColor(Color.BLUE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(10.0)
  if hf and vf:
    pv = sp.addPoints(hf,vf)
    pv.setLineStyle(PointsView.Line.DASH)
    pv.setLineWidth(3.0)
    pv.setLineColor(Color.RED)
  if png and _pngDir: sp.paintToPng(300,6,_pngDir+png+".png")

def plot(f,x1,x2,s1,s2,g,png=None,cv=True,mv=True):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setSize(1083,648)
  sp.setFontSizeForSlide(1.0,1.0)
  sp.setHLabel("Easting (km)")
  sp.setVLabel("Northing (km)")
  sp.setLimits(-241.0,9.0,191.0,301.0)
  sp.addColorBar("Potentiometric level (km)")
  sp.plotPanel.setVInterval(100)
  pv = sp.addPixels(s1,s2,g)
  pv.setClips(0.2,1.2)
  #pv.setColorModel(ColorMap.JET)
  pv.setColorModel(makeTransparentColorModel())
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cv:
    cv = sp.addContours(s1,s2,g)
    cv.setLineColor(Color.BLACK)
    #cv.setContours(Sampling(26,0.04,0.2))
    cv.setContours(Sampling(11,0.1,0.2))
  if mv:
    mv = sp.addPoints(x1,x2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkSize(5)
  if png and _pngDir: sp.paintToPng(300,6,_pngDir+png+".png")

def makeTransparentColorModel():
  cm = ColorMap.getJet()
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

def plot3(f,x1,x2,s1,s2,g):
  f = mul(300.0,f)
  g = mul(300.0,g)
  pg = makePointGroup(x1,x2,f)
  tg = makeTriangleGroup(s1,s2,g)
  world = World()
  world.addChild(pg)
  world.addChild(tg)
  sf = SimpleFrame(world,AxesOrientation.XOUT_YRIGHT_ZUP)
  sf.setSize(1200,750)
  sf.setWorldSphere(-240.0,10.0,90.0,170.0,300.0,240.0)
  ov = sf.getOrbitView()
  ov.setScale(1.5)
  ov.setAzimuthAndElevation(50.0,12.0)
  ov.setTranslate(Vector3(0.070,0.000,0.087))

def makePointGroup(x,y,z):
  n = len(z)
  xyz = zerofloat(3*n)
  copy(n,0,1,x,0,3,xyz)
  copy(n,0,1,y,1,3,xyz)
  copy(n,0,1,z,2,3,xyz)
  pg = PointGroup(3.0,xyz)
  ss = StateSet.forTwoSidedShinySurface(Color.RED)
  pg.setStates(ss)
  return pg

def makeTriangleGroup(sx,sy,sz):
  sz = transpose(sz)
  tg = TriangleGroup(True,sx,sy,sz)
  ss = StateSet.forTwoSidedShinySurface(Color.LIGHT_GRAY)
  tg.setStates(ss)
  return tg

#############################################################################
# Data

def readWolfcamp():
  s = Scanner(FileInputStream("wa.txt"))
  for line in range(6):
    s.nextLine()
  n = 85
  f = zerofloat(n)
  x1 = zerofloat(n)
  x2 = zerofloat(n)
  for i in range(n):
    x1[i] = s.nextFloat()
    x2[i] = s.nextFloat()
    f[i] = s.nextFloat()
  s.close()
  f = mul(0.3048,f)
  x1 = mul(1.609344,x1)
  x2 = mul(1.609344,x2)
  #print "f min =",min(f)," max =",max(f)
  #print "x1 min =",min(x1)," max =",max(x1)
  #print "x2 min =",min(x2)," max =",max(x2)
  return f,x1,x2

#############################################################################
if __name__ == "__main__":
  run(main)
