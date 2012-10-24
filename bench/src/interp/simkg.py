#############################################################################
# Simulation and kriging gridder for smooth covariance function.

from shared import *

s1 = Sampling(431,1.0,-240.0) # Easting (km)
s2 = Sampling(291,1.0,  10.0) # Northing (km)
n1,d1,f1 = s1.count,s1.delta,s1.first
n2,d2,f2 = s2.count,s2.delta,s2.first

def main(args):
  demoTopo()
  demoTensors()
  #demoBlended()
  #demoKriging()

def demoTensors():
  z = makeTopo(s1,s2)
  et = makeTensors(z)
  plot(None,None,None,s1,s2,z,mv=False,et=et)

def demoTopo():
  z = makeTopo(s1,s2)
  plotTopo(z)

def makeTopo(s1,s2):
  z = mul(1.0e4*s1.delta,pow(randfloat(Random(3),n1,n2),8.0))
  ref = RecursiveExponentialFilter(2.0)
  for i in range(5):
    ref.apply(z,z)
  z = sub(z,min(z))
  return z

def makeTensors(z):
  n1,n2 = len(z[0]),len(z)
  et = EigenTensors2(n1,n2)
  for i2 in range(n2):
    i2m = max(i2-1,0)
    i2p = min(i2+1,n2-1)
    for i1 in range(n1):
      i1m = max(i1-1,0)
      i1p = min(i1+1,n1-1)
      z1 = (z[i2][i1p]-z[i2][i1m])/(i1p-i1m)
      z2 = (z[i2p][i1]-z[i2m][i1])/(i2p-i2m)
      den = 1.0+z1*z1+z2*z2
      d11 = (1.0+z2*z2)/den
      d22 = (1.0+z1*z1)/den
      d12 = -z1*z2/den
      et.setTensor(i1,i2,(d11,d12,d22))
  return et

def demoVariogram():
  f,x1,x2 = readWolfcamp()
  ft = fitTrend(f,x1,x2); f = subTrend(ft,f,x1,x2)
  #fa = sum(f)/len(f); f = sub(f,fa)
  hs,vs,hb,vb,hf,vf = makeVariogram(f,x1,x2)
  plotVariogram(hs,vs,None,None,None,None,png="vars")
  plotVariogram(hs,vs,hb,vb,None,None,png="varsb")
  plotVariogram(hs,vs,hb,vb,hf,vf,png="varsbf")

def demoSimple():
  f,x1,x2 = readWolfcamp()
  gsa = gridSimple(f,x1,x2,s1a,s2a);
  gsb = gridSimple(f,x1,x2,s1b,s2b);
  plot(f,x1,x2,s1a,s2a,gsa,"gsa",cv=False)
  plot(f,x1,x2,s1b,s2b,gsb,"gsb",cv=False)

def demoKriging():
  print "demoKriging ..."
  sd,rm = 0.0,80.0
  f,x1,x2 = readWolfcamp()
  f,x1,x2 = gridWolfcamp(f,x1,x2,s1,s2)
  #sd = mul(sd,randfloat(len(f)))
  #gk0 = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm)
  gkip = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=eti,pa=True)
  #gkit = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=eti,pa=False)
  #gkap = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=eta,pa=True)
  #gkat = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=eta,pa=False)
  #gkbp = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=etb,pa=True)
  #gkbt = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=etb,pa=False)
  #gksp = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=ets,pa=True)
  #gkst = gridKriging(f,x1,x2,s1,s2,sigmaD=sd,rangeM=rm,et=ets,pa=False)
  #plot(f,x1,x2,s1,s2,gk0,"gk0")
  #plot(f,x1,x2,s1,s2,gki,"gki")
  #plot(f,x1,x2,s1,s2,gka,"gka")
  #plot(f,x1,x2,s1,s2,gk0,"gkt0",et=eti)
  #plot(f,x1,x2,s1,s2,gkip,"gktip",et=eti)
  #plot(f,x1,x2,s1,s2,gkit,"gktit",et=eti)
  #plot(f,x1,x2,s1,s2,gkap,"gktap",et=eta)
  #plot(f,x1,x2,s1,s2,gkat,"gktat",et=eta)
  #plot(f,x1,x2,s1,s2,gkbp,"gktbp",et=etb)
  #plot(f,x1,x2,s1,s2,gkbt,"gktbt",et=etb)
  #plot(f,x1,x2,s1,s2,gksp,"gktsp",et=ets)
  #plot(f,x1,x2,s1,s2,gkst,"gktst",et=ets)
  #plot3(f,x1,x2,s1,s2,gk0)
  plot3(f,x1,x2,s1,s2,gkip)
  #plot3(f,x1,x2,s1,s2,gkit)
  #plot3(f,x1,x2,s1,s2,gkap)
  #plot3(f,x1,x2,s1,s2,gkat)
  #plot3(f,x1,x2,s1,s2,gkbp)
  #plot3(f,x1,x2,s1,s2,gkbt)
  #plot3(f,x1,x2,s1,s2,gksp)
  #plot3(f,x1,x2,s1,s2,gkst)
  print "done"

def demoBlended():
  f,x1,x2 = readWolfcamp()
  f,x1,x2 = gridWolfcamp(f,x1,x2,s1,s2)
  #gbi = gridBlended(f,x1,x2,s1,s2,0.5,et=eti)
  #gba = gridBlended(f,x1,x2,s1,s2,0.5,et=eta)
  #gbb = gridBlended(f,x1,x2,s1,s2,0.5,et=etb)
  #gbr = gridBlended(f,x1,x2,s1,s2,0.5,et=etr)
  gbs = gridBlended(f,x1,x2,s1,s2,0.5,et=ets)
  #plot(f,x1,x2,s1,s2,gbi,"gbi")
  #plot(f,x1,x2,s1,s2,gba,"gba")
  #plot(f,x1,x2,s1,s2,gbb,"gbb")
  #plot(f,x1,x2,s1,s2,gbr,"gbr")
  #plot(f,x1,x2,s1,s2,gbs,"gbs")
  #plot(f,x1,x2,s1,s2,gbi,"gbit",et=eti)
  #plot(f,x1,x2,s1,s2,gba,"gbat",et=eta)
  #plot(f,x1,x2,s1,s2,gbb,"gbbt",et=etb)
  #plot(f,x1,x2,s1,s2,gbr,"gbrt",et=etr)
  #plot(f,x1,x2,s1,s2,gbs,"gbst",et=ets)
  #plot3(f,x1,x2,s1,s2,gbi)
  #plot3(f,x1,x2,s1,s2,gba)
  #plot3(f,x1,x2,s1,s2,gbb)
  #plot3(f,x1,x2,s1,s2,gbr)
  plot3(f,x1,x2,s1,s2,gbs)

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
  #vf = mul(0.1,sub(1.0,exp(neg(div(mul(hf,hf),2.0*220.0*220.0)))))
  return hs,vs,hb,vb,hf,vf

def gridSimple(f,x1,x2,s1,s2):
  return SimpleGridder2(f,x1,x2).grid(s1,s2);

def gridBlended(f,x1,x2,s1,s2,smooth=0.5,et=None):
  bg = BlendedGridder2(f,x1,x2)
  bg.setSmoothness(smooth)
  if et:
    bg.setTensors(et)
    #d = zerofloat(3); et.getTensor(0,0,d); bg.setTensor(d[0],d[1],d[2])
  return bg.grid(s1,s2);

def gridKriging(f,x1,x2,s1,s2,sigmaD=0.0,rangeM=40.0,et=None,pa=False):
  cm = SmoothCovariance(0.5,1.0,rangeM,2)
  cm.setTensors(et)
  #cm.testSpd(s1.count,s2.count)
  kg = KrigingGridder2(f,x1,x2)
  kg.setDataError(sigmaD)
  kg.setModelCovariance(cm);
  kg.setPolyTrend(0)
  kg.setPaciorek(pa)
  #kg.setIdentityTensors()
  kg.setTensors(et)
  return kg.grid(s1,s2)

#############################################################################
# Plotting

_pngDir = None # directory to use for png files
#_pngDir = "../../png/interp" # directory to use for png files
 
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

def plot(f,x1,x2,s1,s2,g,png=None,cv=True,mv=True,et=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setSize(1083,648)
  sp.setFontSizeForSlide(1.0,1.0)
  sp.setHLabel("Easting (km)")
  sp.setVLabel("Northing (km)")
  sp.setLimits(-241.0,9.0,191.0,301.0)
  sp.addColorBar("Elevation (km)")
  sp.plotPanel.setVInterval(100)
  pv = sp.addPixels(s1,s2,g)
  #pv.setClips(0.2,1.2)
  pv.setColorModel(ColorMap.JET)
  #pv.setColorModel(makeTransparentColorModel())
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cv and not et:
    cv = sp.addContours(s1,s2,g)
    cv.setLineColor(Color.BLACK)
    #cv.setContours(Sampling(26,0.04,0.2))
    cv.setContours(Sampling(21,0.05,0.2))
  if mv:
    mv = sp.addPoints(x1,x2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkSize(5)
  if et:
    tv = TensorsView(s1,s2,et)
    sp.add(tv)
  if png and _pngDir: sp.paintToPng(300,6,_pngDir+png+".png")

def makeTransparentColorModel():
  a = fillfloat(1.0,256); a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a);

def plot3(f,x1,x2,s1,s2,g):
  f = mul(120.0,f)
  g = mul(120.0,g)
  pg = makePointGroup(x1,x2,f)
  tg = makeTriangleGroup(s1,s2,g)
  world = World()
  world.addChild(pg)
  world.addChild(tg)
  sf = SimpleFrame(world,AxesOrientation.XOUT_YRIGHT_ZUP)
  sf.setSize(1200,750)
  sf.setWorldSphere(-240.0,10.0,30.0,170.0,300.0,80.0)
  ov = sf.getOrbitView()
  ov.setScale(1.5)

def plotTopo(s1,s2,z):
  tg = makeTriangleGroup(s1,s2,z)
  world = World()
  world.addChild(tg)
  sf = SimpleFrame(world,AxesOrientation.XOUT_YRIGHT_ZUP)
  sf.setSize(1200,750)
  #sf.setWorldSphere(-240.0,10.0,30.0,170.0,300.0,80.0)
  ov = sf.getOrbitView()
  ov.setScale(2.0)

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

def gridWolfcamp(f,x1,x2,s1,s2):
  return SimpleGridder2.samplesOnGrid(s1,s2,f,x1,x2)

#############################################################################
if __name__ == "__main__":
  run(main)
