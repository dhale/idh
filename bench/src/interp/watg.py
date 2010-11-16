#############################################################################
# Demo tensor-guided interpolation of Wolfcamp Aquifer potentiometric levels

from shared import *

#s1 = Sampling(271,1.0,-150.0) # Easting (miles)
#s2 = Sampling(201,1.0,0.0) # Northing (miles)
s1 = Sampling(431,1.0,-240.0) # Easting (km)
s2 = Sampling(291,1.0,  10.0) # Northing (km)
s1t,s2t = s1,s2
s1a = Sampling(87,5.0,-240.0)
s2a = Sampling(59,5.0,  10.0)
s1b = Sampling(44,10.0,-240.0)
s2b = Sampling(30,10.0,  10.0)
eti,eta,etb,etr,ets = None,None,None,None,None

def main(args):
  makeTensors()
  #demoTimeMaps()
  #demoVariogram()
  demoSimple()
  #demoBlended()
  #demoKriging()

def makeTensors():
  global eti,eta,etb,etr,ets
  eti = ConstantTensors2(s1.count,s2.count,-30,1.0,1.0)
  eta = ConstantTensors2(s1.count,s2.count,-30,0.1,1.0)
  etb = BlockTensors2(s1.count,s2.count)
  etr = RandomTensors2(s1.count,s2.count)
  ets = SineTensors2(s1.count,s2.count,100.0,30.0)

def demoTimeMaps():
  f,x1,x2 = readWolfcamp()
  f,x1,x2 = gridWolfcamp(f,x1,x2,s1,s2)
  #makeTimeMaps(f,x1,x2,eti,"tmi.dat")
  #makeTimeMaps(f,x1,x2,eta,"tma.dat")
  #makeTimeMaps(f,x1,x2,etb,"tmb.dat")
  #makeTimeMaps(f,x1,x2,etr,"tmr.dat")
  #makeTimeMaps(f,x1,x2,ets,"tms.dat")

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
  sigma,delta = 1.0,40.0 # exponential model after linear trend removed
  #sigma,delta = 1.0,320.0 # exponential model after linear trend removed
  #sigma,delta = 1.0,220.0 # gaussian model after average removed
  f,x1,x2 = readWolfcamp()
  f,x1,x2 = gridWolfcamp(f,x1,x2,s1,s2)
  #gk0 = gridKriging(sigma,delta,f,x1,x2,s1,s2)
  gki = gridKrigingTm(sigma,delta,f,x1,x2,"tmi.dat")
  gka = gridKrigingTm(sigma,delta,f,x1,x2,"tma.dat")
  #gkb = gridKrigingTm(sigma,delta,f,x1,x2,"tmb.dat")
  #gkr = gridKrigingTm(sigma,delta,f,x1,x2,"tmr.dat")
  gks = gridKrigingTm(sigma,delta,f,x1,x2,"tms.dat")
  #plot(f,x1,x2,s1,s2,gk0,"gk0")
  plot(f,x1,x2,s1,s2,gki,"gki")
  plot(f,x1,x2,s1,s2,gka,"gka")
  #plot(f,x1,x2,s1,s2,gkb,"gkb")
  #plot(f,x1,x2,s1,s2,gkr,"gkr")
  plot(f,x1,x2,s1,s2,gks,"gks")
  #plot(f,x1,x2,s1,s2,gk0,"gkt0",et=eti)
  plot(f,x1,x2,s1,s2,gki,"gkti",et=eti)
  plot(f,x1,x2,s1,s2,gka,"gkta",et=eta)
  #plot(f,x1,x2,s1,s2,gkr,"gktr",et=etr)
  plot(f,x1,x2,s1,s2,gks,"gkts",et=ets)
  #plot3(f,x1,x2,s1,s2,gki)
  #plot3(f,x1,x2,s1,s2,gka)
  #plot3(f,x1,x2,s1,s2,gkb)
  #plot3(f,x1,x2,s1,s2,gkr)
  #plot3(f,x1,x2,s1,s2,gks)

def demoBlended():
  f,x1,x2 = readWolfcamp()
  f,x1,x2 = gridWolfcamp(f,x1,x2,s1,s2)
  gbi = gridBlended(f,x1,x2,s1,s2,0.5,et=eti);
  gba = gridBlended(f,x1,x2,s1,s2,0.5,et=eta);
  #gbb = gridBlended(f,x1,x2,s1,s2,0.5,et=etb);
  #gbr = gridBlended(f,x1,x2,s1,s2,0.5,et=etr);
  gbs = gridBlended(f,x1,x2,s1,s2,0.5,et=ets);
  plot(f,x1,x2,s1,s2,gbi,"gbi")
  plot(f,x1,x2,s1,s2,gba,"gba")
  #plot(f,x1,x2,s1,s2,gbb,"gbb")
  #plot(f,x1,x2,s1,s2,gbr,"gbr")
  plot(f,x1,x2,s1,s2,gbs,"gbs")
  #plot(f,x1,x2,s1,s2,gbi,"gbit",et=eti)
  #plot(f,x1,x2,s1,s2,gba,"gbat",et=eta)
  #plot(f,x1,x2,s1,s2,gbb,"gbbt",et=etb)
  #plot(f,x1,x2,s1,s2,gbr,"gbrt",et=etr)
  #plot(f,x1,x2,s1,s2,gbs,"gbst",et=ets)
  #plot3(f,x1,x2,s1,s2,gbi)
  #plot3(f,x1,x2,s1,s2,gba)
  #plot3(f,x1,x2,s1,s2,gbb)
  #plot3(f,x1,x2,s1,s2,gbr)
  #plot3(f,x1,x2,s1,s2,gbs)

def makeTimeMaps(f,x1,x2,et,fileName):
  nk = len(f)
  n1,n2 = s1.count,s2.count
  tm = TimeMarker2(n1,n2,et)
  t = fillfloat(-1.0,n1,n2,nk)
  m = zeroint(n1,n2)
  for k in range(nk):
    i1 = s1.indexOfNearest(x1[k])
    i2 = s2.indexOfNearest(x2[k])
    print "k =",k," i1 =",i1," i2 =",i2
    t[k][i2][i1] = 0.0
    tm.apply(t[k],m)
  SimplePlot.asContours(t[1])
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(t)
  aos.close()

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
  return bg.grid(s1,s2);

def gridKrigingTm(sigma,delta,f,x1,x2,tmFileName):
  tm = zerofloat(s1.count,s2.count,len(f))
  ais = ArrayInputStream(tmFileName)
  ais.readFloats(tm)
  ais.close()
  return gridKriging(sigma,delta,f,x1,x2,s1,s2,tm)

#############################################################################
# Tensors

class BlockTensors2(EigenTensors2):
  """
  2D tensors constant within each of 2x2 blocks of samples
  """
  def __init__(self,n1,n2):
    EigenTensors2.__init__(self,n1,n2)
    m1,m2 = n1/2,n2/2
    du,dv = 0.01,1.00
    for i2 in range(n2):
      for i1 in range(n1):
        self.setEigenvalues(i1,i2,du,dv)
        if i1<m1 and i2<m2 or i1>=m1 and i2>=m2:
          self.setEigenvectorU(i1,i2,1.0,0.0)
        else:
          self.setEigenvectorU(i1,i2,0.0,1.0)

class ConstantTensors2(EigenTensors2):
  """
  2D constant tensors
  dip is angle clockwise between eigenvector v and x2 axis (in degrees)
  du is eigenvalue corresponding to eigenvector u
  dv is eigenvalue corresponding to eigenvector v
  """
  def __init__(self,n1,n2,dip,du,dv):
    EigenTensors2.__init__(self,n1,n2)
    a = dip*PI/180.0
    u1 =  cos(a)
    u2 = -sin(a)
    for i2 in range(n2):
      for i1 in range(n1):
        self.setEigenvalues(i1,i2,du,dv)
        self.setEigenvectorU(i1,i2,u1,u2)

class RandomTensors2(EigenTensors2):
  """
  2D random tensors
  """
  def __init__(self,n1,n2):
    EigenTensors2.__init__(self,n1,n2)
    r = Random(3) # 31
    for i2 in range(n2):
      for i1 in range(n1):
        a = 2.0*PI*r.nextFloat()
        u1,u2 = cos(a),sin(a)
        du,dv = 0.01+0.09*r.nextFloat(),0.01+0.99*r.nextFloat()
        self.setEigenvectorU(i1,i2,u1,u2)
        self.setEigenvalues(i1,i2,du,dv)

class SineTensors2(EigenTensors2):
  """
  2D tensors for geodesic distance in sinusoidal topography
  ampmax is max height of sine waves
  dipmax is maximum dip (in degrees) of sine waves
  """
  def __init__(self,n1,n2,ampmax,dipmax):
    EigenTensors2.__init__(self,n1,n2)
    b1 = 3.0*2.0*PI/(n1-1) # 3 cycles horizontally
    b2 = 3.0*2.0*PI/(n2-1) # 3 cycles vertically
    a1 = ampmax
    a2 = atan(dipmax*PI/180.0)/b2
    f = zerofloat(n1,n2)
    for i2 in range(n2):
      for i1 in range(n1):
        s2 = a2*sin(b2*i2)
        f[i2][i1] = a1*cos(b1*(i1+s2))
        e1 = -a1*b1*sin(b1*(i1+s2))
        e2 = e1*a2*b2*cos(b2*i2)
        den = 1.0+e1*e1+e2*e2
        d11 = (1.0+e2*e2)/den
        d22 = (1.0+e1*e1)/den
        d12 = -e1*e2/den
        self.setTensor(i1,i2,(d11,d12,d22))

#############################################################################
# Plotting

#_pngDir = None # directory to use for png files
_pngDir = "png/" # directory to use for png files
 
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
  sp.addColorBar("Potentiometric level (km)")
  sp.plotPanel.setVInterval(100)
  pv = sp.addPixels(s1,s2,g)
  pv.setClips(0.2,1.2)
  #pv.setColorModel(ColorMap.JET)
  pv.setColorModel(makeTransparentColorModel())
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
    tv = TensorsView(s1t,s2t,et)
    sp.add(tv)
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
  fnull = -999.9
  sg = SimpleGridder2(f,x1,x2)
  sg.setNullValue(fnull)
  g = sg.grid(s1,s2);
  f,x1,x2 = sg.getGriddedSamples(fnull,s1,s2,g)
  return f,x1,x2

#############################################################################
if __name__ == "__main__":
  run(main)
