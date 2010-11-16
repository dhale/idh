#############################################################################
# Test for SPD co-variance function with different distance metrics.

from shared import *

s1 = Sampling(501,0.002,0.0);
s2 = Sampling(501,0.002,0.0);

def main(args):
  #demoTensor()
  #demoLinearVelocity()
  demoLinearSloth()
  #demoCovariance()

def demoTensor():
  n1,n2 = s1.count,s2.count
  #et = WavyTensors2(n1,n2,1000.0,45.0)
  et = SinSinTensors2(n1,n2,0.01,0.1)
  #et = RandomTensors2(n1,n2)
  #et = BlockTensors2(n1,n2)
  n = 4
  x1 = [0.0,1.0,0.0,1.0]
  x2 = [0.0,0.0,1.0,1.0]
  #x1,x2 = makeRandomPoints(n)
  tm = makeTimeMaps(x1,x2,s1,s2,et)
  tm = mul(s1.delta,tm)
  for i in range(n):
    plot(None,x1,x2,s1,s2,tm[i],et=et)
  sigma = 1.0
  delta = 8.0
  c = DMatrix(n,n)
  t = DMatrix(n,n)
  for i in range(n):
    for j in range(n):
      j1 = s1.indexOfNearest(x1[j])
      j2 = s2.indexOfNearest(x2[j])
      tj = tm[i][j2][j1]
      c.set(i,j,covExp(sigma,delta,tj))
      #c.set(i,j,covGauss(sigma,delta,tj))
      t.set(i,j,tj)
  print "t =\n",t
  print "c =\n",c
  evd = DMatrixEvd(c)
  dump(evd.getRealEigenvalues())

def demoCovariance():
  #ss0,ss1,ss2 = 1.0,0.0,0.0
  ss0,ss1,ss2 = 1.0,3.0,-0.5
  #v0,v1,v2 = 1.0,5.0,5.0
  #v0,v1,v2 = 1.0,-0.4,-0.4
  v0,v1,v2 = 1.0,4.0,-0.4
  c = DMatrix(4,4)
  t = DMatrix(4,4)
  x1 = [0.0,1.0,0.0,1.0]
  x2 = [0.0,0.0,1.0,1.0]
  sigma = 1.0
  delta = sqrt(2.0)
  for i in range(4):
    for j in range(4):
      tij = Eikonal.timeForLinearSloth(ss0,ss1,ss2,x1[i],x2[i],x1[j],x2[j])
      #tij = Eikonal.timeForLinearVelocity(v0,v1,v2,x1[i],x2[i],x1[j],x2[j])
      #c.set(i,j,covExp(sigma,delta,tij))
      c.set(i,j,covGauss(sigma,delta,tij))
      t.set(i,j,tij)
  print "t =\n",t
  print "c =\n",c
  evd = DMatrixEvd(c)
  dump(evd.getRealEigenvalues())

def demoLinearVelocity():
  #v0,v1,v2 = 1.0,FLT_EPSILON,FLT_EPSILON
  #v0,v1,v2 = 1.0,3.0,3.0
  v0,v1,v2 = 1.0,4.0,-0.4
  for p1 in [0.0,1.0]:
    for p2 in [0.0,1.0]:
      t = Eikonal.timesForLinearVelocity(v0,v1,v2,p1,p2,s1,s2)
      plot(None,None,None,s1,s2,t,"t00")
      print "p1 =",p1,"p2 =",p2,"tmax =",max(t)

def demoLinearSloth():
  #ss0,ss1,ss2 = 1.0,0.0,0.0
  #ss0,ss1,ss2 = 1.0,3.0,3.0
  ss0,ss1,ss2 = 1.0,4.0,4.0
  for p1 in [0.0,1.0]:
    for p2 in [0.0,1.0]:
      t = Eikonal.timesForLinearSloth(ss0,ss1,ss2,p1,p2,s1,s2)
      plot(None,None,None,s1,s2,t,"t00")
      print "p1 =",p1,"p2 =",p2,"tmax =",max(t)

def makeTimeMaps(x1,x2,s1,s2,et):
  nk = len(x1)
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
  return t

def makeRandomPoints(n):
  r = Random(314159)
  x1 = zerofloat(n)
  x2 = zerofloat(n)
  for i in range(n):
    x1[i] = s1.getValue(s1.indexOfNearest(r.nextFloat()))
    x2[i] = s2.getValue(s2.indexOfNearest(r.nextFloat()))
  return x1,x2

def computeTimesForLinearVelocity(v0,v1,v2,p1,p2,s1,s2):
  return Eikonal.solveForLinearVelocity(v0,v1,v2,p1,p2,s1,s2)

def computeTimesForLinearSloth(ss0,ss1,ss2,p1,p2,s1,s2):
  return Eikonal.solveForLinearSloth(ss0,ss1,ss2,p1,p2,s1,s2)

def covExp(sigma,delta,d):
  return sigma*sigma*exp(-d/delta)

def covGauss(sigma,delta,d):
    return sigma*sigma*exp(-0.5*d*d/(delta*delta))

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

class SinSinTensors2(EigenTensors2):
  """
  2D tensors with a low-velocity hump in the middle.
  """
  def __init__(self,n1,n2,vmin,power):
    EigenTensors2.__init__(self,n1,n2)
    s1 = Sampling(n1,1.0/(n1-1),0.0)
    s2 = Sampling(n2,1.0/(n2-1),0.0)
    d = zerofloat(3)
    vmins = vmin*vmin
    va = zerofloat(n1,n2)
    for i2 in range(n2):
      sin2 = sin(PI*s2.getValue(i2))
      for i1 in range(n1):
        sin1 = sin(PI*s1.getValue(i1))
        vs = 1.0-(1.0-vmins)*pow(sin1*sin2,power)
        d[0] = vs
        d[1] = 0.0
        d[2] = vs
        self.setTensor(i1,i2,d)
        va[i2][i1] = vs

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

class WavyTensors2(EigenTensors2):
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

def plot(f,x1,x2,s1,s2,g,png=None,cv=True,mv=True,et=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setSize(1000,800)
  sp.setFontSizeForSlide(1.0,1.0)
  #sp.plotPanel.setVInterval(100)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.JET)
  #pv.setColorModel(makeTransparentColorModel())
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  sp.addColorBar("Time")
  if cv and not et:
    cv = sp.addContours(s1,s2,g)
    cv.setLineColor(Color.BLACK)
  if mv and x1 and x2:
    mv = sp.addPoints(x1,x2)
    mv.setLineStyle(PointsView.Line.NONE)
    mv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    mv.setMarkSize(10)
    mv.setMarkColor(Color.WHITE)
  if et:
    tv = TensorsView(s1,s2,et)
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
