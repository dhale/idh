import sys
#from math import *
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from ldf import BilateralFilter

gauss = BilateralFilter.Type.GAUSS
huber = BilateralFilter.Type.HUBER
tukey = BilateralFilter.Type.TUKEY

#pngDir = None
pngDir = "../../png/blf/"


#############################################################################

def main(args):
  #goRandomBlocks()
  #goRangeFunctions()
  #goApproximation()
  goCtRock()
  #goImage("f3d")
  #goFilter("f3d")
  #goFilter("tpd")
  #goFilter("atw")
  #goFilter("f3d",lim=(2.1,1.4,7.9,1.85))
  #goFilter("f3d",lim=(2.3,0.98,7.7,1.42))
  #goFilterQC("f3d",lim=(2.1,1.4,7.9,1.85))
  #goFilterImpulses()
  #goFilterRandom()
  #goFilterWithGpow()

def goRandomBlocks():
  n1 = 801
  x = makeBlocks(n1)
  plotB(x,x,-12,12,"rbx")
  return
  y = add(x,makeNoiseForBlocks(3141,8.0,3.0,n1))
  plotB(x,y,-12,12,"rbxy")
  sigmaS = 20.0
  yqqd = qqd(y)
  z = zerofloat(n1)
  for scale in [0.01,0.5,1.0,1.5,10,100.0]:
    sigmaX = scale*yqqd
    print "sigmaS =",sigmaS," sigmaX =",sigmaX
    bf = BilateralFilter(sigmaS,sigmaX)
    bf.setType(BilateralFilter.Type.TUKEY)
    bf.apply(y,z)
    png = "rbxz"+str(int(scale*10+0.5))
    if len(png)==1: png = "0"+png
    plotB(x,z,-12,12,png=png)
def plotB(x,y,ymin=0.0,ymax=0.0,png=None):
  sp = SimplePlot()
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setBackground(backgroundColor)
  #sp.setHLabel("Sample index")
  #sp.setVLabel("Amplitude")
  sp.setSize(720,500)
  pv = sp.addPoints(x) 
  pv.setLineWidth(3)
  pv.setLineStyle(PointsView.Line.DOT)
  pv = sp.addPoints(y) 
  pv.setLineWidth(3)
  #pv.setLineColor(Color.RED)
  if ymin<ymax:
    sp.setVLimits(ymin,ymax)
  if png and pngDir:
    sp.paintToPng(720,3.3,pngDir+png+".png")

def goRangeFunctions():
  xmin,xmax,sigma = -3.5,3.5,1.0
  nx = 351
  dx = (xmax-xmin)/(nx-1)
  fx = xmin
  sx = Sampling(nx,dx,fx)
  yg = BilateralFilter.sampleRangeFunction(gauss,sigma,sx)
  yt = BilateralFilter.sampleRangeFunction(tukey,sigma,sx)
  #yh = BilateralFilter.sampleRangeFunction(huber,sigma,sx)
  #x = rampfloat(fx,dx,nx)
  #yg = mul(x,yg)
  #yh = mul(x,yh)
  #yt = mul(x,yt)
  sp = SimplePlot()
  sp.setBackground(backgroundColor)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(700,500)
  sp.setHLabel(" ")
  sp.setVLabel(" ")
  solid = PointsView.Line.SOLID
  dash = PointsView.Line.DASH
  dot = PointsView.Line.DOT
  #pv = sp.addPoints(sx,yh); pv.setLineStyle(dot); pv.setLineWidth(4)
  pv = sp.addPoints(sx,yg); pv.setLineStyle(dash); pv.setLineWidth(4)
  pv = sp.addPoints(sx,yt); pv.setLineStyle(solid); pv.setLineWidth(4)
  if pngDir:
    sp.paintToPng(720,3.3,pngDir+"blftg.png")

def goApproximation():
  factor = 0.67448 # 3rd quartile of standard normal distribution (sigma = 1)
  sigma = 1.401114
  pmin,pmax = -9.0,9.0
  np,dp,fp = 301,0.06,pmin
  nk,dk,fk =  13,1.5,pmin
  sp = Sampling(np,dp,fp)
  sk = Sampling(nk,dk,fk)
  type = BilateralFilter.Type.GAUSS
  rf = BilateralFilter.sampleRangeFunction(type,sigma,sp)
  rg = BilateralFilter.sampleRangeFunction(type,sigma/sqrt(2.0),sp)
  method = CubicInterpolator.Method.MONOTONIC
  rfi = CubicInterpolator(method,np,rampfloat(fp,dp,np),rf)
  rgi = CubicInterpolator(method,np,rampfloat(fp,dp,np),rg)
  rr = zerofloat(np,np)
  ra = zerofloat(np,np)
  r0 = zerofloat(np,np)
  rh = zerofloat(np,np)
  for i in range(np):
    pi = sp.getValue(i)
    for j in range(np):
      pj = sp.getValue(j)
      rr[i][j] = rfi.interpolate(pi-pj)
      for k in range(nk):
        pk = sk.getValue(k)
        rjk = rgi.interpolate(pj-pk)
        rik = rgi.interpolate(pi-pk)
        #hik = hat(dk,pi-pk)
        #ra[i][j] += hik*rjk
        ra[i][j] += rik*rjk
        if k==nk/2:
          r0[i][j] = rjk
          rh[i][j] = rik*rjk
          #rh[i][j] = hik*rjk
  rr = div(rr,max(rr))
  r0 = div(r0,max(r0))
  rh = div(rh,max(rh))
  ra = div(ra,max(ra))
  plotR2(sp,rr,"blfrr")
  plotR2(sp,r0,"blfr0")
  plotR2(sp,rh,"blfrh")
  plotR2(sp,ra,"blfra")
def plotR2(sp,rp,what):
  plot = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  plot.setBackground(backgroundColor)
  plot.addColorBar()
  plot.setFontSizeForSlide(1.0,0.9)
  plot.setSize(740,630)
  plot.setHLabel(" ")
  plot.setVLabel(" ")
  plot.setHInterval(2)
  plot.setVInterval(2)
  pv = plot.addPixels(sp,sp,rp)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(0,1)
  if pngDir:
    plot.paintToPng(720,3.3,pngDir+what+".png")
def hat(dk,dp):
  if dp<-dk: 
    return 0.0
  elif dp<0.0:
    return 1.0+dp/dk
  elif dp<=dk:
    return 1.0-dp/dk
  else: 
    return 0.0

def goFilterQC(what,lim=None):
  x,s1,s2,clip = getImage(what)
  if lim: z = "z"+str(int(lim[1]))
  else: z = "" 
  sigmaR = qqd(x)
  sigmaS = 16.0
  print "sigmaS =",sigmaS," sigmaR =",sigmaR
  doFilterQC(what+"_blfq"+z,sigmaS,sigmaR,True,x,s1,s2,clip,lim)

def doFilterQC(what,sigmaS,sigmaR,useT,x,s1,s2,clip,limits=None):
  if useT:
    t = diffusionTensors(2.0,x)
  else:
    t = None
  sx,sn,sd,tn,td = bilateralFilterQC(sigmaS,sigmaR,t,x)
  nx,dx,fx, = sx.count,sx.delta,sx.first
  print "x sampling: count =",nx," delta =",dx," first =",fx
  for kx in range(nx):
    plot(what,sn[kx],s1,s2,clip=clip,limits=limits)
    plot(what,sd[kx],s1,s2,clip=clip,limits=limits)
    plot(what,tn[kx],s1,s2,clip=clip,limits=limits)
    plot(what,td[kx],s1,s2,clip=clip,limits=limits)

def goCtRock():
  x,s1,s2,clip = getImageCtr()
  sigmaR = 0.5*qqd(x)
  sigmaS = 8.0
  y = bilateralFilter(sigmaS,sigmaR,None,x)
  z = bilateralFilter(sigmaS,100.0*sigmaR,None,x)
  xy = sub(x,y)
  n1,n2 = s1.count,s2.count
  xt = threshold(0.57,0.78,x)
  yt = threshold(0.57,0.78,y)
  plot("ctr_blfx",x,s1,s2,cmin=0.35,cmax=0.99)
  plot("ctr_blfy",y,s1,s2,cmin=0.35,cmax=0.99)
  plot("ctr_blfxy",xy,s1,s2,cmin=-0.1,cmax=0.1)
  plot("ctr_blfz",z,s1,s2,cmin=0.35,cmax=0.99)
  plot("ctr_blfxt",xt,s1,s2,cmin=0.35,cmax=0.99)
  plot("ctr_blfyt",yt,s1,s2,cmin=0.35,cmax=0.99)
  hx = Histogram(flatten(x))
  vmin,vmax,nbin = hx.minValue,hx.maxValue,hx.binCount
  nbin *= 2
  hx = Histogram(flatten(x),vmin,vmax,nbin)
  hy = Histogram(flatten(y),vmin,vmax,nbin)
  plotHistograms("ctr_blfhx",hx)
  plotHistograms("ctr_blfhxy",hx,hy)

def goFilter(what,lim=None):
  x,s1,s2,clip = getImage(what)
  if lim: z = "z"+str(int(lim[1]))
  else: z = "" 
  plot(what+"_input"+z,x,s1,s2,clip=clip,limits=lim)
  sigmaR = qqd(x)
  sigmaS = 16.0
  print "sigmaS =",sigmaS," sigmaR =",sigmaR
  #doFilter(what+"_blf"+z,sigmaS,sigmaR,False,False,x,s1,s2,clip,lim)
  doFilter(what+"_blft"+z,sigmaS,sigmaR,True,False,x,s1,s2,clip,lim)
  #doFilter(what+"_blftc"+z,sigmaS,sigmaR,True,True,x,s1,s2,clip,lim)
  #doFilter(what+"_lsf"+z,sigmaS,100.0*sigmaR,False,False,x,s1,s2,clip,lim)
  #doFilter(what+"_lsft"+z,sigmaS,100.0*sigmaR,True,False,x,s1,s2,clip,lim)
  doFilter(what+"_lsftc"+z,sigmaS,100.0*sigmaR,True,True,x,s1,s2,clip,lim)
  
def doFilter(what,sigmaS,sigmaR,useT,useC,x,s1,s2,clip,limits=None):
  if useT:
    t = diffusionTensors(2.0,x)
    if useC:
      c = coherence(2.0,t,x)
      s = pow(c,1.0/8.0)
      plot(what+"s",s,s1,s2,cmin=0,cmax=1,limits=limits)
      plot(what+"c",c,s1,s2,cmin=0,cmax=1,limits=limits)
      t.scale(mul(c,c))
    plot(what+"e",x,s1,s2,clip=clip,tensors=t,limits=limits)
  else:
    t = None
  y = x
  y = bilateralFilter(sigmaS,sigmaR,t,y)
  print "y min =",min(y)," max =",max(y)
  plot(what,y,s1,s2,clip=clip,limits=limits)
  plot(what+"d",sub(x,y),s1,s2,clip=clip*0.251,limits=limits)

def goFilterImpulses():
  xa,s1,s2,clip = getImage()
  xb = makeImpulses(12,len(xa[0]),len(xa))
  t = diffusionTensors(2.0,xa)
  plot(xa,s1,s2,1.3)
  #plot(xb,s1,s2,0.9)
  for sigmaR in [0.1,1.0,100.0]:
  #for sigmaR in [100.0]:
    y = like(xa)
    bf = BilateralFilter(30.0,sigmaR)
    bf.setType(BilateralFilter.Type.TUKEY)
    bf.applyAB(t,xa,xb,y)
    y = smoothS(y)
    #plot(y,s1,s2,0.5*max(y))
    plot(y,s1,s2,0.1)

def goFilterRandom():
  xa,s1,s2,clip = getImage()
  plot(xa,s1,s2,clip)
  n1,n2 = len(xa[0]),len(xa)
  xb = makeRandom(n1,n2)
  t = diffusionTensors(2.0,xa)
  #plot(xa,s1,s2,1.3)
  #plot(xb,s1,s2,0.5)
  xqqd = qqd(x)
  for sigmaR in [xqqd,100*xqqd]:
    y = like(xa)
    bf = BilateralFilter(30.0,sigmaR)
    bf.setType(BilateralFilter.Type.TUKEY)
    bf.applyAB(t,xa,xb,y)
    plot(y,s1,s2,0.1)
    bf.apply(t,xa,y)
    plot(y,s1,s2,clip)

def goImage(what):
  x,s1,s2,clip = getImage(what)
  plot(what,x,s1,s2,clip=clip)

def structureTensors(sigma,x):
  n1,n2 = len(x[0]),len(x)
  if n1==500 and n2==500:
    lof = LocalOrientFilter(2.0*sigma)
    lof.setGradientSmoothing(1.0)
  else:
    lof = LocalOrientFilter(2.0*sigma,sigma)
    lof.setGradientSmoothing(sigma)
  t = lof.applyForTensors(x)
  return t

def diffusionTensors(sigma,x):
  t = structureTensors(sigma,x) # structure tensors
  t.invertStructure(0.0,4.0) # inverted with ev = 1.0, eu = small
  return t

def bilateralFilter(sigmaS,sigmaX,t,x):
  y = like(x)
  bf = BilateralFilter(sigmaS,sigmaX)
  if t:
    bf.apply(t,x,y)
  else:
    bf.apply(x,y)
  return y

def bilateralFilterQC(sigmaS,sigmaX,t,x):
  y = like(x)
  bf = BilateralFilter(sigmaS,sigmaX)
  if t:
    qc = bf.applyQC(t,x,y)
  else:
    qc = bf.applyQC(x,y)
  return qc.sx,qc.sn,qc.sd,qc.tn,qc.td

def smoothS(x):
  y = like(x)
  lsf = LocalSmoothingFilter()
  lsf.applySmoothS(x,y)
  return y

def smooth(sigma,t,x):
  z = copy(x)
  if t==None:
    rgf = RecursiveGaussianFilter(sigma)
    rgf.apply00(x,z)
  else:
    lsf = LocalSmoothingFilter(0.001,int(10*sigma))
    y = copy(x)
    #lsf.applySmoothS(y,y)
    #lsf.applySmoothL(kmax,y,y)
    lsf.apply(t,0.5*sigma*sigma,y,z)
  return z

def coherence(sigma,t,x):
  n1,n2 = len(x[0]),len(x)
  s = semblance(sigma,t,x) # structure-oriented semblance s
  if n1!=500 or n2!=500:
    s = pow(s,8.0) # make smaller sembances in [0,1] much smaller
  return s

def semblance(sigma,t,x):
  n1,n2 = len(x[0]),len(x)
  if n1==500 and n2==500:
    lsf = LocalSemblanceFilter(2*int(sigma),2*int(sigma))
    return lsf.semblance(LocalSemblanceFilter.Direction2.UV,t,x)
  else:
    lsf = LocalSemblanceFilter(int(sigma),4*int(sigma))
    return lsf.semblance(LocalSemblanceFilter.Direction2.V,t,x)

def threshold(xa,xb,x):
  n1,n2 = len(x[0]),len(x)
  y = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if x[i2][i1]<xa:
        y[i2][i1] = 0.0
      elif x[i2][i1]<xb:
        y[i2][i1] = 0.5
      else:
        y[i2][i1] = 1.0
  return y

def qqd(x):
  return 0.5*(Quantiler.estimate(0.75,x)-Quantiler.estimate(0.25,x))

def like(x):
  n1,n2 = len(x[0]),len(x)
  return zerofloat(n1,n2)

##############################################################################
# images

def getImage(what):
  if what[:3]=="ctr": return getImageCtr()
  if what[:3]=="f3d": return getImageF3d()
  if what[:3]=="tpd": return getImageTpd()
  if what[:3]=="atw": return getImageAtw()
  if what[:3]=="syn": return getImageSyn()

def getImageCtr():
  n1,n2 = 161,161
  d1,d2 = 0.0125,0.0125
  f1,f2 = 0.0000,0.0000
  fileName = "/data/rock/rock500.dat"
  x = readImage(fileName,n1,n2)
  print "rock min =",min(x)," max =",max(x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,1.0

def getImageF3d():
  n1,n2 = 462,951
  d1,d2 = 0.004,0.025
  f1,f2 = 0.004,0.000
  fileName = "/data/seis/f3d/f3d75.dat"
  x = readImage(fileName,n1,n2)
  subset = True
  if subset:
    j1,j2 = 240,0
    n1,n2 = n1-j1,440
    f1,f2 = f1+j1*d1,f2+j2*d2
    x = copy(n1,n2,j1,j2,x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,5.0

def getImageTpd():
  n1,n2 = 251,357
  d1,d2 = 0.004,0.025
  f1,f2 = 0.500,0.000
  fileName = "/data/seis/tp/csm/oldslices/tp73.dat"
  x = readImage(fileName,n1,n2)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,2.0

def getImageAtw():
  n1,n2 = 500,500
  d1,d2 = 0.02,0.02
  f1,f2 = 0.00,0.00
  fileName = "/data/seis/atw/atwj1s.dat"
  x = readImage(fileName,n1,n2)
  x = mul(0.001,x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return x,s1,s2,9.0

def getImageSyn():
  n1,n2 = 251,357
  d1,d2 = 1.0,1.0
  f1,f2 = 0.0,0.0
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  x = zerofloat(n1,n2)
  for i2 in range(n2):
    if i2<n2/2:
      x[i2] = sin(rampfloat(0.00,0.1,n1))
    else:
      x[i2] = sin(rampfloat(3.14,0.1,n1))
    x[i2] = mul(1.4,x[i2])
  return x,s1,s2,1.4

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def makeImpulses(ni,n1,n2):
  x = zerofloat(n1,n2)
  ns = max(n1/ni,n2/ni)
  m1 = (n1-1)/ns
  m2 = (n2-1)/ns
  j1 = (n1-1-(m1-1)*ns)/2
  j2 = (n2-1-(m2-1)*ns)/2
  for i2 in range(j2,n2,ns):
    for i1 in range(j1,n1,ns):
      x[i2][i1] = 1.0
  return x

def makeRandom(n1,n2):
  x = mul(2.0,sub(randfloat(n1,n2),0.5))
  return smooth(1.0,None,x)

def makeBlocks(n1):
  nb = 17
  db = 1.0
  pb = 8.0
  sb = 1.0
  xb = pb
  m1 = 1+n1/nb
  x = zerofloat(n1)
  for i1 in range(n1):
    if (i1+1)%m1==0:
      pb = pb-db
      sb = -sb
    x[i1] = sb*pb
  return x

def makeNoiseForBlocks(seed,scale,sigma,n1):
  r = Random(seed)
  x = mul(2.0*scale,sub(randfloat(r,n1),0.5))
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply0(x,x)
  return x

#############################################################################
# plot

backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent

def plot(what,f,s1,s2,tensors=None,limits=None,
         clip=None,cmin=None,cmax=None,
         cbar="",cbarInterval=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  #sp.setFontSizeForPrint(8.0,240) # print
  sp.setFontSizeForSlide(1.0,0.9) # slide
  if what[:3]=="atw":
    sp.setSize(980,872)
    sp.setHInterval(2.0)
    sp.setVInterval(2.0)
    sp.setHLabel("Inline (km)")
    sp.setVLabel("Crossline (km)")
  elif what[:3]=="ctr":
    sp.setSize(1040,890)
    sp.setHInterval(0.5)
    sp.setVInterval(0.5)
    sp.setHLabel("X (mm)")
    sp.setVLabel("Y (mm)")
  elif what[:3]=="f3d" or what[:3]=="tpd":
    sp.setSize(1040,700)
    if limits:
      sp.setHInterval(1.0)
      sp.setVInterval(0.1)
    else:
      sp.setHInterval(2.0)
      sp.setVInterval(0.2)
    sp.setHLabel("Inline (km)")
    sp.setVLabel("Time (s)")
  if limits:
    sp.setLimits(limits[0],limits[1],limits[2],limits[3])
  if clip:
    cmin = -clip
    cmax =  clip
  if cbar!=None:
    if len(cbar)>0:
      cbar = sp.addColorBar(cbar)
    else:
      cbar = sp.addColorBar()
    if cbarInterval:
      cbar.setInterval(cbarInterval)
    if not cbarInterval and cmin and cmax:
      if cmax-cmin>6.0:
        cbar.setInterval(2.0)
      elif cmax-cmin>=4.0:
        cbar.setInterval(1.0)
      elif cmax-cmin<1.0:
        cbar.setInterval(0.1)
      elif cmax-cmin<=2.0:
        cbar.setInterval(0.2)
    sp.plotPanel.setColorBarWidthMinimum(110)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin!=None and cmax!=None:
    pv.setClips(cmin,cmax)
  if tensors:
    tv = TensorsView(s1,s2,tensors)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(24)
    #tv.setScale(3)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir:
    sp.paintToPng(720,3.3,pngDir+what+".png")

def plotHistograms(what,hx,hy=None):
  sp = SimplePlot()
  sp.setBackground(backgroundColor)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(720,550)
  pv = sp.addPoints(hx.binSampling,hx.densities)
  pv.setLineWidth(3)
  pv.setLineColor(Color.RED)
  if hy:
    pv = sp.addPoints(hy.binSampling,hy.densities)
    pv.setLineWidth(5)
    pv.setLineColor(Color.BLUE)
  sp.setHLabel("Value")
  sp.setVLabel("Density")
  sp.setVLimits(0.0,0.119);
  if pngDir: 
    sp.paintToPng(720,3.3,pngDir+what+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
