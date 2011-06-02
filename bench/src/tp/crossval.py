"""
Cross-validation (leave-one-out)
"""
from tputils import *

setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
logSet = "deep" # use only deep logs
logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.6
#logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.0,3.0
#logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
#logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
vomit = [9,5,24,7,19,6] # velocity logs to leave out (omit)
vlabel = {9:"A",5:"B",24:"C",7:"D",19:"E",6:"F"}
smin,smax = -5.5,5.5
sfile = "tpsz" # seismic image
esfile = "tpets" # eigen-tensors scaled by semblances
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,
#smooth = 50 # half-width of smoothing filter for logs
smooth = 0 # no smoothing

#pngDir = "png/"
pngDir = None

def main(args):
  #goCrossVal()
  goErrors()
  #goDisplay()

def goCrossVal():
  for omit in vomit:
    #crossVal(omit)
    writeLogs(omit)

def crossVal(omit):
  #gridWellLogs(omit)
  #gridBlendedP(omit)
  gridBlendedQ(omit)
  return

def goErrors():
  print "median absolute error"
  print "log","  all "," shal "," deep "
  for omit in vomit:
    errors(omit,"mda")
  print "root-mean-square error"
  print "log","  all "," shal "," deep "
  for omit in vomit:
    errors(omit,"rms")

def goDisplay():
  displayWellPoints(1.0)
  displayWellPoints(1.5)
  for omit in vomit:
    displayLogs(omit)

def errors(ilog,type="mda"):
  #fw,x1w,x2w,x3w = readLog(wellLog(ilog))
  fg,x1g,x2g,x3g = readLog(griddedLog(ilog))
  fs,x1s,x2s,x3s = smoothLog(fg),x1g,x2g,x3g
  #fs,x1s,x2s,x3s = fg,x1g,x2g,x3g
  #fp,x1p,x2p,x3p = readLog(nearestLog(ilog))
  fq,x1q,x2q,x3q = readLog(blendedLog(ilog))
  #fgs,fgd = splitShallowDeep(x1g,fg)
  fss,fsd = splitShallowDeep(x1s,fs)
  fqs,fqd = splitShallowDeep(x1q,fq)
  if type=="mda":
    efunc = mdae
  elif type=="rms":
    efunc = rmse
  aerr = efunc(fs, fq )
  serr = efunc(fss,fqs)
  derr = efunc(fsd,fqd)
  print "%3d %6.3f %6.3f %6.3f" % (ilog,aerr,serr,derr)

def splitShallowDeep(x1,f):
  n = len(x1)
  m = binarySearch(x1,1.3)
  if m<0: m = -m-1
  fs = copy(m,0,f)
  fd = copy(n-m,m,f)
  return fs,fd

def mdae(x,y):
  n = len(x)
  e = abs(sub(x,y))
  quickPartialSort(n/2,e)
  return e[n/2]

def rmse(x,y):
  e = sub(x,y)
  n = len(x)
  return sqrt(sum(mul(e,e))/n)

def writeLogs(ilog):
  fw,x1w,x2w,x3w = readWellLog(logSet,logType,ilog)
  print "log",ilog,": x1 =",x1w[0],"x2 =",x2w[0],"x3 =",x3w[0]
  p = readImage(nearestOmit(ilog))
  q = readImage(blendedOmit(ilog))
  wlg = WellLogGridder(s1,s2,s3,0.0)
  fg,x1g,x2g,x3g = wlg.getGriddedSamples(fw,x1w,x2w,x3w)
  fp = wlg.getGriddedValues(x1g,x2g,x3g,p)
  fq = wlg.getGriddedValues(x1g,x2g,x3g,q)
  writeLog(wellLog(ilog),fw,x1w,x2w,x3w)
  writeLog(griddedLog(ilog),fg,x1g,x2g,x3g)
  writeLog(nearestLog(ilog),fp,x1g,x2g,x3g)
  writeLog(blendedLog(ilog),fq,x1g,x2g,x3g)

def writeLog(fileName,f,x1,x2,x3):
  print "writeLog:",fileName," n =",len(f)
  aos = ArrayOutputStream(getSeismicDir()+fileName+".dat")
  aos.writeInt(len(f))
  aos.writeFloats(f)
  aos.writeFloats(x1)
  aos.writeFloats(x2)
  aos.writeFloats(x3)
  aos.close()

def readLog(fileName):
  ais = ArrayInputStream(getSeismicDir()+fileName+".dat")
  n = ais.readInt()
  f = zerofloat(n)
  x1 = zerofloat(n)
  x2 = zerofloat(n)
  x3 = zerofloat(n)
  ais.readFloats(f)
  ais.readFloats(x1)
  ais.readFloats(x2)
  ais.readFloats(x3)
  ais.close()
  return f,x1,x2,x3

def readWellLog(set,type,index):
  fl,x1l,x2l,x3l = readLogSamples(set,type,smooth)
  return fl[index],x1l[index],x2l[index],x3l[index]

def gridWellLogs(omit=-1):
  print "gridWellLogs:",logType,"without log",omit
  fnull = 0.0
  wlg = WellLogGridder(s1,s2,s3,fnull)
  fl,x1l,x2l,x3l = readLogSamples(logSet,logType,smooth)
  for i in range(len(fl)):
    f,x1,x2,x3 = fl[i],x1l[i],x2l[i],x3l[i]
    if i!=omit:
      wlg.insertWellLog(f,x1,x2,x3)
  g = wlg.getGriddedValues()
  writeImage(griddedOmit(omit),g)
  
def gridBlendedP(omit=-1):
  print "gridBlendedP:",logType,"without log",omit
  e = getEigenTensors()
  bi = BlendedGridder3(e)
  p = readImage(griddedOmit(omit))
  t = bi.gridNearest(0.0,p)
  writeImage(nearestOmit(omit),p)
  writeImage(mintimeOmit(omit),t)

def gridBlendedQ(omit=-1):
  print "gridBlendedQ:",logType,"without log",omit
  e = getEigenTensors()
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(nearestOmit(omit))
  t = readImage(mintimeOmit(omit))
  t = clip(0.0,50.0,t)
  q = copy(p)
  bg.gridBlended(t,p,q)
  writeImage(blendedOmit(omit),q)

def getEigenTensors():
  e = readTensors(esfile)
  return e

def smoothLog(f):
  #sigma = 2.0*s1.delta/0.0001524 # 0.0001524 km = 6 inches
  sigma = 2.5
  lpad = int(3.0*sigma)
  n = len(f)
  npad = lpad+n+lpad
  fpad = zerofloat(npad)
  gpad = zerofloat(npad)
  for i in range(lpad):
    fpad[i] = f[0]
  copy(n,0,f,lpad,fpad)
  for i in range(lpad+n,npad):
    fpad[i] = f[-1]
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply0(fpad,gpad)
  g = copy(n,lpad,gpad)
  return g 

def displayLogs(k):
  fw,x1w,x2w,x3w = readLog(wellLog(k))
  fg,x1g,x2g,x3g = readLog(griddedLog(k))
  fp,x1p,x2p,x3p = readLog(nearestLog(k))
  fq,x1q,x2q,x3q = readLog(blendedLog(k))
  fs,x1s,x2s,x3s = smoothLog(fg),x1g,x2g,x3g
  #fs,x1s,x2s,x3s = fg,x1g,x2g,x3g
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(450,900)
  sp.setFontSizeForSlide(0.5,1.0)
  #sp.setTitle("Log "+str(k))
  sp.plotPanel.setHInterval(1.0)
  sp.plotPanel.setVInterval(0.2)
  sp.setHLimits(2.800,6.200)
  sp.setVLimits(0.575,2.025)
  sp.setHLabel("Velocity (km/s)")
  sp.setVLabel("Depth (km)")
  pv = sp.addPoints(x1w,fw)
  pv.setLineColor(Color.GRAY)
  pv = sp.addPoints(x1s,fs)
  pv.setLineColor(Color.BLACK)
  pv.setLineWidth(4.0)
  #pv = sp.addPoints(x1p,fp)
  #pv.setLineColor(Color.RED)
  #pv.setLineWidth(4.0)
  pv = sp.addPoints(x1q,fq)
  pv.setLineColor(Color.RED)
  pv.setLineWidth(4.0)
  if pngDir:
    png = pngDir+"tpwsqv"+vlabel[k]+".png"
    sp.paintToPng(300,3,png)

def displayWellPoints(x1):
  x2,x3 = getWellIntersections(logSet,logType,x1)
  pv = PointsView(x2,x3)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(8.0)
  pp = PlotPanel(1,1,
                 PlotPanel.Orientation.X1RIGHT_X2UP,
                 PlotPanel.AxesPlacement.NONE)
  pp.setHLimits(s2.first,s2.last)
  pp.setVLimits(s3.first,s3.last)
  pp.addTiledView(pv)
  pf = PlotFrame(pp)
  pf.setSize(600,293)
  pf.setVisible(True)
  if pngDir:
    png = pngDir+"tpwp"+str(int(0.5+10*x1))+".png"
    pf.paintToPng(300,3,png)

def griddedOmit(omit):
  return "tpg"+omitSuffix(omit)
def nearestOmit(omit):
  return "tpp"+omitSuffix(omit)
def mintimeOmit(omit):
  return "tpt"+omitSuffix(omit)
def blendedOmit(omit):
  return "tpq"+omitSuffix(omit)
def wellLog(ilog):
  return "tpw"+ilogSuffix(ilog)
def griddedLog(ilog):
  return "tpg"+ilogSuffix(ilog)
def nearestLog(ilog):
  return "tpp"+ilogSuffix(ilog)
def mintimeLog(ilog):
  return "tpt"+ilogSuffix(ilog)
def blendedLog(ilog):
  return "tpq"+ilogSuffix(ilog)
def omitSuffix(omit): # suffix for image files with a log omitted
  if omit<10:
    return logType[0]+"o0"+str(omit)
  else:
    return logType[0]+"o"+str(omit)
def ilogSuffix(ilog): # suffix for log files with specified index
  if ilog<10:
    return logType[0]+"l0"+str(ilog)
  else:
    return logType[0]+"l"+str(ilog)

#############################################################################
#############################################################################
#############################################################################

def display3d(omit=-1):
  s = readImage(sfile)
  q = readImage(blendedOmit(omit))
  world = World()
  addImage2ToWorld(world,s,q)
  addLogsToWorld(world,logSet,logType)
  #addHorizonToWorld(world,"TensleepASand")
  makeFrame(world)

#############################################################################
run(main)
