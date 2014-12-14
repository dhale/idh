"""
Warping of well logs from Teapot Dome survey.
Author: Dave Hale, Colorado School of Mines
Version: 2013.12.28
"""

from java.awt import *
from java.io import *
from java.lang import *
from java.util import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mesh import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from tp import *
from warp import *

wlw = WellLogWarping()
curve = "v"
logs = None

def main(args):
  global logs
  logs = getLogs("d",curve)
  goShifts()
  #goWarping()
  #goErrors()
  #goResample()
  #goSort()
  #goMesh()

def goShifts():
  sz,fs = resample(logs,curve)
  if curve=="v":
    #fs = [fs[0],fs[4],fs[9],fs[14],fs[17],fs[20]] # deepest 6 velocity logs
    fs = [fs[0],fs[4],fs[9],fs[11],fs[14],fs[17],fs[20]] # 7 velocity logs
  elif curve=="d":
    fs = [fs[ 1],fs[ 2],fs[ 3],fs[ 4],fs[ 7],
          fs[11],fs[21],fs[22],fs[33],fs[35],
          fs[43],fs[48],fs[50],fs[56],fs[66],
          fs[81],fs[88],fs[163]] # deepest 18 density logs
  nk,nl = len(fs[0]),len(fs)
  sl = Sampling(nl,1.0,1.0)
  wlw.setPowError(0.25)
  wlw.setMaxShift(350)
  s = wlw.findShifts(fs)
  cs = [Color.BLACK,
    Color.RED,Color.GREEN,Color.BLUE,
    Color.CYAN,Color.MAGENTA,Color.YELLOW]
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setSize(1200,600)
  for i,si in enumerate(s):
    pv = sp.addPoints(si)
    pv.setLineColor(cs[i%len(cs)])
  """
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(650,550)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar("Shifts (m)")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl),s)
  pv.setClips(-250,250)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  """
  gs = wlw.applyShifts(fs,s)
  s = mul(1000*sz.delta,s) # convert shifts to m
  freplace = 2.0
  if curve=="d":
    freplace = 1.0
  fclips = (2.0,6.0)
  if curve=="d":
    fclips = (2.0,2.8)
  fs = wlw.replaceNulls(fs,freplace)
  gs = wlw.replaceNulls(gs,freplace)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(450,850)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar("Velocity (km/s)")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,sl,fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(450,850)
  sp.setVLabel("Relative geologic time")
  sp.setHLabel("Log index")
  sp.addColorBar("Velocity (km/s)")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,sl,gs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])

def goWarping():
  #pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)] # deepest 6 velocity logs
  pairs = [(4,0),(4,9),(4,14),(4,17),(4,20)] # log 4 is nearest to centroid
  sz,fs = resample(logs,curve)
  wlw.setPowError(0.25)
  #wlw.setPowError(2.00)
  wlw.setMaxShift(300)
  for pair in pairs:
    jf,jg = pair[0],pair[1]
    fi,gj = fs[jf],fs[jg]
    e = wlw.computeErrors(fi,gj)
    nl,nk = len(e[0]),len(e)
    d = wlw.accumulateErrors(e)
    kl = wlw.findWarping(d)
    fk,gk = wlw.applyWarping(kl,fi,gj)
    freplace = 2.0
    fi = wlw.replaceNulls(fi,freplace)
    gj = wlw.replaceNulls(gj,freplace)
    fk = wlw.replaceNulls(fk,freplace)
    gk = wlw.replaceNulls(gk,freplace)
    title = "("+str(jf)+","+str(jg)+")"
    sk = Sampling(2*sz.count-1,0.5*sz.delta,sz.first)
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel("Velocity (km/s)")
      pv = sp.addPoints(sz,fi)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sz,gj)
      pv.setLineColor(Color.RED)
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel("Velocity (km/s)")
      pv = sp.addPoints(sk,fk)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sk,gk)
      pv.setLineColor(Color.RED)

def goErrors():
  # 0 4 9 (11) 14 17 20 # 11 not one of the six deepest velocity logs
  #pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)]
  pairs = [(11,20)] # this pair shows errors shallow
  sz,f = resample(logs,curve)
  wlw.setPowError(0.25)
  wlw.setMaxShift(300)
  for pair in pairs:
    ia,ib = pair[0],pair[1]
    e = wlw.computeErrors(f[ia],f[ib])
    wlw.interpolateOddErrors(e)
    nl,nk = len(e[0]),len(e)
    lmax = (nl-1)/2
    lmin = -lmax
    sl = Sampling(nl,1,lmin)
    sk = Sampling(nk,1,0)
    title = "("+str(ia)+","+str(ib)+")"
    sp = SimplePlot()
    sp.setSize(750,500)
    sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    pv = sp.addPixels(sk,sl,transpose(e))
    pv.setColorModel(cjet)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    d = wlw.accumulateErrors(e)
    wlw.interpolateOddErrors(d)
    sp = SimplePlot()
    sp.setSize(750,500)
    sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    pv = sp.addPixels(sk,sl,transpose(d))
    pv.setColorModel(cjet)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,90)
    kw,lw = wlw.findWarping(d)
    kw = wlw.toFloat(kw)
    lw = wlw.toFloat(lw)
    pv = sp.addPoints(kw,lw)
    pv.setLineColor(Color.WHITE)

def goResample():
  nlog = len(logs)
  sz,f = resample(logs,curve)
  f = wlw.replaceNulls(f,-0.01)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar() #("Velocity (km/s)")
  pv = sp.addPixels(sz,Sampling(nlog),f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ajet)
  if curve=="v":
    pv.setClips(2.0,6.0)
  elif curve=="d":
    pv.setClips(2.0,2.8)
  elif curve=="p":
    pv.setClips(0.0,0.5)
  elif curve=="g":
    pv.setClips(50,250)

def goSort():
  x,y = [],[]
  for index,log in enumerate(logs):
    x.append(log.x2[0])
    y.append(log.x3[0])
  sp = SimplePlot()
  sp.setSize(700,380)
  sp.setHLabel("Crossline (km)")
  sp.setVLabel("Inline (km)")
  sp.setLimits(0.0,0.0,9.0,4.0)
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(6)
  if curve=="v":
    x = [x[0],x[4],x[9],x[14],x[17],x[20]]
    y = [y[0],y[4],y[9],y[14],y[17],y[20]]
    pv = sp.addPoints(x,y)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(10)

def goMesh():
  mesh = TriMesh()
  for i,log in enumerate(logs):
    node = TriMesh.Node(log.x2[0],log.x3[0])
    node.index = i
    mesh.addNode(node)
  sp = SimplePlot()
  sp.setSize(700,380)
  sp.setHLabel("Crossline (km)")
  sp.setVLabel("Inline (km)")
  sp.setLimits(0.0,0.0,9.0,4.0)
  tmv = TriMeshView(mesh)
  tmv.setLineColor(Color.BLACK)
  tmv.setMarkColor(Color.BLACK)
  tmv.setTriBoundsVisible(True)
  sp.add(tmv)

#############################################################################
# utilities

_tpDir = "/data/seis/tpd/"
_csmDir = _tpDir+"csm/"
_wellLogsDir = _csmDir+"welllogs/"
  
def getLogs(set,type):
  fileName = _wellLogsDir+"tpw"+set[0]+".dat"
  wdata = WellLog.Data.readBinary(fileName)
  logs = []
  logt = wdata.getLogsWith(type)
  for log in logt:
    logs.append(log)
  logs = sortLogs(logs)
  return logs

def sortLogs(logs):
  nlog = len(logs)
  xlog = zerodouble(nlog)
  ylog = zerodouble(nlog)
  for i,log in enumerate(logs):
    xlog[i] = log.x2[0]
    ylog[i] = log.x3[0]
  j = WellLogWarping.sortWells(xlog,ylog)
  logt = list(logs)
  for i,log in enumerate(logs):
    logt[i] = logs[j[i]]
  return logt

def resample(logs,curve):
  nlog = len(logs)
  zs = zerofloat(0,nlog)
  fs = zerofloat(0,nlog)
  for i,log in enumerate(logs):
    zs[i] = log.z
    fs[i] = log.getCurve(curve)
  zs = mul(zs,0.0003048) # ft to km
  sz = wlw.getDepthSampling(zs,fs)
  nz,dz,fz,lz = sz.count,sz.delta,sz.first,sz.last
  #print "resample before: nz =",nz," dz =",dz," fz =",fz
  dz = 0.001 # 1 m
  nz = 1+int((lz-fz)/dz)
  sz = Sampling(nz,dz,fz)
  #print "resample  after: nz =",nz," dz =",dz," fz =",fz
  fs = wlw.resampleLogs(sz,zs,fs)
  return sz,fs

def readLogSamples(set,type,smooth=0):
  """ 
  Reads log curves from the specified set that have the specified type.
  set: "s" for shallow, "d" for deep, or "a" for all
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  logs = getLogs(set,type)
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    if smooth: 
      log.smooth(smooth)
    samples = log.getSamples(type)
    if samples:
      f,x1,x2,x3 = samples
      fl.append(f)
      x1l.append(x1)
      x2l.append(x2)
      x3l.append(x3)
  return fl,x1l,x2l,x3l

def getWellIntersections(set,type,x1):
  fileName = _wellLogsDir+"tpw"+set[0]+".dat"
  wdata = WellLog.Data.readBinary(fileName)
  x2,x3 = wdata.getIntersections(type,x1)
  return x2,x3

def fgood(f):
  n = len(f)
  for i in range(n):
    if f[i]!=-999.2500:
      return i
def lgood(f):
  n = len(f)
  for i in range(n):
    if f[n-1-i]!=-999.2500:
      return n-1-i

#############################################################################
# graphics

cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

#############################################################################
# Run the function main on the Swing thread
import sys
from javax.swing import *
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
