import sys
from jarray import *

from java.awt import *
from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

from ldf import *
from joe import *

datadir = "/data/seis/joe/"

# 
nhead=3200 # number of bytes in EBCDIC header
nbhed=400 # number of bytes in binary header
nthed=240 # number of bytes in trace header
n1i=1051 # number of time samples (1st dimension)
d1i=0.004 # time sampling interval
f1i=0.000 # time of first sample
n2i=500 # number of traces in 2nd dimension
d2i=0.005
f2i=0.000
n3i=500 # number of traces in 3rd dimension
d3i=0.005
f3i=0.000

# sampling before windowing to get rid of wierd edges
#k1=750 # index of first sample in subvolume
#n1=251; d1=d1i; f1=f1i+k1*d1 # time sampling in subvolume
#n2=500; d2=d2i; f2=f2i
#n3=500; d3=d3i; f3=f3i

# sampling after windowing
#n1,n2,n3 = 251,400,500
#d1,d2,d3 = d1i,d2i,d3i
#f1,f2,f3 = f1i,f2i,f3i

# sampling after windowing
n1,n2,n3 = 200,200,200
d1,d2,d3 = 1.0,1.0,1.0
f1,f2,f3 = 0.0,0.0,0.0
#s1 = Sampling(n1,d1,f1)
#s2 = Sampling(n2,d2,f2)
#s3 = Sampling(n3,d3,f3)
s1 = Sampling(n1,0.004,3.0)
s2 = Sampling(n2,0.025,0.0)
s3 = Sampling(n3,0.025,0.0)

def main(args):
  #window()
  #x = readFile("x.dat",n1,n2,n3)
  #plot3d(x)
  #xplanes()
  #x = readFile("x.dat",n1,n2,n3)
  #ep = readFile("xep.dat",n1,n2,n3)
  #el = readFile("xel.dat",n1,n2,n3)
  #plotHistogram(el,ep)
  #ep = readFile("yep.dat",n1,n2,n3)
  #el = readFile("yel.dat",n1,n2,n3)
  #plotHistogram(el,ep)
  #plot3ds((x,el,ep))
  #ldf()
  #yplanes()
  #x = readFile("x.dat",n1,n2,n3)
  #y = readFile("y.dat",n1,n2,n3)
  #plot3ds((x,y))
  #ylines()
  #plot3ds((y,el))
  #llf()
  #emask()
  x = readFile("x.dat",n1,n2,n3)
  #exz = readFile("exz.dat",n1,n2,n3)
  #plot3ds((x,exz))
  eyz = readFile("eyz.dat",n1,n2,n3)
  plot3ds((x,eyz))
  #plotAll()
  #writeSlice23()
  return

def writeSlice23():
  x = readFloats23("x.dat",174)
  print "x min/max =",Array.min(x),Array.max(x)
  writeFile("x174.dat",x)

"""
  x clip = 0.0100
 xz clip = 0.0100
exz clip = 0.0030
  y clip = 0.0050
 yz clip = 0.0050
eyz clip = 0.0015
k1 = 174,177
k2 = 48
k3 = 65
"""
def plotAll():
  k1,k2,k3 = 174,48,65
  x = readFile("x.dat",n1,n2,n3)
  eyz = readFile("eyz.dat",n1,n2,n3)
  plot3ds((x,eyz))
  """
  plot3dPli(k1,k2,k3,"xep.dat","xel.dat",(1.0,1.0,1.0),"xe111.png")
  plot3dPli(k1,k2,k3,"yep.dat","yel.dat",(1.0,1.0,1.0),"ye111.png")
  plot3dPli(k1,k2,k3,"xep.dat","xel.dat",(0.9,0.3,0.6),"xe936.png")
  plot3dPli(k1,k2,k3,"yep.dat","yel.dat",(0.9,0.3,0.6),"ye936.png")
  plot3dPanels(k1,k2,k3,  "x.dat",0.0100,ColorMap.GRAY,  "x.png")
  plot3dPanels(k1,k2,k3, "xz.dat",0.0100,ColorMap.GRAY, "xz.png")
  plot3dPanels(k1,k2,k3,"exz.dat",0.0030,ColorMap.GRAY,"exz.png")
  plot3dPanels(k1,k2,k3,  "y.dat",0.0050,ColorMap.GRAY,  "y.png")
  plot3dPanels(k1,k2,k3, "yz.dat",0.0050,ColorMap.GRAY, "yz.png")
  plot3dPanels(k1,k2,k3,"eyz.dat",0.0015,ColorMap.GRAY,"eyz.png")
  """
  """
  xel = readFile("xel.dat",n1,n2,n3)
  xep = readFile("xep.dat",n1,n2,n3)
  plotHistogram(xel,xep,"xelp.png")
  yel = readFile("yel.dat",n1,n2,n3)
  yep = readFile("yep.dat",n1,n2,n3)
  plotHistogram(yel,yep,"yelp.png")
  """

def plotHistogram(el,ep,png=None):
  nl = 101
  np = 101
  h = Array.zerofloat(nl,np)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        il = int(el[i3][i2][i1]*(nl-1)+0.5)
        ip = int(ep[i3][i2][i1]*(np-1)+0.5)
        h[ip][il] += 1.0
  for ip in range(np):
    for il in range(nl):
      h[ip][il] = pow(h[ip][il],0.25)
  sp = SimplePlot()
  sp.setBackground(Color.WHITE)
  sp.setFontSize(30)
  sp.setSize(800,800)
  sp.setHLabel("linearity (%)")
  sp.setVLabel("planarity (%)")
  pv = sp.addPixels(Sampling(nl,1.0,0.0),Sampling(np,1.0,0.0),h)
  pv.setColorModel(ColorMap.JET)
  pv.setPercentiles(2,98)
  if png!=None:
    sp.paintToPng(200,6,png)

def ldf():
  x = readFile("x.dat",n1,n2,n3)
  iu = readFileShorts("xiu.dat",n1,n2,n3)
  sigma = 14
  small = 0.001
  niter = 100
  ldf = LocalDiffusionFilterCg(sigma,small,niter)
  y = Array.zerofloat(n1,n2,n3)
  ldf.applyPlanarKill(None,iu,x,y)
  writeFile("y.dat",y)

def ewarp(e):
  return Array.exp(Array.neg(Array.div(0.3,e)))

def llf():
  iw = readFileShorts("yiw.dat",n1,n2,n3)
  el = readFile("yel.dat",n1,n2,n3)
  il = LocalDiffusionFilter.encodeFractions(ewarp(el))
  sigma = 14
  small = 0.001
  niter = 100
  ldf = LocalDiffusionFilterCg(sigma,small,niter)
  x = readFile("x.dat",n1,n2,n3)
  #y = readFile("y.dat",n1,n2,n3)
  z = Array.zerofloat(n1,n2,n3)
  ldf.applyLinearPass(il,iw,x,z)
  #ldf.applyLinearPass(il,iw,y,z)
  writeFile("xz.dat",z)
  plot3ds((x,z))
  #writeFile("yz.dat",z)
  #plot3ds((y,z))

def emask():
  el = readFile("yel.dat",n1,n2,n3)
  el = ewarp(el)
  xz = readFile("xz.dat",n1,n2,n3)
  exz = Array.mul(el,xz) 
  writeFile("exz.dat",exz)
  #yz = readFile("yz.dat",n1,n2,n3)
  #eyz = Array.mul(el,yz) 
  #writeFile("eyz.dat",eyz)
  #x = readFile("x.dat",n1,n2,n3)
  #y = readFile("y.dat",n1,n2,n3)
  #plot3ds((x,exz))
  #plot3ds((y,eyz))

def window():
  x = readFile("win34.dat",251,500,500)
  #j1,j2,j3 =  0, 50, 0
  #m1,m2,m3 = n1,400,n3
  j1,j2,j3 =   0,150,250
  m1,m2,m3 = 200,200,200
  y = Array.copy(m1,m2,m3,j1,j2,j3,x)
  writeFile("x.dat",y)

def xplanes():
  x = readFile("x.dat",n1,n2,n3)
  lof = LocalOrientFilter(6.0)
  u1 = Array.zerofloat(n1,n2,n3)
  u2 = Array.zerofloat(n1,n2,n3)
  u3 = Array.zerofloat(n1,n2,n3)
  ep = Array.zerofloat(n1,n2,n3)
  el = Array.zerofloat(n1,n2,n3)
  lof.apply(x,
    None,None,
    u1,u2,u3,
    None,None,None,
    None,None,None,
    None,None,None,
    ep,el)
  iu = LocalDiffusionFilter.encodeUnitVectors(u1,u2,u3)
  writeFileShorts("xiu.dat",iu)
  writeFile("xu1.dat",u1)
  writeFile("xu2.dat",u2)
  writeFile("xu3.dat",u3)
  writeFile("xep.dat",ep)
  writeFile("xel.dat",el)

def yplanes():
  y = readFile("y.dat",n1,n2,n3)
  lof = LocalOrientFilter(6.0)
  u1 = Array.zerofloat(n1,n2,n3)
  u2 = Array.zerofloat(n1,n2,n3)
  u3 = Array.zerofloat(n1,n2,n3)
  ep = Array.zerofloat(n1,n2,n3)
  el = Array.zerofloat(n1,n2,n3)
  lof.apply(y,
    None,None,
    u1,u2,u3,
    None,None,None,
    None,None,None,
    None,None,None,
    ep,el)
  iu = LocalDiffusionFilter.encodeUnitVectors(u1,u2,u3)
  writeFileShorts("yiu.dat",iu)
  writeFile("yu1.dat",u1)
  writeFile("yu2.dat",u2)
  writeFile("yu3.dat",u3)
  writeFile("yep.dat",ep)
  writeFile("yel.dat",el)

def ylines():
  x = readFile("y.dat",n1,n2,n3)
  lof = LocalOrientFilter(6.0)
  w1 = Array.zerofloat(n1,n2,n3)
  w2 = Array.zerofloat(n1,n2,n3)
  w3 = Array.zerofloat(n1,n2,n3)
  el = Array.zerofloat(n1,n2,n3)
  lof.apply(x,
    None,None,
    None,None,None,
    None,None,None,
    w1,w2,w3,
    None,None,None,
    None,el)
  iw = LocalDiffusionFilter.encodeUnitVectors(w1,w2,w3)
  writeFileShorts("yiw.dat",iw)
  writeFile("yw1.dat",w1)
  writeFile("yw2.dat",w2)
  writeFile("yw3.dat",w3)
  writeFile("yel.dat",el)

def plot12(i3):
  ais = ArrayInputStream(datadir+"win34.dat")
  x = Array.zerofloat(n1,n2)
  ais.skipBytes(4*n1*n2*i3)
  ais.readFloats(x)
  ais.close()
  print "x min =",Array.min(x)," max =",Array.max(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(x)
  #pv.setPercentiles(1.0,99.0)
  pv.setClips(-1.0e-6,1.0e-6)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)

def plot3ds(xs,k=None):
  world = World()
  for x in xs:
    print "x min =",Array.min(x)," max =",Array.max(x)
    n1 = len(x[0][0])
    n2 = len(x[0])
    n3 = len(x)
    s1 = Sampling(n1)
    s2 = Sampling(n2)
    s3 = Sampling(n3)
    ipg = ImagePanelGroup(s3,s2,s1,SimpleFloat3(x))
    ipg.setPercentiles(1,99)
    #ipg.setPercentiles(0,100)
    world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def plot3d(x,cmin=-0.01,cmax=0.01):
  print "x min =",Array.min(x)," max =",Array.max(x)
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  s3 = Sampling(n3,d3,f3)
  x3 = SimpleFloat3(x)
  ipg = ImagePanelGroup(s3,s2,s1,x3)
  ipg.setClips(cmin,cmax)
  world = World()
  world.addChild(ipg)
  frame = TestFrame(world)
  frame.setVisible(True)

def readFileShorts(file,n1,n2,n3):
  ais = ArrayInputStream(datadir+file)
  x = Array.zeroshort(n1,n2,n3)
  ais.readShorts(x)
  ais.close()
  return x

def writeFileShorts(file,x):
  aos = ArrayOutputStream(datadir+file)
  aos.writeShorts(x)
  aos.close()

def readFile(file,n1,n2,n3):
  ais = ArrayInputStream(datadir+file)
  x = Array.zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  return x

def writeFile(file,x):
  aos = ArrayOutputStream(datadir+file)
  aos.writeFloats(x)
  aos.close()

def readFloats12(file,i3):
  f = readFile(file,n1,n2,n3)
  return slice12(f,i3)

def readFloats13(file,i2):
  f = readFile(file,n1,n2,n3)
  return slice13(f,i2)

def readFloats23(file,i1):
  f = readFile(file,n1,n2,n3)
  return slice23(f,i1)

def slice12(f,i3):
  f12 = Array.copy(f[i3])
  return f12

def slice13(f,i2):
  f13 = Array.zerofloat(n1,n3)
  for i3 in range(n3):
    Array.copy(n1,f[i3][i2],f13[i3])
  return f13

def slice23(f,i1):
  f23 = Array.zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      f23[i3][i2] = f[i3][i2][i1];
  return f23

def makeFloat3(x1,x2,x3):
  x = zeros(3,SimpleFloat3)
  x[0] = SimpleFloat3(x1)
  x[1] = SimpleFloat3(x2)
  x[2] = SimpleFloat3(x3)
  return x

def plot3dPli(k1,k2,k3,epfile,elfile,clips=None,png=None):
  orientation = PlotPanelPixels3.Orientation.X1DOWN
  axesPlacement = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  ep = readFile(epfile,n1,n2,n3)
  el = readFile(elfile,n1,n2,n3)
  ei = Array.sub(1.0,Array.add(ep,el))
  e = makeFloat3(ep,el,ei)
  panel = PlotPanelPixels3(orientation,axesPlacement,s1,s2,s3,e)
  if clips==None:
    clips = (1,1,1)
  panel.setClips(0,0.00,clips[0])
  panel.setClips(1,0.00,clips[1])
  panel.setClips(2,0.00,clips[2])
  #panel.setPercentiles(0,2,98)
  #panel.setPercentiles(1,2,98)
  #panel.setPercentiles(2,2,98)
  print "ep clip min =",panel.getClipMin(0),"  max =",panel.getClipMax(0)
  print "el clip min =",panel.getClipMin(1),"  max =",panel.getClipMax(1)
  print "ei clip min =",panel.getClipMin(2),"  max =",panel.getClipMax(2)
  panel.setSlices(k1,k2,k3)
  panel.setHLabel(0,"inline (km)")
  panel.setHLabel(1,"crossline (km)")
  panel.setVLabel(0,"crossline (km)")
  panel.setVLabel(1,"time (s)")
  panel.setLineColor(Color.WHITE)
  mosaic = panel.getMosaic()
  mosaic.getTileAxisBottom(0).setInterval(2.0)
  mosaic.getTileAxisBottom(1).setInterval(2.0)
  mosaic.getTileAxisLeft(0).setInterval(2.0)
  mosaic.getTileAxisLeft(1).setInterval(0.2)
  return frame(panel,png)

def plot3dPanels(k1,k2,k3,file,clip,cmod=ColorMap.GRAY,png=None):
  orientation = PlotPanelPixels3.Orientation.X1DOWN
  axesPlacement = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  f = readFile(file,n1,n2,n3)
  panel = PlotPanelPixels3(orientation,axesPlacement,s1,s2,s3,f)
  panel.setClips(-clip,clip)
  panel.setSlices(k1,k2,k3)
  panel.setHLabel(0,"inline (km)")
  panel.setHLabel(1,"crossline (km)")
  panel.setVLabel(0,"crossline (km)")
  panel.setVLabel(1,"time (s)")
  if isGray(cmod):
    panel.setLineColor(Color.YELLOW)
  mosaic = panel.getMosaic()
  mosaic.getTileAxisBottom(0).setInterval(2.0)
  mosaic.getTileAxisBottom(1).setInterval(2.0)
  mosaic.getTileAxisLeft(0).setInterval(2.0)
  mosaic.getTileAxisLeft(1).setInterval(0.2)
  return frame(panel,png)

def xplot3dPanels(k1,k2,k3,file,clip,cmod=ColorMap.GRAY,png=None):
  f12 = readFloats12(file,k3)
  f13 = readFloats13(file,k2)
  f23 = readFloats23(file,k1)
  print "plot3dPanels: f12 min =",Array.min(f12),"  max =",Array.max(f12)
  print "plot3dPanels: f13 min =",Array.min(f13),"  max =",Array.max(f13)
  print "plot3dPanels: f23 min =",Array.min(f23),"  max =",Array.max(f23)
  panel = PlotPanel(2,2,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  panel.setHLabel(0,"inline (km)")
  panel.setHLabel(1,"crossline (km)")
  panel.setVLabel(0,"crossline (km)")
  panel.setVLabel(1,"time (s)")
  #panel.setColorBarWidthMinimum(colorBarWidthMinimum)
  mosaic = panel.getMosaic()
  mosaic.getTileAxisBottom(0).setInterval(2.0)
  mosaic.getTileAxisBottom(1).setInterval(2.0)
  mosaic.getTileAxisLeft(0).setInterval(2.0)
  mosaic.getTileAxisLeft(1).setInterval(0.2)
  #mosaic.setWidthElastic(0,100*n2/n3)
  #cb = panel.addColorBar()
  p12 = panel.addPixels(1,0,s1,s2,f12)
  p13 = panel.addPixels(1,1,s1,s3,f13)
  p23 = panel.addPixels(0,0,s2,s3,f23)
  p23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  p12.setClips(-clip,clip)
  p13.setClips(-clip,clip)
  p23.setClips(-clip,clip)
  p12.setColorModel(cmod)
  p13.setColorModel(cmod)
  p23.setColorModel(cmod)
  xa1,xa2,xa3 = s1.getValue(0),s2.getValue(0),s3.getValue(0)
  xk1,xk2,xk3 = s1.getValue(k1),s2.getValue(k2),s3.getValue(k3)
  xb1,xb2,xb3 = s1.getValue(n1-1),s2.getValue(n2-1),s3.getValue(n3-1)
  p12 = panel.addPoints(1,0,((xa1,xb1),(xk1,xk1)),((xk2,xk2),(xa2,xb2)))
  p13 = panel.addPoints(1,1,((xa1,xb1),(xk1,xk1)),((xk3,xk3),(xa3,xb3)))
  p23 = panel.addPoints(0,0,((xa2,xb2),(xk2,xk2)),((xk3,xk3),(xa3,xb3)))
  p23.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  if isGray(cmod):
    p12.setLineColor(Color.YELLOW)
    p13.setLineColor(Color.YELLOW)
    p23.setLineColor(Color.YELLOW)
  return frame(panel,png)

def isGray(icm):
  for p in range(256):
    r = icm.getRed(p)
    g = icm.getGreen(p)
    b = icm.getBlue(p)
    if r!=g or r!=b:
      return False
  return True

#colorBarWidthMinimum = 80
frameWidth=800
frameHeight=800
def frame(panel,png):
  frame = PlotFrame(panel)
  frame.setBackground(Color.WHITE)
  frame.setFontSize(30)
  frame.setSize(frameWidth,frameHeight)
  frame.setVisible(True)
  if png!=None:
    frame.paintToPng(200,6,png)
  return frame

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
