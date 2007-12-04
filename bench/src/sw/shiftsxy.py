import sys
from math import *
from java.awt import *
from java.lang import *
from javax.swing import *
from sw import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

True = 1
False = 0

n1 = 301  
d1 = 0.004
f1 = 3.6
s1 = Sampling(n1,d1,f1)

n2 = 623  
d2 = 0.025
f2 = 0
s2 = Sampling(n2,d2,f2)

n3 = 367
d3 = 0.025
f3 = 0
s3 = Sampling(n3,d3,f3)

datadir = "/data/seis/sw/"

##############################################################################
# Read/write

def readFloats3(file):
  f = Array.zerofloat(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readFloats(f)
  af.close()
  return f

def writeFloats3(file,f):
  af = ArrayFile(datadir+file,"rw")
  af.writeFloats(f)
  af.close()
  return f

def readBytes3(file):
  b = Array.zerobyte(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readBytes(b)
  af.close()
  return b

def writeBytes3(file,lag):
  af = ArrayFile(datadir+file,"rw")
  af.writeBytes(lag)
  af.close()

def readFloats12(file,i3):
  f = readFloats3(file)
  return slice12(f,i3)

def readFloats13(file,i2):
  f = readFloats3(file)
  return slice13(f,i2)

def readFloats23(file,i1):
  f = readFloats3(file)
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

##############################################################################
# Plot
#colorBarWidthMinimum = 100
#frameWidth=1030
#frameHeight=620
colorBarWidthMinimum = 80
frameWidth=1160
frameHeight=768

def frame(panel,png):
  frame = PlotFrame(panel)
  frame.setBackground(Color.WHITE)
  frame.setFont(Font("Arial",Font.PLAIN,30))
  #frame.setFontSize(30)
  frame.setSize(frameWidth,frameHeight)
  frame.setVisible(True)
  if png!=None:
    frame.paintToPng(300,6,png)
  return frame

def plot3d(k1,k2,k3,file,scale,clip,cmod=ColorMap.GRAY,png=None):
  f12 = Array.mul(scale,readFloats12(file,k3))
  f13 = Array.mul(scale,readFloats13(file,k2))
  f23 = Array.mul(scale,readFloats23(file,k1))
  print "plot3d: f12 min =",Array.min(f12),"  max =",Array.max(f12)
  print "plot3d: f13 min =",Array.min(f13),"  max =",Array.max(f13)
  print "plot3d: f23 min =",Array.min(f23),"  max =",Array.max(f23)
  panel = PlotPanel(2,2,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  mosaic = panel.getMosaic()
  mosaic.setWidthElastic(0,100*n2/n3)
  mosaic.getTileAxisBottom(0).setInterval(2.0)
  mosaic.getTileAxisBottom(1).setInterval(2.0)
  mosaic.getTileAxisLeft(0).setInterval(2.0)
  mosaic.getTileAxisLeft(1).setInterval(0.5)
  panel.setColorBarWidthMinimum(colorBarWidthMinimum)
  panel.getMosaic().setWidthElastic(0,100*n2/n3)
  panel.setHLabel(0,"inline (km)")
  panel.setHLabel(1,"crossline (km)")
  panel.setVLabel(0,"crossline (km)")
  panel.setVLabel(1,"time (s)")
  cb = panel.addColorBar()
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

def window(k1,k2,f):
  n1 = len(f[0])
  n2 = len(f)
  g = Array.zerofloat(n1,n2)
  s = -0.5/(12*12)
  for i2 in range(0,n2):
    x2 = i2-k2
    for i1 in range(0,n1):
      x1 = i1-k1
      g[i2][i1] = f[i2][i1]*exp(s*(x1*x1+x2*x2))
  return g

def wplot3d(k1,k2,k3,file,scale,clip,cmod=ColorMap.GRAY,png=None):
  f12 = Array.mul(scale,readFloats12(file,k3))
  f13 = Array.mul(scale,readFloats13(file,k2))
  f23 = Array.mul(scale,readFloats23(file,k1))
  print "plot3d: f12 min =",Array.min(f12),"  max =",Array.max(f12)
  print "plot3d: f13 min =",Array.min(f13),"  max =",Array.max(f13)
  print "plot3d: f23 min =",Array.min(f23),"  max =",Array.max(f23)
  n1 = len(f12[0])
  n2 = len(f12)
  n3 = len(f13)
  n1 = 1+n1/4
  n2 = 1+n2/4
  n3 = 1+n3/4
  j1 = k1-n1/2
  j2 = k2-n2/2
  j3 = k3-n3/2
  k1 = k1-j1
  k2 = k2-j2
  k3 = k3-j3
  s1 = Sampling(n1,d1,f1+j1*d1)
  s2 = Sampling(n2,d2,f2+j2*d2)
  s3 = Sampling(n3,d3,f3+j3*d3)
  f12 = Array.copy(n1,n2,j1,j2,f12)
  f13 = Array.copy(n1,n3,j1,j3,f13)
  f23 = Array.copy(n2,n3,j2,j3,f23)
  f12 = window(k1,k2,f12)
  f13 = window(k1,k3,f13)
  f23 = window(k2,k3,f23)
  panel = PlotPanel(2,2,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  mosaic = panel.getMosaic()
  mosaic.setWidthElastic(0,100*n2/n3)
  mosaic.getTileAxisBottom(0).setInterval(1.0)
  mosaic.getTileAxisBottom(1).setInterval(1.0)
  mosaic.getTileAxisLeft(0).setInterval(1.0)
  mosaic.getTileAxisLeft(1).setInterval(0.1)
  panel.setColorBarWidthMinimum(colorBarWidthMinimum)
  panel.getMosaic().setWidthElastic(0,100*n2/n3)
  panel.setHLabel(0,"inline (km)")
  panel.setHLabel(1,"crossline (km)")
  panel.setVLabel(0,"crossline (km)")
  panel.setVLabel(1,"time (s)")
  cb = panel.addColorBar()
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

def plot3dAll():
  k1 = 201 # = (4.404-3.600)/0.004
  k2 = 293 # = (7.325-0.000)/0.025
  k3 = 170 # = (4.250-0.000)/0.025
  gray = ColorMap.GRAY
  flag = ColorMap.RED_WHITE_BLUE
  #wplot3d(k1,k2,k3,"sw02a.dat",0.002,5.0,gray,"sws02z.png")
  #wplot3d(k1,k2,k3,"sw04a.dat",0.002,5.0,gray,"sws04z.png")
  #wplot3d(k1,k2,k3,"w02.dat",0.002,0.5,gray,"sww02w.png")
  #wplot3d(k1,k2,k3,"w04.dat",0.002,0.5,gray,"sww04w.png")
  #plot3d(k1,k2,k3,"sw02a.dat",0.001,5.0,gray,"sws02.png")
  #plot3d(k1,k2,k3,"sw04a.dat",0.001,5.0,gray,"sws04.png")
  #plot3d(k1,k2,k3,"w02.dat",0.001,0.5,gray,"sww02.png")
  #plot3d(k1,k2,k3,"w04.dat",0.001,0.5,gray,"sww04.png")
  #plot3d(k1,k2,k3,"u1s0.dat",1000*d1,4.5,flag,"swu1s0.png")
  #plot3d(k1,k2,k3,"u2s0.dat",1000*d2,6.5,flag,"swu2s0.png")
  #plot3d(k1,k2,k3,"u3s0.dat",1000*d3,6.5,flag,"swu3s0.png")
  #plot3d(k1,k2,k3,"u1s1.dat",1000*d1,4.5,flag,"swu1s1.png")
  plot3d(k1,k2,k3,"u2s1.dat",1000*d2,6.5,flag,"swu2s1.png")
  plot3d(k1,k2,k3,"u3s1.dat",1000*d3,6.5,flag,"swu3s1.png")
  #plot3d(k1,k2,k3,"e2s1.dat",1000*d2,6.5,flag,"swe2s1.png")
  #plot3d(k1,k2,k3,"e3s1.dat",1000*d3,6.5,flag,"swe3s1.png")

#############################################################################
# Research

def computeDeltaXY():
  dx,dy,dt = d3,d2,d1
  sx,sy,st = s3,s2,s1
  shifts = Shifts(sx,sy,st)
  r = 5.0
  v0 = Array.rampfloat(4.0,0.5*dt,n1)
  deltat = readFloats3("u1s1.dat")
  Array.mul(dt,deltat,deltat) # shifts are in samples
  deltax = shifts.getDeltaX(r,v0,deltat)
  deltay = shifts.getDeltaY(r,v0,deltat)
  Array.mul(1/dx,deltax,deltax) # shifts are in samples
  Array.mul(1/dy,deltay,deltay) # shifts are in samples
  writeFloats3("e2s1.dat",deltay)
  writeFloats3("e3s1.dat",deltax)

def main(args):
  #computeDeltaXY()
  plot3dAll()
  return

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
