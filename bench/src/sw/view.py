import sys
from math import *
from java.awt import *
from java.lang import *
from javax.swing import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

True = 1
False = 0

n1 = 1501  
d1 = 0.004
f1 = 0.0
s1 = Sampling(n1,d1,f1)

n2 = 623  
d2 = 0.025
f2 = 0
s2 = Sampling(n2,d2,f2)

n3 = 367
d3 = 0.025
f3 = 0
s3 = Sampling(n3,d3,f3)

datadir = "/data/seis/sw/all/"
#datadir = "/datc/seis/sw/all/"
#datadir = "/datb/seis/sw/all/"

pngdir = None

##############################################################################
# Read/write

def readFloats3(file):
  f = Array.zerofloat(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readFloats(f)
  af.close()
  return f

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

def readAndScaleSlices(file,scale,i1,i2,i3):
  f = readFloats3(file)
  f12 = Array.mul(scale,slice12(f,i3))
  f13 = Array.mul(scale,slice12(f,i2))
  f23 = Array.mul(scale,slice12(f,i1))
  return f12,f13,f23

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
  frame.setFontSize(30)
  frame.setSize(frameWidth,frameHeight)
  frame.setVisible(True)
  if pngdir!=None and png!=None:
    frame.paintToPng(300,6,pngdir+png+".png")
  return frame

def plot3d(k1,k2,k3,file,scale,clip,cmod=ColorMap.GRAY,png=None):
  f12,f13,f23 = readAndScaleSlices(file,scale,k1,k2,k3)
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

def isGray(icm):
  for p in range(256):
    r = icm.getRed(p)
    g = icm.getGreen(p)
    b = icm.getBlue(p)
    if r!=g or r!=b:
      return False
  return True

def plot3dAll():
  #k1 = 201 # = (4.404-3.600)/0.004
  k1 = 1101 # = (4.404-0.000)/0.004
  k2 =  293 # = (7.325-0.000)/0.025
  k3 =  170 # = (4.250-0.000)/0.025
  gray = ColorMap.GRAY
  flag = ColorMap.RED_WHITE_BLUE
  #plot3d(k1,k2,k3,"s02.dat",0.001,5.0,gray,"sws02")
  #plot3d(k1,k2,k3,"s04.dat",0.001,5.0,gray,"sws04")
  #plot3d(k1,k2,k3,"w02.dat",0.001,0.5,gray,"sww02")
  #plot3d(k1,k2,k3,"w04.dat",0.001,0.5,gray,"sww04")
  #plot3d(k1,k2,k3,"u1s1.dat",1000*d1,4.5,flag,"swu1s1")
  #plot3d(k1,k2,k3,"u2s1.dat",1000*d2,6.5,flag,"swu2s1")
  #plot3d(k1,k2,k3,"u3s1.dat",1000*d3,6.5,flag,"swu3s1")
  #plot3d(k1,k2,k3,"u1s2.dat",1000*d1,4.5,flag,"swu1s2")
  #plot3d(k1,k2,k3,"u2s2.dat",1000*d2,6.5,flag,"swu2s2")
  #plot3d(k1,k2,k3,"u3s2.dat",1000*d3,6.5,flag,"swu3s2")
  plot3d(k1,k2,k3,"u1s3.dat",1000*d1,4.5,flag,"swu1s3")
  plot3d(k1,k2,k3,"u2s3.dat",1000*d2,6.5,flag,"swu2s3")
  plot3d(k1,k2,k3,"u3s3.dat",1000*d3,6.5,flag,"swu3s3")

##############################################################################
# 3-D View

def ipg(file,clip,cmod):
  x = readFloats3(file)
  ipg = ImagePanelGroup(s1,s2,s3,x)
  ipg.setClips(-clip,clip)
  ipg.setColorModel(cmod)
  return ipg

def view3dAll():
  gray = ColorMap.GRAY
  flag = ColorMap.RED_WHITE_BLUE
  world = World()
  #world.addChild(ipg("s02.dat",5000,gray))
  #world.addChild(ipg("u1s2.dat",0.50,flag))
  #world.addChild(ipg("u2s2.dat",0.25,flag))
  world.addChild(ipg("u3s2.dat",0.25,flag))
  #world.addChild(ipg("e2s2.dat",0.50,flag))
  world.addChild(ipg("e3s2.dat",0.25,flag))
  frame = TestFrame(world)
  frame.setVisible(True)

#############################################################################
def main(args):
  #plot3dAll()
  view3dAll()
  return

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
