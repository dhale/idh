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

#n1,d1,f1 = 1501,0.004,0.0
#n1,d1,f1 = 301,0.004,3.6
#n1,d1,f1 = 601,0.004,0.0
n1,d1,f1 = 601,0.004,2.4
n2,d2,f2 = 623,0.025,0.0
n3,d3,f3 = 367,0.025,0.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)

#datadir = "/data/seis/sw/sub0/"
datadir = "/data/seis/sw/sub24/"
#datadir = "/datc/seis/sw/all/"
#datadir = "/datb/seis/sw/all/"

#pngdir = None
pngdir = "png/"

##############################################################################

def readAndScaleVolume(file,scale):
  a = Array.zerofloat(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readFloats(a)
  af.close()
  Array.mul(scale,a,a)
  return a

def readAndScaleSlices(file,scale,i1,i2,i3):
  f12 = Array.zerofloat(n1,n2)
  f13 = Array.zerofloat(n1,n3)
  f23 = Array.zerofloat(n2,n3)
  a = Array.zerofloat(n1)
  af = ArrayFile(datadir+file,"r")
  for j3 in range(n3):
    for j2 in range(n2):
      af.readFloats(a)
      f23[j3][j2] = scale*a[i1]
      if j3==i3:
        Array.mul(scale,a,f12[j2])
      if j2==i2:
        Array.mul(scale,a,f13[j3])
  af.close()
  return f12,f13,f23

##############################################################################

colorBarWidthMinimum = 80
frameWidth = 1160
frameHeight = 768

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

# Problem with this is that time and space samplings are different.
# PlotPanelPixels3 needs a parameter for aspect ratio.
def plot3dNew(k1,k2,k3,file,scale,clip,cmod=ColorMap.GRAY,png=None):
  f = readAndScaleVolume(file,scale)
  print "plot3d: f min =",Array.min(f),"  max =",Array.max(f)
  panel = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,f)
  panel.addColorBar()
  panel.setSlice23(k1)
  panel.setSlice13(k2)
  panel.setSlice12(k3)
  panel.setInterval1(2.0)
  panel.setInterval2(2.0)
  panel.setInterval3(2.0)
  panel.setLabel1("time (s)")
  panel.setLabel2("inline (km)")
  panel.setLabel3("crossline (km)")
  panel.setColorBarWidthMinimum(colorBarWidthMinimum)
  panel.setClips(-clip,clip)
  panel.setColorModel(cmod)
  if isGray(cmod):
    panel.setLineColor(Color.YELLOW)
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
  #k1 = 1101 # = (4.404-0.000)/0.004
  #k1 = 201 # = (4.404-3.600)/0.004
  #k1 = 501 # = (4.404-2.400)/0.004
  k1 = 501 # = (2.404-0.000)/0.004
  k2 = 293 # = (7.325-0.000)/0.025
  k3 = 170 # = (4.250-0.000)/0.025
  gray = ColorMap.GRAY
  flag = ColorMap.RED_WHITE_BLUE
  #plot3d(k1,k2,k3,"s02.dat",0.001,5.0,gray,"sws02")
  #plot3d(k1,k2,k3,"s04.dat",0.001,5.0,gray,"sws04")
  #plot3d(k1,k2,k3,"w02.dat",0.001,0.5,gray,"sww02")
  #plot3d(k1,k2,k3,"w04.dat",0.001,0.5,gray,"sww04")
  #plot3d(k1,k2,k3,"u1s1.dat",1000*d1,3.25,flag,"swu1s1")
  #plot3d(k1,k2,k3,"u2s1.dat",1000*d2,6.50,flag,"swu2s1")
  #plot3d(k1,k2,k3,"u3s1.dat",1000*d3,6.50,flag,"swu3s1")
  #plot3d(k1,k2,k3,"u1s2.dat",1000*d1,3.25,flag,"swu1s2")
  #plot3d(k1,k2,k3,"u2s2.dat",1000*d2,6.50,flag,"swu2s2")
  #plot3d(k1,k2,k3,"u3s2.dat",1000*d3,6.50,flag,"swu3s2")
  #plot3d(k1,k2,k3,"u1s3.dat",1000*d1,3.25,flag,"swu1s3")
  #plot3d(k1,k2,k3,"u2s3.dat",1000*d2,6.50,flag,"swu2s3")
  #plot3d(k1,k2,k3,"u3s3.dat",1000*d3,6.50,flag,"swu3s3")
  plot3d(k1,k2,k3,"e2s3.dat",1000*d2,6.50,flag,"swe2s3")
  plot3d(k1,k2,k3,"e3s3.dat",1000*d3,6.50,flag,"swe3s3")

#############################################################################
def main(args):
  plot3dAll()
  return

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
