"""
Converts 3D seismic images exported from OpendTect to BIG_ENDIAN floats.
In OpendTect, I used the "Simple file" export and specified subsets 
that include no null traces and no sampling information. In other words,
the exported files are assumed to contain only LITTLE_ENDIAN floats. See
the function setGlobals below to find the subsets that I specified in the
OpendTect -> Export -> Simple file dialog box.

More information about these data is available on the OpendTect website.

f3d.dat
Netherlands F3 Block, Offshore North Sea
The dataset provided by OpendTect includes multiple 3D images.
I used this script to convert the "original seismics" image.

pen.dat
Canada Penobscot 3D Survey, Offshore Nova Scotia
The 3D image provided by OpendTect is labeled "PSTM stack agc",
which implies prestack time migration and automatic gain control.

@author: Dave Hale, Colorado School of Mines
@version: 2010.11.05 
"""
from imports import *

#############################################################################
def main(args):
  for name in ["f3d"]: # ["f3d","pen"]:
    setGlobals(name)
    #convertEndian(odtfile,datfile)
    display3d(datfile)
    #displaySlice2(datfile)

def setGlobals(what):
  global n1,n2,n3,d1,d2,d3,f1,f2,f3
  global datDir,datbase,datfile,odtfile
  global clip,scale,scale1
  datDir = "/data/seis/odt/"
  clip = 5000.0
  if what=="f3d":
    n1,n2,n3 = 462,951,591 # time[4,1848], xline[100,690], iline[300,1250]
    d1,d2,d3 = 0.004,0.025,0.025 # s, km, km
    f1,f2,f3 = 0.004,0.000,0.000 # s, km, km
    datbase = "f3d"
    datfile = datDir+"f3d.dat" # flat file of BIG_ENDIAN floats
    odtfile = datDir+"f3dx.dat" # simple file exported from OpendTect
    scale,scale1 = 3.0,8.0
  if what=="pen":
    n1,n2,n3 = 1501,480,456 # time[0,6000] xline[1001,1480], iline[1075,1530]
    d1,d2,d3 = 0.004,0.025,0.012 # s, km, km
    f1,f2,f3 = 0.000,0.000,0.000 # s, km, km
    datbase = "pen"
    datfile = datDir+"pen.dat" # flat file of BIG_ENDIAN floats
    odtfile = datDir+"penx.dat" # simple file exported from OpendTect
    scale,scale1 = 2.0,2.0

def convertEndian(infile,outfile):
  x = zerofloat(n1*n2*n3)
  ais = ArrayInputStream(infile,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(x)
  ais.close()
  aos = ArrayOutputStream(outfile)
  aos.writeFloats(x)
  aos.close()

def display3d(datfile):
  ais = ArrayInputStream(datfile)
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  print "x min =",min(x)," max =",max(x)
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,x)
  ipg.setClips(-clip,clip)
  view = frame.getOrbitView()
  view.setScale(scale)
  view.setAxesScale(1.0,1.0,scale1)
  frame.setSize(1200,900)
  frame.setVisible(True)

def displaySlice2(datfile):
  k2 = 211
  j1 = 25
  m1 = 315
  ais = ArrayInputStream(datfile)
  x = zerofloat(n1,n2,n3)
  ais.readFloats(x)
  ais.close()
  y = zerofloat(m1,n3)
  for i3 in range(n3):
    copy(m1,j1,x[i3][k2],0,y[i3])
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(y)
  pv.setClips(-clip,clip)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 

"""
x3=x   x2=y   x1=z
0.034, 5.275, 0.100
5.412, 5.275, 1.357

  0 211  25
480 211 340 
"""
