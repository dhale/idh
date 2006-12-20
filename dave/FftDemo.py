from math import *
from edu.mines.jtk.util import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.mosaic import *
from dave import *

True = 1
False = 0

def plot(x,y,z):
  panel = PlotPanel(3,1)
  panel.addSequence(0,0,x)
  panel.addSequence(1,0,y)
  panel.addSequence(2,0,z)
  frame = PlotFrame(panel)
  frame.setSize(950,800)
  frame.setVisible(True)
  return frame

lx = 31;  kx = 0
ly = 31;  ky = 0
lz = 11;  kz = -5
sx = max(0, kz)
sy = max(0,-kz)
lpad = lx+max(-kz,lz-1+kz)
nfft = FftReal.nfftSmall(lpad)
print "lpad =",lpad," nfft =",nfft
fft = FftReal(nfft)
x = Array.fillfloat(1.0,lx)
y = Array.reverse(Array.fillfloat(1.0,ly))
z = Array.zerofloat(lz)
xpad = Array.zerofloat(nfft+2)
ypad = Array.zerofloat(nfft+2)
Array.copy(lx,  0,x,sx,xpad)
Array.copy(ly,  0,y,sy,ypad)
fft.realToComplex(-1,xpad,xpad)
fft.realToComplex(-1,ypad,ypad)
xpad = Array.cconj(xpad)
zpad = Array.cmul(xpad,ypad)
fft.complexToReal(1,zpad,zpad)
fft.scale(lz,zpad)
Array.copy(lz,zpad,z)
#Conv.xcor(lx,kx,x,ly,ky,y,lz,kz,z)

frame = plot(x,y,z)
