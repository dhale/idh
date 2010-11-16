import sys
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from lss import *

D22 = LocalDiffusionKernel.Stencil.D22
D33 = LocalDiffusionKernel.Stencil.D33
D71 = LocalDiffusionKernel.Stencil.D71

#############################################################################

def main(args):
  benchStencils()
  #showStencils()

def benchStencils():
  stencils = [D22,D33,D71]
  n1,n2,n3 = 200,200,200
  x = sub(randfloat(n1,n2,n3),0.5)
  x = applyBandPassFilter(x)
  plotH3(x)
  y = zerofloat(n1,n2,n3)
  au,av,aw = 1.0,0.0,0.0
  a1,a2,a3 = 45.0,22.0,0.0
  d = makeConstantTensors(n1,n2,n3,au,av,aw,a1,a2,a3)
  for stencil in stencils:
    ldk = LocalDiffusionKernel(stencil)
    sw = Stopwatch()
    sw.restart()
    while sw.time()<0.1: # warmup
      zero(y)
      ldk.apply(d,x,y)
    n = 0
    sw.restart()
    while sw.time()<1.0:
      zero(y)
      ldk.apply(d,x,y)
      n += 1
    sw.stop()
    rate = int(n*n1*n2*n3*1.0e-3/sw.time())
    sumy = sum(y)
    print "stencil =",stencil," sum =",sumy," rate =",rate,"kilosamples/s"
    plotH3(y)

def showStencils():
  stencils = [D22,D33,D71]
  n1,n2,n3 = 101,101,101
  x = zerofloat(n1,n2,n3); x[n3/2][n2/2][n1/2] = 1.0
  x = applyBandPassFilter(x)
  au,av,aw = 1.0,0.0,0.0
  a1,a2,a3 = 45.0,22.0,0.0
  d = makeConstantTensors(n1,n2,n3,au,av,aw,a1,a2,a3)
  for stencil in stencils:
    ldk = LocalDiffusionKernel(stencil)
    h = zerofloat(n1,n2,n3)
    ldk.apply(d,x,h)
    plotA3(h)

def makeConstantTensors(n1,n2,n3,au,av,aw,a1,a2,a3):
  a1 *= PI/180
  a2 *= PI/180
  a3 *= PI/180
  c1,s1 = cos(a1),sin(a1)
  c2,s2 = cos(a2),sin(a2)
  c3,s3 = cos(a3),sin(a3)
  u1 = c1*c3-c2*s1*s3
  u2 = c2*c3*s1+c1*s3
  w1 = s2*s3
  w2 = -c3*s2
  u1 = fillfloat(u1,n1,n2,n3)
  u2 = fillfloat(u2,n1,n2,n3)
  w1 = fillfloat(w1,n1,n2,n3)
  w2 = fillfloat(w2,n1,n2,n3)
  au = fillfloat(au,n1,n2,n3)
  av = fillfloat(av,n1,n2,n3)
  aw = fillfloat(aw,n1,n2,n3)
  return EigenTensors3(u1,u2,w1,w2,au,av,aw,True)

def applyBandPassFilter(f):
  bpf = BandPassFilter(0.00,0.45,0.10,0.01)
  bpf.setExtrapolation(BandPassFilter.Extrapolation.ZERO_SLOPE)
  g = copy(f)
  bpf.apply(f,g)
  return g
 
#############################################################################
# plot

def plotH3(h):
  h = copy(h)
  print "amplitude min =",min(h)," max =",max(h)
  n1,n2,n3 = len(h[0][0]),len(h[0]),len(h)
  f1,f2,f3 = -(n1-1)/2,-(n2-1)/2,-(n3-1)/2
  s1,s2,s3 = Sampling(n1,1,f1),Sampling(n2,1,f2),Sampling(n3,1,f3)
  sf = SimpleFrame()
  ip = sf.addImagePanels(s1,s2,s3,h)
  #ip.setColorModel(ColorMap.JET)
  sf.orbitView.setScale(3.0)
  sf.setSize(900,900)

def plotA3(h):
  fft = Fft(h)
  fft.setCenter(True)
  fft.setPadding(100)
  ak = cabs(fft.applyForward(h))
  print "amplitude max =",max(ak)
  sk1 = fft.getFrequencySampling1()
  sk2 = fft.getFrequencySampling2()
  sk3 = fft.getFrequencySampling2()
  sf = SimpleFrame()
  ip = sf.addImagePanels(sk1,sk2,sk3,ak)
  ip.setClips(0,1)
  ip.setColorModel(ColorMap.JET)
  sf.orbitView.setScale(3.0)
  sf.setSize(900,900)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
