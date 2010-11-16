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
D24 = LocalDiffusionKernel.Stencil.D24
D33 = LocalDiffusionKernel.Stencil.D33
D71 = LocalDiffusionKernel.Stencil.D71

#############################################################################

def main(args):
  benchStencils()
  #showStencils()

def benchStencils():
  stencils = [D22,D33,D71]
  n1,n2,n3 = 200,200,200
  x = randfloat(n1,n2,n3)
  x = applyBandPassFilter(x)
  plotH3(x)
  y = zerofloat(n1,n2,n3)
  au,av,aw = 1.0,0.0,0.0
  a1,a2,a3 = 25.0,25.0,0.0
  d = makeConstantTensors(n1,n2,n3,au,av,aw,a1,a2,a3)
  for stencil in stencils:
    ldk = LocalDiffusionKernel(stencil)
    sw = Stopwatch()
    sw.restart()
    while sw.time()<1.0: # warmup
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
  n1,n2,n3 = 15,15,15
  x = zerofloat(n1,n2,n3); x[n3/2][n2/2][n1/2] = 1.0
  au,av,aw = 1.0,0.0,0.0
  a1,a2,a3 = 45.0,45.0,0.0
  d = makeConstantTensors(n1,n2,n3,au,av,aw,a1,a2,a3)
  for stencil in stencils:
    ldk = LocalDiffusionKernel(stencil)
    h = zerofloat(n1,n2,n3)
    ldk.apply(d,x,h)
    #plotH3(h)
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

def makeSphereTensors(n1,n2,n3):
  c1,c2,c3 = 0,0,0
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  w1 = zerofloat(n1,n2,n3)
  w2 = zerofloat(n1,n2,n3)
  au = zerofloat(n1,n2,n3)
  av = zerofloat(n1,n2,n3)
  aw = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    x3 = i3-c3
    for i2 in range(n2):
      x2 = i2-c2
      for i1 in range(n1):
        x1 = i1-c1
        xu = x1*x1+x2*x2+x3*x3
        if xu!=0.0:
          xu = 1.0/sqrt(xu)
        else:
          xu = 1.0
        xw = x1*x1+x2*x2
        if xw!=0.0:
          xw = 1.0/sqrt(xw)
        else:
          xw = 1.0
        u1[i3][i2][i1] =  xu*x1
        u2[i3][i2][i1] =  xu*x2
        w1[i3][i2][i1] = -xw*x2
        w2[i3][i2][i1] =  xw*x1 
        au[i3][i2][i1] = 0.0
        av[i3][i2][i1] = 1.0
        aw[i3][i2][i1] = 1.0
  return EigenTensors3(u1,u2,w1,w2,au,av,aw)

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
  ip.setColorModel(ColorMap.JET)
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

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
